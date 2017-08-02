## Script: Analysis of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-12
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(reshape2)

acti_files <- readRDS(".\\Rmd\\Data\\actigraphy_header.rds")
actigraphy_static <- readRDS(".\\Rmd\\Data\\actigraphy_data.rds") 

results <- list()

## ------ Off Wrist Filtering ------

off_wrist <- as.data.frame.table(xtabs(~patient_ID + Off_Wrist + Day, data = actigraphy_static)) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ Off_Wrist, value.var = "Freq") %>%
  rename(Watch_On = `FALSE`, Watch_Off = `TRUE`, Day_number = Day) %>%
  mutate(Percent_Off = 100 * Watch_Off/(Watch_On + Watch_Off),
         Total_Epochs = Watch_On + Watch_Off) %>%
  dplyr::filter(complete.cases(.))

actigraphy <-  merge(actigraphy_static,
                     off_wrist[,c("patient_ID", "Day_number",
                                  "Percent_Off", "Total_Epochs")] %>% 
               rename(Day = Day_number), by = c("patient_ID", "Day")) %>%
  dplyr::filter(Percent_Off < (2 / 24), Total_Epochs >= (720 / 2)) %>%
  arrange(patient_ID, DateTime) #%>%
  # filter(patient_ID != "Rocko_Jones", patient_ID != "Carolyn_Jones",
  #        patient_ID != "Lew_Andrews")

saveRDS(actigraphy, ".\\Rmd\\Data\\actigraphy_filtered.rds")

## ------ Percent Time Off-Wrist ------

percent_off_wrist <- aggregate(cbind(Watch_Off, Total_Epochs) ~ patient_ID,
                               off_wrist, sum) %>%
  mutate(Percent_Off = 100 * Watch_Off / Total_Epochs)

results$percent_off_wrist <- percent_off_wrist

## ------ Daytime activity ratio (DAR) ------
# taken from http://www.neurology.org/content/88/3/268.long
# Daytime is defined as 6:30 AM - 11 PM --> can be modified in 01_ETL.R script

DAR_parent <- aggregate(Activity ~ patient_ID + DAR_Period + Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ DAR_Period, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Day_number", "Day", "Night")) %>%
  mutate(DAR = 100 * (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number"))

DAR_child <- cbind(aggregate(DAR ~ patient_ID, DAR_parent, mean),
                   aggregate(DAR ~ patient_ID, DAR_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                    })[2]) %>%
  set_colnames(c("patient_ID", "DAR", "SEM")) %>%
  mutate(ymax = DAR + SEM,
         ymin = DAR - SEM)

results$DAR <- DAR_child

## ------ Nighttime sleep percentage ------
# `Line` is being used as a dummy variable

NSD_parent <- aggregate(Line ~ patient_ID + Sleep_Smooth + SunPeriod + Day,
                        actigraphy, length) %>%
  reshape2::dcast(., formula = patient_ID + SunPeriod + Day ~ Sleep_Smooth,
                  value.var = "Line") %>%
  mutate(Wake = ifelse(is.na(Wake), 0, Wake),
         Sleep = ifelse(is.na(Sleep), 0, Sleep),
         Percent_Sleep = 100 * (Sleep/(Sleep + Wake))) %>%
  filter(SunPeriod == "Night")

NSD_child <- cbind(aggregate(Percent_Sleep ~ patient_ID, NSD_parent, mean),
                   aggregate(Percent_Sleep ~ patient_ID, NSD_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                   })[2]) %>%
  set_colnames(c("patient_ID", "Percent_Sleep", "SEM")) %>%
  mutate(ymax = Percent_Sleep + SEM,
         ymin = Percent_Sleep - SEM)

results$NSD <- NSD_child

## ------ Percent Sleep/Social Jetlag ------

social_jetlag <- aggregate(Sleep_Smooth ~ patient_ID + Date + Day, actigraphy, table) %>%
  as.data.frame.table(.) %>%
  filter(Var2 == "patient_ID", Freq.Day != 1) %>%
  select(-Var1, -Var2, -Freq.Day) %>%
  set_colnames(c("patient_ID", "Date", "Sleep", "Wake")) %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake)),
         Date = as.POSIXct(strptime(Date, format = "%F")),
         Weekend = weekdays(Date),
         Weekend = factor(ifelse(Weekend %in% c("Saturday", "Sunday"),
                                 "Weekend", "Weekday")))

results$social_jetlag <- social_jetlag

## ------ Bedtime ------

get_time <- function (time) {
  # Converts POSIX to decimal
  # 
  # Args:
  #   time (POSIX)
  # 
  # Returns:
  #   Time as a decimal
  
  if(!any(class(time) == "POSIXt"))
    stop("time must be POSIX")
  
  tim <- as.numeric(time - trunc(time, "days"))
  return(tim)
}

find_bedtime <- function(noon_day, ID, sleep_column, status = "bed", data) {
  # Bedtime/Waketime
  # 
  # Args:
  #   noon_day (int): Which day (noon-noon) we're investiagting
  #   ID (str): Patient name
  #   sleep_column (str): Which column are we running the function on
  #   status (str): "bed" or "wake" for bed- and wake-times
  #   data (df): Data frame to iterate upon
  # 
  # Returns:
  #   Bed- or Wake-time as POSIXct
  
  if (! toupper(status) %in% c("BED", "WAKE"))
    stop("'status' must be either 'bed' or 'wake'")
  
  df <- data %>%
    dplyr::filter(patient_ID == ID, Noon_Day == noon_day)
  
  sleep <- data.frame(unclass(rle(as.character(df[, sleep_column])))) %>%
    mutate(
      Index = cumsum(c(1, lengths[-length(lengths)])),
      Time = df[, "DateTime"][Index]
    )
  
  max_sleep <- which(sleep$lengths == max(sleep$lengths[sleep$values == "Sleep"])) %>%
    .[. %in% which(sleep$values == "Sleep")]
  
  bed_time <- sleep$Time[max_sleep]
  wake_time <- sleep$Time[max_sleep + 1]
  
  if (toupper(status) == "BED") {
    return(bed_time)
  } else {
    return(wake_time)
  }
}

get_bedtime_loop <- function(data, status) {
  # lapply 'find_bedtime' function for 'bedtime' data.frame
  # 
  # Args:
  #   data (data.frame)
  #   status (str): Bed or Wake, same as find_bedtime
  # 
  # Returns:
  #   Vector of POSIX times
  
  lapply(rownames(data), function(ii) {
    ii <- as.integer(ii)
    
    find_bedtime(data[, "Noon_Day"][ii], ID = data[, "patient_ID"][ii],
                 sleep_column = "Sleep_Smooth", status = status,
                 data = actigraphy)
  }) %>%
    do.call("c", .)
}
 
bedtime <- as.data.frame(table(actigraphy$patient_ID, actigraphy$Noon_Day)) %>%
  set_colnames(c("patient_ID", "Noon_Day", "Freq")) %>%
  dplyr::filter(Freq == 720) %>% 
  arrange(patient_ID, Noon_Day) %>%
  mutate(Bedtime = get_bedtime_loop(., 'bed'),
         Waketime = get_bedtime_loop(., 'wake'),
    Day_of_Week = weekdays(Bedtime),
    Weekend = factor(ifelse(Day_of_Week %in% c("Saturday", "Sunday"),
                            "Weekend", "Weekday")),
    Bedtime_Decimal = get_time(Bedtime),
    Waketime_Decimal = get_time(Waketime),
    Bedtime_NoDate = as.POSIXct(strptime(substr(Bedtime, 12, 20), format = "%T")),
    Waketime_NoDate = as.POSIXct(strptime(substr(Waketime, 12, 20), format = "%T"))) %>%
  dplyr::filter(complete.cases(.))

is_noon_time <- aggregate(DateTime ~ patient_ID + Noon_Day, actigraphy, min) %>%
  rename(Prev_Day = DateTime) %>%
  mutate(IsNoon = substr(Prev_Day, 12, 20) == "12:00:00") %>%
  dplyr::filter(IsNoon == T) %>%
  select(-IsNoon)

bedtime %<>%
  merge(., is_noon_time, by = c("patient_ID", "Noon_Day")) %>%
  mutate(Bedtime_from_Noon = as.integer(difftime(Bedtime, Prev_Day, units = "mins")),
         Waketime_from_Noon = as.integer(difftime(Waketime, Prev_Day, units = "mins")))

## ------ Average Bed/Wake Times ------

# Need to calculate bed and wake time from previous noon -- 
average_bedtime_SEM <- aggregate(cbind(Bedtime_from_Noon, Waketime_from_Noon) ~ patient_ID, 
          bedtime, opelr::sem) %>%
  melt(., id.vars = "patient_ID") %>%
  rename(Status = variable, SEM = value) %>%
  mutate(Status = gsub("time_from_Noon", "", Status))

average_bedtime <- aggregate(cbind(Bedtime_from_Noon, Waketime_from_Noon) ~ patient_ID,
                             bedtime, mean) %>%
  melt(., id.vars = "patient_ID") %>%
  rename(Status = variable, Time = value) %>%
  mutate(Status = gsub("time_from_Noon", "", Status)) %>%
  merge(., average_bedtime_SEM, by = c("patient_ID", "Status")) %>%
  mutate(upper = Time + SEM,
         lower = Time - SEM,
         Template = as.POSIXct(strptime("2017-01-01 12:00:00", format = "%F %T")),
         DateTime = Template + seconds(Time * 60),
         upper_Date = DateTime + seconds(SEM * 60),
         lower_Date = DateTime - seconds(SEM * 60))

results$avg_bedtime <- average_bedtime

## ------ Light Adherence ------

clean_adherence <- function(data, column, patient, date) {
  # Recursively smooth light/activity adherence
  # 
  # Args:
  #   data (df): Data frame to iterate upon
  #   column (str): Name of column to call smoothing functions on
  #   patient (str): Patient name/ID
  #   date (str): Adherence per day
  # 
  # Returns:
  #   RLE DF of TRUE and FALSE values
  
  # Convert RLE object to data.frame
  get_rle_df <- function(x) {
    data.frame(Values = rle(as.character(x))$values,
               Lengths = rle(as.character(x))$lengths)
  }
  # Convert RLE data.frame to a vector
  rle_df_to_vector <- function(df) {
    rep(df$Values, df$Lengths)
  }
  
  post_process_light <- function(df) {
    # Copy data for comparison
    df_old <- df
    
    for (ii in 1:nrow(df)) {
      if (is.na(df$Values[ii])) {
        df$Values[ii] <- NA
      } else if (df$Values[ii] == "FALSE") {
        if (ii == 1) next
        
        # 2) Sleep of 1, 3, or 4 minutes was rescored as wake if it preceded
        #    at least 4(2), 10(5), or 15(7) minutes of wake;
        # 3) Sleep of 6 or 10 minutes surrounded by at least 10 or 20 minutes
        #    of wake was rescored as wake.
        if ((df$Lengths[ii] == 1 & df$Lengths[ii - 1] >= 2) | 
            (df$Lengths[ii] == 2 & df$Lengths[ii - 1] >= 6) |
            (df$Lengths[ii] %in% c(3, 4) & df$Lengths[ii - 1] >= 5 & 
             df$Lengths[ii + 1] >= 5) |
            (df$Lengths[ii] == 5 & df$Lengths[ii - 1] >= 10 &
             df$Lengths[ii + 1] >= 10)) {
          df$Values[ii] <- "TRUE"
        }
      }
    }
    
    if (identical(df_old, df)) {
      return(df)
    } else {
      post_process_light(get_rle_df(rle_df_to_vector(df)))
    }
  }
  
  vec <- data[data$patient_ID == patient & data$Date == date, column]
  rle_df <- get_rle_df(vec)
  
  out <- post_process_light(rle_df) %>%
    cbind(., patient, date)
}

patient_dates <- table(actigraphy$patient_ID, actigraphy$Date) %>%
  as.data.frame(.) %>%
  filter(Freq > 0) %>%
  set_colnames(c("patient_ID", "Date", "Epochs")) %>%
  arrange(patient_ID, Date)

light_adherence_parent <- lapply(1:nrow(patient_dates), function(ii) {
  clean_adherence(actigraphy, "Light_15_10k", patient_dates$patient_ID[ii],
                        patient_dates$Date[ii])
}) %>%
  do.call("rbind", .) %>%
  dplyr::filter(complete.cases(.)) %>%
  rename(patient_ID = patient, Date = date)
  
light_adherence <- merge(aggregate(Lengths ~ patient_ID + Date + Values,
                                   light_adherence_parent, sum),
      aggregate(Lengths ~ patient_ID + Date + Values, light_adherence_parent,
                length) %>%
        rename(Bouts = Lengths),
      by = c("patient_ID", "Date", "Values")) %>%
  filter(Values == TRUE) %>%
  arrange(patient_ID, Date) %>%
  mutate(Minutes = 30 * Bouts + 2 * (Lengths - Bouts),
         ABL = Minutes/Bouts)

results$light_adherence <- light_adherence

## ------ Average Light Exposure -------

median_light_exposure <- cbind(aggregate(Light ~ patient_ID + Hour,
                                      actigraphy, median),
                            aggregate(Light ~ patient_ID + Hour,
                                      actigraphy, opelr::sem)[3] %>%
                              rename(SEM = Light)) %>%
  mutate(ymax = Light + SEM,
         ymin = Light - SEM)

results$median_light_exposure <- median_light_exposure

## ------ Sleep Duration vs. Light/Activity per day ------

compare <- merge(social_jetlag, aggregate(Light ~ patient_ID + Date + Day, actigraphy, median),
                 by = c("patient_ID", "Date")) %>%
  merge(., aggregate(Activity ~ patient_ID + Date, actigraphy, median),
        by = c("patient_ID", "Date")) %>%
  mutate(Hours = (Percent_Sleep / 100) * 24)

results$sleep_dur_light_activity_per_day <- compare

## ----- Day Radial Plots -------

radial_individuals <- aggregate(cbind(Light, Activity) ~ patient_ID + Hour, actigraphy, mean) %>%
  melt(., id.vars = c("patient_ID", "Hour")) %>%
  rename(Metric = variable, Value = value) %>%
  arrange(patient_ID, Metric, Hour) %>%
  mutate(Time = as.POSIXct(paste0(sprintf("%02d", Hour), ":00:00"), format = "%T"))

radial_means <- aggregate(Value ~ Metric + Hour, radial_individuals, mean) %>%
  mutate(patient_ID = "mean",
         Time = as.POSIXct(paste0(sprintf("%02d", Hour), ":00:00"), format = "%T")) %>%
  select(one_of(colnames(radial_individuals)))

radial_intermediary <- rbind(radial_individuals, radial_means)

radial_24 <- dplyr::filter(radial_intermediary, Hour == 0) %>%
  mutate(Hour = 24,
         Time = as.POSIXct(paste0(sprintf("%02d", Hour), ":00:00"), format = "%T"))

radial <- rbind(radial_intermediary, radial_24)
  
results$radial_plots <- radial

## ------ Activity EE Metrics ------

### Consolidated Bouts

activity_consolidation_parent <- lapply(1:nrow(patient_dates), function(ii) {
  get_rle_df <- function(x) {
    data.frame(Values = rle(as.character(x))$values,
               Lengths = rle(as.character(x))$lengths)
  }
  x <- actigraphy$Activity_15_2k[actigraphy$patient_ID == patient_dates$patient_ID[ii] &
                            actigraphy$Date == patient_dates$Date[ii]]
  
  cbind(get_rle_df(x), patient_dates$patient_ID[ii], patient_dates$Date[ii])
}) %>%
  do.call("rbind", .) %>%
  dplyr::filter(complete.cases(.)) %>%
  set_colnames(c("Values", "Lengths", "patient_ID", "Date"))

activity_consolidation <- merge(aggregate(Lengths ~ patient_ID + Date + Values,
                                          activity_consolidation_parent, sum),
                         aggregate(Lengths ~ patient_ID + Date + Values, activity_consolidation_parent,
                                   length) %>%
                           rename(Bouts = Lengths),
                         by = c("patient_ID", "Date", "Values")) %>%
  filter(Values == TRUE) %>%
  arrange(patient_ID, Date) %>%
  mutate(Minutes = 30 * Bouts + 2 * (Lengths - Bouts),
         ABL = Minutes/Bouts)

results$activity_consolidation <- activity_consolidation

### 

## ------ Basic Sleep Metrics ------
### TST

TST <- merge(aggregate(Sleep ~ patient_ID, social_jetlag, mean),
             aggregate(Sleep ~ patient_ID, social_jetlag, opelr::sem) %>%
               rename(SEM = Sleep),
             by = "patient_ID") %>%
  mutate(TST_minutes = Sleep * 2,
         SEM = SEM * 2,
         ymax = TST_minutes + SEM,
         ymin = TST_minutes - SEM)

results$TST <- TST

### WASO
ii <- actigraphy[actigraphy$patient_ID == "Sophie_Lee" & actigraphy$Noon_Day == 4, ]

rle(as.character(ii$Sleep_Smooth))
rle(as.character(ii$Sleep_Smooth))


### Sleep Latency



## ------ Save RDS ------

saveRDS(results, ".\\Rmd\\Data\\results.rds")
