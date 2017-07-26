## Script: Analysis of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-12
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(reshape2)

acti_files <- readRDS(".\\Data\\actigraphy_header.rds")
actigraphy <- readRDS(".\\Data\\actigraphy_data.rds") 

results <- list()

## ------ Visualizing a single patient ------

d1 <- dplyr::filter(actigraphy, patient_ID == "Ryan_Opel") %>%
  select(., Time, Light, Activity, Day, DateAbbr) %>%
  mutate(Activity_Scale = Activity * (max(Light) / max(Activity)))

ggplot(d1, mapping = aes(x = Time)) +
  geom_line(aes(y = Light), color = "Orange") +
  geom_line(aes(y = Activity_Scale), color = "Black") +
  facet_grid(Day + DateAbbr ~ .) +
  scale_x_datetime(date_labels = "%I %p", date_minor_breaks = "1 hour") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * (max(d1$Activity) / max(d1$Light)),
                                         name = "Activity")) +
  ylab("Light")

## ------ Off Wrist Filtering ------

off_wrist <- as.data.frame.table(xtabs(~patient_ID + Off_Wrist + Day, data = actigraphy)) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ Off_Wrist, value.var = "Freq") %>%
  rename(Watch_On = `FALSE`, Watch_Off = `TRUE`, Day_number = Day) %>%
  mutate(Percent_Off = 100 * Watch_Off/(Watch_On + Watch_Off),
         Total_Epochs = Watch_On + Watch_Off) %>%
  dplyr::filter(complete.cases(.))

actigraphy %<>% merge(., off_wrist[,c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")] %>% 
               rename(Day = Day_number), by = c("patient_ID", "Day")) %>%
  dplyr::filter(Percent_Off < 10, Total_Epochs >= (720 / 2)) %>%
  arrange(patient_ID, DateTime)

## ------ Percent Time Off-Wrist ------

percent_off_wrist <- aggregate(cbind(Watch_Off, Total_Epochs) ~ patient_ID,
                               off_wrist, sum) %>%
  mutate(Percent_Off = 100 * Watch_Off / Total_Epochs)

ggplot(off_wrist, aes(Day_number, Percent_Off, group = patient_ID, color = patient_ID)) + 
  geom_line(size = 1) +
  facet_grid(patient_ID ~ .)

## ------ Daytime activity ratio (DAR) ------
# taken from http://www.neurology.org/content/88/3/268.long
# Daytime is defined as 7AM - 10 PM --> can be modified in 01_ETL.R script

DAR_parent <- aggregate(Activity ~ patient_ID + DAR_Period + Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ DAR_Period, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Day_number", "Day", "Night")) %>%
  mutate(DAR = 100 * (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  # Use `off_wrist` to verify that watch is on and recording for most of the day
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number")) #%>%
  # dplyr::filter(Percent_Off < 10, Total_Epochs >= (720 / 2))

DAR_child <- cbind(aggregate(DAR ~ patient_ID, DAR_parent, mean),
                   aggregate(DAR ~ patient_ID, DAR_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                    })[2]) %>%
  set_colnames(c("patient_ID", "DAR", "SEM")) %>%
  mutate(ymax = DAR + SEM,
         ymin = DAR - SEM)

ggplot(DAR_child, aes(x = patient_ID, y = DAR)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25) +
  labs(y = "Daytime Activity Ratio")

## ------ Nighttime sleep duration ------

# Need to decide which sleep determinant is the best in previous script
# Again, daytime is defined as 7AM - 10 PM --> can be modified in 01_ETL.R script

# `Line` is being used as a dummy variable
NSD_parent <- aggregate(Line ~ patient_ID + Sleep_Wake + Day + DAR_Period,
                        actigraphy, length) %>%
  reshape2::dcast(., formula = patient_ID + DAR_Period + Day ~ Sleep_Wake,
                  value.var = "Line") %>%
  rename(Day_number = Day) %>%
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number")) %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake)))

NSD_child <- cbind(aggregate(Percent_Sleep ~ patient_ID, NSD_parent, mean),
                   aggregate(Percent_Sleep ~ patient_ID, NSD_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                   })[2]) %>%
  set_colnames(c("patient_ID", "Percent_Sleep", "SEM")) %>%
  mutate(ymax = Percent_Sleep + SEM,
         ymin = Percent_Sleep - SEM)

ggplot(NSD_child, aes(x = patient_ID, y = Percent_Sleep)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25)

## ------ Percent Sleep/Social Jetlag ------

percent_sleep <- lapply(unique(actigraphy$patient_ID), function (ii) {
  df <- actigraphy[actigraphy$patient_ID == ii, ]
  spec <- as.data.frame.table(100 * prop.table(table(df$Date, df$Sleep_Wake), 1))
  sleep <- cbind(spec, factor(ii))
  
  return(sleep)
}) %>% 
  do.call("rbind", .) %>%
  data.frame(.) %>%
  set_colnames(c('Date', 'Stage', 'Sleep', 'patient_ID')) %>%
  filter(complete.cases(.), Stage == "Sleep") %>%
  mutate(Date = as.POSIXct(strptime(Date, format = "%F")),
         Weekend = weekdays(Date),
         Weekend = factor(ifelse(Weekend %in% c("Saturday", "Sunday"),
                                 "Weekend", "Weekday")))

percent_sleep_agg <- merge(aggregate(Sleep ~ patient_ID + Weekend, percent_sleep, mean),
                           aggregate(Sleep ~ patient_ID + Weekend, percent_sleep, opelr::sem) %>%
                             rename(SEM = Sleep),
                           by = c("patient_ID", "Weekend")) %>%
  mutate(ymax = Sleep + SEM,
         ymin = Sleep - SEM)

ggplot(percent_sleep, aes(patient_ID, Sleep)) + 
  geom_bar(aes(fill = Weekend), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 33)

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

# Plotting
ggplot(bedtime, aes(Bedtime_NoDate,  Day_of_Week, color = Weekend)) +
  geom_point() +
  scale_x_datetime(date_breaks = "2 hours",
                   date_labels = "%I:%M %p") +
  labs(x = "Hour") +
  facet_grid(patient_ID ~ .)

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

ggplot(average_bedtime, aes(DateTime, patient_ID, color = Status)) +
  geom_point() +
  scale_x_datetime(date_breaks = "2 hours",
                   date_labels = "%H:%M") +
  geom_errorbarh(aes(xmax = upper_Date, xmin = lower_Date), height = 0.1)

## ------ KNN Sleep Stage Prediction ------

## ------ Light Adherence ------

clean_light_adherence <- function(data, column, patient, date) {
  # Recursively smooth light adherence
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
  clean_light_adherence(actigraphy, "Light_15_10k", patient_dates$patient_ID[ii],
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
  mutate(Minutes = 60 * Bouts + 2 * (Lengths - Bouts))
