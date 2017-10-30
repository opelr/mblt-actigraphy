## Script: Analysis of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-12
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(zoo)

actigraphy <- readRDS(".\\Rmd\\Data\\actigraphy_filtered.rds") 

results <- list()

## ------ Basic Sleep Metrics ------

#' Creating a parent dataframe used to calculate TST, WASO, and SE

sleep_metrics_parent <- xtabs(~ patient_ID + Sleep_Bouts + Sleep_Acti_Smooth + Noon_Day,
                              data = actigraphy) %>%
  as.data.frame %>%
  dcast(., patient_ID + Noon_Day + Sleep_Bouts ~ Sleep_Acti_Smooth,
        value.var = "Freq") %>%
  filter(Sleep_Bouts == "Sleep_Bout", (Sleep + Wake) > 0) %>%
  mutate(Sleep_Minutes = Sleep * 2,
         WASO_Minutes = Wake * 2,
         SE = 100 * Sleep / (Wake + Sleep))

results$TST <- aggregate(Sleep_Minutes ~ patient_ID, sleep_metrics_parent, mean)

results$WASO <- aggregate(WASO_Minutes ~ patient_ID, sleep_metrics_parent, mean)

results$SE <- aggregate(SE ~ patient_ID, sleep_metrics_parent, mean)

## ------ Sleep Latency ------



## ------ Rhythm Stability ------

### Interdaily Stability
calculate_IS <- function(var, patient) {
  #' Calculate Interdaily Stability for a non-parametric actigraphy variable
  #' The extent to which the profiles of indivudal days resemble each other
  #' 
  #' Args:
  #'   var (str): Column name present in your 'actigraphy' dataframe
  #'   patient (int): Patient ID number
  #' 
  #' Returns:
  #'   Data frame containing IS for every unique patient in your actigraphy data frame.
  #'    Values range from 0 (Guassian noise) to 1 (perfect IS).
  #'   
  #' Reference:
  #'   Van Someren, Chronobiology International, 1999
  
  ## Helper fn
  fp <- function(x, y) {
    # Shorthand notation for 'as.formula(paste(.))'
    return(as.formula(paste(x, "~", y)))
  }
  
  pat_unique_days <- unique(actigraphy$Noon_Day[actigraphy$patient_ID == patient])
  pat_seq_days <- seq(pat_unique_days[1], pat_unique_days[length(pat_unique_days)], 1)
  
  IS <- lapply(3:length(pat_unique_days), function (ii) {
    dat <- actigraphy[actigraphy$patient_ID == patient &
                      actigraphy$Noon_Day %in% pat_unique_days[1]:pat_unique_days[ii], ]
    
    ## Create numerator
    hourly_mean <- merge(aggregate(fp(var, "patient_ID + Hour"), dat, mean) %>%
                           rename_(.dots=setNames(names(.), gsub(var, "Hourly_Mean", names(.)))),
                         aggregate(fp(var, "patient_ID"), dat, mean) %>%
                           rename_(.dots=setNames(names(.), gsub(var, "Total_Mean", names(.)))),
                         by = "patient_ID") %>%
      arrange(patient_ID, Hour) %>%
      mutate(Variance = (Hourly_Mean - Total_Mean) ** 2)
    
    numerator <- merge(aggregate(Variance ~ patient_ID, hourly_mean, sum),
                       aggregate(fp(var, "patient_ID"), dat, length) %>%
                         rename_(.dots=setNames(names(.), gsub(var, "n", names(.)))),
                       by = "patient_ID") %>%
      mutate(Numerator = n * Variance)
    
    ## Denominator
    total_mean <- merge(actigraphy[, c("patient_ID", "Epoch", var)] %>%
                          rename_(.dots=setNames(names(.), gsub(var, "Indiv_Data", names(.)))),
                        aggregate(fp(var, "patient_ID"), dat, mean) %>%
                          rename_(.dots=setNames(names(.), gsub(var, "Total_Mean", names(.)))),
                        by = "patient_ID") %>%
      arrange(patient_ID, opelr::numfac(Epoch)) %>%
      mutate(Variance = (Indiv_Data - Total_Mean) ** 2)
    
    denominator <- aggregate(Variance ~ patient_ID, total_mean, sum) %>%
      mutate(Denominator = Variance * 24)
    
    ## Combine
    IS <- merge(numerator[, c("patient_ID", "Numerator")],
                denominator[, c("patient_ID", "Denominator")],
                by = "patient_ID") %>%
      mutate(IS = Numerator / Denominator,
             Noon_Day = pat_unique_days[ii]) %>%
      select(patient_ID, Noon_Day, IS)
  }) %>%
    do.call("rbind", .)
  
  IS_blank <- data.frame(patient_ID = patient,
                         Noon_Day = pat_seq_days)
  
  IS %<>% merge(IS_blank, ., by = c("patient_ID", "Noon_Day"), all.x = T)
  
  return(IS)
}

### Intradaily Variability
calculate_IV <- function(var, patient) {
  #' Calculate Intradaily Variability for a non-parametric actigraphy variable
  #' Quantifies how fragmented the rhythm is relative to the overall variance
  #' 
  #' Args:
  #'   var (str): Column name present in your 'actigraphy' dataframe
  #'   patient (int): Patient ID number
  #' 
  #' Returns:
  #'   Data frame containing IV for every unique patient in your actigraphy data frame.
  #'     Values reach near zero for a perfect sine wave, around 2 for Gaussian noise,
  #'     and may be higher when a definite ultradian component with a period of
  #'     2h is present.
  #'   
  #' Reference:
  #'   Van Someren, Chronobiology International, 1999
  
  ## Helper fn
  fp <- function(x, y) {
    # Shorthand notation for 'as.formula(paste(.))'
    return(as.formula(paste(x, "~", y)))
  }
  
  pat_unique_days <- unique(actigraphy$Noon_Day[actigraphy$patient_ID == patient])
  pat_seq_days <- seq(pat_unique_days[1], pat_unique_days[length(pat_unique_days)], 1)
  
  IV <- lapply(3:length(pat_unique_days), function (ii) {
    dat <- actigraphy[actigraphy$patient_ID == patient &
                      actigraphy$Noon_Day %in% pat_unique_days[1]:pat_unique_days[ii], ]
    
    ## Numerator
    lagged_var <- dat[, c("patient_ID", "Epoch", var)] %>%
      split(., .[, "patient_ID"]) %>%
      lapply(., function(ii) {
        ii[[paste0(var, "_Lag")]] <- lag(ii[, var], 1)
        
        ii %<>%
          mutate(Lag_1_Variance = (.[, var] - .[, paste0(var, "_Lag")]) ** 2)
      }) %>%
      do.call('rbind', .) %>%
      set_rownames(1:nrow(.))
    
    numerator <- merge(aggregate(Lag_1_Variance ~ patient_ID, lagged_var, sum),
                       aggregate(Epoch ~ patient_ID, lagged_var, length),
                       by = "patient_ID") %>%
      mutate(Numerator = Lag_1_Variance * Epoch)
    
    ## Denominator
    total_mean <- merge(dat[, c("patient_ID", "Epoch", var)] %>%
                          rename_(.dots=setNames(names(.), gsub(var, "Indiv_Data", names(.)))),
                        aggregate(fp(var, "patient_ID"), dat, mean) %>%
                          rename_(.dots=setNames(names(.), gsub(var, "Total_Mean", names(.)))),
                        by = "patient_ID") %>%
      arrange(patient_ID, opelr::numfac(Epoch)) %>%
      mutate(Variance = (Indiv_Data - Total_Mean) ** 2)
    
    denominator <- merge(aggregate(Variance ~ patient_ID, total_mean, sum),
                         aggregate(Epoch ~ patient_ID, total_mean, length),
                         by = "patient_ID") %>%
      mutate(Denominator = Variance * (Epoch - 1))
    
    ## Combine
    IV <- merge(numerator[, c("patient_ID", "Numerator")],
                denominator[, c("patient_ID", "Denominator")],
                by = "patient_ID") %>%
      mutate(IV = Numerator / Denominator,
             Noon_Day = pat_unique_days[ii]) %>%
      select(patient_ID, Noon_Day, IV)
  }) %>%
    do.call("rbind", .)
  
  
  IV_blank <- data.frame(patient_ID = patient,
                         Noon_Day = pat_seq_days)
  
  IV %<>% merge(IV_blank, ., by = c("patient_ID", "Noon_Day"), all.x = T)
  
  return(IV)
}

### Domintant Rest Phase Onset
# DRPO <- function() {
#' Represents the clock time at which the 5-hr period of lowest activity in
#' the average 24-hr pattern started
#   
# }

results$interdaily_stab <- lapply(unique(actigraphy$patient_ID), function(ii) {
    calculate_IS("Activity", ii)
  }) %>%
  do.call("rbind", .)

results$intradaily_var <- lapply(unique(actigraphy$patient_ID), function(ii) {
  calculate_IV("Activity", ii)
}) %>%
  do.call("rbind", .)

## ------ Relative Amplitude (RA), M10, and L5 ------

calc_rest_phase <- function(len, fn = 'max') {
  hours <- unique(actigraphy$Hour)[1:(length(unique(actigraphy$Hour)) - (len - 1))]
  
  rest_phase <- lapply(hours, function (i) {
    out <- aggregate(Activity ~ patient_ID + Day,
                     dplyr::filter(actigraphy, Hour >= i, Hour <= i + (len - 1)),
                     sum) %>%
      rename(Activity_Sum = Activity) %>%
      mutate(Hours = paste0(i, "-", i + (len - 1)),
             Activity_Mean = Activity_Sum / len)
  }) %>%
    do.call("rbind", .)
  
  rp2 <- aggregate(Activity_Mean ~ patient_ID + Hours, rest_phase, mean)
  
  rp3 <- aggregate(Activity_Mean ~ patient_ID, rp2, function(j) {f <- get(fn); f(j)}) %>%
    merge(., rp2, by = c("patient_ID", "Activity_Mean"), all_x = T) %>%
    rename_(.dots=setNames(names(.), tolower(gsub("Hours",
                                                  paste0(fn, len, "Hours"),
                                                  names(.)))))
  
  return(rp3)
}

results$RA <- merge(calc_rest_phase(10, 'max') %>%
                      rename(M10_Activity = activity_mean),
                    calc_rest_phase(5, 'min')  %>%
                      rename(L5_Activity = activity_mean),
                    by = "patient_id") %>%
  rename(patient_ID = patient_id) %>%
  mutate(RA = (M10_Activity - L5_Activity) / (M10_Activity + L5_Activity))

## ------ Percent Time Off-Wrist ------

## TODO:\\ This is no longer accurate, because each "Watch_Off" counts is 3 hours

percent_off_wrist <- aggregate(cbind(Watch_Off, Total_Epochs) ~ patient_ID,
                               off_wrist[opelr::numfac(off_wrist$Day) <= 28, ], sum) %>%
  mutate(Percent_Off = 100 * Watch_Off / Total_Epochs)

results$percent_off_wrist <- percent_off_wrist

## ------ Daytime activity ratio (DAR) ------
# taken from http://www.neurology.org/content/88/3/268.long

### Daytime is defined as 6:30 AM - 11 PM --> can be modified in 01_ETL.R script
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

### Daytime is defined as 7:00 AM - 7 PM
DAR_77_parent <- aggregate(Activity ~ patient_ID + DAR_77 + Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ DAR_77, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Day_number", "Day", "Night")) %>%
  mutate(DAR = 100 * (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number"))

DAR_77_child <- cbind(aggregate(DAR ~ patient_ID, DAR_77_parent, mean),
                   aggregate(DAR ~ patient_ID, DAR_77_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                   })[2]) %>%
  set_colnames(c("patient_ID", "DAR", "SEM")) %>%
  mutate(ymax = DAR + SEM,
         ymin = DAR - SEM)

results$DAR_77 <- DAR_77_child

## ------ Nighttime sleep duration ------
# `Line` is being used as a dummy variable

NSD_parent <- aggregate(Line ~ patient_ID + Sleep_Thresh_Smooth + SunPeriod + Day,
                        actigraphy, length) %>%
  reshape2::dcast(., formula = patient_ID + SunPeriod + Day ~ Sleep_Thresh_Smooth,
                  value.var = "Line") %>%
  opelr::replace_na(., 0) %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake))) %>%
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

percent_sleep <- aggregate(Sleep_Thresh_Smooth ~ patient_ID + MBLT_Group + Date + Day, actigraphy, table) %>%
  as.data.frame.table(.) %>%
  filter(Var2 == "patient_ID", Freq.Day != 1) %>%
  select(-Var1, -Var2, -Freq.Day) %>%
  set_colnames(c("patient_ID", "MBLT_Group", "Date", "Sleep", "Wake")) %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake)),
         Date = as.POSIXct(strptime(Date, format = "%F")),
         DayOfWeek = weekdays(Date),
         Weekend = factor(ifelse(DayOfWeek %in% c("Saturday", "Sunday"),
                                 "Weekend", "Weekday"))) %>%
  arrange(patient_ID, Date)

# results$percent_sleep <- percent_sleep

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
                 sleep_column = "Sleep_Thresh_Smooth", status = status,
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
results$patient_dates <- patient_dates

## ------ Median Light Exposure -------

median_light_exposure <- cbind(aggregate(Light ~ patient_ID + Hour,
                                      actigraphy, median),
                            aggregate(Light ~ patient_ID + Hour,
                                      actigraphy, opelr::sem)[3] %>%
                              rename(SEM = Light)) %>%
  mutate(ymax = Light + SEM,
         ymin = Light - SEM)

results$median_light_exposure <- median_light_exposure

## ------ Average Light Exposure during Lightbox Period ------

actigraphy %<>% mutate(Prescribed_Time = Hour %in% 6:10)

avg_light_exposure <- cbind(aggregate(Light ~ patient_ID + Prescribed_Time,
                                         actigraphy, mean),
                               aggregate(Light ~ patient_ID + Prescribed_Time,
                                         actigraphy, opelr::sem)[3] %>%
                                 rename(SEM = Light)) %>%
  mutate(ymax = Light + SEM,
         ymin = Light - SEM)

results$avg_light_exposure_ <- avg_light_exposure

ggplot(avg_light_exposure, aes(patient_ID, Light, group = Prescribed_Time)) +
  geom_bar(aes(fill = Prescribed_Time), stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin, group = Prescribed_Time), 
                width = 0.5, position = position_dodge(0.9)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Light Exposure", x = "")

## ------ Sleep Duration vs. Light/Activity per day ------

compare <- merge(percent_sleep, aggregate(Light ~ patient_ID + Date + Day, actigraphy, median),
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
aggregate(Date ~ patient_ID, patient_dates, length)

## ------ Save RDS ------

saveRDS(results, ".\\Rmd\\Data\\results.rds")
