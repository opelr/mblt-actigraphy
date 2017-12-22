## Script: Analysis of Actiwatch Data
## Project: mblt-actigraphy
## Author: Ryan Opel -- @opelr
## Date: 2017-11-15
## Version: 0.3.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(zoo)

source("Scripts/Functions.R")

actigraphy <- readRDS(".\\Rmd\\Data\\actigraphy_filtered.rds")
results <- list()

## ------ Basic Sleep Metrics ------

#' Creating a parent dataframe used to calculate TST, WASO, and SE

results$sleep_metrics <- xtabs(~ patient_ID + Sleep_Bouts + Sleep_Acti_Smooth + Noon_Day,
                               data = actigraphy) %>%
  as.data.frame %>%
  dcast(., patient_ID + Noon_Day + Sleep_Bouts ~ Sleep_Acti_Smooth,
        value.var = "Freq") %>%
  filter(Sleep_Bouts == "Sleep_Bout", (Sleep + Wake) > 0) %>%
  mutate(Sleep_Minutes = Sleep * 2,
         WASO_Minutes = Wake * 2,
         SE = 100 * Sleep / (Wake + Sleep))

## ------ Sleep Latency ------



## ------ Activity Rhythm Stability ------

### Domintant Rest Phase Onset
# DRPO <- function() {
#' Represents the clock time at which the 5-hr period of lowest activity in
#' the average 24-hr pattern started
#   
# }

IV_IS_combs <- expand.grid(c("Activity", "Light"), #, "Sleep_Thresh_Smooth_Int"
                           c("moving", "expanding"), c(3, 7)) %>% 
  rename(var = Var1, window = Var2, size = Var3) %>%
  arrange(var, window, size)

for (j in 1L:nrow(IV_IS_combs)) {
  var <- IV_IS_combs$var[j] %>% as.character
  window <- IV_IS_combs$window[j] %>% as.character
  size <- IV_IS_combs$size[j]
  
  IV_title <- paste(var, "IV", window, size, sep = "_")
  IS_title <- paste(var, "IS", window, size, sep = "_")
  
  IV <- lapply(unique(actigraphy$patient_ID), function(ii) {
    calculate_IV(var, ii, window, size, actigraphy)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IV", IV_title, names(.))))
  
  IS <- lapply(unique(actigraphy$patient_ID), function(ii) {
    calculate_IS(var, ii, window, size, actigraphy)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IS", IS_title, names(.))))
  
  results[[IV_title]] <- IV
  results[[IS_title]] <- IS
}

## ------ Relative Amplitude (RA), M10, and L5 ------

results$RA <- merge(calc_rest_phase(10, 'max') %>%
                      rename(M10_Activity = Activity_Mean),
                    calc_rest_phase(5, 'min')  %>%
                      rename(L5_Activity = Activity_Mean),
                    by = c("patient_ID", "Noon_Day")) %>%
  mutate(RA = (M10_Activity - L5_Activity) / (M10_Activity + L5_Activity))

## ------ Activity Consolidation ------

results$activity_consolidation <- calc_activity_consolidation("Activity_15_2k") %>%
  select(-Values, -Lengths) %>%
  rename(Activity_Consol_Bouts = Bouts,
         Activity_Consol_Minutes = Minutes,
         Activity_Consol_ABL = ABL)

## ------ Light Adherence ------

results$light_adherence <- calc_light_adherence("Light_15_5k_smooth") %>%
  select(-Values, -Lengths) %>%
  rename(Light_Adherent_Bouts = Bouts,
         Light_Adherent_Minutes = Minutes,
         Light_Adherent_ABL = ABL)

## ------ Daytime activity ratio (DAR) ------
# taken from http://www.neurology.org/content/88/3/268.long

#' DAR == 12-hour period based on individual's self-reported chronotype
#' DAR_Period == 6:30 AM - 11 PM
#' DAR_77 == 7AM - 7PM
#' Additional DAR variables can be added in 01_ETL.R script

results$DAR <- aggregate(Activity ~ patient_ID + DAR + Noon_Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Noon_Day ~ DAR, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Noon_Day", "Day", "Night")) %>%
  mutate(DAR = (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  arrange(patient_ID, Noon_Day)

## ------ Total Activity and Light Exposure ------

results$total_activity_light <- aggregate(cbind(Activity, Light) ~ patient_ID + Noon_Day,
                                  actigraphy, sum)

## ------ Average Light Exposure By Bins ------

# actigraphy$

## ------ Major Metrics by Patient by Day ------

results$patient_catalog <- xtabs(~ patient_ID + Noon_Day + Group + Date + watch_ID, actigraphy) %>%
  as.data.frame %>%
  filter(Freq > 0) %>%
  mutate(Group = factor(Group, levels = c("C", "CCT"),
                        labels = c("Ctrl", "CCT"))) %>%
  select(-Freq) %>%
  arrange(patient_ID, Noon_Day) %>%
  split(., .[, "watch_ID"]) %>%
  lapply(., function(i) i[seq(1, nrow(i), 2), ]) %>%
  do.call("rbind", .) %>%
  arrange(patient_ID, Noon_Day)

major_results <- results$sleep_metrics %>%
          select(-Sleep_Bouts, -Sleep, -Wake) 

for (i in 2:length(results)) {
  major_results <- merge(major_results, results[i][[1]],
                         by = c("patient_ID", "Noon_Day"),
                         all.x = T)
}

major_results %<>% arrange(patient_ID, Noon_Day) %>%
  mutate(Noon_Day = as.numeric(as.character(Noon_Day))) %>%
  merge(., xtabs(~ patient_ID + Group_PTSD + MBLT_Group, actigraphy) %>%
          as.data.frame() %>%
          filter(Freq > 0) %>%
          select(-Freq),
        by = "patient_ID")

results$major_results <- major_results %>%
  arrange(patient_ID)
rm(major_results)

## ------ Save RDS ------

saveRDS(results, ".\\Rmd\\Data\\results.rds")

## ----------------------------------------------------------------------------

## ------ Slopegraphs: First and Last Week Averages ------

FL_results <- results$major_results %>%
  mutate(Week = floor((Noon_Day - 1)/7) + 1) %>%
  aggregate(cbind(Activity_IV_expanding_3, Activity_IS_expanding_3,
                  Activity_IV_moving_3, Activity_IS_moving_3,
                  Activity_IV_expanding_7, Activity_IS_expanding_7,
                  Activity_IV_moving_7, Activity_IS_moving_7, RA, DAR,
                  Sleep_Minutes, WASO_Minutes, SE) ~ patient_ID + Week,
            ., mean, na.rm = T, na.action = NULL)

## Determine Week numbers
weeks <- aggregate(Week ~ patient_ID, FL_results, max)

# FIXME: Need to fix PTSD, light definitions later
FL_results_FL <-  split(FL_results, FL_results$patient_ID) %>%
  lapply(., function (ii) {
    ii %<>%
      filter(Week %in% c(1, weeks$Week[weeks$patient_ID == unique(.[, "patient_ID"])]))
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1L:nrow(.)) %>%
  mutate(Week_Name = factor(ifelse(Week == 1, "First", "Last")),
         PTSD = ifelse(patient_ID %in% c(3, 8, 9, 10, 11, 13, 15, 16), "PTSD", "No PTSD"),
         Light = ifelse(as.numeric(as.character(patient_ID)) < 10, "Control", "Light"))


## First and Second Weeks
FL_results_F2L2 <-  split(FL_results, FL_results$patient_ID) %>%
  lapply(., function (ii) {
    ii %<>%
      filter(Week %in% c(1, 2, weeks$Week[weeks$patient_ID == unique(.[, "patient_ID"])],
                         weeks$Week[weeks$patient_ID == unique(.[, "patient_ID"])] - 1))
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1L:nrow(.)) %>%
  mutate(Week_Name = factor(ifelse(Week %in% c(1, 2), "First", "Last"))) %>%
  select(-Week) %>%
  aggregate(. ~ patient_ID + Week_Name, ., mean) %>%
  mutate(PTSD = ifelse(patient_ID %in% c(3, 8, 9, 10, 11, 13, 15, 16), "PTSD", "No PTSD"),
         Light = ifelse(as.numeric(as.character(patient_ID)) < 10, "Control", "Light")) %>%
  arrange(patient_ID, Week_Name)

write.csv(FL_results_F2L2, "First2_Last2_Weeks.csv")

## Most Light Exposure
max_light_days_7 <- lapply(unique(actigraphy$patient_ID), function(ii) {
  light_week(ii, 7, fn = 'max')
}) %>%
  do.call("rbind", .) %>%
  rename(Max_Days = Days)

min_light_days_7 <- lapply(unique(actigraphy$patient_ID), function(ii) {
  light_week(ii, 7, fn = 'min')
}) %>%
  do.call("rbind", .) %>%
  rename(Min_Days  = Days)

max_light_days_14 <- lapply(unique(actigraphy$patient_ID), function(ii) {
  light_week(ii, 14, fn = 'max')
}) %>%
  do.call("rbind", .) %>%
  rename(Max_Days = Days)

min_light_days_14 <- lapply(unique(actigraphy$patient_ID), function(ii) {
  light_week(ii, 14, fn = 'min')
}) %>%
  do.call("rbind", .) %>%
  rename(Min_Days  = Days)

light_days_7 <- merge(max_light_days_7, min_light_days_7, by = 'patient_ID')
light_days_14 <- merge(max_light_days_14, min_light_days_14, by = 'patient_ID')

FL_results_light_days_7 <- split(results$major_results, results$major_results$patient_ID) %>%
  lapply(., function(ii) {
    pat_ID <- unique(ii$patient_ID)
    max_day_start <- gsub("-(.*)", "", light_days_7$Max_Days[light_days_7$patient_ID == pat_ID]) %>% as.numeric
    max_day_end <- gsub("(.*)-", "", light_days_7$Max_Days[light_days_7$patient_ID == pat_ID]) %>% as.numeric
    min_day_start <- gsub("-(.*)", "", light_days_7$Min_Days[light_days_7$patient_ID == pat_ID]) %>% as.numeric
    min_day_end <- gsub("(.*)-", "", light_days_7$Min_Days[light_days_7$patient_ID == pat_ID]) %>% as.numeric
    
    ii %<>%
      filter(Noon_Day %in% c(seq(max_day_start, max_day_end), seq(min_day_start, min_day_end))) %>%
      mutate(Min_Max_Light_7 = ifelse(Noon_Day %in% seq(min_day_start, min_day_end), "Min", "Max"))
  }) %>% do.call("rbind", .) %>%
  set_rownames(1L:nrow(.)) %>%
  aggregate(cbind(Activity_IV_expanding_3, Activity_IS_expanding_3,
                  Activity_IV_moving_3, Activity_IS_moving_3,
                  Activity_IV_expanding_7, Activity_IS_expanding_7,
                  Activity_IV_moving_7, Activity_IS_moving_7, RA, DAR,
                  Sleep_Minutes, WASO_Minutes, SE) ~ patient_ID + Group_PTSD + MBLT_Group + Min_Max_Light_7,
            ., mean) %>%
  arrange(patient_ID, Min_Max_Light_7)

FL_results_light_days_14 <- split(results$major_results, results$major_results$patient_ID) %>%
  lapply(., function(ii) {
    pat_ID <- unique(ii$patient_ID)
    max_day_start <- gsub("-(.*)", "", light_days_14$Max_Days[light_days_14$patient_ID == pat_ID]) %>% as.numeric
    max_day_end <- gsub("(.*)-", "", light_days_14$Max_Days[light_days_14$patient_ID == pat_ID]) %>% as.numeric
    min_day_start <- gsub("-(.*)", "", light_days_14$Min_Days[light_days_14$patient_ID == pat_ID]) %>% as.numeric
    min_day_end <- gsub("(.*)-", "", light_days_14$Min_Days[light_days_14$patient_ID == pat_ID]) %>% as.numeric
    
    ii %<>%
      filter(Noon_Day %in% c(seq(max_day_start, max_day_end), seq(min_day_start, min_day_end))) %>%
      mutate(Min_Max_Light_14 = ifelse(Noon_Day %in% seq(min_day_start, min_day_end), "Min", "Max"))
  }) %>% do.call("rbind", .) %>%
  set_rownames(1L:nrow(.)) %>%
  aggregate(cbind(Activity_IV_expanding_3, Activity_IS_expanding_3,
                  Activity_IV_moving_3, Activity_IS_moving_3,
                  Activity_IV_expanding_7, Activity_IS_expanding_7,
                  Activity_IV_moving_7, Activity_IS_moving_7, RA, DAR,
                  Sleep_Minutes, WASO_Minutes, SE) ~ patient_ID + Group_PTSD + MBLT_Group + Min_Max_Light_14,
            ., mean) %>%
  arrange(patient_ID, Min_Max_Light_14)

write.csv(FL_results_light_days_7, "Min_Max_Light_7.csv")
write.csv(FL_results_light_days_14, "Min_Max_Light_14.csv")

## ---

results$FL_results_final <- FL_results_final

FL_dt <- data.table::as.data.table(FL_results_final)

data.table::dcast(FL_dt, formula = patient_ID + Week + PTSD + Light ~ .)

# anova(lm())

## ------ First/Middle/Last ------ 

# # Creation
# 
# 
# 
# ## Helper dataframe
# results$FML_days <- merge(aggregate(Noon_Day ~ patient_ID, results$major_results, function(x) min(x) + 2) %>%
#         rename(Min_Day = Noon_Day),
#         aggregate(Noon_Day ~ patient_ID, results$major_results, function(x) find_closest(x, 7)) %>% 
#           rename(Middle_Day = Noon_Day),
#       by = "patient_ID") %>%
#   merge(., aggregate(Noon_Day ~ patient_ID, results$major_results, function(x) find_closest(x, 28) - 1) %>% 
#           rename(Max_Day = Noon_Day),
#         by = "patient_ID") %>%
#   arrange(patient_ID)
# 
# ## FML DF
# results$FML <- split(results$major_results, results$major_results$patient_ID) %>%
#   lapply(., function(ii) {
#     ID <- unique(ii$patient_ID)
#     dayz <- unlist(results$FML_days[results$FML_days$patient_ID == ID, ][-1])
#     ii <- ii[ii$Noon_Day %in% dayz, ]
#     ii$Time_Point <- c("Initial", "Base_End", "Final")
#     return(ii)
# }) %>%
#   do.call("rbind", .) %>%
#   set_rownames(1:nrow(.)) %>%
#   filter(., Time_Point != "Base_End")
# 
# #  Analysis
# 
# FML_vars <- c("Sleep_Minutes", "WASO_Minutes", "SE", "Activity_IS",
#               "Activity_IV", "Light_IS", "Light_IV", "RA")
# 
# lapply(FML_vars, function(j) {
#   form <- as.formula(paste(j, "~ Group * Time_Point"))
#   ANOV <- anova(lm(form, results$FML))
#   
#   sf <- function(num) {
#     signif(num, 3)
#   }
#   
#   data.frame(IV = j,
#              Group_F = sf(ANOV$`F value`[1]),
#              Group_p = sf(ANOV$`Pr(>F)`[1]),
#              Time_F = sf(ANOV$`F value`[2]),
#              Time_p = sf(ANOV$`Pr(>F)`[2]),
#              Int_F = sf(ANOV$`F value`[3]),
#              Int_p = sf(ANOV$`Pr(>F)`[3]))
#   
# }) %>%
#   do.call("rbind", .)
# 
# # Group Plot
# 
# Group_FML_plot <- function(var) {
#   i_1 <- select(results$FML, one_of(var), Group, Time_Point)
#   
#   i_2 <- merge(aggregate(. ~ Group + Time_Point, i_1, mean),
#                aggregate(. ~ Group + Time_Point, i_1, sciplot::se) %>%
#                  set_colnames(c("Group", "Time_Point", "SEM")),
#                by = c("Group", "Time_Point")) %>%
#     mutate(ymax = .[, var] + SEM,
#            ymin = .[, var] - SEM,
#            Time_Point = factor(Time_Point, levels = c("Initial", "Base_End", "Final")))
#   
#   p <- ggplot(i_2, aes_string("Group", var, fill = "Time_Point")) +
#     geom_bar(position = "dodge", stat = "identity", colour='black') + 
#     geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, position=position_dodge(.9))
#   
#   plot(p)
# }
# 
# Group_FML_plotFML_plot(FML_vars[1])
# 
# # Individual Plot
# 
# Indiv_FML_plot <- function(var) {
#   i_1 <- select(results$major_results, patient_ID, one_of(var), Group, Noon_Day) %>%
#     filter(Group == "CCT")
#   
#   p <- ggplot(i_1, aes_string("Noon_Day", var, group = "patient_ID")) +
#     geom_point(aes(color = Group), size = 3) + 
#     geom_line(aes(color = Group), size = 1)
#     
#   
#   plot(p)
# }
# 
# Indiv_FML_plot(FML_vars[1])


## ------------ In Development ------------

## ------ Median Light Exposure -------

#' TODO: This aggregation should also take into account Baseline vs. Test
#'       timepoints once we have subjects with Negative Ion Generators and
#'       lightboxes with built-in baseline periods.

results$median_light_exposure <-
  cbind(aggregate(Light ~ patient_ID + Hour, actigraphy, median),
        aggregate(Light ~ patient_ID + Hour, actigraphy, opelr::sem)[3] %>%
          rename(SEM = Light)) %>%
  mutate(ymax = Light + SEM,
         ymin = Light - SEM) %>%
  arrange(patient_ID, Hour)

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
