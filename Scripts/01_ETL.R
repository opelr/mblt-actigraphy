## Script: Import and Shaping of Actiwatch Data
## Project: mblt-actigraphy
## Author: Ryan Opel -- @opelr
## Date: 2017-11-15
## Version: 0.3.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(maptools)
library(zoo)
library(openxlsx)
library(reshape2)
# devtools::install_github("opelr/opelR")

source("Scripts/Functions.R")

FILE_PATH <- "D:/data/acquired_data/human/h4085/Actigraphy/CSV/"
FILE_MASK <- "(.*)_New_Analysis.csv"

## ------ File Repository ------
#' Build repository of file paths to be read in

acti_files <- data.frame(File = list.files(path = FILE_PATH, pattern = FILE_MASK)) %>%
  mutate(rootpath = normalizePath(paste0(FILE_PATH, File)))

## ------ Patient Catalog ------
#' 

catalog <- read.xlsx("Data/Raw/4085_Events_Master_Current.xlsx", 1) %>%
  select(ID, Group, Actiwatch_ID, Second_Actiwatch_ID, LightPad_Kit)

watch_ids <- melt(select(catalog, -LightPad_Kit), id = c("ID", "Group")) %>%
  mutate(Watch_Assign = factor(variable,
                               levels = c("Actiwatch_ID", "Second_Actiwatch_ID"),
                               labels = c("Primary", "Secondary"))) %>%
  select(-variable) %>%
  rename(watch_ID = value, patient_ID = ID)

light_box <- catalog[, c("ID", "LightPad_Kit")] %>%
  rename(patient_ID = ID) %>%
  mutate(MBLT_Group = ifelse(is.na(LightPad_Kit), F, T)) %>%
  select(-LightPad_Kit)

## ------ PTSD Status from PCL-5 ------

PCL <- read.csv("Data/Raw/4085_Baseline_PCL.csv") %>% 
  select(-redcap_event_name, -pcl5_notes, -pcl5_complete) %>%
  rename(patient_ID = record_id) %>%
  mutate(pcl_sum = rowSums(.[, grepl("pcl5_[0-9]", colnames(.))]))

pcl_b_cols <- colnames(PCL)[grepl('pcl5_[1-5]{1}', colnames(PCL))]
pcl_c_cols <- colnames(PCL)[grepl('pcl5_[6-7]{1}', colnames(PCL))]
pcl_d_cols <- colnames(PCL)[grepl('pcl5_([8-9]{1}|1[0-4]{1})', colnames(PCL))]
pcl_e_cols <- colnames(PCL)[grepl('pcl5_(1[5-9]|20)', colnames(PCL))]

PCL %<>%
  mutate(pcl_b_intrusion = survSums(pcl_b_cols),
         pcl_c_avoidance = survSums(pcl_c_cols),
         pcl_d_mood_cognition = survSums(pcl_d_cols),
         pcl_e_arousal_activity = survSums(pcl_e_cols),
         cluster_ptsd =
           rowSums(.[, pcl_b_cols] >= 2, na.rm = T) >= 1 &
           rowSums(.[, pcl_c_cols] >= 2, na.rm = T) >= 1 &
           rowSums(.[, pcl_d_cols] >= 2, na.rm = T) >= 2 &
           rowSums(.[, pcl_e_cols] >= 2, na.rm = T) >= 2,
         Group_PTSD = factor(ifelse(cluster_ptsd == T & pcl_sum >= 33,
                                    "PTSD", "No_PTSD")))

rm(list = regmatches(ls(), regexpr("pcl_[a-z]_cols", ls())))

## ------ Patient Chronotype from rMEQ ------

rMEQ <- read.csv("Data/Raw/4085_Baseline_rMEQ.csv") %>%
  set_colnames(c("patient_ID", "Event_Name", "Wake_Time", "Refreshed",
                 "Bed_Time", "Best_Time", "Chronotype", "Complete")) %>%
  mutate(DAR_Wake_Time = mapply(Wake_Time, FUN = strip_waketime),
         DAR_Bed_Time = mapply(DAR_Wake_Time, FUN = function(i) waketime_plus(i, 12)),
         DAR_Light_Box_Time = mapply(DAR_Wake_Time, FUN = function(i) waketime_plus(i, 3))) %>%
  select(patient_ID, DAR_Wake_Time, DAR_Bed_Time, DAR_Light_Box_Time)

## ------ Create Data Frames ------

# Append header information onto file list data.frame
acti_files <- do.call("rbind", lapply(acti_files$rootpath, get_actigraphy_headers)) %>%
  cbind(acti_files, .)

# Create main data.frame
actigraphy <- do.call("rbind", lapply(acti_files$rootpath, parse_actigraphy_data)) %>%
  mutate(Date = as.character(Date))

# Get sunrise and sunset times
sunsetDF <- do.call("rbind", lapply(unique(actigraphy$Date), function (x) {
  get_sun_times(lat = 45.5231, long = -122.6765, date = as.character(x))
  })) %>%
  mutate(Date = as.character(Date))

# Merge all helper DFs wtih major DF
actigraphy %<>% merge(., sunsetDF, by = "Date") %>%
  mutate(SunPeriod = factor(ifelse(DateTime > Sunrise & DateTime < Sunset, "Day",
                            "Night"))) %>%
  dplyr::arrange(watch_ID, DateTime) %>%
  mutate(Activity = opelr::numfac(Activity) + 1,
         Light = opelr::numfac(Light) + 1,
         DAR_77 = factor(ifelse(DateTime <= paste(Date, "07:00:00 PDT") |
                                  DateTime >= paste(Date, "18:59:59 PDT"),
                                "Night", "Day")),
         log_Activity = log10(Activity),
         log_Light = log10(Light),
         Day = factor(Day, levels = unique(Day))) %>%
  merge(., watch_ids, by = "watch_ID") %>%
  merge(., light_box, by = "patient_ID") %>%
  merge(., PCL[, c("patient_ID", "Group_PTSD")], by = "patient_ID") %>%
  merge(., rMEQ, by = "patient_ID") %>%
  arrange(patient_ID, DateTime) %>%
  filter(!patient_ID %in% c("15")) %>%
  mutate(DAR = DAR_ifelse(., "DAR_Wake_Time", "DAR_Bed_Time", "Day", "Night"),
         Light_Box_Period = DAR_ifelse(., "DAR_Wake_Time", "DAR_Light_Box_Time", "Light", "Non_Light"),
         ZT_6 = cut(Time, seq(strptime("00:00", "%R"), by = "6 hours", length.out = 5),
                    labels = c("12AM-6AM", "6AM-12PM", "12PM-6PM", "6PM-12AM")),
         ZT_12 = cut(Time, seq(strptime("06:00", "%R"), by = "12 hours", length.out = 3),
                    labels = c("6AM-6PM", "6PM-6AM")),
         ZT_12 = factor(ifelse(is.na(ZT_12), "6PM-6AM", as.character(ZT_12)),
                        levels = c("6AM-6PM", "6PM-6AM")))

rm(list = c("PCL", "sunsetDF", "light_box"))

## ------ Filter out Overlapping Dates ------

#' FIXME: Will be good to put this data.frame into a CSV when we have more CCT
#' patients, but this works fine until January.

overlap_dates <- data.frame(patient_ID = c(10, 11, 13, 14),
                            DateTimePrimary = as.POSIXct("2017-10-13 10:58:00",
                                                         tz = "America/Los_Angeles")) %>%
  mutate(DateTimeSecondary = DateTimePrimary + minutes(2))

actigraphy <- split(actigraphy, actigraphy$patient_ID) %>%
  lapply(., function (j) {
    
    pat_ID <- unique(j$patient_ID)
    
    if (pat_ID %in% overlap_dates$patient_ID) {
      j <- split(j, j$Watch_Assign) %>%
        lapply(., function (ii) {
          if (unique(ii$Watch_Assign) == "Primary") {
            ii <- filter(ii, DateTime <= overlap_dates$DateTimePrimary[overlap_dates$patient_ID == pat_ID])
          } else if (unique(ii$Watch_Assign) == "Secondary") {
            ii <- filter(ii, DateTime >= overlap_dates$DateTimeSecondary[overlap_dates$patient_ID == pat_ID])
          }
        }) %>% 
        do.call("rbind", .)
    }
    return(j)
  }) %>%
  do.call("rbind", .)

## ------ Determining Day Number by Date ------

#' Day number breaks when we give patients multiple watches, as each watch
#' starts at 1. Instead, we'll enumerate the ordered dates for each patient to
#' determine day number.

updated_day_date <- xtabs(~ patient_ID + Date, actigraphy) %>%
  as.data.frame %>%
  filter(Freq > 0) %>%
  arrange(patient_ID, Date) %>%
  mutate(Day = lapply(rle(as.numeric(as.character(.[, "patient_ID"])))$lengths,
                      function(i) seq(0, i - 1, 1)) %>%
           do.call('c', .)) %>%
  select(-Freq)

actigraphy %<>%
  select(-Day) %>%
  merge(., updated_day_date, by = c("patient_ID", "Date"))

updated_day_date %<>%
  mutate(DT = as.POSIXct(strptime(paste(Date, "12:00:00"), "%F %T")))

rm(updated_day_date)

## ------ Lotjonen Parameters, Off-Wrist Detection, and Sleep Smoothing ------

actigraphy <- split(actigraphy, actigraphy$patient_ID) %>%
  lapply(., function(ii) {
    ii %<>% 
      mutate(
        Sleep_Actiware = actiware_sleep(., "Activity", 40),
        Lotjonen_mean = moving_window(., "Activity", 15, 1, FUN = "mean", "center"),
        Lotjonen_sd = moving_window(., "Activity", 17, 1, FUN = "sd", "center"),
        Lotjonen_ln = log(Activity) + 0.1,
        Lotjonen_Counts = Activity > 10,
        Sleep_Thresh = factor(ifelse(Lotjonen_mean < 100, "Sleep", "Wake")),
        Activity_diff = do.call("c", lapply(unique(.[, "watch_ID"]), function(ii) {
          c(NaN, diff(.[.[, "watch_ID"] == ii, "Activity"], 1))})),
        No_Activity_Change_Window = moving_window(., "Activity", 90, 1,
                                                  FUN = "zero_range", "left"),
        No_Activity_Change_Window = na.locf(No_Activity_Change_Window),
        No_Activity_Length = rep(rle(No_Activity_Change_Window)[["lengths"]],
                                 rle(No_Activity_Change_Window)[["lengths"]]),
        Thresh_Activity_Change_Window = moving_window(., "Activity", 90, 1,
                                                      FUN = "threshold_range", "left"),
        Thresh_Activity_Change_Window = na.locf(Thresh_Activity_Change_Window),
        Thresh_Activity_Length = rep(rle(Thresh_Activity_Change_Window)[["lengths"]],
                                 rle(Thresh_Activity_Change_Window)[["lengths"]])) %>%
      mutate(Lotjonen_nat = moving_window(., "Lotjonen_Counts", 23, 1,
                                          FUN = "sum", "center"),
             Off_Wrist = off_wrist_detector(., 1),
             Lotjonen_Sleep = 1.687 + (0.002 * Activity) - (0.034 * Lotjonen_mean) - 
               (0.419 * Lotjonen_nat) + (0.007 * Lotjonen_sd) - (0.127 * Lotjonen_ln),
             Lotjonen_Sleep = factor(ifelse(Lotjonen_Sleep < 0.0, "Wake", "Sleep")),
             Lotjonen_Sleep_2 = 2.457 - (0.004 * Activity) - (0.689 * Lotjonen_nat) -
               (0.007 * Lotjonen_sd) - (0.108 * Lotjonen_ln),
             Lotjonen_Sleep_2 = factor(ifelse(Lotjonen_Sleep_2 < 0.0, "Wake", "Sleep")),
             Noon_Date = as.character(date(DateTime + hours(12))))
    
    # Sleep Smoothing
    ii %<>%
      mutate(
        Sleep_Thresh_Smooth = smooth_sleep(., "Sleep_Thresh"),
        Sleep_Acti_Smooth = smooth_sleep(., "Sleep_Actiware"),
        Lotjonen_Sleep_Smooth = smooth_sleep(., "Lotjonen_Sleep"),
        Lotjonen_Sleep_2_Smooth = smooth_sleep(., "Lotjonen_Sleep_2")
      )
    
    ## Enumerate consolidated sleep/wake bouts
    ii %<>%
      mutate(Sleep_Bouts = get_sleep_bouts(Sleep_Acti_Smooth, 8))
    
    # Calculate Noon-Noon days
    which_day_date <- table(ii[, "Day"], ii[, "Date"]) %>%
      as.data.frame.table(.) %>%
      dplyr::filter(Freq > 0) %>%
      rename(Noon_Day = Var1, Noon_Date = Var2) %>%
      mutate(Noon_Date = as.character(Noon_Date)) %>%
      dplyr::select(-Freq)
    
    ii %<>% merge(., which_day_date, by = "Noon_Date", all.x = T) %>%
      select(-Noon_Date) %>%
      mutate(Noon_Day = opelr::numfac(Noon_Day),
             Noon_Day = Noon_Day - (min(Noon_Day, na.rm = T))) %>%
      arrange(DateTime)

    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.))

## ------ Light Adherence ------

# Count periods where median Light window > 10,000/5,000 lux/m2 across 15/30 epochs
actigraphy <- split(actigraphy, actigraphy$watch_ID) %>%
  lapply(., function(ii) {
    ii %<>% 
      mutate(
        Light_median_15 = moving_window(., column = "Light", 15, 1, FUN = "median", "left"),
        Light_15_5k = Light_median_15 >= 5000,
        Light_15_10k = Light_median_15 >= 10000,
        Light_median_30 = moving_window(., column = "Light", 31, 1, FUN = "median", "left"),
        Light_30_5k = Light_median_30 >= 5000,
        Light_30_10k = Light_median_30 >= 10000)
    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.))

# Smooth Light adherence
actigraphy <- split(actigraphy, actigraphy$watch_ID) %>%
  lapply(., function(ii) {
    ii %<>% 
      mutate(
        Light_15_5k_smooth = clean_adherence(., "Light_15_5k"),
        Light_15_10k_smooth = clean_adherence(., "Light_15_10k"),
        Light_30_5k_smooth = clean_adherence(., "Light_30_5k"),
        Light_30_10k_smooth = clean_adherence(., "Light_30_10k"))
    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.))

## ------ Activity Windowing - EE ------

actigraphy <- split(actigraphy, actigraphy$watch_ID) %>%
  lapply(., function(ii) {
    ii %<>% 
      mutate(
        Activity_median_15 = moving_window(., column = "Activity", 15, 1, FUN = "median", "left"),
        Activity_15_500 = Activity_median_15 >= 500, 
        Activity_15_1k = Activity_median_15 >= 1000, 
        Activity_15_2k = Activity_median_15 >= 2000, 
        Activity_median_30 = moving_window(., column = "Activity", 30, 1, FUN = "median", "left"),
        Activity_30_500 = Activity_median_30 >= 500, 
        Activity_30_1k = Activity_median_30 >= 1000,
        Activity_30_2k = Activity_median_30 >= 2000)
    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.))

## ------ Saving RDS ------

saveRDS(acti_files, ".\\Rmd\\Data\\actigraphy_header.rds")
saveRDS(actigraphy, ".\\Rmd\\Data\\actigraphy_static.rds")

## ------ Off Wrist Filtering ------

#' TODO: How should the overlap of data be handled? Currently we're throwing
#' out two days per patient in the middle of recording, which seems wasteful.
#' I'm thinking it would be good to record the exact time a new watch is
#' placed on the wrist of a CCT patient, and we can cut the tail/head of watch
#' #1/#2 at that point, respectively. 

### Exclude by any 3+ hour window
off_wrist <- xtabs(~ patient_ID + Thresh_Activity_Change_Window + Noon_Day, data = actigraphy) %>%
  as.data.frame.table %>%
  reshape2::dcast(., formula = patient_ID + Noon_Day ~ Thresh_Activity_Change_Window,
                  value.var = "Freq") %>% 
  filter((`FALSE` + `TRUE`) > 0) %>%
  rename(Watch_On = `FALSE`, Watch_Off = `TRUE`) %>%
  mutate(Percent_Off = 100 * Watch_Off/(Watch_On + Watch_Off),
         Total_Epochs = Watch_On + Watch_Off) %>%
  dplyr::filter(complete.cases(.)) %>%
  mutate(Consec_Days = split(., .[, "patient_ID"]) %>%
           lapply(., function(ii) {
             rep(rle(ii$Watch_On == 720)$lengths, rle(ii$Watch_On == 720)$lengths)
           }) %>%
           do.call("c", .),
         At_Least_96_Hours_On = Consec_Days >= 3 & Watch_On == 720)

### Merge
actigraphy <- merge(actigraphy,
                    off_wrist[,c("patient_ID", "Noon_Day", "At_Least_96_Hours_On")],
                    by = c("patient_ID", "Noon_Day")) %>%
  dplyr::filter(At_Least_96_Hours_On == T) %>%
  arrange(patient_ID, DateTime)

## ------ Hand-Added Medical History -----

CCT_medHist <- data.frame(patient_ID = c(10, 11, 12, 13, 14, 15, 16),
                          Sex = c("F", "M", "F", "F", "M", "M", "M"),
                          PTSD_CC = c(T, T, F, T, F, T, T),
                          TBI_CC = c(F, F, F, F, F, T, F),
                          Depression_CC = c(F, F, F, F, T, F, F),
                          Anxiety_CC = c(T, F, F, F, F, F, F),
                          OSA_CC = c(F, F, F, F, F, T, T))

actigraphy %<>% merge(., CCT_medHist, "patient_ID", all.x = T)

## ------ Save RDS ------

saveRDS(actigraphy, ".\\Rmd\\Data\\actigraphy_filtered.rds")
