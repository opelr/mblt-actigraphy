## Script: Import and Shaping of Actiwatch Data
## Project: mblt-actigraphy
## Author: Ryan Opel -- @opelr
## Date: 2018-01-05
## Version: 1.0.0

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
#' Build repository of file paths to be read later by two functions.

acti_files <- data.frame(File = list.files(path = FILE_PATH,
                                           pattern = FILE_MASK)) %>%
  mutate(rootpath = normalizePath(paste0(FILE_PATH, File)))

## ------ Patient Catalog ------
#' Import (or create) a table that links the unique ID for each
#' participant (`patient_ID`), the watch serial # (`Actiwatch_ID`), and
#' a way of discerning if your participant got light therapy.
#' 
#' If patients have multiple watches, separate serial numbers into multiple
#' columns, and melt as done below. This will allow for merging of all data
#' for a single participant, which will be done for the user in a later
#' section.
#' 
#' Output should look like the following (example data):
#' 
#'     >head(catalog, 5)
#'     ID   Group Actiwatch_ID Second_Actiwatch_ID LightPad_Kit
#'      1     MBL     A12694.1                <NA>           19
#'      2     MBL     A12641.1                <NA>            6
#'      3     MBL     A12616.1            A22618.1           17
#'      4 Control     A12622.1                <NA>         <NA>
#'      5     MBL     A12624.1                <NA>           13
#'     
#'     >head(watch_ids,5 )
#'     patient_ID   Group watch_ID Watch_Assign
#'              1     MBL A12694.1      Primary
#'              2     MBL A12641.1      Primary
#'              3     MBL A12616.1      Primary
#'              3     MBL A22618.1    Secondary
#'              4 Control A12622.1      Primary
#'     
#'     >head(light_box, 5)
#'     patient_ID MBLT_Group
#'              1       TRUE
#'              2       TRUE
#'              3       TRUE
#'              4      FALSE
#'              5       TRUE
#' 
#' For Lim Lab: Pull the most recent Events_Master Excel document from the
#'              H4085 folder in the Box account.

catalog <- read.xlsx("Data/Raw/4085_Events_Master_Current.xlsx", 1) %>%
  select(ID, Group, Actiwatch_ID, Second_Actiwatch_ID, LightPad_Kit)

watch_ids <- melt(select(catalog, -LightPad_Kit), id = c("ID", "Group")) %>%
  mutate(Watch_Assign = factor(variable,
                               levels = c("Actiwatch_ID", "Second_Actiwatch_ID"),
                               labels = c("Primary", "Secondary"))) %>%
  select(-variable) %>%
  rename(watch_ID = value, patient_ID = ID)

light_box <- catalog[, c("ID", "Group")] %>%
  rename(patient_ID = ID) %>%
  mutate(MBLT_Group = ifelse(Group == "Control", F, T)) %>%
  select(-Group)

## ------ Create Data Frames ------

#' Append header information onto file list dataframe.
acti_files <- lapply(acti_files$rootpath, get_actigraphy_headers) %>%
  do.call("rbind", .) %>%
  cbind(acti_files, .)

#' Create the main dataframe using the `parse_actigraphy_data` function.
actigraphy <- lapply(acti_files$rootpath, parse_actigraphy_data) %>%
  do.call("rbind", .) %>%
  mutate(Date = as.character(Date))

#' Pulling sunrise and sunset times for Portland (note the long and lat).
#' It would be good to incorporate individual participants' addresses/cities,
#' though I'm worried that this could interfere with PHI.
sunsetDF <- do.call("rbind", lapply(unique(actigraphy$Date), function (x) {
  get_sun_times(lat = 45.5231, long = -122.6765, date = as.character(x))
  })) %>%
  mutate(Date = as.character(Date))

#' Merge all the helper dataframes wtih our main dataframe.
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
  mutate(ZT_6 = cut(Time, seq(strptime("00:00", "%R"), by = "6 hours", length.out = 5),
                    labels = c("12AM-6AM", "6AM-12PM", "12PM-6PM", "6PM-12AM")),
         ZT_12 = cut(Time, seq(strptime("06:00", "%R"), by = "12 hours", length.out = 3),
                    labels = c("6AM-6PM", "6PM-6AM")),
         ZT_12 = factor(ifelse(is.na(ZT_12), "6PM-6AM", as.character(ZT_12)),
                        levels = c("6AM-6PM", "6PM-6AM"))) %>%
  arrange(patient_ID, DateTime)

rm(list = c("sunsetDF", "light_box"))

## ------ Filter out Overlapping Dates ------

#' For patients with multiple watches, we frequently run into the problem of
#' overlapping temporal data, especially when watches are replaced during a
#' long, continuous recording period.
#' 
#' This section solves this problem by filtering data that comes before/after
#' the watch switch. Here, the switch is hard-coded to occur at 11:00 AM,
#' but that wouldn't be hard to change on a per-patient basis.
#' 
#' The date origin at '1899-12-30' is due to importing from Excel -- worth
#' noting depending where your data is coming from.
#' 
#' Again, you'll need a dataframe that indicates the date of change, and
#' the watch serial number.

overlap_dates <- read.xlsx("Data/Raw/4085_Events_Master_Current.xlsx", 1) %>%
  select(ID, Actiwatch.Return.Date, Second_Actiwatch_ID) %>%
  set_colnames(gsub("[.]", "_", colnames(.))) %>%
  rename(patient_ID = ID) %>%
  mutate(Actiwatch_Return_Date = as.Date(Actiwatch_Return_Date,
                                         origin = "1899-12-30"),
         DateTimePrimary = as.POSIXct(paste(Actiwatch_Return_Date, "10:58:00"),
                                      tz = "America/Los_Angeles",
                                      format = "%F %T"),
         DateTimeSecondary = DateTimePrimary + minutes(2)) %>%
  filter(!is.na(Second_Actiwatch_ID)) %>%
  select(-Actiwatch_Return_Date, -Second_Actiwatch_ID)

actigraphy %<>% split(., actigraphy$patient_ID) %>%
  lapply(., function (j) {
    
      pat_ID <- unique(j$patient_ID)
    
    if (pat_ID %in% overlap_dates$patient_ID) {
      j <- split(j, j$Watch_Assign) 
      
      if (all(sapply(j, nrow) > 0)) {
        j %<>%
          lapply(., function (ii) {
            if (unique(ii$Watch_Assign) == "Primary") {
              ii <- filter(ii, DateTime <= overlap_dates$DateTimePrimary[overlap_dates$patient_ID == pat_ID])
            } else if (unique(ii$Watch_Assign) == "Secondary") {
              ii <- filter(ii, DateTime >= overlap_dates$DateTimeSecondary[overlap_dates$patient_ID == pat_ID])
            }
          })
      }
      j %<>% do.call("rbind", .)
    }
    return(j)
  }) %>%
  do.call("rbind", .) %>%
  arrange(patient_ID, DateTime) %>%
  set_rownames(1:nrow(.))

rm(overlap_dates)

## ------ Determining Day Number by Date ------

#' Similar to the prvious section, day number breaks when we give patients
#' multiple watches, as each watch starts at Day 1. Instead, we'll enumerate
#' the ordered dates from each patient to determine day number.

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

rm(updated_day_date)

## ------ Lotjonen Parameters, Off-Wrist Detection, and Sleep Smoothing ------

#' This and the following few sections add a ton of variables to our main
#' dataframe. Most relate to various ways to interpret sleep staging, light
#' exposure, and activity, but there are a few notable exceptions, like
#' `Noon_Day`.

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
      mutate(Noon_Day = opelr::numfac(Noon_Day)) %>%
      #        Noon_Day = Noon_Day - (min(Noon_Day, na.rm = T))) %>%
      arrange(DateTime)

    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.))

## ------ Light Adherence ------

#' Count periods where median Light window > 5,000/10,000 lux/m2 across
#'  15/30 epochs.
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

#' Smooth Light adherence columns generated in previous code block.
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

## ------ Activity Windowing ------

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

#' Save dataframes before we filter out 'off-wrist' days.

saveRDS(acti_files, ".\\Rmd\\Data\\actigraphy_header.rds")
saveRDS(actigraphy, ".\\Rmd\\Data\\actigraphy_static.rds")

## ------ Off Wrist Filtering ------

#' Exclude days that contain any 3+ hour window of inactivity
#' (`Thresh_Activity_Change_Window`).
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

actigraphy  %<>% merge(.,
                    off_wrist[,c("patient_ID", "Noon_Day", "At_Least_96_Hours_On")],
                    by = c("patient_ID", "Noon_Day")) %>%
  dplyr::filter(At_Least_96_Hours_On == T) %>%
  arrange(patient_ID, DateTime)

## ------ Save RDS ------

saveRDS(actigraphy, ".\\Rmd\\Data\\actigraphy_filtered.rds")
