## Script: Import and Shaping of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-12
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(zoo)
# devtools::install_github("opelr/opelR")

FILE_PATH <- ".\\Data\\Raw\\"
FILE_MASK <- "(.*)_New_Analysis.csv"

acti_files <- data.frame(File = list.files(path = FILE_PATH, pattern = FILE_MASK)) %>%
  mutate(rootpath = normalizePath(paste0(getwd(), "\\Data\\Raw\\", File)))

## ------- Major ETL Functions ------

get_actigraphy_headers <- function(path) {
  acti <- read.table(path, header = FALSE, sep = ",", row.names = NULL,
                     col.names = paste0("X", 1:44), fill = T)

  ## Strip header information to data frame
  patient_name <- as.character(acti$X2[acti$X1 == "Identity:"])
  patient_sex <- as.character(acti$X2[acti$X1 == "Gender:"])
  patient_DOB <- as.character(acti$X2[acti$X1 == "Date of Birth:"])
  patient_age <- as.character(acti$X2[acti$X1 == "Age (at start of data collection):"])
  
  recording_date_start <- acti$X2[acti$X1 == "Data Collection Start Date:"]
  recording_date_end <- acti$X2[acti$X1 == "Data Collection End Date:"]
  recording_days_span <- acti$X2[acti$X1 == "Number of Days:"]
  recording_epoch_length <- acti$X2[acti$X1 == "Epoch Length:"]
  
  watch_serial_no <- acti$X2[acti$X1 == "Actiwatch Serial Number:"]
  calibration_activity <- acti$X2[acti$X1 == "Activity Calibration Factor:"]
  calibration_light <- acti$X2[acti$X1 == "White Calibration Factor:"]
  threshold_wake <- acti$X2[acti$X1 == "Wake Threshold Value:"]
  threshold_light <- acti$X2[acti$X1 == "White Light Threshold:"]
  
  header_info <- data.frame(patient_ID = gsub(" ", "_", patient_name)) %>%
    mutate(patient_sex = patient_sex,
           patient_DOB = as.POSIXct(strptime(patient_DOB,
                                             format = "%Y-%m-%d")),
           patient_age = patient_age,
           recording_date_start = as.POSIXct(strptime(recording_date_start,
                                                      format = "%Y-%m-%d")),
           recording_date_end = as.POSIXct(strptime(recording_date_end,
                                                    format = "%Y-%m-%d")),
           recording_days_span = recording_days_span,
           recording_epoch_length = recording_epoch_length,
           watch_serial_no = watch_serial_no,
           calibration_activity = calibration_activity,
           calibration_light = calibration_light,
           threshold_wake = threshold_wake,
           threshold_light = threshold_light)
  
  return(header_info)
}

parse_actigraphy_data <- function(path) {
  ## Pull epoch-by-epoch data
  acti_1 <- read.table(path, header = FALSE, sep = ",", row.names = NULL,
                     col.names = paste0("X", 1:44), fill = T)
  start_row <- which(acti_1$X1 == "-------------------- Epoch-by-Epoch Data -------------------")
  
  acti_2 <- acti_1[start_row:nrow(acti_1), ]
  line_row <- which(acti_2$X1 == "Line")
  
  acti_3 <- acti_2[line_row:nrow(acti_2), ]
  header_cols <- acti_3[1, ]
  header_cols <- as.character(header_cols[header_cols != ""])
  
  acti <- acti_3[2:nrow(acti_3), 1:length(header_cols)] %>%
    set_colnames(as.character(header_cols)) %>%
    set_rownames(1:nrow(.)) %>%
    dplyr::select(., -`Interval Status`, -`S/W Status`) %>%
    rename(., Light = `White Light`, Sleep_Wake = `Sleep/Wake`) %>%
    mutate(Sleep_Wake = factor(Sleep_Wake, levels = 0:1,
                               labels = c("Sleep", "Wake")))
  
  ## Get Name
  nam <- as.character(acti_files$File[acti_files$rootpath == path])
  nam <- strsplit(nam, "_")[[1]][1:2]
  
  acti$patient_ID = paste0(nam, collapse = "_")
  
  ## Date/Time Calculations
  acti$DateTime = as.POSIXct(strptime(x = paste(acti$Date, acti$Time),
                                          format = "%Y-%m-%d %I:%M:%S %p"))
  acti$Time = as.POSIXct(strptime(acti$Time, format = "%I:%M:%S %p"))
  
  acti %<>%
    mutate(Hour = hour(DateTime),
           AM_PM = factor(floor(Hour/12), levels = 0:1,
                          labels = c("AM", "PM")),
           Day_of_Week = weekdays(DateTime),
           Weekend = factor(ifelse(Day_of_Week %in% c("Saturday", "Sunday"),
                            "Weekend", "Weekday")),
           DateAbbr = format(DateTime, format = "%b %d"))
    
  return(acti)
  
}

## ------ Sunset Function ------

sunset <- function(city, state = "OR", date) {
  datetime <- as.POSIXct(strptime(date, format = "%F")) 
  if (is.na(datetime)) {stop("date variable must be in 'yyyy-mm-dd' format")}
  
  link <- paste0("http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=", year(datetime),
                 "&task=0&state=", state, "&place=", city)
  
  sun <- suppressWarnings(readLines(link))
  
  cityVerify <- sun[26]
  if(grepl("Unable", cityVerify)) {stop("City is not in specified state. Please revise search")}
  
  times <- suppressWarnings(lapply(35:65, function(ii) {
    jj <- gsub("             ", "  NA  ", sun[ii]) %>%
      strsplit(., "  ")
    return(jj)
  }) %>%
    unlist(.) %>%
    matrix(., ncol = 13, byrow = T) %>%
    as.data.frame(.) %>%
    set_colnames(c("Day", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) %>%
    reshape2::melt(., "Day") %>%
    set_colnames(c("Day", "Month", "Value")) %>%
    mutate(Sunrise = gsub(" [0-9]{4}", "", Value),
           Sunset = gsub("[0-9]{4} ", "", Value),
           Date = as.POSIXct(strptime(paste0(Month, Day, year(datetime)),
                                      format = "%b%d%Y")),
           Sunrise = as.POSIXct(strptime(paste0(Month, Day, year(datetime),
                                                Sunrise), format = "%b%d%Y%H%M")),
           Sunset = as.POSIXct(strptime(paste0(Month, Day, year(datetime),
                                               Sunset), format = "%b%d%Y%H%M"))) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::select(-Value, -Day, -Month))
  
  out <- times[as.character(times$Date) == date, ]
  return(out)
  
}

## ------ Create Data Frames ------

acti_files <- do.call("rbind", lapply(acti_files$rootpath, get_actigraphy_headers)) %>%
  cbind(acti_files, .)
actigraphy <- do.call("rbind", lapply(acti_files$rootpath, parse_actigraphy_data)) %>%
  mutate(Date = as.character(Date))

sunsetDF <- do.call("rbind", lapply(unique(actigraphy$Date), function (x) {
    sunset(city = "Portland", state = "OR", date = as.character(x))
  })) %>%
  mutate(Date = as.character(Date))

actigraphy %<>% merge(., sunsetDF, by = "Date") %>%
  mutate(SunPeriod = factor(ifelse(DateTime <= Sunrise, "Dawn",
                            ifelse(DateTime > Sunrise & DateTime < Sunset, "Day",
                            ifelse(DateTime >= Sunset, "Night", NA))))) %>%
  dplyr::arrange(., patient_ID, DateTime) %>%
  mutate(Activity = opelr::numfac(Activity) + 1,
         Light = opelr::numfac(Light) + 1,
         DAR_Period = factor(ifelse(DateTime < paste(Date, "07:00:00 PDT") |
                                    DateTime > paste(Date, "21:59:59 PDT"),
                                    "Night", "Day")),
         log_Activity = log10(Activity),
         log_Light = log10(Light),
         Day = factor(Day, levels = unique(Day)))

## ------ Rolling Window & Off-Wrist Functions ------

moving_window <- function(dataframe, column = "Activity", window, step, FUN, align) {
  if (! align %in% c("left", "right", "center")) stop('"align" must be "left", "right", or "center"')
  if (align == "center" & window %% 2 == 0) stop("'window' must be odd when 'align = center'")
  
  df <- dataframe
  
  do.call("c", lapply(unique(df[, "patient_ID"]), function(ii) {
    data <- df[df[, "patient_ID"] == ii, column]
    
    total <- length(data)
    
    if (align == "left") {
      spots <- seq(from = 1, to = total - window + 1, by = step)
    } else if (align == "right") {
      spots <- seq(from = window, to = total, by = step)
    } else {
      spots <- seq(from = ceiling(window/2), to = total - floor(window/2),
                   by = step)
    }
    
    result <- vector(length = length(spots))
    
    for(i in 1:length(spots)){
      if (align == "left") {
        result[i] <- eval(call(FUN, data[spots[i]:(spots[i] + window - 1)]))
      } else if (align == "right") {
        result[i] <- eval(call(FUN, data[(spots[i] - window + 1):spots[i]]))
      } else {
        result[i] <- eval(call(FUN, data[(spots[i] - floor(window/2)):(spots[i] + floor(window/2))]))
      }
    }
    
    if (align == "left") {
      result <- c(result, rep(NA, window - 1))
    } else if (align == "right") {
      result <- c(rep(NA, window - 1), result)
    } else {
      result <- c(rep(NA, floor(window/2)), result, rep(NA, floor(window/2)))
    }
    
    return(result)
  }))
}

zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  ## Tests for equality
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

off_wrist_detector <- function(data, threshold) {
  df <- data
  
  rle_values <- rle(df$No_Activity_Change_Window)$lengths
  
  lst <- lapply(1:length(rle_values), function(ii) {
    
    current_value <- rle_values[ii]
    
    ## grab elapsed rows (ii - 2)
    ## use na.locf to carry down decision
    if (ii == 1) {
      
      out <- rep(df$No_Activity_Change_Window[1], current_value)
      
    } else if (ii == length(rle_values)) {
      
      out <- rep(NA, current_value)
      
    } else {
      
      previous_value <- rle_values[ii - 1]
      next_value <- rle_values[ii + 1]
      
      elapsed_rows <- sum(rle_values[seq(0, ii - 1)]) + 1
      
      change_bool <- df[elapsed_rows, "No_Activity_Change_Window"]
      
      if (change_bool == TRUE & previous_value >= threshold & next_value >= threshold) {
        out <- rep.int(FALSE, current_value)
      } else if (change_bool == FALSE & previous_value >= threshold & next_value >= threshold) {
        out <- rep.int(TRUE, current_value)
      } else {
        out <- rep.int(change_bool, current_value)
      }
    }
    return(out)
  })
  
  return(do.call("c", lst))
}

## ------ Smooth Sleep Function ------

data <- dplyr::filter(actigraphy, patient_ID == "Ryan_Opel")
column <- "Sleep_Wake"

smooth_sleep <- function(data, column) {
  ## Internal functions
  get_rle_df <- function(x) {
    data.frame(Values = rle(as.character(x))$values,
               Lengths = rle(as.character(x))$lengths)
  }
  
  rle_df_to_vector <- function(df) {rep(df$Values, df$Lengths)}
  
  cole_post_process <- function(df) {
    df_old <- df
    
    for (ii in 1:nrow(df)) {
      if (is.na(df$Values[ii])) {
        df$Values[ii] <- NA
      } else if (df$Values[ii] == "Sleep") {
        # 2) Sleep of 1, 3, or 4 minutes was rescored as wake if it preceded at least 4(2), 10(5), or 15(7) minutes of wake;
        # 3) Sleep of 6 or 10 minutes surrounded by at least 10 or 20 minutes of wake was rescored as wake.
        if ((df$Lengths[ii] == 1 & df$Lengths[ii - 1] >= 2) | 
            (df$Lengths[ii] == 2 & df$Lengths[ii - 1] >= 6) |
            (df$Lengths[ii] %in% c(3, 4) & df$Lengths[ii - 1] >= 5 & df$Lengths[ii + 1] >= 5) |
            (df$Lengths[ii] == 5 & df$Lengths[ii - 1] >= 10 & df$Lengths[ii + 1] >= 10)) {
          df$Values[ii] <- "Wake"
        }
      }
    }
    
    if (identical(df_old, df)) {
      return(df)
    } else {
      cole_post_process(get_rle_df(rle_df_to_vector(df)))
    }
  }
  
  
  vec <- data[, column]
  rle_df <- get_rle_df(vec)
  
  # 1) Wake of 1 epoch was rescored as sleep
  rle_df$Values[rle_df$Values == "Wake" & rle_df$Lengths == 1] <- "Sleep"
  rle_df <- get_rle_df(rle_df_to_vector(rle_df))
  

  
  out <- cole_post_process(rle_df) %>%
    rle_df_to_vector(.)
  
  return(out)
}

## ------ Lotjonen Parameters, Off-Wrist Detection, and Sleep Smoothing ------

actigraphy <- split(actigraphy, actigraphy$patient_ID) %>%
  lapply(., function(ii) {
    ii %<>% 
      mutate(
        Lotjonen_mean = moving_window(., "Activity", 15, 1,FUN = "mean", "center"),
        Lotjonen_sd = moving_window(., "Activity", 17, 1, FUN = "sd", "center"),
        Lotjonen_ln = log(Activity) + 0.1,
        Lotjonen_Counts = Activity > 10,
        Sleep = factor(ifelse(Lotjonen_mean < 100, "Sleep", "Wake")),
        Activity_diff = do.call("c", lapply(unique(.[, "patient_ID"]), function(ii) {
          c(NaN, diff(.[.[, "patient_ID"] == ii, "Activity"], 1))})),
        No_Activity_Change_Window = moving_window(., "Activity", 22, 1,
                                                  FUN = "zero_range", "left"),
        No_Activity_Change_Window = na.locf(No_Activity_Change_Window),
        No_Activity_Length = rep(rle(No_Activity_Change_Window)[["lengths"]],
                                 rle(No_Activity_Change_Window)[["lengths"]])) %>%
      mutate(Lotjonen_nat = moving_window(., "Lotjonen_Counts", 23, 1,
                                          FUN = "sum", "center"),
             Off_Wrist = off_wrist_detector(., 30),
             Lotjonen_Sleep = 1.687 + (0.002 * Activity) - (0.034 * Lotjonen_mean) - 
               (0.419 * Lotjonen_nat) + (0.007 * Lotjonen_sd) - (0.127 * Lotjonen_ln),
             Lotjonen_Sleep = factor(ifelse(Lotjonen_Sleep < 0, "Wake", "Sleep")),
             Noon_Date = factor(date(DateTime + hours(12))))
    
    ## Calculate Noon-Noon days
    which_day_date <- table(ii[, "Day"], ii[, "Date"]) %>%
      as.data.frame.table(.) %>%
      dplyr::filter(Freq > 0) %>%
      rename(Noon_Day = Var1, Noon_Date = Var2) %>%
      dplyr::select(-Freq)
    
    ii %<>% merge(., which_day_date, by = "Noon_Date") %>%
      select(-Noon_Date) %>%
      mutate(Sleep_Smooth = smooth_sleep(., "Sleep"),
             Sleep_Wake_Smooth = smooth_sleep(., "Sleep_Wake"),
             Lotjonen_Sleep_Smooth = smooth_sleep(., "Lotjonen_Sleep"))

    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.)) %>%
  dplyr::filter(!is.na(Activity))

## ------ Sleep Agreement ------

dat <- actigraphy[actigraphy$patient_ID == "Ryan_Opel", ] %>%
  mutate(Sleep_Agree = apply(.[, c("Sleep", "Sleep_Wake", "Lotjonen_Sleep")],
                           1, function(ii) {ifelse(any(is.na(ii)), NA, length(unique(unlist(ii))) == 1)}))

View(dat[, c("Activity", "Light", "DateTime", "Sleep", "Sleep_Wake", "Lotjonen_Sleep",
        "Sleep_Agree", "No_Activity_Change_Window", "No_Activity_Length", "Off_Wrist")])

## ------ Light Adherence ------

### Count period where mean Light window > 10,000/5,000 across 15/30 epochs
### How many of these bouts per day

actigraphy <- split(actigraphy, actigraphy$patient_ID) %>%
  lapply(., function(ii) {
    ii %<>% 
      mutate(
        Light_mean_15 = moving_window(., "Light", 15, 1, FUN = "median", "center"),
        Light_15_5k = Light_mean_15 >= 5000, 
        Light_15_10k = Light_mean_15 >= 10000, 
        Light_mean_30 = moving_window(., "Light", 31, 1, FUN = "median", "center"),
        Light_30_5k = Light_mean_30 >= 5000,
        Light_30_10k = Light_mean_30 >= 10000)
    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.))

## ------ Saving RDS ------

saveRDS(acti_files, ".\\Data\\actigraphy_header.rds")
saveRDS(actigraphy, ".\\Data\\actigraphy_data.rds")
