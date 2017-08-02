## Script: Import and Shaping of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-21
## Version: 0.2.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(maptools)
library(zoo)
# devtools::install_github("opelr/opelR")

FILE_PATH <- ".\\Data\\Raw\\"
FILE_MASK <- "(.*)_New_Analysis.csv"

acti_files <- data.frame(File = list.files(path = FILE_PATH, pattern = FILE_MASK)) %>%
  mutate(rootpath = normalizePath(paste0(getwd(), "\\Data\\Raw\\", File)))

## ------- Major ETL Functions ------

get_actigraphy_headers <- function(path) {
  # Strips header information from exported Actiware CSV files
  # 
  # Args:
  #   path (str): Full filepath of an Actiware CSV
  # 
  # Returns:
  #   Data frame containing header information from a single file
  
  # Read CSV to table
  acti <- read.table(path, header = FALSE, sep = ",", row.names = NULL,
                     col.names = paste0("X", 1:44), fill = T)
  
  # Strip header information to temporary variables
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
  
  # Create data frame with stripped variables
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
  # Strips epoch-by-epoch data from exported Actiware CSV files
  # 
  # Args:
  #   path (str): Full filepath of an Actiware CSV
  # 
  # Returns:
  #   Data frame containing epoch-by-epoch data from a single file
  
  # Read entire CSV file to a table, force number of columns
  acti_1 <- read.table(path, header = FALSE, sep = ",", row.names = NULL,
                     col.names = paste0("X", 1:44), fill = T)
  start_row <- which(acti_1$X1 == "-------------------- Epoch-by-Epoch Data -------------------")
  
  # Skip down to Epoch-by-Epoch Data
  acti_2 <- acti_1[start_row:nrow(acti_1), ]
  line_row <- which(acti_2$X1 == "Line")
  
  acti_3 <- acti_2[line_row:nrow(acti_2), ]
  header_cols <- acti_3[1, ]
  header_cols <- as.character(header_cols[header_cols != ""]) %>%
    .[!is.na(.)]
  
  acti <- acti_3[2:nrow(acti_3), 1:length(header_cols)] %>%
    set_colnames(as.character(header_cols)) %>%
    set_rownames(1:nrow(.)) %>%
    dplyr::select(., -`Interval Status`, -`S/W Status`) %>%
    rename(., Light = `White Light`, Sleep_Wake = `Sleep/Wake`) %>%
    mutate(Sleep_Wake = factor(Sleep_Wake, levels = 0:1,
                               labels = c("Sleep", "Wake"))) %>%
    rename(Sleep_Acti = Sleep_Wake)
  
  # Get patient name
  nam <- as.character(acti_files$File[acti_files$rootpath == path])
  nam <- strsplit(nam, "_")[[1]][1:2]
  
  acti$patient_ID = paste0(nam, collapse = "_")
  
  # Date/Time Calculations
  acti$DateTime = as.POSIXct(strptime(x = paste(acti$Date, acti$Time),
                                          format = "%Y-%m-%d %I:%M:%S %p"))
  acti$Time = as.POSIXct(strptime(acti$Time, format = "%I:%M:%S %p"))
  
  acti %<>%
    mutate(ClockTime = as.character(strftime(Time, format = "%T")),
           Hour = hour(DateTime),
           AM_PM = factor(floor(Hour/12), levels = 0:1,
                          labels = c("AM", "PM")),
           Day_of_Week = weekdays(DateTime),
           Weekend = factor(ifelse(Day_of_Week %in% c("Saturday", "Sunday"),
                            "Weekend", "Weekday")),
           DateAbbr = format(DateTime, format = "%b %d"),
           Month = format(DateTime, format = "%B"))
    
  return(acti)
  
}

## ------ Sunset Function ------

get_sun_times <- function(lat, long, date, tz = "America/Los_Angeles") {
  # Gets sunrise and sunset time for a specific date and city (lat, long)
  # 
  # Args:
  #   lat (numeric): Latitude in degrees
  #   long (numeric): Longitude in degrees
  #   date (POSIXct): Date in yyyy-mm-dd format
  # 
  # Returns:
  #   Data frame containing epoch-by-epoch data from a single file
  
  coordin <- matrix(c(long, lat), nrow = 1)
  day <- as.POSIXct(date, tz=tz)
  
  sunrise <- sunriset(coordin, day, direction="sunrise", POSIXct.out=TRUE)
  sunset <- sunriset(coordin, day, direction="sunset", POSIXct.out=TRUE)
  
  data.frame(Date = as.Date(sunrise$time),
             Sunrise = sunrise$time,
             Sunset = sunset$time,
             Day_Length = as.numeric(sunset$time - sunrise$time))
}

## ------ Create Data Frames ------

# Append header information onto file list data.frame
acti_files <- do.call("rbind", lapply(acti_files$rootpath, get_actigraphy_headers)) %>%
  cbind(acti_files, .)

# Create main data.frame
actigraphy <- do.call("rbind", lapply(acti_files$rootpath, parse_actigraphy_data)) %>%
  mutate(Date = as.character(Date))

sunsetDF <- do.call("rbind", lapply(unique(actigraphy$Date), function (x) {
  get_sun_times(lat = 45.5231, long = -122.6765, date = as.character(x))
  })) %>%
  mutate(Date = as.character(Date))

actigraphy %<>% merge(., sunsetDF, by = "Date") %>%
  mutate(SunPeriod = factor(ifelse(DateTime > Sunrise & DateTime < Sunset, "Day",
                            "Night"))) %>%
  dplyr::arrange(., patient_ID, DateTime) %>%
  mutate(Activity = opelr::numfac(Activity) + 1,
         Light = opelr::numfac(Light) + 1,
         DAR_Period = factor(ifelse(DateTime < paste(Date, "06:30:00 PDT") |
                                    DateTime > paste(Date, "22:59:59 PDT"),
                                    "Night", "Day")),
         log_Activity = log10(Activity),
         log_Light = log10(Light),
         Day = factor(Day, levels = unique(Day)))

rm(sunsetDF)

## ------ Rolling Window & Off-Wrist Functions ------

moving_window <- function(dataframe, column = "Activity", window, step, FUN, align) {
  # Moving window wrapper for other functions
  # 
  # Args:
  #   dataframe (df): Data frame to iterate upon
  #   column (str): Which column will we be acting on
  #   window (int): Size of window (number of rows)
  #   step (int): How many rows to move each iteration
  #   FUN (str): Function to utilize
  #   align (str): Look forward ('left'), backward ('right'), or around ('center')
  # 
  # Returns:
  #   Vector of the same length
  
  if (! align %in% c("left", "right", "center"))
    stop('"align" must be "left", "right", or "center"')
  if (align == "center" & window %% 2 == 0)
    stop("'window' must be odd when 'align = center'")
  
  # Make copy of data
  df <- dataframe
  
  do.call("c", lapply(unique(df[, "patient_ID"]), function(ii) {
    data <- df[df[, "patient_ID"] == ii, column]
    
    total <- length(data)
    
    # Sequence window
    if (align == "left") {
      spots <- seq(from = 1, to = total - window + 1, by = step)
    } else if (align == "right") {
      spots <- seq(from = window, to = total, by = step)
    } else {
      spots <- seq(from = ceiling(window/2), to = total - floor(window/2),
                   by = step)
    }
    
    result <- vector(length = length(spots))
    
    # Call function upon the vector window
    for(i in 1:length(spots)){
      if (align == "left") {
        result[i] <- eval(call(FUN, data[spots[i]:(spots[i] + window - 1)]))
      } else if (align == "right") {
        result[i] <- eval(call(FUN, data[(spots[i] - window + 1):spots[i]]))
      } else {
        result[i] <- eval(call(FUN, data[(spots[i] - floor(window/2)):(spots[i] + floor(window/2))]))
      }
    }
    
    # Pad NA's to make vector the same length
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
  # Test for equality in a vector
  # 
  # Args:
  #   x (int): Vector of numbers
  # 
  # Returns:
  #   Boolean value
  
  if (length(x) == 1)
    return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

off_wrist_detector <- function(dataframe, threshold) {
  # Determines wearer status based on the number of non-changing rows in a look-ahead
  # 
  # Args:
  #   dataframe (df): Data frame to iterate upon
  #   threshold (int): Number of non-changing rows before off-wrist is called
  # 
  # Returns:
  #   Boolean value denoting off-wrist status
  
  # Copy data
  df <- dataframe
  
  rle_values <- rle(df$No_Activity_Change_Window)$lengths
  
  # Loop through each bout of non-changing activity
  lst <- lapply(1:length(rle_values), function(ii) {
    
    current_value <- rle_values[ii]
    
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

smooth_sleep <- function(data, column) {
  # Recursively smooth sleep staging
  # 
  # Args:
  #   data (df): Data frame to iterate upon
  #   column (str): Name of column to call smoothing functions on
  # 
  # Returns:
  #   Vector containing 'Sleep' and 'Wake', same length as 'column'
  
  # Convert vector to RLE object to data.frame
  get_rle_df <- function(x) {
    data.frame(Values = rle(as.character(x))$values,
               Lengths = rle(as.character(x))$lengths)
  }
  
  # Convert RLE data.frame to a vector
  rle_df_to_vector <- function(df) {rep(df$Values, df$Lengths)}
  
  # Recursively apply sleep smoothing to RLE DF object
  cole_post_process <- function(df) {
    
    # Copy data for comparison
    df_old <- df
    
    for (ii in 1:nrow(df)) {
      if (is.na(df$Values[ii]) | is.null(df$Values[ii])) {
        df$Values[ii] <- NA
      } else if (df$Values[ii] == "Sleep") {
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
  
  # 1) Wake of 1 or 2 epochs were rescored as sleep
  rle_df$Values[rle_df$Values == "Wake" & rle_df$Lengths %in%  c(1, 2)] <- "Sleep"
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
        Lotjonen_mean = moving_window(., "Activity", 15, 1, FUN = "mean", "center"),
        Lotjonen_sd = moving_window(., "Activity", 17, 1, FUN = "sd", "center"),
        Lotjonen_ln = log(Activity) + 0.1,
        Lotjonen_Counts = Activity > 10,
        Sleep_Thresh = factor(ifelse(Lotjonen_mean < 100, "Sleep", "Wake")),
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
             Lotjonen_Sleep = factor(ifelse(Lotjonen_Sleep < 0.0, "Wake", "Sleep")),
             Lotjonen_Sleep_2 = 2.457 - (0.004 * Activity) - (0.689 * Lotjonen_nat) -
               (0.007 * Lotjonen_sd) - (0.108 * Lotjonen_ln),
             Lotjonen_Sleep_2 = factor(ifelse(Lotjonen_Sleep_2 < 0.0, "Wake", "Sleep")),
             Noon_Date = as.character(date(DateTime + hours(12))))
    
    # Sleep Smoothing
    ii %<>%
      mutate(
        Sleep_Thresh_Smooth = smooth_sleep(., "Sleep_Thresh"),
        Sleep_Acti_Smooth = smooth_sleep(., "Sleep_Acti"),
        Lotjonen_Sleep_Smooth = smooth_sleep(., "Lotjonen_Sleep"),
        Lotjonen_Sleep_2_Smooth = smooth_sleep(., "Lotjonen_Sleep_2")
      )
    
    # Calculate Noon-Noon days
    which_day_date <- table(ii[, "Day"], ii[, "Date"]) %>%
      as.data.frame.table(.) %>%
      dplyr::filter(Freq > 0) %>%
      rename(Noon_Day = Var1, Noon_Date = Var2) %>%
      dplyr::select(-Freq)
    
    ii %<>% merge(., which_day_date, by = "Noon_Date") %>%
      select(-Noon_Date) %>%
      mutate(Noon_Day = opelr::numfac(Noon_Day),
             Noon_Day = Noon_Day - (min(Noon_Day) - 1))

    return(ii)
  }) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.)) %>%
  dplyr::filter(!is.na(Activity))

## ------ Light Adherence ------

# Count periods where median Light window > 10,000/5,000 lux/m2 across 15/30 epochs
actigraphy <- split(actigraphy, actigraphy$patient_ID) %>%
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

## ------ Activity Windowing - EE ------

actigraphy <- split(actigraphy, actigraphy$patient_ID) %>%
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

## ------ Filter Bad Recordings ------

actigraphy %<>%
  filter(!patient_ID %in% c("Matt_Wickham", "Rick_Teutsch"))

## ------ Saving RDS ------

saveRDS(acti_files, ".\\Rmd\\Data\\actigraphy_header.rds")
saveRDS(actigraphy, ".\\Rmd\\Data\\actigraphy_data.rds")
