## Script: Functions for ETL, Analysis, and Plotting of Philips Actiware Data
## Project: mblt-actigraphy
## Author: Ryan Opel -- @opelr
## Date: 2017-11-15
## Version: 0.3.0

# ------------------------ ETL ------------------------

## ------ Major ETL Functions ------

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
  
  # Patient name from rootpath
  patient_name <- gsub("_(.*)", "", gsub("(.*)\\\\", "", path))
  
  
  # Strip header information to temporary variables
  # patient_name <- as.character(acti$X2[acti$X1 == "Identity:"])
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
  header_info <- data.frame(watch_ID = gsub(" ", "_", patient_name)) %>%
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
  acti$watch_ID <- strsplit(nam, "_")[[1]][1]
  
  # acti$watch_ID = paste0(nam, collapse = "_")
  
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

## ------ Survey Data Summation ------

survSums <- function(cols) {
  ifelse(apply(PCL[, cols], 1, FUN = function(x) all(is.na(x))) == TRUE,
         NA, rowSums(PCL[, cols], na.rm = T))
}


## ------ Patient Chronotype from rMEP ------

strip_waketime <- function(x) {
  wake_time <- strsplit(as.character(x), "-")[[1]][1] %>%
    # unlist %>%
    # .[seq(1, length(.) - 1, 2)] %>%
    strptime(., "%I:%M %p", tz = "America/Los_Angeles") %>%
    as.character %>%
    substr(., 12, nchar(.))
  
  if (is.na(wake_time)) {
    return(NA)
  } else {
    return(paste0(wake_time, " PST"))
  }
}

waketime_plus <- function(x, n) {
  if (is.na(x)) {
    return(NA)
  } else {
    
    hours_out <- strsplit(x, ":")[[1]][1] %>%
      as.numeric %>%
      add(., n)
    
    return(gsub(strsplit(x, ":")[[1]][1], hours_out, x))
  }
}

## ------ DAR Calculation ------

DAR_ifelse <- function(df, morning, evening, true_label, false_label) {
  
  factor(ifelse((is.na(df[, morning]) | is.na(df[, evening])), NA, 
                ifelse(df[, "DateTime"] >= as.POSIXct(paste(df[, "Date"], df[, morning])) &
                         df[, "DateTime"] <= as.POSIXct(paste(df[, "Date"], df[, evening])),
                       true_label, false_label)))
}

## ------ Sleep Staging Functions ------

actiware_sleep <- function(df, column = "Activity", threshold) {
  # Replicates Actiware's threshold-based sleep-staging algorithm
  # 
  # Args:
  #   dataframe (df): Data frame to iterate upon; must be single patient to work correctly
  #   column (str): Which column will we be acting on
  #   threshold (int): Lower bound for levels to be considered "Wake"
  # 
  # Returns:
  #   A same-length vector of Sleep and Wake values
  #
  # Example:
  #   actiware_sleep(dplyr::filter(actigraphy, watch_ID == "A26598.1"), threshold = 40))
  
  x <- df[, column]
  
  calc_values <- sapply(2:(length(x) - 1), function(ii) {
    y = (x[ii - 1] * 0.12) + (x[ii] * 0.5)  + (x[ii + 1] * 0.12)
  }) %>%
    c(NA, ., NA)
  
  sleep_stats <- factor(ifelse(calc_values > threshold, "Wake", "Sleep"))
  
  return(sleep_stats)
}

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
  
  do.call("c", lapply(unique(df[, "watch_ID"]), function(ii) {
    data <- df[df[, "watch_ID"] == ii, column]
    
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
  x <- round(x)
  x <- range(x, na.rm = T) / mean(x, na.rm = T)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

threshold_range <- function(x, thresh = 7) {
  #' Test for sub-threshold variance in a vector
  #' 
  #' Args:
  #'   x (int): Vector of numbers
  #'   thresh (int): Maximum difference between the extremes of 'x'
  #' 
  #' Returns:
  #'   Boolean value
  
  if (length(x) == 1)
    return(TRUE)
  
  (max(x, na.rm = T) - min(x, na.rm = T)) <= thresh
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

get_rle_df <- function(x) {
  # Convert vector to RLE object to data.frame
  data.frame(Values = rle(as.character(x))$values,
             Lengths = rle(as.character(x))$lengths)
}

smooth_sleep <- function(data, column) {
  # Recursively smooth sleep staging
  # 
  # Args:
  #   data (df): Data frame to iterate upon
  #   column (str): Name of column to call smoothing functions on
  # 
  # Returns:
  #   Vector containing 'Sleep' and 'Wake', same length as 'column'
  
  # Convert RLE data.frame to a vector
  rle_df_to_vector <- function(df) {rep(df$Values, df$Lengths)}
  
  # Recursively apply sleep smoothing to RLE DF object
  cole_post_process <- function(df) {
    
    # Copy data for comparison
    df_old <- df
    
    for (ii in 1:nrow(df)) {
      if (is.null(df$Values[ii])) {
        df$Values[ii] <- NA
      } else if (is.na(df$Values[ii])) {
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

## ------ Consolidated Sleep/Wake Bouts ------

get_sleep_bouts <- function(x, threshold) {
  # Determine if span of sleep/wake constitutes a consolidated bout
  # 
  # Args:
  #   x (vector): Vector containing values "Sleep" and "Wake"
  #   threshold (int): Number of continous bouts required to constitute a bout
  # 
  # Returns:
  #   Vector containing 'Sleep_Bout' and 'Wake_Bout', same length as 'x'
  
  # Run smoothing logic on RLE DF
  sleep_rle <- get_rle_df(x)
  
  
  sleep_bouts <- sapply(1:nrow(sleep_rle), function (jj) {
    if (sleep_rle$Values[jj] == "Wake" & sleep_rle$Lengths[jj] >= threshold) {
      o <- "Wake_Bout"
    } else if (sleep_rle$Values[jj] == "Sleep" & sleep_rle$Lengths[jj] >= threshold) {
      o <- "Sleep_Bout"
    } else {
      o <- NA
    }
  }) %>%
    na.locf(., na.rm = F)
  
  sleep_rle["Sleep_Bouts"] <- sleep_bouts
  
  output <- rep(sleep_rle$Sleep_Bouts, sleep_rle$Lengths)
  return(output)
}

## ------ Clean Light Adherence ------

clean_adherence <- function(data, column) {
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
      if (is.null(df$Values[ii])) {
        df$Values[ii] <- NA
      } else if (is.na(df$Values[ii])) {
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
  
  rle_df <- get_rle_df(data[, column])
  
  out <- post_process_light(rle_df) %>%
    rle_df_to_vector
  return(out)
}



# ------------------------ Process ------------------------

## ------ Rhythm Stability ------

### Intradaily Variability
calculate_IV <- function(var, patient, window, size, dat_frm) {
  #' Calculate Intradaily Variability for a non-parametric fake_actig variable
  #' Quantifies how fragmented the rhythm is relative to the overall variance
  #' 
  #' Args:
  #'   var (str): Column name present in your 'fake_actig' dataframe
  #'   patient (int): Patient ID number
  #'   window (str): Options are "moving" or "expanding"
  #'   size (int): If window is set to "moving", this is the number of days
  #'               considered. If "expanding", specifies the buffer before
  #'               initial IS is calcualted.
  #'   dat_frm (data.frame): Data frame being used
  #' 
  #' Returns:
  #'   Data frame containing IV for every unique patient in your fake_actig data frame.
  #'   Values reach near zero for a perfect sine wave, around 2 for Gaussian noise,
  #'   and may be higher when a definite ultradian component with a period of
  #'   2h is present.
  #'   
  #' Reference:
  #'   Van Someren, Chronobiology International, 1999
  
  ## Helper fn
  fp <- function(x, y) {
    # Shorthand notation for 'as.formula(paste(.))'
    return(as.formula(paste(x, "~", y)))
  }
  
  ### Local Variables
  dat <- dat_frm
  pat_unique_days <- unique(dat$Noon_Day[dat$patient_ID == patient])
  pat_seq_days <- seq(pat_unique_days[1], pat_unique_days[length(pat_unique_days)], 1)
  
  IV <- lapply(size:length(pat_unique_days), function (ii) {
    
    if (window == "expanding") {
      dat <- dat[dat$patient_ID == patient &
                   dat$Noon_Day %in% pat_unique_days[1]:pat_unique_days[ii], ]
    } else if (window == "moving") {
      dat <- dat[dat$patient_ID == patient &
                   dat$Noon_Day %in% pat_unique_days[ii - size + 1]:pat_unique_days[ii], ]
    }
    
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

### Interdaily Stability
calculate_IS <- function(var, patient, window, size, dat_frm) {
  #' Calculate Interdaily Stability for a non-parametric fake_actig variable
  #' The extent to which the profiles of individual days resemble each other
  #' 
  #' Args:
  #'   var (str): Column name present in your 'fake_actig' dataframe
  #'   patient (int): Patient ID number
  #'   window (str): Options are "moving" or "expanding"
  #'   size (int): If window is set to "moving", this is the number of days
  #'               considered. If "expanding", specifies the buffer before
  #'               initial IS is calcualted.
  #'   dat_frm (data.frame): Data frame being used
  #' 
  #' Returns:
  #'   Data frame containing IS for every unique patient in your fake_actig data frame.
  #'   Values range from 0 (Guassian noise) to 1 (perfect IS).
  #'   
  #' Reference:
  #'   Van Someren, Chronobiology International, 1999
  
  ## Helper fn
  fp <- function(x, y) {
    # Shorthand notation for 'as.formula(paste(.))'
    return(as.formula(paste(x, "~", y)))
  }
  
  ### Local Variables
  dat <- dat_frm
  pat_unique_days <- unique(dat$Noon_Day[dat$patient_ID == patient])
  pat_seq_days <- seq(pat_unique_days[1], pat_unique_days[length(pat_unique_days)], 1)
  
  IS <- lapply(size:length(pat_unique_days), function (ii) {
    
    if (window == "expanding") {
      dat <- dat[dat$patient_ID == patient &
                   dat$Noon_Day %in% pat_unique_days[1]:pat_unique_days[ii], ]
    } else if (window == "moving") {
      dat <- dat[dat$patient_ID == patient &
                   dat$Noon_Day %in% pat_unique_days[ii - size + 1]:pat_unique_days[ii], ]
    }
    
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
    total_mean <- merge(dat[, c("patient_ID", "Epoch", var)] %>%
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

## ------ Relative Amplitude (RA), M10, and L5 ------

calc_rest_phase <- function(len, fn = 'max') {
  hours <- unique(actigraphy$Hour)[1:(length(unique(actigraphy$Hour)) - (len - 1))]
  
  rest_phase <- lapply(hours, function (i) {
    out <- aggregate(Activity ~ patient_ID + Noon_Day,
                     dplyr::filter(actigraphy, Hour >= i, Hour <= i + (len - 1)),
                     sum) %>%
      rename(Activity_Sum = Activity) %>%
      mutate(Hours = paste0(i, "-", i + (len - 1)),
             Activity_Mean = Activity_Sum / len)
  }) %>%
    do.call("rbind", .)
  
  rp2 <- aggregate(Activity_Mean ~ patient_ID + Noon_Day, rest_phase, function(j) {f <- get(fn); f(j)}) %>%
    merge(., rest_phase, by = c("patient_ID", "Noon_Day", "Activity_Mean")) %>%
    rename_(.dots=setNames(names(.), tolower(gsub("Hours",
                                                  paste0(fn, len, "Hours"),
                                                  names(.))))) %>%
    rename(patient_ID = patient_id, Noon_Day = noon_day,
           Activity_Mean = activity_mean) %>%
    arrange(patient_ID, Noon_Day) %>%
    select(-activity_sum)
  
  return(rp2)
}

## ------ Activity Consolidation ------

### Consolidated Bouts
calc_activity_consolidation <- function(column) {
  
  get_rle_df <- function(x) {
    data.frame(Values = rle(as.character(x))$values,
               Lengths = rle(as.character(x))$lengths)
  }
  
  patient_dates <- xtabs(~ patient_ID + Noon_Day, actigraphy) %>%
    as.data.frame %>%
    filter(Freq > 0) %>%
    rename(Epochs = Freq) %>%
    arrange(patient_ID, Noon_Day)
  
  activity_consolidation_parent <- lapply(1:nrow(patient_dates), function(ii) {
    x <- actigraphy[actigraphy$patient_ID == patient_dates$patient_ID[ii] &
                      actigraphy$Noon_Day == patient_dates$Noon_Day[ii], column]
    
    cbind(get_rle_df(x), patient_dates$patient_ID[ii], patient_dates$Noon_Day[ii])
  }) %>%
    do.call("rbind", .) %>%
    dplyr::filter(complete.cases(.)) %>%
    set_colnames(c("Values", "Lengths", "patient_ID", "Noon_Day"))
  
  activity_consolidation <- merge(aggregate(Lengths ~ patient_ID + Noon_Day + Values,
                                            activity_consolidation_parent, sum),
                                  aggregate(Lengths ~ patient_ID + Noon_Day + Values,
                                            activity_consolidation_parent, length) %>%
                                    rename(Bouts = Lengths),
                                  by = c("patient_ID", "Noon_Day", "Values")) %>%
    filter(Values == TRUE) %>%
    arrange(patient_ID, Noon_Day) %>%
    mutate(Minutes = 30 * Bouts + 2 * (Lengths - Bouts),
           ABL = Minutes/Bouts)
  
  return(activity_consolidation)
}

## ------ Light Adherence ------

calc_light_adherence <- function(column) {
  
  get_rle_df <- function(x) {
    data.frame(Values = rle(as.character(x))$values,
               Lengths = rle(as.character(x))$lengths)
  }
  
  patient_dates <- xtabs(~ patient_ID + Noon_Day, actigraphy) %>%
    as.data.frame %>%
    filter(Freq > 0) %>%
    rename(Epochs = Freq) %>%
    arrange(patient_ID, Noon_Day)
  
  light_adherence_parent <- lapply(1:nrow(patient_dates), function(ii) {
    pat_ID <- patient_dates$patient_ID[ii]
    noon_day <- patient_dates$Noon_Day[ii]
    
    dat <- filter(actigraphy, patient_ID == pat_ID, Noon_Day == noon_day) %>%
      select(starts_with(column)) %>%
      unlist %>% 
      get_rle_df %>%
      data.frame(patient_ID = pat_ID, Noon_Day = noon_day, .)
    
    return(dat)
  }) %>%
    do.call("rbind", .) %>%
    dplyr::filter(complete.cases(.))
  
  light_adherence <- merge(aggregate(Lengths ~ patient_ID + Noon_Day + Values,
                                     light_adherence_parent, sum),
                           aggregate(Lengths ~ patient_ID + Noon_Day + Values, light_adherence_parent,
                                     length) %>%
                             rename(Bouts = Lengths),
                           by = c("patient_ID", "Noon_Day", "Values")) %>%
    # mutate(Values = as.logical(as.numeric(Values))) %>%
    filter(Values == TRUE) %>%
    arrange(patient_ID, Noon_Day) %>%
    mutate(Minutes = 30 * Bouts + 2 * (Lengths - Bouts),
           ABL = Minutes/Bouts)
  
  return(light_adherence)
}

## ------ Max/Min Week by Light ------

light_week <- function(pat_ID, len, fn = 'max') {
  df <- filter(results$major_results, patient_ID == pat_ID)
  days <- df$Noon_Day[1:(length(df$Noon_Day) - (len - 1))]
  
  light_days <- lapply(days, function (i) {
    out <- aggregate(Light ~ patient_ID,
                     dplyr::filter(df, Noon_Day >= i, Noon_Day <= i + (len - 1)),
                     sum) %>%
      rename(Light_Sum = Light) %>%
      mutate(Day_Range = paste0(i, "-", i + (len - 1)),
             Light_Mean = Light_Sum / len)
  }) %>%
    do.call("rbind", .)
  
  if (fn == 'max') {
    day_range <- light_days$Day_Range[which.max(light_days$Light_Sum)]
  } else if (fn == 'min') {
    day_range <- light_days$Day_Range[which.min(light_days$Light_Sum)]
  }
  
  return(data.frame(patient_ID = pat_ID, Days = day_range))
}

# ------------------------ Visualize ------------------------ 

## ------ Entire Actogram ------

plot_patient <- function(patient, df) {
  d1 <- dplyr::filter(df, patient_ID == patient) %>%
    select(., Time, Light, Activity, Day, DateAbbr, Sleep_Acti_Smooth, Sleep_Bouts) %>%
    mutate(Activity_Scale = Activity * (max(Light, na.rm = T) / max(Activity, na.rm = T)),
           Log_Light = log10(Light))
  
  ## Rectangle shading
  bouts <- data.frame(Lengths = rle(d1$Sleep_Bouts)$lengths,
                      Values = rle(d1$Sleep_Bouts)$values) %>%
    mutate(x2 = cumsum(Lengths),
           x1 = lag(x2, 1) + 1,
           x1 = ifelse(is.na(x1), 1, x1)) %>%
    filter(Values == "Sleep_Bout") %>%
    select(Lengths, Values, x1, x2)
  
  time_bouts <- data.frame(x1 = d1$Time[bouts$x1],
                           x2 = d1$Time[bouts$x2],
                           y1 = 0, y2 = Inf)
  
  p <- ggplot(d1, mapping = aes(x = Time)) +
    geom_line(aes(y = Light), color = "Orange") +
    geom_line(aes(y = Activity_Scale), color = "Black") +
    # geom_rect(data=time_bouts, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
    # color='blue', alpha=0.2) +
    facet_grid(Day + DateAbbr ~ .) +
    scale_x_datetime(date_labels = "%I %p", date_minor_breaks = "1 hour") + 
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . * (max(d1$Activity, na.rm = T) / max(d1$Light, na.rm = T)),
                                           name = "Activity")) +
    ylab("Light") + 
    labs(title = patient)
  
  plot(p)
}

## ------ Average Day Activity ------

plot_avg_patient_day <- function(patient, df) {
  d1 <- dplyr::filter(df, patient_ID == patient) %>%
    select(., Time, Light, Activity, Day, DateAbbr, Sleep_Acti_Smooth, Sleep_Bouts) %>%
    mutate(Activity_Scale = Activity * (max(Light, na.rm = T) / max(Activity, na.rm = T)),
           Log_Light = log10(Light)) %>%
    aggregate(cbind(Light, Activity, Activity_Scale, Log_Light) ~ Time, ., mean)
  
  p <- ggplot(d1, mapping = aes(x = Time)) +
    geom_line(aes(y = Light), color = "Orange") +
    geom_line(aes(y = Activity_Scale), color = "Black") +
    # geom_rect(data=time_bouts, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
    # color='blue', alpha=0.2) +
    scale_x_datetime(date_labels = "%I %p", date_minor_breaks = "1 hour") + 
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . * (max(d1$Activity, na.rm = T) / max(d1$Light, na.rm = T)),
                                           name = "Activity (Black)")) +
    ylab("Light (Yellow)") + 
    labs(title = paste0("Subject #", patient))
  
  plot(p)
}

## ------ All Patients Avg Activity, Overlayed ------

plot_avg_all_patients  <- function(variable, df) {
  d1 <- dplyr::filter(df) %>%
    mutate(Activity_Scale = Activity * (max(Light, na.rm = T) / max(Activity, na.rm = T)),
           Log_Light = log10(Light),
           patient_ID = factor(patient_ID)) %>%
    aggregate(cbind(Light, Activity, Activity_Scale, Log_Light) ~ patient_ID + Time, ., mean)
  
  p <- ggplot(d1, aes_string("Time", variable, group = "patient_ID",
                             color = "patient_ID")) +
    geom_line() + 
    scale_x_datetime(date_labels = "%I %p", date_minor_breaks = "1 hour")
  
  plot(p)
}

## ------ Slopegraphs ------

slopegraph <- function(variable, dat_frm) {
  p <- ggplot(dat_frm,
              aes_string("Week_Name", variable, group = "patient_ID",
                         color = "patient_ID")) +
    geom_point(aes(shape = Light), size = 3) + 
    geom_line(aes(linetype = Light), size = 1) +
    facet_grid(PTSD ~ .)
  
  plot(p)
}

## ------ Longitudinal Plots ------

longitudinal_plot <- function(var) {
  ggplot(results$major_results, aes_string("Noon_Day", var, colour="patient_ID")) +
    geom_point(aes(shape = MBLT_Group), size = 2) +
    geom_smooth(aes(linetype = MBLT_Group)) +
    facet_grid(Group_PTSD ~ .)
}

# ------------------------ Model/Analyze ------------------------ 

# ------------------------ Helper/Auxiliary Functions ------------------------

## ------ Find Closest ------
find_closest <- function(x, num) {
  #' Find "Price is Right" style number
  x_sort <- x[order(x)]
  i <- findInterval(num, x_sort)
  return(x_sort[i])
}