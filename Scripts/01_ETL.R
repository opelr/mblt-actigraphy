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
         Light = opelr::numfac(Light),
         DAR_Period = factor(ifelse(DateTime < paste(Date, "07:00:00 PDT") |
                                    DateTime > paste(Date, "21:59:59 PDT"),
                                    "Night", "Day")),
         log_Activity = log10(Activity),
         log_Light = log10(Light))

## ------ Rolling Window/Off-Wrist Detection ------

moving_window <- function(data, window, step, FUN, align) {
  if (! align %in% c("left", "right", "center")) stop('"align" must be "left", "right", or "center"')
  if (align == "center" & window %% 2 == 0) stop("'window' must be odd when 'align = center'")
  
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
}

dat <- dplyr::filter(actigraphy, patient_ID == "Carolyn_Jones", Day == "4")

ggplot(dat) +
  # geom_line(aes(x = Time, y = Light), color = "Orange") +
  geom_line(aes(x = Time, y = Activity), color = "Black") +
  scale_x_datetime(date_labels = "%I %p") +
  # scale_y_log10() +
  facet_grid(Day + DateAbbr ~ .)

dat %<>% 
  mutate(
    Rolling_Activity_SD = do.call("c", lapply(unique(dat$patient_ID), function(ii) {
      moving_window(dat$Activity[dat$patient_ID == ii], 7, 1, FUN = "sd", "center")})),
    Rolling_Activity_mean = do.call("c", lapply(unique(dat$patient_ID), function(ii) {
      moving_window(dat$Activity[dat$patient_ID == ii], 15, 1, FUN = "mean", "center")})),
    Sleep = factor(ifelse(Rolling_Activity_mean < 100, "Sleep", "Wake")),
    Activity_diff = do.call("c", lapply(unique(dat$patient_ID), function(ii) {
      c(NaN, diff(dat$Activity[dat$patient_ID == ii], 1))})))

## ------ Lotjonen ------

dat %<>% 
  mutate(
    Lotjonen_mean = do.call("c", lapply(unique(dat$patient_ID), function(ii) {
      moving_window(dat$Activity[dat$patient_ID == ii], 15, 1, FUN = "mean", "center")})),
    Lotjonen_sd = do.call("c", lapply(unique(dat$patient_ID), function(ii) {
      moving_window(dat$Activity[dat$patient_ID == ii], 17, 1, FUN = "sd", "center")})),
    Lotjonen_ln = log(Activity) + 0.1,
    Lotjonen_Counts = Activity > 10,
    Lotjonen_nat = do.call("c", lapply(unique(dat$patient_ID), function(ii) {
      moving_window(dat$Lotjonen_Counts[dat$patient_ID == ii], 23, 1, FUN = "sum", "center")})))

dat[, c("Rolling_Activity_mean", "Rolling_Activity_SD", "Sleep", "Activity_diff")]

## How can I be certain this person isn't just asleep during the day?
## compare to sleep diary? What if they don't complete one?
## If Activity Diff is 0 (consistent?) for j number of epochs, call it watch off

## ------ Saving RDS ------

saveRDS(acti_files, ".\\Data\\actigraphy_header.rds")
saveRDS(actigraphy, ".\\Data\\actigraphy_data.rds")
