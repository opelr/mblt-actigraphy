## Script: Import and Shaping of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-12
## Version: 1.1.0

library(magrittr)
library(tidyverse)
library(lubridate)
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
  
  header_info <- data.frame(patient_name = patient_name) %>%
    mutate(patient_sex = patient_sex) %>%
    mutate(patient_DOB = as.POSIXct(strptime(patient_DOB,
                                             format = "%Y-%m-%d"))) %>%
    mutate(patient_age = patient_age) %>%
    mutate(recording_date_start = as.POSIXct(strptime(recording_date_start,
                                                      format = "%Y-%m-%d"))) %>%
    mutate(recording_date_end = as.POSIXct(strptime(recording_date_end,
                                                    format = "%Y-%m-%d"))) %>%
    mutate(recording_days_span = recording_days_span) %>%
    mutate(recording_epoch_length = recording_epoch_length) %>%
    mutate(watch_serial_no = watch_serial_no) %>%
    mutate(calibration_activity = calibration_activity) %>%
    mutate(calibration_light = calibration_light) %>%
    mutate(threshold_wake = threshold_wake) %>%
    mutate(threshold_light = threshold_light)
  
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
  
  acti$ID = paste0(nam, collapse = "_")
  
  ## Date/Time Calculations
  acti$DateTime = as.POSIXct(strptime(x = paste(acti$Date, acti$Time),
                                          format = "%Y-%m-%d %I:%M:%S %p"))
  acti$Time = as.POSIXct(strptime(acti$Time, format = "%I:%M:%S %p"))
  
  acti %<>%
    mutate(Hour = hour(DateTime)) %>%
    mutate(AM_PM = factor(floor(Hour/12), levels = 0:1,
                          labels = c("AM", "PM"))) %>%
    mutate(Weekend = factor(chron::is.weekend(DateTime),
                            labels = c("Weekend", "Weekday")))
    
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
    mutate(Sunrise = gsub(" [0-9]{4}", "", Value)) %>%
    mutate(Sunset = gsub("[0-9]{4} ", "", Value)) %>%
    mutate(Date = as.POSIXct(strptime(paste0(Month, Day, year(datetime)),
                                      format = "%b%d%Y"))) %>%
    mutate(Sunrise = as.POSIXct(strptime(paste0(Month, Day, year(datetime),
                                                Sunrise), format = "%b%d%Y%H%M"))) %>%
    mutate(Sunset = as.POSIXct(strptime(paste0(Month, Day, year(datetime),
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

actigraphy <- merge(actigraphy, sunsetDF, by = "Date") %>%
  mutate(SunPeriod = factor(ifelse(DateTime <= Sunrise, "Dawn",
                            ifelse(DateTime > Sunrise & DateTime < Sunset, "Day",
                            ifelse(DateTime >= Sunset, "Night", NA))))) %>%
  dplyr::arrange(., ID, DateTime) %>%
  mutate(Activity = opelr::numfac(Activity),
         Light = opelr::numfac(Light),
         DAR_Period = factor(ifelse(DateTime < paste(Date, "07:00:00 PDT") |
                                    DateTime > paste(Date, "21:59:59 PDT"),
                                    "Night", "Day")))

## ------ Saving RDS ------

saveRDS(acti_files, ".\\Data\\actigraphy_header.rds")
saveRDS(actigraphy, ".\\Data\\actigraphy_data.rds")
