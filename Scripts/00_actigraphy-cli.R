## Script: Actigraphy Pipeline Command-Line Interface
## Project: mblt-actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2018-04-18
## Version: 1.0.0

# TODO: If patient is updated, pull their data from the old data frames
#       before appending new data

#' This script acts as a general interface between the bash script `actigraphy_csv_monitor`
#' and the R processing scripts `01_...` - `02_...`.
#'
#' The data-processing pipeline is detailed below:
#'   1. Run the `actigraphy_csv_monitor` bash script on a 24-hour cron job. This
#'      script monitors the directory containing raw actigraphy CSV files. It gets
#'      the names of files modified in the past 24 hours, and passes them to this script;
#'   2. The final products of this pipeline are four RDS objects. When this script is run,
#'      it checks if those files exist. If they do not, it ignores the command-line arguments
#'      and runs through the entire database of CSV files;
#'   3. If those RDS files exist, we use the command-line inputs and filter the `acti_files` object to only include them;
#'   4. Run patients through to the end, normally;
#'   5. Load old versions of RDS objects (assuming they exist), and rbind new patients on;
#'   6. Save RDS object;
#'   7. Repeat process with the other scripts.

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

## ------ CLI Options ------

option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL,
              help="dataset files", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

### Handle CLI Input
input_string <- opt$file

if (!is.null(input_string)) {
  new_files <- data.frame(Files = strsplit(input_string, "---")) %>%
    set_colnames("Files")
}

## ------ File Repository ------
#' Build repository of file paths to be read later by two functions.
#' If CLI arguments were passed in, filter dataframe.

FILE_PATH <- "D:/data/acquired_data/human/h4085/Actigraphy/CSV/"
FILE_MASK <- "(.*)_Bedtime.csv"

acti_files <- data.frame(File = list.files(path = FILE_PATH,
                                           pattern = FILE_MASK)) %>%
  mutate(rootpath = normalizePath(paste0(FILE_PATH, File))) %>%
  filter(!grepl('-', File))

if (exists("new_files")) {
  acti_files <- acti_files %>%
    filter(File %in% new_files$Files)
}

## ------ Load Old RDS Objects ------

#' Do old RDS objects exist?
#' If no: Run through all patients
#' If yes and new_files: Proceed
#' If yes and no new_files: quit

tryRDS <- function(path) {
  if (file.exists(path)) {
    readRDS(path)
  } else {
    return(FALSE)
  }
}

actigraphy_headers_old <- tryRDS("./Rmd/Data/actigraphy_header.rds")
actigraphy_static_old  <- tryRDS("./Rmd/Data/actigraphy_static.rds")
actigraphy_filtered_old <- tryRDS("./Rmd/Data/actigraphy_filtered.rds")
actigraphy_results_old <- tryRDS("./Rmd/Data/actigraphy_results.rds")

## All exist and no new files
if (is.data.frame(actigraphy_headers_old) &
    is.data.frame(actigraphy_static_old) &
    is.data.frame(actigraphy_results_old) &
    (is.null(input_string)) | input_string == " ") {
  message("Message: No new files to add... quitting")
  quit(save="no")
}

## ------ Source ETL Script and Save Concat. RDS files ------

try_save_RDS <- function(new_obj, old_obj, column, path) {
  if (is.logical(old_obj)) {
    saveRDS(new_obj, path)
  } else {
    rows_select <- !(old_obj[, column] %in% unique(new_obj[, column]))
    old_obj <- old_obj[rows_select, ]
    saveRDS(rbind(old_obj, new_obj), path)
  }
}

message("Message: Source... Script #1")
source("Scripts/01_ETL.R")

## Save RDS -- Headers
try_save_RDS(acti_files, actigraphy_headers_old, "File", "./Rmd/Data/actigraphy_header.rds")
## Save RDS -- Static Data
try_save_RDS(actigraphy, actigraphy_static_old, "patient_ID", "./Rmd/Data/actigraphy_static.rds")
message("Message: Saved... Headers & Static DataFrame")

## Validation
validate_actigraphy <- readRDS("./Rmd/Data/actigraphy_static.rds")

assertthat::assert_that(all(table(validate_actigraphy$patient_ID, validate_actigraphy$Noon_Day) <= 720),
                        msg = "Error: Some patients have more than 720 observations per day")

## ------ Off Wrist Filtering ------

#' Exclude days that contain any 3+ hour window of inactivity
#' (`Thresh_Activity_Change_Window`).

# actigraphy <- readRDS("./Rmd/Data/actigraphy_static.rds")

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

## ------ Source Proprietary Additions ------

message("Message: Source... Script #2")
source("Scripts/Proprietary_Additions.R")

## Save RDS
try_save_RDS(actigraphy, actigraphy_filtered_old, "patient_ID", "./Rmd/Data/actigraphy_filtered.rds")
message("Message: Saved... Filtered DataFrame")

## ------ Source Results Script ------

message("Message: Source... Script #3")
source("Scripts/02_Process.R")

# rbind each dataframe in the results list object

if (!is.logical(actigraphy_results_old)) {
  results_bind <- list()
  for (nam in names(results)) {
    new <- results[[nam]]
    old <- actigraphy_results_old[[nam]]
    
    rows_select <- !(old[, "patient_ID"] %in% unique(new[, "patient_ID"]))
    old <- old[rows_select, ]
    
    results_bind[[nam]] <- rbind(old, new)
  }
} else {
  results_bind <- results
}

saveRDS(results_bind, "./Rmd/Data/actigraphy_results.rds")
message("Message: Saved... Results")
message("Finished: goodbye")
