## Script: Inclusion of Lab-Specific Metrics (REDCap, Latency Estimations, etc)
## Project: mblt-actigraphy
## Author: Ryan Opel -- @opelr
## Date: 2018-04-18
## Version: 1.0.1

suppressMessages(library(magrittr))
suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(reshape2))
suppressMessages(library(zoo))

#' In previous versions of this repository, there were many code snippets that
#' were not only a) of exclusive interest to the lab, but b) caused the scripts
#' to break if our proprietary information was not present.
#' 
#' This script collates all those bits and pieces, making the rest of the
#' repository more friendly to collaborators, other users, etc.
#' 
#' This repository is still not perfectly generalized, and probably never will
#' be. Users will continue to change certain file-paths and -masks, however
#' with the changes made by this script/branch, it is my hope that the
#' repository will *at least* not break when run outside of our development
#' environment.

# source("Scripts/Functions.R")
# actigraphy <- readRDS(".\\Rmd\\Data\\actigraphy_filtered.rds")

## ------ PTSD Status from PCL-5 ------

PCL <- read.csv("data/raw/4085_Baseline_PCL.csv") %>% 
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

## ------ Hand-Added Medical History -----

#' FIXME: This doesn't belong in this script
#' TODO: Filter out withdrawn patients ?

CCT_medHist <- data.frame(patient_ID = c(10, 11, 12, 13, 14, 15, 16),
                          Sex = c("F", "M", "F", "F", "M", "M", "M"),
                          PTSD_CCT = c(T, T, F, T, F, T, T),
                          TBI_CCT = c(F, F, F, F, F, T, F),
                          Depression_CCT = c(F, F, F, F, T, F, F),
                          Anxiety_CCT = c(T, F, F, F, F, F, F),
                          OSA_CCT = c(F, F, F, F, F, T, T))

## ------ Patient Chronotype from rMEQ ------
#' Create DAR estimates (in following section) based on published criteria

rMEQ <- read.csv("data/raw/4085_Baseline_rMEQ.csv") %>%
  set_colnames(c("patient_ID", "Event_Name", "Wake_Time", "Refreshed",
                 "Bed_Time", "Best_Time", "Chronotype", "Complete")) %>%
  mutate(DAR_Wake_Time = mapply(Wake_Time, FUN = strip_waketime),
         DAR_Bed_Time = mapply(DAR_Wake_Time, FUN = function(i) waketime_plus(i, 12)),
         DAR_Light_Box_Time = mapply(DAR_Wake_Time, FUN = function(i) waketime_plus(i, 3))) %>%
  select(patient_ID, DAR_Wake_Time, DAR_Bed_Time, DAR_Light_Box_Time)

## ------ Mutate Main Dataframe ------

actigraphy %<>%
  merge(., PCL[, c("patient_ID", "Group_PTSD")], by = "patient_ID", all.x = T) %>%
  merge(., rMEQ, by = "patient_ID", all.x = T) %>%
  mutate(DAR = DAR_ifelse(., "DAR_Wake_Time", "DAR_Bed_Time", "Day", "Night"),
         Light_Box_Period = DAR_ifelse(., "DAR_Wake_Time", "DAR_Light_Box_Time", "Light", "Non_Light")) %>%
  merge(., CCT_medHist, "patient_ID", all.x = T) %>%
  arrange(patient_ID, DateTime)

rm(CCT_medHist, rMEQ)
## ------ Save RDS ------

# saveRDS(actigraphy, "./data/rds/actigraphy_filtered.rds")
