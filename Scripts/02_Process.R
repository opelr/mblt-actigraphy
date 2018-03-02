## Script: Analysis of Actiwatch Data
## Project: mblt-actigraphy
## Author: Ryan Opel -- @opelr
## Date: 2018-01-05
## Version: 1.0.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(zoo)

source("Scripts/Functions.R")

actigraphy <- readRDS(".\\Rmd\\Data\\actigraphy_filtered.rds")
results <- list()

## ------ Basic Sleep Metrics ------

#' Calculating Total Sleep Time (TST), sleep efficiency (SE), and
#' Wake After Sleep Onset (WASO)

results$sleep_metrics <- xtabs(~ patient_ID + Interval + Sleep_Acti_Smooth + Noon_Day,
                               data = actigraphy) %>%
  as.data.frame %>%
  dcast(., patient_ID + Noon_Day + Interval ~ Sleep_Acti_Smooth,
        value.var = "Freq") %>%
  filter(Interval == "Sleeping", (Sleep + Wake) > 0) %>%
  mutate(Sleep_Minutes = Sleep * 2,
         WASO_Minutes = Wake * 2,
         SE = 100 * Sleep / (Wake + Sleep))

## ------ Sleep Latency ------

results$sleep_latency <- xtabs(~ patient_ID + Interval + Noon_Day,
                               data = actigraphy) %>%
  as.data.frame %>%
  filter(Interval %in% c("Falling_Asleep", "Waking_Up")) %>%
  dcast(., patient_ID + Noon_Day ~ Interval, value.var = "Freq") %>%
  mutate(Sleep_Latency_Minutes = Falling_Asleep * 2,
         Wake_Latency_Minutes = Waking_Up * 2) %>%
  select(-Falling_Asleep, -Waking_Up)

## ------ Activity Rhythm Stability ------

#' Calculating IV and IS
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

#' DAR_77 == 7AM - 7PM
#' Additional DAR variables can be added in 01_ETL.R script

results$DAR_77 <- aggregate(Activity ~ patient_ID + DAR_77 + Noon_Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Noon_Day ~ DAR_77, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Noon_Day", "Day", "Night")) %>%
  mutate(DAR_77_Ratio = (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  arrange(patient_ID, Noon_Day)

if ("DAR" %in% colnames(actigraphy)) {
  #' DAR == 12-hour period based on individual's self-reported chronotype
  results$DAR <- aggregate(Activity ~ patient_ID + DAR + Noon_Day, actigraphy, mean) %>%
    reshape2::dcast(., formula = patient_ID + Noon_Day ~ DAR, value.var = "Activity") %>%
    set_colnames(c("patient_ID", "Noon_Day", "Day", "Night")) %>%
    mutate(DAR = (Day/(Day + Night))) %>%
    filter(complete.cases(.)) %>%
    arrange(patient_ID, Noon_Day)
}

## ------ Total Activity and Light Exposure ------

results$total_activity_light <- aggregate(cbind(Activity, Light) ~ patient_ID + Noon_Day,
                                  actigraphy, sum)

## ------ Patient Catalog ------

results$patient_catalog <- xtabs(~ patient_ID + Noon_Day + Group + Date, actigraphy) %>%
  as.data.frame %>%
  filter(Freq > 0) %>%
  mutate(Group = factor(Group, levels = c("C", "CCT"),
                        labels = c("Ctrl", "CCT"))) %>%
  select(-Freq) %>%
  arrange(patient_ID, Noon_Day) %>%
  split(., .[, "patient_ID"]) %>%
  lapply(., function(i) i[seq(1, nrow(i), 2), ]) %>%
  do.call("rbind", .) %>%
  arrange(patient_ID, Noon_Day)

## ------ Major Metrics by Patient by Day ------

major_results <- results$sleep_metrics %>%
          select(-Interval, -Sleep, -Wake) 

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

saveRDS(results, ".\\Rmd\\Data\\actigraphy_results.rds")
