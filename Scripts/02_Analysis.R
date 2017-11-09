## Script: Analysis of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-12
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(zoo)

actigraphy <- readRDS(".\\Rmd\\Data\\actigraphy_filtered.rds") 

results <- list()

## ------ Basic Sleep Metrics ------

#' Creating a parent dataframe used to calculate TST, WASO, and SE

results$sleep_metrics <- xtabs(~ patient_ID + Sleep_Bouts + Sleep_Acti_Smooth + Noon_Day,
                               data = actigraphy) %>%
  as.data.frame %>%
  dcast(., patient_ID + Noon_Day + Sleep_Bouts ~ Sleep_Acti_Smooth,
        value.var = "Freq") %>%
  filter(Sleep_Bouts == "Sleep_Bout", (Sleep + Wake) > 0) %>%
  mutate(Sleep_Minutes = Sleep * 2,
         WASO_Minutes = Wake * 2,
         SE = 100 * Sleep / (Wake + Sleep))

## ------ Sleep Latency ------



## ------ Activity Rhythm Stability ------

### Interdaily Stability
calculate_IS <- function(var, patient, window, size) {
  #' Calculate Interdaily Stability for a non-parametric actigraphy variable
  #' The extent to which the profiles of individual days resemble each other
  #' 
  #' Args:
  #'   var (str): Column name present in your 'actigraphy' dataframe
  #'   patient (int): Patient ID number
  #'   window (str): Options are "moving" or "expanding"
  #'   size (int): If window is set to "moving", this is the number of days
  #'               considered. If "expanding", specifies the buffer before
  #'               initial IS is calcualted.
  #' 
  #' Returns:
  #'   Data frame containing IS for every unique patient in your actigraphy data frame.
  #'   Values range from 0 (Guassian noise) to 1 (perfect IS).
  #'   
  #' Reference:
  #'   Van Someren, Chronobiology International, 1999
  
  ## Helper fn
  fp <- function(x, y) {
    # Shorthand notation for 'as.formula(paste(.))'
    return(as.formula(paste(x, "~", y)))
  }
  
  pat_unique_days <- unique(actigraphy$Noon_Day[actigraphy$patient_ID == patient])
  pat_seq_days <- seq(pat_unique_days[1], pat_unique_days[length(pat_unique_days)], 1)
  
  IS <- lapply(size:length(pat_unique_days), function (ii) {
    
    if (window == "expanding") {
      dat <- actigraphy[actigraphy$patient_ID == patient &
                        actigraphy$Noon_Day %in% pat_unique_days[1]:pat_unique_days[ii], ]
    } else if (window == "moving") {
      dat <- actigraphy[actigraphy$patient_ID == patient &
                          actigraphy$Noon_Day %in% pat_unique_days[ii - size + 1]:pat_unique_days[ii], ]
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
    total_mean <- merge(actigraphy[, c("patient_ID", "Epoch", var)] %>%
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

### Intradaily Variability
calculate_IV <- function(var, patient, window, size) {
  #' Calculate Intradaily Variability for a non-parametric actigraphy variable
  #' Quantifies how fragmented the rhythm is relative to the overall variance
  #' 
  #' Args:
  #'   var (str): Column name present in your 'actigraphy' dataframe
  #'   patient (int): Patient ID number
  #'   window (str): Options are "moving" or "expanding"
  #'   size (int): If window is set to "moving", this is the number of days
  #'               considered. If "expanding", specifies the buffer before
  #'               initial IS is calcualted.
  #' 
  #' Returns:
  #'   Data frame containing IV for every unique patient in your actigraphy data frame.
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
  
  pat_unique_days <- unique(actigraphy$Noon_Day[actigraphy$patient_ID == patient])
  pat_seq_days <- seq(pat_unique_days[1], pat_unique_days[length(pat_unique_days)], 1)
  
  IV <- lapply(size:length(pat_unique_days), function (ii) {
    
    if (window == "expanding") {
      dat <- actigraphy[actigraphy$patient_ID == patient &
                          actigraphy$Noon_Day %in% pat_unique_days[1]:pat_unique_days[ii], ]
    } else if (window == "moving") {
      dat <- actigraphy[actigraphy$patient_ID == patient &
                          actigraphy$Noon_Day %in% pat_unique_days[ii - size + 1]:pat_unique_days[ii], ]
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

### Domintant Rest Phase Onset
# DRPO <- function() {
#' Represents the clock time at which the 5-hr period of lowest activity in
#' the average 24-hr pattern started
#   
# }

IV_IS_combs <- expand.grid(c("Activity", "Light"), c("moving", "expanding"),
                           c(3, 7)) %>% 
  rename(var = Var1, window = Var2, size = Var3) %>%
  arrange(var, window, size)

for (j in 1L:nrow(IV_IS_combs)) {
  var <- IV_IS_combs$var[j] %>% as.character
  window <- IV_IS_combs$window[j] %>% as.character
  size <- IV_IS_combs$size[j]
  
  IV_title <- paste(var, "IV", window, size, sep = "_")
  IS_title <- paste(var, "IS", window, size, sep = "_")
  
  IV <- lapply(unique(actigraphy$patient_ID), function(ii) {
    calculate_IV(var, ii, window, size)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IV", IV_title, names(.))))
  
  IS <- lapply(unique(actigraphy$patient_ID), function(ii) {
    calculate_IS(var, ii, window, size)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IS", IS_title, names(.))))
  
  results[[IV_title]] <- IV
  results[[IS_title]] <- IS
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

results$RA <- merge(calc_rest_phase(10, 'max') %>%
                      rename(M10_Activity = Activity_Mean),
                    calc_rest_phase(5, 'min')  %>%
                      rename(L5_Activity = Activity_Mean),
                    by = c("patient_ID", "Noon_Day")) %>%
  mutate(RA = (M10_Activity - L5_Activity) / (M10_Activity + L5_Activity))

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

results$activity_consolidation <- calc_activity_consolidation("Activity_15_2k") %>%
  select(-Values, -Lengths) %>%
  rename(Activity_Consol_Bouts = Bouts,
         Activity_Consol_Minutes = Minutes,
         Activity_Consol_ABL = ABL)

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

results$light_adherence <- calc_light_adherence("Light_15_5k_smooth") %>%
  select(-Values, -Lengths) %>%
  rename(Light_Adherent_Bouts = Bouts,
         Light_Adherent_Minutes = Minutes,
         Light_Adherent_ABL = ABL)

## ------ Major Metrics by Patient by Day ------

results$patient_catalog <- xtabs(~ patient_ID + Noon_Day + Group, actigraphy) %>%
  as.data.frame %>%
  filter(Freq > 0) %>%
  mutate(Group = factor(Group, levels = c("C", "CCT"),
                        labels = c("Ctrl", "CCT"))) %>%
  select(-Freq)

major_results <-results$sleep_metrics %>%
          select(-Sleep_Bouts, -Sleep, -Wake) 

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

results$major_results <- major_results
rm(major_results)

## ----------------------------------------------------------------------------

## ------ Longitudinal Plots ------

# TODO: Chad, this is the function that generates the plot I sent

longitudinal_plot <- function(var) {
  ggplot(results$major_results, aes_string("Noon_Day", var, colour="patient_ID")) +
    geom_point(aes(shape = MBLT_Group), size = 2) +
    geom_smooth(aes(linetype = MBLT_Group)) +
    facet_grid(Group_PTSD ~ .)
}

longitudinal_plot("Activity_IS_moving_3")

## ------ First/Middle/Last ------ 

# Creation

find_closest <- function(x, num) {
  #' Find "Price is Right" style number
  x_sort <- x[order(x)]
  i <- findInterval(num, x_sort)
  return(x_sort[i])
}

## Helper dataframe
results$FML_days <- merge(aggregate(Noon_Day ~ patient_ID, results$major_results, function(x) min(x) + 2) %>%
        rename(Min_Day = Noon_Day),
        aggregate(Noon_Day ~ patient_ID, results$major_results, function(x) find_closest(x, 7)) %>% 
          rename(Middle_Day = Noon_Day),
      by = "patient_ID") %>%
  merge(., aggregate(Noon_Day ~ patient_ID, results$major_results, function(x) find_closest(x, 28) - 1) %>% 
          rename(Max_Day = Noon_Day),
        by = "patient_ID") %>%
  arrange(patient_ID)

## FML DF
results$FML <- split(results$major_results, results$major_results$patient_ID) %>%
  lapply(., function(ii) {
    ID <- unique(ii$patient_ID)
    dayz <- unlist(results$FML_days[results$FML_days$patient_ID == ID, ][-1])
    ii <- ii[ii$Noon_Day %in% dayz, ]
    ii$Time_Point <- c("Initial", "Base_End", "Final")
    return(ii)
}) %>%
  do.call("rbind", .) %>%
  set_rownames(1:nrow(.)) %>%
  filter(., Time_Point != "Base_End")

#  Analysis

FML_vars <- c("Sleep_Minutes", "WASO_Minutes", "SE", "Activity_IS",
              "Activity_IV", "Light_IS", "Light_IV", "RA")

lapply(FML_vars, function(j) {
  form <- as.formula(paste(j, "~ Group * Time_Point"))
  ANOV <- anova(lm(form, results$FML))
  
  sf <- function(num) {
    signif(num, 3)
  }
  
  data.frame(IV = j,
             Group_F = sf(ANOV$`F value`[1]),
             Group_p = sf(ANOV$`Pr(>F)`[1]),
             Time_F = sf(ANOV$`F value`[2]),
             Time_p = sf(ANOV$`Pr(>F)`[2]),
             Int_F = sf(ANOV$`F value`[3]),
             Int_p = sf(ANOV$`Pr(>F)`[3]))
  
}) %>%
  do.call("rbind", .)

# Group Plot

Group_FML_plot <- function(var) {
  i_1 <- select(results$FML, one_of(var), Group, Time_Point)
  
  i_2 <- merge(aggregate(. ~ Group + Time_Point, i_1, mean),
               aggregate(. ~ Group + Time_Point, i_1, sciplot::se) %>%
                 set_colnames(c("Group", "Time_Point", "SEM")),
               by = c("Group", "Time_Point")) %>%
    mutate(ymax = .[, var] + SEM,
           ymin = .[, var] - SEM,
           Time_Point = factor(Time_Point, levels = c("Initial", "Base_End", "Final")))
  
  p <- ggplot(i_2, aes_string("Group", var, fill = "Time_Point")) +
    geom_bar(position = "dodge", stat = "identity", colour='black') + 
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, position=position_dodge(.9))
  
  plot(p)
}

Group_FML_plotFML_plot(FML_vars[1])

# Individual Plot

Indiv_FML_plot <- function(var) {
  i_1 <- select(results$major_results, patient_ID, one_of(var), Group, Noon_Day) %>%
    filter(Group == "CCT")
  
  p <- ggplot(i_1, aes_string("Noon_Day", var, group = "patient_ID")) +
    geom_point(aes(color = Group), size = 3) + 
    geom_line(aes(color = Group), size = 1)
    
  
  plot(p)
}

Indiv_FML_plot(FML_vars[1])

## ------ Time Series Analysis ------

anova(lm(Sleep_Minutes ~ Group * Noon_Day, major_results))


ggplot(major_results, aes(Noon_Day, Sleep_Minutes, group = patient_ID,
                          color = patient_ID)) +
  geom_point() +
  geom_smooth()

## ----------------------------------------------------------------------------

## ------ Median Light Exposure -------

#' TODO: This aggregation should also take into account Baseline vs. Test
#'       timepoints once we have subjects with Negative Ion Generators and
#'       lightboxes with built-in baseline periods.

results$median_light_exposure <-
  cbind(aggregate(Light ~ patient_ID + Hour, actigraphy, median),
        aggregate(Light ~ patient_ID + Hour, actigraphy, opelr::sem)[3] %>%
          rename(SEM = Light)) %>%
  mutate(ymax = Light + SEM,
         ymin = Light - SEM) %>%
  arrange(patient_ID, Hour)


## ------ Percent Time Off-Wrist ------

## TODO: This is no longer accurate, because each "Watch_Off" counts is 3 hours

percent_off_wrist <- aggregate(cbind(Watch_Off, Total_Epochs) ~ patient_ID,
                               off_wrist[opelr::numfac(off_wrist$Day) <= 28, ], sum) %>%
  mutate(Percent_Off = 100 * Watch_Off / Total_Epochs)

results$percent_off_wrist <- percent_off_wrist

## ------ Daytime activity ratio (DAR) ------
# taken from http://www.neurology.org/content/88/3/268.long

### Daytime is defined as 6:30 AM - 11 PM --> can be modified in 01_ETL.R script
DAR_parent <- aggregate(Activity ~ patient_ID + DAR_Period + Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ DAR_Period, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Day_number", "Day", "Night")) %>%
  mutate(DAR = 100 * (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number"))

DAR <- cbind(aggregate(DAR ~ patient_ID, DAR_parent, mean),
                   aggregate(DAR ~ patient_ID, DAR_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                    })[2]) %>%
  set_colnames(c("patient_ID", "DAR", "SEM")) %>%
  mutate(ymax = DAR + SEM,
         ymin = DAR - SEM)

results$DAR <- DAR

### Daytime is defined as 7:00 AM - 7 PM
DAR_77_parent <- aggregate(Activity ~ patient_ID + DAR_77 + Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ DAR_77, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Day_number", "Day", "Night")) %>%
  mutate(DAR = 100 * (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number"))

DAR_77 <- cbind(aggregate(DAR ~ patient_ID, DAR_77_parent, mean),
                   aggregate(DAR ~ patient_ID, DAR_77_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                   })[2]) %>%
  set_colnames(c("patient_ID", "DAR", "SEM")) %>%
  mutate(ymax = DAR + SEM,
         ymin = DAR - SEM)

results$DAR_77 <- DAR_77

## ------ Nighttime sleep duration ------
# `Line` is being used as a dummy variable

NSD_parent <- aggregate(Line ~ patient_ID + Sleep_Thresh_Smooth + SunPeriod + Day,
                        actigraphy, length) %>%
  reshape2::dcast(., formula = patient_ID + SunPeriod + Day ~ Sleep_Thresh_Smooth,
                  value.var = "Line") %>%
  opelr::replace_na(., 0) %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake))) %>%
  filter(SunPeriod == "Night")

NSD_child <- cbind(aggregate(Percent_Sleep ~ patient_ID, NSD_parent, mean),
                   aggregate(Percent_Sleep ~ patient_ID, NSD_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                   })[2]) %>%
  set_colnames(c("patient_ID", "Percent_Sleep", "SEM")) %>%
  mutate(ymax = Percent_Sleep + SEM,
         ymin = Percent_Sleep - SEM)

results$NSD <- NSD_child

## ------ Percent Sleep/Social Jetlag ------

percent_sleep <- aggregate(Sleep_Thresh_Smooth ~ patient_ID + MBLT_Group + Date + Day, actigraphy, table) %>%
  as.data.frame.table(.) %>%
  filter(Var2 == "patient_ID", Freq.Day != 1) %>%
  select(-Var1, -Var2, -Freq.Day) %>%
  set_colnames(c("patient_ID", "MBLT_Group", "Date", "Sleep", "Wake")) %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake)),
         Date = as.POSIXct(strptime(Date, format = "%F")),
         DayOfWeek = weekdays(Date),
         Weekend = factor(ifelse(DayOfWeek %in% c("Saturday", "Sunday"),
                                 "Weekend", "Weekday"))) %>%
  arrange(patient_ID, Date)

# results$percent_sleep <- percent_sleep

## ------ Bedtime ------

get_time <- function (time) {
  # Converts POSIX to decimal
  # 
  # Args:
  #   time (POSIX)
  # 
  # Returns:
  #   Time as a decimal
  
  if(!any(class(time) == "POSIXt"))
    stop("time must be POSIX")
  
  tim <- as.numeric(time - trunc(time, "days"))
  return(tim)
}

find_bedtime <- function(noon_day, ID, sleep_column, status = "bed", data) {
  # Bedtime/Waketime
  # 
  # Args:
  #   noon_day (int): Which day (noon-noon) we're investiagting
  #   ID (str): Patient name
  #   sleep_column (str): Which column are we running the function on
  #   status (str): "bed" or "wake" for bed- and wake-times
  #   data (df): Data frame to iterate upon
  # 
  # Returns:
  #   Bed- or Wake-time as POSIXct
  
  if (! toupper(status) %in% c("BED", "WAKE"))
    stop("'status' must be either 'bed' or 'wake'")
  
  df <- data %>%
    dplyr::filter(patient_ID == ID, Noon_Day == noon_day)
  
  sleep <- data.frame(unclass(rle(as.character(df[, sleep_column])))) %>%
    mutate(
      Index = cumsum(c(1, lengths[-length(lengths)])),
      Time = df[, "DateTime"][Index]
    )
  
  max_sleep <- which(sleep$lengths == max(sleep$lengths[sleep$values == "Sleep"])) %>%
    .[. %in% which(sleep$values == "Sleep")]
  
  bed_time <- sleep$Time[max_sleep]
  wake_time <- sleep$Time[max_sleep + 1]
  
  if (toupper(status) == "BED") {
    return(bed_time)
  } else {
    return(wake_time)
  }
}

get_bedtime_loop <- function(data, status) {
  # lapply 'find_bedtime' function for 'bedtime' data.frame
  # 
  # Args:
  #   data (data.frame)
  #   status (str): Bed or Wake, same as find_bedtime
  # 
  # Returns:
  #   Vector of POSIX times
  
  lapply(rownames(data), function(ii) {
    ii <- as.integer(ii)
    
    find_bedtime(data[, "Noon_Day"][ii], ID = data[, "patient_ID"][ii],
                 sleep_column = "Sleep_Thresh_Smooth", status = status,
                 data = actigraphy)
  }) %>%
    do.call("c", .)
}
 
bedtime <- as.data.frame(table(actigraphy$patient_ID, actigraphy$Noon_Day)) %>%
  set_colnames(c("patient_ID", "Noon_Day", "Freq")) %>%
  dplyr::filter(Freq == 720) %>% 
  arrange(patient_ID, Noon_Day) %>%
  mutate(Bedtime = get_bedtime_loop(., 'bed'),
         Waketime = get_bedtime_loop(., 'wake'),
    Day_of_Week = weekdays(Bedtime),
    Weekend = factor(ifelse(Day_of_Week %in% c("Saturday", "Sunday"),
                            "Weekend", "Weekday")),
    Bedtime_Decimal = get_time(Bedtime),
    Waketime_Decimal = get_time(Waketime),
    Bedtime_NoDate = as.POSIXct(strptime(substr(Bedtime, 12, 20), format = "%T")),
    Waketime_NoDate = as.POSIXct(strptime(substr(Waketime, 12, 20), format = "%T"))) %>%
  dplyr::filter(complete.cases(.))

is_noon_time <- aggregate(DateTime ~ patient_ID + Noon_Day, actigraphy, min) %>%
  rename(Prev_Day = DateTime) %>%
  mutate(IsNoon = substr(Prev_Day, 12, 20) == "12:00:00") %>%
  dplyr::filter(IsNoon == T) %>%
  select(-IsNoon)

bedtime %<>%
  merge(., is_noon_time, by = c("patient_ID", "Noon_Day")) %>%
  mutate(Bedtime_from_Noon = as.integer(difftime(Bedtime, Prev_Day, units = "mins")),
         Waketime_from_Noon = as.integer(difftime(Waketime, Prev_Day, units = "mins")))

## ------ Average Bed/Wake Times ------

# Need to calculate bed and wake time from previous noon -- 
average_bedtime_SEM <- aggregate(cbind(Bedtime_from_Noon, Waketime_from_Noon) ~ patient_ID, 
          bedtime, opelr::sem) %>%
  melt(., id.vars = "patient_ID") %>%
  rename(Status = variable, SEM = value) %>%
  mutate(Status = gsub("time_from_Noon", "", Status))

average_bedtime <- aggregate(cbind(Bedtime_from_Noon, Waketime_from_Noon) ~ patient_ID,
                             bedtime, mean) %>%
  melt(., id.vars = "patient_ID") %>%
  rename(Status = variable, Time = value) %>%
  mutate(Status = gsub("time_from_Noon", "", Status)) %>%
  merge(., average_bedtime_SEM, by = c("patient_ID", "Status")) %>%
  mutate(upper = Time + SEM,
         lower = Time - SEM,
         Template = as.POSIXct(strptime("2017-01-01 12:00:00", format = "%F %T")),
         DateTime = Template + seconds(Time * 60),
         upper_Date = DateTime + seconds(SEM * 60),
         lower_Date = DateTime - seconds(SEM * 60))

results$avg_bedtime <- average_bedtime

## ------ Sleep Duration vs. Light/Activity per day ------

compare <- merge(percent_sleep, aggregate(Light ~ patient_ID + Date + Day, actigraphy, median),
                 by = c("patient_ID", "Date")) %>%
  merge(., aggregate(Activity ~ patient_ID + Date, actigraphy, median),
        by = c("patient_ID", "Date")) %>%
  mutate(Hours = (Percent_Sleep / 100) * 24)

results$sleep_dur_light_activity_per_day <- compare

## ----- Day Radial Plots -------

radial_individuals <- aggregate(cbind(Light, Activity) ~ patient_ID + Hour, actigraphy, mean) %>%
  melt(., id.vars = c("patient_ID", "Hour")) %>%
  rename(Metric = variable, Value = value) %>%
  arrange(patient_ID, Metric, Hour) %>%
  mutate(Time = as.POSIXct(paste0(sprintf("%02d", Hour), ":00:00"), format = "%T"))

radial_means <- aggregate(Value ~ Metric + Hour, radial_individuals, mean) %>%
  mutate(patient_ID = "mean",
         Time = as.POSIXct(paste0(sprintf("%02d", Hour), ":00:00"), format = "%T")) %>%
  select(one_of(colnames(radial_individuals)))

radial_intermediary <- rbind(radial_individuals, radial_means)

radial_24 <- dplyr::filter(radial_intermediary, Hour == 0) %>%
  mutate(Hour = 24,
         Time = as.POSIXct(paste0(sprintf("%02d", Hour), ":00:00"), format = "%T"))

radial <- rbind(radial_intermediary, radial_24)
  
results$radial_plots <- radial

## ------ Save RDS ------

saveRDS(results, ".\\Rmd\\Data\\results.rds")
