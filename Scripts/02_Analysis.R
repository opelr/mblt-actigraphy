## Script: Analysis of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-07-12
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(lubridate)

acti_files <- readRDS(".\\Data\\actigraphy_header.rds")
actigraphy <- readRDS(".\\Data\\actigraphy_data.rds") 

## ------ Visualizing a single patient ------

dat <- actigraphy[actigraphy$patient_ID == "Jon_Elliott" & actigraphy$Day == "11", ]

dat <- actigraphy[actigraphy$patient_ID == "Carolyn_Jones",]

ggplot(dat) +
  # geom_line(aes(x = Time, y = Light), color = "Orange") + 
  geom_line(aes(x = Time, y = Activity), color = "Black") +
  scale_x_datetime(date_labels = "%I %p") + 
  # scale_y_log10() + 
  facet_grid(Day + DateAbbr ~ .)


## ------ Off Wrist Filtering ------

off_wrist <- as.data.frame.table(xtabs(~patient_ID + Off_Wrist + Day, data = actigraphy)) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ Off_Wrist, value.var = "Freq") %>%
  rename(Watch_On = `FALSE`, Watch_Off = `TRUE`, Day_number = Day) %>%
  mutate(Percent_Off = 100 * Watch_Off/(Watch_On + Watch_Off),
         Total_Epochs = Watch_On + Watch_Off) %>%
  dplyr::filter(complete.cases(.))

actigraphy %<>% merge(., off_wrist[,c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")] %>% 
               rename(Day = Day_number), by = c("patient_ID", "Day")) %>%
  dplyr::filter(Percent_Off < 10, Total_Epochs >= (720 / 2))

## ------ Percent Time Off-Wrist ------

percent_off_wrist <- aggregate(cbind(Watch_Off, Total_Epochs) ~ patient_ID,
                               off_wrist, sum) %>%
  mutate(Percent_Off = 100 * Watch_Off / Total_Epochs)

ggplot(off_wrist, aes(Day_number, Percent_Off, group = patient_ID, color = patient_ID)) + 
  geom_line(size = 1) +
  facet_grid(patient_ID ~ .)

## ------ Daytime activity ratio (DAR) ------
# taken from http://www.neurology.org/content/88/3/268.long
# Daytime is defined as 7AM - 10 PM --> can be modified in 01_ETL.R script

DAR_parent <- aggregate(Activity ~ patient_ID + DAR_Period + Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = patient_ID + Day ~ DAR_Period, value.var = "Activity") %>%
  set_colnames(c("patient_ID", "Day_number", "Day", "Night")) %>%
  mutate(DAR = 100 * (Day/(Day + Night))) %>%
  filter(complete.cases(.)) %>%
  # Use `off_wrist` to verify that watch is on and recording for most of the day
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number")) #%>%
  # dplyr::filter(Percent_Off < 10, Total_Epochs >= (720 / 2))

DAR_child <- cbind(aggregate(DAR ~ patient_ID, DAR_parent, mean),
                   aggregate(DAR ~ patient_ID, DAR_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                    })[2]) %>%
  set_colnames(c("patient_ID", "DAR", "SEM")) %>%
  mutate(ymax = DAR + SEM,
         ymin = DAR - SEM)

ggplot(DAR_child, aes(x = patient_ID, y = DAR)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25)

## ------ Nighttime sleep duration ------

# Need to decide which sleep determinant is the best in previous script
# Again, daytime is defined as 7AM - 10 PM --> can be modified in 01_ETL.R script

# `Line` is being used as a dummy variable
NSD_parent <- aggregate(Line ~ patient_ID + Sleep_Wake + Day + DAR_Period,
                        actigraphy, length) %>%
  reshape2::dcast(., formula = patient_ID + DAR_Period + Day ~ Sleep_Wake,
                  value.var = "Line") %>%
  rename(Day_number = Day) %>%
  merge(., off_wrist[, c("patient_ID", "Day_number", "Percent_Off", "Total_Epochs")],
        by = c("patient_ID", "Day_number")) %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake)))

NSD_child <- cbind(aggregate(Percent_Sleep ~ patient_ID, NSD_parent, mean),
                   aggregate(Percent_Sleep ~ patient_ID, NSD_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                   })[2]) %>%
  set_colnames(c("patient_ID", "Percent_Sleep", "SEM")) %>%
  mutate(ymax = Percent_Sleep + SEM,
         ymin = Percent_Sleep - SEM)

ggplot(NSD_child, aes(x = patient_ID, y = Percent_Sleep)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25)

## ------ Percent Sleep/Social Jetlag ------

percent_sleep <- lapply(unique(actigraphy$patient_ID), function (ii) {
  df <- actigraphy[actigraphy$patient_ID == ii, ]
  spec <- as.data.frame.table(100 * prop.table(table(df$Date, df$Sleep_Wake), 1))
  sleep <- cbind(spec, factor(ii))
  
  return(sleep)
}) %>% 
  do.call("rbind", .) %>%
  data.frame(.) %>%
  set_colnames(c('Date', 'Stage', 'Sleep', 'patient_ID')) %>%
  filter(complete.cases(.)) %>%
  mutate(Date = as.POSIXct(strptime(Date, format = "%F")),
         Weekend = weekdays(Date),
         Weekend = factor(ifelse(Weekend %in% c("Saturday", "Sunday"),
                                 "Weekend", "Weekday")))

percent_sleep_agg <- merge(aggregate(Sleep ~ patient_ID + Stage + Weekend, percent_sleep, mean),
                           aggregate(Sleep ~ patient_ID + Stage + Weekend, percent_sleep, opelr::sem) %>%
                             rename(SEM = Sleep),
                           by = c("patient_ID", "Stage", "Weekend")) %>%
  mutate(ymax = Sleep + SEM,
         ymin = Sleep - SEM) %>%
  dplyr::filter(Stage == "Sleep")

ggplot(percent_sleep, aes(patient_ID, Sleep)) + 
  geom_bar(aes(fill = Weekend), stat = "identity", position = "dodge")

## ------ Bedtime ------

## Converts POSIX to decimal
getTime <- function (x) {
  if(!any(class(x) == "POSIXt")) stop("x must be POSIX")
  tim <- as.numeric(x - trunc(x, "days"))
  return(tim)
}

find_bedtime <- function(day, ID, data, sleep_column) {
  df <- data %>%
    dplyr::filter(patient_ID == ID & (Day == day & AM_PM == "PM") | (Day == day + 1 & AM_PM == "AM"))
  
  sleep <- data.frame(unclass(rle(as.character(df[, sleep_column])))) %>%
    mutate(
      Index = cumsum(c(1, lengths[-length(lengths)])),
      Time = df[, "DateTime"][Index]
    )
}

bedtime <- do.call("rbind", sapply(unique(watch$Day), function (ii) {
  if (nrow(watch[watch$Day == ii & watch$Noon == "Afternoon", ]) > 1) {
    sl <- watch[watch$Day == ii & watch$Noon == "Afternoon", ]

    x <- data.frame(unclass(rle(as.integer(sl$sleep_agreement))))
    x$index <- cumsum(c(1, x$lengths[-length(x$lengths)]))
    x$values <- levels(sl$sleep_agreement)[x$values]
    x$Time <- sl$DateTime[x$index]
    # x <- x[complete.cases(x),]
    x$Hour <- as.integer(substr(x$Time, 12, 13))

    if(any(is.na(x$values))) {
      x <- x[!is.na(x$values), ]
    }

    maxBedtime <- x$Time[x$lengths == max(x$lengths[x$values == "Sleep"])]

    xx <- x[x$values == "Sleep" & x$Hour > 19, ]
    probBedtime <- xx$Time[which(xx$lengths >= 10)][1]

    if(length(maxBedtime) > 1) {maxBedtime <- maxBedtime[length(maxBedtime)]}

    return(data.frame(Day = ii, maxBedtime = maxBedtime, probBedtime = probBedtime,
               ConsecEpochs = max(x$lengths[x$values == "Sleep"])))
  }

}))

bedtime <- do.call("rbind", sapply(unique(actigraphy$Day), function(day, patient_ID) {
  
  if (nrow(actigraphy[actigraphy$Day == day & actigraphy$Noon == "Afternoon" & actigraphy$patient_ID == patient_ID, ]) > 1) {
    sl <-  actigraphy[actigraphy$Day == day & actigraphy$Noon == "Afternoon" & actigraphy$patient_ID == patient_ID, ]
    
    x <- data.frame(unclass(rle(as.integer(sl$sleep_agreement))))
    
    if(sum(is.na(x$values)) > (length(x$values) / 3)) {
      
  } else {
    
    x$index <- cumsum(c(1, x$lengths[-length(x$lengths)]))
    x$values <- levels(sl$sleep_agreement)[x$values]
    x$Time <- sl$DateTime[x$index]
    # x <- x[complete.cases(x),]
    x$Hour <- as.integer(substr(x$Time, 12, 13))
    
    if(any(is.na(x$values))) {
      x <- x[!is.na(x$values), ]
    }
    
    maxBedtime <- x$Time[x$lengths == max(x$lengths[x$values == "Sleep"])]
    
    xx <- x[x$values == "Sleep" & x$Hour > 19, ]
    probBedtime <- xx$Time[which(xx$lengths >= 10)][1]
    
    if(length(maxBedtime) > 1) {maxBedtime <- maxBedtime[length(maxBedtime)]}
    
    return(data.frame(patient_ID = patient_ID, Day = day, maxBedtime = maxBedtime, probBedtime = probBedtime,
                      ConsecEpochs = max(x$lengths[x$values == "Sleep"])))
  }
    
    
}}, ID= unique(actigraphy$patient_ID)))


# ===---===---===---====----====----====----====---


bedtime$Weekend <- factor(ifelse(isWeekday(bedtime$maxBedtime), "Weekday", "Weekend"))

bedtime$maxNoDay <- strptime(substr(bedtime$maxBedtime, 12, 20), "%T")
bedtime$probNoDay <- strptime(substr(bedtime$probBedtime, 12, 20), "%T")
bedtime$dayOfWeek <- factor(weekdays(bedtime$maxBedtime),
                            levels = rev(c('Sunday', 'Monday', 'Tuesday', 'Wednesday',
                                           'Thursday', 'Friday', 'Saturday')))
bedtime$bedDecimal <- getTime(bedtime$maxNoDay)

## Model: Does sleep differ on weekends?
mod <- anova(lm(bedDecimal ~ Weekend, bedtime))
p <- signif(mod$`Pr(>F)`[1], 3)
F.val <- signif(mod$`F value`[1], 3)

## Plotting
ggplot(bedtime, aes(maxNoDay, dayOfWeek, color = Weekend)) +
  geom_point() +
  scale_x_datetime(date_breaks = "1 hour",
                   date_labels = "%I:%M %p") +
  labs(x = "Hour", title = paste0("Miranda's Bedtime, (p=", p, ", F=", F.val, ")"))


## ------ KNN Sleep Stage Prediction ------
