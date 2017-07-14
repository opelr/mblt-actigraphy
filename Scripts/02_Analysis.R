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

## ------ Can I determine when the light was on? ------

dat <- actigraphy[actigraphy$patient_ID == "Ryan_Opel", ]

ggplot(dat[dat$Day == 5, ]) +
  geom_line(aes(x = Time, y = Light), color = "Orange") + 
  geom_line(aes(x = Time, y = Activity), color = "Black") +
  scale_x_datetime(date_labels = "%I %p") + 
  scale_y_log10() + 
  facet_grid(Day + DateAbbr ~ .)

## ------ Sleep-wake assessments ------
# taken from http://www.neurology.org/content/88/3/268.long

### Daytime activity ratio (DAR)

DAR_parent <- aggregate(Activity ~ ID + DAR_Period + Day, actigraphy, mean) %>%
  reshape2::dcast(., formula = ID + Day ~ DAR_Period, value.var = "Activity") %>%
  set_colnames(c("ID", "Day_number", "Day", "Night")) %>%
  mutate(DAR = 100 * (Day/(Day + Night))) %>%
  filter(complete.cases(.))

## how do we deal with infinites?

DAR_child <- cbind(aggregate(DAR ~ ID, DAR_parent, mean),
                   aggregate(DAR ~ ID, DAR_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                    })[2]) %>%
  set_colnames(c("ID", "DAR", "SEM")) %>%
  mutate(ymax = DAR + SEM) %>%
  mutate(ymin = DAR - SEM)

ggplot(DAR_child, aes(x = ID, y = DAR)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25)

### Nighttime sleep duration

NSD_parent <- aggregate(Line ~ ID + Sleep_Wake + Day + DAR_Period,
                        actigraphy, length) %>%
  reshape2::dcast(., formula = ID + DAR_Period + Day ~ Sleep_Wake,
                  value.var = "Line") %>%
  mutate(Percent_Sleep = 100 * (Sleep/(Sleep + Wake)))

NSD_child <- cbind(aggregate(Percent_Sleep ~ ID, NSD_parent, mean),
                   aggregate(Percent_Sleep ~ ID, NSD_parent, function(X) {
                     opelr::sem(X, na.rm = T)
                   })[2]) %>%
  set_colnames(c("ID", "Percent_Sleep", "SEM")) %>%
  mutate(ymax = Percent_Sleep + SEM) %>%
  mutate(ymin = Percent_Sleep - SEM)

ggplot(NSD_child, aes(x = ID, y = Percent_Sleep)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25)

## ------ Bedtime ------

sleep <- data.frame(do.call("rbind", lapply(unique(actigraphy$ID), function (ii) {
  spec <- as.data.frame.table(100 * prop.table(table(actigraphy$Date[actigraphy$ID == ii],
                                                     actigraphy$Sleep_Wake[actigraphy$ID == ii]), 1))
  sleep <- cbind(spec, factor(ii))
  return(sleep)
}))) %>% 
  set_colnames(c('Date', 'Stage', 'Sleep', 'ID')) %>%
  filter(complete.cases(.)) %>%
  mutate(Weekend = factor(chron::is.weekend(Date), levels = c(T,F),
                          labels = c("Weekend", "Weekday")))
  
## ------ KNN Sleep Stage Prediction ------
