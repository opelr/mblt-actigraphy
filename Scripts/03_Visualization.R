## Script: Visualization of Actiwatch Data for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-08-01
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(lubridate)
library(reshape2)

actigraphy_static <- readRDS(".\\Rmd\\Data\\actigraphy_data.rds") 
actigraphy <- readRDS(".\\Rmd\\Data\\actigraphy_filtered.rds") 
results <- readRDS(".\\Rmd\\Data\\results.rds") 

## ------ Visualizing a single patient ------

d1 <- dplyr::filter(actigraphy_static, patient_ID == "Peyton_Teutsch") %>%
  select(., Time, Light, Activity, Day, DateAbbr, Sleep_Thresh_Smooth) %>%
  mutate(Activity_Scale = Activity * (max(Light) / max(Activity)),
         Log_Light = log10(Light))

ggplot(d1, mapping = aes(x = Time)) +
  geom_line(aes(y = Light), color = "Orange") +
  geom_line(aes(y = Activity_Scale), color = "Black") +
  facet_grid(Day + DateAbbr ~ .) +
  scale_x_datetime(date_labels = "%I %p", date_minor_breaks = "1 hour") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * (max(d1$Activity) / max(d1$Light)),
                                         name = "Activity")) +
  ylab("Light")

## ------ Light vs Activity ------

ggplot(d1, aes(Log_Light, Activity)) + 
  geom_jitter() +
  geom_smooth(method = "lm")

## ------ Percent Time Off-Wrist ------

ggplot(results$percent_off_wrist, aes(Day_number, Percent_Off,
                                      group = patient_ID, color = patient_ID)) + 
  geom_line(size = 1) +
  facet_grid(patient_ID ~ .)

## ------ Daytime activity ratio (DAR) ------

ggplot(results$DAR, aes(x = patient_ID, y = DAR)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25) +
  labs(y = "Daytime Activity Ratio (6:30 AM - 11 PM)", x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## ------ Daytime activity ratio (DAR) - 7AM - 7PM ------

ggplot(results$DAR_77, aes(x = patient_ID, y = DAR)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25) +
  labs(y = "Daytime Activity Ratio (7 AM - 7 PM)", x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## ------ Nighttime sleep percentage ------
# `Line` is being used as a dummy variable

ggplot(results$NSD, aes(x = patient_ID, y = Percent_Sleep)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.25) + 
  labs(y = "Percent Time Sleeping (Sunset - Sunrise)")

## ------ Percent Sleep/Social Jetlag ------

ggplot(results$social_jetlag, aes(patient_ID, Percent_Sleep)) + 
  geom_bar(aes(fill = Weekend), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 33)

## ------ Average Bed/Wake Times ------

ggplot(results$avg_bedtime, aes(DateTime, patient_ID, color = Status)) +
  geom_point() +
  scale_x_datetime(date_breaks = "2 hours",
                   date_labels = "%H:%M") +
  geom_errorbarh(aes(xmax = upper_Date, xmin = lower_Date), height = 0.1)

## ------ Light Adherence ------

ggplot(cbind(aggregate(ABL ~ patient_ID, results$light_adherence, mean),
             aggregate(Bouts ~ patient_ID, results$light_adherence, sum)[2]),
       aes(Bouts, ABL)) +
  geom_jitter(aes(color = patient_ID), width = 0.05, height = 0.05) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 4))

## ------ Average Light Exposure -------

ggplot(results$median_light_exposure, aes(Hour, Light, group= patient_ID,
                                       color = patient_ID)) + 
  scale_x_continuous(breaks = seq(0, 24, 4)) + 
  geom_line(stat = "identity", size = 1)

## ------ Sleep Duration vs. Light/Activity per day ------

p <- ggplot(results$sleep_dur_light_activity_per_day, aes(Activity, Hours)) +
  geom_point(aes(color = patient_ID)) +
  geom_smooth(method = "lm") +
  labs(y = "Sleep (hours)")

plotly::ggplotly(p)

## ------ Radial Plot ------

pat <- "Ryan_Opel"

radial <-  dplyr::filter(results$radial_plots, patient_ID == pat | patient_ID == "mean" )

ggplot(radial, aes(Time, Value, color = Metric)) +
  geom_line(data = radial[radial$Metric == "Light" & radial$patient_ID == pat, ],
            color = "darkorange", size = 2) +
  geom_line(data = radial[radial$Metric == "Activity" & radial$patient_ID == pat, ],
            color = "black", size = 2) +
  geom_line(data = radial[radial$Metric == "Light" & radial$patient_ID == "mean", ],
            color = "darkorange", size = 2, alpha=0.2) +
  geom_line(data = radial[radial$Metric == "Activity" & radial$patient_ID == "mean", ],
            color = "black", size = 2, alpha=0.2) +
  scale_x_datetime(date_labels = "%I %p", date_breaks = "3 hours", date_minor_breaks = "1 hour") +
  # scale_y_continuous(limits = c(0, 10 ^ ceiling(log10(max(radial$Value))))) +
  coord_polar() +
  theme_minimal() +
  labs(y = "", x = "")

## ------ Activity Consolidation ------

ggplot(results$activity_consolidation, aes(Bouts, ABL)) +
  geom_point(aes(color = patient_ID)) +
  geom_smooth(method = "lm")

consoli_days <- merge(aggregate(Date ~ patient_ID, results$patient_dates, length) %>%
                             rename(Total_Days = Date),
                      aggregate(Date ~ patient_ID, results$activity_consolidation,
                                length) %>%
                        rename(Consol_days = Date),
                           by = "patient_ID", all.x = T) %>%
  mutate(Percent = 100 * (Consol_days / Total_Days)) %>%
  opelr::replace_na(., 0)

ggplot(consoli_days, aes(Total_Days, Percent)) +
  geom_point(aes(color = patient_ID)) +
  geom_smooth(method = "lm")


ggplot(consoli_days, aes(patient_ID, Percent)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## ------ WASO ------

WASO_agg <- merge(aggregate(WASO ~ patient_ID,  WASO, mean),
                  aggregate(WASO ~ patient_ID,  WASO, opelr::sem) %>%
                    rename(SEM = WASO),
                  by = "patient_ID") %>%
  mutate(ymax = WASO + SEM,
         ymin = WASO - SEM)

ggplot(WASO_agg, aes(patient_ID, WASO)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.5) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "WASO (%)", x = "")

