## Script: Tests -- IV & IS
## Project: mblt-actigraphy
## Author: Ryan Opel -- @opelr
## Date: 2017-11-16
## Version: 0.3.0

library(magrittr)
library(tidyverse)
library(reshape2)

source("Scripts/Functions.R")

## ------ Create Main Data Frame ------
#' This object serves as fake test data

fake_actig <- data.frame(patient_ID = rep(1:3, each = 14*720),
                         Noon_Day = rep(rep(1:14, each = 720), 3),
                         Epoch = rep(seq(1, 720 * 14, 1), 3),
                         Hour = rep(rep(1:24, each = 30), 3 * 14)) %>%
  mutate(Sine_Normal = sin(Epoch / round(720 / (2 * pi))),
         Sign_Sine_Normal = sign(Sine_Normal),
         Gaussian_Noise = rnorm(Epoch),
         Sine_Fuzz_p1 = Sine_Normal + rnorm(Epoch, 0, 0.1),
         Sine_Fuzz_p5 = Sine_Normal + rnorm(Epoch, 0, 0.5),
         Sine_Fuzz_p1_Zero = pmax(Sine_Fuzz_p1, 0),
         Sine_Fuzz_p5_Zero = pmax(Sine_Fuzz_p5, 0),
         Sine_Phase_p875 = sin(Epoch / round(720 / (1.75 * pi))),
         Sine_Phase_p750 = sin(Epoch / round(720 / (1.50 * pi))),
         Sine_Fuzz_p1_Phase_p875 = Sine_Phase_p875 + rnorm(Epoch, 0, 0.1),
         Sine_Fuzz_p5_Phase_p875 = Sine_Phase_p875 + rnorm(Epoch, 0, 0.5),
         Sine_Fuzz_p1_Phase_p875_Zero = pmax(Sine_Fuzz_p1_Phase_p875, 0),
         Sine_Fuzz_p5_Phase_p875_Zero = pmax(Sine_Fuzz_p5_Phase_p875, 0),
         Square_Wave = sign(Sine_Normal),
         Square_Fuzz_p1 = Square_Wave + rnorm(Epoch, 0, 0.1),
         Square_Fuzz_p5 = Square_Wave + rnorm(Epoch, 0, 0.5),
         Square_Phase_p875 = sign(sin(Epoch / round(720 / (1.75 * pi)))),
         Sqaure_Fuzz_p1_Phase_p875 = Square_Phase_p875 + rnorm(Epoch, 0, 0.1),
         Sqaure_Fuzz_p5_Phase_p875 = Square_Phase_p875 + rnorm(Epoch, 0, 0.5))

## ------ Combintations ------

IV_IS_combs <- expand.grid(c("Sine_Normal", "Sign_Sine_Normal", "Gaussian_Noise",
                             "Sine_Fuzz_p1", "Sine_Fuzz_p5",
                             "Sine_Fuzz_p1_Zero", "Sine_Fuzz_p5_Zero",
                             "Sine_Phase_p875", "Sine_Phase_p750",
                             "Sine_Fuzz_p1_Phase_p875", "Sine_Fuzz_p5_Phase_p875",
                             "Sine_Fuzz_p1_Phase_p875_Zero", "Sine_Fuzz_p5_Phase_p875_Zero",
                             "Square_Wave", "Square_Fuzz_p1", "Square_Fuzz_p5",
                             "Square_Phase_p875", "Sqaure_Fuzz_p1_Phase_p875",
                             "Sqaure_Fuzz_p5_Phase_p875"),
                           c("moving", "expanding"),
                           c(3, 7)) %>% 
  rename(var = Var1, window = Var2, size = Var3) %>%
  arrange(var, window, size)

## ------ Test ------

IV_IS_test <- list()

for (j in 1L:nrow(IV_IS_combs)) {
  ## Get function variable strings from data.frame
  var <- IV_IS_combs$var[j] %>% as.character
  window <- IV_IS_combs$window[j] %>% as.character
  size <- IV_IS_combs$size[j]
  
  ## Concatenate strings to generate title
  IV_title <- paste(var, "IV", window, size, sep = "_")
  IS_title <- paste(var, "IS", window, size, sep = "_")
  
  ## Loop to calculate IV and IS
  IV <- lapply(unique(fake_actig$patient_ID), function(ii) {
    calculate_IV(var, ii, window, size, fake_actig)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IV", IV_title, names(.))))
  
  IS <- lapply(unique(fake_actig$patient_ID), function(ii) {
    calculate_IS(var, ii, window, size, fake_actig)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IS", IS_title, names(.))))
  
  ## Append IV & IS for specific variable to list
  IV_IS_test[[IV_title]] <- IV
  IV_IS_test[[IS_title]] <- IS
}

## ------ Merge IV & IS Results Objects ------

IV_IS_test_results <- IV_IS_test[1][[1]]

for (i in 2L:length(IV_IS_test)) {
  IV_IS_test_results <- merge(IV_IS_test_results, IV_IS_test[i][[1]],
        by = c("patient_ID", "Noon_Day"))
}

IV_IS_test_results %<>%
  arrange(patient_ID, Noon_Day) %>%
  mutate(patient_ID = factor(patient_ID))

rm(IV_IS_test)

## ------ Plot fake_actig Variables ------

pdf("Tests/Sine_Noise_IV_IS_Test_Plots.pdf", width = 15, height = 10)

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Normal, main = "Standard Sine Wave"))
text(50, -0.5, "IV = 0\nIS = 1", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Gaussian_Noise, main = "Gaussian Noise"))
text(50, -2.5, "IV = 1\nIS = 0", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p1, main = "Sine Wave + Gaussian Noise (SD = 0.1)"))
text(50, -0.5, "IV = 0.05\nIS = 0.95", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p5, main = "Sine Wave + Gaussian Noise (SD = 0.5)"))
text(50, -1.5, "IV = 0.65\nIS = 0.65", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p1_Zero, main = "Sine Wave + Gaussian Noise (SD = 0.1), Min Zero"))
text(600, 1, "IV = 0.05\nIS = 0.95", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p5_Zero, main = "Sine Wave + Gaussian Noise (SD = 0.5), Min Zero"))
text(600, 1.5, "IV = 0.65\nIS = 0.65", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Phase_p875, main = "Sine Wave at 87.5% Phase"))
text(50, -0.5, "IV = 0\nIS = 0.60", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p1_Phase_p875, main = "Sine Wave at 87.5% Phase + Gaussian Noise (SD = 0.1)"))
text(50, -0.5, "IV = 0.70\nIS = 0.40", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p5_Phase_p875, main = "Sine Wave at 87.5% Phase + Gaussian Noise (SD = 0.5)"))
text(50, -1.5, "IV = 0.70\nIS = 0.40", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p1_Phase_p875_Zero,
          main = "Sine Wave at 87.5% Phase + Gaussian Noise (SD = 0.1), Min Zero"))
text(600, 1, "IV = 0.07\nIS = 0.50", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sine_Fuzz_p5_Phase_p875_Zero,
          main = "Sine Wave at 87.5% Phase + Gaussian Noise (SD = 0.5), Min Zero"))
text(600, 1.5, "IV = 0.90\nIS = 0.35", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Square_Wave,
          main = "Square Wave"))
text(50, -0.5, "IV = 0.01\nIS = 0.95", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Square_Fuzz_p1,
          main = "Square Wave + Gaussian Noise (SD = 0.1)"))
text(50, -0.5, "IV = 0.05\nIS = 0.90", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3), ##
     plot(Epoch, Square_Fuzz_p5,
          main = "Square Wave + Gaussian Noise (SD = 0.5)"))
text(50, -1.5, "IV = 0.40\nIS = 0.75", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Square_Phase_p875,
          main = "Square Wave at 87.5% Phase"))
text(50, -0.5, "IV = 0\nIS = 0.45", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sqaure_Fuzz_p1_Phase_p875,
          main = "Square Wave at 87.5% Phase + Gaussian Noise (SD = 0.1)"))
text(50, -0.5, "IV = 0.03\nIS = 0.45", adj = 0)
abline(v = seq(720, 720 * 2, 720))

with(filter(fake_actig, patient_ID == 1, Noon_Day %in% 1:3),
     plot(Epoch, Sqaure_Fuzz_p5_Phase_p875,
          main = "Square Wave at 87.5% Phase + Gaussian Noise (SD = 0.5)")) ##
text(50, -0.5, "IV = 0.40\nIS = 0.40", adj = 0)
abline(v = seq(720, 720 * 2, 720))

dev.off()

## ------ Plot IV and IS ------

longitudinal_plot <- function(var) {
  p <- ggplot(IV_IS_test_results, aes_string("Noon_Day", var, color = "patient_ID")) +
    geom_point(size = 2) +
    geom_smooth() +
    ylim(0, 2)
  
  plot(p)
}

test_IV_IS_names <- data.frame(beep = rep(c("IV", "IS"), each = nrow(IV_IS_combs)),
                               rbind(IV_IS_combs, IV_IS_combs)) %>%
  mutate(nam = paste(var, beep, window, size, sep = "_")) %>%
  arrange(var, beep, window, size) %>%
  select(nam) %>% 
  unlist

## Plotting Loop
pdf("Tests/test_IV_IS_plots.pdf", width = 15, height = 9)

for (j in 1:length(test_IV_IS_names)) {
  longitudinal_plot(test_IV_IS_names[j])
}

dev.off()
