## Script: Sample Size Estimates/Power Analysis for H4085
## Project: actigraphy -- H4085
## Author: Ryan Opel -- @opelr
## Date: 2017-11-17
## Version: 0.1.0

library(magrittr)
library(tidyverse)
library(pwr2)
library(reshape2)

source("Scripts/Functions.R")

actigraphy <- readRDS(".\\Rmd\\Data\\actigraphy_filtered.rds")
results <- readRDS(".\\Rmd\\Data\\results.rds")

## ------ Sample Size Estimates ------

pwr.2way(a=3, b=3, alpha=0.05, f.A=0.8, f.B=0.4, power=0.8)
