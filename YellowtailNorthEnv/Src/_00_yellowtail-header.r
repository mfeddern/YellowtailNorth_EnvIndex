# Load packages
library(data.table)
library(tidyverse)
library(tidync)
library(ncdf4)
library(readr) # faster writing
library(fst)
library(lubridate) # for working with dates, very handy.
library(here)

# file locations
(home_dir = getwd())
data_raw = "~/GitHub/glorys-download/data-raw/"
results_dir = paste0(home_dir, "/results-yellowtail/")
fig_dir= paste0(home_dir,'/figures-yellowtail/')
data_dir = paste0(home_dir, '/data-yellowtail/') # for processed time series for hake
