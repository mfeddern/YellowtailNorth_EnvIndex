library(tidyverse)
library(ggplot2)
library(readr) # faster writing
library(data.table) # faster writing of large files
library(lubridate)
library(mgcv)
library(MuMIn)
library(corrplot)

source("_00_yellowtail-header.r")
source("Src/Functions-for-envir-index.r")
# set directories ##############################################################

# bring in data ################################################################
# combined fish and environmental drivers file
df = data.frame(read.csv(paste0(data_dir,"DATA_Combined_glorys_yellowtail.csv"), header = T))

# get predictors to create model formula #######################################
envir_data = df %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('year','sd','Y_rec','ZOOpjuv','ZOOben')))
head(envir_data)
dim(envir_data)
data_years = 1994:2014

envir_data <- envir_data[complete.cases(envir_data), ]
M = data.frame(cor(envir_data))
M <- tibble::rownames_to_column(M, "xvar")
M<-M%>%pivot_longer(!xvar, names_to = 'yvar', values_to = 'corr')
uncorr<-M%>%filter(corr<0.4&corr> -0.4)

unique(uncorr$xvar)
corrPlot<-corrplot.mixed(M)




pdf(file = "figures-yellowtail/CorrPlot.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12) # The height of the plot in inches
corrplot.mixed(M)

dev.off()
