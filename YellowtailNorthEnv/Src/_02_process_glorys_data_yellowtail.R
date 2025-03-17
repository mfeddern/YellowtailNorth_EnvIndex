# Load packages
library(tidyverse)
library(tidync)
library(ncdf4)
library(readr) # faster writing
library(data.table) # faster writing of large files
library(fst)
library(lubridate)

#source("_00_yellowtail-header.r")

################################################################################
# DAILY TIME SERIES ############################################################
################################################################################

# load data files ####
temp_daily = data.table(read.csv(paste0("data-yellowtail/glorys-daily-temperature.csv")))
dt = temp_daily
# fix date
# get year
dt[,year:=lubridate::year(date),]
# get month
dt[,month:=lubridate::month(date),]
# lag year
dt[, yr_egg:=ifelse(month %in% c(11,12), year+1,year)]
dt[, yr_pre:=ifelse(month %in% c(7:12), year+1,year)]
# reference temp is 3.5
ref_temp = 3.5

DDpre = dt[ month %in% c(7:12,1:3), list(DDpre = mean(dd_90_180, na.rm=TRUE)), by='yr_pre']
  colnames(DDpre)[1] <- 'year'
DDegg = dt[ month %in% c(11,12), list(DDegg = mean(dd_90_180, na.rm=TRUE)), by='yr_egg']
  colnames(DDegg)[1] <- 'year'
DDlarv = dt[ month %in% c(2,3), list(DDlarv = mean(dd_0_90, na.rm=TRUE)), by='year']
DDpjuv = dt[ month %in% c(4:8), list(DDpjuv = mean(dd_30_130, na.rm=TRUE)), by='year']
DDben = dt[ month %in% c(9:12), list(DDben = mean(dd_180_549, na.rm=TRUE)), by='year']
Tcop = dt[ month %in% c(8:10), list(Tcop = mean(temp_90_180, na.rm=TRUE)), by='year']
Tpart = dt[ month %in% c(1:4), list(Tpart = mean(temp_0_180, na.rm=TRUE)), by='year']
#ZOOpjuv = dt[ month %in% c(4:8), list(ZOOpjuv = mean(zoo_30_130, na.rm=TRUE)), by='year']
#ZOOben = dt[ month %in% c(9:12), list(ZOOben = mean(zoo_180_549, na.rm=TRUE)), by='year']

# combine time daily time series
glorys_ts = left_join(DDpre,DDegg)
glorys_ts = left_join(glorys_ts,DDlarv)
glorys_ts = left_join(glorys_ts,DDpjuv)
glorys_ts = left_join(glorys_ts,DDben)
glorys_ts = left_join(glorys_ts,Tcop)
glorys_ts = left_join(glorys_ts,Tpart)

################################################################################
## MONTHLY TIME SERIES #########################################################
################################################################################

# load and combine data
cst_monthly = data.table(read.csv(paste0("data-yellowtail/glorys-monthly-means-cst.csv")))
lst_monthly = data.table(read.csv(paste0("data-yellowtail/glorys-monthly-means-lst.csv")))
mld_monthly = data.table(read.csv(paste0("data-yellowtail/glorys-monthly-means-mld.csv")))

dt = merge.data.table(cst_monthly, lst_monthly, by='date')
dt = merge.data.table(dt, mld_monthly)
dt[,year:=lubridate::year(mdy(date)),]
dt[,month:=lubridate::month(mdy(date)),]

MLDpart = dt[ month %in% c(1:4), list(MLDpart = mean(mld_0_180, na.rm=TRUE)), by='year']
MLDlarv = dt[ month %in% c(2,3), list(MLDlarv = mean(mld_0_90, na.rm=TRUE)), by='year']
MLDpjuv = dt[ month %in% c(4:8), list(MLDpjuv = mean(mld_30_130, na.rm=TRUE)), by='year']

CSTlarv = dt[ month %in% c(2,3), list(CSTlarv = mean(cst_0_90, na.rm=TRUE)), by='year']
CSTpjuv = dt[ month %in% c(4:8), list(CSTpjuv = mean(cst_30_130, na.rm=TRUE)), by='year']
LSTlarv = dt[ month %in% c(2,3), list(LSTlarv = mean(lst_0_90, na.rm=TRUE)), by='year']
LSTpjuv = dt[ month %in% c(4:8), list(LSTpjuv = mean(lst_30_130, na.rm=TRUE)), by='year']

glorys_ts = left_join(glorys_ts,MLDpart)
glorys_ts = left_join(glorys_ts,MLDlarv)
glorys_ts = left_join(glorys_ts,MLDpjuv)
glorys_ts = left_join(glorys_ts,CSTlarv)
glorys_ts = left_join(glorys_ts,CSTpjuv)
glorys_ts = left_join(glorys_ts,LSTlarv)
glorys_ts = left_join(glorys_ts,LSTpjuv)

fwrite(glorys_ts, paste0(data_dir, "glorys-data-annual-yellowtail.csv"))

#### add other predictors ######################################################

# habitat compression index (area of cold water, sort of)
hci1 = data.table(read.csv( paste0("data-yellowtail/ei_hci_rgn1_M.csv" )))
hci1_larv = hci1[ month %in% c(2:3), list(hci1_larv = mean(data, na.rm=TRUE)), by='year']
hci1_pjuv = hci1[ month %in% c(4:8), list(hci1_pjuv = mean(data, na.rm=TRUE)), by='year']

hci2 = data.table(read.csv( paste0("data-yellowtail/ei_hci_rgn2_M.csv")))
hci2_larv = hci1[ month %in% c(2:3), list(hci2_larv = mean(data, na.rm=TRUE)), by='year']
hci2_pjuv = hci1[ month %in% c(4:8), list(hci2_pjuv = mean(data, na.rm=TRUE)), by='year']

# ocean nino index
oni = data.table(read.csv( paste0(data_dir, "cciea_OC_ONI.csv" )))
oni = oni[-1,]
oni$year = lubridate::year(oni$time)
oni$month = lubridate::month(oni$time)
# year  = substring(oni$time, 1,4)
oni$ONI <- as.numeric(as.character(oni$ONI))
oni[, oni_pre:=ifelse(month %in% c(7:12), year+1,year)]


oni_pre = oni[ month %in% c(1:3, 7:12), list(oni_pre = mean(ONI, na.rm=TRUE)), by='oni_pre']
colnames(oni_pre)[1] <- 'year'
oni_larv = oni[ month %in% c(2:3), list(oni_larv = mean(ONI, na.rm=TRUE)), by='year']
oni_pjuv = oni[ month %in% c(4:8), list(oni_pjuv = mean(ONI, na.rm=TRUE)), by='year']

#PDO 
pdo = data.table(read.csv( paste0(data_dir, "cciea_OC_PDO.csv" )))
pdo = pdo[-1]
pdo$year = lubridate::year(pdo$time)
pdo$month = lubridate::month(pdo$time)
pdo$PDO <- as.numeric(as.character(pdo$PDO))
pdo_larv = pdo[ month %in% c(2:3), list(pdo_larv = mean(PDO, na.rm=TRUE)), by='year']
pdo_pjuv = pdo[ month %in% c(4:8), list(pdo_pjuv = mean(PDO, na.rm=TRUE)), by='year']

#BEUTI
beuti = data.table(read.csv( paste0(data_dir, "oc_beuti_M_FY2024.csv" )))
beuti_39_45 = beuti[beuti$timeseries != "33N",]
beuti_45 = beuti_39_45[beuti_39_45$timeseries != "39N",]
beuti_45$year = lubridate::year(beuti_45$time)
beuti_45$month = lubridate::month(beuti_45$time)
beuti_larv = beuti_45[ month %in% c(2:3), list(oni_larv = mean(index, na.rm=TRUE)), by='year']
beuti_pjuv = beuti_45[ month %in% c(4:8), list(oni_pjuv = mean(index, na.rm=TRUE)), by='year']

# LUSI (annual)
lusi = data.table(read.csv( paste0("data-yellowtail/cciea_OC_LUSI.csv" )))
lusi = lusi[-1,]
lusi$year = lubridate::year(lusi$time)
lusi$lusi <- as.numeric(as.character(lusi$lusi))
lusi_annual = lusi[ year %in% c(1993:2024), list(lusi_annual = mean(lusi, na.rm=TRUE)), by='year']

# Zooplankton
zoo = data.table(read.csv( paste0("data-yellowtail/Zoop_Newport_copepods_FY2024.csv")))
zoo_NBA = zoo[zoo$timeseries != "NorthernBiomassAnomaly",]
ZOOpjuvS = zoo_NBA[ month %in% c(4:8), list(ZOOpjuvS = mean(index, na.rm=TRUE)), by='year']
ZOObenS = zoo_NBA[ month %in% c(9:12), list(ZOObenS = mean(index, na.rm=TRUE)), by='year']

zoo_NBA = zoo[zoo$timeseries != "SouthernBiomassAnomaly",]
ZOOpjuvN = zoo_NBA[ month %in% c(4:8), list(ZOOpjuvN = mean(index, na.rm=TRUE)), by='year']
ZOObenN = zoo_NBA[ month %in% c(9:12), list(ZOObenN = mean(index, na.rm=TRUE)), by='year']

# Environmental Data
env = data.table(read.csv( paste0("data-yellowtail/Processed_Environmental/2024update_Env_Indices.csv")))

glorys_ts = left_join(glorys_ts, hci1_larv)
glorys_ts = left_join(glorys_ts, hci1_pjuv)
glorys_ts = left_join(glorys_ts, hci2_larv)
glorys_ts = left_join(glorys_ts, hci2_pjuv)
glorys_ts = left_join(glorys_ts, oni_pre)
glorys_ts = left_join(glorys_ts, oni_larv)
glorys_ts = left_join(glorys_ts, oni_pjuv)
glorys_ts = left_join(glorys_ts, pdo_larv)
glorys_ts = left_join(glorys_ts, pdo_pjuv)
#glorys_ts = left_join(glorys_ts, beuti_larv)
#glorys_ts = left_join(glorys_ts, beuti_pjuv)
glorys_ts = left_join(glorys_ts, lusi_annual)
glorys_ts = left_join(glorys_ts, ZOOpjuvN)
glorys_ts = left_join(glorys_ts, ZOObenN)
glorys_ts = left_join(glorys_ts, ZOOpjuvS)
glorys_ts = left_join(glorys_ts, ZOObenS)
glorys_ts = left_join(glorys_ts, env)

fwrite(glorys_ts, paste0("2024Env-annual-yellowtail.csv"))

standardized <-data.frame(glorys_ts)%>%
  dplyr::select(!any_of(c('year')))%>% 
  mutate_all(~ scale(.))
standardized_env<-cbind(year=data.frame(glorys_ts)$year,standardized)

fwrite(standardized_env, paste0("2024Env-annual-yellowtail-standardized.csv"))


