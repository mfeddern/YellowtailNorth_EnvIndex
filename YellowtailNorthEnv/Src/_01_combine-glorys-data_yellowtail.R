
source("_00_yellowtail-header.r")

# note on data.table format

# new.file = old.file[sub-setting, calculation, grouping]
# new.file = old.file[bottom_depth <=549 & lat >=40, .(mld = mean(mlotst)), by='date']  
# new.file = old.file[bottom_depth <=549 & lat >=40, list(mld = mean(mlotst)), by='date']  

# add some helper functions ####
if(str_detect(home_dir,"GitHub")==TRUE){
  source("~/GitHub/glorys-download/glorys_helper_functions.r")
  }else{source("~/Documents/glorys-download/glorys_helper_functions.r")}

# bring in bottom depth to add to other files ##################################
# bathymetry from GLORYS
df <- tidync::tidync(paste0(data_raw,'glorys-bathymetry.nc')) %>% 
  hyper_tibble( force = TRUE) %>%
  drop_na() %>% 
  group_by(longitude,latitude)

# head(df)
dt_bathy <- data.table(df)
dt_bathy = dt_bathy %>% rename(lat = latitude, lon = longitude, bottom_depth = deptho)
head(dt_bathy)
tail(dt_bathy)
fwrite(dt_bathy, paste0(data_dir,"bathymetry-m.csv" ))
rm(df)

# distance from shore ##########################################################
dist_shore = data.frame(read.csv(
  "~/GitHub/glorys-download/GLORYS distance to shore 15May2024.csv", header=TRUE))
head(dist_shore)
d_shore = dist_shore[,c("longitude","latitude", 'Distance.to.shore..km.', 'Water_land')]
d_shore = d_shore %>% rename(lat=latitude, lon=longitude, dist_shore_km = Distance.to.shore..km.)
head(d_shore)
fwrite(d_shore, paste0(data_dir,'distance-to-shore.csv'))

################################################################################
# DAILY MIXED-LAYER-DEPTH ######################################################
################################################################################

# Note need separate monthly and daily files below. Keep separate.

# dir()
# combine reanalysis and interim files

df1 <- tidync::tidync(paste0(data_raw, "glorys-daily-MLD-1993-01-01-2021-06-30.nc")) %>% 
  hyper_tibble( force = TRUE) %>%
  drop_na() %>% 
  group_by(longitude,latitude,time)

df2 <- tidync::tidync(paste0(data_raw, "glorys-daily-MLD-2021-07-01-2023-12-31.nc")) %>% 
  hyper_tibble( force = TRUE) %>%
  drop_na() %>% 
  group_by(longitude,latitude,time)

dt_x = rbindlist( list(df1, df2)) # faster; now a data.table
# clean up for space
rm(df1, df2)

# helper function to clean up date #############################################
mld_daily <- data_prep(data.file=dt_x, 
                    bathy_file = dt_bathy, 
                    shore_dist = d_shore)
rm(dt_x)

mld_daily_processed = mld_daily[bottom_depth <=549 & lat >40, .(mld = mean(mlotst)), by='date']  

head(mld_daily_processed)
dim(mld_daily_processed)

# to combine with other files and output
# glorys_daily_means <- mld_daily

# mld processed for any time series
# mld_daily for combining with other files because it retains the lat long

fwrite(mld_daily_processed, paste0(data_dir, "glorys-daily-means-mld.csv")) 

rm(mld_daily) # clean up; might need to retain if subsetting by mld below
rm(mld_daily_processed)

################################################################################
# MONTHY MIXED-LAYER-DEPTH #####################################################
################################################################################

# dir(data_raw)
# for filtering monthly transport below

df1 <- tidync::tidync(paste0(data_raw, "glorys-monthly-MLD-1993-01-01-2021-06-30-1200m.nc")) %>% 
  hyper_tibble( force = TRUE) %>%
  drop_na() %>% 
  group_by(longitude,latitude,time)

df2 <- tidync::tidync(paste0(data_raw, "glorys-monthly-MLD-2021-07-01-2023-12-31-1200m.nc")) %>% 
  hyper_tibble( force = TRUE) %>%
  drop_na() %>% 
  group_by(longitude,latitude,time) 

dt = rbindlist( list(df1, df2)) # faster; now a data.table
# clean up for space
rm(df1, df2)


mld_monthly <- data_prep(data.file=dt, bathy_file = dt_bathy)

rm(dt)
head(mld_monthly)
dim(mld_monthly)

mld_monthly = mld_monthly[lat>40,,]

# use mld_monthly to combine with other files because it contains lat long info

# mld_monthly_processed for time series if included
mld_0_180 = mld_monthly[bottom_depth >= 0 & bottom_depth <=180,.(mld_0_180 = mean(mlotst)), by='date']
mld_0_90  = mld_monthly[bottom_depth >= 0 & bottom_depth <=90, .(mld_0_90 = mean(mlotst)), by='date']
mld_30_130  = mld_monthly[bottom_depth >= 30 & bottom_depth <=130 , .(mld_30_130 = mean(mlotst)), by='date']

# no loop here so don't need the if statement

mld_month_out = left_join(mld_0_180, mld_0_90)
mld_month_out = left_join(mld_month_out, mld_30_130)

fwrite(mld_month_out, paste0(data_dir, "glorys-monthly-means-mld.csv")) 

rm(mld_month_out, mld_0_180, mld_0_90, mld_30_130)

# continue here ####

################################################################################
# DAILY WATER COLUMN TEMPERATURE ###############################################
################################################################################
x = dir(data_raw)
y = grep("watercolumn-temp", x)
x = x[y]
x
# reference values for degree days
# Yamada, J., Kusakari, M. Staging and the time course of embryonic development in kurosoi,Sebastes schlegeli . 
# Environ Biol Fish 30, 103–110 (1991). https://doi.org/10.1007/BF02296881

# Kusakari (1995) Y=−0.00142X3+48.0028 @ 9.8C

DDref = 4 # just a default at present

for(i in 1:length(x)){  # run for each year because the data are large # just output the time daily time series
  print(paste0("i = ",i))
  begin.time = Sys.time()
  print( paste0( "Start time = ",begin.time) )
  dfx <- tidync::tidync(paste0(data_raw,x[i])) %>% 
    hyper_tibble( force = TRUE) %>%
    drop_na() %>% 
    group_by(longitude,latitude,time, depth)
  
  # reduce size immediately
  dfx = dfx[dfx$latitude>40,] # not a data table yet
  #shortcut using data_prep
  (prep.time = Sys.time())
  dt_temp <- data_prep(data.file=dfx, bathy_file = dt_bathy)
  
  rm(dfx)
  
  # dt_temp = dt_temp[lat>40,,]
  
  print( paste0( "Starting averaging at ", Sys.time() ) )
  
  temp_90_180 = dt_temp[bottom_depth >= 90 & bottom_depth <=180 , .(temp_90_180 = mean(thetao)), by='date']
  temp_0_180 = dt_temp[bottom_depth >= 0 & bottom_depth <=180 , .(temp_0_180 = mean(thetao)), by='date']
  temp_0_90 = dt_temp[bottom_depth >= 0 & bottom_depth <=90 , .(temp_0_90 = mean(thetao)), by='date']
  temp_30_130 = dt_temp[bottom_depth >= 30 & bottom_depth <=130 , .(temp_30_130 = mean(thetao)), by='date']
  temp_180_549 = dt_temp[bottom_depth >= 180 & bottom_depth <=549 , .(temp_180_549 = mean(thetao)), by='date']
  
  dd_90_180 = dt_temp[bottom_depth >= 90 & bottom_depth <=180 , .(dd_90_180 = sum(thetao-DDref)), by='date']
  dd_0_180 = dt_temp[bottom_depth >= 0 & bottom_depth <=180 , .(dd_0_180 = sum(thetao-DDref)), by='date']
  dd_0_90 = dt_temp[bottom_depth >= 0 & bottom_depth <=90 , .(dd_0_90 = sum(thetao-DDref)), by='date']
  dd_30_130 = dt_temp[bottom_depth >= 30 & bottom_depth <=130 , .(dd_30_130 = sum(thetao-DDref)), by='date']
  dd_180_549 = dt_temp[bottom_depth >= 180 & bottom_depth <=549 , .(dd_180_549 = sum(thetao-DDref)), by='date']
  
  rm(dt_temp) # clean up for memory space
  
  if(i==1){
    Temp_90_180 = temp_90_180
    Temp_0_180 = temp_0_180
    Temp_0_90 = temp_0_90
    Temp_30_130 = temp_30_130
    Temp_180_549 = temp_180_549
    
    DD_90_180 = dd_90_180
    DD_0_180 = dd_0_180
    DD_0_90 = dd_0_90
    DD_30_130 = dd_30_130
    DD_180_549 = dd_180_549
    
  }else{
    Temp_90_180 = rbindlist( list(Temp_90_180,temp_90_180))
    Temp_0_180 = rbindlist( list(Temp_0_180,temp_0_180))
    Temp_0_90 = rbindlist( list(Temp_0_90,temp_0_90))
    Temp_30_130 = rbindlist( list(Temp_30_130,temp_30_130))
    Temp_180_549 = rbindlist( list(Temp_180_549,temp_180_549))
    
    DD_90_180 = rbindlist( list(DD_90_180,dd_90_180))
    DD_0_180 = rbindlist( list(DD_0_180,dd_0_180))
    DD_0_90 = rbindlist( list(DD_0_90,dd_0_90))
    DD_30_130 = rbindlist( list(DD_30_130,dd_30_130))
    DD_180_549 = rbindlist( list(DD_180_549,dd_180_549))
    
    #remove old files
    rm(temp_90_180,temp_0_180,temp_0_90,temp_30_130,temp_180_549,
       dd_90_180,dd_0_180,dd_0_90,dd_30_130,dd_180_549)
    
  } # end if
  (stop.time = Sys.time())
  print( paste0('Elapsed time = ', (stop.time - begin.time)))
} # end i loop

temp_daily_out = left_join(Temp_90_180,Temp_0_180)
temp_daily_out = left_join(temp_daily_out,Temp_0_90)
temp_daily_out = left_join(temp_daily_out,Temp_30_130)
temp_daily_out = left_join(temp_daily_out,Temp_180_549)
temp_daily_out = left_join(temp_daily_out,DD_90_180)
temp_daily_out = left_join(temp_daily_out,DD_0_180)
temp_daily_out = left_join(temp_daily_out,DD_0_90)
temp_daily_out = left_join(temp_daily_out,DD_30_130)
temp_daily_out = left_join(temp_daily_out,DD_180_549)

# clean up to save memory
fwrite(temp_daily_out, paste0(data_dir, "glorys-daily-temperature.csv"))

rm(temp_daily_out)

################################################################################
# MONTHLY CROSS SHELF TRANSPORT ################################################
################################################################################

x = dir(data_raw)
y = grep("monthly-crossshelf", x)
x = x[y]
x

for(i in 1:length(x)){ 
  print(x[i])
  (t1 = Sys.time())
  dfx <- tidync::tidync(paste0(data_raw,x[i])) %>% 
    hyper_tibble( force = TRUE) %>%
    drop_na() %>% 
    group_by(longitude,latitude,time)
  
  dt_cst = data_prep(data.file=dfx, bathy_file = dt_bathy)
  rm(dfx) # clean up memory for space
  # as in temporary not temperature
  
  # calculations for depth x lat zones = mean daily value
  # subset and average by date later
  
  cst_90_180 = dt_cst[bottom_depth >= 90 & bottom_depth <=180 , .(cst_90_180 = mean(uo)), by='date']
  cst_0_180 = dt_cst[bottom_depth >= 0 & bottom_depth <=180 , .(cst_0_180 = mean(uo)), by='date']
  cst_0_90 = dt_cst[bottom_depth >= 0 & bottom_depth <=90 , .(cst_0_90 = mean(uo)), by='date']
  cst_30_130 = dt_cst[bottom_depth >= 30 & bottom_depth <=130 , .(cst_30_130 = mean(uo)), by='date']
  cst_180_549 = dt_cst[bottom_depth >= 180 & bottom_depth <=549 , .(cst_180_549 = mean(uo)), by='date']
  
  rm(dt_cst)
  
  if(i==1){
    CST_90_180 = cst_90_180
    CST_0_180 = cst_0_180
    CST_0_90 = cst_0_90
    CST_30_130 = cst_30_130
    CST_180_549 = cst_180_549
  }else{
    CST_90_180 = rbindlist( list(CST_90_180,cst_90_180))
    CST_0_180 = rbindlist( list(CST_0_180,cst_0_180))
    CST_0_90 = rbindlist( list(CST_0_90,cst_0_90))
    CST_30_130 = rbindlist( list(CST_30_130,cst_30_130))
    CST_180_549 = rbindlist( list(CST_180_549,cst_180_549))
  } # end if
} # end i loop

cst_monthly = left_join(CST_90_180,CST_0_180)
cst_monthly = left_join(cst_monthly,CST_0_90)
cst_monthly = left_join(cst_monthly,CST_30_130)
cst_monthly = left_join(cst_monthly,CST_180_549)

# clean up to save memory
fwrite(cst_monthly, paste0(data_dir, "glorys-monthly-means-cst.csv")) 
rm(cst_monthly)

################################################################################
# MONTHLY LONGSHORE TRANSPORT ##################################################
################################################################################

# need to go through this and sea surface height

# glorys_daily_means = fread( paste0(data_dir, "glorys-daily-means-3.csv"))
# glorys_daily_means$date = lubridate::as_date(glorys_daily_means$date)

x = dir(data_raw)
y = grep("monthly-longshore" , x)
x = x[y]
x

for(i in 1:length(x)){ 
  print(x[i])
  (t1 = Sys.time())
  dfx <- tidync::tidync(paste0(data_raw,x[i])) %>% 
    hyper_tibble( force = TRUE) %>%
    drop_na() %>% 
    group_by(longitude,latitude,time)
  
  dt_lst = data_prep(data.file=dfx, bathy_file = dt_bathy)
  rm(dfx) # clean up memory for space
  # as in temporary not temperature
  
  # calculations for depth x lat zones = mean daily value
  # subset and average by date later
  
  lst_90_180 = dt_lst[bottom_depth >= 90 & bottom_depth <=180 , .(lst_90_180 = mean(vo)), by='date']
  lst_0_180 = dt_lst[bottom_depth >= 0 & bottom_depth <=180 , .(lst_0_180 = mean(vo)), by='date']
  lst_0_90 = dt_lst[bottom_depth >= 0 & bottom_depth <=90 , .(lst_0_90 = mean(vo)), by='date']
  lst_30_130 = dt_lst[bottom_depth >= 30 & bottom_depth <=130 , .(lst_30_130 = mean(vo)), by='date']
  lst_180_549 = dt_lst[bottom_depth >= 180 & bottom_depth <=549 , .(lst_180_549 = mean(vo)), by='date']
  
  rm(dt_lst)
  
  if(i==1){
    LST_90_180 = lst_90_180
    LST_0_180 = lst_0_180
    LST_0_90 = lst_0_90
    LST_30_130 = lst_30_130
    LST_180_549 = lst_180_549
  }else{
    LST_90_180 = rbindlist( list(LST_90_180,lst_90_180))
    LST_0_180 = rbindlist( list(LST_0_180,lst_0_180))
    LST_0_90 = rbindlist( list(LST_0_90,lst_0_90))
    LST_30_130 = rbindlist( list(LST_30_130,lst_30_130))
    LST_180_549 = rbindlist( list(LST_180_549,lst_180_549))
  } # end if
} # end i loop



lst_monthly = left_join(LST_90_180,LST_0_180)
lst_monthly = left_join(lst_monthly,LST_0_90)
lst_monthly = left_join(lst_monthly,LST_30_130)
lst_monthly = left_join(lst_monthly,LST_180_549)

# clean up to save memory
fwrite(lst_monthly, paste0(data_dir, "glorys-monthly-means-lst.csv")) 
rm(lst_monthly)
