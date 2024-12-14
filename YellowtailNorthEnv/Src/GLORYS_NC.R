library(ncdf4) # package for netcdf manipulation
library(terra)
library(tidyr)
library(dplyr)

nc <- nctidyrnc <- nc_open('DATAFILE_cuti_daily_glorys_1993-2019.nc') # read in the file

#pulling out dimensions for nc
v <- nc$var[[1]]
size <- v$varsize 
dims <- v$ndims

# get lat/long
lat <- ncvar_get(nc, "latitude") #pulling out latitude (we don't have long)

# get times 
raw <- ncvar_get(nc, "time") #pulling out the raw times, these are just indices
tunits<-ncatt_get(nc,"time",attname="units") #getting the units the raw names refer to
dates <- RNetCDF::utcal.nc(tunits$value, raw) #getting the dates split into year/month/day
dates <- as.data.frame(dates[,c("year","month", "day")]) # we are working with daily data so lets ditch hour/minute/second columns

#defining variables and dims for the dataframe
var_name = "CUTI" #defining which varuiable we want from .nc
nt<- length(raw) #number of time steps for each location

#looping over the timesteps and pulling out CUTI obs for all latitudes (columns) for each time step (rows)

CUTI_temp <-data.frame() #temp datafreme to fill in with the loop
for (i in 1:nt) {
  start <- rep(1,dims)     # begin with start=(1,1,...,1)
  start[dims] <- i             # change to start=(1,1,...,i) to read    timestep i
  count <- size                # begin with count=(nx,ny,...,nt), reads entire var
  count[dims] <- 1             # change to count=(nx,ny,...,1) to read 1 tstep
  temp<-ncvar_get(nc, varid = var_name, start = start, count = count) #temp var that has CUTI for each lat for timestep i
  CUTI_temp<-rbind(CUTI_temp, temp) #putting the temp var into the temp dataframe and stiching it together for ever i of the loop so it is saved
}

#checking the output and adding names
dim(CUTI_temp) #checking to make sure we have all the data we think we need...
#just cleaning up the dataframe
colnames(CUTI_temp)<-gsub(" ","",paste(lat,'N'))#setting latitude as our rownames but adding in 'N' so it matches Ellens file
CUTI_output <-cbind(dates,CUTI_temp) #combining our restructured CUTI data with the appropriate dates from the dates dataframe


write.csv(CUTI_output, "CUTI_daily_GLORYS.csv")
