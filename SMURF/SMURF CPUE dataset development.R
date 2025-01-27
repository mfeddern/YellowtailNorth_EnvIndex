
###################################################

# PROJECT:    CREATION OF A SMURF RECRUITMENT INDEX
# AUTHORS:    ALISON WHITMAN

# WHAT DOES THIS SCRIPT DO?
#   (1) COMBINES DATA FROM OSU SMURF TEAM ON SEBASTES FLAVIDUS SETTLEMENT WITH 
#   OCEANOGRPAHY DATA FROM ODFW MOORINGS AT SPECIFIC MARINE RESERVES 

# UPDATED: 01/17/2025 BY A. WHITMAN

###################################################

### LOAD PACKAGES ###
library(easypackages)
libraries("lattice", "Rmisc", "ggplot2", "tidyr", "reshape2", "readxl", 
          "writexl", "zoo", "RODBC", "splitstackshape",
          "lubridate","plyr","dplyr")

# BRING IN DATA -----------------------------------------------------------


setwd("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS")

ocean<-read.csv("oceanography.csv") # ODFW oceanography data, note updated to corrected data on 1/14

mooring_central<-read.csv("mooring_locations_central.csv") # central moorings 
settle_central<-read.csv("SFLA_settlement_central.csv") # central settlement

mooring_south<-read.csv("mooring_locations_south.csv") # south moorings 
settle_south<-read.csv("SFLA_settlement_south.csv") # south settlement


str(settle_central)
str(settle_south)

## COMBINE CENTRAL AND SOUTHERN SETTLEMENT DATASETS

# have to rename some fields
names(settle_central)[names(settle_central) == 'Site.Code'] <- 'Site'
names(settle_south)[names(settle_south) == 'sfla_settlement_rate'] <- 'SFLA_settlement_rate'

# fix the dates
settle_central$Date<-as.POSIXct(settle_central$Date, format = "%m/%d/%Y")
settle_south$Date<-as.POSIXct(settle_south$Date, format = "%m/%d/%y")

settle<-rbind(settle_central,settle_south)
str(settle)

# merge the mooring locations for the sdmTMB code
str(mooring_central)
str(mooring_south)

# have to rename ALL THE fields

oldnames = names(mooring_central)
newnames = names(mooring_south)

mooring_central<-mooring_central %>% 
  rename_at(vars(all_of(oldnames)), ~ newnames)

moorings<-rbind(mooring_central,mooring_south)

## COMBINE WITH OCEANOGRAPHY DATA

# note that the oceanography data is a finer temporal scale than the settlements
# which are only recorded by the date

# Ryan Fields (ODFW- MR program) gave me both a raw file and one with a day-level average 
# for all three instruments

str(ocean)
str(settle)

# will need to match the average oceanography data with the date and location of the 
# settlement areas

# how do I match location between the two datasets? do I need the mooring lat/longs? 
# emails since 01/10/2025 indicate that mooring codes should work 
# (not all moorings have oceanography data associated)
# but i should be able to link them with detailed site code and sampling date

# bring in the odfw moorings information to match the codes correctly 

odfw_moorings<-read.csv("odfw moorings.csv")
head(odfw_moorings)

# also noticed that year does not match the year in the date for the oceanography data 
# ** sent email to Ryan about this on 1/9/25 - using corrected version now (1/14)

# pull out year for the settlement data
settle$year<-as.numeric(format(settle$Date,"%Y"))

summary(settle$year)
summary(ocean$year) 

ocean<-ocean[ocean$year>2010,] # removed early years that we don't need

settle$Site<-as.factor(settle$Site)
summary(settle$Site) # all four sites
settle$Mooring.Code<-as.factor(settle$Mooring.Code)
summary(settle$Mooring.Code)

odfw_moorings$Site.Code<-as.factor(odfw_moorings$Site.Code)
odfw_moorings$SMURF.Mooring.Name<-as.factor(odfw_moorings$SMURF.Mooring.Name)
summary(odfw_moorings$Site.Code)
summary(odfw_moorings$SMURF.Mooring.Name)

# add detailed site code to settlement data via odfw moorings info

settle$Site.Code<-odfw_moorings$Site.Code[match(settle$Mooring.Code,odfw_moorings$SMURF.Mooring.Name)]
summary(settle$Site.Code)
settle$Site.Code<-droplevels(settle$Site.Code)

# now the oceanography data can be linked via the date and the detailed site code 
# i think it would be cleaner to link via julian day, year and site code

# so add julian to settlement data
settle$julian<-as.numeric(format(settle$Date,"%j"))
summary(settle$julian)


# remove the sites we don't need from the oceanography data 
# and see if we have multi-depth sampling days after that
ocean<-ocean[ocean$site %in% c("OR","ORF","RR","RRH"),]

ocean$site<-as.factor(ocean$site)
ocean$ID<-paste0(ocean$year,"_",ocean$julian,"_",ocean$site)

depths<-as.data.frame.matrix(table(ocean$ID,ocean$depth_m))
depths$total<-rowSums(depths)
table(depths$total) # yes we have multiple depths of sampling for each date and site 

# so going to need to restructure the oceanography data to include more than one depth
# not entirely sure how to do this. 

# ryan noted that these depths are not perfect, so should categorize them 

unique(ocean$depth_m)
ocean$depth_m<-as.factor(ocean$depth_m)

# next steps, categorize depth, spread the data by ID ~ depth bin (?) or filter by depth bin
# then link with settlement data

ocean$depth_cat<-ifelse(ocean$depth_m==1.0,"shallow",
                        ifelse(ocean$depth_m==15,"deep","mid"))
summary(as.factor(ocean$depth_cat))

# reorganize to wrap my head around this a bit 

ocean<-ocean[,c(1:6,10:12,7:9)]
names(ocean)

# an easier way might be to filter to mid-water column only and see how much matches up 
# not sure about how to code the idea of keeping at least one measurement, regardless of depth, 
# and then only including the mid-water if i have more than one 

# notes from meeting with Meghan on 1/17 - she will work on a state-space model with the oceanographic 
# data to fill in holes in space and time 
# but I should also try my nested ifelse statements to fill in the holes where I can with the data 
# we do have already from multiple depths 

# SUBSETTING TO MID-WATER DEPTHS ONLY FOR A TEST RUN 

ocean_test<-ocean[ocean$depth_cat=="mid",]
depths_test<-as.data.frame.matrix(table(ocean_test$ID,ocean_test$depth_m))
depths_test$total<-rowSums(depths_test)
table(depths_test$total) # so that excludes all data where we have multiple sampled depths

# COMBINING OCEANOGRAPHIC DATA WITH SETTLEMENT DATA 

str(settle)
str(ocean_test) # note that I'm using the mid-water only depths test dataset

# create the ID for the settlement dataset

settle$ID<-paste0(settle$year,"_",settle$julian,"_",settle$Site.Code)

# add the oceanographic fields I want to the settlement dataframe - use the ID to match

settle_ocean<-merge(settle, ocean_test[,c(8:12)], by = c("ID"), all.x = TRUE)

names(settle_ocean)
head(settle_ocean)
summary(settle_ocean)

# EXPORT DATA! WOOHOO! 
# though I may be back here shortly to see if I can match more oceanographic data...

#write.csv(settle_ocean,"combined_settlement_ocean.csv",row.names = F)
