
###################################################

# PROJECT:    CREATION OF A SMURF RECRUITMENT INDEX
# AUTHORS:    ALISON WHITMAN & MEGAN FEDDERN

# WHAT DOES THIS SCRIPT DO?
#   (1) COMBINES DATA FROM OSU SMURF TEAM ON SEBASTES FLAVIDUS SETTLEMENT WITH 
#   OCEANOGRPAHY DATA FROM ODFW MOORINGS AT SPECIFIC MARINE RESERVES 

# UPDATED: 02/20/2025 BY M. FEDDERN

###################################################

### LOAD PACKAGES ###
library(easypackages)
libraries("lattice", "Rmisc", "ggplot2", "tidyr", "reshape2", "readxl", 
          "writexl", "zoo", "RODBC", "splitstackshape",
          "lubridate","plyr","dplyr","corrplot")
# BRING IN DATA -----------------------------------------------------------


#setwd("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS")

#ocean<-read.csv("oceanography.csv") # ODFW oceanography data, note updated to corrected data on 1/14
ocean<-read.csv('Oceanography_DailyMean.csv')
mooring_central<-read.csv("mooring_locations_central.csv") # central moorings 
settle_central<-read.csv("SFLA_settlement_central.csv")


mooring_south<-read.csv("mooring_locations_south.csv") # south moorings 
settle_south<-read.csv("SFLA_settlement_south.csv") # south settlement


str(settle_central)
str(settle_south)

## COMBINE CENTRAL AND SOUTHERN SETTLEMENT DATASETS

# have to rename some fields
names(settle_central)[names(settle_central) == 'Site.Code'] <- 'Site'
names(settle_south)[names(settle_south) == 'sfla_settlement_rate'] <- 'SFLA_settlement_rate'

# fix the dates

settle_south$Date<-as.POSIXct(settle_south$Date, format = "%m/%d/%y")
settle_central$Date<-as.POSIXct(settle_central$Date, format = "%m/%d/%y")

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
year(settle$Date)

summary(settle$year)
summary(ocean$year) 
unique(settle$year)
unique(ocean$year)
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

# Evaluating the correlation across deths, location

marine_reserves<-unique(ocean$marine_reserve) # create a reserve var
loc<-unique(ocean$location) # create a location var


  for(i in 1:2){
    p<-ggplot(ocean%>%filter(marine_reserve==marine_reserves[i]), 
              aes(x = julian, y = temp_c, group=as.factor(depth_cat), col=as.factor(depth_cat)))+
      geom_line()+ 
      ggtitle(marine_reserves[i])+
      facet_wrap(~year)
    plot(p)
  } 



for(i in 1:4){
  p<-ggplot(ocean%>%filter(location==loc[i]), 
            aes(x = julian, y = temp_c, group=as.factor(depth_cat), col=as.factor(depth_cat)))+
    geom_line()+ 
    ggtitle(loc[i])+
    facet_wrap(~year)
  plot(p)
}

#lets look first at how these correlate based on location
cor_dat_loc<-ocean%>%filter(depth_cat=="mid")%>%
  select(year,julian,location, temp_c)%>%
  pivot_wider(names_from=location, values_from=temp_c)
df.subset <- cor_dat_loc[, names(cor_dat_loc)[(names(cor_dat_loc) %in% loc)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

# based on location - these are very correlated. CF and OR 0.99 and RR and HH 0.96
# We will get singularity issues if we try to fit these separately as the contain
# nearly identical information. Lets start by combining across reserve/control
# we will aggregate the data by taking the mean

ocean_reserve=ocean%>%
  group_by(depth_cat,date,julian,year,marine_reserve)%>%
  dplyr::summarise(temp=mean(na.omit(temp_c)))%>%
  ungroup()
#reserve_plot(ocean_reserve)
# now we have reserve-specific temp data. Next we can look at the correlation at 
# different depths. We really only need one "water column" index, so we do not
# need to retain depth data if it is highly correlated

cor_dat_depth<-ocean_reserve%>%filter(marine_reserve=="Redfish Rocks")%>%
  select(year,julian,depth_cat, temp)%>%
  pivot_wider(names_from=depth_cat, values_from=temp)
depth<-unique(ocean_reserve$depth_cat)
df.subset <- cor_dat_depth[, names(cor_dat_depth)[(names(cor_dat_depth) %in% depth)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

cor_dat_depth<-ocean_reserve%>%filter(marine_reserve=="Otter Rock")%>%
  select(year,julian,depth_cat, temp)%>%
  pivot_wider(names_from=depth_cat, values_from=temp)
depth<-unique(ocean_reserve$depth_cat)
df.subset <- cor_dat_depth[, names(cor_dat_depth)[(names(cor_dat_depth) %in% depth)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

# Our depth data is also insanely correlated. After chatting with Eric W., we 
# think a MARSS model might be over kill. There also is not a lot of shared
# seasonality or variability across years, so lets make a single variable across
# depths

# lets start by getting teh means for depths so we can standardize the data, which
# within each reserve and depth cat, this will help control for having years with 
# different depths 

mean_tab<- ocean_reserve%>%
  dplyr::group_by(marine_reserve, depth_cat)%>%
  dplyr::summarise(mean=mean(na.omit(temp)), sd=sd(na.omit(temp)))

ocean_scaled<- ocean_reserve%>%
  left_join(mean_tab)%>%
  mutate(temp_scale=(temp-mean)/sd)
#reserve_plot(ocean_scaled)
# lets replot this to make sure it looks right...
for(i in 1:2){
  p<-ggplot(ocean_scaled%>%filter(marine_reserve==marine_reserves[i]), 
            aes(x = julian, y = temp_scale, group=as.factor(depth_cat), col=as.factor(depth_cat)))+
    geom_line()+ 
    ggtitle(marine_reserves[i])+
    facet_wrap(~year)
  plot(p)
}

# It seems pretty reasonable to aggregate these, we can retain the mean values and we can back transform/unscale
# if we want to. For the sake of our sanity, lets set the NA values to the mean 
# which will mean it will have minmal impact on our results. We can get fancier if
# the model runs; lets also expand our grid so we have observations for every day
#dummy_jd <-data.frame(julian=seq(1,365,1))
mean_mid<-mean_tab%>%filter(depth_cat=='mid')

temp_ind<- ocean_scaled%>%
  dplyr::group_by(marine_reserve, year, julian)%>%
  dplyr::summarise(temp_index=mean(temp_scale))%>%
  dplyr::ungroup()%>%
  left_join(mean_mid)%>%
  mutate(temp_c_mid = (temp_index*sd)+mean)

#temp_index<-merge(temp_index, dummy_jd,all.y = TRUE)
temp_ind[is.na(temp_index)] <- 0

for(i in 1:2){
  p<-ggplot(temp_ind%>%filter(marine_reserve==marine_reserves[i]), 
            aes(x = julian, y = temp_index))+
    geom_line()+ 
    ggtitle(marine_reserves[i])+
    facet_wrap(~year)
  plot(p)
}

# okay now how should we summarise data for the settlement? I am going to just use a simple
# rolling window. If the var is important and we get the SDMtmb model to fit 
# I can make this better
mean(settle$Sampling.Interval) # mean sampling interval is 15.5 
min(settle$Sampling.Interval) # minn sampling interval is 7 days so lets use a 7 day rolling window

reserve=unique(temp_ind$marine_reserve)
years=unique(temp_ind$year)
rolling_temp<-data.frame()
for(i in 1:2){
  for(j in 1:length(years)){
    temp_dat=temp_ind%>%filter(year==years[j]& marine_reserve==reserve[i])
temp<-rollapply(temp_dat%>%
            select(temp_index), 8, mean,fill=NA,na.rm=T)
temp1<-rollapply(temp_dat%>%
                  select(temp_c_mid), 8, sum,fill=NA,na.rm=T)
temp2<-data.frame(rolling16 =temp, temp1,marine_reserve=reserve[i],year=years[j],julian=temp_dat$julian+4)

rolling_temp=rbind(rolling_temp,temp2)
}
}
rolling_temp3<-rolling_temp%>%dplyr::rename(rolling8d=temp_index, cdd_8=temp_c_mid)
ocean_rolling=temp_ind%>%full_join(rolling_temp)

rolling_temp<-data.frame()
for(i in 1:2){
  for(j in 1:length(years)){
    temp_dat=temp_ind%>%filter(year==years[j]& marine_reserve==reserve[i])
    temp<-rollapply(temp_dat%>%
                      select(temp_index), 16, mean,fill=NA,na.rm=T)
    temp1<-rollapply(temp_dat%>%
                       select(temp_c_mid), 16, sum,fill=NA,na.rm=T)
    temp2<-data.frame(rolling16 =temp, temp1,marine_reserve=reserve[i],year=years[j],julian=temp_dat$julian+8)
    
    rolling_temp=rbind(rolling_temp,temp2)
  }
}
rolling_temp2<-rolling_temp%>%dplyr::rename(rolling16d=temp_index, cdd_16=temp_c_mid)

ocean_rolling=temp_ind%>%
  full_join(rolling_temp3)%>%
  full_join(rolling_temp2)

# an easier way might be to filter to mid-water column only and see how much matches up 
# not sure about how to code the idea of keeping at least one measurement, regardless of depth, 
# and then only including the mid-water if i have more than one 

# notes from meeting with Meghan on 1/17 - she will work on a state-space model with the oceanographic 
# data to fill in holes in space and time 
# but I should also try my nested ifelse statements to fill in the holes where I can with the data 
# we do have already from multiple depths 

# SUBSETTING TO MID-WATER DEPTHS ONLY FOR A TEST RUN 

ocean_full<-ocean%>%select(marine_reserve, site, year, date, julian, ID)%>%
  dplyr::rename(Site.Code=site)%>%distinct()%>%
  dplyr::left_join(ocean_rolling%>%select(marine_reserve, year, julian, temp_index, temp_c_mid, cdd_8, rolling8d, cdd_16, rolling16d))

#ocean_test<-ocean[ocean$depth_cat=="mid",]
#depths_test<-as.data.frame.matrix(table(ocean_test$ID,ocean_test$depth_m))
#depths_test$total<-rowSums(depths_test)
#table(depths_test$total) # so that excludes all data where we have multiple sampled depths

# COMBINING OCEANOGRAPHIC DATA WITH SETTLEMENT DATA 

str(settle)
#str(ocean_test) # note that I'm using the mid-water only depths test dataset
str(ocean_final) # Megan's datasets

# create the ID for the settlement dataset

settle$ID<-paste0(settle$year,"_",settle$julian,"_",settle$Site.Code)

# add the oceanographic fields I want to the settlement dataframe - use the ID to match
#settle_ocean<-merge(settle, ocean_test[,c(8:12)], by = c("ID"), all.x = TRUE)

settle_ocean<-merge(settle, ocean_full, by = c("ID","year","julian","Site.Code"), all.x = TRUE)

sum(is.na(settle_ocean['rolling8d'])) #149 missing observations with rolling windows
sum(is.na(settle_ocean['temp_c_mid'])) #153 missing with regular data
nrow(settle_ocean['rolling8d']) #1376 total - we dont lose much using rolling windows

positive <-settle_ocean%>%filter(SFLA_count>0)
sum(is.na(positive['rolling8d'])) #we lose 4 positive observations - not bad!
sum(is.na(positive['temp_c_mid'])) #we still lost 4
nrow(positive['rolling8d']) #255 total - we dont lose much using rolling windows


names(settle_ocean)
head(settle_ocean)
summary(settle_ocean)

# EXPORT DATA! WOOHOO! 
# though I may be back here shortly to see if I can match more oceanographic data...

write.csv(settle_ocean,"combined_settlement_ocean.csv",row.names = F)


##### Plots of oceanographic data ######

for(i in 1:2){
  p<-ggplot(ocean_full%>%filter(marine_reserve==marine_reserves[i]), 
            aes(x = julian, y = temp_index))+
    geom_line()+ 
    ggtitle(marine_reserves[i])+
    facet_wrap(~year)
  plot(p)
}

for(i in 1:2){
  p<-ggplot(ocean_full%>%filter(marine_reserve==marine_reserves[i]), 
            aes(x = julian, y = rolling16d))+
    geom_line()+ 
    ggtitle(marine_reserves[i])+
    facet_wrap(~year)
  plot(p)
} 

for(i in 1:2){
  p<-ggplot(ocean_full%>%filter(marine_reserve==marine_reserves[i], ), 
            aes(x = julian, y = cdd_16))+
    geom_line()+ 
    ggtitle(marine_reserves[i])+
    facet_wrap(~year)
  plot(p)
}
