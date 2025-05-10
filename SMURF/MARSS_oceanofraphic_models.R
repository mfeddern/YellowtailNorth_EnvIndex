library(MARSS)
library(ggplot2)
library(forecast)
library(dplyr)
library(lubridate)
library(mgcv)
library(tidyr)
library(corrplot)
library(reshape2)
autoplot(resids, plot.type = "acf")
#### Reading in Data ####
data <-read.csv('Oceanography_DailyMean.csv')
head(data)
unique(data$site)
marine_reserves<-unique(data$marine_reserve)
unique(data$location)

####Herring####

yy=data%>%
  #filter(marine_reserve=='Otter Rock'| marine_reserve=='Redfish Rocks')%>% #filtering for reserves that we have SMURF data
  mutate(depth_bin=ifelse(depth_m==1,'surface',ifelse(depth_m==5|depth_m==7.5,"mid","bottom")))%>% # binning 5 and 7.5 depth to a single "mid" depth
 # group_by(depth_bin,date,julian,year,marine_reserve)%>%
  #summarise(temp=mean(temp_c))%>%
  ungroup()

yy$date<-as.Date(yy$date)


ggplot(yy, aes(x = date, y = temp_c, group=as.factor(depth_bin), col=as.factor(depth_bin)))+
  geom_line()+ 
  facet_wrap(~marine_reserve)

for(i in 1:4){
p<-ggplot(yy%>%filter(marine_reserve==marine_reserves[i]), aes(x = julian, y = temp_c, group=as.factor(depth_bin), col=as.factor(depth_bin)))+
    geom_line()+ 
    ggtitle(marine_reserves[i])+
    facet_wrap(~year)
plot(p)

}
ggsave("Figures/CascadeHeadTS.png", height = 8, width = 8)
ggsave("Figures/RedfishRocksTS.png", height = 8, width = 8)
ggsave("Figures/OtterRockTS.png", height = 8, width = 8)
for(i in 1:length(marine_reserves)){
p<-ggplot(yy%>%filter(marine_reserve==marine_reserves[i]), aes(x = julian, y = oxy_mll, group=as.factor(depth_bin), col=as.factor(depth_bin)))+
    geom_line()+ 
    ggtitle(marine_reserves[i])+
    facet_wrap(~year)
plot(p)

}
yy<-yy%>%filter(marine_reserve=="Otter Rock"|
                marine_reserve=="Redfish Rocks")
#### correlation across location ####


cor_dat_loc<-yy%>%filter(depth_bin=="mid")%>%
  select(year,julian,location, temp_c)%>%
  pivot_wider(names_from=location, values_from=temp_c)
locations<-unique(yy$location)
df.subset <- cor_dat_loc[, names(cor_dat_loc)[(names(cor_dat_loc) %in% locations)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

cor_dat_loc<-yy%>%filter(depth_bin=="surface")%>%
  select(year,julian,location, temp_c)%>%
  pivot_wider(names_from=location, values_from=temp_c)
locations<-unique(yy$location)
df.subset <- cor_dat_loc[, names(cor_dat_loc)[(names(cor_dat_loc) %in% locations)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

cor_dat_loc<-yy%>%filter(depth_bin=="bottom")%>%
  select(year,julian,location, temp_c)%>%
  pivot_wider(names_from=location, values_from=temp_c)
locations<-unique(yy$location)
df.subset <- cor_dat_loc[, names(cor_dat_loc)[(names(cor_dat_loc) %in% locations)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

cor_dat_depth<-yy%>%filter(location=="Otter Rock")%>%
  select(year,julian,depth_bin, temp_c)%>%
  pivot_wider(names_from=depth_bin, values_from=temp_c)
depth<-unique(yy$depth_bin)
df.subset <- cor_dat_depth[, names(cor_dat_depth)[(names(cor_dat_depth) %in% depth)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

cor_dat_depth<-yy%>%filter(location=="Redfish Rocks")%>%
  select(year,julian,depth_bin, temp_c)%>%
  pivot_wider(names_from=depth_bin, values_from=temp_c)
depth<-unique(yy$depth_bin)
df.subset <- cor_dat_depth[, names(cor_dat_depth)[(names(cor_dat_depth) %in% depth)]]
corrplotsat<-cor(df.subset, use = "complete.obs")
corr_loc<-corrplot.mixed(corrplotsat)

##### Otter Rock MARSS #####
full_days<-data.frame(date=seq(as.Date("2010/5/1"), as.Date("2024/10/15"), "days"))%>%
  mutate(julian=yday(date))
yyy<- yy%>%filter(marine_reserve=="Otter Rock"&julian>100&julian<275)
merge(yyy%>%select(-julian),full_days)

dates<-seq(as.Date(min(yy$date)), as.Date(max(yy$date)), "days")
years<-lubridate::year(dates)
months<-lubridate::month(dates)
depth_bins<-c('bottom', 'mid','surface')
covariates_df<-data.frame(date=as.Date(rep(dates,3)),
           depth_bin=rep(depth_bins,each=length(dates)),
           year=as.numeric(years), months=months, julian=yday(dates)
)

df<-data.frame(yyy)%>%
  mutate(year2=ifelse(year==2010,1,
         ifelse(year==2011,2,
         ifelse(year==2012,3,
         ifelse(year==2013,4,
         ifelse(year==2014,5,
         ifelse(year==2015,6,
        ifelse(year==2016,7,
         ifelse(year==2017,8,
         ifelse(year==2018,9,
         ifelse(year==2019,10,
        # ifelse(year==2020,11,
         ifelse(year==2021,11,
        ifelse(year==2022,12,
        ifelse(year==2023,13,
         ifelse(year==2024,14,0)))))))))))))))%>%
  select(date, year, year2)

cov <- reshape2::acast(na.omit(df), year ~ date,fun.aggregate=min, value.var = "year2")
cov[!is.finite(cov)] <- 0
ncol(cov)

years<-seq(2011,2024,1)
C3<-matrix(as.character(years), 3, 14, byrow = TRUE)
C1<-matrix(as.character(years), 1, 14, byrow = TRUE)

dim(data.frame(C))
dat <- reshape2::acast(yyy%>%filter(marine_reserve=="Otter Rock"&julian>100&julian<275),fun.aggregate=mean, depth_bin ~ date, value.var = "temp_c")
#the.mean <- apply(dat, 1, mean, na.rm = TRUE)
#the.sigma <- sqrt(apply(dat, 1, var, na.rm = TRUE))
#dat <- (dat - the.mean) * (1/the.sigma)


# Each taxon has unique density-dependence
B <- "equal"
# Assume independent process errors
Q <- "diagonal and unequal"
# We have demeaned the data & are fitting a mean-reverting
# model by estimating a diagonal B, thus
U <- "zero"
# Each obs time series is associated with only one process
Z <- "identity"
# The data are demeaned & fluctuate around a mean
A <- "zero"
# equation
D <- "zero"
d <- "zero"
B = "identity" 
# U = "identity" 
Q = "equalvarcov"
Z = "identity" #matrix(1, 3, 1)
R = "zero"
x0 = "unequal"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, 
    C = C3, c = cov, D = D, d = d)
seas.mod.1 <- MARSS(dat, model = model.list, control = list(maxit = 1500))

coef(seas.mod.1, type = "matrix")$C

states<-seas.mod.1$states
row.names(states)<-depth_bins
colnames(states)<-colnames(dat)
states.t<-data.frame(t(states))
states.t<-states.t%>%
  mutate(dates=as.Date(colnames(dat)))%>%
  pivot_longer(-dates,names_to='depth_bin', values_to = 'temp')%>%
  mutate(julian=yday(dates),Year=lubridate::year(dates))

p<-ggplot(states.t, aes(x = julian, y = temp, group=as.factor(depth_bin), col=as.factor(depth_bin)))+
    geom_line()+ 
    ggtitle("Otter Rock")+
    facet_wrap(~Year)
p

ggsave("Figures/OtterRockFit.png", height = 8, width = 8)


##### Otter Rock MARSS #####

yyy<- yy%>%filter(marine_reserve=="Redfish Rocks"&julian>100&julian<275)
merge(yyy%>%select(-julian),full_days)

dates<-seq(as.Date(min(yy$date)), as.Date(max(yy$date)), "days")
years<-lubridate::year(dates)
months<-lubridate::month(dates)
depth_bins<-c('bottom', 'mid','surface')
covariates_df<-data.frame(date=as.Date(rep(dates,3)),
           depth_bin=rep(depth_bins,each=length(dates)),
           year=as.numeric(years), months=months, julian=yday(dates)
)

df<-data.frame(yyy)%>%
  mutate(year2=ifelse(year==2010,1,
         ifelse(year==2011,2,
         ifelse(year==2012,3,
         ifelse(year==2013,4,
         ifelse(year==2014,5,
         ifelse(year==2015,6,
        ifelse(year==2016,7,
         ifelse(year==2017,8,
         ifelse(year==2018,9,
         ifelse(year==2019,10,
        # ifelse(year==2020,11,
         ifelse(year==2021,11,
        ifelse(year==2022,12,
        ifelse(year==2023,13,
         ifelse(year==2024,14,0)))))))))))))))%>%
  select(date, year, year2)

cov <- reshape2::acast(na.omit(df), year ~ date,fun.aggregate=min, value.var = "year2")
cov[!is.finite(cov)] <- 0
ncol(cov)

years<-c(seq(2010,2019,1),seq(2021,2024,1))
C3<-matrix(as.character(years), 3, 14, byrow = TRUE)
C1<-matrix(as.character(years), 1, 14, byrow = TRUE)

dim(data.frame(C))
dat <- reshape2::acast(yyy%>%filter(marine_reserve=="Redfish Rocks"&julian>100&julian<275),fun.aggregate=mean, depth_bin ~ date, value.var = "temp_c")
#the.mean <- apply(dat, 1, mean, na.rm = TRUE)
#the.sigma <- sqrt(apply(dat, 1, var, na.rm = TRUE))
#dat <- (dat - the.mean) * (1/the.sigma)


# We are not including covariate effects in the obs
# equation
D <- "zero"
d <- "zero"

B = "identity" 
U = "zero"
Q = "equalvarcov"
Z = "identity" #matrix(1, 3, 1)
A = "scaling"
R = "equalvarcov"
x0 = "unequal"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, 
    C = C3, c = cov, D = D, d = d)
seas.mod.1 <- MARSS(dat, model = model.list, control = list(maxit = 1500))

coef(seas.mod.1, type = "matrix")$C

states<-seas.mod.1$states
row.names(states)<-depth_bins
colnames(states)<-colnames(dat)
states.t<-data.frame(t(states))
states.t<-states.t%>%
  mutate(dates=as.Date(colnames(dat)))%>%
  pivot_longer(-dates,names_to='depth_bin', values_to = 'temp')%>%
  mutate(julian=yday(dates),Year=lubridate::year(dates))

p<-ggplot(states.t, aes(x = julian, y = temp, group=as.factor(depth_bin), col=as.factor(depth_bin)))+
    geom_line()+ 
    ggtitle("Redfish Rocks")+
    facet_wrap(~Year)
p

ggsave("Figures/RedfishRocksFit.png", height = 12, width = 12)
