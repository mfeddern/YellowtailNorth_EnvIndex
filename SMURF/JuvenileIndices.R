library(MARSS)
library(ggplot2)
library(forecast)
library(dplyr)
library(lubridate)
library(mgcv)
library(tidyr)
library(corrplot)
library(reshape2)

OCNMS<-read.csv('JuvenileIndexDatasets/Estimated-YOY-trend-coast.csv')%>%
  rename(Index=grand.mean)%>%
  select(year, Index)%>%
  mutate(dataset="OCNMS Dive Survey")

RREAS_North<-read.csv('JuvenileIndexDatasets/YOYGroundfishCPUEPerYr_North.csv')%>%
  rename(year=YEAR, Index=mean_cpue)%>%
  select(year, Index)%>%
  mutate(dataset="Northern YOY Rockfish")

Oceanographic<-read.csv('JuvenileIndexDatasets/OceanographicIndexV1.csv')%>%
  rename(Index=fit)%>%
  select(year, Index)%>%
  mutate(dataset="Oceanographic Index")

SMURF<-read.csv('JuvenileIndexDatasets/index_forSS.csv')%>%
  rename(Index=obs)%>%
  select(year, Index)%>%
  mutate(dataset="SMURF Juvenile Survey")

RREAS_Yellowtail_Coast<-read.csv('JuvenileIndexDatasets/ytail_coastwide_indices.csv')%>%
  rename(Index=est, year=YEAR)%>%
  select(year, Index)%>%
  mutate(dataset="Yellowtail RREAS Coastwide")

dat_long<-bind_rows(OCNMS,RREAS_North)%>%
  bind_rows(Oceanographic)%>%
  bind_rows(SMURF)%>%
  bind_rows(RREAS_Yellowtail_Coast)

dat_long<-dat_long%>%filter(year>2013)%>%
  group_by(dataset)%>%
  mutate(scaled_index=scale(Index))%>%
  rename(Dataset=dataset)

YOYIndex<-ggplot(data=dat_long%>%filter(year>2013)%>%
         group_by(Dataset)%>%
         mutate(scaled_index=scale(Index)),
       aes(group=Dataset,col=Dataset,x=year,y=scaled_index))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0,lty=2)+
  ylab("Standardized Index")+
  xlab("Year")+
  theme_classic()

YOYno2016<-ggplot(data=dat_long%>%filter(year>2016)%>%
         group_by(Dataset)%>%
         mutate(scaled_index=scale(Index)),
       aes(group=Dataset,col=Dataset,x=year,y=scaled_index))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0,lty=2)+
  ylab("Standardized Index (omitting 2016)")+
  xlab("Year")+
  theme_classic()

png("YOYindices.png",width=6,height=3,units="in",res=1200)
YOYIndex
dev.off()

png("YOYindicesNo2016.png",width=6,height=3,units="in",res=1200)
YOYno2016
dev.off()
