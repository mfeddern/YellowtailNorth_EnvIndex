library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(plyr)
library(corrplot)

chl_index <- readRDS('data-yellowtail/Raw_Environmental/chl_indices_WCBTS1996_2023.rds')%>%
  filter(region=='north')%>%
  filter(month==4|month==5|month==6|month==7|month==8)%>%
  group_by(year)%>%
  dplyr::summarise(CHL=mean(chl))%>%
  mutate(CHLpjuv=(CHL-mean(CHL))/sd(CHL))%>%
  select(year, CHLpjuv)

#chl_index <- readRDS('data-yellowtail/Raw_Environmental/chl_indices_WCBTS2025.rds')%>%
#  filter(region=='north')%>%
#  filter(month==4|month==5|month==6|month==7|month==8)%>%
#  group_by(year)%>%
#  dplyr::summarise(CHL=mean(chl))%>%
#  mutate(CHLpjuv=(CHL-mean(CHL))/sd(CHL))%>%
#  select(year, CHLpjuv)

pp_index <- readRDS('data-yellowtail/Raw_Environmental/pp_indices_WCBTS1996_2023.rds')%>%
  filter(region=='north')%>%
  filter(month==4|month==5|month==6|month==7|month==8)%>%
  group_by(year)%>%
  dplyr::summarise(PP=mean(pp))%>%
  mutate(PPpjuv=(PP-mean(PP))/sd(PP))%>%
  select(year, PPpjuv)

#pp_index <- readRDS('data-yellowtail/Raw_Environmental/pp_indices_WCBTS2025.rds')%>%
#  filter(region=='north')%>%
#  filter(month==4|month==5|month==6|month==7|month==8)%>%
#  group_by(year)%>%
#  dplyr::summarise(PP=mean(pp))%>%
#  mutate(PPpjuv=(PP-mean(PP))/sd(PP))%>%
#  select(year, PPpjuv)


Beuti_STI<- read.csv('data-yellowtail/Raw_Environmental/beuti_sti.csv')%>%
  select(X, X45N, X46N, X47N)%>%
  pivot_longer(cols=-X,names_to='location', values_to='BeutiSTI')%>%
  group_by(X)%>%
  dplyr::summarise(BeutiSTI=mean(BeutiSTI))%>%
  mutate(BeutiSTIpjuv=(BeutiSTI-mean(BeutiSTI))/sd(BeutiSTI))%>%
  dplyr::rename(year=X)%>%
  select(year, BeutiSTIpjuv)


BEUTI_TUMI<-read.csv('data-yellowtail/Raw_Environmental/beuti_tumi.csv')%>%
  select(X, X45N, X46N, X47N)%>%
  pivot_longer(cols=-X,names_to='location', values_to='BeutiTUMI')%>%
  group_by(X)%>%
  dplyr::summarise(BeutiTUMI=mean(BeutiTUMI))%>%
  mutate(BeutiTUMIpjuv=(BeutiTUMI-mean(BeutiTUMI))/sd(BeutiTUMI))%>%
  dplyr::rename(year=X)%>%
  select(year, BeutiTUMIpjuv)

CUTI_STI<- read.csv('data-yellowtail/Raw_Environmental/cuti_sti.csv')%>%
  select(X, X45N, X46N, X47N)%>%
  pivot_longer(cols=-X,names_to='location', values_to='CutiSTI')%>%
  group_by(X)%>%
  dplyr::summarise(CutiSTI=mean(CutiSTI))%>%
  mutate(CutiSTIpjuv=(CutiSTI-mean(CutiSTI))/sd(CutiSTI))%>%
  dplyr::rename(year=X)%>%
  select(year, CutiSTIpjuv)


CUTI_TUMI<-read.csv('data-yellowtail/Raw_Environmental/cuti_tumi.csv')%>%
  select(X, X45, X46, X47)%>%
  pivot_longer(cols=-X,names_to='location', values_to='CutiTUMI')%>%
  group_by(X)%>%
  dplyr::summarise(CutiTUMI=mean(CutiTUMI))%>%
  mutate(CutiTUMIpjuv=(CutiTUMI-mean(CutiTUMI))/sd(CutiTUMI))%>%
  dplyr:: rename(year=X)%>%
  select(year, CutiTUMIpjuv)

bakun_sti<- read.csv('data-yellowtail/Raw_Environmental/bakun_sti.csv')%>%
  select(-time)%>%
  filter(lat==45)%>%
  #pivot_wider(names_from = lat,values_from = bakun_sti)%>%
  group_by(year)%>%
  dplyr::summarise(bakun_sti=mean(bakun_sti))%>%
  mutate(bakun_sti=scale(bakun_sti))

Env_Indices<-join_all(list(CUTI_STI, CUTI_TUMI,BEUTI_TUMI, Beuti_STI,chl_index,pp_index,bakun_sti), by = 'year')

#write.csv(Env_Indices, "data-yellowtail/Processed_Environmental/2024updateGLORYs.csv")

write.csv(Env_Indices, "data-yellowtail/Processed_Environmental/2024update_Env_Indices.csv")

glorys<-read.csv('Data/Processed_Environmental/glorys-data-annual-yellowtail_subset.csv')
glorys_full<-glorys%>%
  left_join(Env_Indices)
corplotdat<-cor(na.omit(scale(glorys_full%>%select(-year))))
corrplot(corplotdat, method='number') # colorful number  

pdf(file = "Corr.Plot.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 18) # The height of the plot in inches

corrplot.mixed(corplotdat, lower = 'circle', upper = 'number', order = 'hclust')
dev.off()


CUTI_corr<- read.csv('data-yellowtail/Raw_Environmental/cuti_sti.csv') %>%
select(-X)
corrplot.mixed(cor(na.omit(CUTI_corr)), lower = 'circle', upper = 'number')


combine_data<-read.csv("data-yellowtail/02_DATA_Combined_glorys_yellowtail.csv")%>%
  left_join(bakun_sti)
write.csv(combine_data,"data-yellowtail/02_DATA_Combined_glorys_yellowtail.csv")
write.csv(CUTI_STI%>%left_join(bakun_sti),"data-yellowtail/upwellingts.csv")


cuti<-read.csv('data-yellowtail/Raw_Environmental/cuti_sti.csv')%>%
  select(X, X45N, X46N, X47N)%>%
  pivot_longer(cols=-X,names_to='location', values_to='CutiSTI')%>%
  group_by(X)%>%
  dplyr::summarise(CutiSTI=mean(CutiSTI))
mean(cuti$CutiSTI)
