library(dplyr)
library(tidyr)
library(ggpubr)
library(corrplot)
library(mgcv)
library(DHARMa)
library(mgcViz)
library(gridExtra)
library(ROCR)
library(recdata)
library(predRecruit)
library(dplyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(corrplot)
library(car)
library(gratia)
library(ggpubr)

### Stock assessment data was pulled from Assessment on May 12 2025 (prior to STAR panel)
yrmin <- 1996
yrmax<- 2018
recruitmentdevs<- data.frame(read.csv("data-yellowtail/RecruitmentDeviations2025draft.csv"))%>%
  filter(year>=yrmin & year <=yrmax)
Bio_param<-read.csv("data-yellowtail/BiologicalParameters2025.csv")%>%
  filter(Year>=yrmin & Year <=yrmax)

SpawningBiomassTS <- ggplot(data=Bio_param, aes(y=Total.Biomass.4...mt., x=Year))+
  geom_line()+ theme_classic() +
  labs(x="Year", y="Spawning Biomass (mt)")+
  xlim(c(1995,2020))+
  #geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
SpawningBiomassTS

Age0SpawningBiomass <- ggplot(data=Bio_param, aes(y=Age.0.Recruits..1.000s., x=Total.Biomass.4...mt.))+
  geom_point()+ theme_classic() +
  labs(x="Spawning Biomass (mt)", y="Age-0 (x1000)")+
  xlim(c(70000,115000))+
  #geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
Age0SpawningBiomass

Age0TS <- ggplot(data=Bio_param, aes(y=Age.0.Recruits..1.000s., x=Year))+
  geom_line()+ theme_classic() +
  labs(x="Year", y="Age-0 (x1000)")+
  ylim(c(1995,2020))+
  #geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
Age0TS


recdevTS <-ggplot(data=recruitmentdevs%>%filter(year>=1996&year<=2019), aes(y=Y_rec, x=year))+
  geom_line()+
  geom_ribbon(aes(ymin=Y_rec-sd,ymax=Y_rec+sd),alpha=0.2)+
  theme_classic() +
  labs(y="Recruitment Deviations", x="Year")+
  ylim(c(1995,2020))+
  #geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

Figure1<-ggarrange(SpawningBiomassTS,Age0TS, Age0SpawningBiomass, recdevTS,
           labels = c('A', 'B', 'C', "D"))

png("Figures/Figure1.png",width=6,height=6,units="in",res=1200)
Figure1
dev.off()
