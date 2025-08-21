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
yrmax<- 2025
recruitmentdevs<- data.frame(read.csv("data-yellowtail/RecruitmentDeviations2025draft.csv"))%>%
  filter(year>=yrmin & year <=yrmax)
Bio_param<-read.csv("data-yellowtail/BiologicalParameters2025.csv")#%>%
 # filter(Year>=yrmin & Year <=yrmax)
spawn_rec<-read.csv("data-yellowtail/spwn_rec_curve_non_SMURF.csv")%>%
  filter(Yr>=yrmin & Yr <=yrmax)
curve<- read.csv("data-yellowtail/spwn_rec_curve_non_SMURF_to0.csv")

SpawningBiomassTS <- ggplot(data=Bio_param, aes(y=Spawning.output..trillions.of.eggs., 
                                                x=Year))+
  geom_line()+ theme_classic() +
  labs(x="Year", y="Spawning Output \n (trillions of eggs)")+
  xlim(c(1995,2020))+
  #geom_hline(yintercept=0,lty=2)+
  ylim(c(0,11))+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
SpawningBiomassTS

Age0SpawningBiomass <- ggplot(data=curve, aes(y=Recruitment, x=SSB))+
  geom_line()+
  geom_point(data=Bio_param%>%filter(Year>1945), aes(y=Age.0.Recruits..1.000s., x=Spawning.output..trillions.of.eggs., colour=Year))+ 
  theme_classic() +
  labs(x="Spawning Output (trillions of eggs)", y="Age-0 (x1000)")+
  #geom_hline(yintercept=0,lty=2)+
  ylim(c(0,70000))+xlim(c(0,15))+
  scale_color_gradientn(colours =c("#08519C","#0f85a0","#dd4124"))+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
Age0SpawningBiomass

Age0TS <- ggplot(data=spawn_rec, aes(y=pred_recr, x=Yr))+
  geom_line()+ theme_classic() +
  ylim(c(0,75000))+
  labs(x="Year", y="Age-0 (x1000)")+
  #geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
Age0TS


recdevTS <-ggplot(data=recruitmentdevs%>%filter(year>=1996&year<=2025&Datatreatment=="2025 No_SMURF"), aes(y=Y_rec, x=year))+
  geom_line()+
 geom_ribbon(aes(ymin=Y_rec-sd,ymax=Y_rec+sd),alpha=0.2)+
  theme_classic() +
  labs(y="Recruitment Deviations", x="Year")+
  #geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
recdevTS
Figure1<-ggarrange(SpawningBiomassTS,Age0TS, Age0SpawningBiomass, nrow=3,
           labels = c('A', 'B', 'C', "D"))

png("Figures/Figure1.png",width=4,height=8,units="in",res=1200)
Figure1
dev.off()

pdf(file = "Figures/Figure1.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 8)
Figure1
dev.off()
