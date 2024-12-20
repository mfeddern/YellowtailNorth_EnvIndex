library(bayesdfa)
library(tidyverse)
library(ggplot2)
library(readr) # faster writing
library(data.table) # faster writing of large files
library(lubridate)
library(mgcv)
library(MuMIn)
library(corrplot)
library(car)
library(ggpubr)
library(fst)
library(viridis)
library(ggbeeswarm)
library(gratia)

#### Best LMs ####

##### LOO #####

lmloo<-lm(Y_rec~LSTpjuv+ONIpjuv+CutiSTIpjuv+I(CutiSTIpjuv^2), data=dat)
summary(lmloo)

lmloo.vif<-tableGrob(data.frame(VIF=round(vif(lmloo),2)))

lmloo.predict <- cbind(dat, predict(lmloo, interval = 'confidence'), residuals= residuals(lmloo))

lmloo.plot<- ggplot(lmloo.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=lwr, ymin=upr), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("LM, LOO: Rsquared = ",round(summary(lmloo)$r.squared,2)), subtitle=formula(lmloo),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
lmloo.plot

lmloo.res <- ggplot(lmloo.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggarrange(lmloo.plot, ggarrange(lmloo.res,lmloo.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/lmloo.png", height = 8, width = 10)
##### LFO5 #####

lmlfo5<-lm(Y_rec~DDben+DDpre+HCIpjuv, data=dat)
lmlfo5.vif<-tableGrob(data.frame(VIF=round(vif(lmlfo5),2)))
vif(lmlfo5)
lmlfo5.predict <- cbind(dat, predict(lmlfo5, interval = 'confidence'), residuals= residuals(lmlfo10))

lmlfo5.plot<- ggplot(lmlfo5.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=lwr, ymin=upr), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("LM, LFO5: Rsquared = ",round(summary(lmlfo5)$r.squared,2)), subtitle=formula(lmlfo5),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
lmlfo5.plot

lmlfo5.res <- ggplot(lmlfo5.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  labs(title="Best LM, LFO5", subtitle=formula(lmlfo5),x="Fitted", y="Residuals")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
lmlfo5.res

ggarrange(lmlfo5.plot, ggarrange(lmlfo5.res,lmlfo5.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/lmlfo5.png", height = 8, width = 10)
##### LFO10 #####

lmlfo10<-lm(Y_rec~MLDpjuv+CutiSTIpjuv+I(CutiSTIpjuv^2), data=dat)
summary(lmlfo10)
lmlfo10.vif<-tableGrob(data.frame(VIF=round(vif(lmlfo10),3)))
vif(lmlfo10)
lmlfo10.predict <- cbind(dat, predict(lmlfo10, interval = 'confidence'), residuals= residuals(lmlfo10))

lmlfo10.plot<- ggplot(lmlfo10.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=lwr, ymin=upr), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("LM, LFO10: Rsquared = ",round(summary(lmlfo10)$r.squared,2)),x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
lmlfo10.plot

lmlfo10.res <- ggplot(lmlfo10.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  labs(title="Best LM, LFO10", subtitle=formula(lmlfo10),x="Fitted", y="Residuals")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
lmlfo10.res
ggarrange(lmlfo10.plot, ggarrange(lmlfo10.res,lmlfo10.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/lmlfo10.png", height = 8, width = 10)


#### Best GAMS ####

##### LOO #####

gamloo<-gam(Y_rec~s(HCIlarv,k=3)+s(PDOpjuv,k=3)+s(LSTlarv,k=3), data=dat)
summary(gamloo)

gamloo.vif<-tableGrob(data.frame(VIF=round(vif(lm(Y_rec~HCIlarv+PDOpjuv+LSTlarv, data=dat)
),2)))

gamloo.predict <- cbind(dat, predict.gam(gamloo,se.fit=TRUE), residuals= residuals(gamloo))

gamloo.plot<- ggplot(gamloo.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2*se.fit, ymin=fit-2*se.fit), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("gam, LOO: Deviance Explained = ",round(summary(gamloo)$dev.expl,2)), subtitle=formula(gamloo),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gamloo.plot

gamloo.res <- ggplot(gamloo.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggarrange(gamloo.plot, ggarrange(gamloo.res,gamloo.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/gamloo.png", height = 8, width = 10)
##### LFO5 #####

gamlfo5<-gam(Y_rec~s(HCIpjuv,k=3)+s(Tpart,k=3)+s(LSTlarv,k=3), data=dat)
gamlfo5.vif<-tableGrob(data.frame(VIF=round(vif(lm(Y_rec~HCIpjuv+Tpart+LSTlarv, data=dat)
),2)))
gamlfo5.predict <- cbind(dat, predict.gam(gamlfo5,se.fit=TRUE), residuals= residuals(gamlfo5))

gamlfo5.plot<- ggplot(gamlfo5.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2*se.fit, ymin=fit-2*se.fit), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("gam, LFO5: Deviance Explained = ",round(summary(gamlfo5)$dev.expl,2)), subtitle=formula(gamlfo5),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gamlfo5.plot

gamlfo5.res <- ggplot(gamlfo5.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  labs(title="Best gam, LFO5", subtitle=formula(gamlfo5),x="Fitted", y="Residuals")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gamlfo5.res

ggarrange(gamlfo5.plot, ggarrange(gamlfo5.res,gamlfo5.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/gamlfo5.png", height = 8, width = 10)
##### LFO10 #####

gamlfo10<-gam(Y_rec~s(CutiSTIpjuv,k=3)+s(MLDpjuv,k=3), data=dat)
gamlfo10.vif<-tableGrob(data.frame(VIF=round(vif(lm(Y_rec~MLDpjuv+CutiSTIpjuv, data=dat)
),2)))
gamlfo10.predict <- cbind(dat, predict.gam(gamlfo10,se.fit=TRUE), residuals= residuals(gamlfo10))

gamlfo10.plot<- ggplot(gamlfo10.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2*se.fit, ymin=fit-2*se.fit), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("gam, LFO10: Deviance Explained = ",round(summary(gamlfo10)$dev.expl,2)), subtitle=formula(gamlfo10),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gamlfo10.plot


gamlfo10.res <- ggplot(gamlfo10.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  labs(title="Best gam, LFO10", subtitle=formula(gamlfo10),x="Fitted", y="Residuals")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gamlfo10.res
ggarrange(gamlfo10.plot, ggarrange(gamlfo10.res,gamlfo10.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/gamlfo10.png", height = 8, width = 10)


g2 <- ggplot(resid_df, aes(fitted,residuals)) +
  geom_point() +
  facet_wrap(~ var, scale="free_x") +
  theme_bw() +
  geom_smooth() +
  xlab("Fitted") + ylab("Residuals") +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"))
g2

g3 <- ggplot(resid_df, aes(year,fitted)) +
  geom_line() +
  geom_point(aes(y=Y_Rec))+
  facet_wrap(~ var, scale="free_x") +
  theme_bw() +
  xlab("Year") + ylab("Fitted") +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"))
g3


### Best ###

best.gam<-gam(Y_rec~s(ONIpre,k=3)+s(CutiSTIpjuv,k=3)+s(MLDpart,k=3)+s(LSTlarv,k=3), data=dat)
best.gam.vif<-tableGrob(data.frame(VIF=round(vif(lm(Y_rec~ONIpre+CutiSTIpjuv+MLDpart+LSTlarv, data=dat)
),2)))
best.gam.predict <- cbind(dat, predict.gam(best.gam,se.fit=TRUE), residuals= residuals(best.gam))

best.gam.plot<- ggplot(best.gam.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2*se.fit, ymin=fit-2*se.fit), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("gam, LFO10: Deviance Explained = ",round(summary(best.gam)$dev.expl,2)), subtitle=formula(best.gam),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
best.gam.plot


best.gam.res <- ggplot(best.gam.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  #geom_smooth() +
  labs(title="", subtitle=formula(best.gam),x="Fitted", y="Residuals")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
best.gam.res
ggarrange(best.gam.plot, ggarrange(best.gam.res,best.gam.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/gambest.png", height = 8, width = 10)



best.lm<-lm(Y_rec~ONIpre+CutiSTIpjuv+I(CutiSTIpjuv^2)+HCIlarv+MLDpart, data=dat)
summary(best.lm)
best.lm.vif<-tableGrob(data.frame(VIF=round(vif(lm(Y_rec~ONIpre+CutiSTIpjuv+MLDpart+HCIlarv, data=dat)
),2)))
best.lm.predict <-cbind(dat, predict(best.lm, interval = 'confidence'), residuals= residuals(best.lm))
best.lm.plot<- ggplot(best.lm.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=upr, ymin=lwr), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("LM R squared = ",round(summary(best.lm)$r.sq,2)), subtitle=formula(best.lm),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
best.lm.plot


best.lm.res <- ggplot(best.lm.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  #geom_smooth() +
  labs(title="", subtitle=formula(best.lm),x="Fitted", y="Residuals")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
best.lm.res
ggarrange(best.lm.plot, ggarrange(best.lm.res,best.lm.vif, ncol=2),nrow=2)+ bgcolor("White")

ggsave("figures-yellowtail/lmbest.png", height = 8, width = 10)



df.subset <- envir_data[, names(envir_data)[(names(envir_data) %in% unique(margtop5$cov))]]
corrplotsat<-cor(df.subset)
corrplot.mixed(corrplotsat)
colSums(corrplotsat)

corr2<-cor(envir_data%>%select(ONIpre,CutiSTIpjuv,MLDpart,LSTlarv))
best<-gam(Y_rec~s(ONIpre,k=3)+s(CutiSTIpjuv,k=3)+s(MLDpart,k=3)+s(LSTlarv,k=3), data=dat)
best_lm<-lm(Y_rec~ONIpre+CutiSTIpjuv+LSTlarv+MLDlarv, data=dat)
vif(best_lm)
summary(best_lm)
plot.gam(best)
draw(best, transform = "response", cex = 2)
