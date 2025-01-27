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
library(gridExtra)

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
lmlfo5.predict <- cbind(dat, predict(lmlfo5, interval = 'confidence'), residuals= residuals(lmlfo5))

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


#### Best Overall ####

best.gam<-gam(Y_rec~s(HCIlarv,k=3)+s(CutiSTIpjuv,k=3)+s(MLDlarv,k=3)+s(LSTlarv,k=3), data=dat)
best.gam.vif<-tableGrob(data.frame(VIF=round(vif(lm(Y_rec~HCIlarv+CutiSTIpjuv+MLDlarv+LSTlarv, data=dat)
),2)))
best.gam.predict <- cbind(dat, predict.gam(best.gam,se.fit=TRUE), residuals= residuals(best.gam))

best.gam.plot<- ggplot(best.gam.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2*se.fit, ymin=fit-2*se.fit), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("Best GAM Deviance Explained = ",round(summary(best.gam)$dev.expl,2)), subtitle=formula(best.gam),
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



best.lm<-lm(Y_rec~LSTlarv+I(LSTlarv^2)+CutiSTIpjuv+I(CutiSTIpjuv^2)+HCIlarv+MLDlarv, data=dat)
summary(best.lm)
best.lm.vif<-tableGrob(data.frame(VIF=round(vif(lm(Y_rec~LSTlarv+CutiSTIpjuv+MLDlarv+HCIlarv, data=dat)
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

cuti<- ggplot(dat, aes(CutiSTIpjuv, Y_rec)) + 
  geom_point() + 
  xlab('CutiSTIpjuv') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x + I(x^2), 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(CutiSTIpjuv,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()

lst <-ggplot(dat, aes(LSTlarv, Y_rec)) + 
  geom_point() + 
  xlab('LSTlarv') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x + I(x^2), 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(LSTlarv,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()


hci<-ggplot(dat, aes(HCIlarv, Y_rec)) + 
  geom_point() + 
  xlab('HCIlarv') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x, 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(HCIlarv,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()

mld<-ggplot(dat, aes(MLDpart, Y_rec)) + 
  geom_point() + 
  xlab('MLDpart') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x, 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(MLDpart,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()

ggarrange(hci,lst,cuti,mld, nrow=2, ncol=2)

ggsave("figures-yellowtail/partieleffects.png", width = 11, height =  8)       



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

#### Diagnostic Plots ####

##### Linear Model Diagnostics #####


m1=lm(Y_rec~LSTlarv+I(LSTlarv^2)+CutiSTIpjuv+I(CutiSTIpjuv^2)+HCIlarv+MLDpart, data=dat)
#m1=lmlfo5
(s1 = summary(m1))
capture.output(s1,file="R_Core.Model.txt")

png("figures-yellowtail/lmdiag.png",width=8,height=8,units="in",res=1200)
par(mfrow=c(2,2))
plot(m1)
dev.off()

##### Predict last 5 years ####

predict5 = c()
predict5$year = forecast_years
predict5 = data.frame(predict5)
predict5$pred = NA
head(predict5)

predictfit = c()
predictfit$year = short_years
predictfit = data.frame(predictfit)
predictfit$pred = NA
head(predictfit)

m1 = lm(as.formula(m1), data=dat[dat$year %in% short_years,])
NewData = data.frame(dat[dat$year %in% forecast_years,])
predict5[,'pred'] = predict(m1, newdata = NewData, interval = 'confidence')
predict5
write.table(predict5,paste0(results_dir, 'D_Predict_Last_5.csv'),col.names=TRUE, row.names=FALSE,sep=',')
s1 = summary(m1)

capture.output(s1, file = 'R_Predict_Last_5.txt')
saveRDS(s1,file='R_Predict_Last_5.rds')


predict5<-cbind(predict5, predict(m1, newdata = NewData, interval = 'confidence'))
predictfit<-cbind(predictfit, predict(m1, newdata = dat[dat$year %in% short_years,], interval = 'confidence'))
ggplot(predictfit, aes(x=year, y=fit)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              color='grey', alpha = 0.05) + 
  xlab("Year") + ylab('ln (rec devs)') + 
  geom_point(data=dat,aes(year,Y_rec)) + 
 #  geom_ribbon(data=dat,aes(x=year, y=Y_rec,ymin = Y_rec-sd, ymax = Y_rec+sd), 
 #             color='grey',fill='white', alpha = 0.05) +
  geom_line(data=predict5,aes(x=year, y=fit),col='red')+ 
  geom_point(data=predict5,aes(year,fit),col='red') + 
    geom_ribbon(data=predict5,aes(ymin = lwr, ymax = upr), 
              color='pink',fill='pink', alpha = 0.05)+
  theme_bw()

ggsave("figures-yellowtail/lm5pred.png",
       height=2.0, width = 3.5)

##### GAM Diag #####

m2 = gam(Y_rec~s(LSTlarv,k=3)+s(CutiSTIpjuv,k=3)+s(MLDlarv,k=3)+s(HCIlarv,k=3), data=dat)

png("figures-yellowtail/gamdiag.png",width=8,height=8,units="in",res=1200)
par(mfrow=c(2,2))
gam.check(m2)
dev.off()

dat$cooks.distance <- cooks.distance(m2)
dat$leverage <-influence.gam(m2)
dat$residuals.standard <- scale(residuals.gam(m2))
n <- length(residuals)  # Number of data points
influence  <- dat[(dat$cooks.distance> 4 / n |dat$cooks.distance> 1|dat$cooks.distance> 0.5),]


ggplot(dat, aes(x = leverage, y = residuals.standard)) +
  geom_hline(yintercept = 3, linetype = 'dotted') +
  geom_vline(xintercept = 2* mean(dat$leverage), linetype = 'dotted') +
  geom_point(alpha = 0.2) +
  ylab("Standardized Residuals")+
  xlab("Leverage")+
  geom_point(data = influence, aes(x = leverage, y = residuals.standard), color = 'darkorange') +
theme_minimal() +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5),  # Center plot title
    axis.title = element_text(size = 12),    # Axis title size
    axis.text = element_text(size = 10)      # Axis text size
  ) +
  geom_text(aes(label=ifelse(leverage > 4 / n, year, "")), cex=3, hjust=1.1)
ggsave("figures-yellowtail/gamleverage.png",
       height=2.0, width = 3.5)



predict5 = c()
predict5$year = forecast_years
predict5 = data.frame(predict5)
predict5$pred = NA
head(predict5)

predictfit = c()
predictfit$year = short_years
predictfit = data.frame(predictfit)
predictfit$pred = NA
head(predictfit)

m2 = gam(as.formula(m2), data=dat[dat$year %in% short_years,])
NewData = data.frame(dat[dat$year %in% forecast_years,])
predict5= data.frame(predict.gam(m2, newdata = NewData, se.fit=TRUE))
predict5<-cbind(NewData, predict5)
predictfit<-data.frame(predict.gam(m2, se.fit=TRUE, newdata = dat[dat$year %in% short_years,]))
predictfit<-cbind(predictfit,dat[dat$year %in% short_years,])
ggplot(predictfit, aes(x=year, y=fit)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = fit-2*se.fit, ymax = fit+2*se.fit), 
              color='grey', alpha = 0.05) + 
  xlab("Year") + ylab('ln (rec devs)') + 
  geom_point(data=dat,aes(year,Y_rec)) + 
  geom_line(data=predict5,aes(x=year, y=fit),col='red')+ 
  geom_point(data=predict5,aes(year,fit),col='red') + 
    geom_ribbon(data=predict5,aes(ymin = fit-2*se.fit, ymax = fit+2*se.fit), 
              color='pink',fill='pink', alpha = 0.05)+
  theme_bw()

ggsave("figures-yellowtail/gam5pred.png",
       height=2.0, width = 3.5)



png("figures-yellowtail/gampartial.png",width=8,height=8,units="in",res=1200)
par(mfrow=c(2,2))
plot.gam(m2, residuals=TRUE,se=TRUE)
dev.off()


#### An alternative best model ####
peplot <- function(mod,var,ci=.95, plot_points = "n",
                   xlab=var,ylab=names(mod[12]$model)[1],
                   main="Partial Effect Plot",
                   pe_lty=1,pe_lwd=3,pe_col="black",
                   ci_lty=1,ci_lwd=1,ci_col="black",
                   pch_col="black",pch_ty=19,
                   ylim=c(min(pred[,"lwr"]),max(pred[,"upr"]))){
  modDat <- mod[12]$model
  modDat1 <- modDat[,-1]
  modDat2 <- modDat[,which(names(modDat)!=var)]
  x <- resid(lm(modDat1[,var] ~., data=modDat1[,which(names(modDat1)!=var)]))
  y <- resid(lm(modDat2[,1] ~ ., modDat2[,-1]))
  plot(x,y,type=plot_points,xlab=xlab,ylab=ylab,
       ylim=ylim,col=pch_col,pch=pch_ty,
       main=main)
  part <- lm(y~x)
  wx <- par("usr")[1:2]
  new.x <- seq(wx[1],wx[2],len=100)
  pred <- predict(part, new=data.frame(x=new.x), interval="conf",
                  level = ci)
  lines(new.x,pred[,"fit"],lwd=pe_lwd,lty=pe_lty,col=pe_col)
  lines(new.x,pred[,"lwr"],lwd=ci_lwd,lty=ci_lty,col=ci_col)
  lines(new.x,pred[,"upr"],lwd=ci_lwd,lty=ci_lty,col=ci_col)
}

peplot(best.lm,var="LSTlarv",plot_points="p",ylim=NULL,pch_col="steelblue",ci_lty=2,
       xlab="Horsepower",ylab="Miles per Gallon")
peplot(best.lm,var="I(LSTlarv^2)",plot_points="p",ylim=NULL,pch_col="steelblue",ci_lty=2,
       xlab="Horsepower",ylab="Miles per Gallon")
peplot(best.lm,var="HCIlarv",plot_points="p",ylim=NULL,pch_col="steelblue",ci_lty=2,
       xlab="Horsepower",ylab="Miles per Gallon")


cuti<- ggplot(dat, aes(CutiSTIpjuv, Y_rec)) + 
  geom_point() + 
  xlab('CutiSTIpjuv') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x + I(x^2), 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(CutiSTIpjuv,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()

lst <-ggplot(dat, aes(LSTlarv, Y_rec)) + 
  geom_point() + 
  xlab('LSTlarv') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x + I(x^2), 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(LSTlarv,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()


hci<-ggplot(dat, aes(HCIlarv, Y_rec)) + 
  geom_point() + 
  xlab('HCIlarv') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x, 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(HCIlarv,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()

ggplot(dat, aes(DDlarv, Y_rec)) + 
  geom_point() + 
  xlab('DDlarv') + ylab("ln Rec Devs") + 
  geom_smooth(method = lm, formula=y ~ x, 
              color='black', fill='grey90') +
  ggrepel::geom_text_repel(data=dat,
    aes(DDlarv,Y_rec, label= year), color = "black", size = 3) +
  theme_classic()
