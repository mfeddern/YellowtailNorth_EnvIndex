library(VGAM)
########## BEST MODEL INFO #############################################

###############################
#### set best model here ######
###############################
#### Linear Model Diagnostics ####


m1=lm(Y_rec~ONIpre+CutiSTIpjuv+I(CutiSTIpjuv^2)+HCIlarv+MLDpart, data=dat)
(s1 = summary(m1))
capture.output(s1,file="R_Core.Model.txt")
par(mfrow=c(2,2))
plot(m1)

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
  geom_line(data=predict5,aes(x=year, y=fit))+ 
  geom_point(data=predict5,aes(year,fit),col='red') + 
    geom_ribbon(data=predict5,aes(ymin = lwr, ymax = upr), 
              color='pink', alpha = 0.05)+
  theme_bw()

ggsave("figures-yellowtail/lm5pred.png",
       height=2.0, width = 3.5)












m2 = gam(Y_rec~s(ONIpre,k=3)+s(CutiSTIpjuv,k=3)+s(MLDpart,k=3)+s(LSTlarv,k=3), data=dat)

gam.check(m2)


dat$cooks.distance <- cooks.distance(m2)
dat$leverage <-influence.gam(m2)
dat$residuals.standard <- scale(residuals.gam(m2))
n <- length(residuals)  # Number of data points
influence  <- dat[(dat$cooks.distance> 4 / n |dat$cooks.distance> 1|dat$cooks.distance> 0.5),]


ggplot(dat, aes(x = leverage, y = residuals.standard)) +
  geom_hline(yintercept = 3, linetype = 'dotted') +
  geom_vline(xintercept = 2* mean(dat$leverage), linetype = 'dotted') +
  geom_point(alpha = 0.2) +
  geom_point(data = influence, aes(x = leverage, y = residuals.standard), color = 'orange') +
theme_minimal() +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5),  # Center plot title
    axis.title = element_text(size = 12),    # Axis title size
    axis.text = element_text(size = 10)      # Axis text size
  ) +
  geom_text(aes(label=ifelse(residuals.standard > 2, model, "")), cex=3, hjust=1.1)




