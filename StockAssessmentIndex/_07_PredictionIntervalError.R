library(dplyr)
library(ggplot2)
library(Metrics)

data_years = 1993:2019
results_dir = paste0(main_dir, "/results_",data_years[1],'-',max(data_years),'_predict_',max(df$year),'/')
fig_dir = paste0(main_dir, "/figures_",data_years[1],'-',max(data_years),'_predict_',max(df$year),'/')
mod = dir(results_dir)
mod = dir(results_dir)

mod = mod[grep('Environmental-index-4Assessement', mod)][1]
mod
# load data
Index<-read.csv(paste0(results_dir,mod)) #reading in data that Nick provided
m1 = readRDS(paste0(results_dir, "Best_fit_model_results.rds"))
s1 = summary(m1)
# residualerror <- 1.032 #residual error from Nick's model output from "Best-fit-summary" file
residualerror <- s1$sigma
S<- residualerror

# assign the critical value (only necessary for checking results against prediction interval), 
# this is based on DF from Best-fit-summary
# critical_value <- 2.101 
critical_value = stats::qt(p = 0.975, df = s1$df[2])
critical_value


IndexV2<-Index%>%
  # adding in the standard deviation of the residual error to the SE
  mutate(prediction_error = sqrt(S^2+se^2))%>% 
  # calculating an interval check
  mutate(int_lwr_check = fit-critical_value*prediction_error, 
         int_upr_check = fit+critical_value*prediction_error)%>% 
  # back calculating from prediction interval 
  # to check that the intervals line up as a triple check
  mutate(prederror1 = (fit-pred_int_lwr)/critical_value, 
         prederror2 = (fit-pred_int_upr)/-critical_value) %>% 
  # add in confidence interval just to plot
  mutate(conf_upr=fit+critical_value*se, 
         conf_lwr=fit-critical_value*se) 

# just checking if these are the same from back calculating from the prediction interval

IndexV2 %>% select(prederror1, prederror2, prediction_error) 

plot <- ggplot(data=IndexV2, aes(x=year, y=fit))+
  geom_line()+
  #prediction interval from the prediction error calculated from model
  geom_ribbon(aes(ymin=pred_int_lwr, ymax=pred_int_upr), alpha = 0.2, fill = 'blue')+ 
  #prediction interval from predict()
  geom_ribbon(aes(ymin=int_lwr_check, ymax=int_upr_check), alpha = 0.2, fill = 'red', color='black')+
  #confidence interval
  geom_ribbon(aes(ymin=conf_lwr, ymax=conf_upr), alpha = 0.2, fill = 'green')+
  theme_bw()
plot
png(filename = paste0(fig_dir,"PredictionIntervalsSablefish-",data_years[1],'-',max(data_years),".png"),
    width = 700, height = 500)
plot
dev.off()

# fix name
mod2 = stringr::str_remove(mod, ".csv")
# select columns
IndexV3<-IndexV2%>%select(year, fit, pred_int_lwr, pred_int_upr, se, Y_rec, prediction_error)
# write out index
write.csv(IndexV3,paste0(results_dir, mod2,"_wPredictionError.csv"))
