library(tidyverse)
library(ggplot2)
library(readr) # faster writing
library(data.table) # faster writing of large files
library(lubridate)
library(mgcv)
library(MuMIn)

source("Src/_00_yellowtail-header.r")
source("Src/Functions-for-envir-index.r")
# set directories ##############################################################

# bring in data ################################################################
# combined fish and environmental drivers file
df = data.frame(read.csv(paste0(data_dir,"02_DATA_Combined_glorys_yellowtail.csv"), header = T))

# get predictors to create model formula #######################################
envir_data = df %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('year','sd','Y_rec','ZOOpjuv','ZOOben')))%>%  # drop terms not in model statement
    dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI','ONIlarv')))%>% 
    mutate_all(~ scale(.))
head(envir_data)
dim(envir_data)
data_years = 1994:2014

envir_data2<- envir_data[complete.cases(envir_data), ]
data.pca <- princomp(envir_data2)
summary(data.pca)

# build model - full formula ###################################################
quadratic_terms = c('LSTlarv','ONIpre','ONIlarv', 'CutiSTIpjuv','bakun_sti')

form_dredge = make_dredge_equation(envir_data = envir_data, quadratic_vars = quadratic_terms)
form_dredge
# save out for model testing later
saveRDS(form_dredge, paste0(results_dir,"formula_for_dredge.rds"))

# Final data preparation #######################################################
# remove ZOO terms because they miss two years; model can't run with NAs

data_1 = df %>% 
  dplyr::select(!any_of( c('ZOOpjuv', 'ZOOben')))  %>% 
  filter(year %in% data_years)

######## RUN DREDGE ##############################################################

fit = lm(form_dredge, data=data_1, na.action = na.fail)

# dredge_text
# write as a function to use in model testing later
fit_dredge <- function(fit){
  fit = fit
  mtable = dredge(fit, rank = AICc, m.lim = c(NA,3),
                  subset= # block correlated terms from the model
                    !(CutiSTIpjuv && BeutiSTIpjuv) && 
                    #!(DDben && ZOOben) && 
                    !(DDlarv && DDpjuv) && 
                    !(DDlarv && HCIlarv) && 
                    #!(DDlarv && ONIlarv) && 
                    !(DDlarv && PDOlarv) && 
                    !(DDlarv && Tpart) && 
                    #!(DDlarv && ZOOpjuv) && 
                    !(DDpjuv && HCIlarv) && 
                    !(DDpjuv && HCIpjuv) && 
                    !(DDpjuv && PDOlarv) && 
                    !(DDpjuv && PDOpjuv) && 
                    !(DDpjuv && Tpart) && 
                    #!(DDpjuv && ZOOpjuv) && 
                    !(DDpre && DDegg) && 
                    !(DDpre && DDlarv) && 
                    !(DDpre && HCIlarv) &&
                    !(DDpre && PDOlarv) && 
                    !(DDpre && Tcop) && 
                    !(DDpre && Tpart) && 
                    #!(DDpre && ZOOpjuv) && 
                    !(HCIlarv && PDOlarv) && 
                    !(HCIpjuv && PDOpjuv) && 
                    !(MLDpart && MLDlarv) && 
                    #!(ONIlarv && PDOlarv) && 
                   # !(ONIpre && ONIlarv) && 
                    #!(PDOlarv && ZOOpjuv) && 
                    !(Tpart && HCIlarv) && 
                    !(Tpart && PDOlarv) && 
                    #!(Tpart && ZOOpjuv) &&
                    dc(LSTlarv, I(LSTlarv^2)) &&
                    dc(ONIpre, I(ONIpre^2)) &&
                    dc(bakun_sti, I(bakun_sti^2)) &&
                   # dc(ONIlarv, I(ONIlarv^2)) &&
                    dc(CutiSTIpjuv, I(CutiSTIpjuv^2)),
                  extra = list(R2 = function(x)
                    summary(x)$r.squared[[1]], 
                    F = function(y)
                      summary(y)$fstatistic[[1]]),
                  trace=2 )
  
  return(mtable)
}

# output the dredge function to use in model testing
saveRDS(fit_dredge, paste0(results_dir,"fit_dredge.rds"))

# run dredge and fit models using function ####################################
mtable = fit_dredge(fit)

mtable4 = subset(mtable, delta<4)
print(mtable4, abbrev.names = FALSE)

saveRDS( mtable, file=paste0(results_dir,'Table_dredge-model-fits.rds'))
fwrite(mtable,paste0(results_dir,"Table_dredge-aic-delta4.csv"))
mtable = readRDS(paste0(results_dir,'Table_dredge-model-fits.rds') )

# Refit best fit model #########################################################

best_fit = find_best_fit(mtable, fewest_params=TRUE)
bf_mod = bf_mod_equation(best_fit)
bf_mod
saveRDS(bf_mod, paste0(results_dir, "Best_fit_model.rds")) # save out to use later

best_fit = lm( bf_mod, data=data_1)
summary(best_fit)
capture.output( summary(best_fit), file = paste0(results_dir, "Best-fit-summary.txt"))
capture.output( anova(best_fit), file = paste0(results_dir, "Best-fit-anova.txt"))

# predict model ################################################################

df_index = predict_best_fit(old_data=data_1, new_data=df, bf_mod=bf_mod)
fwrite(df_index, paste0(results_dir, "Environmental-index.csv"))

# plot predicted and rec devs ##################################################

# df_index = na.omit(df_index)

ggplot(df_index, aes(x=year, y=fit)) + 
  geom_ribbon(data=df_index, aes(ymin = fit-1.96*se, ymax = fit+1.96*se), 
              color='grey', alpha = 0.05) + 
  geom_line() + 
  ggtitle("CutiSTIpjuv + I(CutiSTIpjuv^2)")+
  xlab("Year") + ylab('ln (rec devs)') + 
  geom_point(aes(year,Y_rec)) + 
  scale_x_continuous(breaks=seq(1995,2025,5) , minor_breaks = 1993:max(df$year)) + 
  annotate(geom='segment', x=1993, xend=max(df_index$year), y=0, linetype='dotted') +
  theme_bw()

ggsave(paste0(fig_dir, "Predicted-recruitment-index.png"),
       height=2.0, width = 3.5)


mtable


CutiSTIpjuvPRED<-data.frame(CutiSTIpjuv=seq(-2,2,0.1))
mod_bf<- lm(Y_rec ~ CutiSTIpjuv + I(CutiSTIpjuv^2), data_1)
df_index = data.frame(y_rec_pred=predict(mod_bf,newdata =CutiSTIpjuvPRED, interval = 'confidence'), CutiSTIpjuvPRED)

ggplot(df_index, aes(x=CutiSTIpjuv, y=y_rec_pred.fit)) + 
  geom_ribbon(data=df_index, aes(ymin = y_rec_pred.lwr, ymax =y_rec_pred.upr), 
              color='grey', alpha = 0.05) + 
  geom_line() + 
  ggtitle("CutiSTIpjuv + I(CutiSTIpjuv^2)")+
  xlab("Year") + ylab('ln (rec devs)') + 
  geom_point() + 
  theme_bw()
# Refit best fit model #########################################################

best_fit = find_best_fit(mtable, fewest_params=TRUE)
bf_mod = bf_mod_equation(best_fit)
bf_mod
saveRDS(bf_mod, paste0(results_dir, "Best_fit_model.rds")) # save out to use later

best_fit = lm( bf_mod, data=data_1)
summary(best_fit)
capture.output( summary(best_fit), file = paste0(results_dir, "Best-fit-summary.txt"))
capture.output( anova(best_fit), file = paste0(results_dir, "Best-fit-anova.txt"))

# predict model ################################################################

df_index = predict_best_fit(old_data=data_1, new_data=df, bf_mod=bf_mod)
fwrite(df_index, paste0(results_dir, "Environmental-index.csv"))

# plot predicted and rec devs ##################################################

# df_index = na.omit(df_index)

mtable[is.na(mtable)] <- 0
param_names<-colnames(mtable)[2:31]
AICcTab=NA
for(i in 2:31){
  AICweight =NA
  for(j in 1:nrow(mtable)){
  AICweight = rbind(AICweight,ifelse(mtable[j,i]!= 0, mtable[j,38],0))
  }
  AICcTab=cbind(AICcSum,AICweight)
}

AICSum<-AICcTab[2:nrow(AICcTab),2:31]
colnames(AICSum)<-param_names
AICSum<-data.frame(colSums(AICSum))







mtable2<-mtable
sum(mtable$weight)
mtable2[is.na(mtable2)] <- 0
param_names<-colnames(mtable2)[31:2]
AICcTab=NA
for(i in 2:31){
  AICweight =NA
  for(j in 1:nrow(mtable2)){
    AICweight = rbind(AICweight,ifelse(mtable2[j,i]!= 0, mtable2[j,38],0))
  }
  AICcTab=cbind(AICweight,AICcTab)
}


AICSum<-data.frame(AICcTab[2:nrow(AICcTab),1:30])
colnames(AICSum)<-param_names

AICSum2<-data.frame(colSums(AICSum))
AICsummary<-data.frame(Predictor=param_names,AICWeight=round(AICSum2[,1],3))
print(AICsummary)
sum(AICsummary$AICWeight)
m2<-subset(mtable, delta<2)



ncol(AICSum)
ggplot(df_index, aes(x=year, y=fit)) + 
  geom_ribbon(data=df_index, aes(ymin = fit-1.96*se, ymax = fit+1.96*se), 
              color='grey', alpha = 0.05) + 
  geom_line() + 
  ggtitle("MLDlarv")+
  xlab("Year") + ylab('ln (rec devs)') + 
  geom_point(aes(year,Y_rec)) + 
  scale_x_continuous(breaks=seq(1995,2025,5) , minor_breaks = 1993:max(df$year)) + 
  annotate(geom='segment', x=1993, xend=max(df_index$year), y=0, linetype='dotted') +
  theme_bw()

ggsave(paste0(fig_dir, "Predicted-recruitment-index.png"),
       height=2.0, width = 3.5)


################################################################################
# ADDITIONAL INDICES ###########################################################
################################################################################
# combined fish and environmental drivers file
df_addvars = data.frame(read.csv(paste0(data_dir,"additional-indicies-data-yellowtail.csv"), header = T))
# load terms to make equation to include stuff from best-fit model.
pars = readRDS(paste0(results_dir,'parms.rds'))
parsx = readRDS(paste0(results_dir,'parms_x.rds'))

df_main = df[,c('year',pars)]
# df_main = data_1[,c('year',pars)]
df_addvars = full_join(df_addvars, df_main)

# get predictors to create model formula #######################################
dropvars = c('year','Y_rec','sd')
envir_data = df_addvars %>% dplyr::select(!all_of(dropvars))

data_years = 1998:2014

# build model - full formula ###################################################
#quadratic_terms = c('CutiSTIpjuv')

form_dredge_2 = make_dredge_equation(envir_data = envir_data, quadratic_vars = NULL)
form_dredge_2

# save out for model testing later
saveRDS(form_dredge_2, paste0(results_dir,"formula_for_dredge_2.rds"))

# Final data preparation #######################################################
# add in parameters from original best fit model

data_2 = df_addvars %>% filter(year %in% data_years)

######## RUN DREDGE ##############################################################

fit2 = lm(form_dredge_2, data=data_2, na.action = na.fail)

# write as a function to use in model testing later
fit_dredge2 <- function(fit2){
  fit2 = fit2
  mtable = dredge(fit2, rank = AICc, m.lim = c(NA,1),
                  subset= # block correlated terms from the model
                    !(CHLpjuv && PPpjuv) && 
                    !(CHLpjuv && PPpjuv),
                  extra = list(R2 = function(x)
                    summary(x)$r.squared[[1]], 
                    F = function(y)
                      summary(y)$fstatistic[[1]]),
                  trace=2 )
  
  return(mtable)
}

# output the dredge function to use in model testing
saveRDS(fit_dredge2, paste0(results_dir,"fit_dredge_2.rds"))

# rund dredge and fit models using function ####################################
mtable2 = fit_dredge2(fit2)

mtable2_4 = subset(mtable2, delta<4)
print(mtable2_4, abbrev.names = FALSE)

saveRDS( mtable, file=paste0(results_dir,'Table_dredge-model-fits_2.rds'))
fwrite(mtable,paste0(results_dir,"Table_dredge-aic-delta4_2.csv"))

# Refit best fit model #########################################################

best_fit = find_best_fit(mtable2, fewest_params=FALSE)
bf_mod = bf_mod_equation(best_fit, params_suffix = 2)
bf_mod
saveRDS(bf_mod, paste0(results_dir, "Best_fit_model_2.rds")) # save out to use later

best_fit = lm( bf_mod, data=data_2)
best_fit = lm( bf_mod, data=data_2)
summary(best_fit)
capture.output( summary(best_fit), file = paste0(results_dir, "Best-fit-summary_2.txt"))
capture.output( anova(best_fit), file = paste0(results_dir, "Best-fit-anova_2.txt"))

# predict model ################################################################
df_index = predict_best_fit(old_data=data_2, new_data=df_addvars, bf_mod=bf_mod)
fwrite(df_index, paste0(results_dir, "Environmental-index_2.csv"))

# plot predicted and rec devs ##################################################

# df_index = na.omit(df_index)

ggplot(df_index, aes(x=year, y=fit)) + 
  geom_ribbon(data=df_index, aes(ymin = fit-1.96*se, ymax = fit+1.96*se), 
              color='grey', alpha = 0.05) + 
  geom_line() + 
  xlab("Year") + ylab('ln (rec devs)') + 
  geom_point(aes(year,Y_rec)) + 
  scale_x_continuous(breaks=seq(1995,2025,5) , minor_breaks = 1993:max(df$year)) + 
  annotate(geom='segment', x=1993, xend=max(df_index$year), y=0, linetype='dotted') +
  theme_bw()

ggsave(paste0(fig_dir, "Predicted-recruitment-index_2.png"),
       height=2.0, width = 3.5)

