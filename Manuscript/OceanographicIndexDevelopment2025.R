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

##### Reading in the data #####
df <- data.frame(read.csv("data-yellowtail/2024Env-annual-yellowtail-standardized.csv"))
colnames(df)
envir_data = df %>% 
  rename(ONIlarv=oni_larv,ONIpjuv=oni_pjuv,ONIpre=oni_pre,
         PDOlarv=pdo_larv,PDOpjuv=pdo_pjuv,
         LUSI=lusi_annual,
         BeutiTUMI=BeutiTUMIpjuv, BeutiSTI=BeutiSTIpjuv,
         CutiSTI=CutiSTIpjuv,CutiTUMI=CutiTUMIpjuv, 
         BakunSTI=bakun_sti,
         HCIlarv=hci1_larv,HCIpjuv=hci1_pjuv,
         HCI2larv=hci2_larv,HCI2pjuv=hci2_pjuv,
         NCOPpjuv=ZOOpjuvN, NCOPben=ZOObenN,
         SCOPpjuv=ZOOpjuvS, SCOPben=ZOObenS)
to_omit<-c('X',"Year","Year.1","LUSI", "BakunSTI")
envir_data2 =envir_data%>%dplyr::select(!any_of(to_omit))
recruitmentdevs<- data.frame(read.csv("data-yellowtail/RecruitmentDeviations2025draft.csv"))%>%
  filter(Datatreatment=="Expanded PacFIN")
df = envir_data2%>%left_join(recruitmentdevs)
data_years = 1996:2021 ## timeframe over which to stabilize environmental date
colnames(df)
dat = df %>% 
  dplyr::select(!any_of(to_omit))

full_dat<-dat
dat<-full_dat%>%  filter(year %in% data_years)

#### Evaluating Correlations between Covariates #####

threshold <-0.3 #assiging a threshold of correlation to filter out 
envir_data2 <- envir_data2[complete.cases(envir_data2), ] %>% 
  dplyr::select(!any_of( c('year'))) # getting environmental data
M = data.frame(cor(envir_data2)) # generating a correlation matrix
M <- tibble::rownames_to_column(M, "xvar") #putting the rownames in their own column
M<-M%>%pivot_longer(!xvar, names_to = 'yvar', values_to = 'corr') #pivoting to a longer table
uncorr<-M%>%filter(abs(corr)<threshold) #generating the uncorrelated thresholds 
corrvar<-M%>%filter(abs(corr)>threshold)%>%filter(xvar!=yvar) #generating the correlated thresholds
corrvar2<-corrvar%>%select(-corr)
combinations_to_omit<-list() #setting up an empty list to fill with pairs to omit
for(i in 1:nrow(corrvar2)){
  combinations_to_omit[[i]]<-c(as.character(corrvar2[i, ]))
}

#### Generate Null Model ####
predict(gam(Y_rec ~ 1, data=full_dat))

null<-full_dat%>%mutate(sr_null = 0, mean_null = -0.0917625)%>%filter(year<=2019)

sqerror<-(null$sr_null-null$Y_rec)^2
rmse_sr_full <- sqrt(mean(sqerror, na.rm=T))

sqerror<-(null$mean_null-null$Y_rec)^2
rmse_mean_full <- sqrt(mean(sqerror, na.rm=T))

null_lfo<-full_dat%>%mutate(sr_null = 0, mean_null = -0.0917625)%>%filter(year>=2014 & year<2019)

sqerror<-(null_lfo$sr_null-null_lfo$Y_rec)^2
rmse_sr_lfo <- sqrt(mean(sqerror, na.rm=T))
sqerror<-(null_lfo$mean_null-null_lfo$Y_rec)^2
rmse_mean_lfo <- sqrt(mean(sqerror, na.rm=T))

#### Models to Fit ####
covariates <- names(envir_data2)
maxcovar <- 4 #max number of covars for a given model
combinations <- lapply(1:maxcovar, function(i) {
  combn(covariates, i, simplify = FALSE)
})

#setting up all possible combinations
combinations <- unlist(combinations, recursive = FALSE) 
combinations <- unique(lapply(combinations, function(x) sort(x)))
length(combinations) #all possible combinations 15275

# Function to check if a combination is a partial match of any combination to omit
is_partial_match <- function(comb, omit_list) {
  any(sapply(omit_list, function(omit) all(omit %in% comb)))
}

# Remove combinations that are partial matches (but not exact matches) to the combinations to omit
combinations <- combinations[!sapply(combinations, is_partial_match, omit_list = combinations_to_omit)]
# Check the length of remaining combinations
length(combinations) # how many left? only 887 (best model does not change compared to threshold of 0.4 but reduces model by 1/2)

#### LOO ####
models <- list()
results <- data.frame()
predicted <- data.frame()
jstart<-1

#start by dredging all possible combinations of the GAMs

for (i in seq_along(combinations)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  
  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
  predictions <- numeric(nrow(dat))
  n_year <- length(unique(dat$year))
  # Loop over each observation
  for (j in jstart:n_year) {
    train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
    test_index <- j                 # The j-th index
    
    # Fit model on n-1 observations
    gam_model <- gam(as.formula(formula_str),
                     # weights = number_cwt_estimated,
                     data = dat[which(dat$year != unique(dat$year)[j]), ])
    
    # Predict the excluded observation
    predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
  }
  
  # re-fit the model
  gam_model <- gam(as.formula(formula_str),
                   data = dat)
  
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec - predictions)^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # Extract variable names
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Store results with variable names padded to ensure there are always 3 columns
  padded_vars <- c(var_names, rep(NA, 4 - length(var_names)))
  
  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq_full=round(r2,2),
    dev.ex=round(dev.expl,4),
    rmse_imp=(rmse_sr_full-rmse)/rmse_sr_full, 
    #AUC = auc,
    #direction = direction,
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var4 = padded_vars[4]
    
    
  ))
  
  predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[jstart:n_year],
    year=unique(dat$year)[jstart:n_year],
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var3 = padded_vars[4]
    
  ))
  #saving the one step ahead predictions
  print(i)
}
results <- arrange(results,AIC)%>%
  mutate(rankmod=seq_along(ModelID))%>%
  mutate(deltaAIC=AIC-min(AIC))
arrange(results_full,AIC)
arrange(results_full,RMSE)

results_arr_RMSE<- arrange(results,RMSE)
results_arr_AIC<- arrange(results,AIC)

#### Extracting LOO-CV Highest Ranked Model ####
rankmod<-1 #which model rank do you want to look at?
combos= c(results_arr_RMSE_LOO$var1[rankmod], results_arr_RMSE_LOO$var2[rankmod], results_arr_RMSE_LOO$var3[rankmod]) #this needs to be mod depending on which model and selection criteria you want

smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ") #pasting in model structure
formula_str <- paste("Y_rec ~ ", smooth_terms) #formula for selected model
years<-1996:2021
years5<-1996:2019
##### Generating pred error for the training period #####
yearmin<- 2019
yearmax<-2021

train_data<- full_dat%>%filter(year %in% years) 
gam <- gam(as.formula(formula_str),data = train_data)
train_data5<- full_dat%>%filter(year %in% years5) 
gam5 <- gam(as.formula(formula_str),data = train_data5)
summary(gam)
gam.train <- cbind(train_data, predict.gam(newdata=train_data,gam,se.fit=TRUE), residuals= residuals(gam))
Vb <- vcov(gam)
N <- 10000
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)),V = Vb)
Cg <- predict(gam, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)
absDev <- abs(sweep(simDev, 1, gam.train$se.fit, FUN = "/"))
masd <- apply(absDev, 2L, max)
crit2 <- quantile(masd, prob = 0.95, type = 8)
train_dat_pred <- full_dat%>%left_join(gam.train%>%select(year, fit, se.fit))%>%
  mutate(uprP = fit + (crit2 * se.fit),
                                  lwrP = fit - (crit2 * se.fit),
                                  uprCI = fit + (2 * se.fit),
                                  lwrCI = fit - (2 * se.fit))%>%
  mutate(se.p=(uprP-fit)/1.96)%>%
  left_join(predicted%>%filter(ModelID==616)%>%select(year,pred))%>%
  left_join(data.frame(year=seq(yearmin,yearmax),
                    last5 = predict.gam(gam5,newdata=train_dat_pred%>%filter(year>=yearmin&year<=yearmax)))
  )%>%
  mutate(type2=ifelse(type=="Main_RecrDev","Main", "Late"))


#### Plotting Highest Ranked Model ####
# Fit a model that includes multiple variables from the top model
gam.plot<- ggplot(train_dat_pred, aes(year, Y_rec)) +
  geom_point(data=na.omit(train_dat_pred%>%filter(year>1995)%>%select(year, Y_rec, type2)),aes(shape=type2)) +
  geom_point(aes(y=last5, fill = "Jackknife"),pch=21) +
  geom_point(data=train_dat_pred%>%filter(year>=2019),aes(y=pred,fill= "Predicted"), pch=21) +
  geom_point(data=train_dat_pred%>%filter(year<2019),aes(y=pred,fill= "Predicted"), pch=24) +
  geom_line(aes(year, fit)) +
  geom_ribbon(data=train_dat_pred,aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  theme_classic() +
  labs(x="Year", y="ln(Recruitment Deviations)")+
  ylim(c(-2,1.5))+
  labs(fill = "Diagnostic", shape="Recruitment Deviation \n Type")+
  geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.plot

png("Figures/Figure2.png",width=8,height=5,units="in",res=1200)
gam.plot
dev.off()



#### Partial Effects Figures ####
# Use smooth_estimates() to get the smooth estimates and confidence intervals
smooth_data <- smooth_estimates(gam)%>%  # Or specify the smooth term you want to plot
  add_constant(model_constant(gam)) %>% # Add the intercept
  add_confint()%>% # Add the credible interval
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")%>%
  select(-c(.by))
observations<-full_dat%>%
  select(year,Y_rec,combos)%>%
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")
smooth_data <- na.omit(smooth_data )

partial_effects<-ggplot(smooth_data, aes(x = Value, y = .estimate)) +  # Setup the plot with the fitted values
  facet_wrap(~Smooth)+
  geom_line() + # Add estimated smooth
  geom_ribbon(aes(ymax = .upper_ci, ymin = .lower_ci), fill = "black", alpha = 0.2) + # Add credible interval
  geom_point(data = observations, aes(x = Value, y = Y_rec), color = "black") + # Add your data points
  labs(x = "Standardized Oceanographic Conditions", y = "Partial Residual")+ # Add labels
  geom_text(data = observations, aes(x = Value, y = Y_rec,label=year),hjust=0,nudge_x = 0.1)+
  theme_bw()
partial_effects 
png("Figures/LOO Model/ModelPartialsRMSE.png",width=6,height=6,units="in",res=1200)
partial_effects 
dev.off()


smooth_data_lfo <- smooth_estimates(gam_lfo)%>%  # Or specify the smooth term you want to plot
  add_constant(model_constant(gam_lfo)) %>% # Add the intercept
  add_confint()%>% # Add the credible interval
  pivot_longer(c(combos_LFO),names_to = "Smooth",values_to ="Value")%>%
  select(-c(.by))
observations_LFO<-full_dat%>%
  select(year,Y_rec,combos_LFO)%>%
  pivot_longer(c(combos_LFO),names_to = "Smooth",values_to ="Value")
smooth_data_lfo <- na.omit(smooth_data_lfo)
partial_effects_lfo<-ggplot(smooth_data_lfo, aes(x = Value, y = .estimate)) +  # Setup the plot with the fitted values
  facet_wrap(~Smooth)+
  geom_line() + # Add estimated smooth
  geom_ribbon(aes(ymax = .upper_ci, ymin = .lower_ci), fill = "black", alpha = 0.2) + # Add credible interval
  geom_point(data = observations_LFO, aes(x = Value, y = Y_rec), color = "black") + # Add your data points
  labs(x = "Standardized Oceanographic Conditions", y = "Partial Residual")+ # Add labels
  geom_text(data = observations_LFO, aes(x = Value, y = Y_rec,label=year),hjust=0,nudge_x = 0.1)+
  theme_bw()
partial_effects_lfo 
png("Figures/LFO Model/ModelPartialsRMSE.png",width=6,height=6,units="in",res=1200)
partial_effects_lfo
dev.off()

#partial_effects<-ggplot(smooth_data, aes(x = Value, y = .estimate)) +  # Setup the plot with the fitted values
#  facet_wrap(~Smooth)+
#  geom_line() + # Add estimated smooth
#  geom_ribbon(aes(ymax = .upper_ci, ymin = .lower_ci), fill = "black", alpha = 0.2) + # Add credible interval
#  geom_point(data = observations%>%filter(year>2021), aes(x = Value, y = Y_rec), color = "black") + # Add your data points
#  labs(x = "Standardized Oceanographic Conditions", y = "Partial Residual")+ # Add labels
#  geom_text(data = observations%>%filter(year>2021), aes(x = Value, y = Y_rec,label=year),hjust=0,nudge_x = 0.1)+
#  theme_bw()
#partial_effects 
#png("Figures/ModelPartials.png",width=6,height=6,units="in",res=1200)
#partial_effects 
#dev.off()

#ts<-full_dat%>%select(c(combos,year))
#tscov<-ggplot(pivot_longer(ts,col=combos,names_to = 'var',values_to = "val"),aes(x=year,y=val))+
#  geom_line()+
#  geom_point()+
#  facet_wrap(~var)+
#  ylab("Oceanographic Index")+
#  xlab("Year")+
#  geom_hline(yintercept = 0,lty=2)+
#  theme_bw()
#tscov

#png("Figures/OceanographicTS.png",width=4,height=4,units="in",res=1200)
#tscov
#dev.off()
#colnames(full_dat)

tscov<-ggplot(pivot_longer(full_dat%>%select(-c(HCIpjuv, HCIlarv)),cols=-c(year,Datatreatment,sd,type,Y_rec),names_to = 'var',values_to = "val"),aes(x=year,y=val))+
  geom_line()+
  geom_point()+
  facet_wrap(~var)+
  ylab("Oceanographic Index")+
  xlab("Year")+
  geom_hline(yintercept = 0,lty=2)+
  theme_bw()
tscov

png("Figures/Appendix/FullTS.png",width=9,height=6,units="in",res=1200)
tscov
dev.off()

FWts<-envir_data%>%select('SCOPben', 'SCOPpjuv','NCOPpjuv','NCOPben',"CHLpjuv","PPpjuv","year")%>%
pivot_longer(cols=-c(year),names_to = 'var',values_to = "val")
  
FWtsplot<-ggplot(FWts,aes(x=year,y=val))+
  facet_wrap(~var)+
   geom_line()+
  geom_point()+
  ylab("Food Web Indices")+
  xlab("Year")+
  geom_hline(yintercept = 0,lty=2)+
  theme_bw()
png("Figures/FwTS.png",width=6,height=6,units="in",res=1200)
FWtsplot
dev.off()

#### GAM diagnostics ####

#checking for multicollinearity
dfx = envir_data2%>%select(-c(HCIpjuv, HCIlarv)) 
dfx = na.omit(dfx) # na's mess up the cor function
cor_xy = cor(dfx)

graphics.off()
png( paste0("Figures/Appendix/oceanographic-correlations-among-variables.png"),
     units = 'in', res=300, width = 6.5, height=6.5)
# plot command
corrplot::corrplot(cor_xy, hod='circle', type='lower', tl.cex = 0.5)
dev.off()


#checking for multicollinearity
dfx = envir_data2%>%select(-c(HCIpjuv, HCIlarv))
dfx = na.omit(dfx) # na's mess up the cor function
cor_xy = cor(dfx)

graphics.off()
png( paste0("Figures/Appendix/oceanographic-correlations-among-variables.png"),
     units = 'in', res=300, width = 12.5, height=12.5)
# plot command
corrplot::corrplot(cor_xy, hod='circle', method='number',type='lower', tl.cex = 1.5)
dev.off()


#combos= c(results_arr_LFO5_AIC$var1[rankmod], results_arr_LFO5_AIC$var2[rankmod], results_arr_LFO5_AIC$var3[rankmod], results_arr_LFO5_AIC$var4[rankmod])
#linear_terms <- paste(combos, collapse = " + ") #pasting in model structure
#formula_str_lm <- paste("Y_rec ~ ", linear_terms) #formula for selected model
#vifs<-vif(lm(as.formula(formula_str_lm),data = train_data))
#vifs.df<-data.frame(vifs)

#resultsSave<-results_arr_LFO5_AIC%>%
#  select(rankmod, var1, var2, var3, var4,AIC,deltaAIC,RMSE,mape,rsq_full,dev.ex_full)

#write.csv(vifs.df,"VIFs.csv")
#write.csv(M,"Correlation.csv")
#write.csv(resultsSave%>%filter(rankmod<=5),"BestModelsunxpanded.csv")

#looking at qq, resid, and response with all rec devs
rankmod<-1 #which model rank do you want to look at?
combos= c(results_arr_LOO_rmse$var1[rankmod], results_arr_LFO5_rmse$var2[rankmod], results_arr_LFO5_rmse$var3[rankmod], results_arr_LFO5_rmse$var4[rankmod])
combos= c(results_arr_LFO5_rmse$var1[rankmod], results_arr_LFO5_rmse$var2[rankmod], results_arr_LFO5_rmse$var3[rankmod])

smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ") #pasting in model structure
formula_str <- paste("Y_rec ~ ", smooth_terms) #formula for selected model
#years<-1994:2013 #these are the training years, 2019 is a bit critical since LST goes off the rails
years<-1994:2019 #these are the training years, 2019 is a bit critical since LST goes off the rails
train_data<- full_dat%>%filter(year %in% years) 
gam <- gam(as.formula(formula_str),data = train_data)


m2 = gam(as.formula(formula_str),
data = diag_dat)
#png("Figures/gamdiagRMSE.png",width=6,height=6,units="in",res=1200)
#par(mfrow=c(2,2))
#gam.check(m2,pch=19,cex=.3)
#dev.off()

#looking for leverage years
m2<-gam
diag_dat<-full_dat%>%filter(year %in% 1994:2019) 
diag_dat$cooks.distance <- cooks.distance(m2)
diag_dat$leverage <-influence.gam(m2)
diag_dat$residuals.standard <- scale(residuals.gam(m2))
n <- length(diag_dat$residuals.standard )  # Number of data points
influence  <- diag_dat[(diag_dat$cooks.distance> 4 / n |diag_dat$cooks.distance> 1|diag_dat$cooks.distance> 0.5),]


cooks<-ggplot(diag_dat, aes(x = cooks.distance, y = residuals.standard)) +
  geom_hline(yintercept = c(-3,3), linetype = 'dotted') +
  #geom_vline(xintercept = 2* mean(diag_dat$leverage), linetype = 'dotted') +
  #geom_vline(xintercept = (2*3)/length(1994:2021), linetype = 'dotted') +
  geom_point(alpha = 0.2) +
  ylab("Standardized Residuals")+
  xlab("Cook's Distance")+
  ggtitle("Influence")+
  geom_point(data = influence, aes(x = cooks.distance, y = residuals.standard), color = 'darkorange') +
  theme_minimal() +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5),  # Center plot title
    axis.title = element_text(size = 12),    # Axis title size
    axis.text = element_text(size = 10)      # Axis text size
  ) +
  #ylim(c(-4,4))+
  #geom_text(aes(label=ifelse(cooks.distance > 4 / n, year, "")), cex=3, hjust=1.1)
  geom_text(aes(label=year),cex=3, hjust=1.1)

cooks
png("Figures/CooksRMSE.png",width=5,height=5,units="in",res=1200)
cooks
dev.off()

### Boot Strapping ####
fm = gam_model
sfm = summary(fm)
#fmedf = t(data.frame(sfm$edf))
Fstats = t(data.frame(c(sfm$s.table[1:3,3])))
Pstats = t(data.frame(c(sfm$s.table[1:3,4])))
fitted.model  = data.frame(cbind(sfm$dev.expl,  sfm$r.sq, Pstats,Fstats))
data_2 = train_data # duplicate so can bootstrap
results_boot <- data.frame()
for(k in 1:1000) {# refit model 1000 times
  # create new recrutment time series
  # BOOTSTRAP sampling with replacement 
  # BOOTSTRAP RECRUITMENT not RperS 
  print(paste("k = ", k))
  ROW = sample(1:nrow(train_data),replace=TRUE)
  # Reorder Y_rec --- NOT WHoLE ROWS
  data_2$Y_rec = train_data$Y_rec[ROW]   
  fit = gam(as.formula(formula_str), data=data_2)
  s1 = summary(fit)
  Fs = t(data.frame(s1$s.table[,3]))
  Ps = t(data.frame(s1$s.table[,3]))
  r1  =  data.frame(cbind(s1$dev.expl,  s1$r.sq, Ps, Fs))
  results_boot<-rbind(r1,results_boot)
} # end k loop
colnames(results_boot) = c('r2','devex',"P1","P2", "P3","P4","F1","F2","F3","F4")
#colnames(results_boot) = c('r2','devex',"P1","P2", "P3","F1","F2","F3")

results_boot$p = 1-pf(results_boot$F1,results_boot$F2, results_boot$F3, results_boot$F4)

#write.table(results, "R_boot_not.std.csv", sep=',',col.names = TRUE, row.names = FALSE)

# get mean and 95% CLs
mn = apply(results_boot,2,mean)
md = apply(results_boot,2,median)
se = apply(results_boot,2,sd)
# quantile function
ci = apply(results_boot,2, qt)
boot.stats = rbind(mn,md, se, ci)
fitted.model$p = NA

colnames(fitted.model) <- colnames(boot.stats)  
fitted.model$p = 1-pf(fitted.model$F1,fitted.model$F2, fitted.model$F3, fitted.model$F4)
boot.stats2 = data.frame(rbind(fitted.model,boot.stats))
bias = mn - boot.stats2[1,] 
mn.corrected = boot.stats2[1,] - bias
# mn.corrected = 2*boot.stats2[1,] - mn
boot.stats3 = rbind(boot.stats2, bias, mn.corrected)
x = rownames(boot.stats3)
x[1] <- "fitted model"
x[length(x)-1] = 'bias'
x[length(x)] = 'mn.corrected'
boot.results <- cbind(x,boot.stats3)
boot.results
#write.table(boot.results,'Stats_boot_not.std.csv', sep=',',col.names = TRUE, row.names = FALSE)

#### Jackknifing ####
#start by dredging all possible combinations of the GAMs
results_loo <- data.frame()
jstart<-1

predictions <- numeric(nrow(dat))
rsqjack <- data.frame()
n_year <- length(unique(dat$year))

# LOO cv

for (j in jstart:n_year) {
  train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
  test_index <- j                 # The j-th index
  # Fit model on n-1 observations
  gam_model <- gam(as.formula(formula_str),
                   # weights = number_cwt_estimated,
                   data = dat[which(dat$year != unique(dat$year)[j]), ])
  # Predict the excluded observation
  predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
  rsq<-summary(gam_model)$r.sq
  rsqjack<-rbind(rsqjack,rsq)
  }

gam.jack<-cbind(gam_pred%>%filter(year>1993&year<2022),jack=predictions[1:28],r2=rsqjack[1:28,1])
jack<-ggplot(gam.jack, aes(year, fit)) +
  #geom_point(aes(shape=Type)) +
  #geom_point(data=new_dat) +
  geom_point(aes(y=jack),bg='yellow',pch=21,cex=1.5,col="black") +
  geom_point(aes(y=Y_rec), bg="red",pch=23,cex=1)+
  #geom_line(data=predicted, col="red",aes(year, fit))+
  #geom_point(data=predicted, col="red",aes(year, fit),shape=15)+
  #geom_line(aes(year, fit)) +
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=uprP, ymin= lwrP), 
              alpha=0.2)+
  #geom_ribbon(data=predicted, fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
  #            alpha=0.2)+
  theme_bw() +
  xlim(c(1993,2019))+
  #ylim(c(-1.5,2))+
 labs(#title="Jackknife",
       x="", y="")+
  ylim(c(-1.5,1.5))+
  
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

jack


#lfo cv
predictions <- numeric(nrow(dat))
rsqjack <- data.frame()

for (j in jstart:n_year) {
  train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
  test_index <- j                 # The j-th index
  # Fit model on n-1 observations
  gam_model <- gam(as.formula(formula_str_LFO),
                   # weights = number_cwt_estimated,
                   data = dat[which(dat$year != unique(dat$year)[j]), ])
  # Predict the excluded observation
  predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
  rsq<-summary(gam_model)$r.sq
  rsqjack<-rbind(rsqjack,rsq)
}
gam.jack_lfo<-cbind(gam_pred_lfo%>%filter(year>1993&year<2022),jack=predictions[1:28],r2=rsqjack[1:28,1])

jack_lfo<-ggplot(gam.jack_lfo, aes(year, fit)) +
  #geom_point(aes(shape=Type)) +
  #geom_point(data=new_dat) +
  geom_point(aes(y=jack),bg='yellow',pch=21,cex=1.5,col="black") +
  #geom_point(aes(y=Y_rec),bg='white',pch=21,col="black") +
  geom_point(aes(y=Y_rec), bg="red",pch=23,cex=1)+
  #geom_point(data=predicted, col="red",aes(year, fit))+
  #geom_line(data=predicted, col="red",aes(year, fit))+
  #geom_point(data=predicted, col="red",aes(year, fit),shape=15)+
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=uprP, ymin= lwrP), 
              alpha=0.2)+
  #geom_ribbon(data=predicted, fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
  #            alpha=0.2)+
  theme_bw() +
  ylim(c(-1.5,1.5))+
  xlim(c(1993,2019))+
  #ylim(c(-1.5,2))+
 labs(#title="Jackknife",
       x="", y="ln(Recruitment Deviations)")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

jack_lfo




png("Figures/LOO Model/Jack.png",width=5,height=3.5,units="in",res=1200)
jack
dev.off()

png("Figures/LFO Model/Jack.png",width=5,height=3.5,units="in",res=1200)
jack_lfo
dev.off()

png("Figures/Appendix/Jackknife.png",width=9,height=4,units="in",res=1200)
comb<-ggarrange(jack_lfo,jack, ncol = 2, labels=c("A.", "B."))
annotate_figure(comb,bottom = c("Year"))
dev.off()


#### R2 comparison ####

r2jack<-ggplot(gam.jack, aes(r2)) +
  geom_histogram(bins=10,fill="lightgrey", col="black")+
  xlim(c(0,0.8))+
  ylab("Frequency")+
  xlab("R-squared")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

r2jack2<-ggplot(gam.jack, aes(x=year, y=r2)) +
  geom_point(col="black")+
 # xlim(c(0.4,0.8))+
  ylim(c(0,1))+
  ylab("R-squared")+
  geom_hline(yintercept=mean(na.omit(gam.jack$r2)),lty=2)+#0.59
  xlab("Year (removed)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))




r2jacklfo<-ggplot(gam.jack_lfo, aes(r2)) +
  geom_histogram(bins=10,fill="lightgrey", col="black")+
  xlim(c(0,0.8))+
  ylab("Frequency")+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

r2jack2lfo<-ggplot(gam.jack_lfo, aes(x=year, y=r2)) +
  geom_point(col="black")+
  # xlim(c(0.4,0.8))+
  ylim(c(0,1))+
  ylab("R-squared")+
  geom_hline(yintercept=mean(na.omit(gam.jack_lfo$r2)),lty=2)+#0.59
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


png("Figures/Appendix/Jackr2rmse.png",width=6,height=5,units="in",res=1200)
ggarrange(r2jacklfo, r2jack2lfo,r2jack,r2jack2,ncol=2,nrow=2,labels=c("A.", "B.", "C.","D."))
dev.off()


### Nick's Jackknifing ####
# add to file below ###
fm = gam
sfm = summary(fm)
#fmedf = t(data.frame(sfm$edf))
Fstats = t(data.frame(c(sfm$s.table[1:4,3])))
Pstats = t(data.frame(c(sfm$s.table[1:4,4])))
fitted.model  = data.frame(cbind(sfm$dev.expl,  sfm$r.sq, Pstats,Fstats))
data_2 = train_data # duplicate so can bootstrap
results_boot <- data.frame()
rmse_fitted = NA
bf = summary(fm)
Pi = predict.gam(fm)
Oi = train_data$Y_rec
Coeffs = t(data.frame(bf$edf))
rmse_fitted = sqrt((sum(Pi - Oi)^2) / length(Oi))

fitted.model = data.frame(cbind(rmse_fitted,  sfm$r.sq,  sfm$dev.expl,  Fstats, Coeffs))
results_jackk<-data.frame()
predicted<-data.frame()
for(k in 1:nrow(train_data)) {# refit model 30 times dropping one year each time
  # jackknifing here drop one datum and run
  data_2 = train_data[-k,]
  print(k)
  fit = gam(as.formula(formula_str), data=data_2, na.action = na.fail)
  s1 = summary(fit)
  Coeffs = t(data.frame(s1$edf))
  Fs = t(data.frame(s1$s.table[1:4,3]))
  # cross validation here
  data_cross = data.frame(train_data[k,])
  pi = predict.gam(fit, newdata = data_cross)
  oi = train_data[k,'Y_rec']
  rmse = sqrt((sum(pi - oi)^2) / length(oi))
  
  r1  = data.frame(cbind(rmse,  s1$r.sq,  s1$dev.expl, pi,oi, Fs, Coeffs))
  predicted<-rbind(predicted,pi[1])
  results_jackk<-rbind(r1,results_jackk)
} # end k loop
colnames(results_jackk) = c('RMSE','r2','devex',"jack","pred","F1","F2","F3","F4","Coef1","Coef2","Coef3","Coef4")
results_jackk$p = 1-pf(results_jackk$F1,results_jackk$F2, results_jackk$F3, results_jackk$F4)

#write.table(results,'R_jackknife-best-fit.csv', sep=',',col.names = TRUE, row.names = FALSE)

# uses same code as bootstrap version, but not bootstrping
# get mean and 95% CLs

mn = apply(results_jackk,2,mean)
md = apply(results_jackk,2,median)
# quantile function
ci = apply(results_jackk,2, qt)
boot.stats = rbind(mn,md, ci)
fitted.model$p = NA

boot.stats

colnames(fitted.model) <- colnames(boot.stats)
fitted.model$p = 1-pf(fitted.model$F1,fitted.model$F2, fitted.model$F3, fitted.model$F4)
boot.stats2 = data.frame(rbind(fitted.model,boot.stats))
x = rownames(boot.stats2)
x[1] <- "fitted model"
boot.results <- cbind(x,boot.stats2)
boot.final =  boot.results #rbind(boot.results,CL2)
boot.final
#write.table(boot.final,'R_jackknife-best-fit-stats.csv', sep=',',col.names = TRUE, row.names = FALSE)


#### LOO ####


cross_validation <-TRUE
models <- list()
results <- data.frame()
predicted <- data.frame()
jstart<-1

#start by dredging all possible combinations of the GAMs

for (i in seq_along(combinations)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  
  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
  predictions <- numeric(nrow(dat))
  n_year <- length(unique(dat$year))
  # Loop over each observation
  for (j in jstart:n_year) {
    train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
    test_index <- j                 # The j-th index
    
    # Fit model on n-1 observations
    gam_model <- gam(as.formula(formula_str),
                     # weights = number_cwt_estimated,
                     data = dat[which(dat$year != unique(dat$year)[j]), ])
    
    # Predict the excluded observation
    predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
  }
  
  # re-fit the model
  gam_model <- gam(as.formula(formula_str),
                   data = dat)
  
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec - predictions)^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # Extract variable names
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Store results with variable names padded to ensure there are always 3 columns
  padded_vars <- c(var_names, rep(NA, 4 - length(var_names)))
  
  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq_full=round(r2,2),
    dev.ex=round(dev.expl,4),
    #AUC = auc,
    #direction = direction,
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var4 = padded_vars[4]
    
    
  ))
  
  predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[jstart:n_year],
    year=unique(dat$year)[jstart:n_year],
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var3 = padded_vars[4]
    
  ))
  #saving the one step ahead predictions
  print(i)
}

results_arr_LFO5_AIC
results_arr_RMSE_LOO <- arrange(results,RMSE)

results_full<-results_arr_RMSE_LOO%>%select(ModelID, RMSE)%>%
  rename(RMSE_LOO=RMSE)%>%
  left_join(results_arr_LFO5_AIC)
arrange(results_full,AIC)
arrange(results_full,RMSE)

gam.loo<-predicted%>%
  filter(ModelID==509)%>%
  left_join(gam.predict)

gam.jack<-cbind(gam.predict,jack=predictions[1:26],r2=rsqjack[1:26,1])
gam.loo<-ggplot(gam.loo, aes(year, fit)) +
  #geom_point(aes(shape=Type)) +
  #geom_point(data=new_dat) +
  geom_point(aes(y=pred),bg='yellow',pch=21,col="black") +
  #geom_point(data=predicted, col="red",aes(year, fit))+
  #geom_line(data=predicted, col="red",aes(year, fit))+
  #geom_point(data=predicted, col="red",aes(year, fit),shape=15)+
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  #geom_ribbon(data=predicted, fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
  #            alpha=0.2)+
  theme_bw() +
  xlim(c(1993,2019))+
  #ylim(c(-1.5,2))+
  labs(#title="Jackknife",
    x="Year", y="ln(Recruitment Deviations)")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

gam.loo

png("Figures/gam.loo.png",width=5,height=3.5,units="in",res=1200)
gam.loo
dev.off()

arrange(results_arr_LFO5_AIC,RMSE)%>%#ordered results by AIC as alternative
  mutate(rankmod=seq_along(ModelID))%>%
  filter(ModelID==509)
arrange(results_arr_LFO5_AIC,mape)%>%#ordered results by AIC as alternative
  mutate(rankmod=seq_along(ModelID))%>%
  filter(ModelID==509)
arrange(results_arr_LFO5_AIC,desc(rsq_full))%>%#ordered results by AIC as alternative
  mutate(rankmod=seq_along(ModelID))%>%
  filter(ModelID==509)
arrange(results_arr_LFO5_AIC,desc(dev.ex_train))%>%#ordered results by AIC as alternative
  mutate(rankmod=seq_along(ModelID))%>%
  filter(ModelID==509)
arrange(results_full,RMSE_LOO)%>%#ordered results by AIC as alternative
  mutate(rankmod=seq_along(ModelID))%>%
  filter(ModelID==509)


png("Figures/LOO Model/gamdiagRMSE.png",width=6,height=6,units="in",res=1200)
par(mfrow=c(2,2))
gam.check(gam,pch=19,cex=.3)
dev.off()


png("Figures/LFO Model/gamdiagRMSE.png",width=6,height=6,units="in",res=1200)
par(mfrow=c(2,2))
gam.check(gam_lfo,pch=19,cex=.3)
dev.off()
summary(gam_lfo)
summary(gam)
