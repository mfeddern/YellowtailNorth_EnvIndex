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
to_omit<-c('SCOPben', 'SCOPpjuv','NCOPben', 'NCOPpjuv','X',
           "CHLpjuv","PPpjuv","Year","Year.1",
           'HCI2pjuv',"HCI2larv", "LUSI", "BakunSTI")
envir_data2 =envir_data%>%dplyr::select(!any_of(to_omit))
recruitmentdevs<- data.frame(read.csv("data-yellowtail/RecruitmentDeviations2025draft.csv"))%>%
  filter(Datatreatment=="Expanded PacFIN")
df = envir_data2%>%left_join(recruitmentdevs)
data_years = 1994:2021 ## timeframe over which to stabilize environmental date
colnames(df)
dat = df %>% 
  dplyr::select(!any_of(to_omit))

full_dat<-dat
dat<-full_dat%>%  filter(year %in% data_years)

mod <-lm(Y_rec~ONIpjuv, data=full_dat%>%filter(year>1993&year<2019))
lm <- data.frame(predict(mod, full_dat%>%filter(year>=2019),se=TRUE,interval = "prediction"))%>%
  mutate(sefitlwr=fit.fit-1.96*se.fit, sefitupr=fit.fit+1.96*se.fit)

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

#### Model Selection: LFO5 ####

n_pred <- 5 #number of predicting years for CV
train_start<-1 #year to start training dataset
`%notin%` <- Negate(`%in%`) #assigning not in
n_year <-length(unique(dat$year)) #total number of years
n_train <- n_year-n_pred #number of training years
jstart<-n_train+1 #year to start the prediciton dataset
models <- list() #models list to fill with loop
predicted<- data.frame()  #predicted DF to fill with 1 year ahead preditions
results <- data.frame() #results to save

for (i in seq_along(combinations)) { #loop over each covariate combination i
  # k represent the number of parameters / knots estimating function at, should be small
  smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ") #generatind smooth term for each combination
  formula_str <- paste("Y_rec ~ ", smooth_terms) #generating the formula string for each smooth term
  predictions <-  numeric(n_pred ) #number of years for predictions

  # Loop over each observation j
  for (j in jstart:n_year) {
    train_index <- setdiff(train_start:n_train, j)  # Setting up training data, All indices except the j-th
    test_index <- j                 # The j-th index, first year not in the training data
    
    # Fit model on n-j observations
    gam_model <- gam(as.formula(formula_str),
                     data = dat[which(dat$year %notin% unique(dat$year)[j:n_year]), ])
    
    # Predict the excluded observation for the next year j
    predictions[which(dat$year == unique(dat$year)[j])] <- predict.gam(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    # re-fit the model
    
  }
  # re-fit the model from j loop to include j and predict j+1
  
  #fitting model from 1994 - 2013 (training fits)
  gam_model_train <- gam(as.formula(formula_str),
                         #data = dat)
                         data = dat[which(dat$year %notin% unique(dat$year)[jstart:n_year]), ])
  #fitting model from 1994 - 2018 (full model fits)
   gam_model_full <- gam(as.formula(formula_str),
                        data = dat)
  #fitting model from 2013 - 2018 (prediction fits)
  gam_model_pred <- gam(as.formula(formula_str),
                        #data = dat)
                        data = dat[which(dat$year %notin% unique(dat$year)[jstart:n_year]), ])
  #saving what we care about
  percdiff<-(na.omit(predictions[jstart:n_year])-dat$Y_rec[jstart:n_year])/dat$Y_rec[jstart:n_year]
  sumdiff<-sum(abs(percdiff))/n_pred
   mape <-sumdiff*100
   sqerror<-((predictions[jstart:n_year])-dat$Y_rec[jstart:n_year])^2
   rmse <- sqrt(mean(sqerror, na.rm=T))
   
  r2_full<-summary(gam_model_full)$r.sq
  r2_train<-summary(gam_model_train)$r.sq
  r2_pred<-summary(gam_model_pred)$r.sq
  dev.expl_full<-summary(gam_model_full)$dev.expl
  dev.expl_train<-summary(gam_model_train)$dev.expl
  dev.expl_pred<-summary(gam_model_pred)$dev.expl
  # Extract variable names
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Extract variable names
  padded_vars <- c(var_names, rep(NA, maxcovar - length(var_names)))
  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i, #model in loop i
    AIC = AIC(gam_model_full), #AIC for the full 1994 - 2018 time period
    RMSE = round(rmse,3), #RMSE for the predictions (LFO-CV 5 year)
    rsq_full=round(r2_full,2), #rsq for full 1994 - 2018
    rsq_train=round(r2_train,2), #rsq for 1994 - 2013
    rsq_pred=round(r2_pred,2), #rsq for 2014 - 2018
    dev.ex_full=round(dev.expl_full,4), #deviance explained for full 1994 - 2018
    dev.ex_train=round(dev.expl_train,4), #deviance explained for 1994 - 2013
    dev.ex_pred=round(dev.expl_pred,4), #deviance explained for 2014 - 2018
    mape=mape,
    #saving the names of each covariate in the results
    var1 = padded_vars[1], 
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var4 = padded_vars[4]
    
    
  ))
  #saving the one step ahead predictions
  predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[jstart:n_year],
    year=unique(dat$year)[jstart:n_year],
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var3 = padded_vars[4]
    
    
  ))
  
  print(i)
}

results2<-results #just storing the results object in a different object 
#results<-results2%>%filter(rsq_full>0&rsq_train>0&rsq_pred>0)
results_arr_LFO5_rmse <- arrange(results,RMSE)%>% #ordered results by RMSE if you want to use this for selection
  filter(rsq_train>0&rsq_pred>0)%>%
  mutate(rankmodRMSE=seq_along(ModelID))
minRMSE=min(results_arr_LFO5_rmse$RMSE)
results_arr_LFO5_AIC <- arrange(results_arr_LFO5_rmse,AIC)%>%#ordered results by AIC as alternative
  mutate(rankmod=seq_along(ModelID))%>%
  mutate(deltaAIC=AIC-min(AIC),deltaRMSE=RMSE-minRMSE)%>%
  mutate( percRMSE=deltaRMSE/minRMSE)
results_arr_LFO5_devex <- arrange(results,desc(dev.ex_pred)) #ordered results by deviance explained
arrange(results,RMSE) #ordered results by deviance explained
# View the results data frame
print(results_arr_LFO5_rmse) #checkin LFO results
print(results_arr_LFO5_AIC) #checking AIC results
predicted_last5<-predicted # saving predictions

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
    smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ") #setting up equation for combinations
  formula_str <- paste("Y_rec ~ ", smooth_terms) #setting up formula string
  predictions <- numeric(nrow(dat)) #number of predicted years
  n_year <- length(unique(dat$year)) #number of years
  # Loop over each observation
  for (j in jstart:n_year) {
    train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
    test_index <- j                 # The j-th index
    
    # Fit model on n-1 observations
    gam_model <- gam(as.formula(formula_str),
                     data = dat[which(dat$year != unique(dat$year)[j]), ]) #fitting the model but excluding each year
    
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

results_arr_RMSE_LOO <- arrange(results,RMSE)

results_full<-results_arr_RMSE_LOO%>%select(ModelID, RMSE,rmse_imp)%>%
  rename(RMSE_LOO=RMSE)%>%
  left_join(results_arr_LFO5_AIC)
arrange(results_full,AIC)
arrange(results_full,RMSE)

results_arr_RMSE_LOO <- arrange(results_full,RMSE_LOO)
results_arr_AIC<- arrange(results_full,AIC)
results_arr_LFO5<- arrange(results_full,RMSE)

#### Extracting LOO-CV Highest Ranked Model ####
rankmod<-1 #which model rank do you want to look at?
#the combos determine which selection you want based on rankmod
combos= c(results_arr_RMSE_LOO$var1[rankmod], results_arr_RMSE_LOO$var2[rankmod], results_arr_RMSE_LOO$var3[rankmod], results_arr_RMSE_LOO$var4[rankmod]) #this needs to be mod depending on which model and selection criteria you want
#combos= c(results_arr_LFO5_rmse$var1[rankmod], results_arr_LFO5_rmse$var2[rankmod], results_arr_LFO5_rmse$var3[rankmod], results_arr_LFO5_rmse$var4[rankmod])
#combos= c(results_arr_LFO5_rmse$var1[rankmod], results_arr_LFO5_rmse$var2[rankmod], results_arr_LFO5_rmse$var3[rankmod])

smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ") #pasting in model structure
formula_str <- paste("Y_rec ~ ", smooth_terms) #formula for selected model
#years<-1994:2013 #these are the training years
years<-1994:2019 #these are the training years

##### Generating pred error for the training period #####
#prediction intervals are a nighmare for gams...here is a not so clean hack to calculate them and the 
#associated prediction error
train_data<- full_dat%>%filter(year %in% years) 
gam <- gam(as.formula(formula_str),data = train_data)
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
train_dat_pred <- gam.train%>%mutate(uprP = fit + (crit2 * se.fit),
                                  lwrP = fit - (crit2 * se.fit),
                                  uprCI = fit + (2 * se.fit),
                                  lwrCI = fit - (2 * se.fit))%>%
  mutate(se.p=(uprP-fit)/1.96)
                    
#new_years<-2014:2018
new_years<-2020:2024
new_dat<-full_dat%>%filter(year %in% new_years) 
predicted<-data.frame(predict.gam(gam,se.fit=TRUE, newdata=new_dat))%>%
  mutate(year=new_years)
gam.predict <- cbind(new_dat, predict.gam(newdata=new_dat,gam,se.fit=TRUE), residuals= NA)
Vb <- vcov(gam) #extract vcov mat
N <- 10000
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)),V = Vb) #Cholesky factorisation of the covariance matrix
Cg <- predict(gam, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)
absDev <- abs(sweep(simDev, 1, gam.predict$se.fit, FUN = "/"))
masd <- apply(absDev, 2L, max)
crit2 <- quantile(masd, prob = 0.95, type = 8)
pred_dat_pred <- gam.predict%>%mutate(uprP = fit + (crit2 * se.fit),
                                     lwrP = fit - (crit2 * se.fit),
                                     uprCI = fit + (1.96 * se.fit),
                                     lwrCI = fit - (1.96 * se.fit))%>%
  mutate(se.p=(uprP-fit)/1.96)
gam_pred<-train_dat_pred%>%bind_rows(pred_dat_pred)
rownames(gam_pred) <- NULL


write.csv(gam_pred, "OceanographicIndexFullV2.csv")
readr::write_rds(gam_pred, "OceanographicIndexFullV2.rds")

index<-gam_pred%>%select(year, type,fit,se.fit,residuals,uprP, lwrP,uprCI,lwrCI,se.p)

write.csv(index, "OceanographicIndexV2.csv")
readr::write_rds(index, "OceanographicIndexV2.rds")


#### Extracting LFO-CV Highest Ranked Model ####
rankmod<-1 #which model rank do you want to look at?
#the combos determine which selection you want based on rankmod
combos_LFO= c(results_arr_LFO5$var1[rankmod], results_arr_LFO5$var2[rankmod], results_arr_LFO5$var3[rankmod]) #this needs to be mod depending on which model and selection criteria you want
smooth_terms_LFO <- paste("s(", combos_LFO, ", k = 3)", collapse = " + ") #pasting in model structure
formula_str_LFO <- paste("Y_rec ~ ", smooth_terms_LFO) #formula for selected model
gam_lfo <- gam(as.formula(formula_str_LFO),data = train_data)
gam.train.lfo <- cbind(train_data, predict.gam(newdata=train_data,gam_lfo,se.fit=TRUE), residuals= residuals(gam_lfo))
Vb_lfo <- vcov(gam_lfo)
BUdiff_lfo <- rmvn(N, mu = rep(0, nrow(Vb_lfo)),V = Vb_lfo)
Cg_lfo <- predict(gam_lfo, type = "lpmatrix")
simDev_lfo <- Cg_lfo %*% t(BUdiff_lfo)
absDev_lfo <- abs(sweep(simDev_lfo, 1, gam.train.lfo$se.fit, FUN = "/"))
masd_lfo <- apply(absDev_lfo, 2L, max)
crit2_lfo <- quantile(masd_lfo, prob = 0.95, type = 8)
train_dat_pred_lfo <- gam.train.lfo%>%mutate(uprP = fit + (crit2 * se.fit),
                                     lwrP = fit - (crit2 * se.fit),
                                     uprCI = fit + (2 * se.fit),
                                     lwrCI = fit - (2 * se.fit))%>%
  mutate(se.p=(uprP-fit)/1.96)

predicted_lfo<-data.frame(predict.gam(gam_lfo,se.fit=TRUE, newdata=new_dat))%>%
  mutate(year=new_years)
gam.predict.lfo <- cbind(new_dat, predict.gam(newdata=new_dat,gam_lfo,se.fit=TRUE), residuals= NA)
Vb_lfo <- vcov(gam_lfo)
BUdiff_lfo <- rmvn(N, mu = rep(0, nrow(Vb_lfo)),V = Vb_lfo)
Cg_lfo <- predict(gam_lfo, type = "lpmatrix")
simDev_lfo <- Cg_lfo %*% t(BUdiff_lfo)
absDev_lfo <- abs(sweep(simDev_lfo, 1, gam.predict.lfo$se.fit, FUN = "/"))
masd_lfo <- apply(absDev_lfo, 2L, max)
crit2_lfo <- quantile(masd_lfo, prob = 0.95, type = 8)
pred_dat_pred_lfo <- gam.predict.lfo%>%mutate(uprP = fit + (crit2 * se.fit),
                                             lwrP = fit - (crit2 * se.fit),
                                             uprCI = fit + (2 * se.fit),
                                             lwrCI = fit - (2 * se.fit))%>%
  mutate(se.p=(uprP-fit)/1.96)
gam_pred_lfo<-train_dat_pred_lfo%>%bind_rows(pred_dat_pred_lfo)
rownames(gam_pred_lfo) <- NULL
#### Plotting Highest Ranked Models ####
# Fit a model that includes multiple variables from the top model
gam_pred<-gam_pred%>%mutate(type=replace_na(gam_pred$type,"Predicted"))
gam.plot<- ggplot(gam_pred , aes(year, Y_rec)) +
  geom_point(aes(shape=type)) +
  geom_point(data=new_dat) +
    geom_line(aes(year, fit)) +
  #geom_point(data=predicted, col="red",aes(year, fit))+
  geom_line(data=predicted, col="red",aes(year, fit))+
  geom_point(data=predicted, col="red",aes(year, fit),shape=15)+
  geom_ribbon(data=gam_pred%>%filter(year<2021),aes(ymax=uprP, ymin=lwrP), 
              alpha=0.2)+
  geom_ribbon(data=predicted, fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  theme_classic() +
  labs(title="Highest Ranked: LOO-CV & AIC",
       x="Year", y="ln(Recruitment Deviations)")+
  ylim(c(-2,1.5))+
  #labs(title=paste("Deviance Explained = ",round(summary(gam)$dev.expl,2), ", R2 = ",round(summary(gam)$r.sq,2)), subtitle=formula(gam),
  #     x="Year", y="ln(Recruitment Deviations)")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.plot

png("Figures/LOO Model/ModelFit.png",width=6,height=3,units="in",res=1200)
gam.plot
dev.off()

gam_pred_lfo<-gam_pred_lfo%>%mutate(type=replace_na(gam_pred_lfo$type,"Predicted"))
gam.plot.lfo<- ggplot(gam_pred_lfo, aes(year, Y_rec)) +
  geom_point(aes(shape=type)) +
  geom_point(data=new_dat) +
  geom_line(aes(year, fit)) +
  #geom_point(data=predicted, col="red",aes(year, fit))+
  geom_line(data=predicted_lfo, col="red",aes(year, fit))+
  geom_point(data=predicted_lfo, col="red",aes(year, fit),shape=15)+
  geom_ribbon(data=gam_pred_lfo%>%filter(year<2021),aes(ymax=uprP, ymin=lwrP), 
              alpha=0.2)+
  geom_ribbon(data=gam_pred_lfo%>%filter(year>=2020), fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  theme_classic() +
  labs(title="Highest Ranked: LFO-CV",
       x="Year", y="ln(Recruitment Deviations)")+
  ylim(c(-2,1.5))+
  #labs(title=paste("Deviance Explained = ",round(summary(gam)$dev.expl,2), ", R2 = ",round(summary(gam)$r.sq,2)), subtitle=formula(gam),
  #     x="Year", y="ln(Recruitment Deviations)")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

gam.plot.lfo

png("Figures/LFO Model/ModelFit.png",width=6,height=3,units="in",res=1200)
gam.plot.lfo
dev.off()

#### Combined Figrues 3 for the Appendix #####

gam.plot.lfo2<- ggplot(gam_pred_lfo, aes(year, Y_rec)) +
  geom_point(aes(shape=type)) +
  geom_point(data=new_dat) +
  geom_line(aes(year, fit)) +
  #geom_point(data=predicted, col="red",aes(year, fit))+
  geom_line(data=predicted_lfo, col="red",aes(year, fit))+
  geom_point(data=predicted_lfo, col="red",aes(year, fit),shape=15)+
  geom_ribbon(data=gam_pred_lfo%>%filter(year<2021),aes(ymax=uprP, ymin=lwrP), 
              alpha=0.2)+
  geom_ribbon(data=gam_pred_lfo%>%filter(year>=2020), fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  theme_classic() +
  labs(title="Highest Ranked: LFO-CV",
       x="", y="ln(Recruitment Deviations)")+
  ylim(c(-2,1.5))+
  #labs(title=paste("Deviance Explained = ",round(summary(gam)$dev.expl,2), ", R2 = ",round(summary(gam)$r.sq,2)), subtitle=formula(gam),
  #     x="Year", y="ln(Recruitment Deviations)")+
  theme(strip.text = element_text(size=6),
        legend.position="none",
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.plot.lfo2

gam.plot2<- ggplot(gam_pred , aes(year, Y_rec)) +
  geom_point(aes(shape=type)) +
  geom_point(data=new_dat) +
  geom_line(aes(year, fit)) +
  #geom_point(data=predicted, col="red",aes(year, fit))+
  geom_line(data=predicted, col="red",aes(year, fit))+
  geom_point(data=predicted, col="red",aes(year, fit),shape=15)+
  geom_ribbon(data=gam_pred%>%filter(year<2021),aes(ymax=uprP, ymin=lwrP), 
              alpha=0.2)+
  geom_ribbon(data=predicted, fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  theme_classic() +
  labs(title="Highest Ranked: LOO-CV & AIC",
       x="", y="")+
  ylim(c(-2,1.5))+
  #labs(title=paste("Deviance Explained = ",round(summary(gam)$dev.expl,2), ", R2 = ",round(summary(gam)$r.sq,2)), subtitle=formula(gam),
  #     x="Year", y="ln(Recruitment Deviations)")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.plot2

png("Figures/Appendix/ModelFit.png",width=9,height=4,units="in",res=1200)
comb<-ggarrange(gam.plot.lfo2,gam.plot2, ncol = 2,widths = c(0.75, 1), labels=c("A.", "B."))
annotate_figure(comb,bottom = c("Year"))
dev.off()


#gam.res <- ggplot(gam_pred,  aes(fit,residuals)) +
#  geom_point() +
#  #geom_line(aes(year, fit)) +
#  theme_bw() +
#  geom_smooth() +
#  theme(strip.text = element_text(size=6),
#        strip.background = element_rect(fill="white"),
#        plot.title = element_text(hjust = 0.5),
#        plot.subtitle = element_text(hjust = 0.5))
#gam.res
#png("Figures/ModelResiduals-unexpanded.png",width=6,height=5,units="in",res=1200)
#gam.res
#dev.off()

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
png("Figures/LOO Model/ModelPartialsLINESPOINTS.png",width=6,height=6,units="in",res=1200)
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
