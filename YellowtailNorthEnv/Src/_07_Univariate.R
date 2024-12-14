library(dplyr)
library(tidyr)
library(corrplot)
library(mgcv)
library(pROC)
library(DHARMa)
library(mgcViz)
library(gridExtra)
library(ROCR)
library(recdata)
library(predRecruit)


df = data.frame(read.csv("data-yellowtail/DATA_Combined_glorys_yellowtail.csv"))
dfa = data.frame(read.csv("data-yellowtail/dfa_trend.csv"))
data_years = 1994:2014
dat = df %>% 
  dplyr::select(!any_of( c('ZOOpjuv', 'ZOOben')))  %>% 
  filter(year %in% data_years)%>%
  left_join(dfa)
envir_data = dat %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','year','Y_rec','ZOOpjuv','ZOOben')))%>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI','ONIlarv')))%>% 
  mutate_all(~ scale(.)) # drop terms not in model statement

head(envir_data)
dim(envir_data)

summary(gam( BeutiSTIpjuv~ ONIpjuv, data=envir_data))$r.sq
# Vector of marine covariates for univariate GAM loop

covariates <- names(envir_data)#[-which(names(ocean) %in% "year")]

####### LMS ######

#### LMS with Leave One OUt Cross Validation ####

cross_validation <-TRUE
models <- list()
results <- data.frame()
jstart<-1
n_year <-length(unique(dat$year))


for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  
  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
 
    smooth_terms <- paste(covariates[[i]], collapse = " + ")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
    predictions <- numeric(nrow(dat))
    n_year <- length(unique(dat$year))
    # Loop over each observation
    for (j in jstart:n_year) {
      train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index
      
      # Fit model on n-1 observations
      lm_model <- lm(as.formula(formula_str),
                       # weights = number_cwt_estimated,
                       data = dat[which(dat$year != unique(dat$year)[j]), ])
      
      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[j])] <- predict(lm_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    
    
    # re-fit the model
    lm_model <- lm(as.formula(formula_str),
                     data = dat)
  } 
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec - predictions)^2, na.rm=T))
  r2<-summary(lm_model)$r.sq
  dev.expl<-summary(lm_model)$dev.expl
  # Extract variable names
  var_name <- covariates[i]
  
  # Store results
  models[[i]] <- lm_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(lm_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    #AUC = auc,
    #direction = direction,
    var = var_name
    
  ))
  print(i)
}



# filter for frowns only

results_arr <- arrange(results,RMSE )
# View the results data frame
print(results_arr)

# Forloop through model results to plot residuals for each covariate with convex
# fit

for(i in 1:length(models)) {
  if(i==1) {
    resid_df <- data.frame("fitted" = fitted(models[[i]]),
                           "residuals" = residuals(models[[i]]),
                           "var" = results$var[i],
                           "year"=dat$year,
                           "Y_Rec"=dat$Y_rec)
  } else {
    resid_df <- rbind(resid_df, data.frame("fitted" = fitted(models[[i]]),
                                           "residuals" = residuals(models[[i]]),
                                           "var" = results$var[i],
                                           "year"=dat$year,
                                           "Y_Rec"=dat$Y_rec))
  }
}

results_arr_LOO_LM <- arrange(results,RMSE )%>%
  rename(AIC_LOO=AIC,RMSE_LOO=RMSE, rsq_LOO=rsq)%>%
  mutate(model="LM")

model_fits_LOO<-resid_df%>%
  rename(fitted_LOO=fitted, residuals_loo=residuals)
#### lmS with leave future out cross valudation ####

n_pred <- 10
train_start<-1
`%notin%` <- Negate(`%in%`)
n_year <-length(unique(dat$year))
n_train <- n_year-n_pred 
predicted<- data.frame()
kstart<-n_train+1
models <- list()
results <- data.frame()
for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  
  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste(covariates[[i]], collapse = " + ")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
  
  predictions <-  numeric(n_pred )
  
  # Loop over each observation
  for (k in kstart:n_year) {
    train_index <- setdiff(train_start:n_train, k)  # All indices except the k-th
    test_index <- k                 # The k-th index
    
    # Fit model on n-k observations
    lm_model <- lm(as.formula(formula_str),
                     data = dat[which(dat$year %notin% unique(dat$year)[k:n_year]), ])
    
    # Predict the excluded observation
    predictions[which(dat$year == unique(dat$year)[k])] <- predict(lm_model, newdata = dat[which(dat$year == unique(dat$year)[k]), ])
    # re-fit the model
    
  }
  # re-fit the model
  lm_model <- lm(as.formula(formula_str),
                   #data = dat)
                   data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
  lm_model2 <- lm(as.formula(formula_str),
                    data = dat)
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec[kstart:n_year] - predictions[kstart:n_year])^2, na.rm=T))
  r2<-summary(lm_model)$r.sq
  dev.expl<-summary(lm_model)$dev.expl
  # Extract variable names
  var_name <- covariates[i]
  
  # Store results
  models[[i]] <- lm_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(lm_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    #dev.ex=round(dev.expl,4),
    #AUC = auc,
    #direction = direction,
    var = var_name
    
  ))
  predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[kstart:n_year],
    var = var_name
    
  ))
  
  print(i)
}
#results_arr_LFO5_lm <- arrange(results,RMSE )%>%
#  rename(AIC_LFO5=AIC,RMSE_LFO5=RMSE, rsq_LFO5=rsq)%>%
#  mutate(model="LM")
 results_arr_LFO10_lm <- arrange(results,RMSE )%>%
rename(AIC_LFO10=AIC,RMSE_LFO10=RMSE, rsq_LFO10=rsq)%>%
  mutate(model="LM")
# View the results data frame
print(results_arr)

# Forloop through model results to plot residuals for each covariate with convex
# fit
resid_df <- data.frame()
for(i in 1:length(models)) {
  
  resid_df_temp <- data.frame("fits" = fitted(models[[i]]),
                              "residuals" = residuals(models[[i]]),
                              "var" = results$var[i],
                              "year"=dat$year[1:kstart-1],
                              "Y_Rec"=dat$Y_rec[1:kstart-1],
                              "sd"=dat$sd[1:kstart-1])
  resid_df<-rbind(resid_df_temp, resid_df)
}
predicted<-data.frame(predicted, "year"=dat$year[kstart:n_year], "Y_Rec"=dat$Y_rec[kstart:n_year])
#resid_df <- resid_df[complete.cases(resid_df),] # we want no missing values
#resid_df <- resid_df[is.element(resid_df$var, v),]

g2 <- ggplot(resid_df, aes(year,fits)) +
  geom_line() +
  geom_ribbon(aes(x=year, y=fits, ymax=fits+sd, ymin=fits-sd), 
              alpha=0.2)+
  geom_point(data=predicted,aes(x=year,y=Y_Rec),col="red") +
  geom_point(aes(y=Y_Rec))+
  geom_line(data=predicted,aes(x=year,y=pred,col="red"))+
  ggtitle("LM Fits")+
  facet_wrap(~ var, scale="free_x") +
  theme_bw() +
  
  xlab("Year") + ylab("Fitted") +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"))
g2
#ggsave("figures-yellowtail/LMFits10.png", height = 8, width = 10)

  facet_wrap(~ var, scale="free_x") +
  theme_bw() +
  
  xlab("Year") + ylab("Fitted") +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"))
g2 

fitted_univariate_LMS<-rbind(predicted%>%
                               rename(fits=pred)%>%mutate(sd=NA, residuals=NA, Type="predicted")%>%
                               select(-ModelID),resid_df%>%select(fits, var, year, Y_Rec, sd, residuals)%>%
                               mutate(Type="fitted"))%>%
  rename(fitted_LFO=fits, residual_LFO=residuals)%>%
  right_join(model_fits_LOO)
results_full_lm <-results_arr_LFO5_lm%>%
  left_join(results_arr_LOO_LM%>%select(-ModelID))%>%
  left_join(results_arr_LFO10_lm%>%select(-ModelID))

  
write.csv(arrange(results_full_lm,RMSE_LOO),"univariateLM_Results.csv")


write.csv(arrange(results_full,RMSE_LOO),"univariateLM_Results.csv")
write.csv(fitted_univariate_LMS,"univariateLM_Fits.csv")

######## GAMS #######
#### GAMS with Leave One OUt Cross Validation ####
models <- list()
results <- data.frame()
jstart<-1
n_year <-length(unique(dat$year))


for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small

  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", covariates[[i]], ", k = 3)")
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
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    #AUC = auc,
    #direction = direction,
    var = var_name
    
  ))
  print(i)
}



# filter for frowns only

results_arr <- arrange(results,RMSE )
# View the results data frame
print(results_arr)

# Forloop through model results to plot residuals for each covariate with convex
# fit

for(i in 1:length(models)) {
  if(i==1) {
    resid_df <- data.frame("fitted" = fitted(models[[i]]),
                           "residuals" = residuals(models[[i]]),
                           "var" = results$var[i],
                           "year"=dat$year,
                           "Y_Rec"=dat$Y_rec)
  } else {
    resid_df <- rbind(resid_df, data.frame("fitted" = fitted(models[[i]]),
                                           "residuals" = residuals(models[[i]]),
                                           "var" = results$var[i],
                           "year"=dat$year,
                           "Y_Rec"=dat$Y_rec))
  }
}

results_arr_LOO_gam <- arrange(results,RMSE )%>%
  rename(AIC_LOO=AIC,RMSE_LOO=RMSE, rsq_LOO=rsq, dev.ex_LOO=dev.ex)
model_fits_LOO<-resid_df%>%
  rename(fitted_LOO=fitted, residuals_loo=residuals)



#### GAMS with leave future out cross valudation ####

n_pred <- 10
train_start<-1
`%notin%` <- Negate(`%in%`)
n_year <-length(unique(dat$year))
n_train <- n_year-n_pred 
predicted<- data.frame()
kstart<-n_train+1
models <- list()
results <- data.frame()
for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small

  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", covariates[[i]], ", k = 3)")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
predictions <-  numeric(n_pred )

    # Loop over each observation
    for (k in kstart:n_year) {
      train_index <- setdiff(train_start:n_train, k)  # All indices except the k-th
      test_index <- k                 # The k-th index

      # Fit model on n-k observations
      gam_model <- gam(as.formula(formula_str),
                     data = dat[which(dat$year %notin% unique(dat$year)[k:n_year]), ])

      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[k])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[k]), ])
        # re-fit the model

  }
    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     #data = dat)
                     data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
        gam_model2 <- gam(as.formula(formula_str),
                     data = dat)
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec[kstart:n_year] - predictions[kstart:n_year])^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # Extract variable names
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    #AUC = auc,
    #direction = direction,
    var = var_name
    
  ))
    predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[kstart:n_year],
    var = var_name
    
  ))
 
   print(i)
}
  results_arr_LFO5_gam <- arrange(results,RMSE )%>%
  rename(AIC_LFO5=AIC,RMSE_LFO5=RMSE, rsq_LFO5=rsq, dev.ex_LFO5=dev.ex)
 # results_arr_LFO10_gam <- arrange(results,RMSE )%>%
#  rename(AIC_LFO10=AIC,RMSE_LFO10=RMSE, rsq_LFO10=rsq, dev.ex_LFO10=dev.ex)
# View the results data frame
print(results_arr)

# Forloop through model results to plot residuals for each covariate with convex
# fit
resid_df <- data.frame()
for(i in 1:length(models)) {
 
    resid_df_temp <- data.frame("fits" = fitted(models[[i]]),
                           "residuals" = residuals(models[[i]]),
                           "var" = results$var[i],
                           "year"=dat$year[1:kstart-1],
                           "Y_Rec"=dat$Y_rec[1:kstart-1],
                           "sd"=dat$sd[1:kstart-1])
    resid_df<-rbind(resid_df_temp, resid_df)
}
predicted<-data.frame(predicted, "year"=dat$year[kstart:n_year], "Y_Rec"=dat$Y_rec[kstart:n_year])
#resid_df <- resid_df[complete.cases(resid_df),] # we want no missing values
#resid_df <- resid_df[is.element(resid_df$var, v),]

g2 <- ggplot(resid_df, aes(year,fits)) +
  geom_line() +
  geom_ribbon(aes(x=year, y=fits, ymax=fits+sd, ymin=fits-sd), 
              alpha=0.2)+
  geom_point(data=predicted,aes(x=year,y=Y_Rec),col="red") +
  geom_point(aes(y=Y_Rec))+
  geom_line(data=predicted,aes(x=year,y=pred,col="red"))+
  ggtitle("GAM Fits")+
  facet_wrap(~ var, scale="free_x") +
  theme_bw() +
  
  xlab("Year") + ylab("Fitted") +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"))
g2
#ggsave("figures-yellowtail/GAMFits10.png", height = 8, width = 10)




fitted_univariate_GAMS<-rbind(predicted%>%
  rename(fits=pred)%>%mutate(sd=NA, residuals=NA, Type="predicted")%>%
  select(-ModelID),resid_df%>%select(fits, var, year, Y_Rec, sd, residuals)%>%
    mutate(Type="fitted"))%>%
  rename(fitted_LFO=fits, residual_LFO=residuals)%>%
  right_join(model_fits_LOO)
results_full_gam <-results_arr_LFO5_gam%>%
  left_join(results_arr_LOO_gam%>%select(-ModelID))%>%
  left_join(results_arr_LFO10_gam%>%select(-ModelID))%>%
  mutate(model="GAM")

  
write.csv(arrange(results_full_gam,RMSE_LOO),"univariateGAM_Results.csv")
write.csv(fitted_univariate_GAMS,"univariateGAM_Fits.csv")

######## NON_STATIONARY #######
models <- list()
results <- data.frame()
jstart<-1
n_year <-length(unique(dat$year))

for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small

  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(year,by=", covariates[[i]], ", k = 6)", collapse = " + ")
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
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    #AUC = auc,
    #direction = direction,
    var = var_name
    
  ))
  print(i)
}

results_arr <- arrange(results,RMSE )
# View the results data frame
print(results_arr)

# Forloop through model results to plot residuals for each covariate with convex
# fit

for(i in 1:length(models)) {
  if(i==1) {
    resid_df <- data.frame("fitted" = fitted(models[[i]]),
                           "residuals" = residuals(models[[i]]),
                           "var" = results$var[i],
                           "year"=dat$year,
                           "Y_Rec"=dat$Y_rec)
  } else {
    resid_df <- rbind(resid_df, data.frame("fitted" = fitted(models[[i]]),
                                           "residuals" = residuals(models[[i]]),
                                           "var" = results$var[i],
                           "year"=dat$year,
                           "Y_Rec"=dat$Y_rec))
  }
}

results_arr_LOO_ns <- arrange(results,RMSE )%>%
  rename(AIC_LOO=AIC,RMSE_LOO=RMSE, rsq_LOO=rsq, dev.ex_LOO=dev.ex)
model_fits_LOO<-resid_df%>%
  rename(fitted_LOO=fitted, residuals_loo=residuals)

#### GAMS with leave future out cross valudation ####

n_pred <- 10
train_start<-1
`%notin%` <- Negate(`%in%`)
n_year <-length(unique(dat$year))
n_train <- n_year-n_pred 
predicted<- data.frame()
kstart<-n_train+1
models <- list()
results <- data.frame()
for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small

  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(year,by=", covariates[[i]], ", k = 6)", collapse = " + ")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
predictions <-  numeric(n_pred )

    # Loop over each observation
    for (k in kstart:n_year) {
      train_index <- setdiff(train_start:n_train, k)  # All indices except the k-th
      test_index <- k                 # The k-th index

      # Fit model on n-k observations
      gam_model <- gam(as.formula(formula_str),
                     data = dat[which(dat$year %notin% unique(dat$year)[k:n_year]), ])

      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[k])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[k]), ])
        # re-fit the model

  }
    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     #data = dat)
                     data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
        gam_model2 <- gam(as.formula(formula_str),
                     data = dat)
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec[kstart:n_year] - predictions[kstart:n_year])^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # Extract variable names
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    #AUC = auc,
    #direction = direction,
    var = var_name
    
  ))
    predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[kstart:n_year],
    var = var_name
    
  ))
 
   print(i)
}
#  results_arr_LFO5_ns <- arrange(results,RMSE )%>%
#  rename(AIC_LFO5=AIC,RMSE_LFO5=RMSE, rsq_LFO5=rsq, dev.ex_LFO5=dev.ex)
 results_arr_LFO10_ns <- arrange(results,RMSE )%>%
  rename(AIC_LFO10=AIC,RMSE_LFO10=RMSE, rsq_LFO10=rsq, dev.ex_LFO10=dev.ex)
# View the results data frame
print(results_arr)

# Forloop through model results to plot residuals for each covariate with convex
# fit
resid_df <- data.frame()
for(i in 1:length(models)) {
 
    resid_df_temp <- data.frame("fits" = fitted(models[[i]]),
                           "residuals" = residuals(models[[i]]),
                           "var" = results$var[i],
                           "year"=dat$year[1:kstart-1],
                           "Y_Rec"=dat$Y_rec[1:kstart-1],
                           "sd"=dat$sd[1:kstart-1])
    resid_df<-rbind(resid_df_temp, resid_df)
}
predicted<-data.frame(predicted, "year"=dat$year[kstart:n_year], "Y_Rec"=dat$Y_rec[kstart:n_year])
#resid_df <- resid_df[complete.cases(resid_df),] # we want no missing values
#resid_df <- resid_df[is.element(resid_df$var, v),]

fitted_univariate_GAMS<-rbind(predicted%>%
  rename(fits=pred)%>%mutate(sd=NA, residuals=NA, Type="predicted")%>%
  select(-ModelID),resid_df%>%select(fits, var, year, Y_Rec, sd, residuals)%>%
    mutate(Type="fitted"))%>%
  rename(fitted_LFO=fits, residual_LFO=residuals)%>%
  right_join(model_fits_LOO)
results_full_ns <-results_arr_LFO5_ns%>%
  left_join(results_arr_LOO_ns%>%select(-ModelID))%>%
  left_join(results_arr_LFO10_ns%>%select(-ModelID))%>%
  mutate(model="NS")

  
write.csv(arrange(results_full_ns,RMSE_LOO),"univariateNS_Results.csv")
write.csv(fitted_univariate_GAMS,"univariateGAM_Fits.csv")
######## ROLLING WINDOW #######
models <- list()
results <- data.frame()
jstart<-1
n_year <-10
length_window<-10
firstyear<-1993:(2014-length_window)

for(k in jstart:length(firstyear)){
  lastyear=firstyear[k]+length_window
  datwindow <- dat%>%filter(year>=firstyear[k]&year<=lastyear)
  results_temp <- data.frame()
for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small

  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", covariates[[i]], ", k = 3)")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
    predictions <- numeric(nrow(datwindow))
    n_year <- length(unique(datwindow$year))
    # Loop over each observation
    for (j in jstart:n_year) {
      train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index

      # Fit model on n-1 observations
      gam_model <- gam(as.formula(formula_str),
                       # weights = number_cwt_estimated,
                       data = datwindow[which(datwindow$year != unique(datwindow$year)[j]), ])

      # Predict the excluded observation
      predictions[which(datwindow$year == unique(datwindow$year)[j])] <- predict(gam_model, newdata = datwindow[which(datwindow$year == unique(datwindow$year)[j]), ])
    }

    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     data = datwindow)

  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((datwindow$Y_rec - predictions)^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # Extract variable names
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results_temp<- rbind(results_temp,data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    firstyear=firstyear[k],
    lastyear=lastyear,
    #AUC = auc,
    #direction = direction,
    var = var_name
  ))

         
  print(i)
}

       results<- rbind(results,results_temp%>%
                   mutate(min_RMSE=min(RMSE),min_AIC=min(AIC))%>%
                    mutate(RMSE_diff=RMSE-min_RMSE,AIC_diff=AIC-min_AIC)%>%
                    mutate(rel_RMSE=exp(-0.5 * RMSE_diff), rel_AIC=exp(-0.5 * AIC_diff))%>%
                     mutate(sum_RMSE=sum(rel_RMSE), sum_AIC=sum(rel_AIC))%>%
                    mutate(RMSE_weight=rel_RMSE/sum_RMSE, AIC_weight=rel_AIC/sum_AIC))   
}

results_rolling<-results%>%mutate(range=paste(firstyear, "-", lastyear))

rollingplot<- ggplot(results_rolling, aes(as.factor(range), var, fill= AIC_weight)) + 
  xlab("Time Period")+
  scale_fill_gradient(low = "white", high = "Darkgreen") +
  ylab("Oceanographic Conditions")+
    ggtitle("Backwards Selection")+
  geom_tile()+
  theme_bw()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
rollingplot
ggsave("figures-yellowtail/rolling_gam.png", height = 8, width = 13)



#### GAMS with Adding 1 - year at a time ####

models <- list()
results <- data.frame()
jstart<-1
n_year <-10
length_window<-10
firstyear<-1993
#rm(firstyear)


for(k in jstart:11){
  lastyear<- firstyear+10+k
  datwindow <- dat%>%filter(year>=firstyear&year<=lastyear)
  results_temp <- data.frame()
for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small

  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", covariates[[i]], ", k = 3)")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
    predictions <- numeric(nrow(datwindow))
    n_year <- length(unique(datwindow$year))
    # Loop over each observation
    for (j in jstart:n_year) {
      train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index

      # Fit model on n-1 observations
      gam_model <- gam(as.formula(formula_str),
                       # weights = number_cwt_estimated,
                       data = datwindow[which(datwindow$year != unique(datwindow$year)[j]), ])

      # Predict the excluded observation
      predictions[which(datwindow$year == unique(datwindow$year)[j])] <- predict(gam_model, newdata = datwindow[which(datwindow$year == unique(datwindow$year)[j]), ])
    }

    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     data = datwindow)

  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((datwindow$Y_rec - predictions)^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # Extract variable names
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results_temp<- rbind(results_temp,data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    lastyear=lastyear,
    firstyear=firstyear,
    #AUC = auc,
    #direction = direction,
    var = var_name
  ))

         
  print(i)
}

       results<- rbind(results,results_temp%>%
                   mutate(min_RMSE=min(RMSE),min_AIC=min(AIC))%>%
                    mutate(RMSE_diff=RMSE-min_RMSE,AIC_diff=AIC-min_AIC)%>%
                    mutate(rel_RMSE=exp(-0.5 * RMSE_diff), rel_AIC=exp(-0.5 * AIC_diff))%>%
                     mutate(sum_RMSE=sum(rel_RMSE), sum_AIC=sum(rel_AIC))%>%
                    mutate(RMSE_weight=rel_RMSE/sum_RMSE, AIC_weight=rel_AIC/sum_AIC))   
}

results_forewards<-results%>%mutate(range=paste(firstyear, "-", lastyear))

forewards_plot <- ggplot(results_forewards, aes(as.factor(range), var, fill= AIC_weight)) + 
  xlab("Time Period")+
  ggtitle("Foreward Selection")+
  scale_fill_gradient(low = "white", high = "Darkgreen") +
  ylab("Oceanographic Conditions")+
  geom_tile()+
  theme_bw()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
forewards_plot
ggsave("figures-yellowtail/forewards_gam.png", height = 8, width = 13)
#### GAMS with Adding 1 - year at a time in reverse ####

models <- list()
results <- data.frame()
jstart<-1
n_year <-10
length_window<-10
lastyear<-2014
#rm(firstyear)


for(k in jstart:11){
  firstyear<- lastyear-10-k
  datwindow <- dat%>%filter(year>=firstyear&year<=lastyear)
  results_temp <- data.frame()
for (i in seq_along(covariates)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small

  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", covariates[[i]], ", k = 3)")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
    predictions <- numeric(nrow(datwindow))
    n_year <- length(unique(datwindow$year))
    # Loop over each observation
    for (j in jstart:n_year) {
      train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index

      # Fit model on n-1 observations
      gam_model <- gam(as.formula(formula_str),
                       # weights = number_cwt_estimated,
                       data = datwindow[which(datwindow$year != unique(datwindow$year)[j]), ])

      # Predict the excluded observation
      predictions[which(datwindow$year == unique(datwindow$year)[j])] <- predict(gam_model, newdata = datwindow[which(datwindow$year == unique(datwindow$year)[j]), ])
    }

    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     data = datwindow)

  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((datwindow$Y_rec - predictions)^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # Extract variable names
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results_temp<- rbind(results_temp,data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    firstyear=firstyear,
    lastyear=lastyear,
    var = var_name
  ))

         
  print(i)
}

       results<- rbind(results,results_temp%>%
                   mutate(min_RMSE=min(RMSE),min_AIC=min(AIC))%>%
                    mutate(RMSE_diff=RMSE-min_RMSE,AIC_diff=AIC-min_AIC)%>%
                    mutate(rel_RMSE=exp(-0.5 * RMSE_diff), rel_AIC=exp(-0.5 * AIC_diff))%>%
                     mutate(sum_RMSE=sum(rel_RMSE), sum_AIC=sum(rel_AIC))%>%
                    mutate(RMSE_weight=rel_RMSE/sum_RMSE, AIC_weight=rel_AIC/sum_AIC))   
}

results_backwards<-results%>%mutate(range=paste(firstyear, "-", lastyear))

backwards_plot <- ggplot(results_backwards, aes(as.factor(range), var, fill= AIC_weight)) + 
  xlab("Time Period")+
  scale_fill_gradient(low = "white", high = "Darkgreen") +
  ylab("Oceanographic Conditions")+
    ggtitle("Backwards Selection")+
  geom_tile()+
  theme_bw()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
backwards_plot

ggsave("figures-yellowtail/backwards_gam.png", height = 8, width = 13)

