library(dplyr)
library(tidyr)
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

df = data.frame(read.csv("data-yellowtail/02_DATA_Combined_glorys_yellowtail2025_unexpanded.csv"))
dfa = data.frame(read.csv("data-yellowtail/dfa_trend.csv"))
data_years = 1994:2019
dat = df %>% 
  dplyr::select(!any_of( c('ZOOpjuv', 'ZOOben')))  %>% 
  filter(year %in% data_years)%>%
  left_join(dfa)
envir_data = dat %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','year','Y_rec','ZOOpjuv','ZOOben',"dfa")))%>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI','ONIlarv',"X","dfa","Type")))%>% 
  mutate_all(~ scale(.)) # drop terms not in model statement
dat<-cbind(dat%>%select(year, Y_rec, sd),envir_data)

summary(gam( BeutiSTIpjuv~ ONIpjuv, data=envir_data))$r.sq
# Vector of marine covariates for univariate GAM loop

covariates <- names(envir_data)#[-which(names(ocean) %in% "year")]
#### Evaluating Correlations #####
envir_data <- envir_data[complete.cases(envir_data), ]
M = data.frame(cor(envir_data))
M <- tibble::rownames_to_column(M, "xvar")
M<-M%>%pivot_longer(!xvar, names_to = 'yvar', values_to = 'corr')
uncorr<-M%>%filter(abs(corr)<0.4)
corrplotdat<-cor(envir_data)
unique(uncorr$xvar)
corrplot.mixed(corrplotdat)

corrvar<-M%>%filter(abs(corr)>0.4)%>%filter(xvar!=yvar)
corrvar2<-corrvar%>%select(-corr)
combinations_to_omit<-list()

for(i in 1:nrow(corrvar2)){
  combinations_to_omit[[i]]<-c(as.character(corrvar2[i, ]))
}

# One decision: what's the maximum number of covariates that any model can have?
combinations <- lapply(1:4, function(i) {
  combn(covariates, i, simplify = FALSE)
})
combinations <- unlist(combinations, recursive = FALSE)
combinations <- unique(lapply(combinations, function(x) sort(x)))

length(combinations)
colnames(envir_data) 
# Function to check if a combination is a partial match of any combination to omit
is_partial_match <- function(comb, omit_list) {
  any(sapply(omit_list, function(omit) all(omit %in% comb)))
}

# Remove combinations that are partial matches (but not exact matches)
combinations <- combinations[!sapply(combinations, is_partial_match, omit_list = combinations_to_omit)]
# Check the length of remaining combinations
length(combinations)
# Sort each combination and remove duplicates


#### LOO ####


cross_validation <-TRUE
models <- list()
results <- data.frame()
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
  print(i)
}


results_arr_RMSE_LOO <- arrange(results,RMSE)
results_arr_LOO_rmse <- arrange(results,RMSE)
results_arr_LOO_AIC <- arrange(results,AIC)

# View the results data frame
print(results_arr_RMSE_LOO)
print(results_arr_devex)
print(results_arr_LOO_AIC)
print(results_arr_rsq)
# Fit a model that includes multiple variables from the top model
combosLOO = c(results_arr_RMSE_LOO$var1[1], results_arr_RMSE_LOO$var2[1], results_arr_RMSE_LOO$var3[1])
smooth_terms <- paste("s(", combosLOO, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", smooth_terms)
gamloo <- gam(as.formula(formula_str),
                 #family = binomial(),
                 #weights = number_cwt_estimated,
                 data = dat)

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

#### LFO5 ####

n_pred <- 5
train_start<-1
`%notin%` <- Negate(`%in%`)
n_year <-length(unique(dat$year))
n_train <- n_year-n_pred 
predicted<- data.frame()
kstart<-n_train+1
models <- list()
results <- data.frame()
for (i in seq_along(combinations)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  
  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ")
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
  gam_model_train <- gam(as.formula(formula_str),
                         #data = dat)
                         data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
  gam_model_full <- gam(as.formula(formula_str),
                        data = dat)
  gam_model_pred <- gam(as.formula(formula_str),
                        #data = dat)
                        data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec[kstart:n_year] - predictions[kstart:n_year])^2, na.rm=T))
  r2_full<-summary(gam_model_full)$r.sq
  r2_train<-summary(gam_model_train)$r.sq
  r2_pred<-summary(gam_model_pred)$r.sq
  dev.expl_full<-summary(gam_model_full)$dev.expl
  dev.expl_train<-summary(gam_model_train)$dev.expl
  dev.expl_pred<-summary(gam_model_pred)$dev.expl
  # Extract variable names
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Extract variable names
  padded_vars <- c(var_names, rep(NA, 4 - length(var_names)))
  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq_full=round(r2_full,2),
    rsq_train=round(r2_train,2),
    rsq_pred=round(r2_pred,2),
    dev.ex_full=round(dev.expl_full,4),
    dev.ex_train=round(dev.expl_train,4),
    dev.ex_pred=round(dev.expl_pred,4),
    #AUC = auc,
    #direction = direction,
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var4 = padded_vars[4]
    
    
  ))
  predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[kstart:n_year],
    year=unique(dat$year)[kstart:n_year],
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var3 = padded_vars[4]
    
    
  ))
  
  print(i)
}

results_arr_LFO5_rmse <- arrange(results,RMSE)
results_arr_LFO5_AIC <- arrange(results,AIC)
results_arr_LFO5_devex <- arrange(results,desc(dev.ex_pred))
# View the results data frame
print(results_arr_LFO5_rmse)
predicted_last5<-predicted


# Fit a model that includes multiple variables from the top model
combosLOO = c(results_arr_LFO5_rmse$var1[1], results_arr_LFO5_rmse$var2[1], results_arr_LFO5_rmse$var3[1], results_arr_LFO5_rmse$var4[1])
smooth_terms <- paste("s(", combosLOO, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", smooth_terms)
gamloo <- gam(as.formula(formula_str),
              #family = binomial(),
              #weights = number_cwt_estimated,
              data = dat)

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




#### LFO10 ####

n_pred <-10
train_start<-1
`%notin%` <- Negate(`%in%`)
n_year <-length(unique(dat$year))
n_train <- n_year-n_pred 
predicted<- data.frame()
kstart<-n_train+1
models <- list()
results <- data.frame()
for (i in seq_along(combinations)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  
  # smooth on total release isn't needed when we're not working with jack rate
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ")
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
  gam_model_train <- gam(as.formula(formula_str),
                   #data = dat)
                   data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
  gam_model_full <- gam(as.formula(formula_str),
                    data = dat)
  gam_model_pred <- gam(as.formula(formula_str),
                        #data = dat)
                        data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec[kstart:n_year] - predictions[kstart:n_year])^2, na.rm=T))
  r2_full<-summary(gam_model_full)$r.sq
  r2_train<-summary(gam_model_train)$r.sq
  r2_pred<-summary(gam_model_pred)$r.sq
  dev.expl_full<-summary(gam_model_full)$dev.expl
  dev.expl_train<-summary(gam_model_train)$dev.expl
  dev.expl_pred<-summary(gam_model_pred)$dev.expl
  # Extract variable names
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Extract variable names
  padded_vars <- c(var_names, rep(NA, 4 - length(var_names)))
  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,3),
    rsq_full=round(r2_full,2),
    rsq_train=round(r2_train,2),
    rsq_pred=round(r2_pred,2),
    dev.ex_full=round(dev.expl_full,4),
    dev.ex_train=round(dev.expl_train,4),
    dev.ex_pred=round(dev.expl_pred,4),
    #AUC = auc,
    #direction = direction,
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var4 = padded_vars[4]
    
    
  ))
  predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[kstart:n_year],
    year=unique(dat$year)[kstart:n_year],
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var3 = padded_vars[4]
    
    
  ))
  
  print(i)
}

results_arr_LFO10_rmse <- arrange(results,RMSE)
results_arr_LFO10_rmse <- arrange(results,AIC)
results_arr_LFO10_devex <- arrange(results,desc(rsq_pred))
# View the results data frame
print(results_arr_LFO10_rmse)

predicted_last10<-predicted

# Fit a model that includes multiple variables from the top model
combosLOO = c(results_arr_LFO10_rmse$var1[1], results_arr_LFO10_rmse$var2[1])
smooth_terms <- paste("s(", combosLOO, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", smooth_terms)
gamloo <- gam(as.formula(formula_str),
              #family = binomial(),
              #weights = number_cwt_estimated,
              data = dat)

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

#### Deviance explaned ####


ggplot(results_arr_LFO5_rmse%>%mutate(rank=seq(from=1,to=nrow(results_arr_LFO5_rmse),by=1)), aes(x=RMSE, y=dev.ex_pred,col=rank))+
  geom_point()+
  theme_bw()

ggplot(results_arr_LFO5_rmse%>%mutate(rank=seq(from=1,to=nrow(results_arr_LFO5_rmse),by=1)), aes(x=RMSE, y=rsq_pred,col=rank))+
  geom_point()+
  theme_bw()

ggplot(results_arr_LFO10_rmse%>%mutate(rank=seq(from=1,to=nrow(results_arr_LFO10_rmse),by=1)), aes(x=RMSE, y=dev.ex_pred,col=rank))+
  geom_point()+
  theme_bw()

ggplot(results_arr_LFO10_rmse%>%mutate(rank=seq(from=1,to=nrow(results_arr_LFO10_rmse),by=1)), aes(x=RMSE, y=rsq_pred,col=rank))+
  geom_point()+
  theme_bw()

ggplot(results_arr_RMSE_LOO%>%mutate(rank=seq(from=1,to=nrow(results_arr_RMSE_LOO),by=1)), aes(x=RMSE, y=rsq,col=rank))+
  geom_point()+
  theme_bw()

ggplot(results_arr_RMSE_LOO%>%mutate(rank=seq(from=1,to=nrow(results_arr_RMSE_LOO),by=1)), aes(x=RMSE, y=dev.ex,col=rank))+
  geom_point()+
  theme_bw()
#### Top 3 models ####
predictions<-predicted_last5%>%filter(ModelID==results_arr_LFO5_rmse[1,1]|
                                        ModelID==results_arr_LFO10_devex[1,1]|
                                        ModelID==results_arr_RMSE_LOO[1,1]|
                                        ModelID==results_arr_LFO5_AIC[1,1])%>%
  mutate(rankedBest=ifelse(ModelID==results_arr_LFO5_rmse[1,1],"LFO-5",
                           ifelse(ModelID==results_arr_LFO10_rmse[1,1],"rsq_pred",
                                  ifelse(ModelID==results_arr_RMSE_LOO[1,1],"LOO","AIC"))))
ensemble<-predictions%>%group_by(year)%>%summarise(ensemble_index=mean(pred))

ggplot(data=df,aes(y=Y_rec, x=year))+
  geom_ribbon(aes(ymax=Y_rec+sd,ymin=Y_rec-sd,alpha=0.2),fill='gray')+
  geom_point(data=ensemble,aes(x=year,y=ensemble_index))+
  geom_line(data=ensemble,aes(x=year,y=ensemble_index))+
  geom_point(col='grey')+
  geom_line(col='grey')+
  theme_bw() +
  xlim(c(1993,2021))+
  geom_point(data=predictions,aes( x=year, y=pred,col=rankedBest))+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



#### Top 3 models ####
predictions<-predicted_last10%>%filter(ModelID==results_arr_LFO5_rmse[1,1]|
                                            ModelID==results_arr_LFO10_devex[1,1]|
                                            ModelID==results_arr_RMSE_LOO[1,1]|
                                            ModelID==results_arr_LFO5_AIC[1,1])%>%
  mutate(rankedBest=ifelse(ModelID==results_arr_LFO5_rmse[1,1],"LFO-5",
                           ifelse(ModelID==results_arr_LFO10_rmse[1,1],"rsq_pred",
                                  ifelse(ModelID==results_arr_RMSE_LOO[1,1],"LOO","AIC"))))
ensemble<-predictions%>%group_by(year)%>%summarise(ensemble_index=mean(pred))

ggplot(data=df,aes(y=Y_rec, x=year))+
  geom_ribbon(aes(ymax=Y_rec+sd,ymin=Y_rec-sd,alpha=0.2),fill='gray')+
  geom_point(data=ensemble,aes(x=year,y=ensemble_index))+
  geom_line(data=ensemble,aes(x=year,y=ensemble_index))+
  geom_point(col='grey')+
  geom_line(col='grey')+
  theme_bw() +
  xlim(c(1993,2021))+
  geom_point(data=predictions,aes( x=year, y=pred,col=rankedBest))+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


# CREATE TABLE OF RMSE/AIC------------------------------------

table_5 <- dplyr::arrange(marginals, total_rmse) %>%
  dplyr::select(cov, total_rmse, total_aic)%>%
  mutate(cv="LFO_5",model="GAM")

# Fit a model that includes multiple variables from the top model
combos = c(table$cov[1], table$cov[2], table$cov[3], table$cov[4], table$cov[5])
smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", smooth_terms)
gam_model <- gam(as.formula(formula_str),
                 #family = binomial(),
                 #weights = number_cwt_estimated,
                 data = dat)













