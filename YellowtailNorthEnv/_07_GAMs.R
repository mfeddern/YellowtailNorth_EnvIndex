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
data_years = 1994:2014
dat = df %>% 
  dplyr::select(!any_of( c('ZOOpjuv', 'ZOOben')))  %>% 
  filter(year %in% data_years)
envir_data = dat %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','year','Y_rec','ZOOpjuv','ZOOben')))%>%  # drop terms not in model statement
    dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI','ONIlarv')))%>% 
    mutate_all(~ scale(.)) # drop terms not in model statement

head(envir_data)
dim(envir_data)

summary(gam( BeutiSTIpjuv~ ONIpjuv, data=envir_data))$r.sq
# Vector of marine covariates for univariate GAM loop

covariates <- names(envir_data)#[-which(names(ocean) %in% "year")]

#### GAMS with Leave One OUt Cross Validation ####

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
  smooth_terms <- paste("s(", covariates[[i]], ", k = 3)")
  formula_str <- paste("Y_rec ~ ", smooth_terms)

  if(cross_validation) {
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
  } 
  
  else {
    # just fit the model once
    gam_model <- gam(as.formula(formula_str),
                     data = dat)
    predictions <- as.numeric(predict(gam_model, type="response")) # probability
    indx <- sort(c(dat[,covariates[i]][[1]]), index.return=T)
    pred_df <- data.frame(x = seq(min(dat[,covariates[i]]), max(dat[,covariates[i]]), length.out=300))
    names(pred_df) <- covariates[i]
    sec_deriv <- predict(gam_model, newdata=pred_df, type = "terms", derivatives = 2)
    lm_fit <- lm(sec_deriv~seq(1,length(sec_deriv)) + I(seq(1,length(sec_deriv))^2))
    #direction <- "up"
        direction <- "up"
    if(coef(lm_fit)[3] < 0) direction <- "down"
  }

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

results <- arrange(results,RMSE )
# View the results data frame
print(results)

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

resid_df <- resid_df[complete.cases(resid_df),] # we want no missing values
resid_df <- resid_df[is.element(resid_df$var, v),]

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
#### GAMS with leave future out cross valudation ####

n_pred <- 4
train_start<-1

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
                       # weights = number_cwt_estimated,
                     data = dat[which(dat$year %notin% unique(dat$year)[k:n_year]), ])

      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[k])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[k]), ])
        # re-fit the model

  }
    # re-fit the model
    gam_model <- gam(as.formula(formula_str),
                     data = dat[which(dat$year %notin% unique(dat$year)[kstart:n_year]), ])
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
  
results <- arrange(results,RMSE )
# View the results data frame
print(results)

# Forloop through model results to plot residuals for each covariate with convex
# fit

for(i in 1:length(models)) {
  if(i==1) {
    resid_df <- data.frame("fitted" = fitted(models[[i]]),
                           "residuals" = residuals(models[[i]]),
                           "var" = results$var[i],
                           "year"=dat$year[1:kstart-1],
                           "Y_Rec"=dat$Y_rec[1:kstart-1],
                           "sd"=dat$sd[1:kstart-1])
  } else {
    resid_df <- rbind(resid_df, data.frame("fitted" = fitted(models[[i]]),
                                           "residuals" = residuals(models[[i]]),
                                           "var" = results$var[i],
                           "year"=dat$year[1:kstart-1],
                           "Y_Rec"=dat$Y_rec[1:kstart-1],
                           "sd"=dat$sd[1:kstart-1]))
  }
}
predicted<-data.frame(predicted, "year"=dat$year[kstart:n_year], "Y_Rec"=dat$Y_rec[kstart:n_year])
#resid_df <- resid_df[complete.cases(resid_df),] # we want no missing values
#resid_df <- resid_df[is.element(resid_df$var, v),]

g2 <- ggplot(resid_df, aes(year,fitted)) +
  geom_line() +
  geom_point(data=predicted,aes(x=year,y=Y_Rec),col="red") +
  geom_point(aes(y=Y_Rec))+
  geom_line(data=predicted,aes(x=year,y=pred,col="red"))+
  facet_wrap(~ var, scale="free_x") +
  theme_bw() +
  xlab("Year") + ylab("Fitted") +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"))
g2 
