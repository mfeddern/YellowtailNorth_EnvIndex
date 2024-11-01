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

cross_validation <- TRUE

models <- list()
results <- data.frame()
jstart<-1
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
  } else {
    # just fit the model once
    gam_model <- gam(as.formula(formula_str),
                     data = dat)
    predictions <- as.numeric(predict(gam_model, type="response")) # probability

    indx <- sort(c(dat[,covariates[i]][[1]]), index.return=T)
    pred_df <- data.frame(x = seq(min(dat[,covariates[i]]), max(dat[,covariates[i]]), length.out=300))
    names(pred_df) <- covariates[i]
    sec_deriv <- predict(gam_model, newdata=pred_df, type = "terms", derivatives = 2)
    # sign_changes <- diff(sign(sec_deriv))  # This will be nonzero if there is a change in sign
    # inflection_pt <- which(sign_changes!=0)
    # mean_deriv <- c(mean(sec_deriv[1:inflection_pt]), mean(sec_deriv[(inflection_pt+2):length(sec_deriv)]))
    # direction <- "down"
    # if(mean_deriv[1] < mean_deriv[2]) direction <- "up"
    # easiest way is just to evaluate sign
    lm_fit <- lm(sec_deriv~seq(1,length(sec_deriv)) + I(seq(1,length(sec_deriv))^2))
    #direction <- "up"
        direction <- "up"
    if(coef(lm_fit)[3] < 0) direction <- "down"
  }

  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec - predictions)^2, na.rm=T))
  r2<-summary(gam_model)$r.sq
  dev.expl<-summary(gam_model)$dev.expl
  # We can instead calculate AUC -- but this involves expanding our compact data frame
  # for(j in 1:nrow(dat)) { # slow version
  #   lived <- dat[rep(j, dat$total_return[j]),]
  #   lived$actual <- 1
  #   died <- dat[rep(j, dat$not_return[j]),]
  #   died$actual <- 0
  #   all_fish <- rbind(lived, died)
  #
  #   if(j == 1) {
  #     expanded_df <- all_fish
  #   } else {
  #     expanded_df <- rbind(expanded_df, all_fish)
  #   }
  # }
  # expanded_df <- dat %>%
    # Repeat rows by 'total_return' for successes
  #   uncount(weights = total_return, .remove = FALSE) %>%
  #   mutate(actual = 1) %>%
    # Combine with rows repeated by 'not_return' for failures
  #   bind_rows(
  #      dat %>%
  #        uncount(weights = not_return, .remove = FALSE) %>%
  #        mutate(actual = 0)
  #    )
#  expanded_predictions <- as.numeric(predict(gam_model, newdata = expanded_df, type="response"))
#  pred_ROCR <- ROCR::prediction(expanded_predictions, expanded_df$actual) #this is throwing an error, from ROcR package right?
#  auc_ROCR <- ROCR::performance(pred_ROCR, measure = "auc")
#  auc <- auc_ROCR@y.values[[1]]

  # Extract variable names
  var_name <- covariates[i]

  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = round(rmse,2),
    rsq=round(r2,2),
    dev.ex=round(dev.expl,4),
    #AUC = auc,
    #direction = direction,
    var = var_name
    
  ))
  print(i)
}


# filter for frowns only

results <- arrange(results, RMSE)
# View the results data frame
print(results)

# Forloop through model results to plot residuals for each covariate with convex
# fit

for(i in 1:length(models)) {
  if(i==1) {
    resid_df <- data.frame("fitted" = fitted(models[[i]]),
                           "residuals" = residuals(models[[i]]),
                           "var" = results$var[i])
  } else {
    resid_df <- rbind(resid_df, data.frame("fitted" = fitted(models[[i]]),
                                           "residuals" = residuals(models[[i]]),
                                           "var" = results$var[i]))
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


##### Kristin's Code #####
min_year <- 1994# first year of complete covariate data
max_year <- 2014 # maximum year of data to be forecast
envir_data2<-envir_data%>%
  cbind(time=dat$year)%>%
  tidyr::pivot_longer(-time)%>%
  rename(cov=value)%>%
  dplyr::group_by(name) %>%
      dplyr::mutate(est_z = (cov - mean(cov, na.rm=T))/sd(cov, na.rm=T))%>%
  select(-cov)%>%
  rename(cov=est_z)%>%
  tidyr::pivot_wider(names_from="name",
                             values_from = "cov")

dat_2<-df%>%select(year, Y_rec, sd)%>%mutate(Label=paste("Main_RecrDev_",year))%>%
  rename(time=year, dev=Y_rec,Parm_StDev=sd)%>%
  dplyr::filter(time <= max_year, time >= min_year)

predictors<-read.csv("envir_data2.csv")%>%select(-X)%>%
  select(time, DDpre,DDegg,DDlarv,DDpjuv,PDOlarv,PDOpjuv,BeutiTUMIpjuv,
         BeutiSTIpjuv,Tpart)
predictors2<-read.csv("envir_data2.csv")%>%select(-X)%>%
  select(time, ONIpjuv,CutiSTIpjuv,CutiTUMIpjuv,Tcop,DDben,MLDpart,MLDlarv,MLDpjuv,
         CSTlarv,CSTpjuv, LSTlarv,LSTpjuv,HCIlarv,HCIpjuv,ONIpre)

gam_example <- univariate_forecast(dat_2,
  as.data.frame(predictors2),
  model_type = "gam", # can be 'lm' or 'gam'
  n_forecast = 3, # how many years of data to use in testing
  n_years_ahead = 1, # how many time steps ahead to forecast
  max_vars = 3)

gam_example$coefs[[1]]
gam_example$marginal[[1]]

lm_example <- univariate_forecast(dat_2,
  as.data.frame(predictors2),
  model_type = "gam", # can be 'lm' or 'gam'
  n_forecast = 3, # how many years of data to use in testing
  n_years_ahead = 1, # how many time steps ahead to forecast
  max_vars = 3)

gam_example$coefs[[1]]




response <- data.frame(time = 1:40, dev = rnorm(40))

predictors <- matrix(rnorm(4 * 40), ncol = 4)
#' #colnames(predictors) = paste0("X",1:ncol(predictors))
predictors <- as.data.frame(predictors)
predictors$time <- 1:40
lm_example <- univariate_forecast(response,
  predictors,
  model_type = "lm", # can be 'lm' or 'gam'
  n_forecast = 10, # how many years of data to use in testing
  n_years_ahead = 1, # how many time steps ahead to forecast
  max_vars = 3) 

df2 = expand.grid(response = "yellowtail",
                 predictors = c("glorys"),
                 model = c("lm","gam"),
                 n_forecast = 3,
                 years_ahead = 1)


for(i in 1:nrow(df2)) {

  # Load response data
  if(df2$response[i]=="yellowtail") {
    #data("recruit_dev_hake_2021")
    rec_devs = dat_2     # Standardize
    #rec_devs$dev = scale(rec_devs$dev)
  }

  # truncate rec_devs and dat to be < max_year
  rec_devs <- dplyr::filter(rec_devs, time <= max_year, time >= min_year)
  #dat <- dplyr::filter(dat, time <= (max_year-1))
  if(df2$predictors[i]=="glorys") {
    dat =envir_data2 #dplyr::rename(new_envir_data,cov=value,time=year)
    max_vars = 3
    # pivot wider to create matrix (years on rows)
    #dat = tidyr::pivot_wider(dat[,c("name","time","cov")],
     #                        names_from="name",
       #                      values_from = "cov")
  }



  if(df2$model[i] %in% c("lm","gam")) {
      forecast = predRecruit::univariate_forecast(response = rec_devs,
                                   predictors = envir_data2,
                                   model_type = df2$model[i],
                                   n_forecast = df2$n_forecast[i],
                                   n_years_ahead = df2$years_ahead[i],
                                   max_vars = 3)
  }
 
  #print(i)
  #saveRDS(forecast, file = paste0("output/2018/",df2$predictors[i],"_",df2$model[i],"_",df2$years_ahead[i],"step.rds"))
}
