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
dat<-cbind(dat%>%select(year, Y_rec, sd),envir_data)

summary(gam( BeutiSTIpjuv~ ONIpjuv, data=envir_data))$r.sq
# Vector of marine covariates for univariate GAM loop

covariates <- names(envir_data)#[-which(names(ocean) %in% "year")]

#### GAMS with Leave One OUt Cross Validation ####
# One decision: what's the maximum number of covariates that any model can have?
combinations <- lapply(1:3, function(i) {
  combn(covariates, i, simplify = FALSE)
})
combinations <- unlist(combinations, recursive = FALSE)
combinations <- unique(lapply(combinations, function(x) sort(x)))

length(combinations)

combinations_to_omit <- list(
  c("DDlarv", "Tpart"),
  c("CutiSTIpjuv", "BeutiSTIpjuv"),
  c("DDlarv", "DDpjuv"), 
  c("DDlarv", "HCIlarv"), 
  c("DDlarv", "PDOlarv"), 
  c("DDlarv", "Tpart"), 
  c("DDpjuv", "HCIlarv"), 
  c("DDpjuv", "HCIpjuv"), 
  c("DDpjuv", "PDOlarv"), 
  c('DDpjuv', "PDOpjuv"), 
  c("DDpjuv", "Tpart"), 
  c("DDpre", "DDegg"), 
  c("DDpre", "DDlarv"), 
  c("DDpre", "HCIlarv"),
  c("DDpre", "PDOlarv"), 
  c("DDpre", "Tcop"), 
  c("DDpre", "Tpart"), 
  c("HCIlarv", "PDOlarv"), 
  c("HCIpjuv", "PDOpjuv"), 
  c("MLDpart", "MLDlarv"), 
  c("Tpart", "HCIlarv"), 
  c("Tpart", "PDOlarv"),
  c("DDpre", "dfa"),
  c("DDegg" , "dfa"),
  c("DDlarv", "dfa"),
  c("DDben" , "dfa"),
  c("Tcop", "dfa"),
  c("Tpart", "dfa"),
  c("MLDpart", "dfa"),
  c( "CSTpjuv", "dfa"),
  c("LSTpjuv", "dfa"),
  c("HCIlarv", "dfa"),
  c("HCIpjuv", "dfa"),
  c("ONIpre", "dfa"),
  c("ONIpjuv", "dfa"),
  c("PDOlarv", "dfa"),
  c("PDOpjuv", "dfa"),
  c("CutiTUMIpjuv", "dfa"),
  c("BeutiTUMIpjuv", "dfa"),
  c("BeutiSTIpjuv", "dfa")) 

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
  padded_vars <- c(var_names, rep(NA, 3 - length(var_names)))

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
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3]
    
  ))
  print(i)
}


# Create a baseline model with only the intercept

    predictions <- numeric(nrow(dat))
    n_year <- length(unique(dat$year))

    # Loop over each observation
    for (j in 1:n_year) {
      train_index <- setdiff(1:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index

      # Fit model on n-1 observations
      gam_model <- gam(Y_rec ~ 1,
                       data = dat[which(dat$year != unique(dat$year)[j]), ])

      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    }

baseline_rmse  <- sqrt(mean((dat$Y_rec - (predictions))^2, na.rm=T))

# filter for frowns only
results_gam_loo<-results%>%
  mutate(RMSE_improv=RMSE/baseline_rmse, cv="LOO", model="GAM")
results_arr_RMSE <- arrange(results_gam_loo,RMSE )
# View the results data frame
print(results_arr_RMSE)
results_arr_AIC <- arrange(results_gam_loo,AIC )
# View the results data frame
print(results_arr_AIC)

# calculate the baseline RMSE

# we can calculate the marginal improvement for each covariate
results$n_cov <- ifelse(!is.na(results$var1), 1, 0) + ifelse(!is.na(results$var2), 1, 0) +
  ifelse(!is.na(results$var3), 1, 0)
marginals <- data.frame(cov = covariates, "rmse_01" = NA, "rmse_12" = NA, "rmse_23" = NA,
                        "aic_01" = NA, "aic_12" = NA, "aic_23" = NA)
for(i in 1:length(covariates)) {
  sub <- dplyr::filter(results, n_cov == 1,
                       var1 == covariates[i])
  marginals$rmse_01[i] <- sub$RMSE / baseline_rmse
  marginals$aic_01[i] <- AIC(gam_model) - sub$AIC

  # next look at all values of models that have 2 covariates and include this model
  sub1 <- dplyr::filter(results, n_cov == 1)
  sub2 <- dplyr::filter(results, n_cov == 2) %>%
    dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0)) %>%
    dplyr::filter(keep == 1) %>% dplyr::select(-keep)
  # loop over every variable in sub2, and find the simpler model in sub1 that just represents the single covariate
  sub2$rmse_diff <- 0
  sub2$AIC_diff <- 0
  for(j in 1:nrow(sub2)) {
    vars <- sub2[j,c("var1", "var2")]
    vars <- vars[which(vars != covariates[i])]
    indx <- which(sub1$var1 == as.character(vars))
    sub2$rmse_diff[j] <- sub2$RMSE[j] / sub1$RMSE[indx]
    sub2$AIC_diff[j] <- sub1$AIC[indx] - sub2$AIC[j]
  }

  # Finally compare models with 3 covariates to models with 2 covariates
  sub2_all <- dplyr::filter(results, n_cov == 2)
  # Apply a function across the rows to sort the values in var1 and var2
  sorted_names <- t(apply(sub2_all[, c("var1", "var2")], 1, function(x) sort(x)))
  # Replace the original columns with the sorted data
  sub2_all$var1 <- sorted_names[, 1]
  sub2_all$var2 <- sorted_names[, 2]

  sub3 <- dplyr::filter(results, n_cov == 3) %>%
    dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0) + ifelse(var3 == covariates[i],1,0)) %>%
    dplyr::filter(keep == 1) %>% dplyr::select(-keep)
  sub3$rmse_diff <- 0
  sub3$AIC_diff <- 0
  for(j in 1:nrow(sub3)) {
    vars <- sub3[j,c("var1", "var2","var3")]
    vars <- sort(as.character(vars[which(vars != covariates[i])]))
    # find the same in sub2
    indx <- which(paste(sub2_all$var1, sub2_all$var2) == paste(vars, collapse=" "))
    sub3$rmse_diff[j] <- sub3$RMSE[j] / sub2_all$RMSE[indx]
    sub3$AIC_diff[j] <- sub2_all$AIC[indx] - sub3$AIC[j]
  }

  # Fill in summary stats
  marginals$rmse_12[i] <- mean(sub2$rmse_diff)
  marginals$aic_12[i] <- mean(sub2$AIC_diff)
  marginals$rmse_23[i] <- mean(sub3$rmse_diff)
  marginals$aic_23[i] <- mean(sub3$AIC_diff)
}

# Calculate avergages of averages
marginals$total_rmse <- apply(marginals[,c("rmse_01","rmse_12", "rmse_23")], 1, mean)
marginals$total_aic <- apply(marginals[,c("aic_01","aic_12", "aic_23")], 1, mean)

# CREATE TABLE OF RMSE/AIC------------------------------------

gam_loo_table <- dplyr::arrange(marginals, total_rmse)%>%
  dplyr::select(cov, total_rmse, total_aic)%>%
  mutate(cv="LOO",model="GAM")

# Fit a model that includes multiple variables from the top model
combos = c(gam_loo_table$cov[1], gam_loo_table$cov[2], gam_loo_table$cov[3], gam_loo_table$cov[4], gam_loo_table$cov[5])
smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", smooth_terms)
gam_model <- gam(as.formula(formula_str),
                 #family = binomial(),
                 #weights = number_cwt_estimated,
                 data = dat)

# Make a plot of the marginal effects from this best model
the_lm2 <-  lm(Y_rec ~  LSTlarv + ONIpjuv + CutiSTIpjuv + MLDpjuv + Tcop, data=dat)
vif(the_lm2)

partial_effects <- gratia::draw(gam_model, transform = "response", cex = 2)
summary(gam_model)


#### GAMS with leave future out cross valudation ####

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
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Extract variable names
 padded_vars <- c(var_names, rep(NA, 3 - length(var_names)))
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
     var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3]
    
    
  ))
    predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[kstart:n_year],
     var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3]
    
    
  ))
 
   print(i)
}
  
results_arr <- arrange(results,RMSE )
# View the results data frame
print(results_arr)



# Create a baseline model with only the intercept

    predictions <- numeric(nrow(dat))
    n_year <- length(unique(dat$year))

    # Loop over each observation
    for (j in 1:n_year) {
      # Fit model on n-1 observations
      gam_model <- gam(Y_rec ~ 1,
                       data = dat[which(dat$year != unique(dat$year)[j]), ])
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
      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    }

baseline_rmse  <- sqrt(mean((dat$Y_rec - (predictions))^2, na.rm=T))

# filter for frowns only

results_arr_RMSE <- arrange(results,RMSE )
# View the results data frame
print(results_arr_RMSE)
results_arr_AIC <- arrange(results,AIC )
# View the results data frame
print(results_arr_AIC)

results_gam_lfo5<-results%>%
  mutate(RMSE_improv=RMSE/baseline_rmse, cv="LFO_5", model="GAM")

# calculate the baseline RMSE

# we can calculate the marginal improvement for each covariate
results$n_cov <- ifelse(!is.na(results$var1), 1, 0) + ifelse(!is.na(results$var2), 1, 0) +
  ifelse(!is.na(results$var3), 1, 0)
marginals <- data.frame(cov = covariates, "rmse_01" = NA, "rmse_12" = NA, "rmse_23" = NA,
                        "aic_01" = NA, "aic_12" = NA, "aic_23" = NA)
for(i in 1:length(covariates)) {
  sub <- dplyr::filter(results, n_cov == 1,
                       var1 == covariates[i])
  marginals$rmse_01[i] <- sub$RMSE / baseline_rmse
  marginals$aic_01[i] <- AIC(gam_model) - sub$AIC

  # next look at all values of models that have 2 covariates and include this model
  sub1 <- dplyr::filter(results, n_cov == 1)
  sub2 <- dplyr::filter(results, n_cov == 2) %>%
    dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0)) %>%
    dplyr::filter(keep == 1) %>% dplyr::select(-keep)
  # loop over every variable in sub2, and find the simpler model in sub1 that just represents the single covariate
  sub2$rmse_diff <- 0
  sub2$AIC_diff <- 0
  for(j in 1:nrow(sub2)) {
    vars <- sub2[j,c("var1", "var2")]
    vars <- vars[which(vars != covariates[i])]
    indx <- which(sub1$var1 == as.character(vars))
    sub2$rmse_diff[j] <- sub2$RMSE[j] / sub1$RMSE[indx]
    sub2$AIC_diff[j] <- sub1$AIC[indx] - sub2$AIC[j]
  }

  # Finally compare models with 3 covariates to models with 2 covariates
  sub2_all <- dplyr::filter(results, n_cov == 2)
  # Apply a function across the rows to sort the values in var1 and var2
  sorted_names <- t(apply(sub2_all[, c("var1", "var2")], 1, function(x) sort(x)))
  # Replace the original columns with the sorted data
  sub2_all$var1 <- sorted_names[, 1]
  sub2_all$var2 <- sorted_names[, 2]

  sub3 <- dplyr::filter(results, n_cov == 3) %>%
    dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0) + ifelse(var3 == covariates[i],1,0)) %>%
    dplyr::filter(keep == 1) %>% dplyr::select(-keep)
  sub3$rmse_diff <- 0
  sub3$AIC_diff <- 0
  for(j in 1:nrow(sub3)) {
    vars <- sub3[j,c("var1", "var2","var3")]
    vars <- sort(as.character(vars[which(vars != covariates[i])]))
    # find the same in sub2
    indx <- which(paste(sub2_all$var1, sub2_all$var2) == paste(vars, collapse=" "))
    sub3$rmse_diff[j] <- sub3$RMSE[j] / sub2_all$RMSE[indx]
    sub3$AIC_diff[j] <- sub2_all$AIC[indx] - sub3$AIC[j]
  }

  # Fill in summary stats
  marginals$rmse_12[i] <- mean(sub2$rmse_diff)
  marginals$aic_12[i] <- mean(sub2$AIC_diff)
  marginals$rmse_23[i] <- mean(sub3$rmse_diff)
  marginals$aic_23[i] <- mean(sub3$AIC_diff)
}

# Calculate avergages of averages
marginals$total_rmse <- apply(marginals[,c("rmse_01","rmse_12", "rmse_23")], 1, mean)
marginals$total_aic <- apply(marginals[,c("aic_01","aic_12", "aic_23")], 1, mean)

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
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Extract variable names
 padded_vars <- c(var_names, rep(NA, 3 - length(var_names)))
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
     var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3]
    
    
  ))
    predicted <- rbind(predicted, data.frame(
    ModelID = i,
    pred=predictions[kstart:n_year],
     var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3]
    
    
  ))
 
   print(i)
}
  
results_arr <- arrange(results,RMSE )
# View the results data frame
print(results_arr)



# Create a baseline model with only the intercept

    predictions <- numeric(nrow(dat))
    n_year <- length(unique(dat$year))

    # Loop over each observation
    for (j in 1:n_year) {
      # Fit model on n-1 observations
      gam_model <- gam(Y_rec ~ 1,
                       data = dat[which(dat$year != unique(dat$year)[j]), ])
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
      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    }

baseline_rmse  <- sqrt(mean((dat$Y_rec - (predictions))^2, na.rm=T))

# filter for frowns only

results_arr_RMSE <- arrange(results,RMSE )
# View the results data frame
print(results_arr_RMSE)
results_arr_AIC <- arrange(results,AIC )
# View the results data frame
print(results_arr_AIC)

results_gam_lfo10<-results%>%
  mutate(RMSE_improv=RMSE/baseline_rmse, cv="LFO_10", model="GAM")



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

# Make a plot of the marginal effects from this best model
the_lm3 <-  lm(Y_rec ~ LSTlarv + PDOlarv + Tcop, data=dat)
vif(the_lm3)

partial_effects <- gratia::draw(gam_model, transform = "response", cex = 2)
summary(gam_model)

fullresults <- bind_rows(results_gam_lfo10,results_gam_lfo5)%>%
  bind_rows(results_gam_loo)%>%
  mutate(var_all=gsub('NA', '', paste(var1,var2,var3)))
  

write.csv(rbind(gam_loo_table,table_5,table_10),"MarginalImprovementGAM.csv")
write.csv(fullresults,"MultivariateGAMresults.csv")









