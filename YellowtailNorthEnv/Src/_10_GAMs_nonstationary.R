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
# One decision: what's the maximum number of covariates that any model can have?
combinations <- lapply(1:3, function(i) {
  combn(covariates, i, simplify = FALSE)
})
combinations <- unlist(combinations, recursive = FALSE)
# Sort each combination and remove duplicates
combinations <- unique(lapply(combinations, function(x) sort(x)))

cross_validation <-TRUE
models <- list()
results <- data.frame()
jstart<-1
n_year <-length(unique(dat$year))


  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(year,by=", combinations[[i]], ", k = 6)", collapse = " + ")
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



# filter for frowns only

results_arr <- arrange(results,RMSE )
# View the results data frame
print(results_arr)



# Create a baseline model with only the intercept

    predictions <- numeric(nrow(dat))
    n_year <- length(unique(dat$year))

    # Loop over each observation
    for (j in 1:n_year) {
      train_index <- setdiff(1:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index

      # Fit model on n-1 observations
      gam_model <- gam(Y_rec ~ s(year, k=6),
                       data = dat[which(dat$year != unique(dat$year)[j]), ])

      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    }

baseline_rmse  <- sqrt(mean((dat$Y_rec - (predictions))^2, na.rm=T))

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

table <- dplyr::arrange(marginals, total_rmse) %>%
  dplyr::select(cov, total_rmse, total_aic)

# Fit a model that includes multiple variables from the top model
combos = c(table$cov[1], table$cov[2], table$cov[3], table$cov[4], table$cov[5])
smooth_terms <- paste("s(year, by=", combos, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", smooth_terms)
gam_model <- gam(as.formula(formula_str),
                 #family = binomial(),
                 #weights = number_cwt_estimated,
                 data = dat)

# Make a plot of the marginal effects from this best model
the_lm2 <-  lm(Y_rec ~  LSTlarv + ONIpjuv + CutiSTIpjuv + MLDpjuv + Tcop, data=dat)
vif(the_lm2)

partial_effects <- gratia::draw(gam_model, transform = "response", cex = 2)
ggsave(partial_effects, filename = paste0("plots/coho/coho_",run,"_gam_partial_effects.png"), height = 7, width = 10)
summary(gam_model)

