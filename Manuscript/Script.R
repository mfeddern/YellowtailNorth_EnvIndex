# Load required libraries
# Load all necessary libraries at the beginning of the script for clarity.
# Note: Several packages are listed multiple times, consolidate for efficiency.
library(dplyr)
library(tidyr)
library(ggpubr)
library(corrplot)
library(mgcv)
library(DHARMa) # For diagnostics of GL(A)M(M)s
library(mgcViz) # For visualizing GAMs
library(gridExtra) # For arranging multiple plots
library(ROCR) # For evaluating classifier performance
library(recdata) # Assuming this is a custom package, ensure it's available
library(predRecruit) # Assuming this is a custom package, ensure it's available
library(tidyverse) # Includes dplyr, tidyr, ggplot2, etc. (supersedes some individual loads)
library(car) # For companion to applied regression

#### --- 1. Data Loading and Preparation --- ####

# Read in the environmental data
# Use `read_csv` from `readr` package (part of tidyverse) for faster and more efficient data loading.
env_data_raw <- read_csv("data-yellowtail/2024Env-annual-yellowtail-standardized.csv")

# Rename columns for clarity and consistency
env_data_renamed <- env_data_raw %>% 
  rename(
    ONIlarv=oni_larv,  ONIpjuv=oni_pjuv,  ONIpre=oni_pre,
    PDOlarv=pdo_larv, PDOpjuv=pdo_pjuv,
    LUSI=lusi_annual,
    BeutiTUMI=BeutiTUMIpjuv, BeutiSTI=BeutiSTIpjuv,
    CutiSTI= CutiSTIpjuv, CutiTUMI = CutiTUMIpjuv, 
    BakunSTI=  bakun_sti,
    HCIlarv= hci1_larv, HCIpjuv=hci1_pjuv,
    HCI2larv=hci2_larv, HCI2pjuv=hci2_pjuv,
    NCOPpjuv = ZOOpjuvN, NCOPben = ZOObenN,
    SCOPpjuv = ZOOpjuvS, SCOPben = ZOObenS
  )

# Define columns to omit
columns_to_omit <- c('X',"Year","Year.1","LUSI", "BakunSTI", "HCI2larv", "HCI2pjuv", "HCIlarv", "HCIpjuv")

# Select relevant environmental variables
envir_data_cleaned <- env_data_renamed %>%
  select(!all_of(columns_to_omit))

# Read in recruitment deviation data
recruitment_deviations <- read_csv("data-yellowtail/RecruitmentDeviations2025draft.csv") %>%
  filter(Datatreatment == "2025 No_SMURF")

# Join environmental data with recruitment deviations based on year
full_dataset <- envir_data_cleaned %>%
  left_join(recruitment_deviations, by = "year")

# Define the years for analysis
analysis_years <- 1998:2021

# Filter the data for the specified years
analyzed_data <- full_dataset %>% 
  filter(year %in% analysis_years)

#### --- 2. Evaluating Correlations Between Covariates --- ####

# Define a correlation threshold to filter out highly correlated variables
correlation_threshold <- 0.3

# Prepare environmental data for correlation analysis
# Ensure complete cases for correlation calculation and remove the 'year' column
envir_data_for_cor <- envir_data_cleaned %>% 
  filter(complete.cases(.)) %>%
  select(-year)

# Generate a correlation matrix
correlation_matrix <- cor(envir_data_for_cor)

# Convert the correlation matrix to a long format for easier filtering
long_correlation_data <- correlation_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("xvar") %>%
  pivot_longer(!xvar, names_to = "yvar", values_to = "correlation")

# Identify highly correlated variable pairs
correlated_pairs <- long_correlation_data %>%
  filter(abs(correlation) > correlation_threshold, xvar != yvar) %>%
  select(-correlation)

# Create a list of combinations to potentially omit
combinations_to_omit <- lapply(1:nrow(correlated_pairs), function(i) {
  as.character(correlated_pairs[i, ])
})

#### --- 3. Generate Null Models --- ####

# Fit a null GAM model with only an intercept
null_gam_model <- gam(Y_rec ~ 1, data = analyzed_data)

# Calculate null predictions (intercept only)
analyzed_data <- analyzed_data %>%
  mutate(sr_null = 0, mean_null = predict(null_gam_model))

# Calculate Root Mean Squared Error (RMSE) for the null models (full dataset)
rmse_sr_full <- sqrt(mean((analyzed_data$sr_null - analyzed_data$Y_rec)^2, na.rm = TRUE))
rmse_mean_full <- sqrt(mean((analyzed_data$mean_null - analyzed_data$Y_rec)^2, na.rm = TRUE))

# Generate null model predictions for the LOOCV timeframe
null_lfo_data <- analyzed_data %>%
  mutate(sr_null = 0, mean_null = predict(null_gam_model)) %>%
  filter(year >= 2014 & year < 2019)

# Calculate RMSE for the null models (LOOCV timeframe)
rmse_sr_lfo <- sqrt(mean((null_lfo_data$sr_null - null_lfo_data$Y_rec)^2, na.rm = TRUE))
rmse_mean_lfo <- sqrt(mean((null_lfo_data$mean_null - null_lfo_data$Y_rec)^2, na.rm = TRUE))


#### --- 4. Model Specification and Filtering Combinations --- ####

# Get names of available covariates from the cleaned environmental data
available_covariates <- names(envir_data_for_cor)

# Define the maximum number of covariates for a model
max_covariates <- 4 

# Generate all possible combinations of covariates
# Use `combn` to create combinations efficiently
all_combinations <- unlist(lapply(1:max_covariates, function(i) {
  combn(available_covariates, i, simplify = FALSE)
}), recursive = FALSE)

# Ensure unique combinations and sort for consistency
all_combinations <- unique(lapply(all_combinations, sort))

# Function to check if a combination is a partial match to a combination to omit
is_partial_match <- function(combination, omit_list) {
  any(sapply(omit_list, function(omit_pair) all(omit_pair %in% combination)))
}

# Remove combinations containing highly correlated variables
filtered_combinations <- all_combinations[!sapply(all_combinations, is_partial_match, omit_list = combinations_to_omit)]

#### --- 5. Leave-One-Out Cross-Validation (LOOCV) and Model Fitting --- ####

# Initialize lists to store model results
model_list <- list()
model_results <- data.frame()
predicted_values <- data.frame()

# Extract unique years for LOOCV
unique_years <- unique(analyzed_data$year)

# Loop through each combination of covariates to fit and evaluate GAM models
for (i in seq_along(filtered_combinations)) {
  
  # Construct the smooth terms for the GAM formula
  smooth_terms <- paste0("s(", filtered_combinations[[i]], ", k = 3)", collapse = " + ")
  formula_string <- paste("Y_rec ~", smooth_terms)
  
  # Initialize a numeric vector to store LOOCV predictions
  loo_predictions <- numeric(nrow(analyzed_data))
  
  # Perform Leave-One-Out Cross-Validation (LOOCV)
  for (j in seq_along(unique_years)) {
    # Define training and testing data for each fold
    training_data <- analyzed_data %>% filter(year != unique_years[j])
    testing_data <- analyzed_data %>% filter(year == unique_years[j])
    
    # Fit the GAM model on the training data
    # Consider using REML or ML for smoothness selection in `gam` for better performance and stability
    gam_loo_model <- gam(as.formula(formula_string), data = training_data, method = "REML") 
    
    # Predict on the excluded observation (test data)
    loo_predictions[which(analyzed_data$year == unique_years[j])] <- predict(gam_loo_model, newdata = testing_data)
  }
  
  # Refit the full model for R-squared and deviance explained
  full_gam_model <- gam(as.formula(formula_string), data = analyzed_data, method = "REML")
  
  # Get predictions from the full model
  full_model_predictions <- predict(full_gam_model)
  
  # Calculate RMSE for LOOCV and the full model
  rmse_loo <- sqrt(mean((analyzed_data$Y_rec - loo_predictions)^2, na.rm = TRUE))
  rmse_full_model <- sqrt(mean((analyzed_data$Y_rec - full_model_predictions)^2, na.rm = TRUE))
  
  # Extract R-squared and explained deviance
  r_squared <- summary(full_gam_model)$r.sq
  deviance_explained <- summary(full_gam_model)$dev.expl
  
  # Extract covariate names and pad with NA for consistent data frame structure
  variable_names <- filtered_combinations[[i]]
  padded_variables <- c(variable_names, rep(NA, max_covariates - length(variable_names)))
  
  # Store results
  model_list[[i]] <- full_gam_model
  model_results <- bind_rows(model_results, data.frame(
    ModelID = i,
    AIC = AIC(full_gam_model),
    RMSE_loo = round(rmse_loo, 3),
    RMSE = round(rmse_full_model, 3),
    rsq_full = round(r_squared, 2),
    dev.ex = round(deviance_explained, 4),
    rmse_imp = (rmse_sr_full - rmse_full_model) / rmse_sr_full,
    rmse_ratio = rmse_full_model / rmse_sr_full,
    var1 = padded_variables[1],
    var2 = padded_variables[2],
    var3 = padded_variables[3],
    var4 = padded_variables[4]
  ))
  
  # Store one-step-ahead predictions
  predicted_values <- bind_rows(predicted_values, data.frame(
    ModelID = i,
    pred = loo_predictions,
    year = analyzed_data$year,
    var1 = padded_variables[1],
    var2 = padded_variables[2],
    var3 = padded_variables[3],
    var4 = padded_variables[4]
  ))
  
  # Print progress
  print(paste("Model", i, "of", length(filtered_combinations), "completed."))
}

#### --- 6. Results Analysis and Output --- ####

# Arrange results by AIC and add rank and delta AIC
model_results <- model_results %>%
  arrange(AIC) %>%
  mutate(rankmod = seq_along(ModelID)) %>%
  mutate(deltaAIC = AIC - min(AIC))

# Arrange results by RMSE (LOOCV) and AIC
results_arr_rmse <- arrange(model_results, RMSE_loo)
results_arr_aic <- arrange(model_results, AIC)

# Filter for top models based on Delta AIC
summary_top_models <- model_results %>% filter(deltaAIC < 2)

# Write top models to a CSV file
write_csv(summary_top_models, "Figures/No_SMURF_Sensitivity/TopModels.csv")

#### --- 7. Calculating Relative Variable Importance --- ####

# Calculate the number of covariates for each model
model_results$n_cov <- rowSums(!is.na(model_results[, c("var1", "var2", "var3", "var4")]))

# Initialize a data frame to store marginal importance values
marginal_importance <- data.frame(
  cov = available_covariates, 
  rmse_01 = NA, rmse_12 = NA, rmse_23 = NA,
  aic_01 = NA, aic_12 = NA, aic_23 = NA
)

# Filter model results to include only models with up to 3 covariates
filtered_results_for_marginal <- model_results %>% 
  filter(n_cov < 4) %>%
  select(-var4)

# Loop to calculate marginal contributions (RMSE and AIC)
for (i in seq_along(available_covariates)) {
  # Marginal improvement from 0 to 1 covariate
  sub_n1 <- filtered_results_for_marginal %>% 
    filter(n_cov == 1, var1 == available_covariates[i])
  marginal_importance$rmse_01[i] <- (rmse_sr_full - sub_n1$RMSE) / rmse_sr_full
  marginal_importance$aic_01[i] <- AIC(null_gam_model) - sub_n1$AIC
  # This section requires careful implementation to calculate the 'average shared variance' or similar metrics for variable importance in GAMs.
}

#### --- 8. Extracting and Analyzing the Highest Ranked LOO-CV Model --- ####

# Determine which ranked model to use for further analysis (e.g., rank 1 for the best)
selected_rank_num <- 1
selected_model_info <- results_arr_rmse[selected_rank_num, ]
combos= c(results_arr_rmse$var1[selected_rank_num], results_arr_rmse$var2[selected_rank_num], results_arr_rmse$var3[selected_rank_num]) #this needs to be mod depending on which model and selection criteria you want

# Extract the covariates for the selected model.
# Adjust the indexing `var1:varX` based on the actual number of covariates in the selected model
# and filter out any `NA` values.
selected_covariates <- as.character(selected_model_info %>%
                                      select(var1, var2, var3, var4) %>% # Select up to var4, adjust if needed
                                      pivot_longer(everything(), names_to = "var_name", values_to = "covariate") %>%
                                      filter(!is.na(covariate)) %>%
                                      pull(covariate))

# Construct the smooth terms for the GAM formula using the selected covariates.
# a small number of knots (`k=3`).
smooth_terms_selected_model <- paste0("s(", selected_covariates, ", k = 3)", collapse = " + ")
formula_selected_model <- as.formula(paste("Y_rec ~", smooth_terms_selected_model))

# Filter the data to the full fit years and fit the GAM model.
gam_model_full_fit <- gam(formula_selected_model, data = analyzed_data) 

# Print a summary of the selected model.
summary(gam_model_full_fit)

#### --- 9. Predicting the Last 5 Years (Fixed Training Window) --- ####

# Define training and prediction periods.
training_years_fixed <- 1998:2016
prediction_years_fixed <- 2017:2021

# Filter data for the training period.
df_train_fixed <- full_dataset %>% filter(year %in% training_years_fixed)

# Fit the GAM model to the fixed training data.
gam_fixed_train <- gam(formula_selected_model, data = df_train_fixed)

# Predict for the fixed prediction period (last 5 years).
df_predict_fixed <- full_dataset %>% filter(year %in% prediction_years_fixed)
predictions_last5_fixed <- predict.gam(gam_fixed_train, newdata = df_predict_fixed)

# Store these predictions with their corresponding years.
df_predictions_last5_fixed <- data.frame(
  year = prediction_years_fixed,
  pred_last5_fixed = predictions_last5_fixed
)

#### --- 10. Predicting the Last 2 Years (Fixed Training Window) --- ####

# Define a training period ending two years before the end of the full dataset
training_years_last2 <- 1998:2019
prediction_years_last2 <- 2020:2021

# Filter data for the adjusted training period.
df_train_last2 <- full_dataset %>% filter(year %in% training_years_last2)

# Fit the GAM model to this training data.
gam_fixed_train_last2 <- gam(formula_selected_model, data = df_train_last2)

# Predict for the last two years.
df_predict_last2 <- full_dataset %>% filter(year %in% prediction_years_last2)
predictions_last2_fixed <- predict.gam(gam_fixed_train_last2, newdata = df_predict_last2)

# Store these predictions.
df_predictions_last2_fixed <- data.frame(
  year = prediction_years_last2,
  pred_last2_fixed = predictions_last2_fixed
)

#### --- 11. Full Dataset Predictions and Combining Data --- ####

# Generate predictions and their standard errors for the entire analysis period (1998-2021)
# using the model fitted to the full data (gam_model_full_fit).
gam_predictions_full_data <- predict(gam_model_full_fit, se.fit = TRUE)

# Create a data frame for the full dataset predictions and standard errors.
df_gam_train_full <- data.frame(
  year = analyzed_data$year, # Assuming df_analysis has years 1998-2021
  fit = gam_predictions_full_data$fit,
  se.fit = gam_predictions_full_data$se.fit
)

# Combine all prediction data into a single data frame.
# This assumes `predicted` (from the LOOCV loop in the main script) is available and contains `ModelID` and `pred` columns.
# It also assumes `predicted` already contains LOOCV predictions for all years in `df_analysis`.
df_combined_predictions <- analyzed_data%>% # Start with the main analysis data (contains Y_rec)
  left_join(df_gam_train_full %>% select(year, fit, se.fit), by = "year") %>% # Add full model predictions
  left_join(predicted_values %>% filter(ModelID == selected_model_info$ModelID) %>% # Add LOOCV predictions for the selected model
              select(year, pred), by = "year") %>%
  rename(pred_loo_selected_model = pred) %>% # Rename the LOOCV predictions for clarity
  left_join(df_predictions_last5_fixed, by = "year") %>% # Add fixed window predictions (last 5 years)
  left_join(df_predictions_last2_fixed, by = "year") # Add fixed window predictions (last 2 years)

#### --- 12. One-Step-Ahead Predictions for the Last 5 Years (Rolling Window) --- ####

# Initialize a data frame to store the one-step-ahead predictions.
df_pred_one_step_ahead <- data.frame()

# Define the starting year for the one-step-ahead prediction sequence.
# Assuming the "last 5 years" for this one-step-ahead process are 2017-2021,
# and the training starts before 2017. The loop simulates rolling forecasts.
one_step_ahead_years <- 2017:2021

for (current_year_pred in one_step_ahead_years) {
  # Define the training data as all years before the current prediction year.
  training_data_roll <- analyzed_data %>%
    filter(year < current_year_pred)
  
  # Fit the GAM model to the rolling training data.
  gam_roll_model <- gam(formula_selected_model, data = training_data_roll, method = "REML")
  
  # Define the data for which to make the prediction (the single next year).
  newdata_roll <-analyzed_data %>%
    filter(year == current_year_pred)
  
  # Make the one-step-ahead prediction.
  prediction_roll <- predict(gam_roll_model, newdata = newdata_roll)
  
  # Store the prediction and the year.
  df_pred_one_step_ahead <- bind_rows(df_pred_one_step_ahead,
                                      data.frame(year = current_year_pred,
                                                 pred_one_step = prediction_roll))
}

# Combine the one-step-ahead predictions with the main combined prediction data frame.
df_combined_predictions <- df_combined_predictions %>%
  left_join(df_pred_one_step_ahead, by = "year")


#### --- 13. Model Diagnostics --- ####
###### --- 1. Cook's Distance --- ####

analyzed_data$cooks.distance <- cooks.distance(gam_model_full_fit)
analyzed_data$leverage <-influence.gam(gam_model_full_fit)
analyzed_data$residuals.standard <- scale(residuals.gam(gam_model_full_fit))
n <- length(analyzed_data$residuals.standard )  # Number of data points
influence  <- analyzed_data[(analyzed_data$cooks.distance> 4 / n |analyzed_data$cooks.distance> 1|analyzed_data$cooks.distance> 0.5),]

###### --- 2. Leave-one-out resampling on Deviance Explained --- ####

loo_dexex <- data.frame()
jstart<-1
predictions <- numeric(nrow(analyzed_data))
dev.expjack <- data.frame()
n_year <- length(unique(analyzed_data$year))
for (j in jstart:n_year) {
  train_index <- setdiff(jstart:n_year, j)  # All indices except the j-th
  test_index <- j                 # The j-th index
  # Fit model on n-1 observations
  gam_model <- gam(as.formula(formula_selected_model),
                   # weights = number_cwt_estimated,
                   data = analyzed_data[which(analyzed_data$year != unique(analyzed_data$year)[j]), ])
  # Predict the excluded observation
  predictions[which(analyzed_data$year == unique(analyzed_data$year)[j])] <- predict(gam_model, newdata = analyzed_data[which(analyzed_data$year == unique(analyzed_data$year)[j]), ])
  dev.exp<-summary(gam_model)$dev.exp
  dev.expjack<-rbind(dev.expjack,dev.exp)
}

gam.jack<-cbind(analyzed_data,jack=predictions[1:24],dev.exp=dev.expjack[1:24,1])

###### --- 3. Dispersion Testing --- ####

simulationOutput <- simulateResiduals(fittedModel = gam_model_full_fit)
testDispersion(simulationOutput)
simulationOutput$scaledResiduals

plotResiduals(simulationOutput$scaledResiduals, analyzed_data$CutiSTI)
plotResiduals(simulationOutput$scaledResiduals, analyzed_data$ONIpjuv)
plotResiduals(simulationOutput$scaledResiduals, analyzed_data$LSTpjuv)

#### --- 14. Generating Main Text Figures --- ####
##### --- 1. Generating Figure 1 --- ####

yrmin <- 1996
yrmax<- 2025
recruitmentdevs<- data.frame(read.csv("data-yellowtail/RecruitmentDeviations2025draft.csv"))%>%
  filter(year>=yrmin & year <=yrmax)
Bio_param<-read.csv("data-yellowtail/BiologicalParameters2025.csv")#%>%
# filter(Year>=yrmin & Year <=yrmax)
spawn_rec<-read.csv("data-yellowtail/spwn_rec_curve_non_SMURF.csv")%>%
  filter(Yr>=yrmin & Yr <=yrmax)
curve<- read.csv("data-yellowtail/spwn_rec_curve_non_SMURF_to0.csv")

SpawningBiomassTS <- ggplot(data=Bio_param, aes(y=Spawning.output..trillions.of.eggs., 
                                                x=Year))+
  geom_line()+ theme_classic() +
  labs(x="Year", y="Spawning Output \n (trillions of eggs)")+
  xlim(c(1995,2020))+
  #geom_hline(yintercept=0,lty=2)+
  ylim(c(0,11))+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
SpawningBiomassTS

Age0SpawningBiomass <- ggplot(data=curve, aes(y=Recruitment, x=SSB))+
  geom_line()+
  geom_point(data=Bio_param%>%filter(Year>1945), aes(y=Age.0.Recruits..1.000s., x=Spawning.output..trillions.of.eggs., colour=Year))+ 
  theme_classic() +
  labs(x="Spawning Output (trillions of eggs)", y="Age-0 (x1000)")+
  #geom_hline(yintercept=0,lty=2)+
  ylim(c(0,70000))+xlim(c(0,15))+
  scale_color_gradientn(colours =c("#08519C","#0f85a0","#dd4124"))+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
Age0SpawningBiomass

Age0TS <- ggplot(data=spawn_rec, aes(y=pred_recr, x=Yr))+
  geom_line()+ theme_classic() +
  ylim(c(0,75000))+
  labs(x="Year", y="Age-0 (x1000)")+
  #geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
Age0TS

Figure1<-ggarrange(SpawningBiomassTS,Age0TS, Age0SpawningBiomass, nrow=3,
                   labels = c('A', 'B', 'C'))
png("Figures/Figure1.png",width=4,height=8,units="in",res=1200)
Figure1
dev.off()

pdf(file = "Figures/Figure1.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 8)
Figure1
dev.off()

##### --- 2. Generating Figure 3 --- ####
tsdata<-pivot_longer(env_data_renamed,cols=-c(year, X, Year, Year.1),names_to = 'var',values_to = "val")%>%
  mutate(type = ifelse(var=='CHLpjuv'|var=='PPpjuv'|var=="NCOPben"|var=="NCOPpjuv"|var=="SCOPben"|var=="SCOPpjuv", "Foodweb", "Oceanographic"))

ts_foodweb <- ggplot(tsdata%>%filter(type=="Foodweb"),aes(x=year,y=val))+
  geom_line(lwd=0.9, col="#7cae00" )+
  # geom_point()+
  facet_wrap(~var, ncol=6)+
  ylab("")+
  xlab("")+
  ggtitle("Foodweb Conditions")+
  geom_hline(yintercept = 0,lty=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

ts_ocean <-ggplot(tsdata%>%filter(type=="Oceanographic"),aes(x=year,y=val))+
  geom_line(lwd=1, col ="#0f85a0" )+
  #geom_point()+
  facet_wrap(~var, ncol=6)+
  ylab("")+
  xlab("")+
  ggtitle("Oceanographic Conditions")+
  geom_hline(yintercept = 0,lty=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))
arranged_ts<-ggarrange(ts_foodweb,ts_ocean, ncol = 1, heights=c(1.5,4.25))

full_ts<-annotate_figure(arranged_ts,
                         bottom = text_grob("Year"),
                         left = text_grob("Standardized Index", rot=90)
)

png("Figures/Figure3.png",width=10,height=8,units="in",res=1200)
full_ts
dev.off()

pdf(file ="Figures/Figure3.pdf",width=10,height=8)
full_ts
dev.off()

##### --- 3. Generating Figure 4 --- ####
dfx = envir_data_for_cor #%>%select(-c(HCIpjuv, HCIlarv, HCI2larv, HCI2pjuv)) 
dfx = na.omit(dfx) # na's mess up the cor function
cor_xy = cor(dfx)
uncorr_sorted  <- long_correlation_data %>%
  filter(abs(correlation) < correlation_threshold, xvar != yvar) %>%
  select(-correlation) %>%
  mutate(key=paste0(pmin(xvar, yvar), pmax(xvar,yvar), sep=""))%>%
  distinct(key, .keep_all=TRUE)%>%
  select(-key)
#unique(lapply(uncorr, function(x) sort(x)))
r<-data.frame()
#r <- rbind(c("DDpre", "MLDpjuv","DDpre","MLDpjuv"))
for(i in 1:length(uncorr_sorted$xvar)){
  temp<-c(uncorr_sorted$xvar[i], uncorr_sorted$yvar[i], uncorr_sorted$xvar[i], uncorr_sorted$yvar[i])
  r<-rbind(r, temp) 
}

graphics.off()
png( paste0("Figures/Figure4.png"),
     units = 'in', res=300, width = 8, height=8)
corrplot::corrplot(cor_xy, method='ellipse',  type='lower', tl.cex = 0.75)%>% 
  corrRect(namesMat = r)
dev.off()

pdf( paste0("Figures/Figure4.pdf"),
     width = 8, height=8)
# plot command
corrplot::corrplot(cor_xy, method='ellipse',  type='lower', tl.cex = 0.75)%>% 
  corrRect(namesMat = r)
dev.off()

##### --- 4. Generating Figure 5 --- ####
sizept<-2.5

rec_panela<-ggplot(df_combined_predictions%>%filter(year  %in% 1998:2021), aes(x=year)) +
  geom_line(data=df_combined_predictions,aes(y=fit)) +
  geom_ribbon(data=df_combined_predictions,aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  geom_point(aes(y=Y_rec),size=sizept) +
  geom_errorbar(aes(ymin=Y_rec-sd,ymax=Y_rec+sd),width = 0)+
  theme_classic() +
  labs(x="", y="")+
  ylim(c(-1.5,1.5))+
  xlim(c(1998,2022))

rec_panelc<-ggplot(df_combined_predictions%>%filter(year  %in% 1998:2022), aes(x=year)) +
  geom_line(data=df_combined_predictions,aes(y=fit)) +
  geom_ribbon(data=df_combined_predictions,aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  geom_point(data=df_combined_predictions,aes(y=pred_last5_fixed),pch=21,fill='#dd4124',size=sizept) +
  theme_classic() +
  labs(x="", y="")+
  # ylim(c(-1.5,1.5))+
  xlim(c(1998,2022)) 


rec_paneld<-ggplot(df_combined_predictions%>%filter(year  %in% 1998:2023), aes(x=year)) +
  geom_line(data=df_combined_predictions,aes(y=fit)) +
  geom_ribbon(data=df_combined_predictions,aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  geom_point(data=df_combined_predictions,aes(y=pred_one_step), pch=21,fill="#edd746",size=sizept) +
  theme_classic() +
  labs(x="", y="")+
  ylim(c(-1.5,1.5))+
  xlim(c(1998,2023)) 

rec_panelb<-ggplot(df_combined_predictions%>%filter(year  %in% 1998:2023), aes(x=year)) +
  geom_line(data=df_combined_predictions,aes(y=fit)) +
  geom_ribbon(data=df_combined_predictions,aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  geom_point(data=df_combined_predictions,aes(y=pred_loo_selected_model), pch=21,fill='#0f85a0',size=sizept) +
  theme_classic() +
  labs(x="", y="")+
  ylim(c(-1.5,1.5))+
  xlim(c(1998,2023)) 

multi<-ggarrange(rec_panela,rec_panelb,rec_panelc,rec_paneld,labels = c("A.", "B.", "C.","D."))

png("Figures/Figure5.png",width=8,height=7,units="in",res=1200)
annotate_figure(multi,left="ln(recruitment devaitions)", bottom="Year")
dev.off()

pdf("Figures/Figure5.pdf",width=8,height=7)
annotate_figure(multi,left="ln(recruitment devaitions)", bottom="Year")
dev.off()

##### --- 5. Generating Figure 6 --- ####

smooth_data <- smooth_estimates(gam_model_full_fit)%>%  # Or specify the smooth term you want to plot
  add_constant(model_constant(gam_model_full_fit)) %>% # Add the intercept
  add_confint()%>% # Add the credible interval
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")%>%
  select(-c(.by))
observations<-analyzed_data%>%
  select(year,Y_rec,combos)%>%
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")
smooth_data <- na.omit(smooth_data )

partial_effects<-ggplot(smooth_data, aes(x = Value, y = .estimate)) +  # Setup the plot with the fitted values
  facet_wrap(~Smooth,ncol=1)+
  geom_line() + # Add estimated smooth
  xlim(c(-2,3))+
  geom_ribbon(aes(ymax = .upper_ci, ymin = .lower_ci), fill = "black", alpha = 0.2) + # Add credible interval
  geom_point(data = observations, aes(x = Value, y = Y_rec, col=year)) + # Add your data points
  labs(x = "Standardized Ecological Conditions", y = "Partial Residual")+ # Add labels
  geom_text(data = observations, aes(x = Value,col=year, y = Y_rec,label=year),hjust=0,nudge_x = 0.1)+
  scale_color_gradient(low="#0f85a0", high="#dd4124" )+
  theme_classic()
partial_effects 

png("Figures/Figure6.png",width=6,height=8,units="in",res=1200)
partial_effects 
dev.off()

pdf("Figures/Figure6.pdf",width=6,height=8)
partial_effects 
dev.off()

##### --- 6. Generating Figure 7 --- ####

gam_loo_table <- dplyr::arrange(marginals, rmse_23)%>%
  dplyr::select(cov, rmse_23)%>%
  mutate(cv="LOO",model="GAM")%>%
  left_join(included)
gam_loo_table[is.na(gam_loo_table)] <-"No"
cols<- c('#dd4124',"#edd746",'#7cae00','#0f85a0')

marginal <- ggplot(gam_loo_table, aes(x = reorder(cov,rmse_23), y = rmse_23, fill=Included)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the axes to make a horizontal bar graph
  labs(y = "Marginal Improvement RMSE", x = "Predictor")+
  scale_fill_manual(values=c('grey',cols[3]))+
  theme_classic()
marginal 

# creating a .png for google docs
png("Figures/Figure7.png",width=5.5,height=5,units="in",res=1200)
marginal
dev.off()

#creating a .pdf for submission
pdf(file = "Figures/Figure7.pdf",   # The directory you want to save the file in
    width = 5.5, # The width of the plot in inches
    height = 5)
marginal
dev.off()


cooks<-ggplot(analyzed_data, aes(x = cooks.distance, y = residuals.standard)) +
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
png("Figures/No_SMURF_Sensitivity/CooksRMSE.png",width=5,height=5,units="in",res=1200)
cooks
dev.off()



#### --- 15. Generating SI Figures --- ####
###### --- S1. Figure S1 --- ####
env_data_unst <- read_csv("data-yellowtail/2024Env-annual-yellowtail_UNSTANDARDIZED.csv")%>% 
  rename(
    ONIlarv=oni_larv,  ONIpjuv=oni_pjuv,  ONIpre=oni_pre,
    PDOlarv=pdo_larv, PDOpjuv=pdo_pjuv,
    LUSI=lusi_annual,
    BeutiTUMI=BeutiTUMIpjuv, BeutiSTI=BeutiSTIpjuv,
    CutiSTI= CutiSTIpjuv, CutiTUMI = CutiTUMIpjuv, 
    BakunSTI=  bakun_sti,
    HCIlarv= hci1_larv, HCIpjuv=hci1_pjuv,
    HCI2larv=hci2_larv, HCI2pjuv=hci2_pjuv,
    NCOPpjuv = ZOOpjuvN, NCOPben = ZOObenN,
    SCOPpjuv = ZOOpjuvS, SCOPben = ZOObenS
  )

tsdata_unst<-pivot_longer(env_data_unst,cols=-c(year, X, Year, Year.1),names_to = 'var',values_to = "val")%>%
  mutate(type = ifelse(var=='CHLpjuv'|var=='PPpjuv'|var=="NCOPben"|var=="NCOPpjuv"|var=="SCOPben"|var=="SCOPpjuv", "Foodweb", "Oceanographic"))

ts_foodweb <- ggplot(tsdata%>%filter(type=="Foodweb"),aes(x=year,y=val))+
  geom_line(lwd=0.9, col="#7cae00" )+
  # geom_point()+
  facet_wrap(~var, ncol=6, scales="free_y")+
  ylab("")+
  xlab("")+
  ggtitle("Foodweb Conditions")+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

ts_ocean <-ggplot(tsdata%>%filter(type=="Oceanographic"),aes(x=year,y=val))+
  geom_line(lwd=1, col ="#0f85a0" )+
  #geom_point()+
  facet_wrap(~var, ncol=6, scales = "free_y")+
  ylab("")+
  xlab("")+
  ggtitle("Oceanographic Conditions")+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))
arranged_ts<-ggarrange(ts_foodweb,ts_ocean, ncol = 1, heights=c(1.5,4.25))

full_ts<-annotate_figure(arranged_ts,
                         bottom = text_grob("Year"),
                         left = text_grob("Standardized Index", rot=90)
)

png("Figures/FigureS1.png",width=10,height=8,units="in",res=1200)
full_ts
dev.off()

pdf(file ="Figures/FigureS1.pdf",width=10,height=8)
full_ts
dev.off()

###### --- S2. Figure S2 --- ####

rec_foodweb <- ggplot(tsdata%>%
                        left_join(analyzed_data%>%select(year, Y_rec))%>%
                        filter(type=="Foodweb"),aes(x=val,y=Y_rec))+
  geom_smooth(method='gam',formula = y ~ s(x,  k=3), col='#7cae00')+
  geom_point()+
  facet_wrap(~var, ncol=6)+
  ylab("")+
  xlab("")+
  ggtitle("Foodweb Conditions")+
  
  geom_hline(yintercept = 0,lty=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))


rec_ocean <-ggplot(tsdata%>%
                     left_join(analyzed_data%>%select(year, Y_rec))%>%
                     filter(type=="Oceanographic"),aes(x=val,y=Y_rec))+
  geom_smooth(method='gam',formula = y ~ s(x,  k=3), col="#0f85a0")+
  
  geom_point()+
  facet_wrap(~var, ncol=6)+
  ylab("")+
  xlab("")+
  ggtitle("Oceanographic Conditions")+
  geom_hline(yintercept = 0,lty=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))
arranged_rec<-ggarrange(rec_foodweb,rec_ocean, ncol = 1, heights=c(1.5,4.25))

full_rec<-annotate_figure(arranged_rec,
                          bottom = text_grob("Standardized Index"),
                          left = text_grob("Recruitment Deviations", rot=90)
)
full_rec
png("Figures/FigureS2.png",width=10,height=8,units="in",res=1200)
full_rec
dev.off()

###### --- S3. Figure S3 --- ####
png("Figures/FigureS3.png",width=10,height=8,units="in",res=1200)
plot(simulationOutput)
dev.off()

###### --- S4. Figure S4 --- ####
png("Figures/FigureS4.png",width=10,height=8,units="in",res=1200)
testDispersion(simulationOutput)
dev.off()

###### --- S5. Figure S5 --- ####
cooks<-ggplot(analyzed_data, aes(x = cooks.distance, y = residuals.standard)) +
  geom_hline(yintercept = c(-3,3), linetype = 'dotted') +
  geom_point(alpha = 0.2) +
  ylab("Standardized Residuals")+
  xlab("Cook's Distance")+
  ggtitle("Influence")+
  geom_point(data=influence,aes(x = cooks.distance, y = residuals.standard), color = 'darkorange') +
  theme_minimal() +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5),  # Center plot title
    axis.title = element_text(size = 12),    # Axis title size
    axis.text = element_text(size = 10)      # Axis text size
  ) +
  #ylim(c(-4,4))+
  #geom_text(aes(label=ifelse(cooks.distance > 4 / n, year, "")), cex=3, hjust=1.1)
  geom_text(aes(label=year),cex=3, hjust=1.1)

png("Figures/FigureS5.png",width=5,height=5,units="in",res=1200)
cooks
dev.off()

###### --- S6. Figure S6 --- ####

dev.expjack<-ggplot(gam.jack, aes(dev.exp)) +
  geom_histogram(bins=40,fill="lightgrey", col="black")+
  xlim(c(0.4,0.8))+
  ylab("Frequency")+
  xlab("Deviance Explained")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.expjack
dev.expjack2<-ggplot(gam.jack, aes(x=year, y=dev.exp)) +
  geom_point(col="black")+
  # xlim(c(0.4,0.8))+
  ylim(c(0,1))+
  ylab("Deviance Explained")+
  geom_hline(yintercept=mean(na.omit(gam.jack$dev.exp)),lty=2)+#0.59
  xlab("Year (removed)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.expjack2
png("Figures/FigureS6.png",width=4,height=6,units="in",res=1200)
ggarrange(dev.expjack,dev.expjack2, ncol=1)
dev.off()

###### --- S7. Figure S7 --- ####
gam_model_short_fit <- gam(formula_selected_model, data = analyzed_data%>%filter(year<=2018)) 
smooth_data <- smooth_estimates(gam_model_short_fit)%>%  # Or specify the smooth term you want to plot
  add_constant(model_constant(gam_model_short_fit)) %>% # Add the intercept
  add_confint()%>% # Add the credible interval
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")%>%
  select(-c(.by))
observations<-analyzed_data%>%
  select(year,Y_rec,combos)%>%
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")
smooth_data <- na.omit(smooth_data )

partial_effects<-ggplot(smooth_data, aes(x = Value, y = .estimate)) +  # Setup the plot with the fitted values
  facet_wrap(~Smooth,ncol=1)+
  geom_line() + # Add estimated smooth
  xlim(c(-2,3))+
  geom_ribbon(aes(ymax = .upper_ci, ymin = .lower_ci), fill = "black", alpha = 0.2) + # Add credible interval
  geom_point(data = observations, aes(x = Value, y = Y_rec, col=year)) + # Add your data points
  labs(x = "Standardized Ecological Conditions", y = "Partial Residual")+ # Add labels
  geom_text(data = observations, aes(x = Value,col=year, y = Y_rec,label=year),hjust=0,nudge_x = 0.1)+
  scale_color_gradient(low="#0f85a0", high="#dd4124" )+
  theme_classic()
partial_effects 

png("Figures/FigureS7.png",width=6,height=8,units="in",res=1200)
partial_effects 
dev.off()
