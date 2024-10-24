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
df = data.frame(read.csv(paste0(data_dir,"DATA_Combined_glorys_yellowtail.csv"), header = T))

# get predictors to create model formula #######################################
envir_data = df %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('year','sd','Y_rec','ZOOpjuv','ZOOben')))%>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI')))%>% 
  mutate_all(~ scale(.))
head(envir_data)
dim(envir_data)
data_years = 1994:2014

envir_data2<- envir_data[complete.cases(envir_data), ]
data.pca <- princomp(envir_data2)
summary(data.pca)

# build model - full formula ###################################################
quadratic_terms = c('LSTlarv','ONIpre','ONIlarv', 'CutiSTIpjuv')

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
                    !(DDlarv && ONIlarv) && 
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
                    !(ONIlarv && PDOlarv) && 
                    !(ONIpre && ONIlarv) && 
                    #!(PDOlarv && ZOOpjuv) && 
                    !(Tpart && HCIlarv) && 
                    !(Tpart && PDOlarv) && 
                    #!(Tpart && ZOOpjuv) &&
                    dc(LSTlarv, I(LSTlarv^2)) &&
                    dc(ONIpre, I(ONIpre^2)) &&
                    dc(ONIlarv, I(ONIlarv^2)) &&
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

# Assess relative AIC support #########################################################

mtable4<-mtable
sum(mtable$weight)
mtable4[is.na(mtable4)] <- 0
param_names<-colnames(mtable2)[31:2]
AICcTab=NA
for(i in 2:31){
  AICweight =NA
  for(j in 1:nrow(mtable4)){
    AICweight = rbind(AICweight,ifelse(mtable4[j,i]!= 0, mtable4[j,37],0))
  }
  AICcTab=cbind(AICweight,AICcTab)
}


AICSum<-AICcTab[2:nrow(AICcTab),1:30]
colnames(AICSum)<-param_names

AICSum2<-data.frame(colSums(AICSum))
AICsummary<-data.frame(Predictor=param_names,AICWeight=round(AICSum2[,1],3))
print(AICsummary)
AICsummary10 <-AICsummary %>%
  arrange(desc(AICWeight)) %>%  # arrange in descending order
  slice(2:12)

which_names <- which(names(data_1) %in% AICsummary10$Predictor)
covs10AIC <- data_1[which_names]

covs10AIC <- covs10AIC[complete.cases(covs10AIC), ]
M = cor(covs10AIC)
colSums(abs(M))
corrplot.mixed(M, order = 'AOE')


M2 = cor(covs10AIC%>%select(!MLDpart&!ONIlarv&!PDOpjuv&!LSTlarv&!MLDlarv))
colSums(abs(M2))
corrplot.mixed(M2, order = 'AOE')

sat<-lm(Y_rec~CutiSTIpjuv+CutiSTIpjuv^2+DDben+Tcop+ONIpre+ONIpjuv, dat=data_1)
vif(sat)
# Refit LMs #########################################################
data_2<-data_1%>%select(CutiSTIpjuv,DDben,Tcop,ONIpre,ONIpjuv)
form_dredge = make_dredge_equation(envir_data = data_2, quadratic_vars = "CutiSTIpjuv")
form_dredge
fit = lm(form_dredge, data=data_1, na.action = na.fail)

# dredge_text
# write as a function to use in model testing later
fit_dredge <- function(fit){
  fit = fit
  mtable = dredge(fit, rank = AICc, m.lim = c(NA,5),
                  
                  subset=  dc(CutiSTIpjuv, I(CutiSTIpjuv^2)),
                  extra = list(R2 = function(x)
                    summary(x)$r.squared[[1]], 
                    F = function(y)
                      summary(y)$fstatistic[[1]]),
                  trace=2 )
  
  return(mtable)
}

mtable = fit_dredge(fit)

mtable4 = subset(mtable, delta<4)
print(mtable4, abbrev.names = FALSE)





# Testing out Marginal improvement ############################################
out = readRDS("parameteric_univariate.rds")
out = dplyr::filter(out, years_ahead == 1)

# add number of predictors
out$n_vars <- 3
out$n_vars[which(is.na(out$name3))] = 2
out$n_vars[which(is.na(out$name2))] = 1

df = expand.grid(spp = unique(out$spp),
                 model = unique(out$model),
                 cov = unique(out$name1))
d <- dplyr::group_by(out, name1) %>%
  dplyr::summarise(predictors = predictors[1]) %>%
  dplyr::rename(cov = name1)

df <- dplyr::left_join(df, d)

df$rmse_improve <- NA

for(i in 1:nrow(df)) {
  # for each row of the dataframe, subset original data to those with
  # the same model / species being predicted / 3 variables used
  sub_3vars_incl <- dplyr::filter(out, n_vars == 3,
                                  spp == df$spp[i],
                                  model == df$model[i],
                                  predictors == df$predictors[i])
  # and further filter to those models definitely including the predictor in 'cov'
  sub_3vars_incl <- dplyr::mutate(sub_3vars_incl,
                                  match1 = ifelse(name1==df$cov[i], 1, 0),
                                  match2 = ifelse(name2==df$cov[i], 1, 0),
                                  match3 = ifelse(name3==df$cov[i], 1, 0),
                                  all_matches = match1 + match2 + match3) %>%
    dplyr::filter(all_matches == 1) %>%
    dplyr::select(-match1, -match2, -match3, -all_matches)
  
  # re-order the time series names in name1/name2/name3 so 'cov' is last
  for(ii in 1:nrow(sub_3vars_incl)) {
    spp <- sub_3vars_incl[ii,c("name1","name2","name3")]
    spp_no_cov <- sort(as.character(spp[which(spp != df$cov[i])]))
    sub_3vars_incl[ii,c("name1")] <- spp_no_cov[1]
    sub_3vars_incl[ii,c("name2")] <- spp_no_cov[2]
    sub_3vars_incl[ii,c("name3")] <- df$cov[i]
  }
  
  # Next pull out the 2-covariate models with out the 'cov' predictor included
  sub_2vars_incl <- dplyr::filter(out, n_vars == 2,
                                  spp == df$spp[i],
                                  model == df$model[i],
                                  predictors == df$predictors[i])
  sub_2vars_incl <- dplyr::mutate(sub_2vars_incl,
                                  match1 = ifelse(name1==df$cov[i], 1, 0),
                                  match2 = ifelse(name2==df$cov[i], 1, 0),
                                  all_matches = match1 + match2) %>%
    dplyr::filter(all_matches == 0) %>%
    dplyr::select(-match1, -match2, -all_matches)
  # re-order the time series names in name1/name2
  for(ii in 1:nrow(sub_2vars_incl)) {
    spp <- sub_2vars_incl[ii,c("name1","name2")]
    sub_2vars_incl[ii,c("name1")] <- sort(as.character(spp))[1]
    sub_2vars_incl[ii,c("name2")] <- sort(as.character(spp))[2]
  }
  sub_2vars_incl <- dplyr::rename(sub_2vars_incl, rmse2 = rmse)
  sub_2vars_incl$model_id <- seq(1, nrow(sub_2vars_incl))
  # for each model in the dataframe with all three models, we need to match up the
  # corresponding 2-parameter model without the 'cov' predictor included
  sub_3vars_incl <- dplyr::left_join(sub_3vars_incl,
                                     sub_2vars_incl[,c("name1","name2","rmse2", "model_id")])
  
  # what's the percent improvement in rmse, e.g. 0.33 = 33% improvement  -- e.g. (1500 - 1000) / 1500
  sub_3vars_incl$pct_rmse <- (sub_3vars_incl$rmse2 - sub_3vars_incl$rmse) / sub_3vars_incl$rmse2
  # raw_mean <- mean(sub_3vars_incl$pct_rmse)
  
  # group by model id, calculate average rmse improvement for each id, then average across ids
  # mean of means is the same here -- so order of operations doesn't really matter, raw mean above works too
  mean_rmse <- dplyr::group_by(sub_3vars_incl, model_id) %>%
    dplyr::summarize(mean_rmse = mean(pct_rmse))
  df$rmse_improve[i] <- mean(mean_rmse$mean_rmse)
}

# all kidns of ways to sort this, e.g. by biggest to smallest rmse improvement by model and species?
df <- dplyr::filter(df, !is.na(df$rmse_improve))
df%>%filter(spp=='Yellowtail_rockfish_North')
