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
library(gratia)

##### Reading in the data #####

df = data.frame(read.csv("data-yellowtail/02_DATA_Combined_glorys_yellowtail2025.csv"))
dfa = data.frame(read.csv("data-yellowtail/dfa_trend.csv"))
data_years = 1994:2018
dat = df %>% 
  dplyr::select(!any_of( c('ZOOpjuv', 'ZOOben')))  %>% 
  filter(year %in% data_years)%>%
  left_join(dfa)
envir_data = dat %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','year','Y_rec','ZOOpjuv','ZOOben',"dfa")))%>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI','ONIlarv',"X","dfa")))%>% 
  mutate_all(~ scale(.)) # drop terms not in model statement
dat<-cbind(dat%>%select(year, Y_rec, sd),envir_data)

#### Evaluating Correlations between Covariates #####

threshold <-0.4 #assiging a threshold of correlation to filter out 
envir_data <- envir_data[complete.cases(envir_data), ] # getting environmental data
M = data.frame(cor(envir_data)) # generating a correlation matrix
M <- tibble::rownames_to_column(M, "xvar") #putting the rownames in their own column
M<-M%>%pivot_longer(!xvar, names_to = 'yvar', values_to = 'corr') #pivoting to a longer table
uncorr<-M%>%filter(abs(corr)<threshold) #generating the uncorrelated thresholds 
corrvar<-M%>%filter(abs(corr)>threshold)%>%filter(xvar!=yvar) #generating the correlated thresholds
corrvar2<-corrvar%>%select(-corr)
combinations_to_omit<-list() #setting up an empty list to fill with pairs to omit
for(i in 1:nrow(corrvar2)){
  combinations_to_omit[[i]]<-c(as.character(corrvar2[i, ]))
}

#### Models to Fit ####

maxcovar <- 4 #max number of covars for a given model
combinations <- lapply(1:maxcovar, function(i) {
  combn(covariates, i, simplify = FALSE)
})

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
length(combinations) # how many left? only 1900

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

for (i in seq_along(combinations)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  #smooth_terms <- paste("s( total_release, k =3) + s(", covariates[[i]], ", k = 3)")
  smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ")
  formula_str <- paste("Y_rec ~ ", smooth_terms)
  predictions <-  numeric(n_pred )

  # Loop over each observation
  for (j in jstart:n_year) {
    train_index <- setdiff(train_start:n_train, j)  # All indices except the j-th
    test_index <- j                 # The j-th index
    
    # Fit model on n-j observations
    gam_model <- gam(as.formula(formula_str),
                     data = dat[which(dat$year %notin% unique(dat$year)[j:n_year]), ])
    
    # Predict the excluded observation
    predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    # re-fit the model
    
  }
  # re-fit the model
  gam_model_train <- gam(as.formula(formula_str),
                         #data = dat)
                         data = dat[which(dat$year %notin% unique(dat$year)[jstart:n_year]), ])
  gam_model_full <- gam(as.formula(formula_str),
                        data = dat)
  gam_model_pred <- gam(as.formula(formula_str),
                        #data = dat)
                        data = dat[which(dat$year %notin% unique(dat$year)[jstart:n_year]), ])
  # keep in mind RMSE is weird for binomial things
  rmse <- sqrt(mean((dat$Y_rec[jstart:n_year] - predictions[jstart:n_year])^2, na.rm=T))
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
    pred=predictions[jstart:n_year],
    year=unique(dat$year)[jstart:n_year],
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3],
    var3 = padded_vars[4]
    
    
  ))
  
  print(i)
}
results2<-results
results<-results2%>%filter(rsq_full>0&rsq_train>0&rsq_pred>0)
results_arr_LFO5_rmse <- arrange(results,RMSE) #ordered results by RMSE
results_arr_LFO5_AIC <- arrange(results,AIC) #ordered results by AIC
results_arr_LFO5_devex <- arrange(results,desc(dev.ex_pred)) #ordered results by deviance explained
# View the results data frame
print(results_arr_LFO5_rmse)
predicted_last5<-predicted # saving predictions

#### Predicting for recent years #####
# Fit a model that includes multiple variables from the top model
rankmod<-1
combos= c(results_arr_LFO5_rmse$var1[rankmod], results_arr_LFO5_rmse$var2[rankmod], results_arr_LFO5_rmse$var3[rankmod], results_arr_LFO5_rmse$var4[rankmod])
smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ")
formula_str <- paste("Y_rec ~ ", smooth_terms)
gam <- gam(as.formula(formula_str),
           #family = binomial(),
           #weights = number_cwt_estimated,
           data = full_dat%>%filter(year>1993&year<2021))
summary(gam)

new_years<-2019:2023
years<-1994:2023

envir_data = df%>%  filter(year %in% years) %>% 
  dplyr::select(!any_of(c('sd','year','Y_rec','ZOOpjuv', 'ZOOben','LUSI','ONIlarv',"X")))%>%  # drop terms that are not covariates
  mutate_all(~ scale(.)) # scale covariates anbout a mean of 0 and sd of 1
full_dat<-cbind(df%>%  filter(year %in% years)%>%select(year, Y_rec, sd),envir_data)
new_dat<-full_dat%>%filter(year %in% new_years) 
old_data<- dat
predicted<-data.frame(predict.gam(gam,se.fit=TRUE, newdata=new_dat))%>%
  mutate(year=new_years)
#plot.gam(gam)


gam.predict <- cbind(full_dat%>%filter(year>1993&year<2021), predict.gam(gam,se.fit=TRUE), residuals= residuals(gam))

gam.plot<- ggplot(gam.predict, aes(year, Y_rec)) +
  geom_point() +
  geom_point(data=new_dat) +
  geom_line(data=predicted, col="red",aes(year, fit))+
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  geom_ribbon(data=predicted, fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  theme_bw() +
  labs(title=paste("gam, LFO-5: Deviance Explained = ",round(summary(gam)$dev.expl,2)), subtitle=formula(gam),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.plot
plot.gam(gam)



gamloo.res <- ggplot(gamloo.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



#### Plotting Recent Years ####
draw(gam)&
  theme_bw()
ts<-full_dat%>%select(c(combos,year))
ggplot(pivot_longer(ts,col=combos,names_to = 'var',values_to = "val"),aes(x=year,y=val))+
  geom_line()+
  geom_point()+
  facet_wrap(~var)+
  geom_hline(yintercept = 0,lty=2)+
  theme_bw()
