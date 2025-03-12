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

df = data.frame(read.csv("data-yellowtail/02_DATA_Combined_glorys_yellowtail2025_unexpanded.csv"))
dfa = data.frame(read.csv("data-yellowtail/dfa_trend.csv"))
data_years = 1994:2021 ## timeframe over which to stabilize environmental date
dat = df %>% 
  dplyr::select(!any_of( c('ZOOpjuv', 'ZOOben')))  %>%
  left_join(dfa)
envir_data = dat %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','year','Y_rec','ZOOpjuv','ZOOben','LUSI','ONIlarv',"dfa","Type")))%>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben',"X","dfa")))%>% 
  mutate_all(~ scale(.)) # drop terms not in model statement
full_dat<-cbind(dat%>%select(year, Y_rec, sd, Type),envir_data)
dat<-full_dat%>%  filter(year %in% data_years)

#### Evaluating Correlations between Covariates #####

threshold <-0.3 #assiging a threshold of correlation to filter out 
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
covariates <- names(envir_data)
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
  
  #fitting model from 1994 - 2013
  gam_model_train <- gam(as.formula(formula_str),
                         #data = dat)
                         data = dat[which(dat$year %notin% unique(dat$year)[jstart:n_year]), ])
  #fitting model from 1994 - 2018
   gam_model_full <- gam(as.formula(formula_str),
                        data = dat)
  #fitting model from 2013 - 2018
  gam_model_pred <- gam(as.formula(formula_str),
                        #data = dat)
                        data = dat[which(dat$year %notin% unique(dat$year)[jstart:n_year]), ])
  #saving what we care
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

# View the results data frame
print(results_arr_LFO5_rmse) #checkin LFO results
print(results_arr_LFO5_AIC) #checking AIC results
predicted_last5<-predicted # saving predictions


#### Predicting for recent years ####
# Fit a model that includes multiple variables from the top model
rankmod<-1 #which model rank do you want to look at?
#the combos determine which selection you want based on rankmod
combos= c(results_arr_LFO5_AIC$var1[rankmod], results_arr_LFO5_AIC$var2[rankmod], results_arr_LFO5_AIC$var3[rankmod], results_arr_LFO5_AIC$var4[rankmod])
smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ") #pasting in model structure
formula_str <- paste("Y_rec ~ ", smooth_terms) #formula for selected model
years<-1994:2019 #these are the training years, 2019 is a bit critical since LST goes off the rails
train_data<- full_dat%>%filter(year %in% years) 
gam <- gam(as.formula(formula_str),data = train_data)
summary(gam)
gam.predict <- cbind(train_data, predict.gam(newdata=train_data,gam,se.fit=TRUE), residuals= residuals(gam))
draw(gam)&
  theme_bw()
new_years<-2019:2024
new_dat<-full_dat%>%filter(year %in% new_years) 

predicted<-data.frame(predict.gam(gam,se.fit=TRUE, newdata=new_dat))%>%
  mutate(year=new_years)

gam.plot<- ggplot(gam.predict, aes(year, Y_rec)) +
  geom_point(aes(shape=Type)) +
  geom_point(data=new_dat) +
  #geom_point(data=predicted, col="red",aes(year, fit))+
  geom_line(data=predicted, col="red",aes(year, fit))+
  geom_point(data=predicted, col="red",aes(year, fit),shape=15)+
  geom_line(aes(year, fit)) +
  geom_ribbon(aes(ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  geom_ribbon(data=predicted, fill="red",aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  theme_bw() +
  #ylim(c(-1.5,2))+
  labs(title=paste("Deviance Explained = ",round(summary(gam)$dev.expl,2), ", R2 = ",round(summary(gam)$r.sq,2)), subtitle=formula(gam),
       x="Year", y="Recruitment Deviations")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.plot

png("Figures/ModelFit.png",width=10,height=6,units="in",res=1200)
gam.plot
dev.off()

gam.res <- ggplot(gam.predict,  aes(fit,residuals)) +
  geom_point() +
  #geom_line(aes(year, fit)) +
  theme_bw() +
  geom_smooth() +
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.res
png("Figures/ModelResiduals-unexpanded.png",width=6,height=5,units="in",res=1200)
gam.res
dev.off()

gam_rel<-draw(gam)&
  theme_bw()
gam_rel

png("Figures/ModelPartials-unexpanded.png",width=8,height=8,units="in",res=1200)
draw(gam)&
  theme_bw()
dev.off()

ts<-full_dat%>%select(c(combos,year))
tscov<-ggplot(pivot_longer(ts,col=combos,names_to = 'var',values_to = "val"),aes(x=year,y=val))+
  geom_line()+
  geom_point()+
  facet_wrap(~var)+
  geom_hline(yintercept = 0,lty=2)+
  theme_bw()
tscov

png("Figures/OceanographicTS.png",width=8,height=8,units="in",res=1200)
tscov
dev.off()



#### GAM diagnostics ####

#checking for multicollinearity
data.frame(cor(na.omit(full_dat %>% dplyr::select(combos)))) # drop terms not in model statement
M<-cor(na.omit(full_dat %>% dplyr::select(combos))) # drop terms not in model statement
corrplot.mixed(M)
png("Figures/CorrPlot.png",width=8,height=8,units="in",res=1200)
corrplot.mixed(M)
dev.off()


combos= c(results_arr_LFO5_AIC$var1[rankmod], results_arr_LFO5_AIC$var2[rankmod], results_arr_LFO5_AIC$var3[rankmod], results_arr_LFO5_AIC$var4[rankmod])
linear_terms <- paste(combos, collapse = " + ") #pasting in model structure
formula_str_lm <- paste("Y_rec ~ ", linear_terms) #formula for selected model
vifs<-vif(lm(as.formula(formula_str_lm),data = train_data))
vifs.df<-data.frame(vifs)

results_arr%>%
  select(rankmod, rankmod_loo,var1, var2, var3, var4,AIC,deltaAIC,RMSE_LOO,rsq_full,dev.ex_full)%>%

write.csv(vifs.df,"VIFs.csv")
write.csv(M,"Correlation.csv")

#looking at qq, resid, and response with all rec devs
diag_dat<-full_dat%>%filter(year %in% 1994:2021) 
m2 = gam(as.formula(formula_str),
data = diag_dat)
png("Figures/gamdiag-unexpanded.png",width=6,height=6,units="in",res=1200)
par(mfrow=c(2,2))
gam.check(m2,pch=19,cex=.3)
dev.off()
plotDiagnostics(m2)
#looking for leverage years
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
png("Figures/Cooks.png",width=5,height=5,units="in",res=1200)
cooks
dev.off()

#### Predictive Capacity Validation ####


mn = apply(results_arr_LFO5_AIC%>%select(RMSE),2,mean)
md = apply(results_arr_LFO5_AIC%>%select(RMSE),2,median)
# quantile function
qt = function(x){quantile(x,c(0.025,0.50,0.975))}

ci = apply(results_arr_LFO5_AIC%>%select(mape),2, qt)

ci_RMSE = apply(results_arr_LFO5_AIC%>%select(RMSE),2, qt)
rmse_plot<-ggplot(results_arr_LFO5_AIC%>%select(RMSE), aes(RMSE)) +
  geom_histogram(bins=10,fill="lightgrey", col="black")+
  #xlim(c(0.4,0.8))+
  ylab("Frequency")+
  geom_vline(xintercept=c(ci_RMSE), lty=2,lwd=1)+
  geom_vline(xintercept=results_arr_LFO5_AIC$RMSE[1], col='red',lwd=1)+
  xlab("RMSE")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ci_MAPE = apply(results_arr_LFO5_AIC%>%select(mape),2, qt)
mape_plot<-ggplot(results_arr_LFO5_AIC%>%select(mape), aes(mape)) +
  geom_histogram(bins=10,fill="lightgrey", col="black")+
  #xlim(c(0.4,0.8))+
  ylab("Frequency")+
  geom_vline(xintercept=c(ci_MAPE), lty=2,lwd=1)+
  geom_vline(xintercept=results_arr_LFO5_AIC$mape[1], col='red',lwd=1)+
  xlab("MAPE")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ci_R2 = apply(results_arr_LFO5_AIC%>%select(rsq_pred),2, qt)
r2_plot<-ggplot(results_arr_LFO5_AIC%>%select(rsq_pred), aes(rsq_pred)) +
  geom_histogram(bins=10,fill="lightgrey", col="black")+
  #xlim(c(0.4,0.8))+
  ylab("Frequency")+
  geom_vline(xintercept=c(ci_R2), lty=2,lwd=1)+
  geom_vline(xintercept=results_arr_LFO5_AIC$rsq_pred[1], col='red',lwd=1)+
  xlab("R-squared (predicted time period)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ci_dev = apply(results_arr_LFO5_AIC%>%select(dev.ex_pred),2, qt)
dev_plot<-ggplot(results_arr_LFO5_AIC%>%select(dev.ex_pred), aes(dev.ex_pred)) +
  geom_histogram(bins=10,fill="lightgrey", col="black")+
  #xlim(c(0.4,0.8))+
  ylab("Frequency")+
  geom_vline(xintercept=c(ci_dev), lty=2,lwd=1)+
  geom_vline(xintercept=results_arr_LFO5_AIC$dev.ex_pred[1], col='red',lwd=1)+
  xlab("Deviance Explained (predicted time period)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

predcap<-ggarrange(rmse_plot,mape_plot,r2_plot,dev_plot)
png("Figures/PredCap.png",width=8,height=6,units="in",res=1200)
predcap
dev.off()
### Boot Strapping ####
fm = gam
sfm = summary(fm)
#fmedf = t(data.frame(sfm$edf))
Fstats = t(data.frame(c(sfm$s.table[1:4,3])))
Pstats = t(data.frame(c(sfm$s.table[1:4,4])))
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
colnames(results_boot) = c('r2','devex',"P1","P2",
                                "P3","P4","F1","F2","F3","F4")
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
  rsq<-summary(gam_model)$r.sq
  rsqjack<-rbind(rsqjack,rsq)
  }




gam.jack<-cbind(gam.predict,jack=predictions[1:26],r2=rsqjack[1:26,1])
ggplot(gam.jack, aes(year, Y_rec)) +
  #geom_point(aes(shape=Type)) +
  geom_point(data=new_dat) +
  geom_point(aes(y=jack),bg='yellow',pch=21,col="black") +
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
  labs(title="Jackknife",
       x="Year", y="ln(Recruitment Deviations)")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggplot(gam.jack, aes(r2)) +
  geom_histogram(bins=10,fill="lightgrey", col="black")+
  xlim(c(0.4,0.8))+
  ylab("Frequency")+
  xlab("R-squared")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(gam.jack, aes(x=year, y=r2)) +
  geom_point(col="black")+
 # xlim(c(0.4,0.8))+
  ylim(c(0,1))+
  ylab("R-squared")+
  geom_hline(yintercept=0.59,lty=2)+
  xlab("Year (removed)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

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
