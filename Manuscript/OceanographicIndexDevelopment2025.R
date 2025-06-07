library(dplyr)
library(tidyr)
library(ggpubr)
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
library(ggpubr)

##### Reading in the data #####
df <- data.frame(read.csv("data-yellowtail/2024Env-annual-yellowtail-standardized.csv"))
colnames(df)
envir_data = df %>% 
  rename(ONIlarv=oni_larv,ONIpjuv=oni_pjuv,ONIpre=oni_pre,
         PDOlarv=pdo_larv,PDOpjuv=pdo_pjuv,
         LUSI=lusi_annual,
         BeutiTUMI=BeutiTUMIpjuv, BeutiSTI=BeutiSTIpjuv,
         CutiSTI=CutiSTIpjuv,CutiTUMI=CutiTUMIpjuv, 
         BakunSTI=bakun_sti,
         HCIlarv=hci1_larv,HCIpjuv=hci1_pjuv,
         HCI2larv=hci2_larv,HCI2pjuv=hci2_pjuv,
         NCOPpjuv=ZOOpjuvN, NCOPben=ZOObenN,
         SCOPpjuv=ZOOpjuvS, SCOPben=ZOObenS)
to_omit<-c('X',"Year","Year.1","LUSI", "BakunSTI", "HCI2larv", "HCI2pjuv", "HCIlarv", "HCIpjuv")
envir_data2 =envir_data%>%dplyr::select(!any_of(to_omit))
recruitmentdevs<- data.frame(read.csv("data-yellowtail/RecruitmentDeviations2025draft.csv"))%>%
  filter(Datatreatment=="Expanded PacFIN")
df = envir_data2%>%left_join(recruitmentdevs)
data_years = 1998:2021 ## timeframe over which to stabilize environmental date
colnames(df)
dat = df %>% 
  dplyr::select(!any_of(to_omit))

full_dat<-dat
dat<-full_dat%>%  filter(year %in% data_years)

#### Evaluating Correlations between Covariates #####

threshold <-0.3 #assiging a threshold of correlation to filter out 
envir_data2 <- envir_data2[complete.cases(envir_data2), ] %>% 
  dplyr::select(!any_of( c('year'))) # getting environmental data
M = data.frame(cor(envir_data2)) # generating a correlation matrix
M <- tibble::rownames_to_column(M, "xvar") #putting the rownames in their own column
M<-M%>%pivot_longer(!xvar, names_to = 'yvar', values_to = 'corr') #pivoting to a longer table
uncorr<-M%>%filter(abs(corr)<threshold) #generating the uncorrelated thresholds 
corrvar<-M%>%filter(abs(corr)>threshold)%>%filter(xvar!=yvar) #generating the correlated thresholds
corrvar2<-corrvar%>%select(-corr)
combinations_to_omit<-list() #setting up an empty list to fill with pairs to omit
for(i in 1:nrow(corrvar2)){
  combinations_to_omit[[i]]<-c(as.character(corrvar2[i, ]))
}

#### Generate Null Model ####
predict(gam(Y_rec ~ 1, data=dat))

null<-dat%>%mutate(sr_null = 0, mean_null = -0.0917625)

sqerror<-(null$sr_null-null$Y_rec)^2
rmse_sr_full <- sqrt(mean(sqerror, na.rm=T))

sqerror<-(null$mean_null-null$Y_rec)^2
rmse_mean_full <- sqrt(mean(sqerror, na.rm=T))

null_lfo<-full_dat%>%mutate(sr_null = 0, mean_null = -0.0917625)%>%filter(year>=2014 & year<2019)

sqerror<-(null_lfo$sr_null-null_lfo$Y_rec)^2
rmse_sr_lfo <- sqrt(mean(sqerror, na.rm=T))
sqerror<-(null_lfo$mean_null-null_lfo$Y_rec)^2
rmse_mean_lfo <- sqrt(mean(sqerror, na.rm=T))

#### Models to Fit ####
covariates <- names(envir_data2)
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

#### LOO ####
models <- list()
results <- data.frame()
predicted <- data.frame()
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
  predict_mod <- predict(gam_model)
  # keep in mind RMSE is weird for binomial things
  rmse_loo <- sqrt(mean((dat$Y_rec - predictions)^2, na.rm=T))
  rmse <- sqrt(mean((dat$Y_rec -  predict_mod)^2, na.rm=T))
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
    RMSE_loo = round(rmse_loo,3),
    RMSE = round(rmse,3),
    rsq_full=round(r2,2),
    dev.ex=round(dev.expl,4),
    rmse_imp=(rmse_sr_full-rmse)/rmse_sr_full, 
    rmse_ratio=rmse/rmse_sr_full,
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
  #saving the one step ahead predictions
  print(i)
}
results <- arrange(results,AIC)%>%
  mutate(rankmod=seq_along(ModelID))%>%
  mutate(deltaAIC=AIC-min(AIC))
results_arr_RMSE<- arrange(results,RMSE_loo)
results_arr_AIC<- arrange(results,AIC)
arrange(results,desc(rmse_imp))
results_arr_AIC%>%filter(rmse_ratio<1)
arrange(results,rmse_ratio)
### Calculating Relative Variable Importance ###

# Create a baseline model with only the intercept

baseline_rmse <-rmse_sr_full

# we can calculate the marginal improvement for each covariate
results_arr_RMSE$n_cov <- ifelse(!is.na(results_arr_RMSE$var1), 1, 0) + ifelse(!is.na(results_arr_RMSE$var2), 1, 0) +
  ifelse(!is.na(results_arr_RMSE$var3), 1, 0) + ifelse(!is.na(results_arr_RMSE$var4), 1, 0)
marginals <- data.frame(cov = covariates, "rmse_01" = NA, "rmse_12" = NA, "rmse_23" = NA,
                        "aic_01" = NA, "aic_12" = NA, "aic_23" = NA)
results<-results_arr_RMSE%>%filter(n_cov<4)%>%select(-var4)
results%>%filter(n_cov==1)
for(i in 1:length(covariates)) {
  sub <- dplyr::filter(results, n_cov == 1,
                       var1 == covariates[i])
   marginals$rmse_01[i] <- (baseline_rmse-sub$RMSE) / baseline_rmse
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
    sub2$rmse_diff[j] <- (sub1$RMSE[indx]-sub2$RMSE[j]) / sub1$RMSE[indx]
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
    sub3$rmse_diff[j] <- (sub2_all$RMSE[indx]-sub3$RMSE[j]) / sub2_all$RMSE[indx]
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

included <- data.frame(cov=c("CutiSTI","LSTpjuv","ONIpjuv","MLDpart","DDegg"),
           Included = c("Yes"))

gam_loo_table <- dplyr::arrange(marginals, total_rmse)%>%
  dplyr::select(cov, total_rmse, total_aic)%>%
  mutate(cv="LOO",model="GAM")%>%
  left_join(included)
gam_loo_table[is.na(gam_loo_table)] <-"No"
cols<- c('#dd4124',"#edd746",'#7cae00','#0f85a0')

marginal <- ggplot(gam_loo_table, aes(x = reorder(cov,total_rmse), y = total_rmse, fill=Included)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the axes to make a horizontal bar graph
  labs(y = "Marginal Improvement RMSE", x = "Predictor")+
  scale_fill_manual(values=c('grey',cols[3]))+
  theme_classic()
marginal 

rmse_ratio <- ggplot(results_arr_RMSE, aes(x = n_cov, y = rmse_ratio)) +
  geom_point(position = "jitter",  alpha=0.4) +
  geom_hline(yintercept=1, lty=2)+
  labs(x = "Number of Predictors", y = expression("RMSE/RMSE"["null"]),col="Number of \n Predictors")+
  #scale_colour_manual(values=cols)+
  ylim(c(0.25,1))+
  theme_classic()
rmse_ratio


png("Figures/Figure7.png",width=8,height=5,units="in",res=1200)
ggarrange(rmse_ratio, marginal, labels = c("A.", "B."), widths = c(0.75,1))
dev.off()
#### Extracting LOO-CV Highest Ranked Model ####
rankmod<-1 #which model rank do you want to look at?
combos= c(results_arr_RMSE$var1[rankmod], results_arr_RMSE$var2[rankmod], results_arr_RMSE$var3[rankmod], results_arr_RMSE$var4[rankmod]) #this needs to be mod depending on which model and selection criteria you want

smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ") #pasting in model structure
formula_str <- paste("Y_rec ~ ", smooth_terms) #formula for selected model
years<-1998:2021
gam_data<- full_dat%>%filter(year %in% years) 
gam <- gam(as.formula(formula_str),data = gam_data)

#### Predicting the last 5 years ####

yearstrain<- 1998:2016
yearspredict<- 2017:2021

train_data<- full_dat%>%filter(year %in% yearstrain) 
gam_train <- gam(as.formula(formula_str),data = train_data)
gam5 <- gam(as.formula(formula_str),data = train_data)
predlast5<-predict.gam(gam5,newdata=full_dat%>%filter(year %in% yearspredict))
gam2 <- gam(as.formula(formula_str),data = full_dat%>%filter(year %in% 1998:2019))
predlast2<-predict.gam(gam2,newdata=full_dat%>%filter(year %in% 2020:2021))
#gam.train<-data.frame(year=train_data$year,predict(gam_train,se=TRUE))

gam3 <- gam(as.formula(formula_str),data = full_dat%>%filter(year %in% 1998:2021))
gam.train<-data.frame(year=full_dat%>%filter(year %in% 1998:2021)%>%select(year),predict(gam3,se=TRUE))

train_dat_pred <- full_dat%>%left_join(gam.train%>%select(year, fit, se.fit))%>%
  left_join(predicted%>%filter(ModelID==616)%>%select(year,pred))%>%
  left_join(data.frame(year=yearspredict,
                       last5 =predlast5))%>%
  left_join(data.frame(year=2020:2021,
                       last2 =predlast2))%>%
  mutate(type2=ifelse(type=="Main_RecrDev","Main", "Late"))


#### Last 5 1-year at a time #####
#start by dredging all possible combinations of the GAMs
pred5<- data.frame()
startval<-1
j1<- 2017
nyear <- 5
jstart<- 21
jend=jstart+nyear

  for (j in jstart:jend) {
    train_index <- 1998:2017  # All indices except the j-th
    test_index <- j                # The j-th index
    
    # Fit model on n-1 observations
    gam_model <- gam(as.formula(formula_str),
                     # weights = number_cwt_estimated,
                     data = dat[which(dat$year %in% unique(dat$year)[startval:(j-1)]), ])
    newdata<-dat[which(dat$year == unique(dat$year)[j]), ]
    # Predict the excluded observation
    pred5 <- rbind(pred5,cbind(predict(gam_model, newdata = newdata), newdata$year))
  }
colnames(pred5)<-c("pred5", "year")
train_dat_pred<-left_join(train_dat_pred,pred5)
#### Plotting Highest Ranked Model ####
# Fit a model that includes multiple variables from the top model
sizept=1.75
gam.plot<- ggplot(train_dat_pred%>%filter(year  %in% 1998:2021), aes(x=year)) +
  geom_line(data=train_dat_pred,aes(y=fit)) +
  geom_ribbon(data=train_dat_pred,aes(x=year,y=fit,ymax=fit+2.1*se.fit, ymin=fit-2.1*se.fit), 
              alpha=0.2)+
  geom_point(aes(y=Y_rec,fill="Recruitment \n Deviation",shape=type2),size=sizept) +
  geom_point(aes(y=pred,fill= "Leave-one-out",shape=type2), pch=21,size=sizept ) +
 # geom_point(aes(y=last2, fill = "Last 2",shape=type2), pch=24,size=sizept) +
  geom_point(data=train_dat_pred%>%filter(year>2018),aes(y=last5, fill = "Last 5",shape=type2), pch=24,size=sizept) +
  geom_point(data=train_dat_pred%>%filter(year<2019),aes(y=last5, fill = "Last 5",shape=type2), pch=21,size=sizept) +
  geom_point(data=train_dat_pred%>%filter(year>2018),aes(y=pred5, fill = "One year ahead",shape=type2), pch=24,size=sizept) +
  geom_point(data=train_dat_pred%>%filter(year<2019),aes(y=pred5, fill = "One year ahead",shape=type2), pch=21,size=sizept) +
  scale_shape_manual(name="Recruitment Deviation \n Type",values=c(24,21))+
  scale_fill_manual(name = "Diagnostic",
                    values = c("Last 5" = '#dd4124',
                               "Leave-one-out" = '#0f85a0',
                               "One year ahead" = "#edd746",
                              "Recruitment \n Deviation"= "black"
                              ), 
                    labels = c("Last 5","Leave-one-out","One year ahead", "Recruitment \n Deviation"))+
  theme_classic() +
  labs(x="Year", y="ln(Recruitment Deviations)")+
  ylim(c(-2,1.5))+
  xlim(c(1998,2022))+
  #labs(fill = "Diagnostic", shape="Recruitment Deviation \n Type")+
  geom_hline(yintercept=0,lty=2)+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
gam.plot

png("Figures/Figure4.png",width=8,height=5,units="in",res=1200)
gam.plot
dev.off()



#### Partial Effects Figures ####
# Use smooth_estimates() to get the smooth estimates and confidence intervals
gam <- gam(as.formula(formula_str),data = gam_data%>%filter(year<2019))
smooth_data <- smooth_estimates(gam)%>%  # Or specify the smooth term you want to plot
  add_constant(model_constant(gam)) %>% # Add the intercept
  add_confint()%>% # Add the credible interval
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")%>%
  select(-c(.by))
observations<-full_dat%>%
  select(year,Y_rec,combos)%>%
  pivot_longer(c(combos),names_to = "Smooth",values_to ="Value")
smooth_data <- na.omit(smooth_data )

partial_effects<-ggplot(smooth_data, aes(x = Value, y = .estimate)) +  # Setup the plot with the fitted values
  facet_wrap(~Smooth)+
  geom_line() + # Add estimated smooth
  xlim(c(-2,3))+
  geom_ribbon(aes(ymax = .upper_ci, ymin = .lower_ci), fill = "black", alpha = 0.2) + # Add credible interval
  geom_point(data = observations%>%filter(year %in% years), aes(x = Value, y = Y_rec), color = "black") + # Add your data points
  labs(x = "Standardized Ecological Conditions", y = "Partial Residual")+ # Add labels
  geom_text(data = observations%>%filter(year %in% years), aes(x = Value, y = Y_rec,label=year),hjust=0,nudge_x = 0.1)+
  theme_classic()
partial_effects 
png("Figures/Figure5.png",width=8,height=8,units="in",res=1200)
partial_effects 
dev.off()

tsdata<-pivot_longer(full_dat,cols=-c(year,Datatreatment,type,sd,Y_rec),names_to = 'var',values_to = "val")%>%
  select(year, var, val,Y_rec)%>%
  mutate(type = ifelse(var=='CHLpjuv'|var=='PPpjuv'|var=="NCOPben"|var=="NCOPpjuv"|var=="SCOPben"|var=="SCOPpjuv", "Foodweb", "Oceanographic"))

ts_foodweb <- ggplot(tsdata%>%filter(type=="Foodweb"),aes(x=year,y=val))+
  geom_line()+
  geom_point()+
  facet_wrap(~var, ncol=6)+
  ylab("")+
  xlab("")+
  ggtitle("Foodweb Conditions")+
  geom_hline(yintercept = 0,lty=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))


ts_ocean <-ggplot(tsdata%>%filter(type=="Oceanographic"),aes(x=year,y=val))+
  geom_line()+
  geom_point()+
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

png("Figures/Figure2.png",width=10,height=8,units="in",res=1200)
full_ts
dev.off()


#### GAM diagnostics ####

#checking for multicollinearity
dfx = envir_data2#%>%select(-c(HCIpjuv, HCIlarv, HCI2larv, HCI2pjuv)) 
dfx = na.omit(dfx) # na's mess up the cor function
cor_xy = cor(dfx)

graphics.off()
png( paste0("Figures/FigureA1.png"),
     units = 'in', res=300, width = 5, height=5)
# plot command
corrplot::corrplot(cor_xy, hod='circle', type='lower', tl.cex = 0.5)
dev.off()

m2<-gam
diag_dat<-full_dat%>%filter(year %in% 1998:2021) 
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
png("Figures/CooksRMSE.png",width=5,height=5,units="in",res=1200)
cooks
dev.off()

#### Jackknifing ####
#start by dredging all possible combinations of the GAMs
results_loo <- data.frame()
jstart<-1

predictions <- numeric(nrow(dat))
rsqjack <- data.frame()
n_year <- length(unique(dat$year))

# LOO cv
dat<-gam_data
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

gam.jack<-cbind(gam_data,jack=predictions[1:24],r2=rsqjack[1:24,1])

#### R2 comparison ####

r2jack<-ggplot(gam.jack, aes(r2)) +
  geom_histogram(bins=40,fill="lightgrey", col="black")+
  xlim(c(0,0.8))+
  ylab("Frequency")+
  xlab("R-squared")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

r2jack2<-ggplot(gam.jack, aes(x=year, y=r2)) +
  geom_point(col="black")+
 # xlim(c(0.4,0.8))+
  ylim(c(0,1))+
  ylab("R-squared")+
  geom_hline(yintercept=mean(na.omit(gam.jack$r2)),lty=2)+#0.59
  xlab("Year (removed)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

png("Figures/Figure5.png",width=4,height=6,units="in",res=1200)
ggarrange(r2jack,r2jack2, ncol=1)
dev.off()


png("Figures/A3?.png",width=6,height=6,units="in",res=1200)
par(mfrow=c(2,2))
gam.check(gam,pch=19,cex=.3)
dev.off()

rec_foodweb <- ggplot(tsdata%>%filter(type=="Foodweb"),aes(x=val,y=Y_rec))+
  geom_smooth(method='gam',formula = y ~ s(x,  k=3), col='#7cae00')+
  geom_point()+
  facet_wrap(~var, ncol=6)+
  ylab("")+
  xlab("")+
  ggtitle("Foodweb Conditions")+
  
  geom_hline(yintercept = 0,lty=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))


rec_ocean <-ggplot(tsdata%>%filter(type=="Oceanographic"),aes(x=val,y=Y_rec))+
  geom_smooth(method='gam',formula = y ~ s(x,  k=3), col='#7cae00')+
  
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
png("Figures/A2.png",width=10,height=8,units="in",res=1200)
full_rec
dev.off()


##### Generating pred error for the training period #####
yearmin<- 2019
yearmax<-2021

train_data<- full_dat%>%filter(year %in% years) 
gam <- gam(as.formula(formula_str),data = train_data)
train_data5<- full_dat%>%filter(year %in% years5) 
gam5 <- gam(as.formula(formula_str),data = train_data5)
summary(gam)
gam.train <- cbind(train_data, predict.gam(newdata=train_data,gam,se.fit=TRUE), residuals= residuals(gam))
Vb <- vcov(gam)
N <- 10000
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)),V = Vb)
Cg <- predict(gam, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)
absDev <- abs(sweep(simDev, 1, gam.train$se.fit, FUN = "/"))
masd <- apply(absDev, 2L, max)
crit2 <- quantile(masd, prob = 0.95, type = 8)
train_dat_pred <- full_dat%>%left_join(gam.train%>%select(year, fit, se.fit))%>%
  mutate(uprP = fit + (crit2 * se.fit),
         lwrP = fit - (crit2 * se.fit),
         uprCI = fit + (2 * se.fit),
         lwrCI = fit - (2 * se.fit))%>%
  mutate(se.p=(uprP-fit)/1.96)%>%
  left_join(predicted%>%filter(ModelID==616)%>%select(year,pred))%>%
  left_join(data.frame(year=seq(yearmin,yearmax),
                       last5 = predict.gam(gam5,newdata=train_dat_pred%>%filter(year>=yearmin&year<=yearmax)))
  )%>%
  mutate(type2=ifelse(type=="Main_RecrDev","Main", "Late"))
