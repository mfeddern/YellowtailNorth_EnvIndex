library(tidyverse)
library(ggplot2)
library(readr) # faster writing
library(data.table) # faster writing of large files
library(lubridate)
library(mgcv)
library(MuMIn)
library(tidyverse)
library(ggplot2)
library(readr) # faster writing
library(data.table) # faster writing of large files
library(lubridate)
library(mgcv)
library(MuMIn)
library(corrplot)
library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
library(stringr)
library(ggpubr)
library(car)
source("Src/_00_yellowtail-header.r")
source("Src/Functions-for-envir-index.r")
# set directories ##############################################################

# bring in data ################################################################
# combined fish and environmental drivers file
df = data.frame(read.csv("data-yellowtail/02_DATA_Combined_glorys_yellowtail.csv", header = T))

# get predictors to create model formula #######################################
envir_data = df %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('year','sd','Y_rec','ZOOpjuv','ZOOben')))
head(envir_data)
dim(envir_data)
data_years = 1994:2014

envir_data <- envir_data[complete.cases(envir_data), ]
M = data.frame(cor(envir_data))
M <- tibble::rownames_to_column(M, "xvar")
M<-M%>%pivot_longer(!xvar, names_to = 'yvar', values_to = 'corr')
uncorr<-M%>%filter(corr<0.4&corr> -0.4)
corrplotdat<-cor(envir_data)
unique(uncorr$xvar)
corrplot.mixed(corrplotdat)
data_1 = df %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','ZOOpjuv','ZOOben')))
head(envir_data)
dim(envir_data)
data_years = 1994:2014

  the_lm <-  lm(Y_rec ~ DDpre + DDegg + DDlarv + DDpjuv + DDben + Tcop + Tpart + 
    MLDpart + MLDlarv + MLDpjuv + CSTlarv + CSTpjuv + LSTlarv + 
    LSTpjuv + HCIlarv + HCIpjuv + ONIpre + ONIlarv + ONIpjuv + 
    PDOlarv + PDOpjuv + CutiSTIpjuv + CutiTUMIpjuv + BeutiTUMIpjuv+bakun_sti, data=data_1)
  #Organize SCC biology data
dat<- df %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI')))
dat_red <- dat[complete.cases(dat), ]
dat_scale <-  dat_red%>%
    dplyr::select(!any_of(c('year'))) %>% 
    mutate_all(~ scale(.))%>%
  cbind(year=dat_red$year)
dat<-dat_scale
ax<- 20
ti<-24
wid <- 28
plot_trends2 <- function(rotated_modelfit,
                         years = NULL,
                         highlight_outliers = FALSE,
                         threshold = 0.01) {
  rotated <- rotated_modelfit
  df <- dfa_trends(rotated, years = years)
  
  # make faceted ribbon plot of trends
  p1 <- ggplot(df, aes_string(x = "time", y = "estimate")) +
    geom_ribbon(aes_string(ymin = "lower", ymax = "upper"), alpha = 0.4) +
    geom_line() +
    facet_wrap("trend_number") +
    xlab("Time") +
    ylab("")+
    
    theme(axis.text=element_text(size=ax),
          axis.title=element_text(size=ti,face="bold"))+
    theme_bw()
  
  if (highlight_outliers) {
    swans <- find_swans(rotated, threshold = threshold)
    df$outliers <- swans$below_threshold
    p1 <- p1 + geom_point(data = df[which(df$outliers), ], color = "red")
  }
  
  p1
}


plot_loadings2 <- function(rotated_modelfit,
                           names = NULL,
                           facet = TRUE,
                           violin = TRUE,
                           conf_level = 0.95,
                           threshold = NULL) {
  v <- dfa_loadings(rotated_modelfit,
                    summary = FALSE,
                    names = names,
                    conf_level = conf_level
  )
  df <- dfa_loadings(rotated_modelfit,
                     summary = TRUE,
                     names = names,
                     conf_level = conf_level
  )
  
  # filter values below threshold
  if (!is.null(threshold)) {
    df <- df[df$prob_diff0 >= threshold, ]
    v <- v[v$prob_diff0 >= threshold, ]
  }
  
  if (!violin) {
    p1 <- ggplot(df, aes_string(
      x = "name", y = "median", col = "trend",
      alpha = "prob_diff0"
    )) +
      geom_point(size = 3, position = position_dodge(0.3)) +
      geom_errorbar(aes_string(ymin = "lower", ymax = "upper"),
                    position = position_dodge(0.3), width = 0
      ) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() +
      xlab("Time Series") +
      ylab("Loading")+
      guides(fill="none", alpha='none')+
      scale_x_discrete(labels = function(x) str_wrap(x, width = wid))+
      theme(legend.position="none")+
      theme_bw()
  }
  
  if (violin) {
    p1 <- ggplot(v, aes_string(
      x = "name", y = "loading", fill = "trend",
      alpha = "prob_diff0"
    )) +
      geom_violin(color = NA) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() +
      xlab("Time Series") +
      ylab("Loading")+
      theme(axis.text=element_text(size=ax),
            axis.title=element_text(size=ti,face="bold"))+
      guides(fill="none", alpha='none')+
      scale_x_discrete(labels = function(x) str_wrap(x, width = wid))+
      theme_bw()
  }
  
  if (facet) {
    p1 <- p1 + facet_wrap(~trend, scales = "free_x")
  }
  
  p1
}


#dat<-dat%>%select(c(year,n4,n3,n1))


#### DFA FUll ####
y1calcofi <- 1994
y2calcofi<- 2021
n1 <- names(dat)[2:26]
dat.calcofi<-dat%>%select(c(year,n1)) #just setting an order so we know how to set the variance index for each survey
remelt = melt(dat.calcofi,id.vars = "year")
names(remelt)<-c("year","variable","value")
Y <- dcast(remelt, variable ~ year)
names = Y$variable
Y = as.matrix(Y[,-which(names(Y) == "variable")])

n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = c(FALSE),
                       var_index = c("survey"), num_trends = 1,
                       elpd_loo = TRUE, se_elpd_loo=NA)

varIndx = c(rep(1,length(n1)))
fit.mod.calcofi1 = fit_dfa(y = Y,
                          num_trends = 1,
                          iter=n_iter,
                          varIndx = NULL,
                          chains=n_chains, estimate_nu=model_df$est_nu[1],
                          estimate_trend_ma = model_df$estimate_trend_ma[1],
                          estimate_trend_ar = model_df$estimate_trend_ar[1],
                          estimate_process_sigma = FALSE,
                          seed=123)

fit.mod.calcofi2 = fit_dfa(y = Y,
                           num_trends = 2,
                           iter=n_iter,
                           varIndx = NULL,
                           chains=n_chains, estimate_nu=model_df$est_nu[1],
                           estimate_trend_ma = model_df$estimate_trend_ma[1],
                           estimate_trend_ar = model_df$estimate_trend_ar[1],
                           estimate_process_sigma = FALSE,
                           seed=123)


fit.mod.calcofi3 = fit_dfa(y = Y,
                           num_trends = 3,
                           iter=n_iter,
                           varIndx = NULL,
                           chains=n_chains, estimate_nu=model_df$est_nu[1],
                           estimate_trend_ma = model_df$estimate_trend_ma[1],
                           estimate_trend_ar = model_df$estimate_trend_ar[1],
                           estimate_process_sigma = FALSE,
                           seed=123)

is_converged(fit.mod.calcofi1)
is_converged(fit.mod.calcofi2)
is_converged(fit.mod.calcofi3)

bayesdfa::loo(fit.mod.calcofi1)
bayesdfa::loo(fit.mod.calcofi2)
bayesdfa::loo(fit.mod.calcofi3)

fit.mod.calcofi<-fit.mod.calcofi1
namescalcofi<-data.frame(names)
pars = rstan::extract(fit.mod.calcofi$model)
write_rds(fit.mod.calcofi,'data-yellowtail/fit.mod.calcofi.rds')


#### DFA REDUCED ####
y1calcofi <- 1994
y2calcofi<- 2021
n1 <- names(dat%>%select("DDben","Tcop","MLDpart","MLDlarv","MLDpjuv","ONIpre","Tpart","ONIpre","ONIlarv","ONIpjuv","CutiSTIpjuv","CutiTUMIpjuv"))
dat.calcofi<-dat%>%select(c(year,n1)) #just setting an order so we know how to set the variance index for each survey
remelt = melt(dat.calcofi,id.vars = "year")
names(remelt)<-c("year","variable","value")
Y <- dcast(remelt, variable ~ year)
names = Y$variable
Y = as.matrix(Y[,-which(names(Y) == "variable")])

n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear
var <- c("glorys","glorys", "glorys","glorys","glorys","ONI","glorys","ONI", "ONI", "roms", "roms")
model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = FALSE,
                        num_trends = 1,var_index = var,
                       elpd_loo = TRUE, se_elpd_loo=NA)

varIndx = c(rep(1,length(n1)))
varIndx2 <- c("glorys","glorys", "glorys","glorys","glorys","ONI","glorys","ONI", "ONI", "roms", "roms")
varIndx2 <- c(1,1, 1,1,1,2,1,2, 2, 3,3)

fit.mod.calcofi1 = fit_dfa(y = Y,
                          num_trends = 1,
                          iter=n_iter,
                          varIndx = varIndx,
                          chains=n_chains, estimate_nu=model_df$est_nu[1],
                          estimate_trend_ma = model_df$estimate_trend_ma[1],
                          estimate_trend_ar = model_df$estimate_trend_ar[1],
                          seed=123)

fit.mod.calcofi2 = fit_dfa(y = Y,
                           num_trends = 2,
                           iter=n_iter,
                           varIndx =varIndx,
                           chains=n_chains, estimate_nu=model_df$est_nu[1],
                           estimate_trend_ma = model_df$estimate_trend_ma[1],
                           estimate_trend_ar = model_df$estimate_trend_ar[1],
                           seed=123)


fit.mod.calcofi3 = fit_dfa(y = Y,
                           num_trends = 3,
                           iter=n_iter,
                           varIndx = varIndx,
                           chains=n_chains, estimate_nu=model_df$est_nu[1],
                           estimate_trend_ma = model_df$estimate_trend_ma[1],
                           estimate_trend_ar = model_df$estimate_trend_ar[1],
                           seed=123)

is_converged(fit.mod.calcofi1)
is_converged(fit.mod.calcofi2)
is_converged(fit.mod.calcofi3)

bayesdfa::loo(fit.mod.calcofi1)
bayesdfa::loo(fit.mod.calcofi2)
bayesdfa::loo(fit.mod.calcofi3)

fit.mod.calcofi<-fit.mod.calcofi3
namescalcofi<-data.frame(names)
pars = rstan::extract(fit.mod.calcofi$model)
#write_rds(fit.mod.calcofi,'data-yellowtail/fit.mod.calcofi.rds')

r.calcofi <- rotate_trends(fit.mod.calcofi)
p.calcofi <- plot_trends2(r.calcofi,years =dat.calcofi$year)
p.calcofi
l.calcofi <- plot_loadings2(r.calcofi,names=namescalcofi$names)
l.calcofi
is_converged(fit.mod.calcofi)
summary(fit.mod.calcofi)
arranged <- ggarrange(p.calcofi, l.calcofi,ncol = 2, nrow = 1,
                      heights=c(1,1.25,1.5))
arranged

png(filename="figures-yellowtail/DFATrend.png")
plot(arranged)
dev.off()



quadratic_terms = c('LSTlarv','ONIpre','ONIlarv', 'CutiSTIpjuv', 'bakun_sti')

form_dredge = make_dredge_equation(envir_data = envir_data, quadratic_vars = quadratic_terms)
form_dredge


sf_dredge_fun <- function(last_year){
  print(last_year)
 the_data <- data_1%>% filter(year <= last_year) %>%
   dplyr::select(Y_rec, DDpre, DDegg, DDlarv, DDpjuv, DDben, Tcop, Tpart, 
    MLDpart, MLDlarv, MLDpjuv, CSTlarv, CSTpjuv, LSTlarv, 
    LSTpjuv, HCIlarv, HCIpjuv, ONIpre, ONIlarv, ONIpjuv, 
    PDOlarv, PDOpjuv, CutiSTIpjuv, CutiTUMIpjuv, BeutiTUMIpjuv, 
    BeutiSTIpjuv)
 
  the_lm <-  lm(Y_rec ~ DDpre + DDegg + DDlarv + DDpjuv + DDben + Tcop + Tpart + 
    MLDpart + MLDlarv + MLDpjuv + CSTlarv + CSTpjuv + LSTlarv + 
    LSTpjuv + HCIlarv + HCIpjuv + ONIpre + ONIlarv + ONIpjuv + 
    PDOlarv + PDOpjuv + CutiSTIpjuv + CutiTUMIpjuv + BeutiTUMIpjuv + 
    BeutiSTIpjuv + I(LSTlarv^2) + I(ONIpre^2) + I(ONIlarv^2) + 
    I(CutiSTIpjuv^2), data=the_data, na.action = "na.fail")
  
  the_lm_d <- dredge(the_lm, rank = "AICc", m.lim = c(0, 5), 
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
                    dc(CutiSTIpjuv, I(CutiSTIpjuv^2)))
  aic_tibble <- as_tibble(the_lm_d)
  aic_tibble[["last_year"]] <- last_year

  important_variables <- list(sw(the_lm_d), last_year)
 
  aic_and_important_variables <- list(aic_tibble, important_variables)
 
  return(aic_and_important_variables)
}


sf_last_return_year_tibble <- tibble(last_year = 2004:2014) #minimum of 10 years of data

sf_dredge_over_time <- sf_last_return_year_tibble %>% pmap(sf_dredge_fun)

#Indexing: [[A]][[B]][[C]]. A = Final return year (index, not number). B = The dredge table [1] or importance table [2]. C = Final return year (number)

sf_dredge_out <- NULL
for(i in 1:length(sf_last_return_year_tibble$last_year)){
  new_row <- as_tibble(t(as.data.frame(sf_dredge_over_time[[i]][[2]][[1]])))
  new_row$last_year <- sf_dredge_over_time[[i]][[2]][[2]]
  sf_dredge_out <- bind_rows(sf_dredge_out, new_row)
}

sf_dredge_out_long_0 <- sf_dredge_out %>% pivot_longer(!last_year, names_to = "variable", values_to = "value")# %>% left_join(sac_nice_indicators)
sf_dredge_out_long_0$indicator <- str_sub(sf_dredge_out_long_0$variable, 1, -3)

sf_dredge_out_long_0_FULL <- sf_dredge_out_long_0
write_rds(sf_dredge_out_long_01, "timevarying_dredge_singlecov.rds")
write_rds(sf_dredge_out_long_0, "timevarying_dredge.rds")
p <- ggplot(sf_dredge_out_long_0%>%filter(last_year <= 2014), aes(as.factor(last_year), variable, fill= value)) + 
  geom_tile()
p 

p <- ggplot(sf_dredge_out_long_01%>%filter(last_year <= 2014), aes(as.factor(last_year), variable, fill= value)) + 
  geom_tile()
p 

sf_dredge_out_long_0%>%filter(last_year==2014&value>0.15)

the_lm <-  lm(Y_rec ~  LSTlarv + ONIpjuv + ONIlarv + 
    ONIpre  + CutiSTIpjuv + MLDpjuv + DDben, data=data_1)
vif(the_lm)

the_lm <-  lm(Y_rec ~  LSTlarv + ONIpjuv  +
    ONIpre  + CutiSTIpjuv + MLDpjuv + DDben+ CutiSTIpjuv, data=data_1)
vif(the_lm)

lm_sat <-data_1%>%select(LSTlarv, ONIpjuv,ONIlarv,ONIpre,CutiSTIpjuv, MLDpjuv, DDben)
corrplotsat<-cor(lm_sat)
corrplot.mixed(corrplotsat)

lm_sat <-data_1%>%select(LSTlarv, ONIpjuv,ONIpre,CutiSTIpjuv, MLDpjuv, DDben)
corrplotsat<-cor(lm_sat)
corrplot.mixed(corrplotsat)