library(tidyverse)
library(ggplot2)
library(readr) # faster writing
library(data.table) # faster writing of large files
library(lubridate)
library(mgcv)
library(MuMIn)
library(corrplot)

source("Src/_00_yellowtail-header.r")
source("Src/Functions-for-envir-index.r")
# set directories ##############################################################

# bring in data ################################################################
# combined fish and environmental drivers file
df = data.frame(read.csv(paste0(data_dir,"DATA_Combined_glorys_yellowtail.csv"), header = T))

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

unique(uncorr$xvar)
corrPlot<-corrplot.mixed(M)




pdf(file = "figures-yellowtail/CorrPlot.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12) # The height of the plot in inches
corrplot.mixed(M)

dev.off()



library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
library(ggplot2)
library(stringr)
library(ggpubr)
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

for(i in 1:ncol(dat)){
  if(names(dat)[i] %in% ids){
    dat[,i] <- log(dat[,i])
  }
}

#dat<-dat%>%select(c(year,n4,n3,n1))


#### CALCOFI####
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
