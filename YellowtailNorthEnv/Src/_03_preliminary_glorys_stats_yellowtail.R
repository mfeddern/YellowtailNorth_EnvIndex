library(tidyverse)
library(ggplot2)
library(readr)
library(data.table) 
library(lubridate)
library(corrplot)
library(mgcv)

source("_00_yellowtail-header.r")
source("~/GitHub/glorys-download/glorys_helper_functions.r")

df = fread(paste0(data_dir, "glorys-data-annual-yellowtail.csv"))
head(df)

# correlations among variables #################################################
dfx = df %>% dplyr::select( DDpre:LSTpjuv )
dfx = na.omit(dfx)
cor_xy = cor(dfx)

graphics.off()
png( paste0(fig_dir,"glorys-correlations-among-variables.png"),
     units = 'in', res=300, width = 6.5, height=6.5)
corrplot::corrplot(cor_xy,  method='number', type='lower', number.cex = 0.6 )
dev.off()

dredge_text = correlated_variables_dredge(dfx) 
dredge_text
saveRDS(dredge_text,file = paste0(results_dir, 'dredge_correlation_text.rds'))


# shortened for missing data ###################################################

dfx = df %>% dplyr::select( !year )
dfx = na.omit(dfx)
cor_xy = cor(dfx)

graphics.off()
png( paste0(fig_dir,"glorys-correlations-among-variables-short.png"),
     units = 'in', res=300, width = 6.5, height=6.5)
corrplot::corrplot(cor_xy,  method='number', type='lower', number.cex = 0.6 )
dev.off()

dredge_text_short = correlated_variables_dredge(dfx) 
dredge_text_short
saveRDS(dredge_text_short,file = paste0(results_dir, 'dredge_correlation_text_short.rds'))

# time series of each parameter ################################################
dflong = pivot_longer(df, cols = 2:ncol(df), names_to = "driver", values_to = "index")

# use to order plots
dflong$driver = factor(dflong$driver, 
                      levels = c("DDpre", "DDegg", "DDlarv", 
                                 "DDpjuv", "DDben", "CSTlarv", 
                                 "CSTpjuv", "LSTlarv", "LSTpjuv",
                                 "MLDlarv", "MLDpart", "MLDpjuv",
                                 "Tcop", "Tpart", "hci1_larv",
                                 "hci1_pjuv", "hci2_larv", "hci2_pjuv",
                                 "oni_pre", "oni_larv", "oni_pjuv",
                                 "pdo_larv", "pdo_pjuv", "lusi_annual"))
                                  
ggplot( dflong, aes(x=year, y=index)) +
  geom_line() + xlab('Year') + ylab("Index") +
  facet_wrap(facets='driver', ncol=3, scales = 'free_y') + 
  scale_x_continuous(breaks = seq(2000,2020,10), minor_breaks = 1993:2023) +
  theme_bw() + theme(strip.background = element_blank())

ggsave(paste0(fig_dir, "glorys_time_series.png"), height=7, width=6)


# non-linear-relationships ###################################################
fish = readRDS(paste0(data_dir,'yellowtail_recruit.rds'))
fish = fish[,c("Yr", "Value","Parm_StDev")]
colnames(fish) = c('year','Y_rec','sd')

df = left_join(df,fish)
cn = colnames(dfx)
dfz = data.frame(df %>% filter(year %in% 1993:2018))

fwrite(df,  paste0(data_dir,"DATA_Combined_glorys_yellowtail.csv"))

for( i in 1:length(cn)){
  Y = dfz[,cn[i]]
  lmx = lm( dfz$Y_rec ~ Y)
  lmp = lm( dfz$Y_rec ~ Y + I(Y^2) )
  gmx =  gam(dfz$Y_rec ~ s(Y) )
  
  a = data.frame('Driver' = cn[i],
                 'Linear.AIC' = AIC(lmx),
                 'Polynomial.AIC' = AIC(lmp),
                 'GAM.AIC' = AIC(gmx))
  
  mina = min(a[,2:4])
  a$Best = ifelse(mina == a[,2], 'Linear', ifelse(mina==a[,3], 'Poly', 'GAM') )
  a$Linear.vs.Poly = ifelse( a$Linear.AIC < a$Polynomial.AIC, 'Linear', 'Poly')
  if(i==1){A = a}else{A = rbind(A,a)  }
}

A$LvPaic = round(A$Linear.AIC - A$Polynomial.AIC,3)
A

fwrite(A, paste0(fig_dir,"Prelim-non-linear-comparison.csv"))

dfx = left_join(dflong,fish)


ggplot(dfx, aes(x=index, y=Y_rec)) +
  geom_point(size=1) + xlab("GLORYS index") + ylab("Recruitment deviations") +
  facet_wrap(facets='driver', ncol=3, scales = 'free') +
  theme_bw() + theme(strip.background = element_blank())
ggsave( paste0(fig_dir, "Correlations recuitment vs driver.png"), height=7, width=6)

