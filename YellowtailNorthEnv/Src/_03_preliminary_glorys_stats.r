source("_00_yellowtail-header.r")
library(corrplot)
library(mgcv)

# main_file = getwd()
# data_loc = paste0(main_file,"/data/")
# results_loc = paste0(main_file, "/results/")
# fig_loc = paste0(main_file,'/figures/')

df = fread(paste0(data_loc, "DATA-glorys-4petrale.csv"))

# correlations among variables #################################################
dfx = df %>% dplyr::select( !year )
cor_xy = cor(dfx)

graphics.off()
png( paste0(fig_loc,"glorys-correlations-among-variables.png"),
     units = 'in', res=300, width = 6.5, height=6.5)
corrplot::corrplot(cor_xy,  method='number', type='lower', number.cex = 0.6 )
dev.off()


# time series of each parameter ################################################

dflong = pivot_longer(df, cols = 2:ncol(df), names_to = "driver", values_to = "index")

dflong$driver = factor(dflong$driver, 
                       levels = c("DDpre","Tpre.a","Tpre.b", 
                                  "DDegg1", "LSTegg1","CSTegg1" , "MLDegg",
                                  "DDegg2","LSTegg2" ,"CSTegg2",
                                  "DDlarv","LSTlarv", "CSTlarv",
                                  "DDpjuv","LSTpjuv","CSTpjuv" ,
                                  "DDbjuv.a","LSTbjuv.a", "CSTbjuv.a",
                                  "DDbjuv.b",  "LSTbjuv.b", "CSTbjuv.b"))


graphics.off()
png (paste0(fig_loc, "glorys_time_series.png"), units = 'in',res=300, height=7, width=6)

ggplot( dflong, aes(x=year, y=index)) +
  geom_line() + xlab('Year') + ylab("Index") +
  facet_wrap(facets='driver', ncol=3, scales = 'free_y') + 
  scale_x_continuous(breaks = seq(2000,2020,10), minor_breaks = 1993:2023) +
  theme_bw() + theme(strip.background = element_blank())
dev.off()

# non-linear-relationships ? ###################################################

fish = fread( paste0(data_loc,"recdevs_2023.csv"))
fish = fish[,c("Year", "Value","Parm_StDev")]
colnames(fish) = c('year','recdev','sd')

df = left_join(df,fish)
cn = colnames(dfx)
dfz = data.frame(df %>% filter(year %in% 1993:2018))

fwrite(df,  paste0(data_loc,"DATA_Combined_glorys_petrale.csv"))

for( i in 1:length(cn)){
  Y = dfz[,cn[i]]
  lmx = lm( dfz$recdev ~ Y)
  lmp = lm( dfz$recdev ~ Y + I(Y^2) )
  gmx =  gam(dfz$recdev ~ s(Y) )
  
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

fwrite(A, paste0(fig_loc,"Prelim-non-linear-comparison.csv"))

dflong = left_join(dflong, fish)

graphics.off()
png( paste0(fig_loc, "Recruitment deviations vs glorys drivers.png"),
     units = 'in',res=300, width=6, height=7)

ggplot(dflong, aes(x=index, y=recdev)) + 
  geom_point(size=1) + xlab("GLORYS index") + ylab("Recruitment deviations") +
  facet_wrap(facets='driver', ncol=4, scales = 'free') + 
  theme_bw() + theme(strip.background = element_blank())
dev.off()






