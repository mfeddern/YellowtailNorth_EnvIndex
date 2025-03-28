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
library(ggpubr)
library(car)
library(forcats)


univariate_ns<-data.frame(read.csv("results-yellowtail/univariateNS_Results.csv"))
univariate_gam<-data.frame(read.csv("results-yellowtail/univariateGAM_Results.csv"))
univariate_lm<-data.frame(read.csv("results-yellowtail/univariateLM_Results.csv"))
mult_gam<-data.frame(read.csv("results-yellowtail/MultivariateGAMresults.csv"))
marginal_gam<-data.frame(read.csv("results-yellowtail/MarginalImprovementGAM.csv"))
marginal_lm<-data.frame(read.csv("results-yellowtail/MarginalImprovementLM.csv"))
marginal_NS<-data.frame(read.csv("results-yellowtail/MarginalImprovementNS.csv"))
df = data.frame(read.csv("data-yellowtail/DATA_Combined_glorys_yellowtail.csv"))
dfa = data.frame(read.csv("data-yellowtail/dfa_trend.csv"))
data_years = 1994:2014
dat = df %>% 
  dplyr::select(!any_of( c('ZOOpjuv', 'ZOOben')))  %>% 
  filter(year %in% data_years)%>%
  left_join(dfa)
envir_data = dat %>%  # drop terms not in model statement
  dplyr::select(!any_of(c('sd','year','Y_rec','ZOOpjuv','ZOOben')))%>%  # drop terms not in model statement
    dplyr::select(!any_of(c('sd','Y_rec','ZOOpjuv','ZOOben','LUSI','ONIlarv')))%>% 
    mutate_all(~ scale(.)) # drop terms not in model statement
data_types<-data.frame(cov=c(colnames(envir_data),"CutiSTIpjuv +I(CutiSTIpjuv^2)"), type=c("Temperature", "Temperature", "Temperature","Temperature","Temperature","Temperature","Temperature",
                  "Transport", "Transport", "Transport", "Transport", "Transport", "Transport", "Transport",
                  "Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Upwelling",
                  "Upwelling", "Upwelling", "Upwelling","Ran", "Temperature","Upwelling"))



head(envir_data)
dim(envir_data)
ggplot(data=dat, aes(x=year, y=Y_rec))+
  geom_point()+
  geom_line()+
  xlim(c(1990,2016))+
  ylab("Recruitment Deviations")+
  xlab("Year")+
  ggtitle("Northern Yellowtail Rockfish")+
  geom_ribbon(aes(ymin=Y_rec-sd,ymax=Y_rec+sd), col = 'grey', alpha=0.2)+
 theme_classic()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/RecDev_TimeSeries.png", height = 6, width = 8)


envir_long<-envir_data%>%
  select(-X)%>%
  cbind(year=dat$year)%>%
  pivot_longer(c(-year),names_to = "variable",values_to = "values")
  
envir_long<-left_join(envir_long,data.frame(variable=unique(envir_long$variable),
           type=c("Temperature", "Temperature", "Temperature","Temperature","Temperature","Temperature","Temperature",
                  "Transport", "Transport", "Transport", "Transport", "Transport", "Transport", "Transport",
                  "Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Upwelling",
                  "Upwelling", "Upwelling", "Upwelling", "Temperature")))

ggplot(data=envir_long, aes(x=year, y=values, group=type))+
  facet_wrap(~variable)+
  geom_point(aes(col=type))+
  geom_line(aes(col=type))+
  xlim(c(1990,2016))+
  ylab("Standardized Anomalies")+
  xlab("Year")+
  ggtitle("Northern Yellowtail Rockfish")+
 theme_classic()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/Env_TimeSeries.png", height = 8, width = 8)

temp<-ggplot(data=envir_long%>%filter(type=="Temperature"), aes(x=year, y=values, group=variable))+
  #facet_wrap(~type)+
  geom_point(aes(col=variable))+
  geom_line(aes(col=variable))+
  xlim(c(1990,2016))+
  ylab("")+
  xlab("")+
  ggtitle("Temperature Environmental Indices")+
 theme_classic()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))

upwelling<-ggplot(data=envir_long%>%filter(type=="Upwelling"), aes(x=year, y=values, group=variable))+
  #facet_wrap(~type)+
  geom_point(aes(col=variable))+
  geom_line(aes(col=variable))+
  xlim(c(1990,2016))+
  ylab("Standardized Anomalies")+
  xlab("")+
  ggtitle("Upwelling Environmental Indices")+
 theme_classic()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))

transport<-ggplot(data=envir_long%>%filter(type=="Transport"), aes(x=year, y=values, group=variable))+
  #facet_wrap(~type)+
  geom_point(aes(col=variable))+
  geom_line(aes(col=variable))+
  xlim(c(1990,2016))+
  ylab("")+
  xlab("Year")+
  ggtitle("Transport Environmental Indices")+
 theme_classic()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))

ggarrange(temp, upwelling, transport, nrow=3)
ggsave("figures-yellowtail/EnvironmnetalIndices.png", height = 10, width = 8)

marginal<-rbind(marginal_gam,marginal_lm,marginal_NS)
ggplot(data=marginal, aes(x=total_rmse, y=cov, group=cv))+
  facet_wrap(~model)+
  geom_point(aes(col=cv),size=2.5)+
  xlim(c(0.5,1.5))+
  ylab("Oceanographic Conditions")+
  xlab("Marginal Mean Improvement RMSE")+
 # ggtitle("Generalized Additive Models")+
 theme_classic()+
   geom_vline(xintercept = 1, 
             color = "grey", linetype = "dashed", size = 1) +
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))

marginal2<-rbind(marginal_gam,marginal_lm,marginal_NS)%>%
  mutate(RMSE_perc=1-total_rmse)%>%
#  filter(RMSE_perc>0)%>%
  left_join(data_types)

ggplot(marginal2%>%
  group_by(cv,model)%>%
  mutate(name = fct_reorder(cov, desc(RMSE_perc))%>%ungroup()), aes(x=cov, y=RMSE_perc,fill=type)) +
    geom_bar(stat="identity",  alpha=.6, width=.4) +
  facet_grid(model~cv)+  
  coord_flip() +
    ylab("Percent Change in RMSE") +
  xlab("Oceanogrpahic Conditions") +
    theme_bw()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/MarginalImprovementv1.png", height = 8, width = 10)

marginal3<-rbind(marginal_gam,marginal_lm)%>%
  mutate(RMSE_perc=1-total_rmse)%>%
  #filter(RMSE_perc>0)%>%
  left_join(data_types)

ggplot(marginal3%>%
  group_by(cv,model)%>%
  mutate(name = fct_reorder(cov, desc(RMSE_perc)))%>%ungroup(), aes(x=cov, y=RMSE_perc,fill=type,group=type)) +
    geom_bar(stat="identity",  alpha=.6, width=.4) +
  facet_grid(model~cv)+  
  coord_flip() +
    ylab("Percent Change in RMSE") +
  xlab("Oceanogrpahic Conditions") +
    theme_bw()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/MarginalImprovementv2.png", height = 12, width = 10)

ggplot(marginal2%>%filter(model=="NS")%>%
  group_by(cv,model)%>%
  mutate(name = fct_reorder(cov, desc(RMSE_perc)))%>%ungroup(), aes(x=cov, y=RMSE_perc,fill=type,group=type)) +
    geom_bar(stat="identity",  alpha=.6, width=.4) +
  facet_grid(model~cv)+  
  coord_flip() +
    ylab("Percent Change in RMSE") +
  xlab("Oceanogrpahic Conditions") +
    theme_bw()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/MarginalImprovementNS.png", height = 12, width = 10)

margtop5<-marginal3%>%
  filter(RMSE_perc>0)%>%
  arrange(cv,model,desc(RMSE_perc))%>%
  group_by(cv,model)%>% 
  slice(1:5)

ggplot(margtop5%>%
  group_by(cv,model)%>%
  mutate(name = fct_reorder(cov, desc(RMSE_perc)))%>%ungroup(), aes(x=cov, y=RMSE_perc,fill=type,group=type)) +
    geom_bar(stat="identity",  alpha=.6, width=.4) +
  facet_grid(model~cv)+  
  coord_flip() +
    ylab("Percent Change in RMSE") +
  xlab("Oceanogrpahic Conditions") +
    theme_bw()+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))

ggsave("figures-yellowtail/MarginalImprovementTop5.png", height = 12, width = 10)
write.csv(margtop5,"results-yellowtail/Top5MarginalImprovement.csv")



ggplot(univariate_lm%>%select(var,RMSE_LFO10_rel, RMSE_LOO_rel,RMSE_LFO5_rel)%>%
  rename(LFO10=RMSE_LFO10_rel, LOO=RMSE_LOO_rel,LFO5=RMSE_LFO5_rel)%>%
  pivot_longer(-var,names_to='Model', values_to = 'value')%>%group_by(Model)%>%mutate(rank=rank(value, ties.method="max"))%>%
    mutate(top5=ifelse(rank>=23,1,0))%>%ungroup(), 
  aes(x=var, y=value,fill=top5)) +
    geom_bar(stat="identity",  alpha=.9, width=.4) +
  facet_grid(~Model,scales = "free_x")+  
  coord_flip() +
  scale_fill_gradient(low ="azure2",high = "darkcyan")+
    ylab("Percent Change in RMSE") +
  xlab("Oceanogrpahic Conditions") +
    theme_bw()+
  geom_hline(yintercept=0,lty=2)+
  ggtitle("Linear Models")+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/UnivariateLmTop5.png", height = 8, width = 10)

ggplot(univariate_gam%>%select(var,RMSE_LFO10_rel, RMSE_LOO_rel,RMSE_LFO5_rel)%>%
  rename(LFO10=RMSE_LFO10_rel, LOO=RMSE_LOO_rel,LFO5=RMSE_LFO5_rel)%>%
  pivot_longer(-var,names_to='Model', values_to = 'value')%>%group_by(Model)%>%mutate(rank=rank(value, ties.method="max"))%>%
    mutate(top5=ifelse(rank>=22,1,0))%>%ungroup(), 
  aes(x=var, y=value,fill=top5)) +
    geom_bar(stat="identity",  alpha=.9, width=.4) +
  facet_grid(~Model,scales = "free_x")+  
  coord_flip() +
    ggtitle("Generalized Additive Models")+
  scale_fill_gradient(low ="azure2",high = "darkcyan")+
    ylab("Percent Change in RMSE") +
  xlab("Oceanogrpahic Conditions") +
    theme_bw()+
  geom_hline(yintercept=0,lty=2)+
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/UnivariateGamTop5.png", height = 8, width = 10)

univariate_full<-univariate_gam%>%bind_rows(univariate_lm)

ggplot(data=univariate_full%>%select(var,model,rsq_LFO10, rsq_LOO,rsq_LFO5)%>%
  rename(LFO10=rsq_LFO10, LOO=rsq_LOO,LFO5=rsq_LFO5)%>%
  pivot_longer(-c(var,model),names_to='cv', values_to = 'value'), 
  aes(x=value, y=var, group=cv))+
  facet_wrap(~model)+
  geom_point(aes(col=cv),size=2.5)+
  xlim(c(-0.25,0.75))+
  ylab("Oceanographic Conditions")+
  xlab("R squared")+
 # ggtitle("Generalized Additive Models")+
 theme_classic()+
   geom_vline(xintercept = 1, 
             color = "grey", linetype = "dashed", size = 1) +
  theme(axis.text = element_text(size = 11),plot.title = element_text(hjust = 0.5))
ggsave("figures-yellowtail/UnivariateRsquared.png", height = 6, width = 10)
