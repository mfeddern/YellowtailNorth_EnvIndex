
#########################################################################
### Run the combined SMURF and ODFW oceanography data to get an index of abundance
### Yellowtail rockfish assessment 2025
### Original script - Melissa Monk (SWFSC), edited by Alison Whitman (ODFW)
#########################################################################
# updated by A. Whitman 01/14/2025


rm(list = ls(all = TRUE))
graphics.off()

library(sdmTMB)
library(tmbstan)
library(ggeffects)
library(MuMIn)
library(here)
library(glue)
library(tidyr)
library(dplyr)
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(bayesplot)
library(grid)
library(devtools)
library(ggeffects)
library(tidybayes)
library(gridExtra)
library(fitdistrplus)

#species and area identifiers - eventually put in function
pacfinSpecies <- "YTRF"
speciesName <- "Yellowtail Rockfish"
modelArea = "oregon"
indexName <-  "SMURF"
modelName <- "full"

# loading helper functions 

#dir<-file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS/Raw Index Scripts")
#setwd(dir)
#list.files()
#source("helper_functions.R")
#source("diagnostics.R")
#source("do_diagnostics.R")
#source("format_hkl_data.R")
#source("format_index.R")
#source("get_index.R")
#source("match.f.R")
#source("plot_betas.R")
#source("plot_index.R")
#source("refactor.R")

# load data
#dir <- file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS")
#setwd(dir)

#load("data_for_GLM.RData")
dat<-read.csv("combined_settlement_ocean.csv")

# subset to species of interest  - unnecessary here
#dat<-dat[dat$Common_Name==speciesName,]

# set dir for full model
#dir <- file.path("C:/Users/daubleal/OneDrive - Oregon/Desktop/2025 Assesssment Cycle/Index_SMURFS/yellowtail_oregon_SMURF_full")
#setwd(dir)

# explore the data 

summary(dat)
names(dat)
#dat$CPUE<-dat$Counts/dat$Total_Drift_Effort
dat$year
## note using their calculated settlement rate rather calculating my own ##

dat$CPUE<-dat$SFLA_settlement_rate

# add some categorical factors to explore 
dat$region<-ifelse(dat$Site %in% c("HH","RR"),"South","North")
dat$treatment<-ifelse(dat$Site %in% c("HH","CF"),"Comparison","Reserve")
dat$Date<-as.POSIXct(dat$Date, format = "%Y-%m-%d")
dat$month<-as.numeric(format(dat$Date,"%m"))
table(dat$month)
summary(dat$julian)
dat$season<-ifelse(dat$julian<196,"early","late") # used July 15 as the divider 
table(dat$season)

#Look at the data
pos <- dat[dat$SFLA_count>0,]

# 255 positives out of 1376 = 18% not bad! 

with(pos, table(Site)) # all pretty even across sites
with(pos, table(treatment)) 
with(pos, table(region)) 
with(pos, table(season))

# start here

ggplot(pos, aes(Site, CPUE)) + # 
  geom_boxplot()

ggplot(pos, aes(treatment, CPUE)) + 
  geom_boxplot()

ggplot(pos, aes(region, CPUE)) + # oh wow, higher in south...
  geom_boxplot()

ggplot(pos, aes(season, CPUE)) + # more in the early season 
  geom_boxplot()

ggplot(pos, aes(temp_c, CPUE)) + # primarily in the 8 to 10 deg range, could bin these
  geom_point(alpha = 0.5)

with(pos, table(year, month)) # thin on samples in August and Sept
with(pos, table(month))
with(dat, table(month))

# fix or add anything 
dat <- dat %>%
  #rename(Year = year) %>%
  dplyr::rename(Effort = Sampling.Interval) %>%
  mutate(logEffort = log(Effort)) %>% 
  #create temp bins for drill
  mutate(temp_bin = cut(temp_c, breaks=c(6,7,8,9,10,11,12)))
  # remove samples - not going to remove anything at this time
  #filter(Ave_Depth > 10) %>%
  #filter(Month != 7)

# pos <- subset(dat, Counts>0)
# 
# with(pos,table(depth_bin))
# 
# ggplot(pos, aes(depth_bin, CPUE)) + 
#   geom_point(alpha = .5)
# 
# ggplot(dat, aes(rock_bin, CPUE)) + 
#   geom_point(alpha = .5)
# 
# with(pos,table(season))
# 
# ggplot(pos, aes(season, CPUE)) + 
#   geom_point(alpha = .5)
# 
# # subset to MR or CA only 
# dat <- dat %>%
#   filter(Treatment == "CA")

# reset dir if not running full model
#dir <- file.path("C:/Users/daubleal/Desktop/2023 Assessment Cycle/Index_MARINE RESERVES/black_oregon_marres_hnl_CAonly")
#setwd(dir)

# define my covars
covars <- c("year", "region", "treatment", "season") # full model 
#covars <- c("year", "season" ,"Site","depth_bin")

#Ensure columns named appropriately and covariates are factors
dat <- dat %>%
  #filter(Treatment == "MR") %>% # testing one with just the MR sites
  mutate_at(covars, as.factor) 

# model selection 

model.full <- MASS::glm.nb(
  SFLA_count ~ year + region + treatment + season + offset(logEffort),
  data = dat,
  na.action = "na.fail") # I'm getting errors using "na.fail" here (I also tried changing it but then it wouldn't work in the dredge function below) 
summary(model.full)
anova(model.full)
#use ggpredict to get an estimate of the logEffort for sdmTMB predictions
#MuMIn will fit all models and then rank them by AICc
model.suite <- MuMIn::dredge(model.full,
                             rank = "AICc", 
                             fixed= c("offset(logEffort)", "year"))

#Create model selection dataframe for the document
Model_selection <- as.data.frame(model.suite) %>%
  dplyr::select(-weight)
Model_selection

# sdmTMB model 

#set the grid
grid <- expand.grid(
  year = unique(dat$year),
  #Month = levels(dat$Month)[1],
  region = levels(dat$region)[1],
  season = levels(dat$season)[1],
  #Site = levels(dat$Site)[1],
  treatment = levels(dat$treatment)[1]
#  temp_bin = levels(dat$temp_bin)[1]
  #rock_bin = levels(dat$rock_bin)[1]
)

fit.nb <- sdmTMB(
  SFLA_count ~  as.factor(year)+region+season+treatment,
  data = dat,
  offset = dat$logEffort,
  time = "year",
  spatial="off",
  spatiotemporal = "off",
  family = nbinom2(link = "log"),
  control = sdmTMBcontrol(newton_loops = 1)) #documentation states sometimes aids convergence?

#}
fit.nb$sdreport
#Get diagnostics and index for SS
do_diagnostics(
  dir = file.path(dir), 
  fit = fit.nb,
  plot_resids = F)

calc_index(
  dir = file.path(dir), 
  fit = fit.nb,
  grid = grid)

#-------------------------------------------------------------------------------
#Format data filtering table and the model selection table for document

# will need to modify my filter script to get the dataFilters to work
#View(dataFilters)

dataFilters <- dataFilters %>%
  rowwise() %>%
  filter(!all(is.na(across((everything()))))) %>%
  ungroup() %>%
  rename(`Positive Samples` = Positive_Samples)
dataFilters <- data.frame(lapply(dataFilters, as.character), stringsasFactors = FALSE)
#write.csv(dataFilters, file = file.path(dir, "data_filters.csv"), row.names = FALSE)

#View(Model_selection)
#format table for the document
out <- Model_selection %>%
  dplyr::select(-`(Intercept)`) %>%
  mutate_at(vars(covars,"year","offset(logEffort)"), as.character) %>%
  mutate(across(c("logLik","AICc","delta"), round, 1)) %>%
  # replace_na(list(district = "Excluded",                      # fix these later
  #                 targetSpecies = "Excluded", month = "Excluded")) %>% # fix later
  mutate_at(c(covars,"year","offset(logEffort)"), 
            funs(stringr::str_replace(.,"\\+","Included"))) %>%
  rename(`Effort offset` = `offset(logEffort)`, 
         `log-likelihood` = logLik) %>%
  rename_with(stringr::str_to_title,-AICc)
View(out)
#write.csv(out, file = file.path(dir,  "model_selection.csv"), row.names = FALSE)

#summary of trips and  percent pos per year
summaries <- dat %>%
  group_by(year) %>%
  summarise(tripsWithTarget = sum(Counts>0),
            tripsWOTarget = sum(Counts==0)) %>%
  mutate(totalTrips = tripsWithTarget+tripsWOTarget,
         percentpos = round(tripsWithTarget/(tripsWithTarget+tripsWOTarget),2)) 
View(summaries)
#write.csv(summaries,
 #       file.path(dir,  "percent_pos.csv"),
  #     row.names=FALSE)



# STAR comparison 

# compare the combined, MR only and CA only to each other 

# reset wd
dir <- file.path("C:/Users/daubleal/Desktop/2023 Assessment Cycle/Index_MARINE RESERVES/STAR_treatment comparison")
setwd(dir)

# pull in all three
ref_model<-read.csv("index_forSS_refmodel.csv")
CAonly<-read.csv("index_forSS_CAonly.csv")
MRonly<-read.csv("index_forSS_MRonly.csv")

# standardize 
ref_model$std_Ref<-ref_model$obs/mean(ref_model$obs)
CAonly$std_CA<-CAonly$obs/mean(CAonly$obs)
MRonly$std_MR<-MRonly$obs/mean(MRonly$obs)

# combine and plot
df1 <- data.frame(Year = ref_model$year, Variable = ref_model$std_Ref)
df2 <- data.frame(Year = CAonly$year, Variable = CAonly$std_CA)
df3 <- data.frame(Year = MRonly$year, Variable = MRonly$std_MR)
df4 <- df1 %>%  mutate(Index = 'Combined') %>%
  bind_rows(df2 %>%
              mutate(Index = 'CA only')) %>%
  bind_rows(df3 %>% 
              mutate(Index = 'MR only'))
df4$Year<-as.factor(df4$Year)

ggplot(df4,aes(y = Variable,x = Year,color = Index)) + 
  geom_line(aes(group = Index))+
  geom_point(size = 1.5) +
  labs(y = "Standardized Index Value")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x = element_blank())
#ggsave("Marine Reserves_treatment index comparison.png",width = 10,height = 5)





