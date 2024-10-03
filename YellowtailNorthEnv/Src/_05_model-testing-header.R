
# directories ##################################################################
setwd("~/GitHub/Env-Index-Devel/")
home_dir = getwd()
code_dir = "~/GitHub/Environmental-index-code"
data_dir = paste0(home_dir,"/data-yellowtail/")
test_dir = paste0(home_dir,"/model_testing-yellowtail/")
results_dir = paste0(home_dir, "/results-yellowtail/")
fig_dir = paste0(home_dir,"/figures-yellowtail/")
fit_dir = paste0(home_dir, "/results-yellowtail/")

# source("C:/Users/tolim/Documents/GitHub/Environmental-index-code/Functions-for-envir-index.r")
# code_dir = "C:/Users/tolim/Documents/GitHub/Environmental-index-code/"

source("~/GitHub/Environmental-index-code/Functions-for-envir-index.r")

# for processed time series for hake
# Set years for testing etc ####################################################

dir.create(test_dir)

data_years = 1994:2014
short_years = 1994:2009 # fit to here. In sample data
forecast_years = 2010:2014 # predict to here. Out of sample data

# bring in data ################################################################
# names in all other files is df0. Do not rename.

df0 = data.frame(read.csv(paste0(data_dir,"DATA_Combined_glorys_yellowtail.csv"), header = T))
df0 = df0[df0$year %in% data_years,]

# set for later plotting.

# model info from dredge file ################################################## 
# rec = "ln_dev" # column name for recruitment in the file
fit_dredge=readRDS(paste0(fit_dir, "fit_dredge.rds"))
form_dredge = readRDS("formula_for_dredge.rds")
form_best = readRDS(paste0(results_dir, "Best_fit_model.rds"))
parms   = readRDS(paste0(results_dir,"parms.rds")) # just column headers for included predictors
parms_x = readRDS(paste0(results_dir,"parms_x.rds")) # all terms for partial resid plots
# just label y-axes
Recruitment_Type = "ln(rec devs)"

# run diagnostics

source(paste0(code_dir, "/_04_model-testing-wrapper.r"))

