
print(code_dir)

# run short version excluding last 5 years #####################################

# refit leaving out the last five years and project ############################

# out.file = paste0(home_dir,"/results-short")
# dir.create(out.file)
# setwd(out.file)
# source(paste0(code_dir,'/_02_dredge-short.r'))

################################################################################

# includes Retro calculations ##################################################
out.file = test_dir
# dir.create(out.file)
setwd(out.file)

# diagnostics ##################################################################
  source(paste0(code_dir,'/01_model_diagnostics.r'))

# bootstrap recruitments #######################################################
# get estimates if recruitment is random
  print("Running bootstrap-1")
  out.file = paste0(test_dir,"02_bootstrap-1")
  dir.create(out.file)
  setwd(out.file)
  source(paste0(code_dir,'/02_bootstrap-1.r'))  

# bootstrap rows for bias estimates ############################################ 
# get bias estimates
  print("Running boostrap-2-bias")
  out.file = paste0(test_dir,"02-bootstrap-2-bias")
  dir.create(out.file)
  setwd(out.file)
  source(paste0(code_dir,'/02-bootstrap-2-bias.r'))
  
# jackknife best fit model by year #############################################
# impact of individual years on best model fit
  # includes cross validation with rmse for out of sample years
  print("Jackknifing")
  out.file = paste0(test_dir,"/03_jackknife-best-fit")
  dir.create(out.file)
  setwd(out.file)
  source(paste0(code_dir,'/03_jackknife-best-fit.r'))
  
# # resample by year with error ###################################################
#   # need to tinker with this file based on form of equation Y vs log(Y)
#   print("resample with error")
#   out.file = paste0(test_dir,"/04_resample-w-error")
#   dir.create(out.file)
#   setwd(out.file)
#   source(paste0(code_dir,'/04_resample-w-error.R'))
  
# Jackknife - by year and refit #########################################################
  # dredge formula saved earlier when running dredge
  # form_dredge
  print ("Jackknife and dredge")
  out.file = paste0(test_dir,"05_jackknife-dredge-refit")
  dir.create(out.file)
  setwd(out.file)
  source(paste0(code_dir,'/05_jackknife-dredge-refit.r'))
  
# # jackknife - by year - re Dredge w/error  #####################################

  #   # need to tinker with this file to set dredge settings
#   # need to tinker with this file based on form of equation Y vs log(Y)
#   print ("dredge with error")
#   out.file = paste0(test_dir,"/06_dredge-refit-w-error")
#   dir.create(out.file)
#   setwd(out.file)
#   source(paste0(code_dir,'/06_dredge-refit-w-error.r'))
#   
# leave out five year chunks and refit moving forward one year at a time ######
  
  print("Leave out 5-year blocks and redredge")
  out.file = paste0(test_dir,"07_dredge-5block")
  dir.create(out.file)
  setwd(out.file)
  source(paste0(code_dir,'/07_dredge-block5.r'))
  
# redredge year by year over the last five years predicting one year ahead #####
  print("Dredge and predict one year ahead")
  out.file = paste0(test_dir,"08_dredge-next-year")
  dir.create(out.file)
  setwd(out.file)
  source(paste0(code_dir,'/08_dredge-next-year.r'))

# figure and table plotting ####################################################  

# needs attention based on available data
# source(paste0(code_dir,'/Figure_1&2.r'))

print("printing fig 3")
source(paste0(code_dir,'/Figure_3.r'))
print("printing fig 4")
source(paste0(code_dir,'/Figure_4.r'))
print("printing fig 5")
source(paste0(code_dir,'/Figure_5.r'))
print("printing table 3")
source(paste0(code_dir,'/Table_3.r'))

