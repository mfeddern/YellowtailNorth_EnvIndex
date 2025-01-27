# various functions
# source("~/GitHub/Environmental-index-code/Functions-for-envir-index.r")

# build model equation for dredge ##############################################

make_dredge_equation <- function(envir_data, quadratic_vars=NULL, Y = "Y_rec"){
  # rm(base_formula)
  for(j in 1:ncol(envir_data)){# build equation
    x1 = colnames(envir_data)[j]
    if(j==ncol(envir_data)){x2 = x1}else{x2 = paste0(x1," + ")}
    if(j==1){base_formula = x2}else{base_formula = paste0(base_formula,x2)}
    }
  if(!is.null(quadratic_vars)){
    qv = paste0("I(",quadratic_vars,"^2)")
    for(q in 1:length(qv)){
      if(q==length(qv)){q2 = qv[q]}else{q2 = paste0(qv[q]," + ")}
      if(q==1){qv_form = q2}else{qv_form = paste0(qv_form,q2)}
    }
    base_formula = paste0(base_formula," + ", qv_form)
  }
  formula_out = as.formula(paste0(Y," ~ ", base_formula))
  return(formula_out)
}

# find best-fit model from mtable, dredge output ###############################
find_best_fit <- function(mtable, fewest_params = TRUE){
  bfit = subset(mtable,delta<2.0)
  # fewest parameters
  if(fewest_params==TRUE){
    min_df = min(bfit$df)
    bfit = subset(bfit, df == min_df)}
  bfit
  # lowest aic of the remaining models
  bfit = bfit[1,]
  return(bfit)
}  

# re-write best fit model equation for later use ###############################

bf_mod_equation <- function(best_fit, params_suffix = NA){
  params = colnames(coef(best_fit))
  params = params[-1]
  # save out parameters for later
  # model terms for for partial residual plots
  
  x = grep("I\\(", params)
  if(length(x) >0){par_cols = params[-x]}else{par_cols = params}
  
  if(is.na(params_suffix)){
    saveRDS(params, file=paste0(results_dir,'parms_x.rds'))
    saveRDS(par_cols, file=paste0(results_dir,'parms.rds'))}else{
      saveRDS(params, file=paste0(results_dir,'parms_x_',params_suffix,'.rds'))
      saveRDS(par_cols, file=paste0(results_dir,'parms_',params_suffix,'.rds'))  
    }
  ###########
  npar = length(params)
  # make model statement
  if(npar==1){
    form =  paste0("Y_rec ~ ",params[1])}else{
      form = paste0("Y_rec ~ ",params[1])
      for(k in 2:npar){
        form = paste0(form,"+",params[k])
      }
    }
  if(str_detect(form,"NA+")==TRUE){form = "Y_rec ~ 1"}
  
  bf_mod = as.formula(form)
  return(bf_mod)
}

# get model predictions and put in data.frame ##################################

predict_best_fit <- function(old_data, new_data, bf_mod){
  best_fit = lm( bf_mod, data=old_data)
  p1 = predict(best_fit, se.fit = T, newdata = new_data)
  df_pred = data.frame(cbind(new_data, fit = p1$fit, se = p1$se.fit))
  df_index = df_pred %>% dplyr::select(year, fit, se, Y_rec)
  return(df_index)  
}
