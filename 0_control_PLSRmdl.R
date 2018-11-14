# ================================================================================ 
# Control file for PLSR HUC models processing
# PLSR code for processing CFSv2 variables for predicting 
#
# Created by S. Baker, June 2018
# ================================================================================ 

## === libraries
suppressMessages(library(lubridate))
suppressMessages(library(ncdf4))
suppressMessages(library(abind))
suppressMessages(library(OceanView))
suppressMessages(library(pls))
suppressMessages(library(dplyr))
# suppressMessages(detach(package:plyr)) # ensure group_by works


### === INPUT prediction: hru, time, variable, lead (if not input)
if (!exists('preds'))     { preds    = c(8, 11)   } # predictors 
if (!exists('pred_var'))  { pred_var = 'prate'    } # variable to predict 'tmp2m' or 'prate' 
if (!exists('wk'))        { wk       = '2_3'      } # week to analyze 
if (!exists('pred_mon'))  { pred_mon = 7          } # month to predict in CV
if (!exists('hru'))       { hru      = '1401'     } # '1802' '1304''1401' '1709'

## predicted month input
mon_nm           = month.abb[pred_mon]   # months for model data (labeling purposes)

### === Domain info
cut_domain       = T           # T/F domain cut?
cut_num          = 1           # label of version of cut spatial domain for variables
res_domain       = 2           # domain resolution x by x degree
l                = 1           # days of lag. If 0, no lag

# ------------------------------------------------------------------------------

# ## input using list
# input_list =  list(P.2_3.Jan = list(base = c('2_3', 'prate',1), preds_ls = list(c(8,11), c(11, 4, 8), c(11,6,7,8), c(10,9,8), c(8,11,10))), 
#                    T.2_3.Jan = list(base = c('2_3', 'tmp2m',1), preds_ls = list(c(8,10), c(10, 2, 4), c(10,5))),
#                    P.2_3.Jul = list(base = c('2_3', 'prate',7), preds_ls = list(c(8,11), c(8,11,10), c(8,6,7,5), c(8,6,7,10), c(8,6,7,11))), 
#                    T.2_3.Jul = list(base = c('2_3', 'tmp2m',7), preds_ls = list(c(8,10), c(8,10,1), c(8,10,2,5), c(8,10,9))),
#                    
#                    P.3_4.Jan = list(base = c('3_4', 'prate',1), preds_ls = list(c(8,11), c(8,4,11),c(8,6,7,11))), 
#                    T.3_4.Jan = list(base = c('3_4', 'tmp2m',1), preds_ls = list(c(8,10), c(8,10,11), c(8,10,9))),
#                    P.3_4.Jul = list(base = c('3_4', 'prate',7), preds_ls = list(c(8,11), c(8,11,4), c(8,11,10))), 
#                    T.3_4.Jul = list(base = c('3_4', 'tmp2m',7), preds_ls = list(c(8,10), c(8,10,5), c(8,10,2)))
# )
# 
# ## CHANGE THESE
# pred_input = input_list[[6]]
# i.pred = 3 # start with 2 since already processed 1 for each list
# 
# 
# wk = pred_input$base[1]
# pred_var = pred_input$base[2]
# pred_mon = as.numeric(pred_input$base[3])
# 
# pred_num = length(pred_input$preds_ls)
# preds = pred_input$preds_ls[[i.pred]]

# ------------------------------------------------------------------------------


### === Directories
base_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/'

## raw data dirs
cfs_dir = paste0(base_dir, 'cfs_annual_new/')
TP_data = paste0(base_dir, '2wk_avg/')

## plsr dirs
out_model_dir = paste0(base_dir, 'plsr_input/X.Input.Processed/')
allVar_model_dir = paste0(base_dir, 'plsr_input/AllvarModels/')
allHRU_model_dir = paste0(base_dir, 'plsr_input/AllhruModels/')
allHRU_multiVar_dir = paste0(base_dir, 'plsr_input/Allhru.MultiVar_PLSR/')
model_stepwisePLSR_dir = paste0(base_dir, 'plsr_input/plsrVarModels_stepwise')
model_pcGLM_dir = paste0(base_dir, 'plsr_input/plsrVarModels_pcGLM')
model_plsrbase_dir = paste0(base_dir, 'plsr_input/plsrVarModels')
plot_dir = paste0(base_dir, 'plsr_input/plots')


### === Variable Info
var_ls = c ('hgt', 'q2m', 'prs', 'slp', 'pwt', 'uwnd', 'vwnd', 'sst', 'olr','tmp', 'prt')
long_nms = c('500 mb Geopotential height', 'Specific Humid 2m', 'Surface Pressure', 
             'Sea Level Pressure', 'Precipitable Water', 'Zonal Winds (850 mb)',
             'Meridional Winds (850 mb)', 'Sea Surface Temperature', 'Outgoing Longwave Radiation',
             'Surface Temperature 2m', 'Surface Prate')
# mon_nm    = ifelse(pred_mon == 7, 'JJA', 'DJF')       # months for model data (labeling purposes)

### === PLSR CV Info
Ntime_sub        = 31          # Define number of timesteps that each sub-period (i.e. validation) will have
Ncomp            = 2           # number of components to predict from

### === Set domain - will trim domain if not NA
lat_domain = matrix(c(25,80, -20,70, -20,70, -20,30, -20,70,  #hgt, q2m, prs, slp, pwat
                      0,80,   0,80, -20,70,  -20,20, 24,53,  24,53), #uwnd, vwnd, sst, olr, tmp, prt
                    ncol = 2, byrow = T)
lon_domain = matrix(c(100,340, 100,340, 100,340, 100,340, 100,340, #hgt, q2m, prs, slp, pwat
                      100,340, 100,340, 100,360, 100,340, 235,293, 235,293), #uwnd, vwnd, sst, olr, tmp, prt
                    ncol = 2, byrow = T)

# ---------------------------------------------------------------------------------------------------------

## Replace NAs with full domain
for (i in 1:length(var_ls)) {
  if (is.na(lat_domain[i,1])) {lat_domain[i,1] = -90}
  if (is.na(lat_domain[i,2])) {lat_domain[i,2] = 90}
  if (is.na(lon_domain[i,1])) {lon_domain[i,1] = 0}
  if (is.na(lon_domain[i,2])) {lon_domain[i,2] = 360}
}

## hru name mapping file
if (exists('hru') & length(hru) == 1) {
  hru_tb = read.table('/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/huc4_id.txt')
  hru_nm = hru_tb[which(hru_tb[,1] == hru),2]
  
  file_allVars           = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'_Allvars.rds')
}

## get month numbers for surrounding months
mon_in           = c(pred_mon-1,pred_mon,pred_mon+1)   # month to analyze
for(i in 1:length(mon_in)) { 
  if(mon_in[i] > 12) { mon_in[i] = mon_in[i] - 12 }
  if(mon_in[i] < 1) { mon_in[i] = mon_in[i] + 12 }
}


### === File Names
if (cut_domain) {
  file_data              = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.cutdomain', cut_num,'.rds')
  file_lag_in            = paste0('plsr.lagged.in.mon_', mon_nm, '_lead_', wk,'.cutdomain', cut_num,'_lag.',l,'.rds')
  file_all               = paste0('CVresults_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'_Allvars.rds')
  file_all_multiVar      = paste0('CVresults_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'_multiVars_',  
                                  paste0(preds, collapse = "."), '.rds')
  file_all_multiVarNoCV  = paste0('noCVresults_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'_multiVars_',  
                                  paste0(preds, collapse = "."), '.rds')
  file_all_extremes      = paste0('CVresults_extremes_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'_multiVars_',  
                                  paste0(preds, collapse = "."), '.rds')
  file_input_multiVar    = paste0('Input_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'_multiVars_',
                                  paste0(preds, collapse = "."), '.rds')
} else {
  file_data = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.rds')
  file_allVars         = paste0('CVresults',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_lag.',l,'_Allvars.rds')
}

pred_var_nm = ifelse(pred_var == 'prate', 'Precipitation', 'Temperature')
gsub('_', '-', wk)
title = paste(month.abb[pred_mon], 'Week', gsub('_', '-', wk),  
              ifelse(pred_var == 'prate', 'Precipitation', 'Temperature'), 'Forecast')
title_compact = paste(month.abb[pred_mon], 'Wk', gsub('_', '-', wk),  
                      ifelse(pred_var == 'prate', 'Prec.', 'Temp.'), 'Fcst')


### ==== Load Correlations for Data

## a function that returns the position of n-th largest
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
maxV <- function(n) function(x) sort(x, decreasing = TRUE)[n]

## Process output of 2a functions
process_cv_output <- function(file) {
  
  ls_in = readRDS(file)
  t = lapply(ls_in, function(x) lapply(x, identity))
  dims = c(length(t), length(t[[1]]), length(t[[1]][[1]]))
  data_all <- list()
  
  ## arrange PLSR data from parallelized output
  for (i in 2:dims[1]) { data_all[[t[[i]]$hru]] = t[[i]]  }
  for (i in 2:dims[2]) { data_all[[t[[1]][[i]]$hru]] = t[[1]][[i]]  }
  for (i in 1:dims[3]) { data_all[[t[[1]][[1]][[i]]$hru]] = t[[1]][[1]][[i]]  }
  
  return(data_all)

}

# ## see top 3 correlated predictors
# setwd(allHRU_model_dir)
# file_cor = paste0('CVcors_allVars_wk', wk, '_mon.', pred_mon, '_',pred_var,'_processed.rds')
# if (!file.exists(file_cor)){
#   ## print error
#   print('Run script 3_All_PLSR_processResults_2a.R to see highest correlated variables')
# } else {
#   cor = readRDS(file_cor)
#   
#   cor_max <- NULL
#   for (i in 1:3) {
#     cor_max = rbind.data.frame(cor_max, 
#                                cbind.data.frame(hru = cor$hru, 
#                                                 var = colnames(cor[c(4:14)])[apply(cor[c(4:14)], 1, maxn(i))],
#                                                 val = apply(cor[c(4:14)], 1, maxV(i)),
#                                                 cfs   = cor$cfs))
#   }
#   
#   # see highest three correlated gridded variables
#   hru_i = hru
#   cor_vars = dplyr::filter(cor_max, hru == as.numeric(hru_i))
#   cor_vars = na.omit(cor_vars[order(cor_vars),])
#   print(cor_vars)
#   print(c(which(var_ls == as.character(cor_vars$var)[1]),
#           which(var_ls == as.character(cor_vars$var)[2]),
#           which(var_ls == as.character(cor_vars$var)[3])))
#   
# }


