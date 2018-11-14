# ================================================================================ 
# Step 3: process PLSR CV results
#
# PLSR code for processing CFSv2 variables for predicting 
# -- Process CV output
# -- calc cor for each month and season for all hrus
#
#   Created by S. Baker 
#
# ================================================================================ 
rm(list=ls())

### === Directories
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/plotLoadings.R')
source('~/s2s/analysis/scripts/cfsv2_analysis/plot_function.r')
source("/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/0_control_PLSRmdl.R")
library(dplyr)
# library(grid)
# library(gridExtra)

## Loop through wk lead, predictor, month
pred_base = c(8,4) # uses two preds max not including T/P

## Files
df_allMon_file = paste0('summMon_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_allSeas_file = paste0('summSeas_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_statMon_file = paste0('summStatsMon_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_statSeas_file = paste0('summStatsSeas_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')


#####   ===========================================================   ######
##                   Read data & Calculate Correlation                    ##
#####   ===========================================================   ######

### === Corr. by month - loop through each predictand, wk period, and month === ###
df_allMon <- NULL
setwd(allHRU_multiVar_dir)
# pred_var = 'prate'; wk = '2_3'; pred_mon = 1
if (!file.exists(df_allMon_file)){
  for (pred_var in c('prate', 'tmp2m')) {
    for (wk in c('2_3', '3_4')) {
      for (pred_mon in c(1:12)) {
        
        ## Organize Predictors
        pred_2 = ifelse(pred_var == 'prate', 11, 10)
        if (length(pred_base) == 2) {
          preds = c(pred_base[1], pred_2, pred_base[2])
        } else {
          preds = c(pred_base[1], pred_2)
        }
        
        ## files to read
        file_all_multiVar = paste0('CVresults_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', 
                                   cut_num,'_lag.',l,'_multiVars_', paste0(preds, collapse = "."),'.rds')
        
        ### ==== Load PLSR Model Data in usable format
        data_all = process_cv_output(file_all_multiVar)
        hru_id = names(data_all)

        ### ==== Calculate correlations
        df_cor <- NULL
        for (hru_i in hru_id) {
          nI = which(hru_i == hru_id)
          
          ## extract data and convert to numeric
          df_i = as.data.frame(data_all[[hru_i]]$df)
          
          ## calc correlations and find max cor
          cor_pred = cor(df_i$Y_plsr, df_i$Y_nld)
          cor_cfs = cor(df_i$Y_cfs, df_i$Y_nld)
          
          # add to matrix
          df_cor = rbind.data.frame(df_cor, 
                                    cbind.data.frame(hru = hru_i, cor_plsr = cor_pred, cor_cfs = cor_cfs))
        }
        
        ## format df (name, to numeric)
        df_cor$inc_cor = df_cor$cor_plsr - df_cor$cor_cfs
        df_cor$val_max <- ifelse(df_cor$inc_cor > 0, df_cor$cor_plsr, df_cor$cor_cfs)
        df_cor$mon_pred = pred_mon; df_cor$pred_var = pred_var; df_cor$wk = wk
        
        ## collect data - cors for other leads
        df_allMon = rbind(df_allMon, df_cor)
      }
    }
  }
  
  ## save df summarizing cor at every lead/mon/var
  setwd(allHRU_multiVar_dir)
  saveRDS(df_allMon, file = df_allMon_file)
}


### === Corr. by season - loop through each predictand, wk period, and month === ###
df_allSeas <- NULL
setwd(allHRU_multiVar_dir)
mon_seas = c(1, 4, 7, 10)
seasons = c('DJF', 'MAM', 'JJA', 'SON')

if (!file.exists(df_allSeas_file)){
  for (pred_var in c('prate', 'tmp2m')) {
    
    ## Organize Predictors
    pred_2 = ifelse(pred_var == 'prate', 11, 10)
    if (length(pred_base) == 2) {
      preds = c(pred_base[1], pred_2, pred_base[2])
    } else {
      preds = c(pred_base[1], pred_2)
    }
    
    for (wk in c('2_3', '3_4')) {
      for (seas_i in mon_seas) {
        
        ## get season
        seas = c(seas_i - 1, seas_i, seas_i + 1)
        seas = ifelse(seas < 1, seas + 12, seas)
        seas_nm = seasons[which(mon_seas == seas_i)]
        
        ## files to read
        file_all_multiVar = paste0('CVresults_',pred_var, '_allHRU_wk', wk, '_mon.', seas,'_cut.', 
                                   cut_num,'_lag.',l,'_multiVars_', paste0(preds, collapse = "."),'.rds')
        
        ### ==== Load PLSR Model Data in usable forma
        data_all.1 = process_cv_output(file_all_multiVar[1])
        data_all.2 = process_cv_output(file_all_multiVar[2])
        data_all.3 = process_cv_output(file_all_multiVar[3])
        
        hru_id.1 = names(data_all.1)
        hru_id.2 = names(data_all.2)
        hru_id.3 = names(data_all.3)
        
        ### ==== Calculate correlations
        for (hru_i in 1:length(hru_id.1)) {
          nI.1 = hru_id.1[hru_i]; nI.2 = hru_id.2[hru_i]; nI.3 = hru_id.3[hru_i]
          
          if (nI.1 != nI.2 | nI.1 != nI.3) { print('Error -- HRU indices not equal')}
          
          ## extract data and convert to numeric
          df_i = rbind(as.data.frame(data_all.1[[nI.1]]$df),
                       as.data.frame(data_all.2[[nI.2]]$df), 
                       as.data.frame(data_all.3[[nI.3]]$df))
          
          ## calc correlations and find max cor
          cor_pred = cor(df_i$Y_plsr, df_i$Y_nld)
          cor_cfs = cor(df_i$Y_cfs, df_i$Y_nld)
          
          # add to matrix
          df_allSeas = rbind.data.frame(df_allSeas,
                                      cbind.data.frame(hru = nI.1, pred_var = pred_var, wk = wk, seas = seas_nm,
                                                       cor_plsr = cor_pred, cor_cfs = cor_cfs))
        }
      }
    }
  }
  
  ## organize
  df_allSeas$inc_cor = df_allSeas$cor_plsr - df_allSeas$cor_cfs
  df_allSeas$val_max <- ifelse(df_allSeas$inc_cor > 0, df_allSeas$cor_plsr, df_allSeas$cor_cfs)
  df_allSeas$seas = as.character(df_allSeas$seas)
  
  ## save df summarizing cor at every lead/mon/var
  setwd(allHRU_multiVar_dir)
  saveRDS(df_allSeas, file = df_allSeas_file)
}


##  === Read data from files
## Month - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_allMon = readRDS(file = df_allMon_file)

## Season - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_allSeas = readRDS(file = df_allSeas_file)



#####   ===========================================================   ######
##                      Calculate overal Statistics                       ##
#####   ===========================================================   ######


#### ===== Monthly - Caclulate summary statistics 
df_statMon <- NULL
var_i = 'tmp2m'; wk_i = '2_3'; mon_i = 4
setwd(allHRU_multiVar_dir)
if (!file.exists(df_statMon_file)) {
  ## Loop through to collect plots / summary statistics
  for (var_i in c('prate', 'tmp2m')) {
    for (wk_i in c('2_3', '3_4')) {
      for (mon_i in 1:12) {
        df_i = filter(df_allMon, wk == wk_i, pred_var == var_i, mon_pred == mon_i) 
        
        ## increase statistic
        sum_cor = sum(df_i$inc_cor[which(df_i$inc_cor > 0)])
        num_cor = length(df_i$inc_cor[which(df_i$inc_cor > 0)])
        sum_stat = c(var_i, wk_i, mon_i, sum_cor, num_cor, sum_cor/num_cor)
        df_statMon = rbind(df_statMon, sum_stat)
      }
    }
  }
  
  ## organize df
  df_statMon = as.data.frame(df_statMon)
  colnames(df_statMon) <- c('var', 'wk', 'mon', 'sum_cor', 'num_cor', 'sum.num')
  cols.num = c('mon', 'sum_cor', 'num_cor', 'sum.num')
  df_statMon[cols.num] <- sapply(df_statMon[cols.num], function(x) as.numeric(as.character(x)))
  df_statMon[c('var', 'wk')] <- sapply(df_statMon[c('var', 'wk')], function(x) as.character(x))
  df_statMon$mon = month.abb[df_statMon$mon]
  
  ## save file
  setwd(allHRU_multiVar_dir)
  saveRDS(df_statMon, file = df_statMon_file)
  
} else {
  df_statMon = readRDS(df_statMon_file)
}


#### ===== Season - Caclulate summary statistics 
df_statSeas <- NULL
setwd(allHRU_multiVar_dir)
if (!file.exists(df_statSeas_file)) {
  ## Loop through to collect plots / summary statistics
  for (var_i in c('prate', 'tmp2m')) {
    for (wk_i in c('2_3', '3_4')) {
      for (seas_i in mon_seas) {
        
        ## get season
        seas = c(seas_i - 1, seas_i, seas_i + 1)
        seas = ifelse(seas < 1, seas + 12, seas)
        seas_nm = seasons[which(mon_seas == seas_i)]
        
        df_i = filter(df_allSeas, wk == wk_i, pred_var == var_i, seas == seas_nm) 
        
        ## increase statistic
        sum_cor = sum(df_i$inc_cor[which(df_i$inc_cor > 0)])
        num_cor = length(df_i$inc_cor[which(df_i$inc_cor > 0)])
        sum_stat = c(var_i, wk_i, seas_nm, sum_cor, num_cor, sum_cor/num_cor)
        df_statSeas = rbind(df_statSeas, sum_stat)
      }
    }
  }
  
  ## organize df
  df_statSeas = as.data.frame(df_statSeas)
  colnames(df_statSeas) <- c('var', 'wk', 'seas', 'sum_cor', 'num_cor', 'sum.num')
  cols.num = c('sum_cor', 'num_cor', 'sum.num')
  df_statSeas[cols.num] <- sapply(df_statSeas[cols.num], function(x) as.numeric(as.character(x)))
  df_statSeas[c('var', 'wk', 'seas')] <- sapply(df_statSeas[c('var', 'wk', 'seas')], function(x) as.character(x))
  
  ## save file
  setwd(allHRU_multiVar_dir)
  saveRDS(df_statSeas, file = df_statSeas_file)
  
} else {
  df_statSeas = readRDS(df_statSeas_file)
}



#####   ===========================================================   ######
##                      Analyze data processed above                      ##
#####   ===========================================================   ######
pred_base_comp = c(paste0(c(8,4), collapse = "."), 
                   paste0(c(8), collapse = "."))


## Monthly Analysis
setwd(allHRU_multiVar_dir)
df_statMon_1 = readRDS(paste0('summStatsMon_CVresults_allLeadsMonsVars_predBase.', pred_base_comp[1], '.rds'))
df_statMon_2 = readRDS(paste0('summStatsMon_CVresults_allLeadsMonsVars_predBase.', pred_base_comp[2], '.rds'))

df_statMon = left_join(df_statMon_1, df_statMon_2, by = c('var', 'wk', 'mon'), 
                       suffix = c(pred_base_comp[1], pred_base_comp[2]))

df_statMon$num_max = apply(df_statMon[paste0('num_cor', pred_base_comp)],1, which.max)
df_statMon$sum.num_max = apply(df_statMon[paste0('sum.num', pred_base_comp)],1, which.max)
df_statMon$best_preds = ifelse(df_statMon$num_max == df_statMon$sum.num_max, df_statMon$num_max, df_statMon$sum.num_max)

## summarize results
table(df_statMon$best_preds) # 2 is the best here (just SST useds as pred)

# save 
setwd(plot_dir)
saveRDS(df_statMon, file = paste0('summStatsMon_CVresults_allLeadsMonsVars_predBase.Compare.rds'))


## Seasonal Analysis
setwd(allHRU_multiVar_dir)
df_statMon_1 = readRDS(paste0('summStatsSeas_CVresults_allLeadsMonsVars_predBase.', pred_base_comp[1], '.rds'))
df_statMon_2 = readRDS(paste0('summStatsSeas_CVresults_allLeadsMonsVars_predBase.', pred_base_comp[2], '.rds'))

df_statMon = left_join(df_statMon_1, df_statMon_2, by = c('var', 'wk', 'seas'), 
                       suffix = c(pred_base_comp[1], pred_base_comp[2]))

df_statMon$num_max = apply(df_statMon[paste0('num_cor', pred_base_comp)],1, which.max)
df_statMon$sum.num_max = apply(df_statMon[paste0('sum.num', pred_base_comp)],1, which.max)
df_statMon$best_preds = ifelse(df_statMon$num_max == df_statMon$sum.num_max, df_statMon$num_max, df_statMon$sum.num_max)

## summarize results
table(df_statMon$best_preds) # 2 is the best here (just SST useds as pred)

# save 
setwd(plot_dir)
saveRDS(df_statMon, file = paste0('summStatsSeas_CVresults_allLeadsMonsVars_predBase.Compare.rds'))
