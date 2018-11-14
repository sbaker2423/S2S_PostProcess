# ================================================================================ 
# Step 3: process PLSR CV results
#
# PLSR code for processing CFSv2 variables for predicting 
# -- Process CV output or non CV output
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
library(tidyverse)
library(hydroGOF)

## Loop through wk lead, predictor, month
pred_base = c(8) # uses two preds max not including T/P
is_CV = TRUE

### === Files Names
if (is_CV) {
  nm_CV = 'CVresults'
} else {
  nm_CV = 'noCVresults'
}
df_anomMon_file = paste0('summMon_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_anomSeas_file = paste0('summSeas_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
stat_anomMon_file = paste0('summStatsMon_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
stat_anomSeas_file = paste0('summStatsSeas_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')

## Directories
dir_scripts = '/home/hydrofcst/s2s/scripts/'

### === load cfsv2 doy avg data
setwd(dir_scripts)
doy_avg_cfs = readRDS('cfsv2_doyAVG.rds') # doy of fcsted day 1 of bi-weekly period
doy_avg_cfs$lead = as.character(doy_avg_cfs$lead)
doy_avg_cfs$var_name = as.character(doy_avg_cfs$var_name)
doy_avg_cfs$hru = as.character(as.numeric(doy_avg_cfs$hru)) # remove leading 0
colnames(doy_avg_cfs)[3] <- c('Y_cfs_doyAvg')

### === load nldas doy avg data
doy_avg_nld = readRDS('nldas_doyAVG.rds') # doy of fcsted day 1 of bi-weekly period
doy_avg_nld$lead = as.character(doy_avg_nld$lead)
doy_avg_nld$var_name = as.character(doy_avg_nld$var_name)
doy_avg_nld$hru = as.character(as.numeric(doy_avg_nld$hru)) # remove leading 0
colnames(doy_avg_nld)[3] <- c('Y_nld_doyAvg')

## combine
doy_avg = left_join(doy_avg_cfs, doy_avg_nld, by = c('hru', 'doy', 'lead', 'var_name'))


#####   ===========================================================   ######
##                   Read data & Calculate Correlation                    ##
#####   ===========================================================   ######

### === Corr. by MONTH - loop through each predictand, wk period, and month === ###
df_anomMon <- NULL
setwd(allHRU_multiVar_dir)
# pred_var = 'tmp2m'; wk = '3_4'; pred_mon = 7
if (!file.exists(df_anomMon_file)){
  for (pred_var in c('prate', 'tmp2m')) {
    for (wk in c('1_2', '2_3', '3_4')) {
      # d = ifelse(wk == '1_2', 0, ifelse(wk == '2_3', 7, 14))
      
      for (pred_mon in c(1:12)) {
        
        ## Organize Predictors
        pred_2 = ifelse(pred_var == 'prate', 11, 10)
        if (length(pred_base) == 2) {
          preds = c(pred_base[1], pred_2, pred_base[2])
        } else {
          preds = c(pred_base[1], pred_2)
        }
        
        ## files to read
        file_all_multiVar = paste0(nm_CV, '_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', 
                                   cut_num,'_lag.',l,'_multiVars_', paste0(preds, collapse = "."),'.rds')
        
        ### ==== Load PLSR Model Data in usable format
        data_all = process_cv_output(file_all_multiVar)
        hru_id = names(data_all)

        ### ==== Calculate correlations
        df_cor <- NULL
        for (hru_i in hru_id) {

          ## extract data and convert to numeric
          df_i = as.data.frame(data_all[[hru_i]]$df)
          df_i$date = as.Date(df_i$date, '%Y%m%d')
          df_i$doy = as.numeric(strftime(df_i$date, format = '%j')) #forecasted date (day 1)
          # df_i$doy = as.numeric(strftime(df_i$date + d, format = '%j')) #forecasted date (day 1)
          
          df_i$doy = ifelse(df_i$doy > 365, 365, df_i$doy) # set any DOY = 366 to 365
          
          df_doyAvg_i = filter(doy_avg, hru == hru_i, lead == wk, var_name == pred_var)
        
          ## combine
          df_i = left_join(df_i[,c('doy', 'Y_plsr', 'Y_cfs', 'Y_nld')], 
                           df_doyAvg_i[,c('doy', 'Y_nld_doyAvg', 'Y_cfs_doyAvg')], by = c('doy'))
          
          ## calculate anomalies and residuals
          df_i$anom_cfs = df_i$Y_cfs - df_i$Y_cfs_doyAvg
          df_i$anom_plsr = df_i$Y_plsr - df_i$Y_cfs_doyAvg
          df_i$anom_nld = df_i$Y_nld - df_i$Y_nld_doyAvg
          df_i$res_cfs = df_i$Y_cfs - df_i$Y_nld
          df_i$res_plsr = df_i$Y_plsr - df_i$Y_nld

          df_cor = rbind.data.frame(df_cor, 
                                    cbind.data.frame(hru = hru_i, 
                                                     cor_plsr = cor(df_i$Y_plsr, df_i$Y_nld), 
                                                     cor_cfs = cor(df_i$Y_cfs, df_i$Y_nld), 
                                                     ACC_plsr = cor(df_i$anom_plsr, df_i$anom_nld), 
                                                     ACC_cfs = cor(df_i$anom_cfs, df_i$anom_nld), 
                                                     Err_plsr = mean(df_i$res_plsr),
                                                     Err_cfs = mean(df_i$res_cfs),
                                                     AbErr_plsr = mean(abs(df_i$res_plsr)),
                                                     AbErr_cfs = mean(abs(df_i$res_cfs)),
                                                     pbias_plsr = pbias(df_i$Y_plsr, df_i$Y_nld),
                                                     pbias_cfs = pbias(df_i$Y_cfs, df_i$Y_nld)
                                    ))
        }
        
        ## format df (name, to numeric)
        df_cor$inc_cor = df_cor$cor_plsr - df_cor$cor_cfs
        df_cor$inc_ACC = df_cor$ACC_plsr - df_cor$ACC_cfs
        df_cor$cor_max <- ifelse(df_cor$inc_cor > 0, df_cor$cor_plsr, df_cor$cor_cfs)
        df_cor$ACC_max <- ifelse(df_cor$inc_ACC > 0, df_cor$ACC_plsr, df_cor$ACC_cfs)
        df_cor$mon_pred = pred_mon; df_cor$pred_var = pred_var; df_cor$wk = wk
        
        ## collect data - cors for other leads
        df_anomMon = rbind(df_anomMon, df_cor)
      }
    }
  }
  
  ## save df summarizing cor at every lead/mon/var
  setwd(allHRU_multiVar_dir)
  saveRDS(df_anomMon, file = df_anomMon_file)
}


### === Anom. Corr. by SEASON - loop through each predictand, wk period, and month === ###
df_anomSeas <- NULL
setwd(allHRU_multiVar_dir)
mon_seas = c(1, 4, 7, 10)
seasons = c('DJF', 'MAM', 'JJA', 'SON')
pred_var = 'tmp2m'; wk = '1_2'; seas_i = 1
if (!file.exists(df_anomSeas_file)){
  for (pred_var in c('prate', 'tmp2m')) {
    
    ## Organize Predictors
    pred_2 = ifelse(pred_var == 'prate', 11, 10)
    if (length(pred_base) == 2) {
      preds = c(pred_base[1], pred_2, pred_base[2])
    } else {
      preds = c(pred_base[1], pred_2)
    }
    
    for (wk in c('1_2', '2_3', '3_4')) {
      d = ifelse(wk == '1_2', 0, ifelse(wk == '2_3', 7, 14))
      for (seas_i in mon_seas) {
        
        ## get season
        seas = c(seas_i - 1, seas_i, seas_i + 1)
        seas = ifelse(seas < 1, seas + 12, seas)
        seas_nm = seasons[which(mon_seas == seas_i)]
        
        ## files to read
        file_all_multiVar = paste0(nm_CV, '_',pred_var, '_allHRU_wk', wk, '_mon.', seas,'_cut.', 
                                   cut_num,'_lag.',l,'_multiVars_', paste0(preds, collapse = "."),'.rds')
        
        ### ==== Load PLSR Model Data in usable forma
        data_all.1 = process_cv_output(file_all_multiVar[1])
        data_all.2 = process_cv_output(file_all_multiVar[2])
        data_all.3 = process_cv_output(file_all_multiVar[3])
        
        hru_id.1 = names(data_all.1)
        hru_id.2 = names(data_all.2)
        hru_id.3 = names(data_all.3)
        # any False - hru's not matching
        !any(hru_id.1 == hru_id.2 | hru_id.1 == hru_id.3)
        
        ### ==== Calculate correlations
        hru_i = hru_id.1[1]
        for (hru_i in hru_id.1) {
         
          ## extract data and convert to numeric
          df_i = rbind(as.data.frame(data_all.1[[hru_i]]$df),
                       as.data.frame(data_all.2[[hru_i]]$df), 
                       as.data.frame(data_all.3[[hru_i]]$df))
          
          ## extract data and convert to numeric
          df_i$date = as.Date(df_i$date, '%Y%m%d') # starts 1999-01-16 ...
          ## should check what 'date' is!!!!
          df_i$doy = as.numeric(strftime(df_i$date, format = '%j')) #forecasted date (day 1)
          # df_i$doy = as.numeric(strftime(df_i$date + d, format = '%j')) #forecasted date (day 1)
          
          df_i$doy = ifelse(df_i$doy > 365, 365, df_i$doy) # set any DOY = 366 to 365
          
          df_doyAvg_i = filter(doy_avg, hru == hru_i, lead == wk, var_name == pred_var)
          
          ## combine
          df_i = left_join(df_i[,c('doy', 'Y_plsr', 'Y_cfs', 'Y_nld')], 
                           df_doyAvg_i[,c('doy', 'Y_nld_doyAvg', 'Y_cfs_doyAvg')], by = c('doy'))
          
          ## calculate anomalies and residuals
          df_i$anom_cfs = df_i$Y_cfs - df_i$Y_cfs_doyAvg
          df_i$anom_plsr = df_i$Y_plsr - df_i$Y_cfs_doyAvg
          df_i$anom_nld = df_i$Y_nld - df_i$Y_nld_doyAvg
          df_i$res_cfs = df_i$Y_cfs - df_i$Y_nld
          df_i$res_plsr = df_i$Y_plsr - df_i$Y_nld
          
          df_anomSeas = rbind.data.frame(df_anomSeas, 
                                        cbind.data.frame(hru = hru_i, pred_var = pred_var, wk = wk, seas = seas_nm,
                                                         cor_plsr = cor(df_i$Y_plsr, df_i$Y_nld), 
                                                         cor_cfs = cor(df_i$Y_cfs, df_i$Y_nld), 
                                                         ACC_plsr = cor(df_i$anom_plsr, df_i$anom_nld), 
                                                         ACC_cfs = cor(df_i$anom_cfs, df_i$anom_nld), 
                                                         Err_plsr = mean(df_i$res_plsr),
                                                         Err_cfs = mean(df_i$res_cfs),
                                                         AbErr_plsr = mean(abs(df_i$res_plsr)),
                                                         AbErr_cfs = mean(abs(df_i$res_cfs)),
                                                         pbias_plsr = pbias(df_i$Y_plsr, df_i$Y_nld),
                                                         pbias_cfs = pbias(df_i$Y_cfs, df_i$Y_nld)
                                        ))
        }
      }
    }
  }
  
  ## organize
  df_anomSeas$inc_cor = df_anomSeas$cor_plsr - df_anomSeas$cor_cfs
  df_anomSeas$inc_ACC = df_anomSeas$ACC_plsr - df_anomSeas$ACC_cfs
  df_anomSeas$cor_max <- ifelse(df_anomSeas$inc_cor > 0, df_anomSeas$cor_plsr, df_anomSeas$cor_cfs)
  df_anomSeas$ACC_max <- ifelse(df_anomSeas$inc_ACC > 0, df_anomSeas$ACC_plsr, df_anomSeas$ACC_cfs)
  df_anomSeas$seas = as.character(df_anomSeas$seas)
  
  ## save df summarizing cor at every lead/mon/var
  setwd(allHRU_multiVar_dir)
  saveRDS(df_anomSeas, file = df_anomSeas_file)
}


##  === Read data from files
## Month - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_anomMon = readRDS(file = df_anomMon_file)

## Season - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_anomSeas = readRDS(file = df_anomSeas_file)



#####   ===========================================================   ######
##                 Calculate overall CONUS Statistics                     ##
#####   ===========================================================   ######


#### ===== Monthly - Caclulate summary statistics 
stat_anomMon <- NULL
var_i = 'tmp2m'; wk_i = '2_3'; mon_i = 4
setwd(allHRU_multiVar_dir)
if (!file.exists(stat_anomMon_file)) {
  ## Loop through to collect plots / summary statistics
  for (var_i in c('prate', 'tmp2m')) {
    for (wk_i in c('1_2', '2_3', '3_4')) {
      for (mon_i in 1:12) {
        df_i = filter(df_anomMon, wk == wk_i, pred_var == var_i, mon_pred == mon_i) 
        
        ## increase statistic
        sum_cor = sum(df_i$inc_cor[which(df_i$inc_cor > 0)])
        num_cor = length(df_i$inc_cor[which(df_i$inc_cor > 0)])
        sum_ACC = sum(df_i$inc_ACC[which(df_i$inc_ACC > 0)])
        num_ACC = length(df_i$inc_ACC[which(df_i$inc_ACC > 0)])
        sum_stat = c(var_i, wk_i, mon_i, 
                     sum_cor, num_cor, sum_cor/num_cor, 
                     sum_ACC, num_ACC, sum_ACC/num_ACC)
        stat_anomMon = rbind(stat_anomMon, sum_stat)
      }
    }
  }
  
  ## organize df
  stat_anomMon = as.data.frame(stat_anomMon)
  colnames(stat_anomMon) <- c('var', 'wk', 'mon', 'sum_cor', 'num_cor', 'sum.num_cor', 
                              'sum_ACC', 'num_ACC', 'sum.num_ACC')
  cols.num = c('mon', 'sum_cor', 'num_cor', 'sum.num_cor', 'sum_ACC', 'num_ACC', 'sum.num_ACC')
  stat_anomMon[cols.num] <- sapply(stat_anomMon[cols.num], function(x) as.numeric(as.character(x)))
  stat_anomMon[c('var', 'wk')] <- sapply(stat_anomMon[c('var', 'wk')], function(x) as.character(x))
  stat_anomMon$mon = month.abb[stat_anomMon$mon]
  
  ## save file
  setwd(allHRU_multiVar_dir)
  saveRDS(stat_anomMon, file = stat_anomMon_file)
  
} else {
  stat_anomMon = readRDS(stat_anomMon_file)
}


#### ===== Season - Caclulate summary statistics 
stat_anomSeas <- NULL
setwd(allHRU_multiVar_dir)
if (!file.exists(stat_anomSeas_file)) {
  ## Loop through to collect plots / summary statistics
  for (var_i in c('prate', 'tmp2m')) {
    for (wk_i in c('1_2', '2_3', '3_4')) {
      for (seas_i in mon_seas) {
        
        ## get season
        seas = c(seas_i - 1, seas_i, seas_i + 1)
        seas = ifelse(seas < 1, seas + 12, seas)
        seas_nm = seasons[which(mon_seas == seas_i)]
        
        df_i = filter(df_anomSeas, wk == wk_i, pred_var == var_i, seas == seas_nm) 
        
        ## increase statistic
        sum_cor = sum(df_i$inc_cor[which(df_i$inc_cor > 0)])
        num_cor = length(df_i$inc_cor[which(df_i$inc_cor > 0)])
        sum_ACC = sum(df_i$inc_ACC[which(df_i$inc_ACC > 0)])
        num_ACC = length(df_i$inc_ACC[which(df_i$inc_ACC > 0)])
        sum_stat = c(var_i, wk_i, seas_nm, 
                     sum_cor, num_cor, sum_cor/num_cor, 
                     sum_ACC, num_ACC, sum_ACC/num_ACC)
        stat_anomSeas = rbind(stat_anomSeas, sum_stat)
        
      }
    }
  }
  
  ## organize df
  stat_anomSeas = as.data.frame(stat_anomSeas)
  colnames(stat_anomSeas) <- c('var', 'wk', 'seas', 'sum_cor', 'num_cor', 'sum.num_cor', 
                              'sum_ACC', 'num_ACC', 'sum.num_ACC')
  cols.num = c('sum_cor', 'num_cor', 'sum.num_cor', 'sum_ACC', 'num_ACC', 'sum.num_ACC')
  stat_anomSeas[cols.num] <- sapply(stat_anomSeas[cols.num], function(x) as.numeric(as.character(x)))
  stat_anomSeas[c('var', 'wk', 'seas')] <- sapply(stat_anomSeas[c('var', 'wk', 'seas')], function(x) as.character(x))
  
  ## save file
  setwd(allHRU_multiVar_dir)
  saveRDS(stat_anomSeas, file = stat_anomSeas_file)
  
} else {
  stat_anomSeas = readRDS(stat_anomSeas_file)
}



#####   ===========================================================   ######
##                             Compare predictors                         ##
#####   ===========================================================   ######
pred_base_comp = c(paste0(c(8,4), collapse = "."), 
                   paste0(c(8), collapse = "."))


## Monthly Analysis
setwd(allHRU_multiVar_dir)
stat_anomMon_file = paste0('summStatsMon_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', pred_base_comp, '.rds')
stat_anomMon_1 = readRDS(stat_anomMon_file[1])
stat_anomMon_2 = readRDS(stat_anomMon_file[2])

stat_anomMon = left_join(stat_anomMon_1, stat_anomMon_2, by = c('var', 'wk', 'mon'), 
                       suffix = c(pred_base_comp[1], pred_base_comp[2]))
# cor
stat_anomMon$num_max_cor = apply(stat_anomMon[paste0('num_cor', pred_base_comp)],1, which.max)
stat_anomMon$sum.num_max_cor = apply(stat_anomMon[paste0('sum.num_cor', pred_base_comp)],1, which.max)
stat_anomMon$best_preds_cor = ifelse(stat_anomMon$num_max_cor == stat_anomMon$sum.num_max_cor, stat_anomMon$num_max_cor, stat_anomMon$sum.num_max_cor)
# ACC
stat_anomMon$num_max_ACC = apply(stat_anomMon[paste0('num_ACC', pred_base_comp)],1, which.max)
stat_anomMon$sum.num_max_ACC = apply(stat_anomMon[paste0('sum.num_ACC', pred_base_comp)],1, which.max)
stat_anomMon$best_preds_ACC = ifelse(stat_anomMon$num_max_ACC == stat_anomMon$sum.num_max_ACC, stat_anomMon$num_max_ACC, stat_anomMon$sum.num_max_ACC)

## summarize results
table(stat_anomMon$best_preds_cor) # 2 is the best here (just SST useds as pred)
table(stat_anomMon$best_preds_ACC)

# save 
setwd(plot_dir)
saveRDS(stat_anomMon, file = paste0('summStatsMon_', nm_CV, '_allLeadsMonsVars_predBase.Compare.rds'))


## Seasonal Analysis
setwd(allHRU_multiVar_dir)
stat_anomSeas_fnames = paste0('summStatsSeas_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', pred_base_comp, '.rds')

stat_anomSeas_1 = readRDS(stat_anomSeas_fnames[1])
stat_anomSeas_2 = readRDS(stat_anomSeas_fnames[2])

stat_anomSeas = left_join(stat_anomSeas_1, stat_anomSeas_2, by = c('var', 'wk', 'seas'), 
                       suffix = c(pred_base_comp[1], pred_base_comp[2]))
stat_anomSeas[is.na(stat_anomSeas)] <- 0 # remove nans that cause errors

# cor
stat_anomSeas$num_max_cor = apply(stat_anomSeas[paste0('num_cor', pred_base_comp)],1, which.max)
stat_anomSeas$sum.num_max_cor = apply(stat_anomSeas[paste0('sum.num_cor', pred_base_comp)],1, which.max)
stat_anomSeas$best_preds_cor = ifelse(stat_anomSeas$num_max_cor == stat_anomSeas$sum.num_max_cor, stat_anomSeas$num_max_cor, stat_anomSeas$sum.num_max_cor)
# ACC
stat_anomSeas$num_max_ACC = apply(stat_anomSeas[paste0('num_ACC', pred_base_comp)],1, which.max)
stat_anomSeas$sum.num_max_ACC = apply(stat_anomSeas[paste0('sum.num_ACC', pred_base_comp)],1, which.max)
stat_anomSeas$best_preds_ACC = ifelse(stat_anomSeas$num_max_ACC == stat_anomSeas$sum.num_max_ACC, stat_anomSeas$num_max_ACC, stat_anomSeas$sum.num_max_ACC)

## summarize results
table(stat_anomSeas$best_preds_cor) # 2 is the best here (just SST useds as pred)
table(stat_anomSeas$best_preds_ACC)

# save 
setwd(plot_dir)
saveRDS(stat_anomSeas, file = paste0('summStats_Anom_Seas_', nm_CV, '_allLeadsMonsVars_predBase.Compare.rds'))

