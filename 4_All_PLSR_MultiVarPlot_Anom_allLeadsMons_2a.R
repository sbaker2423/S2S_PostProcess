# ================================================================================ 
# Step 4: process PLSR CV results
#
# PLSR code for processing CFSv2 variables for predicting 
# -- Process CV output
# -- plot cor for all hrus
# -- plot for each HRU: Mean Loadings & Timeseries 
#
#   Created by S. Baker 
#
# ================================================================================ 
rm(list=ls())

### === Directories
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/plotLoadings.R')
source('~/s2s/analysis/scripts/cfsv2_analysis/plot_function.r')
source("/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/0_control_PLSRmdl.R")
detach('package:tidyverse', unload = T)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

## Inputs: predictors, plot cors, plot hru loadings
pred_base = c(8) # uses two preds max not including T/P
hru_input   = c('1802', '1304', '1306', '1401', '1402', '1709', '303', '1014', '1302', '1704', '1002', '601', '514')
save_cor_plots = T
save_hru_plots = T
save_loadings = F # only when new hucs or in Vars
is_CV = T

##  === Read data from files
nm_CV = ifelse(is_CV, 'CVresults', 'noCVresults')
df_anomMon_file = paste0('summMon_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_anomSeas_file = paste0('summSeas_Anom_', nm_CV, '_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')

## Month - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_anomMon = readRDS(file = df_anomMon_file)

## Season - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_anomSeas = readRDS(file = df_anomSeas_file)


#####   ===========================================================   ######
##                      Plot MONTHLY Correlation Figs                     ##
#####   ===========================================================   ######
stat_type_v = c('ACC', 'AbErr')  #c('cor', 'ACC', 'AbErr') #'Err', 
time_p_v = c('Seas', 'Mon')   # c('Mon', 'Seas')
pal = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
color_fun <- colorRampPalette(pal, space = "rgb") # for discrete

## plot and save
if (save_cor_plots) {
  setwd(plot_dir)
  
  ## === Time period loop vec
  time_p = time_p_v[1]
  for(time_p in time_p_v) {
    
    ## get dataframe for plotting
    if(time_p == 'Mon') {
      df_anom = df_anomMon
      plot_set = c(16, 4, 100)
    } else {
      df_anom = df_anomSeas
      plot_set = c(7, 4, 200)
    }
    
    ## === Statistic loop vec
    stat_type = stat_type_v[2]
    for (stat_type in stat_type_v) {
      
      df_anom$best_fcst = ifelse(df_anom$ACC_cfs > df_anom$ACC_plsr, 'cfs', 'plsr')
      
      ## === Set bins and colors based on plotted statistic
      if (stat_type == 'ACC' | stat_type == 'cor') {
        
        ##  === Correlation Increase Plot - set bins/colors
        bins_Increase = c(seq(from = 0, to = 0.3, by = 0.05), 0.4, 0.5, 0.7, 0.9)
        vec = df_anom[[paste0('inc_',stat_type)]]
        df_anom$bins_inc = cut(vec, breaks = bins_Increase) 
        df_anom$bins_inc <- addNA(df_anom$bins_inc)
        levels(df_anom$bins_inc)[is.na(levels(df_anom$bins_inc))] <- "None"
        ## color palette and plot
        # pal = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
        # color_fun <- colorRampPalette(pal, space = "rgb") # for discrete
        cust_pal_inc = c(color_fun(length(levels(df_anom$bins_inc)) -1 ), 'grey85')  # add none color
        lab_inc = gsub(',',' - ', levels(df_anom$bins_inc)) # replace legend labels
        
        ##  === Correlation Highest and Raw Plots - set bins/colors
        bin_set = c(-0.2, seq(0, 0.5, by = 0.05), 0.6, 0.7, 0.8, 1)
        vec = df_anom[[paste0(stat_type, '_max')]]
        df_anom$bins_base = cut(vec, breaks = bin_set)
        n_col = length(levels(df_anom$bins_base))
        cust_pal_base = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
        lab_base = gsub(',',' - ', levels(df_anom$bins_base)) # replace legend labels
        
        legend_title = c('ACC', 'ACC', 'ACC Inc.')
      } else if (stat_type == 'AbErr') {
        
        ## process statistic
        df_stat = cbind.data.frame(plsr_val = df_anom[[paste0(stat_type, '_plsr')]], 
                                   cfs_val  =  df_anom[[paste0(stat_type, '_cfs')]],
                                   best_mdl = df_anom$best_fcst,
                                   decr_val = NA)
        df_stat$best_val = ifelse(df_stat$best_mdl == 'plsr', df_stat$plsr_val, df_stat$cfs_val)
        df_stat$decr_val = ifelse(df_stat$best_mdl == 'plsr', 
                                  df_stat$cfs_val - df_stat$best_val, df_stat$decr_val) # decrease in MAE
        # df_stat$decr_val = df_stat$cfs_val - df_stat$best_val 
        # df_stat$decr_val = ifelse(df_stat$decr_val == 0, NULL, df_stat$decr_val)
        df_anom[[paste0(stat_type, '_max')]] <- df_stat$best_val
        df_anom[[paste0('inc_', stat_type)]] <- df_stat$decr_val
        
        ##  === Increase Plot
        bins_Increase = seq(from = -0.25, to = 2.5, by = 0.25)
        df_anom$bins_inc = cut(df_stat$decr_val, breaks = bins_Increase) 
        df_anom$bins_inc <- addNA(df_anom$bins_inc)
        levels(df_anom$bins_inc)[is.na(levels(df_anom$bins_inc))] <- "None"
        pal_mae_inc = rev(c("#5E4FA2", "#4470B1", "#3B92B8", "#59B4AA", "#7ECBA4", "#A6DBA4", "#CAE99D", "#E8F69C", "#F7FCB3", "#FEE391"))
        color_fun <- colorRampPalette(pal_mae_inc, space = "rgb") # for discrete
        cust_pal_inc = c(color_fun(length(levels(df_anom$bins_inc)) -1 ), 'grey85')  # add none color
        lab_inc = gsub(',',' - ', levels(df_anom$bins_inc)) # replace legend labels
        
        ##  === Correlation Highest and Raw Plots - set bins/colors
        bin_set = seq(0, 4.5, by = 0.5)
        df_anom$bins_base = cut(df_stat$best_val, breaks = bin_set)
        n_col = length(levels(df_anom$bins_base))
        pal_mae = c("#FFFFBF", "#FDD380", "#F88D51", "#DC494C", "#9E0142")
        color_fun <- colorRampPalette(pal_mae, space = "rgb") # for discrete
        cust_pal_base = color_fun(length(levels(df_anom$bins_base))) #rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
        lab_base = gsub(',',' - ', levels(df_anom$bins_base)) # replace legend labels
        
        legend_title = c('MAE', 'MAE', 'Decr. in MAE')
      }
      
      ## === PLOT - Loop through to collect plots / summary statistics
      # var_i = 'prate'; wk_i = '3_4'
      for (var_i in c('prate', 'tmp2m')) {
        for (wk_i in c('1_2', '2_3', '3_4')) {
          
          ## Filter and plot name
          df_i = filter(df_anom, wk == wk_i, pred_var == var_i)
          if (time_p == 'Mon') {
            df_i$mon_pred = factor(month.abb[df_i$mon_pred], levels = month.abb)
          } else {
            df_i$seas = factor(df_i$seas, levels = c('DJF', 'MAM', 'JJA', 'SON'))
          }
          
          ## Plot raw CFSv2 correlation
          df_i$bins = cut(df_i[[paste0(stat_type, '_cfs')]], breaks = bin_set)
          rawStat = plot_map(var_i, df_i, paste0(stat_type, '_cfs'), type = 'discrete', lab = lab_base, 
                            cust_pal = cust_pal_base, legend_title = legend_title[1])  +
            theme(legend.position="none")
          
          ## Plot PLSR correlation
          df_i$bins = cut(df_i[[paste0(stat_type, '_plsr')]], breaks = bin_set)
          plsrStat = plot_map(var_i, df_i, paste0(stat_type, '_plsr'), type = 'discrete', lab = lab_base, 
                             cust_pal = cust_pal_base, legend_title = legend_title[1])  +
            theme(legend.position="none")
          
          ## Plot max correlation
          df_i$bins = df_i$bins_base
          maxStat = plot_map(var_i, df_i, paste0(stat_type, '_max'), type = 'discrete', lab = lab_base, 
                            cust_pal = cust_pal_base, legend_title = legend_title[2]) +
            theme(legend.position="none")
          
          ## Plot increase in correlation
          df_i$bins = df_i$bins_inc
          statInc = plot_map(var_i, df_i, paste0('inc_', stat_type), type = 'discrete', lab = lab_inc, 
                            cust_pal = cust_pal_inc, legend_title = legend_title[3]) +
            theme(legend.position="none")
          
          ## print and save figures
          if (time_p == 'Mon') {
            statInc = statInc + facet_grid(mon_pred ~ .)
            rawStat = rawStat + facet_grid(mon_pred ~ .)
            maxStat = maxStat + facet_grid(mon_pred ~ .)
            plsrStat = plsrStat + facet_grid(mon_pred ~ .)
          } else {
            statInc = statInc + facet_grid(seas ~ .)
            rawStat = rawStat + facet_grid(seas ~ .)
            maxStat = maxStat + facet_grid(seas ~ .)
            plsrStat = plsrStat + facet_grid(seas ~ .)
          }
          
          ## Save images
          img_nm_base = paste0(stat_type,'_plsrResults_', nm_CV, '_', time_p,'_wk', wk_i, '_', var_i,'_multiVar_', paste0(pred_base, collapse = "."), '.png')
          ggsave(paste0('inc', img_nm_base), statInc, height = plot_set[1], width = plot_set[2], dpi = plot_set[3])
          ggsave(paste0('raw', img_nm_base), rawStat, height = plot_set[1], width = plot_set[2], dpi = plot_set[3])
          ggsave(paste0('plsr', img_nm_base), plsrStat, height = plot_set[1], width = plot_set[2], dpi = plot_set[3])
          ggsave(paste0('max', img_nm_base), maxStat, height = plot_set[1], width = plot_set[2], dpi = plot_set[3])
          
          ## For saving legends - need to have plot with legend
          # ggsave('legend_sample_MAE_base.png', maxStat, width = 5, dp = plot_set[3])
          # ggsave('legend_sample_MAE_inc.png', statInc, width = 5, dp = plot_set[3])
        }
      }
      
      ## === Histogram plots of skill or error metrics
      if (time_p == 'Mon') {
        df_boxplot = df_anom %>% select(wk, mon_pred, paste0(stat_type, '_plsr'), paste0(stat_type, '_cfs'), paste0(stat_type, '_max'))
      } else {
        df_boxplot = df_anom %>% select(wk, seas, paste0(stat_type, '_plsr'), paste0(stat_type, '_cfs'), paste0(stat_type, '_max'))
      }
      colnames(df_boxplot)[c(2:5)] <- c('time_per', 'PLSR', 'CFSv2', 'Best Model')
      df_boxplot = tidyr::gather(df_boxplot, 'Model', 'stat', 3:5)
      
      ## plot
      g <- ggplot(df_boxplot, aes(x = wk, y = stat, fill = Model)) +
        geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed') +
        theme_bw() +
        geom_boxplot() + labs(y = legend_title[1], x = 'Bi-weekly Period') + 
        facet_grid(time_per ~ .)
      ## save
      ggsave(paste0('boxplot_summary_', stat_type, '.png'), g, height = 8, width = 4, dpi = 200)
    }
  }
}

#####   ===========================================================   ######
##             Plot individual HUC's - Loading and timeseries             ##
#####   ===========================================================   ######

### ==== Load PLSR input data / grid data
var_i = 'tmp2m'; wk_i = '1_2'; pred_mon_i = 1

for (var_i in c('prate', 'tmp2m')) {
  for (wk_i in c('2_3', '3_4')) {
    for (pred_mon_i in c(1, 7)) {
      
      ## Organize Predictors
      pred_2 = ifelse(var_i == 'prate', 11, 10)
      if (length(pred_base) == 2) {
        preds = c(pred_base[1], pred_2, pred_base[2])
      } else {
        preds = c(pred_base[1], pred_2)
      }
      
      
      ### ==== Load and process PLSR Model Data to a usable format
      setwd(allHRU_multiVar_dir)
      file_input_multiVar    = paste0('Input_',var_i, '_allHRU_wk', wk_i, '_mon.', pred_mon_i,'_cut.', cut_num,'_lag.',l,'_multiVars_',
                                      paste0(preds, collapse = "."), '.rds')
      file_all_multiVar      = paste0('CVresults_',var_i, '_allHRU_wk', wk_i, '_mon.', pred_mon_i,'_cut.', cut_num,'_lag.',l,'_multiVars_',  
                                      paste0(preds, collapse = "."), '.rds')
      f.out = readRDS(file = file_input_multiVar)
      data_all = process_cv_output(file_all_multiVar)
      
      pred_var_nm = ifelse(var_i == 'prate', 'Precipitation', 'Temperature')
      wk_nm = ifelse(wk_i == '2_3', '2-3', ifelse(wk_i == '3_4', '3-4', '1-2'))
      
      ### === ONE HRU: Loop through to plot loadings and timeseries ==== ###
      setwd(paste0(plot_dir, '/huc_individual_results'))
      hru_n = '1014'
      if (save_hru_plots) {
        for (hru_n in hru_input) {
          hru_nm = hru_tb[which(hru_tb[,1] == hru_n),2]
          
          ## model performance info
          df_n = filter(df_anomMon, wk == wk_i, pred_var == var_i, mon_pred == pred_mon_i, hru == hru_n)
          
          ## set bins
          max = 35 #round(max(ld, abs(min(ld))), -1)
          dif = round(2*(max)/10)
          bin_var = seq(-max, max, dif)
          
          ### === Plot Loadings
          if (save_loadings) {
            plot_ls = list(); n = 0
            for (i in 1:length(preds)) {
              for (comp_i in 1:Ncomp) {
                
                ## variable information
                n = n + 1
                ipred = preds[i]
                var = var_ls[ipred]
                var_long = long_nms[ipred]
                
                ## extract grid
                grid_p = (f.out$grid[ipred, , ])
                plot_grid = na.omit(data.frame(x = grid_p[1,], y = grid_p[2,]))
                
                ## loadings - loading for igroup vs. average loading of all igroups
                igroup = 1
                ind = c((f.out$ind_pred[i]+1):f.out$ind_pred[(i+1)])
                ld = na.omit(data_all[[hru_n]]$load[igroup, ind, comp_i])
                ld_avg = na.omit(apply(data_all[[hru_n]]$load[1:length(data_all[[hru_n]]$load[,1,1]), ind, comp_i], 2, mean))
                
                df_plot = cbind.data.frame(var = ld_avg, # ld_avg
                                           lat_p = plot_grid$x, lon_p = plot_grid$y)
                
                plot_ls[[n]] =
                  plotLoadings(df_plot, var, #max = max, min = -max,
                               bin_vec = bin_var, title_in = paste0(var_long, ' - Comp. ', comp_i))  +
                  theme(legend.position = 'none')
              }
            }
            
            ## plot & save
            title = paste('Mean Loadings -', wk_nm, 'wk Forecast of', month.abb[pred_mon_i], pred_var_nm, '-', hru_nm)
            g <- gridExtra::grid.arrange(grobs = plot_ls, nrow = length(preds),top = textGrob(title, gp=gpar(fontsize=15,font=8)))
            ggsave(paste0('plsrLoadings_wk', wk_i, '_mon.', pred_mon_i, '_',var_i,'_multiVar_huc', hru_n, '_', paste0(preds, collapse = "."), '.png'),
                   g, height = 6, width = 11, dpi = 100)
            # ggsave('loading_legend_ex.png', p, height = 5, dpi = 200) # save legend
          }
          
          ### ==== One HRU: 1:1 plots, timeseries ==== ###
          
          ## extract data and convert to numeric
          df_i = as.data.frame(data_all[[hru_n]]$df)
          if (var_i == 'tmp2m') { df_i[,3:5] = df_i[,3:5] - 273.15 }
          paste0(var_ls[preds], collapse = ', ')
          ymax = max(df_i$Y_nld, df_i$Y_cfs, df_i$Y_plsr)
          ymin = ifelse(var_i == 'tmp2m', min(df_i$Y_nld, df_i$Y_cfs, df_i$Y_plsr), 0)
          
          ## plot two 1:1 plots
          g1 <- ggplot(df_i, aes(Y_nld, Y_plsr)) + geom_abline(color = 'red') + geom_point(size = 0.5) + coord_equal() +
            xlim(ymin,ymax) + ylim(ymin,ymax) + ylab(paste('PLSR', pred_var_nm)) + xlab(paste('NLDAS', pred_var_nm))
          g2 <- ggplot(df_i, aes(Y_nld, Y_cfs)) + geom_abline(color = 'red') + geom_point(size = 0.5) + coord_equal() +
            xlim(ymin,ymax) + ylim(ymin,ymax) + ylab(paste('CFSv2', pred_var_nm)) + xlab(paste('NLDAS', pred_var_nm))
          
          ## plot timeseries
          dplot = reshape2::melt(df_i[,c(1:5)], id.vars = c('date','ind'), na.rm = T)
          dplot$variable = factor(dplot$variable, c('Y_cfs', 'Y_nld', 'Y_plsr'))
          g3 = ggplot(dplot, aes(x = ind, y = value, group = variable, colour = variable)) +
            geom_line() +
            # xlab('Index of January Day (1999 - 2010)') +
            xlab('Day Index') +
            ylab(paste0(pred_var_nm)) +
            # ylab('Precipitation (mm/day per bi-weekly period)') +
            scale_color_manual(values = c("light blue", "black", 'red')) +
            scale_linetype_manual(values = c('dotted', 'solid', 'solid'))
          
          ## export and save
          g <- gridExtra::grid.arrange(grobs = list(g1, g2), nrow = 1,
                                       top = paste(wk_nm, 'wk Forecast of', month.abb[pred_mon_i], pred_var_nm, '-', hru_nm, 
                                                   '(', paste(var_ls[preds], collapse = ', '), ')'),
                                       bottom = paste(' ACC: Raw =', round(df_n$ACC_cfs, 3), '| PLSR = ', round(df_n$ACC_plsr, 3), 
                                                      '\nMAE: Raw =', round(df_n$AbErr_cfs, 3), '| PLSR = ', round(df_n$AbErr_plsr, 3)))
          ggsave(paste0('plsrQQplot_wk', wk_i, '_mon.', pred_mon_i, '_',var_i,'_multiVar_huc', hru_n, '_', paste0(preds, collapse = "."), '.png'),
                 g, height = 4, width = 6, dpi = 100)
          g <- gridExtra::grid.arrange(grobs = list(g1, g2, g3), layout_matrix = rbind(c(1,2), c(3, 3)),
                                       top = paste(wk_nm, 'wk Forecast of', month.abb[pred_mon_i], pred_var_nm, '-', hru_nm, '(', paste(var_ls[preds], collapse = ', '), ')'))
          ggsave(paste0('plsrTimeSeries_wk', wk_i, '_mon.', pred_mon_i, '_',var_i,'_multiVar_huc', hru_n, '_', paste0(preds, collapse = "."), '.png'),
                 g, height = 7, width = 11, dpi = 100)
          
          
        }
      }
    }
  }
}
