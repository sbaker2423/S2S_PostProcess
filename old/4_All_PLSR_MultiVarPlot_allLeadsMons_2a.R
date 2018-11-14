# ================================================================================ 
# Step 3: process PLSR CV results
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
library(dplyr)
library(grid)
library(gridExtra)

## Loop through wk lead, predictor, month
pred_base = c(8) # uses two preds max not including T/P
save_plots = T

### INPUT HRU to look at loadings
# hru_input   = c('1802', '1304', '1306', '1401', '1402', '1709', '303', '1014')
# output_hru = F

## Files
df_allMon_file = paste0('summMon_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_allSeas_file = paste0('summSeas_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_statMon_file = paste0('summStatsMon_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')
df_statSeas_file = paste0('summStatsSeas_CVresults_allLeadsMonsVars_predBase.', paste0(pred_base, collapse = "."), '.rds')


##  === Read data from files
## Month - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_allMon = readRDS(file = df_allMon_file)

## Season - read df summarizing cor at every lead/mon/var
setwd(allHRU_multiVar_dir)
df_allSeas = readRDS(file = df_allSeas_file)


#####   ===========================================================   ######
##                      Plot MONTHLY Correlation Figs                     ##
#####   ===========================================================   ######


## plot and save
if (save_plots) {
  setwd(plot_dir)
  
  ##  === BINS - Increase Plot
  maxCor_in = 0.6; bin_size = 0.05
  df_allMon$bins_inc = cut(df_allMon$inc_cor, breaks = seq(from = 0, to = maxCor_in, by = bin_size)) 
  df_allMon$bins_inc <- addNA(df_allMon$bins_inc)
  levels(df_allMon$bins_inc)[is.na(levels(df_allMon$bins_inc))] <- "None"
  ## color palette and plot
  # n_col = length(levels(df_allMon$bins_inc)) -1 #number of colors needed
  pal = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
  color_fun <- colorRampPalette(pal, space = "rgb") # for discrete
  cust_pal_inc = c(color_fun(length(levels(df_allMon$bins_inc)) -1 ), 'grey85')  # add none color
  lab_inc = gsub(',',' - ', levels(df_allMon$bins_inc)) # replace legend labels
  
  ##  === BINS - Highest and Raw Correlation Plots
  bin_set = c(-0.2, seq(0, 0.5, by = 0.05), 0.6, 0.7, 0.8, 1)
  df_allMon$bins_base = cut(df_allMon$val_max, breaks = bin_set)
  n_col = length(levels(df_allMon$bins_base))
  cust_pal_base = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
  lab_base = gsub(',',' - ', levels(df_allMon$bins_base)) # replace legend labels
  
  ## Loop through to collect plots / summary statistics
  for (var_i in c('prate', 'tmp2m')) {
    for (wk_i in c('2_3', '3_4')) {
      
      ## Filter and plot name
      df_i = filter(df_allMon, wk == wk_i, pred_var == var_i)
      df_i$mon_pred = factor(month.abb[df_i$mon_pred], levels = month.abb)
      
      ## Plot raw CFSv2 correlation
      df_i$bins = cut(df_i$cor_cfs, breaks = bin_set)
      rawCor = plot_map(var_i, df_i, 'cor_cfs', type = 'discrete', lab = lab_base, 
                        cust_pal = cust_pal_base, legend_title = 'Correlation')  +
        theme(legend.position="none")
      
      ## Plot max correlation
      df_i$bins = df_i$bins_base
      maxCor = plot_map(var_i, df_i, 'val_max', type = 'discrete', lab = lab_base, 
                        cust_pal = cust_pal_base, legend_title = 'Correlation') +
        theme(legend.position="none")
      
      ## Plot increase in correlation
      df_i$bins = df_i$bins_inc
      corInc = plot_map(var_i, df_i, 'inc_cor', type = 'discrete', lab = lab_inc, 
                        cust_pal = cust_pal_inc, legend_title = 'Cor. Increase') +
        theme(legend.position="none")
      
      # trying to add title...not working
      # corInc = ggplotGrob(corInc)
      # library(gtable)
      # z <- gtable_add_rows(corInc, unit(2, 'cm'), 3)
      # z <- gtable_add_grob( z, 
      #                                                 list(rectGrob(gp = gpar(col = NA, fill = gray(0.9))),
      #                                                      textGrob("Inc. in Correlation", gp = gpar(col = gray(0)))),
      #                                                 2, 2, 3, 6, name = paste(runif(2)))

      ## print and save figures
      corInc = corInc + facet_grid(mon_pred ~ .)
      ggsave(paste0('incCor_plsrResults_Mon_wk', wk_i, '_', var_i,'_multiVar_', paste0(pred_base, collapse = "."), '.png'),
             corInc, height = 16, width = 4, dpi = 100)
      rawCor = rawCor + facet_grid(mon_pred ~ .)
      ggsave(paste0('rawCor_plsrResults_Mon_wk', wk_i, '_', var_i,'_multiVar_', paste0(pred_base, collapse = "."), '.png'),
             rawCor, height = 16, width = 4, dpi = 100)
      maxCor = maxCor + facet_grid(mon_pred ~ .)
      ggsave(paste0('maxCor_plsrResults_Mon_wk', wk_i, '_', var_i,'_multiVar_', paste0(pred_base, collapse = "."), '.png'),
             maxCor, height = 16, width = 4, dpi = 100)
      
    }
  }
}


#####   ===========================================================   ######
##                      Plot SEASON Correlation Figs                      ##
#####   ===========================================================   ######


## plot and save
if (save_plots) {
  setwd(plot_dir)
  
  ##  === BINS - Increase Plot
  maxCor_in = 0.6; bin_size = 0.05
  df_allSeas$bins_inc = cut(df_allSeas$inc_cor, breaks = seq(from = 0, to = maxCor_in, by = bin_size)) 
  df_allSeas$bins_inc <- addNA(df_allSeas$bins_inc)
  levels(df_allSeas$bins_inc)[is.na(levels(df_allSeas$bins_inc))] <- "None"
  ## color palette and plot
  # n_col = length(levels(df_allSeas$bins_inc)) -1 #number of colors needed
  pal = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
  color_fun <- colorRampPalette(pal, space = "rgb") # for discrete
  cust_pal_inc = c(color_fun(length(levels(df_allSeas$bins_inc)) -1 ), 'grey85')  # add none color
  lab_inc = gsub(',',' - ', levels(df_allSeas$bins_inc)) # replace legend labels
  
  ##  === BINS - Highest and Raw Correlation Plots
  bin_set = c(-0.2, seq(0, 0.5, by = 0.05), 0.6, 0.7, 0.8, 1)
  df_allSeas$bins_base = cut(df_allSeas$val_max, breaks = bin_set)
  n_col = length(levels(df_allSeas$bins_base))
  cust_pal_base = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
  lab_base = gsub(',',' - ', levels(df_allSeas$bins_base)) # replace legend labels
  
  ## Loop through to collect plots / summary statistics
  for (var_i in c('prate', 'tmp2m')) {
    for (wk_i in c('2_3', '3_4')) {
      
      ## Filter and plot name
      df_i = filter(df_allSeas, wk == wk_i, pred_var == var_i)
      df_i$seas = factor(df_i$seas, levels = c('DJF', 'MAM', 'JJA', 'SON'))
      
      ## Plot raw CFSv2 correlation
      df_i$bins = cut(df_i$cor_cfs, breaks = bin_set)
      rawCor = plot_map(var_i, df_i, 'cor_cfs', type = 'discrete', lab = lab_base, 
                        cust_pal = cust_pal_base, legend_title = 'Correlation')  +
        theme(legend.position="none")
      
      ## Plot max correlation
      df_i$bins = df_i$bins_base
      maxCor = plot_map(var_i, df_i, 'val_max', type = 'discrete', lab = lab_base, 
                        cust_pal = cust_pal_base, legend_title = 'Correlation') +
        theme(legend.position="none")
      
      ## Plot increase in correlation
      df_i$bins = df_i$bins_inc
      corInc = plot_map(var_i, df_i, 'inc_cor', type = 'discrete', lab = lab_inc, 
                        cust_pal = cust_pal_inc, legend_title = 'Cor. Increase') +
        theme(legend.position="none")
      
      ## print and save figures
      corInc = corInc + facet_grid(seas ~ .)
      ggsave(paste0('incCor_plsrResults_Seas_wk', wk_i, '_', var_i,'_multiVar_', paste0(pred_base, collapse = "."), '.png'),
             corInc, height = 7, width = 4, dpi = 200)
      rawCor = rawCor + facet_grid(seas ~ .)
      ggsave(paste0('rawCor_plsrResults_Seas_wk', wk_i, '_', var_i,'_multiVar_', paste0(pred_base, collapse = "."), '.png'),
             rawCor, height = 7, width = 4, dpi = 200)
      maxCor = maxCor + facet_grid(seas ~ .)
      ggsave(paste0('maxCor_plsrResults_Seas_wk', wk_i, '_', var_i,'_multiVar_', paste0(pred_base, collapse = "."), '.png'),
             maxCor, height = 7, width = 4, dpi = 200)
      
    }
  }
}


#####   ===========================================================   ######
##                      OLDDDDDDDDDDDDDDDDDDDDDDDDDD                      ##
#####   ===========================================================   ######




# pred_var_nm = ifelse(pred_var == 'prate', 'Precipitation', 'Temperature')
# wk_nm = ifelse(wk == '2_3', '2-3', ifelse(wk == '3_4', '3-4', '1-2'))
# 
# ## save plot
# title = paste(month.abb[pred_mon], 'Week', gsub('_', '-', wk),  
#               pred_var_nm, 'Forecast - Predictors:', paste(var_ls[preds], collapse = ', '))
# g <-gridExtra::grid.arrange(p3,p2,p1,  nrow = 2, top = textGrob(title, gp=gpar(fontsize=16,font=8)))
# ggsave(paste0('plsrResults_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_', paste0(preds, collapse = "."), '.png'),
#        g, height = 8, width = 11, dpi = 100)
# 
# ## only increase figure
# title = paste('Predictors:', paste(var_ls[preds], collapse = ', '))
# g <-gridExtra::grid.arrange(p1 + theme(plot.title = element_blank()),  nrow = 1, top = textGrob(title, gp=gpar(fontsize=15,font=8)))
# ggsave(paste0('plsrResults_corIncrease_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_', paste0(preds, collapse = "."), '.png'), 
#        g, height = 7/2.25, width = 11/2, dpi = 100)
# 
# 
# ### ==== Load PLSR input data / grid data
# pred_var = 'prate'; wk = '2_3'; pred_mon = 1
# preds = c(8,4)
# 
# setwd(allHRU_multiVar_dir)
# file_input_multiVar    = paste0('Input_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'_multiVars_',
#                                 paste0(preds, collapse = "."), '.rds')
# f.out = readRDS(file = file_input_multiVar)
# 
# ### === ONE HRU: Loop through to plot loadings and timeseries ==== ###
# setwd(paste0(plot_dir, '/huc_individual_results'))
# if (output_hru) {
#   for (hru_n in hru_input) {
#     hru_nm = hru_tb[which(hru_tb[,1] == hru_n),2]
#     
#     ## set bins
#     max = 35 #round(max(ld, abs(min(ld))), -1)
#     dif = round(2*(max)/10)
#     bin_var = seq(-max, max, dif)
#     
#     ### === Plot Loadings
#     plot_ls = list(); n = 0
#     for (i in 1:length(preds)) {
#       for (comp_i in 1:Ncomp) {
#         
#         ## variable information
#         n = n + 1
#         ipred = preds[i]
#         var = var_ls[ipred]
#         var_long = long_nms[ipred]
#         
#         
#         
#         ## extract grid
#         grid_p = (f.out$grid[ipred, , ])
#         plot_grid = na.omit(data.frame(x = grid_p[1,], y = grid_p[2,]))
#         
#         ## loadings - loading for igroup vs. average loading of all igroups
#         igroup = 1
#         ind = c((f.out$ind_pred[i]+1):f.out$ind_pred[(i+1)])
#         ld = na.omit(data_all[[hru_n]]$load[igroup, ind, comp_i])
#         ld_avg = na.omit(apply(data_all[[hru_n]]$load[1:length(data_all[[hru_n]]$load[,1,1]), ind, comp_i], 2, mean))
#         
#         df_plot = cbind.data.frame(var = ld_avg, # ld_avg
#                                    lat_p = plot_grid$x, lon_p = plot_grid$y)
#         
#         plot_ls[[n]] =
#           plotLoadings(df_plot, var, #max = max, min = -max,
#                        bin_vec = bin_var, title_in = paste0(var_long, ' - Comp. ', comp_i))  +
#           theme(legend.position = 'none')
#       }
#     }
#     
#     ## plot & save
#     title = paste('Mean Loadings -', wk_nm, 'wk Forecast of', month.abb[pred_mon], pred_var_nm, '-', hru_nm)
#     g <- gridExtra::grid.arrange(grobs = plot_ls, nrow = length(preds),top = textGrob(title, gp=gpar(fontsize=15,font=8)))
#     ggsave(paste0('plsrLoadings_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_huc', hru_n, '_', paste0(preds, collapse = "."), '.png'),
#            g, height = 6, width = 11, dpi = 100)
#     
#     
#     ### ==== One HRU: 1:1 plots, timeseries ==== ###
#     
#     ## extract data and convert to numeric
#     df_i = as.data.frame(data_all[[hru_n]]$df)
#     if (pred_var == 'tmp2m') { df_i[,3:5] = df_i[,3:5] - 273.15 }
#     paste0(var_ls[preds], collapse = ', ')
#     ymax = max(df_i$Y_nld, df_i$Y_cfs, df_i$Y_plsr)
#     ymin = ifelse(pred_var == 'tmp2m', min(df_i$Y_nld, df_i$Y_cfs, df_i$Y_plsr), 0)
#     
#     ## plot two 1:1 plots
#     g1 <- ggplot(df_i, aes(Y_nld, Y_plsr)) + geom_abline(color = 'red') + geom_point(size = 0.5) + coord_equal() +
#       xlim(ymin,ymax) + ylim(ymin,ymax) + ylab(paste('PLSR', pred_var_nm)) + xlab(paste('NLDAS', pred_var_nm))
#     g2 <- ggplot(df_i, aes(Y_nld, Y_cfs)) + geom_abline(color = 'red') + geom_point(size = 0.5) + coord_equal() +
#       xlim(ymin,ymax) + ylim(ymin,ymax) + ylab(paste('CFSv2', pred_var_nm)) + xlab(paste('NLDAS', pred_var_nm))
#     
#     ## plot timeseries
#     dplot = reshape2::melt(df_i[,c(1:5)], id.vars = c('date','ind'), na.rm = T)
#     dplot$variable = factor(dplot$variable, c('Y_cfs', 'Y_nld', 'Y_plsr'))
#     g3 = ggplot(dplot, aes(x = ind, y = value, group = variable, colour = variable)) +
#       geom_line() +
#       # xlab('Index of January Day (1999 - 2010)') +
#       xlab('Day Index') +
#       ylab(paste0(pred_var_nm)) +
#       # ylab('Precipitation (mm/day per bi-weekly period)') +
#       scale_color_manual(values = c("light blue", "black", 'red')) +
#       scale_linetype_manual(values = c('dotted', 'solid', 'solid'))
#     
#     ## export and save
#     g <- gridExtra::grid.arrange(grobs = list(g1, g2, g3), layout_matrix = rbind(c(1,2), c(3, 3)),
#                                  top = paste(wk_nm, 'wk Forecast of', month.abb[pred_mon], pred_var_nm, '-', hru_nm, '(', paste(var_ls[preds], collapse = ', '), ')'))
#     ggsave(paste0('plsrTimeSeries_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_huc', hru_n, '_', paste0(preds, collapse = "."), '.png'),
#            g, height = 7, width = 11, dpi = 100)
#     
#     
#   }
# }
# #     }
# #   }
# # }