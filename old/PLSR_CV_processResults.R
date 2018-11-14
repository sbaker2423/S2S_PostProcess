# ================================================================================ 
# PLSR code for processing CFSv2 variables for predicting 
# -- Process CV output
#   Created by S. Baker 
#
# ================================================================================ 
rm(list=ls())

source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/ploLoadings.R')
source("/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/0_control_PLSRmdl.R")

library(dplyr)

## --- Input Data
preds            = c(5,8) 

### === File Names
if (cut_domain) {
  file_mdlResults         = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.',cut_num,'_lag.',l,
                                   '_', paste(preds, collapse = '.'),'.rds')
} else {
  file_mdlResults         = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'.rds')
}

# ## --- Directories
# TP_data = '/home/sabaker/s2s/analysis/files/cfsv2_files/2wk_avg/'
# out_model_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/plsr_input/'
# source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/ploLoadings.R')
# 
# 
# ## --- Input Data
# preds            = c(5,8) # q2m, sst, prate
# cut_domain       = T
# hru              ='1304'  #'1802' '1304'
# wk               = '2_3' #'2_3'
# pred_var         = 'prate'
# pred_mon         = 1 #month to predict in CV
# # mon_nm           = 'DFJ
# # Ncomp            = 2
# if (cut_domain) {
#   file_out = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.', paste(preds, collapse = '.'),'.rds')
# } else {
#   file_out = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'.rds')
# }
# long_nms = c('500 mb Geopotential height', 'Specific Humid 2m', 'Surface Pressure', 
#              'Sea Level Pressure', 'Precipitable Water', 'Zonal Winds (850 mb)',
#              'Meridional Winds (850 mb)', 'Sea Surface Temperature', 
#              'Surface Temperature 2m', 'Surface Prate')

### ==== Load Data
setwd(out_model_dir)
ls_in = readRDS(file_mdlResults)
nm = names(ls_in)
ind_pred = ls_in$ind_pred
var_nms = ls_in$vars_nms
df = ls_in$df

## calculate residual and correlation for cfsv2 and model
mean(abs(df$res_cfs))
mean(abs(df$res_plsr))
cor(df$Y_nld, df$Y_cfs)
cor(df$Y_nld, df$Y_plsr)

### ==== 1:1 plot
par(mfrow=c(1,2))
ymax = max(df$Y_nld, df$Y_cfs, df$Y_plsr)
plot(df$Y_nld, df$Y_plsr, xlim = c(0,ymax), ylim = c(0,ymax)); abline(0,1) 
plot(df$Y_nld, df$Y_cfs, xlim = c(0,ymax), ylim = c(0,ymax)); abline(0,1) 

### ==== plot timeseries
library(tidyr)
# test = df[,c(1:5)]
# test2 = tidyr::gather(test, fcst, c(Y_plsr, Y_nld, Y_cfs))
test = rbind(
             cbind.data.frame(df[,c(1,2)], val = df[,4], fcst = 'CFSv2'),
             cbind.data.frame(df[,c(1,2)], val = df[,5], fcst = 'NLDAS'),
             cbind.data.frame(df[,c(1,2)], val = df[,3], fcst = 'PLSR'))
ggplot(test, aes(x = ind, y = val, group = fcst, colour = fcst)) +
  geom_line() +
  xlab('Index of January Day (1999 - 2010)') +
  ylab('Precipitation (mm/day per bi-weekly period)') +
  scale_color_manual(values = c("light blue", "black", 'red')) +
  scale_linetype_manual(values = c('dotted', 'solid', 'solid'))

### ==== input
comp_plot  = 1     # EOF
# igroup = 1
# ld = ls_in$load[igroup,,]
ld = apply(ls_in$load[1:length(ls_in$load[,1,1]), ,], c(2,3), mean)


### === Plot gridded Loadings
## set bins
max = round(max(ld, abs(min(ld))), -1)
dif = round(2*(max)/10)
bin_var = seq(-max, max, dif)
#bin_var = c(-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20)


### === plot loadings
plot_ls = list()
for (var_i in 1:length(var_nms)) {
  comp = ld[,comp_plot]
  var = var_nms[var_i]
  var_long = long_nms[preds[var_i]]
  grid_p = (ls_in$grid[preds[var_i], , ])
  beg = ind_pred[var_i] + 1
  end = ind_pred[(var_i+1)]
  # if (var_i == 'sst') {
  #   plot_grid = data.frame(x = grid_sst[,1], y = grid_sst[,2])
  # } else {
  plot_grid = na.omit(data.frame(x = grid_p[1,],
                                 y = grid_p[2,]))
  # }
  
  df_plot = cbind.data.frame(var = comp[(beg:end)], lat_p = plot_grid$x, 
                             lon_p = plot_grid$y)
  
  
  plot_ls[[var]] = plotLoadings(df_plot, var, #max = max, min = -max,
                                  bin_vec = bin_var,
                                  title_in = paste0(var_long, ' Mean Loadings - Component ', comp_plot))  
}

## plot
if (length(var_nms) == 2 ) {
  gridExtra::grid.arrange(plot_ls[[var_nms[1]]] + theme(legend.position = 'none'), 
                          plot_ls[[var_nms[2]]]  + theme(legend.position = 'none'), nrow = 2)
} else if (length(var_nms) == 3 ) {
  gridExtra::grid.arrange(plot_ls[[var_nms[1]]] + theme(legend.position = 'none'), 
                          plot_ls[[var_nms[2]]]  + theme(legend.position = 'none'), 
                          plot_ls[[var_nms[3]]]  + theme(legend.position = 'none'), nrow = 2)
} else if (ength(var_nms) == 4) {
  gridExtra::grid.arrange(plot_ls[[var_nms[1]]] + theme(legend.position = 'none'), 
                          plot_ls[[var_nms[2]]]  + theme(legend.position = 'none'), 
                          plot_ls[[var_nms[3]]]  + theme(legend.position = 'none'),
                          plot_ls[[var_nms[4]]]  + theme(legend.position = 'none'), nrow = 2)
}


# gridExtra::grid.arrange(plot_ls[['hgt']] + theme(legend.position = 'none'), 
#                         plot_ls[['q2m']]  + theme(legend.position = 'none'), 
#                         plot_ls[['prs']]  + theme(legend.position = 'none'), 
#                         plot_ls[['slp']]  + theme(legend.position = 'none'), nrow = 2)
# gridExtra::grid.arrange(plot_ls[['pwt']] + theme(legend.position = 'none'), 
#                         plot_ls[['sst']]  + theme(legend.position = 'none'), 
#                         plot_ls[['uwnd']]  + theme(legend.position = 'none'), 
#                         plot_ls[['vwnd']] + theme(legend.position = 'none'), nrow = 2)
# plot_ls[['tmp']]
# plot_ls[['prt']]

## summary
summary(plsrmod)
plot(RMSEP(plsrmod))
plot(plsrmod, ncomp = 3, asp = 1, line = TRUE) # because scaled?
plot(plsrmod, plottype = "scores", comps = 1:3)
explvar(plsrmod)
