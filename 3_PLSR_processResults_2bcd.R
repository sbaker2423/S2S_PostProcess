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
preds            = c(2,10) 
script           = '2b' #2b (base), 2c (stepwise), 2d (pc glm)
stepwise         = T

### === File Names
if (cut_domain & script ==  '2c') {
  setwd(model_stepwisePLSR_dir)
  file_mdlResults         = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.',cut_num,'_lag.',l,
                                   '_', paste(preds, collapse = '.'),'_stepwise.rds')
} else if (cut_domain & script ==  '2b') {
  setwd(model_plsrbase_dir)
  file_mdlResults         = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.',cut_num,'_lag.',l,
                                   '_', paste(preds, collapse = '.'),'_base.rds')
} else if (cut_domain & script ==  '2d') {
  setwd(model_pcGLM_dir)
  file_mdlResults         = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.',cut_num,'_lag.',l,
                                   '_', paste(preds, collapse = '.'),'_pcGLM.rds')
}

### ==== Load Data
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
test = rbind(
             cbind.data.frame(ind = df$ind, val = df$Y_cfs, fcst = 'CFSv2'),
             cbind.data.frame(ind = df$ind, val = df$Y_nld, fcst = 'NLDAS'),
             cbind.data.frame(ind = df$ind, val = df$Y_plsr, fcst = 'PLSR'))

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
# if (stepwise) {
#   ld = apply(ls_in$load[1:length(ls_in$load[,1,1]), ,], c(2,3), mean)
# } else {
  ld = apply(ls_in$load[1:length(ls_in$load[,1,1]), ,], c(2,3), mean)
# }



### === Plot gridded Loadings
## set bins
max = round(max(na.omit(ld), abs(min(na.omit(ld)))), -1)
dif = round(2*(max)/10)
bin_var = seq(-max, max, dif)
#bin_var = c(-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20)


### === plot loadings
plot_ls = list()
for (var_i in 1:length(preds)) {
  if(stepwise) {
    comp_i = na.omit(ld[,var_i])
  } else {
    comp = ld[,comp_plot]
    beg = ind_pred[var_i] + 1
    end = ind_pred[(var_i+1)]
    comp_i = comp[(beg:end)]
  }
  
  var = var_nms[var_i]
  var_long = long_nms[preds[var_i]]
  grid_p = (ls_in$grid[preds[var_i], , ])

  # if (var_i == 'sst') {
  #   plot_grid = data.frame(x = grid_sst[,1], y = grid_sst[,2])
  # } else {
  plot_grid = na.omit(data.frame(x = grid_p[1,],
                                 y = grid_p[2,]))
  # }
  
  df_plot = cbind.data.frame(var = comp_i, lat_p = plot_grid$x, 
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
                          plot_ls[[var_nms[3]]]  + theme(legend.position = 'none'), nrow = 3)
} else if (length(var_nms) == 4) {
  gridExtra::grid.arrange(plot_ls[[var_nms[1]]] + theme(legend.position = 'none'), 
                          plot_ls[[var_nms[2]]]  + theme(legend.position = 'none'), 
                          plot_ls[[var_nms[3]]]  + theme(legend.position = 'none'),
                          plot_ls[[var_nms[4]]]  + theme(legend.position = 'none'), nrow = 2)
}
