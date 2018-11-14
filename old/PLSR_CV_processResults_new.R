# ================================================================================ 
# Step 3: process PLSR CV results
#
# PLSR code for processing CFSv2 variables for predicting 
# -- Process CV output
#   Created by S. Baker 
#
# ================================================================================ 
rm(list=ls())

## --- Directories
TP_data = '/home/sabaker/s2s/analysis/files/cfsv2_files/2wk_avg/'
out_model_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/plsr_input/'
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/ploLoadings.R')

library(dplyr)

## --- Input Data
preds            = c(5,8) # q2m, sst, prate
cut_domain       = T
cut_num          = 4
l                = 1           # days of lag. If 0, no lag
hru              ='1304'  #'1802' '1304'
wk               = '2_3' #'2_3'
pred_var         = 'prate'
pred_mon         = 1 #month to predict in CV
# mon_nm           = 'DFJ
# Ncomp            = 2

if (cut_domain) {
  #file_out = paste0('CVresultsAllvars_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'.rds')
  file_out         = paste0('CVresultsAllvars_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',l,'.rds')
} else {
  file_out = paste0('CVresultsAllvars_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'.rds')
}
long_nms = c('500 mb Geopotential height', 'Specific Humid 2m', 'Surface Pressure', 
             'Sea Level Pressure', 'Precipitable Water', 'Zonal Winds (850 mb)',
             'Meridional Winds (850 mb)', 'Sea Surface Temperature', 
             'Surface Temperature 2m', 'Surface Prate')

### ==== Load Data
setwd(out_model_dir)
ls_in = readRDS(file_out)
nm = names(ls_in)
ind_pred = ls_in$ind_pred
var_nms = ls_in$vars_nms
df = ls_in$df

## calculate residual and correlation for cfsv2 and model
df_grp = group_by(df, pred)
t = cbind.data.frame(cor.cfs = summarize(df_grp, cor(Y_cfs,Y_nld)),
      cor.plsr = summarize(df_grp, cor(Y_plsr,Y_nld)),
      res.cfs = summarize(df_grp, mean(abs(res_cfs))),
      res.plsr = summarize(df_grp, mean(abs(res_plsr))))
# mean(abs(df$res_cfs))
# mean(abs(df$res_plsr))
# cor(df$Y_nld, df$Y_cfs)
# cor(df$Y_nld, df$Y_plsr)

### ==== 1:1 plot
par(mfrow=c(3,4))
ymax = max(df$Y_nld, df$Y_cfs, df$Y_plsr)
for (i in 1:length(var_nms)) {
  p = filter(df, pred == i)
  plot(p$Y_nld, p$Y_plsr, xlim = c(0,ymax), ylim = c(0,ymax)); abline(0,1) 
}
plot(p$Y_nld, p$Y_cfs, xlim = c(0,ymax), ylim = c(0,ymax)); abline(0,1) 

### ==== plot timeseries
# library(tidyr)
# test = df[,c(1:5)]
# test2 = tidyr::gather(test, fcst, c(Y_plsr, Y_nld, Y_cfs))
# test = rbind(
#              cbind.data.frame(df[1:(nrow(df)/length(var_nms)),c(2,3)], val = df[,5], fcst = 'CFSv2'),
#              cbind.data.frame(df[1:(nrow(df)/length(var_nms)),c(2,3)], val = df[,6], fcst = 'NLDAS'),
#              cbind.data.frame(df[,c(2,3)], val = df[,3], fcst = df[,1]))
# ggplot(test, aes(x = ind, y = val, group = fcst, colour = fcst)) +
#   geom_line() +
#   xlab('Index of January Day (1999 - 2010)') +
#   ylab('Precipitation (mm/day per bi-weekly period)') +
#   scale_color_manual(values=c("light blue", "black", 'red')) +
#   scale_linetype_manual(value = c('dotted', 'solid', 'solid'))

### ==== input
comp_plot  = 1     # EOF
igroup = 1
ipred = 1
ld = na.omit(ls_in$load[ipred,igroup,,])
ld_avg = na.omit(apply(ls_in$load[ipred,1:length(ls_in$load[1,,1,1]), ,], 2, mean))


### === Plot gridded Loadings
## set bins
max = round(max(ld, abs(min(ld))), -1)
dif = round(2*(max)/10)
bin_var = seq(-max, max, dif)
#bin_var = c(-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20)


### === plot loadings
plot_ls = list()
for (ipred in 1:length(var_nms)) {
  ld = na.omit(ls_in$load[ipred,igroup,,])
  ld_avg = na.omit(apply(ls_in$load[
    ipred,1:length(ls_in$load[1,,1,1]), ,], 2, mean))
  
  
  #comp = ld[,comp_plot]
  var = var_nms[ipred]
  var_long = long_nms[ipred]
  grid_p = (ls_in$grid[ipred, , ])
  #beg = ind_pred[ipred] + 1
  #end = ind_pred[(ipred+1)]
  # if (ipred == 'sst') {
  #   plot_grid = data.frame(x = grid_sst[,1], y = grid_sst[,2])
  # } else {
  plot_grid = na.omit(data.frame(x = grid_p[1,],
                                 y = grid_p[2,]))
  # }
  
  df_plot = cbind.data.frame(var = ld_avg, lat_p = plot_grid$x, 
                             lon_p = plot_grid$y)
  
  
  plot_ls[[var]] = plotLoadings(df_plot, var, #max = max, min = -max,
                                  bin_vec = bin_var,
                                  title_in = paste0(var_long, ' Mean Loadings - Component ', comp_plot))  
}

## plot
gridExtra::grid.arrange(plot_ls[['hgt']] + theme(legend.position = 'none'), 
                        plot_ls[['q2m']]  + theme(legend.position = 'none'), 
                        plot_ls[['prs']]  + theme(legend.position = 'none'), 
                        plot_ls[['slp']]  + theme(legend.position = 'none'), nrow = 2)
gridExtra::grid.arrange(plot_ls[['pwt']] + theme(legend.position = 'none'), 
                        plot_ls[['sst']]  + theme(legend.position = 'none'), 
                        plot_ls[['uwnd']]  + theme(legend.position = 'none'), 
                        plot_ls[['vwnd']] + theme(legend.position = 'none'), nrow = 2)
gridExtra::grid.arrange(plot_ls[['tmp']] + theme(legend.position = 'none'),
                        plot_ls[['prt']] + theme(legend.position = 'none'),nrow = 2)
plot_ls[['prt']]

## summary
summary(plsrmod)
plot(RMSEP(plsrmod))
plot(plsrmod, ncomp = 3, asp = 1, line = TRUE) # because scaled?
plot(plsrmod, plottype = "scores", comps = 1:3)
explvar(plsrmod)
