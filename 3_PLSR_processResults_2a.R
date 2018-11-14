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
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/ploLoadings.R')
source("/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/0_control_PLSRmdl.R")

library(dplyr)


### ==== Load Data
setwd(allVar_model_dir)
ls_in = readRDS(file_allVars)
nm = names(ls_in)
ind_pred = ls_in$ind_pred
var_nms = ls_in$vars_nms
df = ls_in$df

## summarize statistics
df_grp = group_by(df, pred) 
# dplyr::summarize(df_grp, c = cor(Y_cfs,Y_nld))
t = cbind.data.frame(cor.cfs = dplyr::summarize(df_grp, cor(Y_cfs,Y_nld)),
                     cor.plsr = dplyr::summarize(df_grp, cor(Y_plsr,Y_nld))[,2],
                     res.cfs = dplyr::summarize(df_grp, mean(abs(res_cfs)))[,2],
                     res.plsr = dplyr::summarize(df_grp, mean(abs(res_plsr)))[,2])
colnames(t) <- c('pred', 'cor.cfs', 'cor.plsr', 'res.cfs', 'res.plsr')
round(t, digits = 3)

### ==== 1:1 plot
par(mfrow=c(3,4))
ymax = max(df$Y_nld, df$Y_cfs, df$Y_plsr)
for (i in 1:length(var_nms)) {
  p = filter(df, pred == i)
  plot(p$Y_nld, p$Y_plsr, xlim = c(0,ymax), ylim = c(0,ymax), ylab = var_nms[i]); abline(0,1) 
}
plot(p$Y_nld, p$Y_cfs, xlim = c(0,ymax), ylim = c(0,ymax), ylab = 'Raw'); abline(0,1) 


### ==== plot timeseries
g <- list()
test = df[,c(1,3:6)]
test2 = reshape2::melt(test, id.vars = c('pred', 'ind'), na.rm = T)
for (i in 1:length(var_nms)) {
  dplot = filter(test2, pred == i)
  dplot$variable = factor(dplot$variable, c('Y_cfs', 'Y_nld', 'Y_plsr'))
  g[[var_nms[i]]] = ggplot(dplot, aes(x = ind, y = value, group = variable, colour = variable)) +
    geom_line() +
    # xlab('Index of January Day (1999 - 2010)') +
    xlab('Day Index') +
    ylab(paste0('Precip - ', var_nms[i])) +
    # ylab('Precipitation (mm/day per bi-weekly period)') +
    scale_color_manual(values = c("light blue", "black", 'red')) +
    scale_linetype_manual(values = c('dotted', 'solid', 'solid')) +
    theme(legend.position = 'none')
}

## plot
gridExtra::grid.arrange(g[['hgt']] + theme(legend.position = 'none'), 
                        g[['q2m']]  + theme(legend.position = 'none'), 
                        g[['prs']]  + theme(legend.position = 'none'), 
                        g[['slp']]  + theme(legend.position = 'none'), 
                        g[['pwt']] + theme(legend.position = 'none'), 
                        g[['sst']]  + theme(legend.position = 'none'), 
                        g[['uwnd']]  + theme(legend.position = 'none'), 
                        g[['vwnd']] + theme(legend.position = 'none'), 
                        g[['olr']] + theme(legend.position = 'none'),
                        g[['tmp']] + theme(legend.position = 'none'),
                        g[['prt']] + theme(legend.position = 'none'),nrow = 6)


# gridExtra::grid.arrange(g[['hgt']] + theme(legend.position = 'none'), 
#                         g[['q2m']]  + theme(legend.position = 'none'), 
#                         g[['prs']]  + theme(legend.position = 'none'), 
#                         g[['slp']]  + theme(legend.position = 'none'), nrow = 2)
# gridExtra::grid.arrange(g[['pwt']] + theme(legend.position = 'none'), 
#                         g[['sst']]  + theme(legend.position = 'none'), 
#                         g[['uwnd']]  + theme(legend.position = 'none'), 
#                         g[['vwnd']] + theme(legend.position = 'none'), nrow = 2)
# gridExtra::grid.arrange(g[['tmp']] + theme(legend.position = 'none'),
#                         g[['prt']] + theme(legend.position = 'none'),nrow = 2)



### ==== Plot loadings
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
                        plot_ls[['slp']]  + theme(legend.position = 'none'), 
                        plot_ls[['pwt']] + theme(legend.position = 'none'), 
                        plot_ls[['sst']]  + theme(legend.position = 'none'), 
                        plot_ls[['uwnd']]  + theme(legend.position = 'none'), 
                        plot_ls[['vwnd']] + theme(legend.position = 'none'), 
                        plot_ls[['olr']] + theme(legend.position = 'none'),
                        plot_ls[['tmp']] + theme(legend.position = 'none'),
                        plot_ls[['prt']] + theme(legend.position = 'none'),nrow = 6)



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