# ================================================================================ 
# Analyze climate indices forecasted for extreme precipatation events
#   Created by S. Baker 
# ================================================================================ 
rm(list=ls())

## --- Directories
TP_data = '/home/sabaker/s2s/analysis/files/cfsv2_files/2wk_avg/'
out_model_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/plsr_input/'
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/ploLoadings.R')


## --- Input Data
#beg.time         = Sys.time()
preds            = 1:10 # q2m, sst, prate
cut_domain       = T
cut_num          = 1
hru              = '1304'  #'1802' '1304'
wk               = '2_3' #'2_3'
pred_var         = 'prate'
pred_mon         = 1 #month to predict in CV
mon_nm           = 'DFJ'
# Ntime_sub        = 30       # Define number of timesteps that each sub-period (i.e. validation) will have
# Ncomp            = 2
# flagPress        = 0 #flag to indicate if the PRESS statistic is used (1) or not (0) to predict. If zero, only one component is used
if (cut_domain) {
  file_in = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.cutdomain', cut_num,'.rds')
  # file_out         = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'_cut.', paste(preds, collapse = '.'),'.rds')
} else {
  file_in = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.rds')
  # file_out         = paste0('CVresults_',pred_var, '_',hru,'_wk', wk, '_mon.', pred_mon,'.rds')
}


## --- libraries
library(lubridate)
library(ncdf4)
library(pls)
library(dplyr)


##### ===== Load CFSv2 Data 
setwd(out_model_dir)
load(file = file_in)
nms = names(data_out)

## extract objects
pred_cfs = data_out[[nms[1]]]
grid = data_out[[nms[2]]]
var_nms = c(data_out[[nms[3]]], 'hru_cfs')
date_cfs = data_out[[nms[4]]]
rm(data_out)


##### ===== Load NLDAS Data 
setwd(TP_data)
file_nld   = nc_open(paste0('nldas.1_2.', pred_var, '.nc')) 

# Extract hru coordinates for NLDAS
hru_v   = ncvar_get(file_nld, "hru")  # hru
hruInd = which(hru_v == hru)

## read time and find overlapping index
time = ncvar_get(file_nld, "time")
date_nld = as.POSIXct((time - 7) * 86400, 
                      origin = '1980-01-01', tz = 'utc') #timeseries of day 1 of obs
ind_overlap = inner_join(data.frame(date = date_cfs), 
                         data.frame(date = date_nld, ind = 1:length(date_nld)), by = 'date')

## read variable and convert units
nld_df   = ncvar_get(file_nld, pred_var)[hruInd, unique(ind_overlap$ind)]
if (pred_var == 'prate') { nld_df = nld_df * 24 } ## convert mm/hr to mm/d
nc_close(file_nld)

##### ===== Load CFSv2 T & P prediction on HRU
file_pred   = nc_open(paste0('cfsv2.', wk, '.', pred_var, '.nc')) 
time_pred = ncvar_get(file_pred, "time")
hru_v   = ncvar_get(file_pred, "hru")  # hru
hruInd = which(hru_v == hru)

## get date and mearge
date_pred = as.POSIXct(as.Date(as.POSIXct(time_pred - 6*86400, 
                                          origin = '1970-01-01', tz = 'utc'))) #timeseries of day 1 in cst
df_date = data.frame(date = date_pred, ind = 1:length(date_pred))
ind_overlap = left_join(ind_overlap, df_date, by = 'date')
cfs_pred_df   = ncvar_get(file_pred, pred_var)[hruInd, ind_overlap$ind.y]
if (pred_var == 'prate') { cfs_pred_df = cfs_pred_df * 86400 } ## convert kg/m^2/s to mm/d
cfs_pred_df = data.frame(date = ind_overlap$date, var = cfs_pred_df)
cfs_pred_df = cfs_pred_df %>% group_by(date) %>% summarise(var= mean(var)) # daily average

## Scale NLDAS variable
nldasVar      = scale(nld_df)
Ymean         = mean(nld_df)     # Mean nldasVar during training period
Ysd           = sd(nld_df)       # Standard deviation of nldasVar during training period


##### =====  get groups for validation?
Ntime           = length(pred_cfs[1,1,1,])
date_v          = unique(ind_overlap$date)
date_v          = data.frame(date = date_v, mon = month(date_v), ind = 1:length(date_v))
date_val        = filter(date_v, mon == pred_mon)
Nval            = nrow(date_val)
# Ngroups         = ceiling(Nval/Ntime_sub)
# cutgCroup        = sapply(1:Ngroups, function(j) c(rep(j,Ntime_sub)))[1:Nval] #not random
# xvalindex       = split(date_val$ind, cutgroup) # validate on these indexes
# Nsegments       = length(xvalindex)-1

##### ===== Arrange data 
## extract high prate events (scaled)
exceed_var = 2
mon_var = cbind.data.frame(var = nldasVar, ind = 1:length(nldasVar))[date_val$ind,]
var_ex = filter(mon_var, var > exceed_var)

## calculate month climatology for each variable
plot_ls <- NULL
for(ipred in preds){ # Start loop over selected predictors
  Ny = length(na.omit(unique(grid[ipred, 1,])))
  Nx = length(na.omit(unique(grid[ipred, 2,])))
  
  # t = which(!is.na(pred_cfs[ipred,1:Nx,1:Ny,1]), TRUE)
  # pred_i = pred_i[,t]
  
  Ymean         = mean(na.omit(as.vector(pred_cfs[ipred,1:Nx,1:Ny,])))
  Ysd           = sd(na.omit(as.vector(pred_cfs[ipred,1:Nx,1:Ny,])))
  
  # calc scaled climatology, extreme and anomaly for grid cell
  anom = array(0,dim=c(Nx*Ny))
  for(ilon in 1:Nx){
    for(ilat in 1:Ny){
      pred_clim  = (mean(pred_cfs[ipred,ilon,ilat,date_val$ind]) - Ymean ) / Ysd
      pred_ex  = (mean(pred_cfs[ipred,ilon,ilat,var_ex$ind]) - Ymean ) / Ysd
      anom[(Ny * (ilon - 1) + ilat)] = pred_ex - pred_clim
      #pred_i[,(Ny * (ilon - 1) + ilat)]  = pred_cfs[ipred,ilon,ilat,]
      
    }
  }
  
  ## plot
  # get grid and df
  grid_p = (grid[ipred, , ])
  plot_grid = na.omit(data.frame(x = grid_p[1,],
                                 y = grid_p[2,]))
  df_plot = na.omit(cbind.data.frame(var = anom, lat_p = plot_grid$x, 
                             lon_p = plot_grid$y))
  #set bounds
  max = max(df_plot$var, abs(min(df_plot$var)))
  max_r = round(max, 0.5)
  max = ifelse((max_r - max) > 0, max, max_r + 0.1)
  dif = (2*(max)/10)
  bin_var = seq(-max, max, dif)
  
  plot_ls[[ipred]] = plotLoadings(df_plot, var, #max = max, min = -max,
                                bin_vec = bin_var,
                                title_in = paste0('Anomalies of ', var_nms[ipred], 
                                                  ' for prate > ', exceed_var))  
}



gridExtra::grid.arrange(plot_ls[[1]]  + theme(legend.position = 'none'), 
                        plot_ls[[2]]  + theme(legend.position = 'none'), 
                        plot_ls[[3]]  + theme(legend.position = 'none'), 
                        plot_ls[[4]]  + theme(legend.position = 'none'), nrow = 2)
gridExtra::grid.arrange(plot_ls[[5]]  + theme(legend.position = 'none'), 
                        plot_ls[[6]]  + theme(legend.position = 'none'), 
                        plot_ls[[7]]  + theme(legend.position = 'none'), 
                        plot_ls[[8]] + theme(legend.position = 'none'), 
                        plot_ls[[9]]  + theme(legend.position = 'none'), 
                        plot_ls[[10]] + theme(legend.position = 'none'), nrow = 3)


gridExtra::grid.arrange(plot_ls[[1]], 
                        plot_ls[[2]], 
                        plot_ls[[3]], 
                        plot_ls[[4]], nrow = 2)
gridExtra::grid.arrange(plot_ls[[5]], 
                        plot_ls[[6]], 
                        plot_ls[[7]],
                        plot_ls[[8]], nrow = 2)
plot_ls[[9]]
plot_ls[[10]]
