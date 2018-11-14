# ================================================================================ 
# Step 1: process CFSv2 variables to correct format
#         - read from netCDFs, regrid, extract desired time
# PLSR code for processing CFSv2 variables for predicting 
# NLDAS temperature and precipitation for 202 HUCs in CONUS domain
#   Created by S. Baker 
#
# Number of CFSv2 predictor variables
# (1)  Geopotential height (z500_f)
# (2)  Specific Humid 2m (q2m_f)
# (3)  Surface Pressure (pressfc_f)
# (4)  Sea Level Pressure (meanSLP)
# (5)  precipitable water (pwat_f)
# (6)  Zonal winds (U_wnd850_f)
# (7)  Meridional winds (V_wnd850_f)
# (8)  Sea Surface Temperature (ocnsst_h)
# (9)  Outgoing Longwave Radiation (..._f)
# (10) Surface Temperature (tmp2m_f)
# (11) Surface Prate (prate_f)
#
#   Adapted from P. Mendoza script (May 2016) and Andy Wood, NCAR
#   from PLSR R code (A. Wood) for cross-validated application of PCR/PLSR
# ================================================================================ 
rm(list=ls())

# ### === Input Data
source("/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/0_control_PLSRmdl.R")

### === libraries
library(lubridate)
library(ncdf4)
library(abind)
library(OceanView)
library(pls)
library(dplyr)


##### ===== Load Predictors & Predictanad ===== #####

## Predictors - Open netCDF files
setwd(cfs_dir)
file_hgt   = nc_open(paste0('cfsv2.', wk, '.z500_f.nc')) 
file_q2m   = nc_open(paste0('cfsv2.', wk, '.q2m_f.nc'))
file_prs   = nc_open(paste0('cfsv2.', wk, '.pressfc_f.nc'))
file_slp   = nc_open(paste0('cfsv2.', wk, '.prmsl_f.nc'))
file_pwt   = nc_open(paste0('cfsv2.', wk, '.pwat_f.nc'))
file_uwnd  = nc_open(paste0('cfsv2.', wk, '.U_wnd850_f.nc'))
file_vwnd  = nc_open(paste0('cfsv2.', wk, '.V_wnd850_f.nc'))
file_sst   = nc_open(paste0('cfsv2.', wk, '.ocnsst_h.nc'))
file_olr   = nc_open(paste0('cfsv2.', wk, '.ulwtoa_f.nc'))
file_tmp   = nc_open(paste0('cfsv2.', wk, '.tmp2m_f.nc'))
file_prt   = nc_open(paste0('cfsv2.', wk, '.prate_f.nc'))


##### ===== Get Lat & Lon Index ===== #####

## Extract lat/lon coordinates - different for variables
lat_in <- lon_in <- list()
lat_in[[var_ls[1]]]  = ncvar_get(file_hgt, "latitude"); lon_in[[var_ls[1]]]  = ncvar_get(file_hgt, "longitude") 
lat_in[[var_ls[2]]]  = ncvar_get(file_q2m, "latitude"); lon_in[[var_ls[2]]]  = ncvar_get(file_q2m, "longitude") 
lat_in[[var_ls[3]]]  = ncvar_get(file_prs, "latitude"); lon_in[[var_ls[3]]]  = ncvar_get(file_prs, "longitude") 
lat_in[[var_ls[4]]]  = ncvar_get(file_slp, "latitude"); lon_in[[var_ls[4]]]  = ncvar_get(file_slp, "longitude") 
lat_in[[var_ls[5]]]  = ncvar_get(file_pwt, "latitude"); lon_in[[var_ls[5]]]  = ncvar_get(file_pwt, "longitude") 
lat_in[[var_ls[6]]]  = ncvar_get(file_uwnd, "latitude"); lon_in[[var_ls[6]]]  = ncvar_get(file_uwnd, "longitude")
lat_in[[var_ls[7]]]  = ncvar_get(file_vwnd, "latitude"); lon_in[[var_ls[7]]]  = ncvar_get(file_vwnd, "longitude") 
lat_in[[var_ls[8]]]  = ncvar_get(file_sst, "latitude"); lon_in[[var_ls[8]]]  = ncvar_get(file_sst, "longitude")
lat_in[[var_ls[9]]]  = ncvar_get(file_sst, "latitude"); lon_in[[var_ls[9]]]  = ncvar_get(file_olr, "longitude")
lat_in[[var_ls[10]]]  = ncvar_get(file_tmp, "lat"); lon_in[[var_ls[10]]]  = ncvar_get(file_tmp, "lon")
lat_in[[var_ls[11]]] = ncvar_get(file_prt, "lat"); lon_in[[var_ls[11]]] = ncvar_get(file_prt, "lon")

## Extract lon/lat indices of interest
latInd_all  <- lonInd_all  <- Nx <- Ny <- NULL
lat_Nmax = length(unlist(lat_in[which.max(lapply(lat_in, nrow))]))
lon_Nmax = length(unlist(lon_in[which.max(lapply(lon_in, nrow))]))

for (var_i in var_ls) { 
  n_var = which(var_ls == var_i)
  
  if (cut_domain) {
    latInd = which(lat_in[[var_i]] >= lat_domain[n_var,1]  & lat_in[[var_i]] <= lat_domain[n_var,2])
    lonInd = which(lon_in[[var_i]] >= lon_domain[n_var,1]  & lon_in[[var_i]] <= lon_domain[n_var,2])

  } else {
    latInd = which(latitude>=-90 & latitude<=90)
    lonInd = which(longitude>=0 & longitude<=360)
  }
  # set to max length and bind df
  latInd = c(latInd, rep(NA, lat_Nmax - length(latInd)))
  lonInd = c(lonInd, rep(NA, lon_Nmax - length(lonInd)))
  latInd_all = cbind(latInd_all, latInd)
  lonInd_all = cbind(lonInd_all, lonInd)
  
  # get length of indicies
  Ny = c(Ny, length(na.omit(latInd)))
  Nx = c(Nx, length(na.omit(lonInd)))
}


##### ===== Get Time Index ===== #####

time = list()
time[[var_ls[1]]] = ncvar_get(file_hgt, "time") 
time[[var_ls[2]]] = ncvar_get(file_q2m, "time") 
time[[var_ls[3]]] = ncvar_get(file_prs, "time") 
time[[var_ls[4]]] = ncvar_get(file_slp, "time") 
time[[var_ls[5]]] = ncvar_get(file_pwt, "time") 
time[[var_ls[6]]] = ncvar_get(file_uwnd, "time") 
time[[var_ls[7]]] = ncvar_get(file_vwnd, "time") 
time[[var_ls[8]]] = ncvar_get(file_sst, "time") 
time[[var_ls[9]]] = ncvar_get(file_olr, "time")
time[[var_ls[10]]] = ncvar_get(file_tmp, "time")
time[[var_ls[11]]] = ncvar_get(file_prt, "time")

## Extract and process date
date = list()
for (i in 1:length(var_ls)) {
  if (var_ls[i] == 'nld.prate' | var_ls[i] == 'nld.tmp2m') {
    date[[var_ls[i]]] = as.POSIXct((time[[var_ls[i]]] - 7) * 86400, 
                                   origin = '1980-01-01', tz = 'utc') #timeseries of day 1 of obs
  } else {
    date[[var_ls[i]]] = as.POSIXct(as.Date(as.POSIXct(time[[var_ls[i]]] - 6*86400, 
                                                      origin = '1970-01-01', tz = 'utc'))) #timeseries of day 1 in fcst (fcsted date)
  }
}

## Define useful dimensions
#avg_fcst_day     = 7
# fcst_date = as.POSIXct(as.Date(as.POSIXct(time, origin = '1970-01-01', tz = 'utc') - (avg_fcst_day - 1) * 86400)) #timeseries of forecast date
# fcst_dt_nld = date_nld - avg_fcst_day*86400 #timeseries of forecast date

## Get index for overlapping times
for (i in 1:length(var_ls)) {
  d_i = data.frame(date = date[[var_ls[i]]], ind = 1:length(date[[var_ls[i]]]))
  if (i == 1) {
    ind_overlap = d_i
  } else {
    ind_overlap = inner_join(ind_overlap, d_i, by = 'date')
  }
}
colnames(ind_overlap) <- c('date', paste0('ind.', var_ls))


## Get month index in all time
indTime <- NULL
for (mon_i in 1:length(mon_in)) {
  ind_i = filter(ind_overlap, month(ind_overlap$date) == mon_in[mon_i])
  indTime = rbind(indTime, ind_i)
}
Ntime = nrow(indTime)
time = indTime$date


##### ===== Read variables based on lat, lon, time indicies ===== #####

## Read Variables
vars_df = list()
vars_df[['hgt']]  = ncvar_get(file_hgt, "HGT_500mb")[na.omit(lonInd_all[,1]),na.omit(latInd_all[,1]),indTime$ind.hgt] 
vars_df[['q2m']]  = ncvar_get(file_q2m, "SPFH_2maboveground")[na.omit(lonInd_all[,2]),na.omit(latInd_all[,2]),indTime$ind.q2m]
vars_df[['prs']]  = ncvar_get(file_prs, "PRES_surface")[na.omit(lonInd_all[,3]),na.omit(latInd_all[,3]),indTime$ind.prs]
vars_df[['slp']]  = ncvar_get(file_slp, "PRMSL_meansealevel")[na.omit(lonInd_all[,4]),na.omit(latInd_all[,4]),indTime$ind.slp]
vars_df[['pwt']]  = ncvar_get(file_pwt, "PWAT_entireatmosphere_consideredasasinglelayer_")[na.omit(lonInd_all[,5]),na.omit(latInd_all[,5]),indTime$ind.pwt]
vars_df[['uwnd']]  = ncvar_get(file_uwnd, "UGRD_850mb")[na.omit(lonInd_all[,6]),na.omit(latInd_all[,6]),indTime$ind.uwnd]
vars_df[['vwnd']]  = ncvar_get(file_vwnd, "VGRD_850mb")[na.omit(lonInd_all[,7]),na.omit(latInd_all[,7]),indTime$ind.vwnd]
vars_df[['sst']]  = ncvar_get(file_sst, "POT_5mbelowsealevel")[na.omit(lonInd_all[,8]),na.omit(latInd_all[,8]),indTime$ind.sst]
vars_df[['olr']]  = ncvar_get(file_olr, "ULWRF_topofatmosphere")[na.omit(lonInd_all[,9]),na.omit(latInd_all[,9]),indTime$ind.olr]
vars_df[['tmp']]  = ncvar_get(file_tmp, "TMP_2maboveground")[na.omit(lonInd_all[,10]),na.omit(latInd_all[,10]),indTime$ind.tmp]
vars_df[['prt']]  = ncvar_get(file_prt, "PRATE_surface")[na.omit(lonInd_all[,11]),na.omit(latInd_all[,11]),indTime$ind.prt]

## close netCDFs
nc_close(file_hgt); nc_close(file_q2m); nc_close(file_prs); nc_close(file_slp); nc_close(file_olr)
nc_close(file_pwt); nc_close(file_uwnd); nc_close(file_vwnd); nc_close(file_sst)


##### ===== Remap to match other coordinates ===== #####
# coarsen grid and cut domain
lat_ls <- lon_ls <- vars_new <- list()
for (var_i in var_ls) { 
  n_var = which(var_ls == var_i)
  lat_in_n = (lat_in[[var_i]])[na.omit(latInd_all[,n_var])]
  lon_in_n = (lon_in[[var_i]])[na.omit(lonInd_all[,n_var])]
  
  lat_out = seq(ceiling(min(lat_in_n)), floor(max(lat_in_n)), by = res_domain)
  lon_out = seq(ceiling(min(lon_in_n)), floor(max(lon_in_n)), by = res_domain)
  Ny[n_var] = length(na.omit(lat_out)); Nx[n_var] = length(na.omit(lon_out))
  
  lat_ls[[var_i]] = lat_out
  lon_ls[[var_i]] = lon_out
  vars_new[[var_i]] = remap(vars_df[[var_i]], x=lon_in_n, y=lat_in_n, z=1:Ntime,
              xto=lon_out, yto=lat_out, zto=1:Ntime)$var
}


##### ===== Reorganize predictor data as matrix with Ntime ===== #####
pred_cfs = array(NA,dim=c(length(var_ls),max(Nx),max(Ny),Ntime))

## creat array of predictors
for (var_i in var_ls) {
  n_var = which(var_ls == var_i)
  
  data = vars_new[[var_i]]
  # loop through lat & lon
  for(ilat in 1:Ny[n_var]){
    for(ilon in 1:Nx[n_var]){
      pred_cfs[n_var,ilon,ilat,]  = data[ilon,ilat,]
    }
  }
}


##### ===== Store grid for each variable ===== #####
grid <- array(NA, dim = c(length(var_ls), 2, max(Nx)*max(Ny)))

## create grid for each variable
for (var_i in var_ls) {
  n_var = which(var_ls == var_i)
  latitude = lat_ls[[var_i]]
  longitude = lon_ls[[var_i]]
  
  for(ilon in 1:Nx[n_var]){
    for(ilat in 1:Ny[n_var]){
      
      grid[n_var,1,(Ny[n_var] * (ilon - 1) + ilat)]  = latitude[ilat]
      grid[n_var,2,(Ny[n_var] * (ilon - 1) + ilat)] = longitude[ilon]
    }
  }
}


## save data for use in PLSR
data_out = list(pred_cfs = pred_cfs, grid = grid, vars = var_ls, dates = time)
setwd(out_model_dir)
rm(list=setdiff(ls(), c('data_out', 'file_data', 'l', 'file_lag_in')))
save(data_out, file = file_data)


###### ====== Average lagged variables loop through time ====== ######
# takes a long time (~20 mins)

# ##### ===== Load CFSv2 Data 
# setwd(out_model_dir)
# load(file = file_data)
nms = names(data_out)

## extract objects
pred_cfs = data_out[[nms[1]]]
grid = data_out[[nms[2]]]
var_nms = c(data_out[[nms[3]]], 'hru_cfs')
date_cfs = data_out[[nms[4]]]
rm(data_out)

## get lengths of vars 
Ntime           = length(pred_cfs[1,1,1,])
Nlat            = length(pred_cfs[1,1,,1])
Nlon            = length(pred_cfs[1,,1,1])
Nvar            = length(pred_cfs[,1,1,1])

## get overlapping indexes for lagged analysis
t_ind <- NULL
for (itime in 1:Ntime) {
  t_i = data.frame(matrix(which(as.Date(date_cfs) >= (as.Date(date_cfs[itime]) - l) & as.Date(date_cfs) <= (as.Date(date_cfs[itime]))),
                          nrow = 1))
  colnames(t_i) <- paste0('X', seq(1:length(t_i)))
  t_ind = rbind.fill(t_ind, t_i)
}


### === average lagged variables loop through time
# takes a long time (~20 mins)

## check if analysis is lagged and lagged data is saved yet
if (l > 0) {
  if (!file.exists(file_lag_in)){
    pred_avg = array(NA,dim=c(Nvar,Nlon,Nlat,Ntime))
    
    ## loop through time
    for (itime in 1:Ntime) {
      pred_avg[,,,itime]  = apply(pred_cfs[,,, na.omit(as.numeric(t_ind[itime,])) ], c(1,2,3),mean)
    }
    saveRDS(pred_avg, file = file_lag_in )
    
    pred_cfs = pred_avg
    rm(pred_avg)
  } else {
    pred_cfs = readRDS(file_lag_in)
  }
}