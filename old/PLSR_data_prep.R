# ================================================================================ 
# Step 1: process CRSv2 variables
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
# (9)  Surface Temperature (tmp2m_f)
# (10) Surface Prate (prate_f)
#
#
#   Adapted from P. Mendoza script (May 2016) and Andy Wood, NCAR
#   from PLSR R code (A. Wood) for cross-validated application of PCR/PLSR
# ================================================================================ 
rm(list=ls())

## --- Directories
cfs_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/cfs_annual_new/'
#TP_data = '/home/sabaker/s2s/analysis/files/cfsv2_files/2wk_avg/'
out_model_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/plsr_input/'

## --- Input Data
wk               = '2_3'  #'2_3' '3_4'
mon_in           = c(12,1,2) # month to analyze
mon_nm           = 'DFJ'
cut_domain       = T
cut_num          = 3
var_ls = c ('hgt', 'q2m', 'prs', 'slp', 'pwt', 'uwnd', 'vwnd', 'sst', 'tmp', 'prt')

## Set domain - will trim domain if not NA
lat_domain = matrix(c(25,80, -20,70, -20,30, -20,30, -20,70, #hgt, q2m, prs, slp, pwat
                      0,80, 0,80, -20,55, NA,NA, NA,NA), #uwnd, vwnd, sst, tmp, prt
                    ncol = 2, byrow = T)
lon_domain = matrix(c(100,320, NA,NA, NA,NA, NA,NA, 100,320, #hgt, q2m, prs, slp, pwat
                      NA,NA, NA,NA, NA,NA, NA,NA, NA,NA), #uwnd, vwnd, sst, tmp, prt
                    ncol = 2, byrow = T)
# replace NAs with full domain
for (i in 1:10) {
  if (is.na(lat_domain[i,1])) {lat_domain[i,1] = -90}
  if (is.na(lat_domain[i,2])) {lat_domain[i,2] = 90}
  if (is.na(lon_domain[i,1])) {lon_domain[i,1] = 0}
  if (is.na(lon_domain[i,2])) {lon_domain[i,2] = 360}
}


## --- libraries
library(lubridate)
library(ncdf4)
library(abind)
library(OceanView)
library(pls)
library(dplyr)


##### ===== Load Predictors & Predictanad ===== #####

## PredictorsOpen netCDF files
setwd(cfs_dir)
file_hgt   = nc_open(paste0('cfsv2.', wk, '.z500_f.nc')) 
file_q2m   = nc_open(paste0('cfsv2.', wk, '.q2m_f.nc'))
file_prs   = nc_open(paste0('cfsv2.', wk, '.pressfc_f.nc'))
file_slp   = nc_open(paste0('cfsv2.', wk, '.prmsl_f.nc'))
file_pwt   = nc_open(paste0('cfsv2.', wk, '.pwat_f.nc'))
file_uwnd  = nc_open(paste0('cfsv2.', wk, '.U_wnd850_f.nc'))
file_vwnd  = nc_open(paste0('cfsv2.', wk, '.V_wnd850_f.nc'))
file_sst   = nc_open(paste0('cfsv2.', wk, '.ocnsst_h.nc'))
file_tmp   = nc_open(paste0('cfsv2.', wk, '.tmp2m_f.nc'))
file_prt   = nc_open(paste0('cfsv2.', wk, '.prate_f.nc'))


##### ===== Get Lat & Lon Index ===== #####

# Extract lat/lon coordinates - different for variables
latitude_in       = ncvar_get(file_hgt, "latitude" ) 
longitude_in      = ncvar_get(file_hgt, "longitude" )
latitude_odd_in   = ncvar_get(file_q2m, "latitude" )  #q2m, prs, pwat
longitude_odd_in  = ncvar_get(file_q2m, "longitude" ) 
latitude_sst_in   = ncvar_get(file_sst, "latitude" ) 
longitude_sst_in  = ncvar_get(file_sst, "longitude" ) 
latitude_PT_in   = ncvar_get(file_prt, "lat" )
longitude_PT_in  = ncvar_get(file_prt, "lon" )
# latitude_odd2_in   = ncvar_get(file_pwt, "latitude" ) 
# longitude_odd2_in  = ncvar_get(file_pwt, "longitude" ) 

## Extract lon/lat indices of interest
latInd_all  <- lonInd_all  <- Nx <- Ny<- NULL
lat_max = max(length(latitude_in), length(latitude_odd_in), length(latitude_PT_in), length(latitude_sst_in))
lon_max = max(length(longitude_in), length(longitude_odd_in), length(longitude_PT_in), length(longitude_sst_in))
for (var_i in 1:length(var_ls)) {
  
  if (cut_domain) {
    ## hgt, slp, uwnd, vwnd
    if (var_i == 1 | var_i == 4 | var_i == 6 | var_i == 7) {
      latInd = which(latitude_in >= lat_domain[var_i,1]  & latitude_in <= lat_domain[var_i,2])
      lonInd = which(longitude_in >= lon_domain[var_i,1] & longitude_in <= lon_domain[var_i,2])
    }
    ## q2m, prs, pwat
    if (var_i == 2 | var_i == 3 | var_i == 5) {
      latInd = which(latitude_odd_in >= lat_domain[var_i,1]  & latitude_odd_in  <= lat_domain[var_i,2])
      lonInd = which(longitude_odd_in >= lon_domain[var_i,1] & longitude_odd_in <= lon_domain[var_i,2])
    }
    ## sst
    if (var_i == 8) {
      latInd = which(latitude_sst_in >= lat_domain[var_i,1]  & latitude_sst_in  <= lat_domain[var_i,2])
      lonInd = which(longitude_sst_in >= lon_domain[var_i,1] & longitude_sst_in <= lon_domain[var_i,2])
    }
    if (var_i == 9 | var_i == 10) {
      latInd = which(latitude_PT_in >= lat_domain[var_i,1]  & latitude_PT_in  <= lat_domain[var_i,2])
      lonInd = which(longitude_PT_in >= lon_domain[var_i,1] & longitude_PT_in <= lon_domain[var_i,2])
    }
  } else {
    latInd = which(latitude>=-90 & latitude<=90)
    lonInd = which(longitude>=0 & longitude<=360)
  }
  # set to max length and bind df
  latInd = c(latInd, rep(NA, lat_max - length(latInd)))
  lonInd = c(lonInd, rep(NA, lon_max - length(lonInd)))
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
time[[var_ls[9]]] = ncvar_get(file_tmp, "time")
time[[var_ls[10]]] = ncvar_get(file_prt, "time")

# Extract and process date
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
hgt_df   = ncvar_get(file_hgt, "HGT_500mb")[na.omit(lonInd_all[,1]),na.omit(latInd_all[,1]),indTime$ind.hgt] 
q2m_df   = ncvar_get(file_q2m, "SPFH_2maboveground")[na.omit(lonInd_all[,2]),na.omit(latInd_all[,2]),indTime$ind.q2m]
prs_df   = ncvar_get(file_prs, "PRES_surface")[na.omit(lonInd_all[,3]),na.omit(latInd_all[,3]),indTime$ind.prs]
slp_df   = ncvar_get(file_slp, "PRMSL_meansealevel")[na.omit(lonInd_all[,4]),na.omit(latInd_all[,4]),indTime$ind.slp]
pwt_df   = ncvar_get(file_pwt, "PWAT_entireatmosphere_consideredasasinglelayer_")[na.omit(lonInd_all[,5]),na.omit(latInd_all[,5]),indTime$ind.pwt]
uwnd_df   = ncvar_get(file_uwnd, "UGRD_850mb")[na.omit(lonInd_all[,6]),na.omit(latInd_all[,6]),indTime$ind.uwnd]
vwnd_df   = ncvar_get(file_vwnd, "VGRD_850mb")[na.omit(lonInd_all[,7]),na.omit(latInd_all[,7]),indTime$ind.vwnd]
sst_df   = ncvar_get(file_sst, "POT_5mbelowsealevel")[na.omit(lonInd_all[,8]),na.omit(latInd_all[,8]),indTime$ind.sst]
tmp_df   = ncvar_get(file_tmp, "TMP_2maboveground")[na.omit(lonInd_all[,9]),na.omit(latInd_all[,9]),indTime$ind.tmp]
prt_df   = ncvar_get(file_prt, "PRATE_surface")[na.omit(lonInd_all[,10]),na.omit(latInd_all[,10]),indTime$ind.prt]

## close netCDFs
nc_close(file_hgt); nc_close(file_q2m); nc_close(file_prs); nc_close(file_slp)
nc_close(file_pwt); nc_close(file_uwnd); nc_close(file_vwnd); nc_close(file_sst)


##### ===== Remap to match other coordinates ===== #####

## SST - coarsen grid
longitude_sst_in = longitude_sst_in[na.omit(lonInd_all[,8])]; latitude_sst_in = latitude_sst_in[na.omit(latInd_all[,8])]
longitude_sst = unique(round(longitude_sst_in)[1:(length(longitude_sst_in)-1)])
latitude_sst = unique(round(latitude_sst_in)[2:(length(latitude_sst_in)-1)])
Ny[8] = length(na.omit(latitude_sst)); Nx[8] = length(na.omit(longitude_sst))
sst = remap(sst_df, x=longitude_sst_in, y=latitude_sst_in, z=1:Ntime,
            xto=longitude_sst, yto=latitude_sst, zto=1:Ntime)$var

## Precip and Temp - coarsen grid & smaller domain
longitude_PT = unique(round(longitude_PT_in)[1:(length(longitude_PT_in))])
latitude_PT = unique(round(latitude_PT_in)[2:(length(latitude_PT_in)-1)])
Ny[9] = Ny[10] = length(na.omit(latitude_PT)); Nx[9] = Nx[10] = length(na.omit(longitude_PT))
tmp = remap(tmp_df, x=longitude_PT_in, y=latitude_PT_in, z=1:Ntime,
            xto=longitude_PT, yto=latitude_PT, zto=1:Ntime)$var
prt = remap(prt_df, x=longitude_PT_in, y=latitude_PT_in, z=1:Ntime,
            xto=longitude_PT, yto=latitude_PT, zto=1:Ntime)$var

## q2m, pwat, pressfc - change grid resolution
# longitude_odd = longitude[117:236]; latitude_odd = latitude[68:136]
longitude_odd_in_q2m = longitude_odd_in[na.omit(lonInd_all[,2])]; latitude_odd_in_q2m = latitude_odd_in[na.omit(latInd_all[,2])]
longitude_odd_q2m = unique(round(longitude_odd_in_q2m)[2:(length(longitude_odd_in_q2m))])
latitude_odd_q2m = unique(round(latitude_odd_in_q2m)[2:(length(latitude_odd_in_q2m)-1)])
Ny[2] = length(na.omit(latitude_odd_q2m)); Nx[2]= length(na.omit(longitude_odd_q2m))
q2m = remap(q2m_df, x=longitude_odd_in_q2m, y=latitude_odd_in_q2m, z=1:Ntime,
            xto=longitude_odd_q2m, yto=latitude_odd_q2m, zto=1:Ntime)$var

longitude_odd_in_pwt = longitude_odd_in[na.omit(lonInd_all[,5])]; latitude_odd_in_pwt = latitude_odd_in[na.omit(latInd_all[,5])]
longitude_odd_pwt = unique(round(longitude_odd_in_pwt)[2:(length(longitude_odd_in_pwt)-1)])
latitude_odd_pwt = unique(round(latitude_odd_in_pwt)[2:(length(latitude_odd_in_pwt)-1)])
Ny[5] = length(na.omit(latitude_odd_pwt)); Nx[5]= length(na.omit(longitude_odd_pwt))
pwt = remap(pwt_df, x=longitude_odd_in_pwt, y=latitude_odd_in_pwt, z=1:Ntime,
            xto=longitude_odd_pwt, yto=latitude_odd_pwt, zto=1:Ntime)$var

longitude_odd_in_prs = longitude_odd_in[na.omit(lonInd_all[,3])]; latitude_odd_in_prs = latitude_odd_in[na.omit(latInd_all[,3])]
longitude_odd_prs = unique(round(longitude_odd_in_prs)[2:(length(longitude_odd_in_prs))])
latitude_odd_prs = unique(round(latitude_odd_in_prs)[2:(length(latitude_odd_in_prs)-1)])
Ny[3] = length(na.omit(latitude_odd_prs)); Nx[3]= length(na.omit(longitude_odd_prs))
prs = remap(prs_df, x=longitude_odd_in_prs, y=latitude_odd_in_prs, z=1:Ntime,
            xto=longitude_odd_prs, yto=latitude_odd_prs, zto=1:Ntime)$var


# pwat = remap(pwt_df, x=longitude_odd_in, y=latitude_odd_in, z=1:Ntime,
#             xto=longitude_odd, yto=latitude_odd, zto=1:Ntime)$var
# prs = remap(prs_df, x=longitude_odd_in, y=latitude_odd_in, z=1:Ntime,
#             xto=longitude_odd, yto=latitude_odd, zto=1:Ntime)$var


##### ===== Reorganize predictor data as matrix with Ntime ===== #####
data_ls = list(hgt_df, q2m, prs, slp_df, pwt, uwnd_df, vwnd_df, sst, tmp, prt)
pred_cfs = array(NA,dim=c(length(var_ls),max(Nx),max(Ny),Ntime))
## creat array of predictors
for (var_i in 1:length(var_ls)) {
  data = data_ls[[var_i]]
  # loop through lat & lon
  for(ilat in 1:Ny[var_i]){
    for(ilon in 1:Nx[var_i]){
      pred_cfs[var_i,ilon,ilat,]  = data[ilon,ilat,]
    }
  }
}


##### ===== Store grid for each variable ===== #####
grid <- array(NA, dim = c(length(var_ls), 2, max(Nx)*max(Ny)))
lat_ls = list(latitude_in[na.omit(latInd_all[,1])], latitude_odd_q2m, latitude_odd_prs, latitude_in[na.omit(latInd_all[,4])], 
              latitude_odd_pwt, latitude_in[na.omit(latInd_all[,6])], latitude_in[na.omit(latInd_all[,7])], 
              latitude_sst, latitude_PT, latitude_PT)
lon_ls = list(longitude_in[na.omit(lonInd_all[,1])], longitude_odd_q2m, longitude_odd_prs, longitude_in[na.omit(lonInd_all[,4])], 
              longitude_odd_pwt, longitude_in[na.omit(lonInd_all[,6])], longitude_in[na.omit(lonInd_all[,7])], 
              longitude_sst, longitude_PT, longitude_PT)

## create grid for each variable
for (var_i in 1:length(var_ls)) {
  latitude = lat_ls[[var_i]]
  longitude = lon_ls[[var_i]]
  
  for(ilon in 1:Nx[var_i]){
    for(ilat in 1:Ny[var_i]){
      
      grid[var_i,1,(Ny[var_i] * (ilon - 1) + ilat)]  = latitude[ilat]
      grid[var_i,2,(Ny[var_i] * (ilon - 1) + ilat)] = longitude[ilon]
    }
  }
}


## save data for use in PLSR
data_out = list(pred_cfs = pred_cfs, grid = grid, vars = var_ls, dates = time)
setwd(out_model_dir)
if (cut_domain) {
  file = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.cutdomain', cut_num,'.rds')
} else {
  file = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.rds')
}
rm(list=setdiff(ls(), c('data_out', 'file')))
save(data_out, file = file)


# ## Remove raw data matrices
# rm(hgt_df, q2m_df, prs_df, slp_df, pwt_df, uwnd_df, vwnd_df, sst_df, sst) 
# 
# for(ilat in 1:Ny){
#   for(ilon in 1:Nx){
#     grid[(Nx * (ilat - 1) + ilon),1] = latitude[ilat]
#     grid[(Nx * (ilat - 1) + ilon),2] = longitude[ilon]
#   }
# }
# for(ilat in 1:Ny_sst){
#   for(ilon in 1:Nx_sst){
#     grid[(Nx * (ilat - 1) + ilon),3] = latitude_sst[ilat]
#     grid[(Nx * (ilat - 1) + ilon),4] = longitude_sst[ilon]
#   }
# }
# for(ilat in 1:Ny_PT){
#   for(ilon in 1:Nx_PT){
#     grid[(Nx * (ilat - 1) + ilon),5] = latitude_PT[ilat]
#     grid[(Nx * (ilat - 1) + ilon),6] = longitude_PT[ilon]
#     #to undo: pred_i[,(Nx * (ilat - 1) + ilon)]  = pred_cfs[ipred,ilon,ilat,]
#   }
# }