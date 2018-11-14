# ================================================================================ 
# PLSR code for processing CFSv2 variables for predicting 
# NLDAS temperature and precipitation for 202 HUCs in CONUS domain
#   Created by S. Baker 
#
# Number of CFSv2 predictor variables
# (1) Geopotential height (z500_f)
# (2) Specific Humid 2m (q2m_f)
# (3) Surface Pressure (pressfc_f)
# (4) Sea Level Pressure (meanSLP)
# (5) precipitable water (pwat_f)
# (6) Zonal winds (U_wnd850_f)
# (7) Meridional winds (V_wnd850_f)
# (8) Sea Surface Temperature (ocnsst_h)
#
#
#   Adapted from P. Mendoza script (May 2016) and Andy Wood, NCAR
#   from PLSR R code (A. Wood) for cross-validated application of PCR/PLSR
# ================================================================================ 
rm(list=ls())

## --- Directories
cfs_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/cfs_annual_new/'
TP_data = '/home/sabaker/s2s/analysis/files/cfsv2_files/2wk_avg/'

## --- Input Data
wk = '2_3'
avg_fcst_day = 7
hru = '1304'

## --- old inputs
#flagflex         = 0       # Flag to define if predictand is fixed to a seasonal value (0) or not (1)
MasterDir        = '/hydro_d2/pmendoza/overtheloop' # Define master directory (NCAR)
mon_in           = 1 # month to analyze
Ntime_sub        = 14       # Define number of timesteps that each sub-period (i.e. validation) will have


## --- libraries
library(lubridate)
library(ncdf4)
library(abind)
library(OceanView)
library(pls)
library(dplyr)

## --- custom R functions
# FcnDir = paste(MasterDir,'/Rscripts/R_functions/',sep='') # Location of function code
# source(paste(FcnDir, 'PlotMapRegGrid.R',sep=''))           # Plot correlation maps
# source(paste(FcnDir, 'PlotTwoClustersMap.R',sep=''))       # Plot two clusters
# source(paste(FcnDir, 'CreateLandMask.R',sep=''))           # Returns matrix with land mask
# source(paste(FcnDir, 'PredictPLSRCrossVal.R',sep=''))      # Returns matrix with land mask
# source(paste(FcnDir, 'BoxplotEnsForecastsv2.R', sep = '')) # Boxplots w/ ensembles
# source(paste(FcnDir, 'ScatterPlotsCorBias.R', sep = ''))   # Boxplots with ensembles
# Define long variable names for title of plots
# ReanalNames   = c('Geopotential Height 500mb')
# ReanalShort   = c('z500')


##### ===== Part 1: Load Predictors & Predictanad ===== #####

## Input number of predictors
Npred = 2 

## Predictand - Open netCDF files
setwd(TP_data)
file_nldprate   = nc_open(paste0('nldas.1_2.prate.nc')) 

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

################### add cfsv2 T and P

# Extract lat/lon coordinates and time! differen!!!!
latitude       = ncvar_get(file_hgt, "latitude" )  # Latitude
longitude      = ncvar_get(file_hgt, "longitude" )  # Longitude
latitude_odd   = ncvar_get(file_q2m, "latitude" )  # Latitude
longitude_odd  = ncvar_get(file_q2m, "longitude" )  # Longitude
latitude_sst   = ncvar_get(file_sst, "latitude" )  # Latitude
longitude_sst  = ncvar_get(file_sst, "longitude" )  # Longitude

# Extract lon/lat indices (does nothing now, could be useful later)
minlat     = min(latitude); maxlat    = max(latitude)
minlon     = min(longitude); maxlon   = max(longitude)
latindex   = which(latitude>=minlat & latitude<=maxlat)
lonindex   = which(longitude>=minlon & longitude<=maxlon)
Ny         = length(latindex)
Nx         = length(lonindex)

# Extract hru coordinates for NLDAS
hru_v   = ncvar_get(file_nldprate, "hru")  # hru
hruInd = which(hru_v == hru)

# times - data are missing
# time_hgt = ncvar_get(file_hgt, "time") 
# time_q2m = ncvar_get(file_q2m, "time") 
# time_prs = ncvar_get(file_prs, "time") 
# time_slp = ncvar_get(file_slp, "time") 
# time_pwt = ncvar_get(file_pwt, "time") 
# time_uwnd = ncvar_get(file_uwnd, "time") 
# time_vwnd = ncvar_get(file_vwnd, "time") 
# time_sst = ncvar_get(file_sst, "time") 
# time_nld = ncvar_get(file_nldprate, "time")
var_ls = c ('nld', 'hgt', 'q2m', 'prs', 'slp', 'pwt', 'uwnd', 'vwnd', 'sst')
time = list()
time[[var_ls[1]]] = ncvar_get(file_nldprate, "time")
time[[var_ls[2]]] = ncvar_get(file_hgt, "time") 
time[[var_ls[3]]] = ncvar_get(file_q2m, "time") 
time[[var_ls[4]]] = ncvar_get(file_prs, "time") 
time[[var_ls[5]]] = ncvar_get(file_slp, "time") 
time[[var_ls[6]]] = ncvar_get(file_pwt, "time") 
time[[var_ls[7]]] = ncvar_get(file_uwnd, "time") 
time[[var_ls[8]]] = ncvar_get(file_vwnd, "time") 
time[[var_ls[9]]] = ncvar_get(file_sst, "time") 

# Extract and process time
date = list()
for (i in 1:length(var_ls)) {
  if (var_ls[i] == 'nld') {
    date[[var_ls[i]]] = as.POSIXct((time[[var_ls[i]]] - 7) * 86400, 
                                   origin = '1980-01-01', tz = 'utc') #timeseries of day 1 of obs
  } else {
    date[[var_ls[i]]] = as.POSIXct(as.Date(as.POSIXct(time[[var_ls[i]]] - 6*86400, 
                                                      origin = '1970-01-01', tz = 'utc'))) #timeseries of day 1 in fcst (fcsted date)
  }
}
# date_nld = as.POSIXct((time_nld - 7) * 86400, origin = '1980-01-01', tz = 'utc') #timeseries of day 1 of obs
# date_hgt = as.POSIXct(as.Date(as.POSIXct(time_hgt - 6*86400, origin = '1970-01-01', tz = 'utc'))) #timeseries of day 1 in fcst (fcsted date)
# date_q2m = as.POSIXct(as.Date(as.POSIXct(time_q2m - 6*86400, origin = '1970-01-01', tz = 'utc'))) #timeseries of day 1 in fcst (fcsted date)
# date_prs = as.POSIXct(as.Date(as.POSIXct(time_prs - 6*86400, origin = '1970-01-01', tz = 'utc')))
# date_slp = as.POSIXct(as.Date(as.POSIXct(time_slp - 6*86400, origin = '1970-01-01', tz = 'utc')))
# date_pwt = as.POSIXct(as.Date(as.POSIXct(time_pwt - 6*86400, origin = '1970-01-01', tz = 'utc')))
# date_uwnd = as.POSIXct(as.Date(as.POSIXct(time_uwnd - 6*86400, origin = '1970-01-01', tz = 'utc')))
# date_vwnd = as.POSIXct(as.Date(as.POSIXct(time_vwnd - 6*86400, origin = '1970-01-01', tz = 'utc')))
# date_sst = as.POSIXct(as.Date(as.POSIXct(time_sst - 6*86400, origin = '1970-01-01', tz = 'utc')))

## Define useful dimensions
# fcst_date = as.POSIXct(as.Date(as.POSIXct(time, origin = '1970-01-01', tz = 'utc') - (avg_fcst_day - 1) * 86400)) #timeseries of forecast date
# fcst_dt_nld = date_nld - avg_fcst_day*86400 #timeseries of forecast date

## Get index for each variable for month 
# (might want to do this in loop to keep more times if fewer variables processed)
for (i in 1:length(var_ls)) {
  d_i = data.frame(date = date[[var_ls[i]]], ind = 1:length(date[[var_ls[i]]]))
  if (i == 1) {
    ind_overlap = d_i
  } else {
    ind_overlap = inner_join(ind_overlap, d_i, by = 'date')
  }
}
colnames(ind_overlap) <- c('date', paste0('ind.', var_ls))
# d_nld = data.frame(date = date_nld, ind.nld = 1:length(date_nld))
# d_hgt = data.frame(date = date_hgt, ind.hgt = 1:length(date_hgt))
# d_q2m = data.frame(date = date_q2m, ind.q2m = 1:length(date_q2m))
# d_prs = data.frame(date = date_prs, ind.prs = 1:length(date_prs))
# d_slp = data.frame(date = date_slp, ind.slp = 1:length(date_slp))
# d_pwt = data.frame(date = date_pwt, ind.pwt = 1:length(date_pwt))
# d_uwnd = data.frame(date = date_uwnd, ind.uwnd = 1:length(date_uwnd))
# d_vwnd = data.frame(date = date_vwnd, ind.vwnd = 1:length(date_vwnd))
# d_sst = data.frame(date = date_sst, ind.sst = 1:length(date_sst))

## join to get overlapping indicies
# ind_overlap = inner_join(d_hgt, d_q2m, by = 'date')
# ind_overlap = inner_join(ind_overlap, d_nld, by = 'date')
# rm(d_hgt, d_q2m, d_nld)

## get month index in all time
# indTime = which(month(ind_overlap$date) == mon_in)
indTime = filter(ind_overlap, month(ind_overlap$date) == mon_in)
# indTime_nld = which(month(date_nld[ind_overlap$ind.nld]) == mon_in)
# indTime_hgt = which(month(date_hgt[ind_overlap$ind.hgt]) == mon_in)
# indTime_q2m = which(month(date_q2m[ind_overlap$ind.q2m]) == mon_in)
ntime = nrow(indTime)

## Read Variables
nld_df   = ncvar_get(file_nldprate, "prate")[hruInd, indTime$ind.nld] * 24 ## convert mm/hr to mm/d
hgt_df   = ncvar_get(file_hgt, "HGT_500mb")[lonindex,latindex,indTime$ind.hgt] 
q2m_df   = ncvar_get(file_q2m, "SPFH_2maboveground")[lonindex,latindex,indTime$ind.q2m]
prs_df   = ncvar_get(file_prs, "PRES_surface")[lonindex,latindex,indTime$ind.prs]
slp_df   = ncvar_get(file_slp, "PRMSL_meansealevel")[lonindex,latindex,indTime$ind.slp]
pwt_df   = ncvar_get(file_pwt, "PWAT_entireatmosphere_consideredasasinglelayer_")[lonindex,latindex,indTime$ind.pwt]
uwnd_df   = ncvar_get(file_uwnd, "UGRD_850mb")[lonindex,latindex,indTime$ind.uwnd]
vwnd_df   = ncvar_get(file_vwnd, "VGRD_850mb")[lonindex,latindex,indTime$ind.vwnd]
sst_df   = ncvar_get(file_sst, "POT_5mbelowsealevel")[lonindex,latindex,indTime$ind.sst]

# close netCDFs
nc_close(file_nldprate); nc_close(file_hgt); nc_close(file_q2m); nc_close(file_prs); nc_close(file_slp)
nc_close(file_pwt); nc_close(file_uwnd); nc_close(file_vwnd); nc_close(file_sst)



## --- Reorganize predictor data as matrix with ntime
pred_cfs = array(0,dim=c(Npred,Nx,Ny,ntime))

for(ilat in 1:Ny){
  for(ilon in 1:Nx){
    pred_cfs[1,ilon,ilat,]  = hgt_df[ilon,ilat,]
    pred_cfs[2,ilon,ilat,]  = q2m_df[ilon,ilat,]
    # pred_cfs[3,ilon,ilat,]  = hgt_df[ilon,ilat,]
    # pred_cfs[4,ilon,ilat,]  = hgt_df[ilon,ilat,]
    # pred_cfs[5,ilon,ilat,]  = hgt_df[ilon,ilat,]
    # pred_cfs[6,ilon,ilat,]  = hgt_df[ilon,ilat,]
    # pred_cfs[7,ilon,ilat,]  = hgt_df[ilon,ilat,]
    # pred_cfs[8,ilon,ilat,]  = hgt_df[ilon,ilat,]
  }
}

## Remove raw data matrices
rm(hgt_df, q2m_df) #, prs_df, slp_df, pwat_df, uwnd_df, vwnd_df, sst_df) 



##### ===== Part 3: pre-process ===== #####
Ncomp            = 4  # number of components

## store grid
grid <- matrix(NA, nrow = Nx*Ny, ncol = 4)
for(ilat in 1:Ny){
  for(ilon in 1:Nx){
    grid[(Nx * (ilat - 1) + ilon),1] = latitude[ilat]
    grid[(Nx * (ilat - 1) + ilon),2] = longitude[ilon]
    grid[(Nx * (ilat - 1) + ilon),3] = latitude_odd[ilat]
    grid[(Nx * (ilat - 1) + ilon),4] = longitude_odd[ilon]
    #to undo: pred_i[,(Nx * (ilat - 1) + ilon)]  = pred_cfs[ipred,ilon,ilat,]
  }
}
grid = data.frame(lat = grid[,1], lon = grid[,2],
                  lat_odd = grid[,3], lon_odd = grid[,4])

#### ? get groups for validation?
Ntime           = length(indTime_nld)
Ngroups         = ceiling(Ntime/Ntime_sub)
cutgroup        = sapply(1:Ngroups, function(j) c(rep(j,Ntime_sub)))[1:Ntime] #not random
xvalindex       = split(1:Ntime, cutgroup) # validate on these indexes
Nsegments       = length(xvalindex)-1

# Scale NLDAS variable
nldasVar      = scale(nld_df)
Ymean         = mean(nld_df)     # Mean nldasVar during training period
Ysd           = sd(nld_df)       # Standard deviation of nldasVar during training period


##### ===== Part 4: Do PLSR over the entire field ===== #####

selpred  = nrow(pred_cfs) #c(1) 
# Scores   = array(NA,dim=c(Ncomp,Nlead,Ngroups,Ntime))

## loop - arrange predictors
Alldata <- NULL
for(ipred in 1:selpred){ # Start loop over selected predictors
  
  # get predictor in correct format
  pred_i = array(0,dim=c(Ntime,Nx*Ny))
  for(ilat in 1:Ny){
    for(ilon in 1:Nx){
      pred_i[,(Nx * (ilat - 1) + ilon)]  = pred_cfs[ipred,ilon,ilat,]
    }
  }
  Alldata   <- cbind(Alldata, pred_i) 
  
} # End loop over predictors

# scale - necessary?
ArrayCFSR        = scale(Alldata)
# ArrayCFSR        = Alldata


igroup =1
# for(igroup in 1:Ngroups){ # Start loop over groups

# Define predictand
nldasMod      = nldasVar[-xvalindex[[igroup]]]        # Observed nldasVar for the training period
Ypredictand   = nldasMod    # scaled above-Column vector containing the predictand to be used ([(Ngroups-1)*Ntime_sub] x 1)
#Ymean         = mean(nldasMod)     # Mean nldasVar during training period
#Ysd           = sd(nldasMod)       # Standard deviation of nldasVar during training period
Yfull         = nldasVar
Yverif        = nldasVar[xvalindex[[igroup]]]     

# Extract predictors (separate training from verification dataset)
Xtrain      = ArrayCFSR[-xvalindex[[igroup]],]
Xverif      = ArrayCFSR[ xvalindex[[igroup]],]
Xfull       = ArrayCFSR[ ,]

# Fit PLSR model and predict!
# from from 'PredictPLSRCrossVal.R' function
mydf        = data.frame(Ypredictand)
mydf$X      = Xtrain
plsrmod     = plsr(Ypredictand ~ X, ncomp = Ncomp, data = mydf, method='simpls', validation = 'CV')
Yhat_ver    = predict(plsrmod, ncomp = Ncomp, Xverif, se.fit=T)
Yhat        = predict(plsrmod, ncomp = Ncomp, Xfull, se.fit=T)

## example does cross validate on training data...(vinette)

res = Yhat[,1,1] - Yfull # residulas
mean(abs(res))
plot(Yfull, Yhat) # negative values?
abline(0,1)

## unscaled data
# Yhat_unscaled = Yhat * attr(nldasVar, 'scaled:scale') + attr(nldasVar, 'scaled:center')
Yhat_unscaled = Yhat * Ysd + Ymean
Yhat_unscaled = ifelse(Yhat_unscaled <0, 0, Yhat_unscaled) # remove negatives
res_unscaled = Yhat_unscaled - nld_df
plot(nld_df, Yhat_unscaled); abline(0,1)


## get gridded coefficients
comp1 = plsrmod$coefficients[,1,1]

P  = scores(plsrmod)             # CCA X component canonical variables
# ie PC time series
#  cols=cvs (time dim by min(cols(X,Y)
U  = loadings(plsrmod)	      # Eigen Vectors (ie EOFs)

# use Y variance explained instead of eigenvalues for weights
# W = plsrmod$Xvar   # eigenval. (X variance)


##
comp1_unscaled = data.frame(comp1 *  attr(ArrayCFSR, 'scaled:scale') + 
                         attr(ArrayCFSR, 'scaled:center'))
# comp1_z500 = data.frame(comp1 *  attr(ArrayCFSR, 'scaled:scale') + 
#                           attr(ArrayCFSR, 'scaled:center'))[1:(Nx*Ny),]
# comp1_x = data.frame(comp1 *  attr(ArrayCFSR, 'scaled:scale') + 
#                        attr(ArrayCFSR, 'scaled:center'))[(Nx*Ny+1):2*(Nx*Ny),]
# df_plot = cbind.data.frame(comp1_z500, comp1_x, grid)
df_plot = cbind.data.frame(z500 = comp1_unscaled[1:(Nx*Ny),],
                           q2m = comp1_unscaled[((Nx*Ny+1):(2*(Nx*Ny))),],
                           grid)



### === plotting things
library(ggplot2)
#library(maps)
#library(RColorBrewer)
library(data.table)
# Lon / Lat labels
ewbrks = seq(0, 360, 20)
nsbrks = seq(-90, 90, 15)
# ewbrks = seq(70, 355, 20)
# nsbrks = seq(-60, 90, 15)
ewlbls = unlist(lapply(ewbrks, function(x) ifelse(x > 0, paste0(x, 'ºE'), 	ifelse(x < 0, paste0(abs(x), 'ºW'), x))))
nslbls = unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), 	'ºS'), ifelse(x > 0, paste0(x, 'ºN') ,x))))

## MAP - correct projection, contries
worldmap = map_data('world2')
setnames(worldmap, c('X','Y', 'PID', 'POS', 'region', 'subregoin'))

## bins - z500
dif = round((max(df_plot$z500) - min(df_plot$z500))/12, -1)
min_bin = floor(min(df_plot$z500)/10) *10
max_bin = ceiling(max(df_plot$z500)/10) *10
dif = round(ceiling((max_bin - min_bin)*1/11), -1)
bin_z500 = seq(min_bin, max_bin, dif)
if (max(bin_z500) < max_bin) {bin_z500 = c(bin_z500, max(bin_z500)+dif)}
cust_pal_z500 <- rev(brewer.pal(n=length(bin_z500), name="RdYlBu"))
## set up bins and color palette
if(max_bin > 1000) {dig= 4} else {dig=3}
df_plot$bins = cut(df_plot$z500, breaks = bin_z500, dig.lab = dig)
cust_pal = cust_pal_z500
# set lat and lon for figure
df_plot$lat_p = df_plot$lat
df_plot$lon_p = df_plot$lon

## bins - q2m
scale_i = 1/10000
dif = round((max(df_plot$q2m) - min(df_plot$q2m))/11, 3)
min_bin = floor(min(df_plot$q2m)/scale_i) *scale_i
max_bin = ceiling(max(df_plot$q2m)/scale_i) *scale_i
#dif = round(ceiling((max_bin - min_bin)*1/11), -1)
bin_q2m = seq(min_bin, max_bin, dif)
if (max(bin_q2m) < max_bin) {bin_q2m = c(bin_q2m, max(bin_q2m)+dif)}
cust_pal_q2m <- rev(brewer.pal(n=length(bin_q2m), name="RdYlBu"))
## set up bins and color palette
if(max_bin > 1000) {dig= 4} else {dig=3}
df_plot$bins = cut(df_plot$q2m, breaks = bin_q2m, dig.lab = dig)
cust_pal = cust_pal_q2m
# set lat and lon for figure
df_plot$lat_p = df_plot$lat_odd
df_plot$lon_p = df_plot$lon_odd

labs = levels(df_plot$bins)

## plot with world map and discrete colors
ggplot() + geom_raster(data = df_plot, aes(x = lon_p, y = lat_p, fill = bins)) +
  scale_fill_manual(values = cust_pal, name = ' ', drop = F, labels = labs,
                    guide = guide_legend(ncol = 1, label.position = 'right', label.hjust = 0,
                                         title.position = 'top', title.hjust = 0.5, reverse = T,
                                         keywidth = 1, keyheight = 1)) +
  theme_bw() +
  coord_equal() + # ratio of x and y...coord equal preserves the relative sizes of each
  geom_polygon(data = worldmap, aes(X,Y,group=PID), color = 'black', size = 0.4, fill = NA) +
  xlab('') +  ylab('') +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0), limits = c(min(df_plot$lon_p), max(df_plot$lon_p))) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0), limits = c(min(df_plot$lat_p), max(df_plot$lat_p)))  
  # ggtitle(paste0('Correlation Coefficient of ', unlist(strsplit(df_nms[i,2], "[_]"))[1], 
  #                ' & ', var_y, ' - ', month.abb[as.numeric(df_nms[i,4])], ' ', df_nms[i,3], ' bi-weekly period ', title_nm))


##
summary(plsrmod)
plot(RMSEP(plsrmod))
plot(plsrmod, ncomp = 3, asp = 1, line = TRUE) # because scaled?
plot(plsrmod, plottype = "scores", comps = 1:3)
explvar(plsrmod)
plot(plsrmod, "loadings", comps = 1:2, legendpos = "topleft")
predict(plsrmod, ncomp = 3, newdata = data.frame(Xverif))
RMSEP(plsrmod, newdata = Xverif)
plot(plsrmod, plottype = "coef", ncomp=1:2, legendpos = "bottomleft")
predplot(plsrmod, ncomp = Ncomp, newdata = data.frame(t(Xverif)))


# Yaux        = PredictPLSRCrossVal(Xtrain,Xverif,Ypredictand,length(xvalindex[[igroup]]),
#                                   Ncomp,Nsegments,0,Ncomp)
Yvar        = Yaux$Yvar
Xvar        = Yaux$explvar
Scores[,init,igroup,-xvalindex[[igroup]]] = t(Yaux$Strain)
Scores[,init,igroup, xvalindex[[igroup]]] = t(Yaux$Sverif)

# } # End loop over groups

# Close plot with loadings
#dev.off()


# coef_df <- array(NA,dim=c(selpred,Nx,Ny))
# for(ipred in 1:selpred){ # Start loop over selected predictors
#   
#   # get predictor in correct format
#   pred_i = array(0,dim=c(Nx,Ny))
#   for(ilat in 1:Ny){
#     for(ilon in 1:Nx){
#       pred_i[ilon,ilat] = comp1[(Nx * (ilat - 1) + ilon) + Nx * Ny * (ipred - 1)]
#     }
#   }
#   coef_df[ipred,,] = pred_i
# } # End loop over predictors
