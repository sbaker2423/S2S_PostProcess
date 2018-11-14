# ================================================================================ 
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
TP_data = '/home/sabaker/s2s/analysis/files/cfsv2_files/2wk_avg/'
out_model_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/modeledVars/'

## --- Input Data
wk = '2_3'
avg_fcst_day = 7
hru = '1304'
mon_in           = 1 # month to analyze
Ntime_sub        = 14       # Define number of timesteps that each sub-period (i.e. validation) will have


## --- libraries
library(lubridate)
library(ncdf4)
library(abind)
library(OceanView)
library(pls)
library(dplyr)

# Define long variable names for title of plots
# ReanalNames   = c('Geopotential Height 500mb')
# ReanalShort   = c('z500')


##### ===== Part 1: Load Predictors & Predictanad ===== #####

## Input number of predictors
Npred = 8 

## Predictand - Open netCDF files
setwd(TP_data)
file_nldprate   = nc_open(paste0('nldas.1_2.prate.nc')) 
file_nldtmp2m   = nc_open(paste0('nldas.1_2.tmp2m.nc')) 

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

# Extract lat/lon coordinates - different for variables
latitude       = ncvar_get(file_hgt, "latitude" ) 
longitude      = ncvar_get(file_hgt, "longitude" )
latitude_odd_in   = ncvar_get(file_q2m, "latitude" ) 
longitude_odd_in  = ncvar_get(file_q2m, "longitude" ) 
latitude_sst_in   = ncvar_get(file_sst, "latitude" ) 
longitude_sst_in  = ncvar_get(file_sst, "longitude" ) 
latitude_PT_in   = ncvar_get(file_prt, "lat" )
longitude_PT_in  = ncvar_get(file_prt, "lon" )
# latitude_odd2_in   = ncvar_get(file_pwt, "latitude" ) 
# longitude_odd2_in  = ncvar_get(file_pwt, "longitude" ) 

# Extract lon/lat indices (does nothing now, could be useful later)
Ny         = length(latitude)
Nx         = length(longitude)


# Extract hru coordinates for NLDAS
hru_v   = ncvar_get(file_nldprate, "hru")  # hru
hruInd = which(hru_v == hru)

# Times - put into list
# var_ls = c ('nld.prate', 'nld.tmp2m', 'hgt', 'q2m', 'prs', 'slp', 'pwt', 'uwnd', 'vwnd', 'sst')
var_ls = c ('nld.prate', 'nld.tmp2m', 'hgt', 'q2m', 'prs', 'slp', 'pwt', 'uwnd', 'vwnd', 'sst', 'tmp', 'prt')
time = list()
time[[var_ls[1]]] = ncvar_get(file_nldprate, "time")
time[[var_ls[2]]] = ncvar_get(file_nldprate, "time")
time[[var_ls[3]]] = ncvar_get(file_hgt, "time") 
time[[var_ls[4]]] = ncvar_get(file_q2m, "time") 
time[[var_ls[5]]] = ncvar_get(file_prs, "time") 
time[[var_ls[6]]] = ncvar_get(file_slp, "time") 
time[[var_ls[7]]] = ncvar_get(file_pwt, "time") 
time[[var_ls[8]]] = ncvar_get(file_uwnd, "time") 
time[[var_ls[9]]] = ncvar_get(file_vwnd, "time") 
time[[var_ls[10]]] = ncvar_get(file_sst, "time") 
time[[var_ls[11]]] = ncvar_get(file_tmp, "time")
time[[var_ls[12]]] = ncvar_get(file_prt, "time")

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
indTime = filter(ind_overlap, month(ind_overlap$date) == mon_in)
Ntime = nrow(indTime)

## Read Variables
nld.prt_df   = ncvar_get(file_nldprate, "prate")[hruInd, indTime$ind.nld] * 24 ## convert mm/hr to mm/d
nld.tmp_df   = ncvar_get(file_nldtmp2m, "tmp2m")[hruInd, indTime$ind.nld] 
hgt_df   = ncvar_get(file_hgt, "HGT_500mb")[,,indTime$ind.hgt] 
q2m_df   = ncvar_get(file_q2m, "SPFH_2maboveground")[,,indTime$ind.q2m]
prs_df   = ncvar_get(file_prs, "PRES_surface")[,,indTime$ind.prs]
slp_df   = ncvar_get(file_slp, "PRMSL_meansealevel")[,,indTime$ind.slp]
pwt_df   = ncvar_get(file_pwt, "PWAT_entireatmosphere_consideredasasinglelayer_")[,,indTime$ind.pwt]
uwnd_df   = ncvar_get(file_uwnd, "UGRD_850mb")[,,indTime$ind.uwnd]
vwnd_df   = ncvar_get(file_vwnd, "VGRD_850mb")[,,indTime$ind.vwnd]
sst_df   = ncvar_get(file_sst, "POT_5mbelowsealevel")[,,indTime$ind.sst]
tmp_df   = ncvar_get(file_tmp, "TMP_2maboveground")[,,indTime$ind.tmp]
prt_df   = ncvar_get(file_prt, "PRATE_surface")[,,indTime$ind.prt]

## close netCDFs
nc_close(file_nldprate); nc_close(file_hgt); nc_close(file_q2m); nc_close(file_prs); nc_close(file_slp)
nc_close(file_nldtmp2m); nc_close(file_pwt); nc_close(file_uwnd); nc_close(file_vwnd); nc_close(file_sst)


## === Remap to match other coordinates
## SST - coarsen grid
longitude_sst = longitude[2:Nx]; latitude_sst = latitude[2:Ny]
sst = remap(sst_df, x=longitude_sst_in, y=latitude_sst_in, z=1:Ntime,
            xto=longitude_sst, yto=latitude_sst, zto=1:Ntime)$var
## Precip and Temp - coarsen grid & smaller domain
longitude_PT = longitude[117:236]; latitude_PT = latitude[68:136]
tmp = remap(tmp_df, x=longitude_PT_in, y=latitude_PT_in, z=1:Ntime,
            xto=longitude_PT, yto=latitude_PT, zto=1:Ntime)$var
prt = remap(prt_df, x=longitude_PT_in, y=latitude_PT_in, z=1:Ntime,
            xto=longitude_PT, yto=latitude_PT, zto=1:Ntime)$var
## q2m, (pwat, pressfc)????? - change grid resolution
# longitude_odd = longitude[117:236]; latitude_odd = latitude[68:136]
q2m = remap(q2m_df, x=longitude_odd_in, y=latitude_odd_in, z=1:Ntime,
            xto=longitude, yto=latitude, zto=1:Ntime)$var

# Extract lon/lat indices - new indices
Ny_odd     = length(latitude_odd)
Nx_odd     = length(longitude_odd)
Ny_sst     = length(latitude_sst)
Nx_sst     = length(longitude_sst)
Ny_PT     = length(latitude_PT)
Nx_PT     = length(longitude_PT)






## --- Reorganize predictor data as matrix with Ntime
pred_cfs = array(NA,dim=c(length(var_ls),Nx,Ny,Ntime))
## NLDAS hru variables
for(i in 1:hru){
  pred_cfs[3,ilon,ilat,]  = hgt_df[ilon,ilat,]
}
## normal CFSv2 variables
for(ilat in 1:Ny){
  for(ilon in 1:Nx){
    pred_cfs[3,ilon,ilat,]  = hgt_df[ilon,ilat,]
    pred_cfs[4,ilon,ilat,]  = q2m_df[ilon,ilat,]
    # pred_cfs[5,ilon,ilat,]  = prs_df[ilon,ilat,]
    pred_cfs[6,ilon,ilat,]  = slp_df[ilon,ilat,]
    # pred_cfs[7,ilon,ilat,]  = pwt_df[ilon,ilat,]
    pred_cfs[8,ilon,ilat,]  = uwnd_df[ilon,ilat,]
    pred_cfs[9,ilon,ilat,]  = vwnd_df[ilon,ilat,]
  }
}
## sst grid
for(ilat in 1:Ny_sst){
  for(ilon in 1:Nx_sst){
    pred_cfs[10,ilon,ilat,]  = sst[ilon,ilat,]
  }
}
## precip & temp grid
for(ilat in 1:Ny_PT){
  for(ilon in 1:Nx_PT){
    pred_cfs[11,ilon,ilat,]  = tmp[ilon,ilat,]
    pred_cfs[12,ilon,ilat,]  = prt[ilon,ilat,]
  }
}

## store grid
grid <- matrix(NA, nrow = Nx*Ny, ncol = 6)
for(ilat in 1:Ny){
  for(ilon in 1:Nx){
    grid[(Nx * (ilat - 1) + ilon),1] = latitude[ilat]
    grid[(Nx * (ilat - 1) + ilon),2] = longitude[ilon]
    grid[(Nx * (ilat - 1) + ilon),3] = latitude_sst[ilat]
    grid[(Nx * (ilat - 1) + ilon),4] = longitude_sst[ilon]
    grid[(Nx * (ilat - 1) + ilon),5] = latitude_PT[ilat]
    grid[(Nx * (ilat - 1) + ilon),6] = longitude_PT[ilon]
    #to undo: pred_i[,(Nx * (ilat - 1) + ilon)]  = pred_cfs[ipred,ilon,ilat,]
  }
}



test = list(pred_cfs, grid)



## Remove raw data matrices
rm(hgt_df, q2m_df, prs_df, slp_df, pwt_df, uwnd_df, vwnd_df, sst_df, sst) 

##### ===== Part 3: pre-process ===== #####
Ncomp            = 4  # number of components





#### ? get groups for validation?
Ngroups         = ceiling(Ntime/Ntime_sub)
cutgroup        = sapply(1:Ngroups, function(j) c(rep(j,Ntime_sub)))[1:Ntime] #not random
xvalindex       = split(1:Ntime, cutgroup) # validate on these indexes
Nsegments       = length(xvalindex)-1

# Scale NLDAS variable
nldasVar      = scale(nld_df)
Ymean         = mean(nld_df)     # Mean nldasVar during training period
Ysd           = sd(nld_df)       # Standard deviation of nldasVar during training period


##### ===== Part 4: Do PLSR over the entire field ===== #####

selpred  = nrow(pred_cfs)
# Scores   = array(NA,dim=c(Ncomp,Nlead,Ngroups,Ntime))

## loop - arrange predictors
Alldata <- NULL
for(ipred in 1:8){ # Start loop over selected predictors
  
  pred_i = array(0,dim=c(Ntime,Nx*Ny))
    for(ilat in 1:Ny){
      for(ilon in 1:Nx){
        pred_i[,(Nx * (ilat - 1) + ilon)]  = pred_cfs[ipred,ilon,ilat,]
      }
    }
  
  # For SST - get predictor in correct format
  if (ipred == 8) {
    ## sst has NA's that need to be removed
    t = which(!is.na(pred_i[1,]), TRUE)
    pred_i = pred_i[,t]
    grid_sst = grid[t,c(1:2)] #save grid
  }

  # bind data
  Alldata   <- cbind(Alldata, pred_i) 
} 

### need to remove NAss



# scale - necessary?
ArrayCFSR        = scale(Alldata) #produces NAs

## loop over training/verification groups
igroup =1
# for(igroup in 1:Ngroups){ # Start loop over groups

# Define predictand
nldasMod      = nldasVar[-xvalindex[[igroup]]]        # Observed nldasVar for the training period
Ypredictand   = nldasMod    # scaled above-Column vector containing the predictand to be used ([(Ngroups-1)*Ntime_sub] x 1)
#Ymean         = mean(nldasMod)     # Mean nldasVar during training period
#Ysd           = sd(nldasMod)       # Standard deviation of nldasVar during training period
#Yfull         = nldasVar
Yverif        = nldasVar[xvalindex[[igroup]]]     

# Extract predictors (separate training from verification dataset)
Xtrain      = ArrayCFSR[-xvalindex[[igroup]],]
Xverif      = ArrayCFSR[ xvalindex[[igroup]],]
#Xfull       = ArrayCFSR[ ,]

# Fit PLSR model and predict!
# from from 'PredictPLSRCrossVal.R' function
mydf        = data.frame(Ypredictand)
mydf$X      = Xtrain

# perform PLSR - with no NAs (SST)
plsrmod     = plsr(Ypredictand ~ X, ncomp = Ncomp, data = mydf, method='simpls', validation = 'CV')
# save
setwd(out_model_dir)
save(plsrmod, file = paste0('prate.',hru,'.plsr_model.rda'))
test = load(paste0('prate.',hru,'.plsr_model.rda'))

## predict Y
Yhat_ver    = predict(plsrmod, ncomp = Ncomp, Xverif, se.fit=T)
Yhat        = predict(plsrmod, ncomp = Ncomp, ArrayCFSR[ ,], se.fit=T)

## example does cross validate on training data...(vinette)

res = Yhat[,1,1] - nldasVar # residulas
mean(abs(res))
plot(nldasVar, Yhat) # negative values?
abline(0,1)

## unscaled data
# Yhat_unscaled = Yhat * attr(nldasVar, 'scaled:scale') + attr(nldasVar, 'scaled:center')
Yhat_unscaled = Yhat * Ysd + Ymean
Yhat_unscaled = ifelse(Yhat_unscaled <0, 0, Yhat_unscaled) # remove negatives
res_unscaled = Yhat_unscaled - nld_df
plot(nld_df, Yhat_unscaled); abline(0,1)


P  = scores(plsrmod)             # CCA X component canonical variables
# ie PC time series
#  cols=cvs (time dim by min(cols(X,Y)
U  = loadings(plsrmod)	      # Eigen Vectors (ie EOFs)

# use Y variance explained instead of eigenvalues for weights
# W = plsrmod$Xvar   # eigenval. (X variance)


##

# comp1_z500 = data.frame(comp1 *  attr(ArrayCFSR, 'scaled:scale') + 
#                           attr(ArrayCFSR, 'scaled:center'))[1:(Nx*Ny),]
# comp1_x = data.frame(comp1 *  attr(ArrayCFSR, 'scaled:scale') + 
#                        attr(ArrayCFSR, 'scaled:center'))[(Nx*Ny+1):2*(Nx*Ny),]
# df_plot = cbind.data.frame(comp1_z500, comp1_x, grid)
# df_plot = cbind.data.frame(z500 = comp1_unscaled[1:(Nx*Ny),],
#                            q2m = comp1_unscaled[((Nx*Ny+1):(2*(Nx*Ny))),],
#                            grid)

### === Plot gridded coefficients
var_nm     = 'sst' # variable
comp_plot  = 1     # EOF
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/ploLoadings.R')
# set range
max = max(plsrmod$loadings)
min = min(plsrmod$loadings)

## plot coefficients
# comp = plsrmod$coefficients[,1,comp_plot]
# comp = data.frame(comp *  attr(ArrayCFSR, 'scaled:scale') + 
#                              attr(ArrayCFSR, 'scaled:center'))
# df_plot = cbind.data.frame(var = comp[(((var-1)*(Nx*Ny)+1):(var*(Nx*Ny))),], grid)

### === plot loadings
comp = plsrmod$loadings[,comp_plot]
var = which(var_ls == var_nm) - 1
beg = ((var-1)*(Nx*Ny)+1)
end = (var*(Nx*Ny))
if (var_nm == 'sst') {
  end = length(comp)
  df_plot = cbind.data.frame(var = comp[(beg:end)], lat_p = grid_sst[,1], 
                             lon_p = grid_sst[,2])
} else if (any(var_nm == c('q2m', 'prs', 'pwt'))) {
  df_plot = cbind.data.frame(var = comp[(beg:end)], lat_p = grid$lat_odd, 
                             lon_p = grid$lon_odd)
} else {
  df_plot = cbind.data.frame(var = comp[(beg:end)], lat_p = grid$lat, 
                             lon_p = grid$lon)
}

plotLoadings(df_plot, var_nm, max = max, min = min,
             bin_vec = c(-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20))



# # create bins of equal widths
# scale_max = 10^(ceiling(log10(max)))/10
# min_bin = floor(min/scale_max) *scale_max
# max_bin = ceiling(max/scale_max) *scale_max
# dif = round((max_bin - min_bin)/10)
# bin_var = seq(min_bin, max_bin, dif)
# if (max(bin_var) < max_bin) {bin_var = c(bin_var, max(bin_var)+dif)}
# if(max_bin > 1000) {dig= 4} else {dig=3} # set sig digs
# # over write bins
# bin_var = c(-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20)
# 
# # create color palette
# cust_pal <- rev(brewer.pal(n=length(bin_var)-1, name="RdYlBu"))
# 
# # add bins to df
# df_plot$bins = cut(df_plot$var, breaks = bin_var, dig.lab = dig)

#labs = levels(df_plot$bins)

# ## plot with world map and discrete colors
# ggplot() + geom_raster(data = df_plot, aes(x = lon_p, y = lat_p, fill = bins)) +
#   scale_fill_manual(values = cust_pal, name = 'Loadings', drop = F, labels = labs,
#                     guide = guide_legend(ncol = 1, label.position = 'right', label.hjust = 0,
#                                          title.position = 'top', title.hjust = 0.5, reverse = T,
#                                          keywidth = 1, keyheight = 1)) +
#   theme_bw() +
#   coord_equal() + # ratio of x and y...coord equal preserves the relative sizes of each
#   geom_polygon(data = worldmap, aes(X,Y,group=PID), color = 'black', size = 0.4, fill = NA) +
#   xlab('') +  ylab('') +
#   scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0), limits = c(min(df_plot$lon_p), max(df_plot$lon_p))) +
#   scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0), limits = c(min(df_plot$lat_p), max(df_plot$lat_p)))  
  # ggtitle(paste0('Correlation Coefficient of ', unlist(strsplit(df_nms[i,2], "[_]"))[1], 
  #                ' & ', var_y, ' - ', month.abb[as.numeric(df_nms[i,4])], ' ', df_nms[i,3], ' bi-weekly period ', title_nm))


# ## bins - z500
# dif = round((max(df_plot$var) - min(df_plot$var))/12, -1)
# min_bin = floor(min(df_plot$var)/10) *10
# max_bin = ceiling(max(df_plot$var)/10) *10
# dif = round(ceiling((max_bin - min_bin)*1/11), -1)
# bin_var = seq(min_bin, max_bin, dif)
# if (max(bin_var) < max_bin) {bin_var = c(bin_var, max(bin_var)+dif)}
# cust_pal_var <- rev(brewer.pal(n=length(bin_var), name="RdYlBu"))
# ## set up bins and color palette
# if(max_bin > 1000) {dig= 4} else {dig=3}
# df_plot$bins = cut(df_plot$var, breaks = bin_var, dig.lab = dig)
# cust_pal = cust_pal_var
# # set lat and lon for figure
# df_plot$lat_p = df_plot$lat
# df_plot$lon_p = df_plot$lon
# 
# ## bins - q2m
# scale_i = 1/10000
# dif = round((max(df_plot$q2m) - min(df_plot$q2m))/11, 3)
# min_bin = floor(min(df_plot$q2m)/scale_i) *scale_i
# max_bin = ceiling(max(df_plot$q2m)/scale_i) *scale_i
# #dif = round(ceiling((max_bin - min_bin)*1/11), -1)
# bin_q2m = seq(min_bin, max_bin, dif)
# if (max(bin_q2m) < max_bin) {bin_q2m = c(bin_q2m, max(bin_q2m)+dif)}
# cust_pal_q2m <- rev(brewer.pal(n=length(bin_q2m), name="RdYlBu"))
# ## set up bins and color palette
# if(max_bin > 1000) {dig= 4} else {dig=3}
# df_plot$bins = cut(df_plot$q2m, breaks = bin_q2m, dig.lab = dig)
# cust_pal = cust_pal_q2m
# # set lat and lon for figure
# df_plot$lat_p = df_plot$lat_odd
# df_plot$lon_p = df_plot$lon_odd

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
