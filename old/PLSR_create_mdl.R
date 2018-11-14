# ================================================================================ 
# PLSR code for processing CFSv2 variables for predicting 
# NLDAS temperature and precipitation for 202 HUCs in CONUS domain
#   Created by S. Baker 
#
#
#
#   Adapted from P. Mendoza script (May 2016) and Andy Wood, NCAR
#   from PLSR R code (A. Wood) for cross-validated application of PCR/PLSR
# ================================================================================ 
rm(list=ls())

## --- Directories
TP_data = '/home/sabaker/s2s/analysis/files/cfsv2_files/2wk_avg/'
out_model_dir = '/home/sabaker/s2s/analysis/files/cfsv2_files/plsr_input/'

## --- Input Data
hru              = '1802'
wk               = '3_4' #'2_3'
pred_var         = 'prate'
mon_in           = 1 # month to analyze
Ntime_sub        = 14       # Define number of timesteps that each sub-period (i.e. validation) will have


## --- libraries
library(lubridate)
library(ncdf4)
library(abind)
library(OceanView)
library(pls)
library(dplyr)
library(ggplot2)

# Define long variable names for title of plots
# ReanalNames   = c('Geopotential Height 500mb')
# ReanalShort   = c('z500')


##### ===== Load CFSv2 Data 
setwd(out_model_dir)
load(file = paste0('plsr.input.mon_', mon_in, '_lead_', wk,'.rds'))
nms = names(data_out)

## extract objects
pred_cfs = data_out[[nms[1]]]
grid = data_out[[nms[2]]]
var_nms = data_out[[nms[3]]]
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
nld_df   = ncvar_get(file_nld, pred_var)[hruInd, ind_overlap$ind]
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

##### =====  get groups for validation?
Ntime           = length(pred_cfs[1,1,1,])
Ngroups         = ceiling(Ntime/Ntime_sub)
cutgroup        = sapply(1:Ngroups, function(j) c(rep(j,Ntime_sub)))[1:Ntime] #not random
xvalindex       = split(1:Ntime, cutgroup) # validate on these indexes
Nsegments       = length(xvalindex)-1

## Scale NLDAS variable
nldasVar      = scale(nld_df)
Ymean         = mean(nld_df)     # Mean nldasVar during training period
Ysd           = sd(nld_df)       # Standard deviation of nldasVar during training period


##### ===== Arrange data for PLSR
# Scores   = array(NA,dim=c(Ncomp,Nlead,Ngroups,Ntime))

## loop - arrange predictors
preds = c(1,2,4,6:10) #two variables arent in correct format yet
grid_in = c(1,1,1,1,1,1,1,2,3,3) #types of grids (based on native grid)
# var_i = var_nms[preds]
Alldata <- ind_pred <- NULL
for(ipred in preds){ # Start loop over selected predictors
  
  Ny = length(na.omit(unique(grid[ ,(2 * grid_in[ipred] - 1) ])))
  Nx = length(na.omit(unique(grid[ ,(2 * grid_in[ipred]) ])))
  
  
  pred_i = array(0,dim=c(Ntime,Nx*Ny))
  for(ilat in 1:Ny){
    for(ilon in 1:Nx){
      pred_i[,(Nx * (ilat - 1) + ilon)]  = pred_cfs[ipred,ilon,ilat,]
    }
  }
  
  # remove NAs (from SST or any other vars)
  t = which(!is.na(pred_i[1,]), TRUE)
  pred_i = pred_i[,t]

  # get beginning of data
  ind_pred = c(ind_pred, ncol(Alldata)) #last ind of var
  
  # For SST - get predictor in correct format
  if (var_nms[ipred] == 'sst') {
    df = na.omit(data.frame(x = grid[ ,3 ],  y = grid[ ,4]))
    ## sst has NA's that need to be removed
    grid_sst = df[t,] #save grid
  }
  
  # bind data
  Alldata   <- cbind(Alldata, pred_i) 
} 
#rm(pred_cfs)

# scale - necessary?
ArrayCFSR        = scale(Alldata) #produces NAs

##### ===== Predicat - loop over training/verification groups
igroup =1
# for(igroup in 1:Ngroups){ # Start loop over groups

## Define predictand
nldasMod      = nldasVar[-xvalindex[[igroup]]]        # Observed nldasVar for the training period
Ypredictand   = nldasMod    # scaled above-Column vector containing the predictand to be used ([(Ngroups-1)*Ntime_sub] x 1)
Yverif        = nldasVar[xvalindex[[igroup]]]     

## Extract predictors (separate training from verification dataset)
Xtrain      = ArrayCFSR[-xvalindex[[igroup]],]
Xverif      = ArrayCFSR[ xvalindex[[igroup]],]
#Xfull       = ArrayCFSR[ ,]

## Fit PLSR model and predict! - from 'PredictPLSRCrossVal.R' function
mydf        = data.frame(Ypredictand)
mydf$X      = Xtrain

# perform PLSR - with no NAs (SST)
Ncomp = 3
plsrmod     = plsr(Ypredictand ~ X, ncomp = Ncomp, data = mydf, method='simpls', validation = 'CV')


## save
# setwd(out_model_dir)
# save(plsrmod, file = paste0('prate.',hru,'.plsr_model.rda'))
# test = load(paste0('prate.',hru,'.plsr_model.rda'))

## predict Y
Yhat_ver    = predict(plsrmod, ncomp = Ncomp, Xverif, se.fit=T)
Yhat        = predict(plsrmod, ncomp = Ncomp, ArrayCFSR[ ,], se.fit=T)




##### ===== organize data

## Yhat - unscaled data
# Yhat_unscaled = Yhat * attr(nldasVar, 'scaled:scale') + attr(nldasVar, 'scaled:center')
Yhat_unscaled = Yhat * Ysd + Ymean
Yhat_unscaled = ifelse(Yhat_unscaled <0, 0, Yhat_unscaled) # remove negatives
res_unscaled = Yhat_unscaled - nld_df

## CFSv2 predicted Y
Y_cfs = cfs_pred_df$var
res_cfs = Y_cfs - nld_df
df = data.frame(date = cfs_pred_df$date, ind = 1:Ntime,
                Y_plsr = Yhat_unscaled[,1,1], Y_cfs = Y_cfs,
                Y_nld = nld_df)

##### ===== Verification
## calculate residual and correlation for cfsv2 and model
mean(abs(res_cfs))
mean(abs(res_unscaled))
cor(df$Y_nld, df$Y_cfs)
cor(df$Y_nld, df$Y_plsr)

## 1:1 plot
par(mfrow=c(1,2))
ymax = max(df$Y_nld, df$Y_cfs, df$Y_plsr)
plot(df$Y_nld, df$Y_plsr, xlim = c(0,ymax), ylim = c(0,ymax)); abline(0,1) 
plot(df$Y_nld, df$Y_cfs, xlim = c(0,ymax), ylim = c(0,ymax)); abline(0,1) 

## plot timeseries
ggplot(df, aes(x = ind, y = Y_nld)) +
  geom_line() +
  geom_line(aes(x = ind, y = Y_cfs), color = 'blue', linetype = 'dotted') +
  geom_line(aes(x = ind, y = Y_plsr), color = 'red')


P  = scores(plsrmod)             # CCA X component canonical variables
# ie PC time series
#  cols=cvs (time dim by min(cols(X,Y)
#U  = loadings(plsrmod)	      # Eigen Vectors (ie EOFs)

# use Y variance explained instead of eigenvalues for weights
W = plsrmod$Xvar   # eigenval. (X variance)




### === Plot gridded coefficients
ind_pred = c(0, ind_pred, ncol(Alldata))
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/ploLoadings.R')
# set range
max = max(plsrmod$loadings)
min = min(plsrmod$loadings)

### === plot loadings
var_nm     = 'prt' # variable
comp_plot  = 2     # EOF
plot_ls = list()
for (var_nm in var_nms[preds]) {
  comp = plsrmod$loadings[,comp_plot]
  var = which(var_nms[preds] == var_nm)
  grid_p = grid_in[preds]
  beg = ind_pred[var] + 1
  end = ind_pred[(var+1)]
  if (var_nm == 'sst') {
    plot_grid = data.frame(x = grid_sst[,1], y = grid_sst[,2])
  } else {
    plot_grid = na.omit(data.frame(x = grid[ ,(2 * grid_p[var] - 1) ],
                                 y = grid[ ,(2 * grid_p[var]) ]))
  }
  
  df_plot = cbind.data.frame(var = comp[(beg:end)], lat_p = plot_grid$x, 
                             lon_p = plot_grid$y)
  
  
  plot_ls[[var_nm]] = plotLoadings(df_plot, var_nm, max = max, min = min,
               bin_vec = c(-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20),
               title_in = paste0('Loadings - ', var_nm, ' - Component ', comp_plot))  
}

#test = var_nms[preds]
gridExtra::grid.arrange(plot_ls[['hgt']] + theme(legend.position = 'none'), 
                        plot_ls[['q2m']]  + theme(legend.position = 'none'), 
                        plot_ls[['slp']]  + theme(legend.position = 'none'), 
                        plot_ls[['uwnd']]  + theme(legend.position = 'none'), nrow = 2)
gridExtra::grid.arrange(plot_ls[['vwnd']] + theme(legend.position = 'none'), 
                        plot_ls[['sst']]  + theme(legend.position = 'none'), 
                        plot_ls[['tmp']]  + theme(legend.position = 'none'), 
                        plot_ls[['prt']] + theme(legend.position = 'none'), nrow = 2)
plot_ls[['hgt']]


## summary
summary(plsrmod)
plot(RMSEP(plsrmod))
plot(plsrmod, ncomp = 3, asp = 1, line = TRUE) # because scaled?
plot(plsrmod, plottype = "scores", comps = 1:3)
explvar(plsrmod)
#plot(plsrmod, "loadings", comps = 1:2, legendpos = "topleft")
#predict(plsrmod, ncomp = 3, newdata = data.frame(Xverif))
#RMSEP(plsrmod, newdata = Xverif)
#plot(plsrmod, plottype = "coef", ncomp=1:2, legendpos = "bottomleft")
#predplot(plsrmod, ncomp = Ncomp, newdata = data.frame(t(Xverif)))


# Yaux        = PredictPLSRCrossVal(Xtrain,Xverif,Ypredictand,length(xvalindex[[igroup]]),
#                                   Ncomp,Nsegments,0,Ncomp)
Yvar        = Yaux$Yvar
Xvar        = Yaux$explvar
Scores[,init,igroup,-xvalindex[[igroup]]] = t(Yaux$Strain)
Scores[,init,igroup, xvalindex[[igroup]]] = t(Yaux$Sverif)

# } # End loop over groups

# Close plot with loadings
#dev.off()

