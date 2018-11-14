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
preds            = c(1,2,8,10) # q2m, sst, prate
cut_domain       = T
hru              ='1304'  #'1802'
wk               = '2_3' #'2_3'
pred_var         = 'prate'
pred_mon         = 1 #month to predict in CV
mon_nm           = 'DFJ'
Ntime_sub        = 30       # Define number of timesteps that each sub-period (i.e. validation) will have
Ncomp            = 2
flagPress        = 0 #flag to indicate if the PRESS statistic is used (1) or not (0) to predict. If zero, only one component is used
if (cut_domain) {
  file_in = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.cutdomain.rds')
  file_out         = paste0('CVresults_',hru,'_wk', wk, '_mon.', pred_mon,'.cutdomain.rds')
} else {
  file_in = paste0('plsr.input.mon_', mon_nm, '_lead_', wk,'.rds')
  file_out         = paste0('CVresults_',hru,'_wk', wk, '_mon.', pred_mon,'.cutdomain.rds')
}


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
Ngroups         = ceiling(Nval/Ntime_sub)
cutgroup        = sapply(1:Ngroups, function(j) c(rep(j,Ntime_sub)))[1:Nval] #not random
xvalindex       = split(date_val$ind, cutgroup) # validate on these indexes
Nsegments       = length(xvalindex)-1

##### ===== Arrange data for PLSR
#preds = c(1,2,4,6:10) #two variables arent in correct format yet
# preds = c(1,2,8,10) # q2m, sst, prate

Alldata <- ind_pred <- NULL
for(ipred in preds){ # Start loop over selected predictors

  # if hru, do different than gridded vars
  if (ipred == 11) {
    pred_i = array(cfs_pred_df$var, dim = c(Ntime,1))
  } else {
    
    Ny = length(na.omit(unique(grid[ipred, 1,])))
    Nx = length(na.omit(unique(grid[ipred, 2,])))
  
    pred_i = array(0,dim=c(Ntime,Nx*Ny))

    # for(ilat in 1:Ny){
    #   for(ilon in 1:Nx){
    #     pred_i[,(Nx * (ilat - 1) + ilon)]  = pred_cfs[ipred,ilon,ilat,]
    #   }
    # }
    for(ilon in 1:Nx){
      for(ilat in 1:Ny){
        pred_i[,(Ny * (ilon - 1) + ilat)]  = pred_cfs[ipred,ilon,ilat,]
      }
    }
  }
  
  # remove NAs (from SST or any other vars)
  t = which(!is.na(pred_i[1,]), TRUE)
  pred_i = pred_i[,t]

  # get beginning of data
  ind_pred = c(ind_pred, ncol(Alldata)) #last ind of var
  
  # For SST - get predictor in correct format
  if (var_nms[ipred] == 'sst') {
    df = na.omit(data.frame(y= grid[ipred ,1, ],  x = grid[ipred ,2,]))
    ## sst has NA's that need to be removed
    grid_sst = df[t,] #save grid
  }
  
  # bind data
  Alldata   <- cbind(Alldata, pred_i) 
} 

# scale - necessary?
ArrayCFSR        = scale(Alldata) #produces NAs
ind_pred         = c(0, ind_pred, ncol(Alldata))



##### ===== Predicat - loop over training/verification groups
df_cv <- NULL
load = array(NA, dim = c(Ngroups, length(ArrayCFSR[1,]), Ncomp))
beg.time = Sys.time()
for(igroup in 1:Ngroups){ # Start loop over groups
  
  ## Define predictand
  Ypredictand   = nldasVar[-xvalindex[[igroup]]]        # Observed nldasVar for the training period
  #Ypredictand   = nldasMod    # scaled above-Column vector containing the predictand to be used ([(Ngroups-1)*Ntime_sub] x 1)
  Yverif        = nldasVar[xvalindex[[igroup]]]     
  
  ## Extract predictors (separate training from verification dataset)
  Xtrain      = ArrayCFSR[-xvalindex[[igroup]],]
  Xverif      = ArrayCFSR[ xvalindex[[igroup]],]
  
  ## Fit PLSR model and predict! - from 'PredictPLSRCrossVal.R' function
  mydf        = data.frame(Ypredictand)
  mydf$X      = Xtrain
  
  ## perform PLSR - with no NAs (SST)
  # apply PLSR using function from 'PredictPLSRCrossVal.R' function
  # plsrmod     = plsr(Ypredictand ~ X, ncomp = Ncomp, data = mydf) #, method='simpls', validation = 'CV'
  plsrmod     = plsr(Ypredictand ~ X, ncomp = Ncomp, data = mydf, 
                     model = TRUE, x = TRUE, y = TRUE, #validation = 'CV',
                     segments=Nsegments, segment.type='consecutive', method='simpls')

  ## perform PRESS
  # if (flagpress == 0){
  #   minpress = Ncomp
  # }
  # if (flagpress == 1){
  #   val         = crossval(plsrmod, segment.type='consecutive')  # Perform cross validation
  #   indval      = which(is.nan(val$validation$PRESS)==F)  # Take out all NaN's
  #   minpress    = which(val$validation$PRESS[indval] == min(val$validation$PRESS[indval]))
  # }
  #B           = coef(plsrmod,ncomp=Ncomp,intercept=TRUE)
  #Yfitman     = B[1,1,1] +  Xtrain %*% B[-1,1,1]

  load[igroup,,]        = loadings(plsrmod)[,1:Ncomp]
  # Strain      = scores(plsrmod)[,1:Ncomp]
  # Sverif      = predict(plsrmod, ncomp = Ncomp, Xverif, se.fit=T, type='scores')
  # Xvar        = explvar(plsrmod)
  # Ycumvar     = drop(R2(plsrmod, estimate='train', intercept='FALSE')$val)*100 #variance

  ## predict Y
  Yhat_ver    = predict(plsrmod, ncomp = Ncomp, Xverif, se.fit=T)
  #Yhat        = predict(plsrmod, ncomp = Ncomp, ArrayCFSR[ ,], se.fit=T)
  
  df_cv = rbind(df_cv, cbind(Yverif, Yhat_ver))
  
}
end.time = Sys.time() - beg.time

##### ===== organize data
df = data.frame(date = date_val$date,
                ind = 1:Nval,
                Y_plsr = df_cv[,2] * Ysd + Ymean, 
                Y_cfs = left_join(date_val, cfs_pred_df, by = 'date')$var,
                Y_nld = df_cv[,1] * Ysd + Ymean)
df$Y_plsr = ifelse(df$Y_plsr <0, 0, df$Y_plsr) # remove negatives
df$res_plsr = df$Y_plsr - df$Y_nld
df$res_cfs = df$Y_cfs - df$Y_nld

## replace sst grid
grid_n = array(NA, dim = c(2, length(grid[1,1,])))
grid_n[1, ] = c(grid_sst[,1], rep(NA, (length(grid[1,1,]) - nrow(grid_sst))))
grid_n[2, ] = c(grid_sst[,2], rep(NA, (length(grid[1,1,]) - nrow(grid_sst))))
grid[which(var_nms == 'sst'),,] = grid_n


## save results
# setwd(out_model_dir)
# data = list(df = df, load = load, run.time = end.time, vars = var_nms[preds], 
#             grid = grid)
# saveRDS(data, file = file_out)


##### ===== Verification
# setwd(out_model_dir)
# ls = readRDS(file_out)
# df = ls$df

## calculate residual and correlation for cfsv2 and model
mean(abs(df$res_cfs))
mean(abs(df$res_plsr))
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


#P  = scores(plsrmod)             # CCA X component canonical variables
# ie PC time series
#  cols=cvs (time dim by min(cols(X,Y)
#U  = loadings(plsrmod)	      # Eigen Vectors (ie EOFs)

# use Y variance explained instead of eigenvalues for weights
# W = plsrmod$Xvar   # eigenval. (X variance)




### === Plot gridded coefficients
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/ploLoadings.R')

## input
var_nm     = 'prt' # variable
comp_plot  = 1     # EOF
igroup = 1
#ld = ls$load[igroup,,]
ld = load[igroup,,]

## set bins
max = round(max(ld, abs(min(ld))), -1)
dif = round(2*(max)/10)
bin_var = seq(-max, max, dif)
#if (max(bin_var) < max) {bin_var = c(bin_var, max(bin_var)+dif)}
#bin_var = c(-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20)


### === plot loadings
plot_ls = list()
for (var_i in var_nms[preds]) {
  comp = ld[,comp_plot]
  var = which(var_nms[preds] == var_i)
  grid_p = (grid[which(var_nms == var_i), , ])
  beg = ind_pred[var] + 1
  end = ind_pred[(var+1)]
  # if (var_i == 'sst') {
  #   plot_grid = data.frame(x = grid_sst[,1], y = grid_sst[,2])
  # } else {
    plot_grid = na.omit(data.frame(x = grid_p[1,],
                                   y = grid_p[2,]))
  # }
  
  df_plot = cbind.data.frame(var = comp[(beg:end)], lat_p = plot_grid$x, 
                             lon_p = plot_grid$y)
  
  
  plot_ls[[var_i]] = plotLoadings(df_plot, var_i, #max = max, min = -max,
               bin_vec = bin_var,
               title_in = paste0('Loadings - ', var_i, ' - Component ', comp_plot))  
}

## plot
gridExtra::grid.arrange(plot_ls[[var_nms[preds][1]]] + theme(legend.position = 'none'), 
                        plot_ls[[var_nms[preds][2]]]  + theme(legend.position = 'none'), 
                        plot_ls[[var_nms[preds][3]]]  + theme(legend.position = 'none'),
                        plot_ls[[var_nms[preds][4]]]  + theme(legend.position = 'none'), nrow = 2)

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

