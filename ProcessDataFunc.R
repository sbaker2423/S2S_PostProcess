ReadProcInputData <- function(preds, pred_cfs, grid, var_nms, date_cfs, 
                              outformat = 'data.frame', extremes_nldas = F) {
  
  print(paste('processing', var_nms[preds]))
  
  ## check if analysis is lagged and lagged data is saved yet
  setwd(out_model_dir)
  if (l > 0) {
    if (!file.exists(file_lag_in)){
      stop("Run '1_PLSR_data_prep.R to created lagged data file.", call.=FALSE)
    } else {
      pred_cfs = readRDS(file_lag_in)
    }
  }
  
  ##### ===== Load NLDAS Data 
  if (!extremes_nldas) {
    setwd(TP_data)
    file_nld   = nc_open(paste0('nldas.1_2.', pred_var, '.nc')) 
    
    # Extract hru coordinates for NLDAS
    hru_v_nld   = ncvar_get(file_nld, "hru")  # hru
    if (hru == 'all') {
      hruInd = 1:length(hru_v_nld)
    } else {
      hruInd = which(hru_v_nld == hru)
    }
    
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
  
    } else {
    
    ## extremes analysis
    library(zoo)
    setwd('~/s2s/nldas/files')
    nroll = 3
    file_extm = paste0(paste0('NLDAS_', pred_var, '_3dayMax_1999-2011_', nroll, 'dayMax.rds'))
    
    if (!file.exists(file_extm)) {
      
      file_nld   = nc_open(paste0('NLDAS_', pred_var, '_daily_1999-2011.nc')) 
      hru_v_nld   = ncvar_get(file_nld, "hru")  # hru
      hruInd = 1:length(hru_v_nld)
      
      ## read time and find overlapping index
      time = ncvar_get(file_nld, "time")
      date_nld = as.POSIXct((time) * 86400, 
                            origin = '1980-01-01', tz = 'utc') #timeseries of day 1 of obs
      
      ## take 3 day average - max for bi-weekly period
      nld_df   = ncvar_get(file_nld, pred_var)[hruInd, ]
      if (pred_var == 'prate') { nld_df = nld_df * 24 } ## convert mm/hr to mm/d
      l = zoo::zoo(t(nld_df), date_nld)
      nld_3mean = zoo::rollmean(l, nroll)
      nld_3dMax = zoo::rollmax(nld_3mean, 14)
      
      ## change date to fcst date (1)
      date_nld = date_nld[1:(length(date_nld)- (14 + nroll -2))]
      
      ## arrange data
      test = as.data.frame(nld_3dMax)
      colnames(test) <- hru_v_nld
      nld_df = cbind.data.frame(date = date_nld, test)
      
      ## save file and close nc
      nc_close(file_nld)
      saveRDS(nld_df, file = file_extm)
      
    } else {
      nld_df = readRDS(file_extm)
    }
    
    ## organize NLDAS df
    date_nld = nld_df$date
    nld_df = nld_df[, 2:ncol(nld_df)]
    hru_v_nld   = colnames(nld_df)
    
    ## Extract hru coordinates for NLDAS
    if (hru == 'all') {
      hruInd = 1:length(hru_v_nld)
    } else {
      hruInd = which(hru_v_nld == hru)
      hru_v_nld = hru
    }
    
    ind_overlap = inner_join(data.frame(date = date_cfs), 
                             data.frame(date = date_nld, ind = 1:length(date_nld)), by = 'date')
    
    ## read variable and convert units
    nld_df   = nld_df[ind_overlap$ind, hruInd]
  }
  
  ##### ===== Load CFSv2 T & P prediction on HRU
  setwd(TP_data)
  file_pred   = nc_open(paste0('cfsv2.', wk, '.', pred_var, '.nc')) 
  time_pred = ncvar_get(file_pred, "time")
  hru_v_cfs   = ncvar_get(file_pred, "hru")  # hru
  if (hru == 'all') {
    hruInd = 1:length(hru_v_cfs)
  } else {
    hruInd = which(hru_v_cfs == hru)
  }
  
  ## get date and mearge
  date_pred = as.POSIXct(as.Date(as.POSIXct(time_pred - 6*86400, 
                                            origin = '1970-01-01', tz = 'utc'))) #timeseries of day 1 in cst
  df_date = data.frame(date = date_pred, ind = 1:length(date_pred))
  ind_overlap = left_join(ind_overlap, df_date, by = 'date')
  cfs_pred_df   = ncvar_get(file_pred, pred_var)[hruInd, ind_overlap$ind.y]
  if (pred_var == 'prate') { cfs_pred_df = cfs_pred_df * 86400 } ## convert kg/m^2/s to mm/d
  cfs_pred_df = data.frame(date = ind_overlap$date, var = t(cfs_pred_df))
  cfs_pred_df = cfs_pred_df %>% dplyr::group_by(date) %>% 
    summarise_at(c(2:ncol(cfs_pred_df)), funs(mean(., na.rm=TRUE))) # daily average
  # cfs_pred_df <- aggregate(cfs_pred_df[,2:ncol(cfs_pred_df)], by =list(cfs_pred_df[,1]), FUN = mean, na.rm = TRUE)
  
  ## Scale NLDAS variable if only one hru considered
  if (hru != 'all') {
    nldasVar      = scale(nld_df)
    Ymean         = mean(nld_df)     # Mean nldasVar during training period
    Ysd           = sd(nld_df)       # Standard deviation of nldasVar during training period
  } else { # no scaling
    nldasVar      = nld_df
    Ymean <- Ysd <- NA
  }
  
  
  ##### =====  get groups for validation
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
  print('arranging data...')
  Alldata <- ind_pred <- NULL
  if (outformat == 'list') { Alldata <- list() }
  
  for(ipred in preds){ # Start loop over selected predictors
    
    # if hru, do different than gridded vars
    if (var_nms[ipred] == 'hru') {
      pred_i = array(cfs_pred_df$var, dim = c(Ntime,1))
    } else {
      
      Ny = length(na.omit(unique(grid[ipred, 1,])))
      Nx = length(na.omit(unique(grid[ipred, 2,])))
      
      pred_i = array(0,dim=c(Ntime,Nx*Ny))
      
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
    if (outformat == 'data.frame') {
      ind_pred = c(ind_pred, ncol(Alldata)) #last ind of var
    }
    
    # For SST - get predictor in correct format (sst has NA's that need to be removed)
    if (var_nms[ipred] == 'sst') {
      df = na.omit(data.frame(y= grid[ipred ,1, ],  x = grid[ipred ,2,]))
      grid_sst = df[t,] #save grid
    }
    
    # collect data
    if (outformat == 'list') {
      Alldata[[ipred]]   <- scale(pred_i)
    } else if (outformat == 'data.frame') {
      Alldata   <- cbind(Alldata, pred_i) 
    } else { stop('Output format must be specified as data.frame or list') }
  } 
  
  ## scale output data for use in plsr
  if (outformat == 'data.frame') {
    Alldata        = scale(Alldata) #produces NAs if NAs
    ind_pred       = c(0, ind_pred, ncol(Alldata))
  }
  
  
  ## replace sst grid
  if (any(var_nms[preds] == 'sst')) {
    grid_n = array(NA, dim = c(2, length(grid[1,1,])))
    grid_n[1, ] = c(grid_sst[,1], rep(NA, (length(grid[1,1,]) - nrow(grid_sst))))
    grid_n[2, ] = c(grid_sst[,2], rep(NA, (length(grid[1,1,]) - nrow(grid_sst))))
    grid[which(var_nms == 'sst'),,] = grid_n
  }
  
  rlist <- list(Alldata = Alldata, ind_pred = ind_pred, nldasVar = nldasVar, 
                Ymean = Ymean, Ysd = Ysd, date_v = date_v, cfs_pred_df = cfs_pred_df,
                grid = grid, hru_v_nld = hru_v_nld, hru_v_cfs = hru_v_cfs,
                preds = preds)
  
  return(rlist)
}




