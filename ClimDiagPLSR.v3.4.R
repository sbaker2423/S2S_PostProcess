##########  SEAS. FORECASTING WITH PARTIAL LEAST SQUARES REGRESSION (PLSR)  #########################

# Created by P. Mendoza (May 24, 2016) under the direction of Andy Wood, NCAR
#   adapts PLSR R code from A. Wood for cross-validated application of PCR/PLSR

# This script reads CFSR variable fields, computes antecedent 3-month seasonal averages 
#   for several initialization times and generates and ensemble seas. flow forecast prediction.

# Updates:
#   20160718 PM - now saves Ncomp score vectors and the matrix with Nyears x Npoints
#   20180328 AW - reorganizing code
# =====================================================================================

##########  PART 1: Load libraries and functions, set directories and define some variables #########

##### Part 1.1: Read control file and load libraries and custom functions #####

filename  = '/hydro_d2/pmendoza/overtheloop/Rscripts/Climate_Diagnostics/control_climdiag.R' ## what is this?
source(filename)

# --- libraries
library(lubridate)
library(ncdf4)
library(abind)
library(OceanView)
library(pls)

# --- custom R functions
FcnDir = paste(MasterDir,'/Rscripts/R_functions/',sep='') # Location of function code
source(paste(FcnDir, 'PlotMapRegGrid.R',sep=''))           # Plot correlation maps
source(paste(FcnDir, 'PlotTwoClustersMap.R',sep=''))       # Plot two clusters
source(paste(FcnDir, 'CreateLandMask.R',sep=''))           # Returns matrix with land mask
source(paste(FcnDir, 'PredictPLSRCrossVal.R',sep=''))      # Returns matrix with land mask
source(paste(FcnDir, 'BoxplotEnsForecastsv2.R', sep = '')) # Boxplots w/ ensembles
source(paste(FcnDir, 'ScatterPlotsCorBias.R', sep = ''))   # Boxplots with ensembles


# --- Start loop over predictand periods
#for( mminit in 4:5){

#Allinit      <- rbind( c(4,7), c(4,9), c(5,7), c(5,5), c(5,9) )
#ipredictand  <- Allinit[mminit,]

#for(ibas in c(5,6,7,8,45)){# Start loop over basins
#ibas         = 6
islog        = 2
Mincorthresh = 0.4
Nens         = 500
Ncomp        = 2
Months2      = c('J','F','M','A','M','J','J','A','S','O','N','D')
MaxComp      = 2

##### Part 1.2: Define other useful variables, call libraries and functions, etc.

# Define size of cells for square to look around max/min correlation points
if (ReaCode      == 'NCEP_NCAR'){  Ncells = 4   }
if (ReaCode      == 'CFSR')     {  
  if(newdx == 1){  Ncells = 20 }
  if(newdx == 2){  Ncells = 10 }
}
# Number of predictor reanalysis variables
# (1) Surface Air temperature (SAT)
# (2) Geopotential height (Z700)
# (3) Precipitable Water (PW)
# (4) Sea Level Pressure (SLP)
# (5) Zonal winds (ZW)
# (6) Meridional winds (MW)
# (7) Sea Surface Temperature (SST)
Npred            = 7

# Set directories
Datadir          = paste(MasterDir,'/Data/',sep='') # Directory with all data
Readir           = paste(Datadir,ReaCode,'/',sep='') # Directory with reanalysis data
RefDir           = paste(Datadir,'CFSv2/1982-2009/9monthReforecasts/netcdf/',sep='') # Directory with CFSv2
Datadirsac       = '/d3/hydrofcst/overtheloop/data/output/fcst/wsf/' # Location of retro SAC simulations
start_reinit     = '1980/1/1'
end_reinit       = '2016/4/1'
## Define level index for geopotential height, zonal winds and meridional winds (case NCEP/NCAR)
# 17 Pressure levels (mb): 1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10
# index goes from 1 (1000 mb) to 17 (10 mb)
levels           = c(1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10)
lindex           = 4  # Index indicating the level pressure to be used to extract predictors

## Define working directory
if(flagflex == 0){
  Datadirpred     = paste(MasterDir,'/Results/TestBasins/Statistical_models/PLSR/SST_GPH.',Ncomp,'comp/',
                          Months2[ipredictand[1]],Months2[ipredictand[2]],sep="")
}
if(flagflex == 1){
  Datadirpred     = paste(MasterDir,'/Results/TestBasins/Statistical_models/PLSR/SST_GPH.',Ncomp,'comp/t',
                          Months2[mmax],sep="")
}

setwd(Datadirpred)

## Define number of months prior the beginning of forecast period to be considered for predictor
## selection (e.g. if it is six and the forecast is for April-July, the preceding season goes
## from October)
if(flagflex == 0){
  if( ipredictand[1] == 4 ){
    monthsprior = 8
    Leadt       = c('1001','1101','1201','0101','0201','0301','0401')
    Leadm       = c(10,11,12,1,2,3,4)  # Month for lead time forecasts
    LeadLegend  = c('Oct 1','Nov 1','Dec 1','Jan 1','Feb 1','Mar 1','Apr 1')
    initfcst    = Leadm
  }
  if( ipredictand[1] == 5 ){ # In case we have a flexible predictand, we compute predictors through May 1  
    monthsprior = 9
    Leadt       = c('1001','1101','1201','0101','0201','0301','0401','0501')
    Leadt       = c('10/1','11/1','12/1','1/1','2/1','3/1','4/1','5/1')
    Leadm       = c(10,11,12,1,2,3,4,5)  # Month for lead time forecasts
    LeadLegend  = c('Oct 1','Nov 1','Dec 1','Jan 1','Feb 1','Mar 1','Apr 1','May 1')
    initfcst    = Leadm
  }
}
if( flagflex == 1){ # In case we have a flexible predictand, we compute predictors through May 1  
  monthsprior = 9
  Leadt       = c('1001','1101','1201','0101','0201','0301','0401','0501','0601','0701')
  Leadt       = c('10/1','11/1','12/1','1/1','2/1','3/1','4/1','5/1')
  Leadm       = c(10,11,12,1,2,3,4,5,6,7)  # Month for lead time forecasts
  LeadLegend  = c('Oct 1','Nov 1','Dec 1','Jan 1','Feb 1','Mar 1','Apr 1','May 1','Jun 1','Jul 1')
  initfcst    = Leadm
}

# Define some reanalysis variables
if (ReaCode == 'NCEP_NCAR'){
  Start_reanalysis = 1948 # NCEP/NCAR reanalysis starts in Jan 1948
  End_reanalysis   = max(Years_analysis) # Ending year
  Start_sst        = 1949 # This dataset starts in January/1949
  End_sst          = max(Years_analysis) # Ending year
}
if (ReaCode == 'CFSR'){
  Start_reanalysis = 1979 # CFSR reanalysis starts in Jan 1948
  End_reanalysis   = max(Years_analysis) # Ending year
  Start_sst        = 1979 # This dataset starts in January/1979
  End_sst          = max(Years_analysis) # Ending year
}
# Define start and end of reforecast datasets
Start_reforecast = 1982 # CFSv2 starts in 1982
End_reforecast   = 2009 # Ending year
Years_ref        = Start_reforecast:End_reforecast
Nyears_ref       = length(Years_ref)
# Define start and end of actual reforecast assessment
Start_evref      = 1983 # Reforecast evaluation period starts in Jan 1983
End_evref        = 2009 # Ending year
Years_evref      = Start_evref:End_evref
Nyears_evref     = length(Years_evref)

# Define months and days
MonthsAcr   = c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
Months      = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
if( flagflex == 0){   Nseas    = length(ipredictand[1]:ipredictand[2])    }
if( flagflex == 1){   Nseas    = length(c(Leadm[1]:12,1:mmax))  }


# Define number of lead times to be used to generate hindcasts
Nlead       = ifelse(flagflex == 0, length(initfcst), length(Leadm) )
Nmonths     = 12

# Define names for reanalysis variables
Nvar          = 6  # Number of reanalysis variables
# Define long variable names for title of plots
ReanalNames   = c('Surface air temperature',
                  paste('Geopotential height (',toString(levels[lindex]),' mb)',sep=''),
                  'Precipitable water','Sea level pressure',
                  paste('Zonal winds (',toString(levels[lindex]),' mb)',sep=''),
                  paste('Meridional winds (',toString(levels[lindex]),' mb)',sep=''),
                  'Sea surface temperature')
# Define short variable names for files
ReanalShort   = c('SAT',paste('GPH',toString(levels[lindex]),'mb',sep=''),
                  'PW','SLP',
                  paste('ZW',toString(levels[lindex]),'mb',sep=''),
                  paste('MW',toString(levels[lindex]),'mb',sep=''),'SST')

# Define long variable names for title of plots
ClimNames     = c('Atlantic Multidecadal Oscillation',
                  'North Atlantic Oscillation ',
                  'East Central Tropical Pacific SST',
                  'Pacific Decadal Oscillation',
                  'Pacific North American index',
                  'Southern Oscillation Index',
                  'Extreme Eastern Tropical Pacific SST',
                  'Eastern Tropical Pacific SST',
                  'Central Tropical Pacific SST',
                  'Multivariate ENSO Index',
                  'Western Pacific Index',
                  'Tropical Northern Atlantic Index')
# Define short variable names for files
ClimShort     = c('AMO','NAO','NINO3.4','PDO','PNA','SOI',
                  'NINO 1+2','NINO3','NINO4','MEI','WP','TNA')
ClimWeb       = c('amon.us','nao','nina34','pdo','pna','soi',
                  'nina1','nina3','nina4','mei','wp','tna')
Nind          = length(ClimShort)
Start_clim    = 1951


# Define variable names for CFSv2
VarNames      = c('TMP_P8_L1_GGA0_avg1m','HGT_P8_L100_GLL0_avg1m',
                  'PWAT_P8_L200_GGA0_avg1m','PRMSL_P8_L101_GLL0_avg1m',
                  'UGRD_P8_L100_GLL0_avg1m','VGRD_P8_L100_GLL0_avg1m',
                  'POT_P8_L160_GLL0_avg1m')
FilePref      = c('tmpsfc','z700','pwat','prmsl','wnd700','wnd700','oceansst','watr')
Nadv          = 10




##########  PART 2: Load runoff data for case study basins  #################################

# Load list with basin datasets
BasinFile        = paste(Datadir,'streamflow/TestBasins.Mly.AF.rsav',sep='')  # Location of R array containing flow data
load(BasinFile)   # Load data

# Read alements in the list
BasinCode     = TestBasinsFlow$BasinCode          # 8-digit code for basins
FlowMonthly   = TestBasinsFlow$FlowMonthly        # Mean monthly runoff, organized (Jan to Dec)
YearsFlows    = TestBasinsFlow$Years              # Years containing monthly flow data
BasinAcr      = TestBasinsFlow$BasinAcr           # Years containing monthly flow data
rm(TestBasinsFlow)




##########   PART 3: Read reanalysis data  #########################################


##### Part 3.0: Download the data

if (ReaCode == 'NCEP_NCAR'){
  
  # Download NCEP/NCAR data
  Readir  = paste(Readir,'Download/',sep='') # Set directory where downloaded data will be saved
  # Download surface air temperature
  download.file('http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc',
                destfile = paste(Readir,'air_mon_mean.nc',sep=''),method='auto',quiet = FALSE,
                mode="wb", cacheOK = TRUE)
  # Download geopotential height
  download.file('http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis.derived/pressure/hgt.mon.mean.nc',
                destfile = paste(Readir,'hgt_mon_mean.nc',sep=''),method='auto',quiet = FALSE,
                mode="wb", cacheOK = TRUE)
  # Download precipitable water
  download.file('http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis.derived/surface/pr_wtr.mon.mean.nc',
                destfile = paste(Readir,'pr_wtr_mon_mean.nc',sep=''),method='auto',quiet = FALSE,
                mode="wb", cacheOK = TRUE)
  # Download sea level pressure
  download.file('http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis.derived/surface/slp.mon.mean.nc',
                destfile = paste(Readir,'slp_mon_mean.nc',sep=''),method='auto',quiet = FALSE,
                mode="wb", cacheOK = TRUE)
  # Download zonal winds
  download.file('http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.mean.nc',
                destfile = paste(Readir,'uwnd_mon_mean.nc',sep=''),method='auto',quiet = FALSE,
                mode="wb", cacheOK = TRUE)
  # Download meridional winds
  download.file('http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis.derived/pressure/vwnd.mon.mean.nc',
                destfile = paste(Readir,'vwnd_mon_mean.nc',sep=''),method='auto',quiet = FALSE,
                mode="wb", cacheOK = TRUE)
}



##### Part 3.1: Read reanalysis data (first six variables...SST is read in the next sub-section)

# Open netCDF files
file_sat   = nc_open(paste(Readir,'air_mon_mean.nc',sep='')) # Surface air temperature
file_hgt   = nc_open(paste(Readir,'hgt_mon_mean.nc',sep='')) # Geopotential height
file_prw   = nc_open(paste(Readir,'pr_wtr_mon_mean.nc',sep='')) # Precipitable water
file_slp   = nc_open(paste(Readir,'slp_mon_mean.nc',sep='')) # Sea level pressure
file_uwnd  = nc_open(paste(Readir,'uwnd_mon_mean.nc',sep='')) # Zonal winds
file_vwnd  = nc_open(paste(Readir,'vwnd_mon_mean.nc',sep='')) # Meridional winds

if (ReaCode == 'CFSR'){ # If the dataset is CFSR, data must be updated until 2014
  file_sat2   = nc_open(paste(Readir,'air_mon_mean_2.nc',sep='')) # Surface air temperature
  file_hgt2   = nc_open(paste(Readir,'hgt_mon_mean_2.nc',sep='')) # Geopotential height
  file_prw2   = nc_open(paste(Readir,'pr_wtr_mon_mean_2.nc',sep='')) # Precipitable water
  file_slp2   = nc_open(paste(Readir,'slp_mon_mean_2.nc',sep='')) # Sea level pressure
  file_uwnd2  = nc_open(paste(Readir,'uwnd_mon_mean_2.nc',sep='')) # Zonal winds
  file_vwnd2  = nc_open(paste(Readir,'vwnd_mon_mean_2.nc',sep='')) # Meridional winds
}

# Extract lat/lon coordinates and time
time       = ncvar_get(file_prw, "time")  # Time
latitude   = ncvar_get(file_prw, "lat" )  # Latitude
longitude  = ncvar_get(file_prw, "lon" )  # Longitude

if (ReaCode == 'CFSR'){ # If the dataset is CFSR, products must be remapped to a lower resolution grid
  lonorig        = longitude
  latorig        = latitude
  latorigindex   = which(latorig>=minlat & latorig<=maxlat)
  lonorigindex   = which(lonorig>=minlon & lonorig<=maxlon)
  longitude      = seq( from = min(lonorig), to = max(lonorig), by = newdx)
  latitude       = seq( from = max(latorig),  to = min(latorig),  by = -newdy)
}

# Extract lon/lat indices of interest
latindex   = which(latitude>=minlat & latitude<=maxlat)
lonindex   = which(longitude>=minlon & longitude<=maxlon)
Ny         = length(latindex)
Nx         = length(lonindex)

# Define useful dimensions
Nmonths_re = length(time)
if (ReaCode == 'CFSR'){
  addmonths  = ( End_reanalysis - 2010 )*12
  Nmonths_re = Nmonths_re + addmonths 
}
Nyears_re  = floor(Nmonths_re/Nmonths)
Years_re   = Start_reanalysis:End_reanalysis

# Extract data
if (ReaCode == 'NCEP_NCAR'){
  # Get variables from netCDF files
  sat_surf   = ncvar_get(file_sat, "air") # Surface air temperature (in C)
  hgt_surf   = ncvar_get(file_hgt, "hgt")[,,lindex,] # Geopotential height (17 levels)
  prw_surf   = ncvar_get(file_prw, "pr_wtr") # Precipitable water
  slp_surf   = ncvar_get(file_slp, "slp") # Sea level pressure (mb)
  uwnd_surf  = ncvar_get(file_uwnd, "uwnd")[,,lindex,] # Zonal winds (17 levels)
  vwnd_surf  = ncvar_get(file_vwnd, "vwnd")[,,lindex,] # Meridional winds (17 levels)
}
if (ReaCode == 'CFSR'){
  # Surface air temperature
  sat_surf   = ncvar_get(file_sat, "TMP_L104_Avg")[lonorigindex,latorigindex,] -273.15 # convert from K to C
  sat_surf2  = ncvar_get(file_sat2, "TMP_L104_Avg")[lonorigindex,latorigindex,1:addmonths] -273.15 # convert from K to C
  sat_surf   = abind(sat_surf,sat_surf2,along=3)
  rm(sat_surf2)
  sat_surf  = remap(sat_surf, x=lonorig[lonorigindex],y=latorig[latorigindex],z=1:Nmonths_re,
                    xto=longitude[lonindex],yto=latitude[latindex],zto=1:Nmonths_re)$var
  # Geopotential height
  hgt_surf   = ncvar_get(file_hgt, "HGT_L100_Avg")[lonorigindex,latorigindex,] 
  hgt_surf2  = ncvar_get(file_hgt2, "HGT_L100_Avg")[lonorigindex,latorigindex,1:addmonths] 
  hgt_surf   = abind(hgt_surf,hgt_surf2,along=3)
  rm(hgt_surf2)
  hgt_surf  = remap(hgt_surf, x=lonorig[lonorigindex],y=latorig[latorigindex],z=1:Nmonths_re,
                    xto=longitude[lonindex],yto=latitude[latindex],zto=1:Nmonths_re)$var
  # Precipitable water
  prw_surf   = ncvar_get(file_prw, "P_WAT_L200_Avg")[lonorigindex,latorigindex,]
  prw_surf2  = ncvar_get(file_prw2, "P_WAT_L200_Avg")[lonorigindex,latorigindex,1:addmonths]
  prw_surf   = abind(prw_surf,prw_surf2,along=3)
  rm(prw_surf2)
  prw_surf  = remap(prw_surf, x=lonorig[lonorigindex],y=latorig[latorigindex],z=1:Nmonths_re,
                    xto=longitude[lonindex],yto=latitude[latindex],zto=1:Nmonths_re)$var
  # Sea Level Pressure
  slp_surf   = ncvar_get(file_slp, "PRES_L101_Avg")[lonorigindex,latorigindex,]*0.01 # convert from Pa to mb
  slp_surf2  = ncvar_get(file_slp2, "PRES_L101_Avg")[lonorigindex,latorigindex,1:addmonths]*0.01 # convert from Pa to mb
  slp_surf   = abind(slp_surf,slp_surf2,along=3)
  rm(slp_surf2)
  slp_surf  = remap(slp_surf, x=lonorig[lonorigindex],y=latorig[latorigindex],z=1:Nmonths_re,
                    xto=longitude[lonindex],yto=latitude[latindex],zto=1:Nmonths_re)$var
  # U-wind speed
  uwnd_surf  = ncvar_get(file_uwnd, "U_GRD_L100_Avg")[lonorigindex,latorigindex,]
  uwnd_surf2 = ncvar_get(file_uwnd2, "U_GRD_L100_Avg")[lonorigindex,latorigindex,1:addmonths]
  uwnd_surf  = abind(uwnd_surf,uwnd_surf2,along=3)
  rm(uwnd_surf2)
  uwnd_surf  = remap(uwnd_surf, x=lonorig[lonorigindex],y=latorig[latorigindex],z=1:Nmonths_re,
                    xto=longitude[lonindex],yto=latitude[latindex],zto=1:Nmonths_re)$var
  # V-wind speed
  vwnd_surf  = ncvar_get(file_vwnd, "V_GRD_L100_Avg")[lonorigindex,latorigindex,]
  vwnd_surf2 = ncvar_get(file_vwnd2, "V_GRD_L100_Avg")[lonorigindex,latorigindex,1:addmonths]
  vwnd_surf  = abind(vwnd_surf,vwnd_surf2,along=3)
  rm(vwnd_surf2)
  vwnd_surf  = remap(vwnd_surf, x=lonorig[lonorigindex],y=latorig[latorigindex],z=1:Nmonths_re,
                     xto=longitude[lonindex],yto=latitude[latindex],zto=1:Nmonths_re)$var}

# Close files
nc_close(file_sat)
nc_close(file_hgt)
nc_close(file_prw)
nc_close(file_slp)
nc_close(file_uwnd)
nc_close(file_vwnd)
rm(file_sat,file_hgt,file_prw,file_slp,file_uwnd,file_vwnd)
if (ReaCode == 'CFSR'){
  nc_close(file_sat2)
  nc_close(file_hgt2)
  nc_close(file_prw2)
  nc_close(file_slp2)
  nc_close(file_uwnd2)
  nc_close(file_vwnd2)
  rm(file_sat2,file_hgt2,file_prw2,file_slp2,file_uwnd2,file_vwnd2)
}
rm(lonorig,latorig,latorigindex,lonorigindex)

# --- Reorganize reanalysis data as matrix with nyears x nmonths
ReanalMonthly = array(0,dim=c(Npred,Nx,Ny,Nyears_re,Nmonths))

for(ilat in 1:Ny){
  for(ilon in 1:Nx){
    ReanalMonthly[1,ilon,ilat,,]  = matrix(sat_surf[ilon,ilat,],ncol=12,byrow=T)[1:Nyears_re,]
    ReanalMonthly[2,ilon,ilat,,]  = matrix(hgt_surf[ilon,ilat,],ncol=12,byrow=T)[1:Nyears_re,]
    ReanalMonthly[3,ilon,ilat,,]  = matrix(prw_surf[ilon,ilat,],ncol=12,byrow=T)[1:Nyears_re,]
    ReanalMonthly[4,ilon,ilat,,]  = matrix(slp_surf[ilon,ilat,],ncol=12,byrow=T)[1:Nyears_re,]
    ReanalMonthly[5,ilon,ilat,,]  = matrix(uwnd_surf[ilon,ilat,],ncol=12,byrow=T)[1:Nyears_re,]
    ReanalMonthly[6,ilon,ilat,,]  = matrix(vwnd_surf[ilon,ilat,],ncol=12,byrow=T)[1:Nyears_re,]
  }
}

## Remove raw data matrices
rm(sat_surf,hgt_surf,prw_surf,slp_surf,uwnd_surf,vwnd_surf)


#### Part 3.2: Read SST 

if (ReaCode == 'NCEP_NCAR'){
  download.file('http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Diagnostic/.surface/.temp/data.nc',
                destfile = paste(Readir,'sst_mon_mean.nc',sep=''),method='auto',quiet = FALSE, mode="wb", cacheOK = TRUE)
}

file_sst  = nc_open(paste(Readir,'sst_mon_mean.nc',sep='')) # Surface air temperature
if (ReaCode == 'CFSR'){ # If the dataset is CFSR, data must be updated until 2014
  file_sst2   = nc_open(paste(Readir,'sst_mon_mean_2.nc',sep='')) # Surface air temperature
}

# Extract data
if (ReaCode == 'NCEP_NCAR'){ # from NCEP/NCAR (this is pointed as a diagnostic surface temperature)
  time_sst  = ncvar_get(file_sst, "T")  # Time
  lat_sst   = ncvar_get(file_sst, "Y")  # Latitude
  lon_sst   = ncvar_get(file_sst, "X")  # Longitude
  sst       = ncvar_get(file_sst, "temp") - 273.15 # Surface air temperature in C
}
if (ReaCode == 'CFSR'){
  time_sst       = ncvar_get(file_sst, "time") # Time
  lon_sst_orig   = ncvar_get(file_sst, "lon")  # Longitude
  lat_sst_orig   = ncvar_get(file_sst, "lat")  # Latitude
  lonorigindex   = which(lon_sst_orig>=minlon & lon_sst_orig<=maxlon)
  latorigindex   = which(lat_sst_orig>=minlat & lat_sst_orig<=maxlat)
  lon_sst        = seq( from = min(lon_sst_orig),  to = max(lon_sst_orig),  by =  newdx)
  lat_sst        = seq( from = max(lat_sst_orig),  to = min(lat_sst_orig),  by = -newdy)
  sst            = ncvar_get(file_sst, "TMP_L1_Avg")[lonorigindex,latorigindex,] - 273.15 # Surface air temperature in C
  sst2           = ncvar_get(file_sst2, "TMP_L1_Avg")[lonorigindex,latorigindex,1:addmonths] - 273.15 # Surface air temperature in C
  sst            = abind(sst,sst2,along=3)
  rm(sst2)

}

# Close files
nc_close(file_sst)
rm(file_sst)
if (ReaCode == 'CFSR'){ 
  nc_close(file_sst2)
  rm(file_sst2)
}

# Extract lon/lat indices of interest
latsst      = which(lat_sst>=minlat & lat_sst<=maxlat)
lonsst      = which(lon_sst>=minlon & lon_sst<=maxlon)
Ny_sst      = length(latsst)
Nx_sst      = length(lonsst)

if (ReaCode == 'CFSR'){ # If the dataset is CFSR, products must be remapped to a lower resolution grid
  sst = remap(sst, x=lon_sst_orig[lonorigindex], y=lat_sst_orig[latorigindex], z=1:Nmonths_re,
              xto=lon_sst[lonsst], yto=lat_sst[latsst], zto=1:Nmonths_re)$var
}


# Define useful dimensions
Nmonths_sst = length(time_sst)
if (ReaCode == 'CFSR'){
  addmonths  = ( End_reanalysis - 2010 )*12
  Nmonths_sst = Nmonths_sst + addmonths 
}
Nyears_sst  = floor(Nmonths_sst/Nmonths)
Years_sst   = Start_sst:End_sst


for(ilat in 1:Ny_sst){
  for(ilon in 1:Nx_sst){
    ReanalMonthly[7,ilon,ilat,,] = matrix(sst[ilon,ilat,],ncol=12,byrow=T)[1:Nyears_sst,]
  }
}
# Remove raw data arrays
rm(sst,lon_sst_orig,lat_sst_orig,lonorigindex,latorigindex)

filename       = paste(RefDir,FilePref[7],'.ensm.',MonthsAcr[1],'.cfsv2.data.nc',sep="")
file_id        = nc_open(filename)
auxvar         = ncvar_get(file_id, VarNames[7]) # Read variable
lon_sst_orig   = ncvar_get(file_id, "lon_0")  # Longitude
lat_sst_orig   = ncvar_get(file_id, "lat_0")  # Latitude
lonorigindex   = which(lon_sst_orig>=minlon & lon_sst_orig<=maxlon)
latorigindex   = which(lat_sst_orig>=minlat & lat_sst_orig<=maxlat)
lon_sst        = seq( from = min(lon_sst_orig),  to = max(lon_sst_orig),  by =  newdx)
lat_sst        = seq( from = max(lat_sst_orig),  to = min(lat_sst_orig),  by = -newdy)
auxvar         = remap(auxvar[lonorigindex,latorigindex,1,1], x=lon_sst_orig[lonorigindex],
                       y=lat_sst_orig[latorigindex], z=1,xto=lon_sst[lonsst], yto=lat_sst[latsst], zto=1)$var

# Obtain land maskfrom CFSv2 and correct SST field from CFSR
indexland = which(is.na(auxvar),arr.ind=TRUE)

# Set values over land as NAs
for( iyear in 1:Nyears_sst ){
  for( imonth in 1:Nmonths){
    ReanalMonthly[7,,,iyear,imonth][indexland] = NA
  }
}

## Remove garbage
rm(indexland,filename,file_id,auxvar,lon_sst_orig,
   lat_sst_orig,lonorigindex,latorigindex,lon_sst,lat_sst)


#### Part 3.3: Read climate indices

# Define variables (In this case, we also consider six climate variables!)
ClimMonthly   = array(0,dim=c(Nind,Nyears_re,Nmonths))

for(mvar in 1:Nind){
  
  if(mvar != 4){
    indexurl    = paste('https://www.esrl.noaa.gov/psd/data/correlation/',ClimWeb[mvar],
                            '.data',sep="")
    auxvector   = matrix(scan(indexurl, nmax = 2), ncol=2, byrow=T)
    auxmatrix   = matrix(scan(indexurl, skip=1, nlines=length(auxvector[1]:auxvector[2])),
                         ncol=Nmonths+1, byrow=T)
    indexclim   = which(auxmatrix[,1] %in% Years_re)
    ClimMonthly[mvar,,]  = auxmatrix[indexclim,2:13]
  }
  if(mvar == 4){
    indexurl    = url('http://research.jisao.washington.edu/pdo/PDO.latest')
    auxfile            = readLines(indexurl)
    CurrentYear        = as.numeric(format(Sys.time(),'%Y'))
    startline          = which(auxfile == 'YEAR     JAN    FEB    MAR    APR    MAY    JUN    JUL    AUG    SEP    OCT    NOV    DEC' ) + 2
    auxmatrix          = matrix(scan(indexurl, what ='list',skip=startline-1,nlines=length(1900:CurrentYear)),
                                ncol=13, byrow=T)
    close(indexurl)  # Close conection with website
    #auxmatrix[,1]      = 1900:CurrentYear
    auxmatrix[,1]      = 1900:(1900+nrow(auxmatrix)-1)
    auxmatrix          = matrix(as.numeric(auxmatrix),ncol=13,byrow=F)
    if( length(which(auxmatrix[nrow(auxmatrix),] == 1900)) >= 1 ){
      auxmatrix[nrow(auxmatrix),which(auxmatrix[nrow(auxmatrix),] == 1900):13] = -99.9
    }
    indexclim          = which(auxmatrix[,1] %in% Years_re)
    ClimMonthly[mvar,,]= auxmatrix[indexclim,2:13]
    
    if( length(which(auxmatrix[nrow(auxmatrix),] == 1900)) >= 1 ){
      auxmatrix[nrow(auxmatrix),which(auxmatrix[nrow(auxmatrix),] == 1900):13] = -99.9
    }
    
  }
  
}
## remove garbage
rm(indexurl,indexclim,auxvector,auxmatrix,auxfile)




#for(ibas in c(6,5,7,8,45)){ # Start loop over basins
###############  PART 4: Compute predictand ##########################################  
# 
# WARNING: From this section, all the analyses focuses on one particular basin
# A key part is the selection of a basin index!!!

### Part 4.1: Define dimensions and arrays to store predictors

# Define initial and final months for the preceding season (reanalysis predictors)
Rean_pred_period   = array(0,dim=c(Npred,Nlead,2))
Clim_pred_period   = array(0,dim=c(Nind, Nlead,2))


##### Part 4.2: Process seasonal values for predictand
#
# NOTE: predictands are computed differently depending on flagflex

# Define starting year for computing predictand
Start_predictand = ifelse( ipredictand[1] >= monthsprior,  min(Years_analysis), min(Years_analysis)+1 )

Nyears          = length(Start_predictand:max(Years_analysis))
Ngroups         = ceiling(Nyears/Nsyears)
cutgroup        = sapply(1:Ngroups, function(j) c(rep(j,Nsyears)))[1:Nyears]
xvalindex       = split(1:Nyears, cutgroup)
Nsegments       = length(xvalindex)-1
init_predictand = which(YearsFlows == Start_predictand) # Index in YearsFlows where analysis period start

# Define array that will contain predictors from reanalysis
Reanalysis_predictors = array(0,dim=c(Npred,Nlead,Ngroups,Nyears))

# Define array that will climate indices
Runoff                = array(NA,dim=c(Nyears,Nlead)) # Initialize vector with seasonal predictand
Ypredmonth            = array(NA,dim=c(Nyears,Nseas)) # Initialize vector with monthly predictand

# --- loop over lead times 
for(ilead in 1:Nlead){
  
  if( flagflex == 1){
    ipredictand[1] = Leadm[ilead]
    ipredictand[2] = mmax
  } # Redefine the first month of the predictand season
  
  # Compute predictand (e.g. April-July runoff)
  # NOTE THAT THE FOLLOWING LINES ONLY WILL WORK IF ipredictand[1] = 4, as defined in control_climdiag
  Runoff[,ilead] = sapply(1:Nyears, function(iyear) ifelse ( ipredictand[1]>ipredictand[2],
                        sum( c(FlowMonthly[ibas,init_predictand+iyear-2,ipredictand[1]:12],
                        c(FlowMonthly[ibas,init_predictand+iyear-1,1:ipredictand[2]]) ), na.rm=TRUE),
                        sum( FlowMonthly[ibas,init_predictand+iyear-1,
                        ipredictand[1]:ipredictand[2]], na.rm=TRUE) ) ) / (10^6)
  
}

## Clean up
#rm(FlowMonthly)



#####################  PART 5: Do PLSR over the entire field #################################  

# Start the reanalysis period in october
selpred  = c(2,7) # For now just use SSTs and GPH as in Sagarika et al. (2014,2015)
ClimInd  = array(0,dim=c(Nind,Nlead,Nyears))
AllCFSR  = array(NA,dim=c(Nlead,Nyears,9945))
Scores   = array(NA,dim=c(Ncomp,Nlead,Ngroups,Nyears))

for (init in 1:Nlead){ # Start loop over init. times

  # Use prior 3-month averages except for the first month
  # get start and end month of 3 month average
  Rean_pred_period[,init,1]  = Clim_pred_period[,init,1] = ifelse(Leadm[init] - 3 < 1, Leadm[init] + 9,  Leadm[init] - 3)
  Rean_pred_period[,init,2]  = Clim_pred_period[,init,2] = ifelse(Leadm[init] - 1 < 1, Leadm[init] + 11, Leadm[init] - 1)
  
  for(ind in 1:Nind){ # Start loop over climate indices
    
    # Define starting year for computing reanalysis predictors
    if (Clim_pred_period[ind,init,1]-monthsprior >= 0 ) {
      Start_predictor =  Start_predictand -1
    }
    if (Clim_pred_period[ind,init,1]-monthsprior < 0 & Clim_pred_period[ind,init,1]>=ipredictand[1] ) {
      Start_predictor =  Start_predictand -1
    }
    if (Clim_pred_period[ind,init,1]-monthsprior < 0 & Clim_pred_period[ind,init,1]<ipredictand[1]) {
      Start_predictor =  Start_predictand
    }
    init_predictor  = which(Years_re == Start_predictor) # Index in Years_re where analysis period start
    
    # get predictor for each year? in analysis period
    if( Clim_pred_period[ind,init,1]>Clim_pred_period[ind,init,2] ){
      ClimInd[ind,init,] = sapply(1:Nyears,
                               function(iyear) mean( c( ClimMonthly[ind,init_predictor+iyear-1,
                               Clim_pred_period[ind,init,1]:12],ClimMonthly[ind,init_predictor+iyear,
                               1:Clim_pred_period[ind,init,2]]) ) )
    }
    
    if( Clim_pred_period[ind,init,1]<=Clim_pred_period[ind,init,2] ){
      ClimInd[ind,init,] = sapply(1:Nyears,
                               function(iyear) mean( c( ClimMonthly[ind,init_predictor+iyear-1,
                               Clim_pred_period[ind,init,1]:Clim_pred_period[ind,init,2]]) ) )
    }
    
  } # End loop over selected climate indices

  
  for(ipred in selpred){ # Start loop over CFSR selected predictors
    
    # Define starting year for computing reanalysis predictors
    if (Rean_pred_period[ipred,init,1]-monthsprior >= 0 ) {
      Start_predictor =  Start_predictand -1
    }
    if (Rean_pred_period[ipred,init,1]-monthsprior < 0 & Rean_pred_period[ipred,init,1]>=ipredictand[1] ) {
      Start_predictor =  Start_predictand -1
    }
    if (Rean_pred_period[ipred,init,1]-monthsprior < 0 & Rean_pred_period[ipred,init,1]<ipredictand[1]) {
      Start_predictor =  Start_predictand
    }
    init_predictor  = which(Years_re == Start_predictor) # Index in Years_re where analysis period start
    
    # Compute seasonal values for reanalysis predictors
    if( Rean_pred_period[ipred,init,1]>Rean_pred_period[ipred,init,2] ){
      
      ArrayPred     = sapply(1:Nyears,
                             function(iyear) apply( abind( ReanalMonthly[ipred,,,
                             init_predictor+iyear-1,Rean_pred_period[ipred,init,1]:12],
                             ReanalMonthly[ipred,,,init_predictor+iyear,
                             1:Rean_pred_period[ipred,init,2]],along=3 ) , 1:2, mean) )

    }# Close condition for month_init > month_end
    
    if( Rean_pred_period[ipred,init,1]<=Rean_pred_period[ipred,init,2] ){
      
      ArrayPred     = sapply(1:Nyears, 
                             function(iyear)  apply(ReanalMonthly[ipred,,,init_predictor+iyear-1,
                             Rean_pred_period[ipred,init,1]:Rean_pred_period[ipred,init,2]],1:2,mean ) )
      
    }# Close condition for month_init <= month_end     
    
    # Step 1: Filter out data that fields containing NAs
    ArrayPred     = t(ArrayPred)
    ValidCells    = which(!is.na(ArrayPred[1,]))
    ArrayPred     = ArrayPred[,ValidCells]
    # Step 2: concatenate with the rest of predictors
    if (ipred == min(selpred)) {  Alldata   <- ArrayPred                  }
    if (ipred  > min(selpred)) {  Alldata   <- cbind(Alldata,ArrayPred)   }
    
    #Alldata   <- ArrayPred
    
  } # End loop over CFSR selected predictors
  
  ArrayCFSR        = Alldata
  ArrayCFSR        = scale(ArrayCFSR)
  AllCFSR[init,,]  = Alldata # num [1:10, 1:35, 1:9945]
  
  for(igroup in 1:Ngroups){ # Start loop over groups
    
    # Define predictand
    Robs          = Runoff[-xvalindex[[igroup]],init]        # Observed runoff for the training period
    if(islog == 1) {  Robsmod = Robs       }           # Flag # 1: streamflow in normal space
    if(islog == 2) {  # Flag # 2: streamflow in log space (Note that -Inf values must be corrected!)
      Robsmod = log(Robs)
      if( length(which(Robsmod==-Inf)==TRUE) > 0 ){
        isinf   = which(Robsmod==-Inf)
        Robsmod[isinf] = min(Robsmod[-isinf])
      }
    }
    Ypredictand   = scale(Robsmod)    # Column vector containing the predictand to be used ([(Ngroups-1)*Nsyears] x 1)
    Ymean         = mean(Robsmod)     # Mean runoff during training period
    Ysd           = sd(Robsmod)       # Standard deviation of runoff during training period
    
    # Process predictors
    # Filter points by correlation
    #indcor        = which(abs(cor(Alldata[-xvalindex[[igroup]],],Ypredictand)) >= Mincorthresh)
    #ArrayCFSR     = Alldata[,indcor]
    
    # Extract predictors (separate training from verification dataset)
    Xtrain      = ArrayCFSR[-xvalindex[[igroup]],]
    Xverif      = ArrayCFSR[ xvalindex[[igroup]],]
    if(length(xvalindex[[igroup]])==1){  Xverif  = as.vector(Xverif)   }

    # apply PLSR using function from 'PredictPLSRCrossVal.R' function
    Yaux        = PredictPLSRCrossVal(Xtrain,Xverif,Ypredictand,length(xvalindex[[igroup]]),
                                      Ncomp,Nsegments,0,Ncomp) #Nsegments?
    Yvar        = Yaux$Yvar #variance
    Xvar        = Yaux$explvar
    Scores[,init,igroup,-xvalindex[[igroup]]] = t(Yaux$Strain)
    Scores[,init,igroup, xvalindex[[igroup]]] = t(Yaux$Sverif)
    
  } # End loop over groups
  
  # Close plot with loadings
  #dev.off()
    
    
  #} # End loop over predictors
  print(paste(BasinAcr[ibas],LeadLegend[init],'is ready'))

} # End loop over init. times




# #############################  PART 6: Save data (predictand and predictors) ##############################  
# 
# # Save data for later reading!
# if ( flagflex == 0 ){ # Define output file containing time series with predictors (fixed predictand)
#   # Create list to save all predictor and predictand indices
#   AllPLSRResults = list( SeasonalRunoff = Runoff[,1], FlowUnits = 'MAF', 
#                          YearsRunoff = Start_predictand:max(Years_analysis),Fcstlead=initfcst,
#                          Nlead = Nlead, Yearsmod = Years_SAC, xvalindex = xvalindex,
#                          ReanalNames = ReanalNames, ReanalShort = ReanalShort,
#                          Scores = Scores, AllCFSR = AllCFSR, Ncomp = Ncomp,
#                          ClimNames = ClimNames, ClimShort = ClimShort, ClimInd = ClimInd,
#                          Ngroups = Ngroups, Nsyears = Nsyears)
#   flname     =  paste(BasinAcr[ibas],'.',ReaCode,'.',toString(Ncross),'xval.L',toString(Nsyears),
#                       'out.',Months2[ipredictand[1]],Months2[ipredictand[2]],'.',
#                       toString(min(Years_SAC)),'-',toString(max(Years_SAC)),'.',Ncomp,'comp.rsav',sep='')
# }
# if ( flagflex == 1 & flagres == 0){ # Define output file containing time series with predictors (flexible predictand)
#   # Create list to save all predictor and predictand indices
#   AllPLSRResults = list( SeasonalRunoff = Runoff, FlowUnits = 'MAF', Ncomp = Ncomp,
#                          YearsRunoff = Start_predictand:max(Years_analysis),Fcstlead=initfcst,
#                          Nlead = Nlead, Yearsmod = Years_SAC, xvalindex = xvalindex,
#                          ReanalNames = ReanalNames, ReanalShort = ReanalShort,
#                          Scores = Scores, AllCFSR = AllCFSR, Ncomp = Ncomp,
#                          ClimNames = ClimNames, ClimShort = ClimShort, ClimInd = ClimInd,
#                          Ngroups = Ngroups, Nsyears = Nsyears)
#   flname     =  paste(BasinAcr[ibas],'.',ReaCode,'.',toString(Ncross),'xval.L',toString(Nsyears),
#                       'out.t',Months2[ipredictand[2]],'.',
#                       toString(min(Years_SAC)),'-',toString(max(Years_SAC)),'.',Ncomp,'comp.rsav',sep='')
# }
# 
# save(AllPLSRResults, file=flname)

############################################## END OF SCRIPT ##############################################

print(paste(BasinAcr[ibas],'is ready'))

#}# End loop over basins

#}# End loop over predictand periods
