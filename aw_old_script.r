#  pKNN -- predictive analogs
#  by A. Wood, 2011, adapting original KNN code from S. Gangopadhyay
#  apply KNN but via PLSR rather than PCA
#  instead of operating on a single feature vector, this version 
#  generates feature vectors from the climatology operates on each, and 
#    also generates analogs for them

# cross-validation version

library(pls)

loc = "dalle"  # loc for prediction (also dalle grvu1 gjnc2 bffu1 glda3)
Rstr = "20"    # r string for input file name

# --- various sets of inputs ---

# full predictors, single predictand
#X=read.table("./clim_ndx/clim_ndx.I_dI.60-10.Nov1.tbl", header=T) # Nov 1
#X=read.table("./clim_ndx/clim_ndx.I_dI.60-10.Dec1.tbl", header=T)  # Dec 1
#q_fl   <- sprintf("./obs_q/%s.AJvol.61-11", loc)   
#outfl  <- sprintf("./out_plsr/plsr.pknn.%s.txt3", loc)     # nov 1
#alogfl <- sprintf("./out_plsr/analogs.pknn2.%s.txt3", loc)
#wgtfl  <- sprintf("./out_plsr/weights.pknn2.%s.txt3", loc)

# DEC 1, reduced predictors, single predictand
in_fl  <- sprintf("./clim_ndx/%s/clim_ndx.60-10.Dec1.r_gt%s.%s", loc,Rstr,loc) 
q_fl   <- sprintf("./obs_q/%s.AJvol.61-11", loc)   # out3
q2_fl  <- sprintf("./obs_q/%s.Nvol.60-10", loc)    # nov flows, predictor
bothfl <- sprintf("./clim_ndx/%s/obs_clim_ndx.Dec1.r_gt%s.%s", loc,Rstr,loc) 
outfl  <- sprintf("./out_plsr/plsr.pls_knn.dec1.%s", loc)     
alogfl <- sprintf("./out_plsr/analogs.pls_knn.dec1.%s", loc)
wgtfl  <- sprintf("./out_plsr/weights.pls_knn.dec1.%s", loc)

# read predictor and predictand tables, make formula
X  = read.table(in_fl, header=T) 
Q  = read.table(q2_fl, header=T)   # fall flows
X2 = data.frame(Q,X)              # add to predictors
Y  = read.table(q_fl,header=T)
Z  = data.frame(Y,X2)   # all data, maybe not used
form <- formula(paste(colnames(Y)," ~ ",paste(colnames(X2),collapse=" + ")))

# ---- settings -----
K     = 17   # perhaps, roughly 1/3 of possible choices (+1 for self-match)
# K must be less than records in training set (all-ntest)
H     = 0    # number of X records to assess (0 for all, less for debugging)
SYR   = 1961 # start year of Y, useful for adjusting a'log indexes

#---- Input Notes ----
#  K - number of nearest neighbor analogs to return 
#  X - Data matrix of historical values used for similarity calculation
#  Y - Data vector of historical predictand values for identifying feature
#      matrices in analog selection
#      assume no missing values
# ---------------------------------

# --- NOT USED -- data analysis on all data first? ---
# prefiltering?
#sX = scale(X2)                 # standardize (mean=0; var=1) 
#sY = scale(Y)                  # could save mu/sigma, see xtra code at end
corr_all = cor(X2,Y)
corrwgt_all = (corr_all^2)/sum(corr_all^2)


# initialize some settings  --------
nrec = dim(X2)[1]                     # number of time elements (records) 
nvar = dim(X2)[2]                     # number of EOFs
dimY = dim(Y)[2]                      # length of predictand feature (# pred)
D = vector()                          # initialize distance vector
wgt = vector()                        # init. weight vector
ayear = vector()                      # initialize alog year vector
if(H==0){ H = nrec }                  # reset # recs to find analogs for
all_ndx = matrix(1:(H*K),ncol=K)      # initialize summary of all indices
all_wgt = matrix(1:(H*K),ncol=K)      # initialize summary of all weights
YR <- seq(length=nt, from=SYR, by=1)  # initialize year list

# open output files
fp1 <- file(alogfl, "w") 
fp2 <- file(wgtfl, "w") 

# diagnostic on all data here?  

# =============================================================
#   LOO evaluation:  analyze the historical data matrix leaving out 
#   one record at a time
# =============================================================
nset=0
avgRsq=0
while (nset < nrec) {  # do one less set than number of records
  nset = nset+1
  print(sprintf("doing PLSR with LOO=rec %d of %d",nset,H))
  
  Ztrn  = Z[-nset,]     # remove test data from training matrices
  YRtrn = YR[-nset]    # ... and year list for training set
  
  sZtrn   = scale(Ztrn)                    # standardize (mean=0; var=1) 
  mu_z    = attr(sZtrn,"scaled:center")    # store mean values of X[1:nvar]
  sigma_z = attr(sZtrn,"scaled:scale")     #  ... and std. devs
  
  # ----- calculate PLSR components in training set -----------
  
  zz = plsr(form, data=data.frame(sZtrn))
  P  = scores(zz)             # CCA X component canonical variables
  # ie PC time series
  #  cols=cvs (time dim by min(cols(X,Y)
  U  = loadings(zz)	      # Eigen Vectors (ie EOFs)
  
  # use Y variance explained instead of eigenvalues for weights
  #  W = zz$Xvar   # eigenval. (X variance)
  W <- drop(R2(zz, estimate = "train", intercept = FALSE)$val)
  # var expl in predictand Y by each X canon. variable
  trW = W[nvar]				#trace of [diag(W)]
  Wcum <- W        # cumulative expl var
  for(c in nvar:2) {
    W[c] = W[c]-W[c-1]  # expl var by component
  }
  
  # ---- picking number of features to use ---
  
  # -- use crossvalidation PRESS (minimum)
  zzcv <- crossval(zz)
  nret <- match(min(zzcv$validation$PRESS), zzcv$validation$PRESS)
  # -- save some stats
  tmpcor = cor(sZtrn[,1],zzcv$validation$pred[,,nret])^2
  avgRsq = avgRsq + tmpcor  # get avg Xval correlation
  print(sprintf("  Rsq=%.3f with nret components = %d",tmpcor,nret))
  
  # nret = nvar, 1           # or use all or 1
  #  -- or set cutoff on var expl. in Y
  #  cutoff = 0.90	# retain PCs to explain ~90% of the total variance expl
  #  p = Wcum/Wcum[nvar]  # cumulative fractional of tot var. expl. in Y by PCs
  #  nret = min(length(p[p<=cutoff])+1,nvar)  
  
  
  # ====== now test LOO record ===========
  # set record as F - "feature vector" for which similarity is sought in sXtrn.
  # reproject test record to the eigenvector space and calculate similarity
  
  Ztest = as.numeric((Z[nset,]-mu_z)/sigma_z)  # standardize testrec
  F = Ztest[2:(nvar+1)]           # carve out X part
  
  pF = F%*%U    # projected test record onto training variates
  # gives vector w/ len of predictand ncomps, a score for each
  
  # for each component c, pF is same as sum(F*U[,c])
  #   can weight this projection by correlations of the component with 
  #   the target, i.e., 
  corwgt <- (cor(sXtrn,sYtrn)^2)/sum(cor(sXtrn,sYtrn)^2)
  
  # distance calculation loop through
  for (i in 1:(nrec-1)){
    
    # loop through number of PCs retained (nret)
    d = 0.0
    for (iret in 1:nret){
      diff = (P[i,iret]-pF[iret])        # difference between PCs and pF
      weight = W[iret]/sum(W[1:nret])    # bi-square weighting
      # for others, see code at end
      d = d + weight*(diff*diff)
    }
    D[i] = sqrt(d)  # distance measure
  }
  
  ndx = order(D)[1:K]  # sort distances and return order, low to high
  # indexes are for training array
  
  # now can transform distance measure into weighting --
  # do inv. dist wgt. over K analogs that will be returned, normalized to 1
  sumw = 0
  for (j in 1:K){         
    wgt[j+1] = 1/D[ndx[j]]
    sumw = sumw + wgt[j+1]
  }
  wgt = wgt/sumw   # normalize
  wgt[1] = YR[nset]  # list year as first col for output processing
  
  ayear[1] = YR[nset]  # list year as first col...
  for (j in 1:K){      # adjust index by year (works for year a'logs)
    ayear[j+1] = YRtrn[ndx[j]]
    all_ndx[nset,j] = YRtrn[ndx[j]]-SYR+1  # store indices of a'log yrs
    all_wgt[nset,j] = wgt[j+1]
  }
  
  # append results for on test record to file
  write(ayear, file=fp1, ncolumns=K+1, append = TRUE)
  write(wgt, file=fp2, ncolumns=K+1, append = TRUE)
  
} # end while loop
print(sprintf("Mean Rsq=%.3f",avgRsq/nrec))

close(fp1)
close(fp2)

# ========================================
# ============== analysis ================
# ========================================
# make weighted analog ensemble means for different ensemble sizes

print("checking analog ensemble predictability")
enscor = vector()
for (NSET in 1:K){
  fY <- Y*0  # initialize object matching dimensions of Y
  for (y in 1:H){
    for (m in 1:dimY){
      tmpwgt = 0
      for (j in 1:NSET){      # top NSET
        fY[y,m] = fY[y,m] + Y[(all_ndx[y,j]),m] * all_wgt[y,j]
        tmpwgt = tmpwgt + all_wgt[y,j]
      }
      fY[y,m] = fY[y,m]/tmpwgt
    }
  }
  enscor[NSET] = cor(Y[,1],fY[,1])
}
enscor

# also make vectors based on all data to output
# could also store vector weights in Xvalidation for first feature

# now write out zz object with various components (must be a better way)
# NOTE: this will show the results last training set used only -- 

#write("explained var in Y",outfl)
#write.table(W,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("explained var in X",outfl)
#write.table(zz$Xvar,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("xcoef",outfl,append=TRUE)
#write.table(zz$xcoef,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("ycoef",outfl,append=TRUE)
#write.table(zz$ycoef,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("xscores",outfl,append=TRUE)
#write.table(zz$scores$xscores,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("yscores",outfl,append=TRUE)
#write.table(zz$scores$yscores,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("corr.X.xscores",outfl,append=TRUE)
#write.table(zz$scores$corr.X.xscores,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("corr.Y.xscores",outfl,append=TRUE)
#write.table(zz$scores$corr.Y.xscores,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("corr.X.yscores",outfl,append=TRUE)
#write.table(zz$scores$corr.X.yscores,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("corr.Y.yscores",outfl,append=TRUE)
#write.table(zz$scores$corr.Y.yscores,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)


# ============== more writes ==============
skill <- cor(Y,fY)^2
write(sprintf("r-sq skill of %d analog ens mean",10),outfl,append=TRUE)
write.table(skill,outfl,append=TRUE)
write(" ",outfl,append=TRUE)

write("forecasts (with total)",outfl,append=TRUE)
write.table(fY,outfl,append=TRUE)
write(" ",outfl,append=TRUE)




# -------- other stuff, not used -----------


#number of principal components retained (in case of many)
#retain PCs to explain approximately 90% of the total variance
#cutoff=1.0			      # just use all...
#p=cumsum(W/trW)                       # % cumul. variance explained by PCs
#nret=length(p[p<=cutoff])             #number of PCs retained


# other outputs

#write("predictor means in standardization",outfl,append=TRUE)
#write.table(mu_x,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("predictor sigmas in standardization",outfl,append=TRUE)
#write.table(sigma_x,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("predictand means in standardization",outfl,append=TRUE)
#write.table(mu_y,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)

#write("predictand sigmas in standardization",outfl,append=TRUE)
#write.table(sigma_y,outfl,append=TRUE)
#write(" ",outfl,append=TRUE)


# other weighting strategies for distance metric
#d=d+(W[iret]/trW)*(diff*diff)   # other
#d=d+(diff*diff)		 # Equal square weighting..