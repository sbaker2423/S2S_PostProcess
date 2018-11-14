# ================================================================================ 
# Step 3: process PLSR CV results
#
# PLSR code for processing CFSv2 variables for predicting 
# -- Process CV output
#   Created by S. Baker 
#
# ================================================================================ 
rm(list=ls())

## --- Directories
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/ploLoadings.R')
source('~/s2s/analysis/scripts/cfsv2_analysis/plot_function.r')
source("/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/0_control_PLSRmdl.R")

library(dplyr)


### ==== Load & Calculate Correlations for Data
setwd(allHRU_model_dir)
file_cor = paste0('CVcors_allVars_wk', wk, '_mon.', pred_mon, '_',pred_var,'_processed.rds')
if (!file.exists(file_cor)){
  
  ### ==== Arrange Data in usable format
  ls_in = readRDS(file_all)
  t = lapply(ls_in, function(x) lapply(x, identity))
  dims = c(length(t), length(t[[1]]), length(t[[1]][[1]]))
  data_all <- list()
  
  for (i in 2:dims[1]) { data_all[[t[[i]]$hru]] = t[[i]]  }
  for (i in 2:dims[2]) { data_all[[t[[1]][[i]]$hru]] = t[[1]][[i]]  }
  for (i in 1:dims[3]) { data_all[[t[[1]][[1]][[i]]$hru]] = t[[1]][[1]][[i]]  }
  hru_id = names(data_all)
  rm(ls_in, t)
  
  ### ==== Calculate correlations
  hru_id = names(data_all)
  cor = matrix(NA, nrow = length(hru_id), ncol = 15)
  cols = c(1, 3:8)  
  for (hru_i in hru_id) {
    nI = which(hru_i == hru_id)
    
    ## extract data and convert to numeric
    df_i = as.data.frame(data_all[[hru_i]]$df)
    colnames(df_i) <- c('pred', 'date', 'ind', 'Y_plsr', 'Y_cfs', 'Y_nld', 'res_plsr', 'res_cfs')
    df_i[,cols] = apply(df_i[,cols], 2, function(x) as.numeric(as.character(x))) # to numeric
    
    ## calc correlations and find max cor
    df_grp = group_by(df_i, pred)
    cor_pred = df_grp %>% summarize(cor(Y_plsr, Y_nld))
    cor_cfs = df_grp %>% summarize(cor(Y_cfs, Y_nld))
    cor_TF = any(cor_pred[,2] > cor_cfs[,2])
    cor_var = ifelse(cor_TF, which.max(as.numeric(unlist(cor_pred[,2]))), NA) # which var max cor
    cor_val = ifelse(cor_TF, max(as.numeric(unlist(cor_pred[,2]))), NA) # which var max cor
    
    # add to matrix
    cor[nI,] = as.matrix(cbind(hru_i, cor_var, cor_val, t(cor_pred[,2]), cor_cfs[1,2] ))
    
  }
  
  ## format df (name, to numeric)
  colnames(cor) <- c('hru', 'var_max', 'val_max', var_ls, 'cfs')
  cor = data.frame(cor)
  cols = 2:ncol(cor)
  cor[,cols] = apply(cor[,cols], 2, function(x) as.numeric(as.character(x))) # to numeric
  cor$inc_cor = cor$val_max - cor$cfs
  
  ## save data - very slow - reading is longest part
  saveRDS(cor, file = file_cor)

} else {
  cor = readRDS(file_cor)
}


## set bases for plotting
#maxCor_in = round(max(na.omit(cor$inc_cor)), digits = 1)
bin_size = 0.05
maxCor_in = 0.6



#### ===== Plot Increase in correlation from a single variable
## arrange df and change factors of NA
df_plot = cor[c('hru','inc_cor', 'val_max', 'var_max', 'cfs')]
df_plot$inc_cor[is.na(df_plot$inc_cor)] <- 0
df_plot$val_max[is.na(df_plot$val_max)] <- df_plot$cfs[is.na(df_plot$val_max)]
df_plot$bins = cut(df_plot$inc_cor, breaks = seq(from = 0, to = maxCor_in, by = bin_size)) 
df_plot$bins <- addNA(df_plot$bins)
levels(df_plot$bins)[is.na(levels(df_plot$bins))] <- "None"

## color palette and plot
n_col = length(levels(df_plot$bins)) -1 #number of colors needed
pal = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
color_test <- colorRampPalette(pal, space = "rgb") #for discrete
cust_pal = color_test(n_col)
cust_pal = c(cust_pal, 'grey85') #add none color
lab = gsub(',',' - ', levels(df_plot$bins)) #replace legend labels

p1 = plot_map('prate', df_plot, 'inc_cor', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Cor. Increase') +
  ggtitle(paste('Correlation Increase from Raw CFSv2'))


#### ===== Plot Raw and PLSR Correlation
# max_cor = 0.75#round(max(na.omit(df_plot$val_max)), 1) #0.65 #
# min_cor = -0.1 #round(min(na.omit(df_plot$val_max)), digits = 1) #0 #
bin_set = c(-0.2, seq(0, 0.5, by = 0.05), 0.6, 0.7, 0.8, 1)
df_plot$bins = cut(df_plot$val_max, breaks = bin_set) 
n_col = length(levels(df_plot$bins))
cust_pal = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
lab = gsub(',',' - ', levels(df_plot$bins)) #replace legend labels

# color_test <- colorRampPalette(c('#2166ac','#f7f7f7', '#b2182b'),space = "rgb") #for discrete
# cust_pal = color_test(n_col)
p2 = plot_map('prate', df_plot, 'val_max', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Correlation')+
  ggtitle(paste('Max. Correlation of CFSv2 (PLSR or Raw)'))

#### ===== Plot Raw CFSv2 Correlation
df_plot$bins = cut(df_plot$cfs, breaks = bin_set) 
p3 = plot_map('prate', df_plot, 'cfs', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Correlation') +
  ggtitle(paste('Correlation of Raw CFSv2'))

#### ===== Plot highest correlation variable from PLSR
df_plot$bins = factor(var_ls[df_plot$var_max], levels = c(var_ls))
lab = c(var_ls, 'None')#levels(df_plot$bins)#
#cust_pal = rev(colorRampPalette(brewer.pal((11),"Spectral"))(length(lab)))
df_plot$bins <- addNA(df_plot$bins)
levels(df_plot$bins)[is.na(levels(df_plot$bins))] <- "None"
cust_pal = c('#fb9a99','#a6cee3','#fdbf6f','#ff7f00', '#b2df8a','#cab2d6','#6a3d9a', '#1f78b4','#ffff99','#e31a1c','#33a02c')
cust_pal = c(cust_pal, 'grey85') #add none color

p4 = plot_map('prate', df_plot, 'var_max', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Variables')+ 
  ggtitle(paste('Highest Correlated CFSv2 Variables'))


## save plot
setwd(plot_dir)
g <-gridExtra::grid.arrange(p3,p1,p2, p4, nrow = 2, top = title)
ggsave(paste0('CVresults_allVars_wk', wk, '_mon.', pred_mon, '_',pred_var,'.png'), 
       g, height = 7, width = 11, dpi = 100)
