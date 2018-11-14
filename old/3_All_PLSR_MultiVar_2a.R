# ================================================================================ 
# Step 3: process PLSR CV results
#
# PLSR code for processing CFSv2 variables for predicting 
# -- Process CV output
# -- plot cor for all hrus
# -- plot for each HRU: Mean Loadings & Timeseries 
#
#   Created by S. Baker 
#
# ================================================================================ 
rm(list=ls())

### === Directories
source('~/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/plotLoadings.R')
source('~/s2s/analysis/scripts/cfsv2_analysis/plot_function.r')
source("/home/sabaker/s2s/analysis/scripts/cfsv2_analysis/post_process/plsr/0_control_PLSRmdl.R")
library(dplyr)
library(grid)

### INPUT HRU to look at loadings
hru_input   = c('1802', '1304', '1401', '1709', '303', '1014')
output_hru = TRUE

input_list =  list(P.2_3.Jan = list(base = c('2_3', 'prate',1), preds_ls = list(c(8,11), c(11, 4, 8), c(11,6,7,8), c(10,9,8), c(8,11,10))),
                   T.2_3.Jan = list(base = c('2_3', 'tmp2m',1), preds_ls = list(c(8,10), c(10, 2, 4), c(10,5))),
                   P.2_3.Jul = list(base = c('2_3', 'prate',7), preds_ls = list(c(8,11), c(8,11,10), c(8,6,7,5), c(8,6,7,10), c(8,6,7,11))),
                   T.2_3.Jul = list(base = c('2_3', 'tmp2m',7), preds_ls = list(c(8,10), c(8,10,1), c(8,10,2,5), c(8,10,9))),

                   P.3_4.Jan = list(base = c('3_4', 'prate',1), preds_ls = list(c(8,11), c(8,4,11),c(8,6,7,11))),
                   T.3_4.Jan = list(base = c('3_4', 'tmp2m',1), preds_ls = list(c(8,10), c(8,10,11), c(8,10,9))),
                   P.3_4.Jul = list(base = c('3_4', 'prate',7), preds_ls = list(c(8,11), c(8,11,4), c(8,11,10))),
                   T.3_4.Jul = list(base = c('3_4', 'tmp2m',7), preds_ls = list(c(8,10), c(8,10,5), c(8,10,2)))
)

## Loop through wk lead, predictor, month
for (l in 1:2) {
  
  pred_input = input_list[[l]]
  wk = pred_input$base[1]
  pred_var = pred_input$base[2]
  pred_mon = as.numeric(pred_input$base[3])
  
  pred_num = length(pred_input$preds_ls)
  
  ## loop through different predictor sets
  for (k in 1:pred_num) {
    
    preds = pred_input$preds_ls[[k]]
    
    pred_var_nm = ifelse(pred_var == 'prate', 'Precipitation', 'Temperature')
    wk_nm = ifelse(wk == '2_3', '2-3', ifelse(wk == '3_4', '3-4', '1-2'))
    lag = 1
    
    ## files to read
    file_all_multiVar      = paste0('CVresults_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',lag,'_multiVars_', paste0(preds, collapse = "."),'.rds')
    file_input_multiVar    = paste0('Input_',pred_var, '_allHRU_wk', wk, '_mon.', pred_mon,'_cut.', cut_num,'_lag.',lag,'_multiVars_',
                                    paste0(preds, collapse = "."), '.rds')
    
    ### ==== Load PLSR input data
    setwd(allHRU_multiVar_dir)
    f.out = readRDS(file = file_input_multiVar)
    
    ### ==== Load PLSR Model Data in usable format
    setwd(allHRU_multiVar_dir)
    ls_in = readRDS(file_all_multiVar)
    t = lapply(ls_in, function(x) lapply(x, identity))
    dims = c(length(t), length(t[[1]]), length(t[[1]][[1]]))
    data_all <- list()
    
    ## arrange PLSR data from parallelized output
    for (i in 2:dims[1]) { data_all[[t[[i]]$hru]] = t[[i]]  }
    for (i in 2:dims[2]) { data_all[[t[[1]][[i]]$hru]] = t[[1]][[i]]  }
    for (i in 1:dims[3]) { data_all[[t[[1]][[1]][[i]]$hru]] = t[[1]][[1]][[i]]  }
    hru_id = names(data_all)
    rm(ls_in, t)
    
    ### ==== Calculate correlations
    df_cor <- NULL
    # hru_i = hru_id[1]
    for (hru_i in hru_id) {
      nI = which(hru_i == hru_id)
      
      ## extract data and convert to numeric
      df_i = as.data.frame(data_all[[hru_i]]$df)
      # colnames(df_i) <- c('date', 'ind', 'Y_plsr', 'Y_cfs', 'Y_nld', 'res_plsr', 'res_cfs')
      
      ## calc correlations and find max cor
      cor_pred = cor(df_i$Y_plsr, df_i$Y_nld)
      cor_cfs = cor(df_i$Y_cfs, df_i$Y_nld)
      # cor_TF = cor_pred > cor_cfs
      
      # add to matrix
      df_cor = rbind.data.frame(df_cor, 
                                cbind.data.frame(hru = hru_i, cor_plsr = cor_pred, cor_cfs = cor_cfs))
    }
    
    ## format df (name, to numeric)
    df_cor$inc_cor = df_cor$cor_plsr - df_cor$cor_cfs
    df_cor$val_max <- ifelse(df_cor$inc_cor > 0, df_cor$cor_plsr, df_cor$cor_cfs)
    
    
    #### ===== Plot Correlation Increase over Raw
    ## set bases for plotting
    maxCor_in = 0.6 # round(max(na.omit(df_cor$inc_cor)), digits = 1)
    bin_size = 0.05
    
    ## set up bins
    df_cor$bins = cut(df_cor$inc_cor, breaks = seq(from = 0, to = maxCor_in, by = bin_size)) 
    df_cor$bins <- addNA(df_cor$bins)
    levels(df_cor$bins)[is.na(levels(df_cor$bins))] <- "None"
    
    ## color palette and plot
    n_col = length(levels(df_cor$bins)) -1 #number of colors needed
    pal = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
    color_fun <- colorRampPalette(pal, space = "rgb") # for discrete
    cust_pal = c(color_fun(n_col), 'grey85')  # add none color
    lab = gsub(',',' - ', levels(df_cor$bins)) # replace legend labels
    
    p1 = plot_map('prate', df_cor, 'inc_cor', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Cor. Increase') +
      ggtitle(paste('Correlation Increase from Raw CFSv2'))
    
    
    #### ===== Plot highest correlation (Raw or PLSR)
    bin_set = c(-0.2, seq(0, 0.5, by = 0.05), 0.6, 0.7, 0.8, 1)
    df_cor$bins = cut(df_cor$val_max, breaks = bin_set)
    n_col = length(levels(df_cor$bins))
    cust_pal = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
    lab = gsub(',',' - ', levels(df_cor$bins)) # replace legend labels
    
    p2 = plot_map('prate', df_cor, 'val_max', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Correlation')+
      ggtitle(paste('Max. Correlation of CFSv2 (PLSR or Raw)'))
    
    #### ===== Plot Raw CFSv2 Correlation
    df_cor$bins = cut(df_cor$cor_cfs, breaks = bin_set)
    p3 = plot_map('prate', df_cor, 'cor_cfs', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Correlation') +
      ggtitle(paste('Correlation of Raw CFSv2')) #+ theme(legend.position="none")
    
    ## save plot
    setwd(plot_dir)
    title = paste(month.abb[pred_mon], 'Week', gsub('_', '-', wk),  
                  pred_var_nm, 'Forecast - Predictors:', paste(var_ls[preds], collapse = ', '))
    g <-gridExtra::grid.arrange(p3,p2,p1,  nrow = 2, top = textGrob(title, gp=gpar(fontsize=16,font=8)))
    ggsave(paste0('plsrResults_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_', paste0(preds, collapse = "."), '.png'),
           g, height = 8, width = 11, dpi = 100)
    
    ## only increase figure
    title = paste('Predictors:', paste(var_ls[preds], collapse = ', '))
    g <-gridExtra::grid.arrange(p1 + theme(plot.title = element_blank()),  nrow = 1, top = textGrob(title, gp=gpar(fontsize=15,font=8)))
    ggsave(paste0('plsrResults_corIncrease_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_', paste0(preds, collapse = "."), '.png'), 
           g, height = 7/2.25, width = 11/2, dpi = 100)
    
    
    ### === ONE HRU: Loop through to plot loadings and timeseries ==== ###
    setwd(paste0(plot_dir, '/huc_individual_results'))
    if (output_hru) {
      for (hru_n in hru_input) {
        hru_nm = hru_tb[which(hru_tb[,1] == hru_n),2]
        
        ## set bins
        max = 35 #round(max(ld, abs(min(ld))), -1)
        dif = round(2*(max)/10)
        bin_var = seq(-max, max, dif)
        
        ### === Plot Loadings
        plot_ls = list(); n = 0
        for (i in 1:length(preds)) {
          for (comp_i in 1:Ncomp) {
            
            ## variable information
            n = n + 1
            ipred = preds[i]
            var = var_ls[ipred]
            var_long = long_nms[ipred]
            
            ## extract grid
            grid_p = (f.out$grid[ipred, , ])
            plot_grid = na.omit(data.frame(x = grid_p[1,], y = grid_p[2,]))
            
            ## loadings - loading for igroup vs. average loading of all igroups
            igroup = 1
            ind = c((f.out$ind_pred[i]+1):f.out$ind_pred[(i+1)])
            ld = na.omit(data_all[[hru_n]]$load[igroup, ind, comp_i])
            ld_avg = na.omit(apply(data_all[[hru_n]]$load[1:length(data_all[[hru_n]]$load[,1,1]), ind, comp_i], 2, mean))
            
            df_plot = cbind.data.frame(var = ld_avg, # ld_avg
                                       lat_p = plot_grid$x, lon_p = plot_grid$y)
            
            plot_ls[[n]] =
              plotLoadings(df_plot, var, #max = max, min = -max,
                           bin_vec = bin_var, title_in = paste0(var_long, ' - Comp. ', comp_i))  +
              theme(legend.position = 'none')
          }
        }
        
        ## plot & save
        title = paste('Mean Loadings -', wk_nm, 'wk Forecast of', month.abb[pred_mon], pred_var_nm, '-', hru_nm)
        g <- gridExtra::grid.arrange(grobs = plot_ls, nrow = length(preds),top = textGrob(title, gp=gpar(fontsize=15,font=8)))
        ggsave(paste0('plsrLoadings_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_huc', hru_n, '_', paste0(preds, collapse = "."), '.png'),
               g, height = 6, width = 11, dpi = 100)
        
        
        ### ==== One HRU: 1:1 plots, timeseries ==== ###
        
        ## extract data and convert to numeric
        df_i = as.data.frame(data_all[[hru_n]]$df)
        if (pred_var == 'tmp2m') { df_i[,3:5] = df_i[,3:5] - 273.15 }
        paste0(var_ls[preds], collapse = ', ')
        ymax = max(df_i$Y_nld, df_i$Y_cfs, df_i$Y_plsr)
        ymin = ifelse(pred_var == 'tmp2m', min(df_i$Y_nld, df_i$Y_cfs, df_i$Y_plsr), 0)
        
        ## plot two 1:1 plots
        g1 <- ggplot(df_i, aes(Y_nld, Y_plsr)) + geom_abline(color = 'red') + geom_point(size = 0.5) + coord_equal() +
          xlim(ymin,ymax) + ylim(ymin,ymax) + ylab(paste('PLSR', pred_var_nm)) + xlab(paste('NLDAS', pred_var_nm))
        g2 <- ggplot(df_i, aes(Y_nld, Y_cfs)) + geom_abline(color = 'red') + geom_point(size = 0.5) + coord_equal() +
          xlim(ymin,ymax) + ylim(ymin,ymax) + ylab(paste('CFSv2', pred_var_nm)) + xlab(paste('NLDAS', pred_var_nm))
        
        ## plot timeseries
        dplot = reshape2::melt(df_i[,c(1:5)], id.vars = c('date','ind'), na.rm = T)
        dplot$variable = factor(dplot$variable, c('Y_cfs', 'Y_nld', 'Y_plsr'))
        g3 = ggplot(dplot, aes(x = ind, y = value, group = variable, colour = variable)) +
          geom_line() +
          # xlab('Index of January Day (1999 - 2010)') +
          xlab('Day Index') +
          ylab(paste0(pred_var_nm)) +
          # ylab('Precipitation (mm/day per bi-weekly period)') +
          scale_color_manual(values = c("light blue", "black", 'red')) +
          scale_linetype_manual(values = c('dotted', 'solid', 'solid'))
        
        ## export and save
        g <- gridExtra::grid.arrange(grobs = list(g1, g2, g3), layout_matrix = rbind(c(1,2), c(3, 3)),
                                     top = paste(wk_nm, 'wk Forecast of', month.abb[pred_mon], pred_var_nm, '-', hru_nm, '(', paste(var_ls[preds], collapse = ', '), ')'))
        ggsave(paste0('plsrTimeSeries_wk', wk, '_mon.', pred_mon, '_',pred_var,'_multiVar_huc', hru_n, '_', paste0(preds, collapse = "."), '.png'),
               g, height = 7, width = 11, dpi = 100)
        
        
      }
    }
  }
}
