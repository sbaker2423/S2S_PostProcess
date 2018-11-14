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

##### ===== Load Libraries ===== #####
library(ggplot2)
library(dplyr) 
library(data.table)
library(PBSmapping)
library(maps)


## Read in huc4_dt rstudio file for use
load('/home/sabaker/s2s/analysis/files/plots/HUC4_df.RData')
huc4_dt$hru = huc4_dt$HUC4
huc4_dt$hru = as.numeric(as.character(huc4_dt$hru))

## a function that returns the position of n-th largest
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
maxV <- function(n) function(x) sort(x, decreasing = TRUE)[n]


## setup inputs for looping purposes
wk_v               = c('2_3', '3_4')   # week to analyze
pred_var_v         = c('prate', 'tmp2m')        # variable to predict 'tmp2m'
pred_mon_v         = c(1, 7)                    # month to predict in CV
mon_nm_v          = c('DJF', 'JJA')            # months for model data (labeling purposes)


### === setup world map variables
## commands for making pretty spatial Lon / Lat labels
ewbrks = seq(-180, 180, 10)
nsbrks = seq(-50, 50, 10)
ewlbls = unlist(lapply(ewbrks, function(x) ifelse(x > 0, paste0(x, '°E'), 	ifelse(x < 0, paste0(abs(x), '°W'), x))))
nslbls = unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), 	'°S'), ifelse(x > 0, paste0(x, '°N') ,x))))

## cut border polygon from world map to CONUS domain
worldmap = map_data('world')
setnames(worldmap, c('X','Y', 'PID', 'POS', 'region', 'subregoin'))
worldmap = clipPolys(worldmap, xlim = c(-126, -66), ylim = c(24, 54))

w <- v <- p <- 1
## Loop through bi-weekly periods
for (wk in wk_v) {
  ## Loop through predictand
  for (pred_var in pred_var_v) {
    ## Loop through predictand month
    for (pred_mon in pred_mon_v) {
      
      
      
      # ## set values
      # wk = wk_v[w]
      # pred_var = pred_var_v[v]
      # pred_mon  = pred_mon_v[p]#; mon_nm = mon_nm_v[p]
      
      ### ==== Load Correlations for Data
      setwd(allHRU_model_dir)
      file_cor = paste0('CVcors_allVars_wk', wk, '_mon.', pred_mon, '_',pred_var,'_processed.rds')
      if (!file.exists(file_cor)){
        
        ## print error
        print('ERROR: run script 3_All_PLSR_processResults_2a.R to create file')
        
      } else {
        cor = readRDS(file_cor)
      }
      
      title = paste(month.abb[pred_mon], 'Week', gsub('_', '-', wk),  
                    ifelse(pred_var == 'prate', 'Precipitation', 'Temperature'), 'Forecast')
      
      ## find 1st, 2nd, 3rd highest correlated variables
      test = cor[c(4:14)]
      
      cor_max <- NULL
      for (i in 1:3) {
        cor_max = rbind.data.frame(cor_max, 
                                   cbind.data.frame(hru = cor$hru, 
                                                    var = colnames(test)[apply(test, 1, maxn(i))],
                                                    val = apply(test, 1, maxV(i)),
                                                    cfs   = cor$cfs))
      }
      
      
      #### ===== Plot 3 highest correlation variable from PLSR
      
      cor_max$bins = factor(cor_max$var, levels = c(var_ls))
      lab = c(var_ls, 'None')#levels(df_cormax$bins)#
      #cust_pal = rev(colorRampPalette(brewer.pal((11),"Spectral"))(length(lab)))
      cor_max$bins <- addNA(cor_max$bins)
      levels(cor_max$bins)[is.na(levels(cor_max$bins))] <- "None"
      cust_pal = c('#fb9a99','#a6cee3','#fdbf6f','#ff7f00', '#b2df8a','#cab2d6','#6a3d9a', '#1f78b4','#ffff99','#e31a1c','#33a02c')
      cust_pal = c(cust_pal, 'grey85') #add none color
      
      bin_set = c(-1, 0, 1)
      cor_max$inc = cor_max$val - cor_max$cfs
      cor_max$bins2 = cut(cor_max$inc, breaks = bin_set) 
      #n_col = length(levels(cor_max$bins2))
      #cust_pal2 = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
      
      ## change data type of hru column based on matrix
      if (!is.numeric(cor_max$hru)) {
        cor_max$hru = as.numeric(as.character(cor_max$hru))
      }
      filter(cor_max, hru == 1401)
      
      p_ls <- list()
      for (var_i in var_ls) {
        
        df_i = filter(cor_max, var == var_i)
        
        # p_ls[[var_i]] = plot_map('prate', df_i, 'var', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Variables') + 
        #   theme(legend.position="none",
        #         plot.title = element_text(size = 10)) +
        #   ggtitle(paste(long_nms[[which(var_i == var_ls)]]))
        
        ## join matricies by HUC4 id 
        var_huc = left_join(huc4_dt, df_i, by = "hru")
        
        #cust_pal = color_test(nlevels(var_huc$bins))
        p_ls[[var_i]] = ggplot() + 
          # add boarders
          geom_polygon(data = worldmap, aes(X,Y,group=PID), color = NA, size = 0.5, fill = 'grey85') + 
          # HUC4 polygon setup
          geom_polygon(data = var_huc, aes(x = long, y = lat, group = hru, fill = bins, alpha = bins2)) + 
          # fill color scale, legend appearance and title
          scale_fill_manual(values = cust_pal, drop = F, labels = lab,
                            guide = guide_legend(ncol = 1, label.position = 'right', label.hjust = 0,
                                                 title.position = 'top', title.hjust = 0.5, reverse = T,
                                                 keywidth = 1, keyheight = 1)) + 
          scale_alpha_manual(values = c(0.30, 1)) +
          # Draw HUC polygons
          geom_path(data = var_huc, aes(x = long, y = lat, group = hru), color = 'grey73', size = 0.15) + #grey60
          # ggplot has built-in themes
          theme_bw() + 
          # legend location...theme can take in other appearance / layout arguments - NOTE, if used with a built-in theme, e.g. theme_bw() it needs to be called after to not get overridden
          theme(legend.position = 'right') +
          # ratio of x and y...coord equal preserves the relative sizes of each
          coord_equal() +
          # label on the x axis
          xlab('') +  
          # label on the y axis
          ylab('') + 
          # state borders
          geom_polygon(data=map_data("state"), aes(x=long, y=lat, group=group), color="grey42", size = 0.3, fill=NA, linetype = 'dashed')+
          # tic mark locations and labels on the x axis
          scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
          # tic mark locations and labels on the y axis 
          scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) + 
          # set margins (wider on right)
          theme(plot.margin=unit(c(0.05,0.01,0.05,0.01),'cm')) + #top, right, bottom, left 
          theme(legend.position="none",
                plot.title = element_text(size = 10)) +
          ggtitle(paste(long_nms[[which(var_i == var_ls)]]))
        
        
      }
      
      
      g <-gridExtra::grid.arrange(p_ls[[var_ls[1]]], p_ls[[2]], p_ls[[3]], p_ls[[4]], p_ls[[5]], 
                                  p_ls[[6]], p_ls[[7]], p_ls[[8]], p_ls[[9]], p_ls[[10]], p_ls[[11]], nrow = 4, top = title)
      
      
      ## save plot
      setwd(plot_dir)
      ggsave(paste0('CorTop3Vars_wk', wk, '_mon.', pred_mon, '_',pred_var,'.png'),
             g, height = 7, width = 11, dpi = 150)
      
      
      
    }
  }
}


# bin_set = c(-0.2, seq(0, 0.5, by = 0.05), 0.6, 0.7, 0.8, 1)
# df_i$inc = df_i$val - df_i$cfs
# df_i$bins2 = cut(df_i$inc, breaks = bin_set) 
# n_col = length(levels(df_i$bins2))
# cust_pal2 = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
# 
# ##### ===== Load Libraries ===== #####
# library(ggplot2)
# library(dplyr) 
# library(data.table)
# library(PBSmapping)
# library(maps)
# 
# ##### ===== Read HUC4 data and combine with forecast ===== #####
# 
# ## Read in huc4_dt rstudio file for use
# load('/home/sabaker/s2s/analysis/files/plots/HUC4_df.RData')
# huc4_dt$hru = huc4_dt$HUC4
# huc4_dt$hru = as.numeric(as.character(huc4_dt$hru))
# 
# 
# ## change data type of hru column based on matrix
# if (!is.numeric(df_i$hru)) {
#   df_i$hru = as.numeric(as.character(df_i$hru))
# }
# 
# ## join matricies by HUC4 id 
# var_huc = left_join(huc4_dt, df_i, by = "hru")
# 
# 
# ##### ====== Plotting ====== #####
# 
# ## commands for making pretty spatial Lon / Lat labels
# ewbrks = seq(-180, 180, 10)
# nsbrks = seq(-50, 50, 10)
# ewlbls = unlist(lapply(ewbrks, function(x) ifelse(x > 0, paste0(x, '°E'), 	ifelse(x < 0, paste0(abs(x), '°W'), x))))
# nslbls = unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), 	'°S'), ifelse(x > 0, paste0(x, '°N') ,x))))
# 
# 
# 
# ## cut border polygon from world map to CONUS domain
# worldmap = map_data('world')
# setnames(worldmap, c('X','Y', 'PID', 'POS', 'region', 'subregoin'))
# worldmap = clipPolys(worldmap, xlim = c(-126, -66), ylim = c(24, 54))
# 
# #cust_pal = color_test(nlevels(var_huc$bins))
# p = ggplot() + 
#   # add boarders
#   geom_polygon(data = worldmap, aes(X,Y,group=PID), color = NA, size = 0.5, fill = 'grey85') + 
#   # HUC4 polygon setup
#   geom_polygon(data = var_huc, aes(x = long, y = lat, group = hru, fill = bins, alpha = bins2)) + 
#   # fill color scale, legend appearance and title
#   scale_fill_manual(values = cust_pal, drop = F, labels = lab,
#                     guide = guide_legend(ncol = 1, label.position = 'right', label.hjust = 0,
#                                          title.position = 'top', title.hjust = 0.5, reverse = T,
#                                          keywidth = 1, keyheight = 1)) + 
#   scale_alpha_manual(values = c(0.1, seq(0.35, by = 0.05, length.out = 14))) +
#   # Draw HUC polygons
#   geom_path(data = var_huc, aes(x = long, y = lat, group = hru), color = 'grey73', size = 0.15) + #grey60
#   # ggplot has built-in themes
#   theme_bw() + 
#   # legend location...theme can take in other appearance / layout arguments - NOTE, if used with a built-in theme, e.g. theme_bw() it needs to be called after to not get overridden
#   theme(legend.position = 'right') +
#   # ratio of x and y...coord equal preserves the relative sizes of each
#   coord_equal() +
#   # label on the x axis
#   xlab('') +  
#   # label on the y axis
#   ylab('') + 
#   # state borders
#   geom_polygon(data=map_data("state"), aes(x=long, y=lat, group=group), color="grey42", size = 0.3, fill=NA, linetype = 'dashed')+
#   # tic mark locations and labels on the x axis
#   scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
#   # tic mark locations and labels on the y axis 
#   scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) + 
#   # set margins (wider on right)
#   theme(plot.margin=unit(c(0.15,0.5,0.15,0.05),'cm')) #top, right, bottom, left 
# 
# 
# 
# 


# library(grid)
# library(ggplot2)
# library(gridExtra)
# 
# grid_arrange_shared_legend <- function(...) {
#   plots <- list(...)
#   g <- ggplotGrob(plots[[1]] + theme(legend.position="right"))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   grid.arrange(
#     do.call(arrangeGrob, lapply(plots, function(x)
#       x + theme(legend.position="none"))),
#     legend,
#     ncol = 1,
#     heights = unit.c(unit(1, "npc") - lheight, lheight))
# }
# 
# 
# grid_arrange_shared_legend(p_ls[[var_ls[1]]], p_ls[[2]], p_ls[[3]], p_ls[[4]], p_ls[[5]], 
#                            p_ls[[6]], p_ls[[7]], p_ls[[8]], p_ls[[9]], p_ls[[10]], p_ls[[11]], nrow = 6)
# 
# p4 = plot_map('prate', df_i, 'var', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Variables')+ 
#   ggtitle(paste('Highest Correlated CFSv2 Variables'))
# 


# ## set bases for plotting
# #maxCor_in = round(max(na.omit(cor$inc_cor)), digits = 1)
# bin_size = 0.05
# maxCor_in = 0.6
# 
# 
# 
# #### ===== Plot Increase in correlation from a single variable
# ## arrange df and change factors of NA
# df_plot = cor[c('hru','inc_cor', 'val_max', 'var_max', 'cfs')]
# df_plot$inc_cor[is.na(df_plot$inc_cor)] <- 0
# df_plot$val_max[is.na(df_plot$val_max)] <- df_plot$cfs[is.na(df_plot$val_max)]
# df_plot$bins = cut(df_plot$inc_cor, breaks = seq(from = 0, to = maxCor_in, by = bin_size)) 
# df_plot$bins <- addNA(df_plot$bins)
# levels(df_plot$bins)[is.na(levels(df_plot$bins))] <- "None"
# 
# ## color palette and plot
# n_col = length(levels(df_plot$bins)) -1 #number of colors needed
# pal = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
# color_test <- colorRampPalette(pal, space = "rgb") #for discrete
# cust_pal = color_test(n_col)
# cust_pal = c(cust_pal, 'grey85') #add none color
# lab = gsub(',',' - ', levels(df_plot$bins)) #replace legend labels
# 
# p1 = plot_map('prate', df_plot, 'inc_cor', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Cor. Increase') +
#   ggtitle(paste('Correlation Increase from Raw CFSv2'))
# 
# 
# #### ===== Plot Raw and PLSR Correlation
# # max_cor = 0.75#round(max(na.omit(df_plot$val_max)), 1) #0.65 #
# # min_cor = -0.1 #round(min(na.omit(df_plot$val_max)), digits = 1) #0 #
# bin_set = c(-0.2, seq(0, 0.5, by = 0.05), 0.6, 0.7, 0.8, 1)
# df_plot$bins = cut(df_plot$val_max, breaks = bin_set) 
# n_col = length(levels(df_plot$bins))
# cust_pal = rev(colorRampPalette(brewer.pal((11),"Spectral"))(n_col))
# lab = gsub(',',' - ', levels(df_plot$bins)) #replace legend labels
# 
# # color_test <- colorRampPalette(c('#2166ac','#f7f7f7', '#b2182b'),space = "rgb") #for discrete
# # cust_pal = color_test(n_col)
# p2 = plot_map('prate', df_plot, 'val_max', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Correlation')+
#   ggtitle(paste('Max. Correlation of CFSv2 (PLSR or Raw)'))
# 
# #### ===== Plot Raw CFSv2 Correlation
# df_plot$bins = cut(df_plot$cfs, breaks = bin_set) 
# p3 = plot_map('prate', df_plot, 'cfs', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Correlation') +
#   ggtitle(paste('Correlation of Raw CFSv2'))
# 
# #### ===== Plot highest correlation variable from PLSR
# df_plot$bins = factor(var_ls[df_plot$var_max], levels = c(var_ls))
# lab = c(var_ls, 'None')#levels(df_plot$bins)#
# #cust_pal = rev(colorRampPalette(brewer.pal((11),"Spectral"))(length(lab)))
# df_plot$bins <- addNA(df_plot$bins)
# levels(df_plot$bins)[is.na(levels(df_plot$bins))] <- "None"
# cust_pal = c('#fb9a99','#a6cee3','#fdbf6f','#ff7f00', '#b2df8a','#cab2d6','#6a3d9a', '#1f78b4','#ffff99','#e31a1c','#33a02c')
# cust_pal = c(cust_pal, 'grey85') #add none color
# 
# p4 = plot_map('prate', df_plot, 'var_max', type = 'discrete', lab = lab, cust_pal = cust_pal, legend_title = 'Variables')+ 
#   ggtitle(paste('Highest Correlated CFSv2 Variables'))


# ## save plot
# setwd(plot_dir)
# g <-gridExtra::grid.arrange(p3,p1,p2, p4, nrow = 2, top = title)
# ggsave(paste0('CVresults_allVars_wk', wk, '_mon.', pred_mon, '_',pred_var,'.png'), 
#        g, height = 7, width = 11, dpi = 100)
