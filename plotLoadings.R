### plot loadings from plsr

# df = must be df with var, lat, lon columns
# va_nm = used to set lat and lon, must input

### === plotting things
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(maptools)
library(mapproj)

plotLoadings <- function(df, var_nm, max = NA, min = NA, bin_vec = NA, title_in = NA) {

  # Lon / Lat labels
  ewbrks = seq(0, 360, 40)
  nsbrks = seq(-90, 90, 15)
  ewlbls = unlist(lapply(ewbrks, function(x) ifelse(x > 0, paste0(x, 'ºE'), 	ifelse(x < 0, paste0(abs(x), 'ºW'), x))))
  nslbls = unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste0(abs(x), 	'ºS'), ifelse(x > 0, paste0(x, 'ºN') ,x))))
  
  ## MAP - correct projection, contries
  worldmap = map_data('world2')
  setnames(worldmap, c('X','Y', 'PID', 'POS', 'region', 'subregoin'))
  
  ## Map - option 2
  # data("wrld_simpl")
  # world = fortify(wrld_simpl)
  # Same plot, but restrict the view instead of removing points
  # allowing the complete render to happen
  # ggplot(world, mapping = aes(x = long, y = lat, group = group)) +
  #   geom_polygon(fill = NA, colour = "black") +
  #   coord_cartesian(xlim = c(-125, -30), ylim = c(-60, 35))
  
  # require(maps)
  # world_map = data.frame(map(plot=FALSE)[c("x","y")])
  # names(world_map) = c("lon","lat")
  # world_map = within(world_map, {
  #   lon = ifelse(lon < 0, lon + 360, lon)
  # })
  mp1 <- fortify(map(fill=TRUE, plot=FALSE))
  mp2 <- mp1
  mp2$long <- mp2$long + 360
  mp2$group <- mp2$group + max(mp2$group) + 1
  mp <- rbind(mp1, mp2)
  # ggplot(aes(x = long, y = lat, group = group), data = mp) + 
  #   geom_path() + 
  #   scale_x_continuous(limits = c(70, 340))
  
  # set range
  max = max(df$var)
  min = min(df$var)
  
  # # create bins of equal widths
  # scale_max = 10^(ceiling(log10(max)))/10
  # min_bin = floor(min/scale_max) *scale_max
  # max_bin = ceiling(max/scale_max) *scale_max
  # dif = round((max_bin - min_bin)/10)
  # bin_var = seq(min_bin, max_bin, dif)
  # if (max(bin_var) < max_bin) {bin_var = c(bin_var, max(bin_var)+dif)}
  if(max > 1000) {dig= 4} else {dig=3} # set sig digs
  
  # input bins
  if (!anyNA(bin_vec)) {
    bin_var = bin_vec
  }
  
  # create color palette
  cust_pal <- rev(brewer.pal(n=length(bin_var)-1, name="RdYlBu"))
  
  # add bins to df
  df$bins = cut(df$var, breaks = bin_var, dig.lab = dig)
  
  labs = levels(df$bins)
  
  ## plot with world map and discrete colors
  # ggplot() + 
  # ggplot(world, mapping = aes(x = long, y = lat, group = group)) +
  #   geom_polygon(fill = NA, colour = "black") +
  # coord_cartesian(xlim = cc(min(df$lon_p), max(df$lon_p)), ylim = c(min(df$lat_p), max(df$lat_p))) +
  ggplot() +
  # ggplot(aes(x = long, y = lat, group = group), data = mp) + 
    # geom_path(aes(x = long, y = lat, group = group), data = mp) + 
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0), limits = c(min(df$lon_p), max(df$lon_p))) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0), limits = c(min(df$lat_p), max(df$lat_p))) +
    geom_raster(data = df, aes(x = lon_p, y = lat_p, fill = bins), inherit.aes = FALSE) +

    scale_fill_manual(values = cust_pal, name = 'Loadings', drop = F, labels = labs,
                      guide = guide_legend(ncol = 1, label.position = 'right', label.hjust = 0,
                                           title.position = 'top', title.hjust = 0.5, reverse = T,
                                           keywidth = 1, keyheight = 1)) +
    theme_bw() +
    geom_path(aes(x = long, y = lat, group = group), data = mp) + 
    
    coord_equal() + # ratio of x and y...coord equal preserves the relative sizes of each
    # geom_polygon(data = worldmap, aes(X,Y,group=PID), color = 'black', size = 0.4, fill = NA) +
    xlab('') +  ylab('') +
    ggtitle(title_in) +
    theme(plot.title = element_text(size=10)) 
    # scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0), limits = c(min(df$lon_p), max(df$lon_p))) +
    # scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0), limits = c(min(df$lat_p), max(df$lat_p)))  
  # ggtitle(paste0('Correlation Coefficient of ', unlist(strsplit(df_nms[i,2], "[_]"))[1], 
  #                ' & ', var_y, ' - ', month.abb[as.numeric(df_nms[i,4])], ' ', df_nms[i,3], ' bi-weekly period ', title_nm))
  
  
}
