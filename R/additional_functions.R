##%######################################################%##
#                                                          #
####                Additional function                 ####
#                                                          #
##%######################################################%##


# Function to binarize extrapolation maps for different thresholds
th_ext <- function(x, v){
  extra_l <- list()
  for(i in 1:length(v)){
    extra_l[[i]] <- x<=v[i]
  }
  names(extra_l) <- v
  extra_l <- stack(x, stack(extra_l))
  return(extra_l)
}

# Function to calculates standard deviation of grid cell for a raster set 
sd_raster <- function(x){
  sdx <- x %>% as.data.frame() %>% na.omit %>% apply(.,1, sd)
  sdxf <- data.frame(coordinates(x), sd=values(x[[1]])) %>% na.omit
  sdxf$sd <- sdx
  gridded(sdxf) <- ~x+y
  sdxf <- stack(sdxf)
  return(sdxf)
}

# Function to transform a raster in a data.frame
rtodf <- function(x){
  x <- data.frame(coordinates(x), val=values(x))
  x <- na.omit(x)
  return(x)
}


# Fucntion to plot raster with continuous values in ggplot2
plot_cont <- function(r, title, cont=T){
  
  r2 <- rtodf(r)
  if(cont==F){
    r2[,'val'] <- as.factor(r2[,'val'])
  }
  ggplot(r2 , aes(x,y))+
    geom_raster(aes(fill=val))+
    coord_equal()+
    # layer_spatial(cntry, fill='transparent', col='gray')+
    # annotation_custom(tableGrob(perf,  theme = mytheme), xmin=-155, xmax=-140, ymin=-60, ymax=-40)+
    theme_minimal()+
    theme(legend.position = c(0.1,0.1), legend.title = element_blank(), axis.title = element_blank(),
          legend.direction = 'horizontal', axis.text = element_blank(), panel.grid = element_blank())+
    labs(subtitle =title)
}

