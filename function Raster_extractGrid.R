# FUNCION PARA GENERAR RASTERS Y MAPAS A PARTIR DE POLYGON DATA FRAME
# 
###########################################################################
require(raster)

Raster_extractGrid <- function (Grid, VARs, pixel_size, plot){
  
  raster_maps<- list()
  
  for (i in 1:length(VARs)){
    print(VARs[i])
    
    sizex <- (extent(Grid)[2]-extent(Grid)[1])/pixel_size
    sizey <- (extent(Grid)[4]-extent(Grid)[3])/pixel_size
    r <- raster(ncols=sizex, nrows=sizey, crs=crs(Grid))
    extent(r) <- extent(Grid)

    rast_tmp <- rasterize(Grid,r,field=VARs[i], fun="last", progress="text", mask=F, background=-9999) #background=-9999
    #ra <- aggregate(x, fact=2)  ## By default aggregates using mean, but see fun=
    
    
    #writeGDAL(as(rast_hab, "SpatialGridDataFrame"), paste(VARs[i],".asc",sep=""), driver = "AAIGrid", mvFlag=-9999)
    raster_maps[[i]] <- rast_tmp
    rast_tmp2 <- rast_tmp
    rast_tmp2@data@values <- ifelse (rast_tmp@data@values==-9999, NA, rast_tmp@data@values)
    
    if (plot==T){
      png(paste(VARs[i],".png",sep=""), width = 900, height = 900)
      print(spplot(rast_tmp2, main=VARs[i]))
      dev.off()
    }
  }
  
  raster_stack <- stack(raster_maps)
  names(raster_stack) <- VARs
  raster_stack
}
