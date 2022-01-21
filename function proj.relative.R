library(sp)
library(raster)
library(rgdal)
library(maptools)


rel.proj.raster <- function (raster, plot, rotation) {
  
  if (projection(plot) != projection(raster)) plot <- spTransform(plot, CRS = projection(raster))
  
  plot.corners <- plot@polygons[[1]]@Polygons[[1]]@coords
  
  x_min <- min(plot.corners[,2])[[1]]
  x_max <- max(plot.corners[,2])[[1]]
  y_min <- min(plot.corners[,1])[[1]]
  y_max <- max(plot.corners[,1])[[1]]
  x_maxy <- plot.corners[,1][plot.corners[,2] == x_max][[1]]
  y_minx <- plot.corners[,2][plot.corners[,1] == y_min][[1]]
  
  a = sqrt((plot.corners[1,1]-plot.corners[2,1])^2+(plot.corners[1,2]-plot.corners[2,2])^2)  #[Side 1 Length]
  b = sqrt((plot.corners[2,1]-plot.corners[3,1])^2+(plot.corners[2,2]-plot.corners[3,2])^2)  #[Side 2 Length]
  
  angle = atan((x_max - y_minx)/ (x_maxy - y_min)) #in radians
  
  deg2rad <- function(deg) {(deg * pi) / (180)}
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  
  if (rotation == 1) angle <- deg2rad(90) - angle 

  # Using an azimuthal projection makes it easy to rotate a raster 
  #https://gis.stackexchange.com/questions/183175/rotating-90-using-two-point-equidistant-projection-with-proj4
  
  rotate <- function(x, angle = 0, resolution = res(x)) {
    y <- x; crs(y) <- "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0"
    projectRaster(y, res=resolution, 
                  crs=paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -angle))
  }
  
  if (attributes(class(raster@file))$package == "raster") {
    raster.r <- rotate(focal(raster, w=matrix(1, 3, 3), mean), rad2deg(angle))
  }
  
  im.raster<- as.im.RasterLayer(focal(raster.r, w=matrix(1, 3, 3), mean)) #smooth a little bit data
  im.raster.r <- rotate.im(im.raster, angle = rotation * deg2rad(angle), centre="centroid")
  im.raster.r <- shift.im(im.raster.r, origin = "midpoint")
  raster.r <- as(im.raster.r, "RasterLayer")
  
  if (rotation == 1) e <- extent(-b/2, b/2, -a/2, a/2)
  if (rotation == -1) e <- extent(-a/2, a/2, -b/2, b/2)
  raster.r <- crop(raster.r, e)
  
  #ras.MDTr <- shift(ras.MDTr, dx=0, dy=0)
  
  #raster.r <- spTransform(MDS, paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -angle))
  return(raster.r)
}



proj.new <- ETRS89
if (projection(plot) != projection(proj.new)) plot <- spTransform(plot, CRS = proj.new)

plot.corners <- plot@polygons[[1]]@Polygons[[1]]@coords

x_min <- min(plot.corners[,2])[[1]]
x_max <- max(plot.corners[,2])[[1]]
y_min <- min(plot.corners[,1])[[1]]
y_max <- max(plot.corners[,1])[[1]]
x_maxy <- plot.corners[,1][plot.corners[,2] == x_max][[1]]
y_minx <- plot.corners[,2][plot.corners[,1] == y_min][[1]]

a = sqrt((plot.corners[1,1]-plot.corners[2,1])^2+(plot.corners[1,2]-plot.corners[2,2])^2)  #[Side 1 Length]
b = sqrt((plot.corners[2,1]-plot.corners[3,1])^2+(plot.corners[2,2]-plot.corners[3,2])^2)  #[Side 2 Length]

angle = atan((x_max - y_minx)/ (x_maxy - y_min)) #in radians

deg2rad <- function(deg) {(deg * pi) / (180)}
rad2deg <- function(rad) {(rad * 180) / (pi)}

if (rotation == 1) angle <- deg2rad(90) - angle

angle <- rad2deg(angle)

plot.sf <- st_as_sf (plot)
plot.sf2 <- sf::st_geometry(plot.sf)

rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
tran = function(geo, ang, center) (geo - center) * rot(ang * pi / 180) + center
center <- st_centroid(st_union(plot.sf))
grd_rot <- tran(plot.sf2, -angle, center)

grid.sf <- st_make_grid(grd_rot, cellsize = c(2,2), 
                        crs = st_crs(grd_rot)$proj4string, what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))
grid.sf2 <- sf::st_geometry(grid.sf)

kk <- as(grid.sf2, "Spatial")
plot(kk)



rotate <- function(x, angle = 0) {
  y <- x; crs(y) <- "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0"
  spTransform(y, CRS = paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -angle))
}


tran = function(geo, ang, center) (geo - center) * rot(ang * pi / 180) + center

plot.r <- rotate(plot, angle)
plot(plot.r)
makegrid()

#https://www.neonscience.org/field-data-polygons-centroids


centroids <- Points

a <- vector('list', length(2))

# loop through each centroid value and create a polygon
# this is where we match the ID to the new plot coordinates
for (i in 1:nrow(centroids)) {  # for each for in object centroids
  a[[i]]<-Polygons(list(Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), CuadCode[i]) 
  # make it an Polygon object with the Plot_ID from object ID
}
