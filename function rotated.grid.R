
make_grid_rot <- function (plot, cellsize, rotation) {
  
  library(sf)
  
  #e <- as(extent(plot), "SpatialPolygons") %>% 
  #  st_as_sf()
  
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
  
  rotang = rad2deg(angle) - 0.5
  
  #https://stackoverflow.com/questions/51282724/creating-a-regular-polygon-grid-over-a-spatial-extent-rotated-by-a-given-angle
  
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  tran = function(geo, ang, center) (geo - center) * rot(ang * pi / 180) + center
  
  inpoly <- st_as_sf(plot) %>% sf::st_geometry()
  center <- st_centroid(st_union(inpoly))
  grd <- sf::st_make_grid(tran(inpoly, -rotang, center), cellsize = cellsize)
  grd_rot <- tran(grd, rotang, center) %>% st_sf('geometry' = ., data.frame('ID' = 1:length(.)))
  grid.plot <- as(grd_rot, "Spatial")
 
  # add projection information to the empty grid
  proj4string(grid.plot) <- projection(plot)
  
  return(grid.plot)
}

# Aqui estaba probando convertir el plot a raster con una resolucion determinada y hacer el grid
# ras <- raster(plot)
# res(ras) <- c(2, 2)
# ras[] <- 0
# 
# # Project the raster
# projection(ras) <- projection(plot)
# 
# st_grid <- rasterToPoints(e, spatial = TRUE)
# gridded(st_grid) <- TRUE
# st_grid <- as(st_grid, "SpatialPixels")
# b <- bbox(st_grid[1])