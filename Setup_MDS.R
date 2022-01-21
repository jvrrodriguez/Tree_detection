library(sp)
library(raster)
library(rgdal)
library(shapefiles)
library(spatstat)
library(maptools)

all.plot = list.files("~/Documentos/Datos NO publicados/BioIntForest/Data/GIS/",pattern = "^parcela.*\\poli.shp$")
all.plot = gsub(".shp","" , all.plot ,ignore.case = TRUE)

all.mdt = list.files("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR/",pattern = "^MDTParcela.*\\.tif") 
all.mds = list.files("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR/",pattern = "^MDS-MDTParcela.*\\.tif") 

rotation <- c(1, -1, 1, -1, -1, -1, -1, -1, -1) #rotate according each north bearing

for (j in 1:length(all.mds)) {

  # Importing files ---------------------------------------------------------
  
  ETRS89 <- CRS("+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  MDT <- raster(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR/", all.mdt[j]), short=TRUE)
  MDS <- raster(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR/", all.mds[j]), short=TRUE)
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/GIS")
  plot   <- readOGR(dsn=".",layer=all.plot[j])
  
  if (projection(MDT) != as.character(ETRS89)) MDT <- spTransform(MDT, ETRS89)
  if (projection(MDS) != as.character(ETRS89)) MDS <- spTransform(MDS, ETRS89)
  if (projection(plot) != as.character(ETRS89)) plot <- spTransform(plot, CRS = ETRS89) #2, 3, 5, 7, 8, 9 son WGS84
  
  source("~/Documentos/Datos NO publicados/BioIntForest/Analysis/function proj.relative.R")
  ras.MDTr <- rel.proj.raster(MDT, plot, rotation[j])
  ras.MDSr <- rel.proj.raster(MDS, plot, rotation[j])
  plot(ras.MDTr)
  

  im.MDTr <- shift.im(as.im.RasterLayer(ras.MDTr), origin = "bottomleft")
  im.MDSr <- shift.im(as.im.RasterLayer(ras.MDSr), origin = "bottomleft")
  im.MDTS <- listof(MDT = im.MDTr, MDS = im.MDSr)
  
  ras.MDTS <- brick(ras.MDTr, ras.MDSr)
  
  png(paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/MDS-MDT_A", j ,".png"), width = 1000, height = 1000)
  plot(im.MDTS, main = paste0("Plot A", j, " // Size: ", round(a, 2), " x ", round(b, 2), "m"))
  dev.off()
  
  save(im.MDTS, ras.MDTS, file = paste0("~/Documentos/Datos NO publicados/BioIntForest/Data/Environ_Maps/MDS-MDT_A",j,".RData"))
   
}
