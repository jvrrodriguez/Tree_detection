rm(list = ls())

# Add required packages 
require(sp)
require(rgdal)
require(spatstat)
library(maptools)
library(spatstat.geom)

library(raster)
library(sp)
library(gstat)
library(rgeos)

data_LidarT.all <- readOGR("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR/AIDFOREST-ASPURZ-PLOTS-2021/trees.shp")
data_LidarT.boxesall <- readOGR("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR/AIDFOREST-ASPURZ-PLOTS-2021/BOXES.shp")

save.out <- TRUE

#la parcela 7 no termina de coger las marcas en ppp_1999

for (i in 1:9){


  # Preparación de los datos ------------------------------------------------

  setwd("~/Documentos/Datos NO publicados/BioIntForest/Tree_detection_Results")
  if (save.out == T) pdf(paste0("Plot",i,"_Match_trees.pdf"), width=7, height=5, compress = TRUE)
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Inventarios_floristicos")
  traits_1999 <- read.csv("Medidas pinos 1999.csv", sep=";", dec=",")
  traits_1999 <- traits_1999[traits_1999$Parcela == paste0("A",i),]
  traits_1999$Height <- rowMeans(cbind(traits_1999$P.MeanHt, traits_1999$A.MeanHt), na.rm=T)
  
  traits_2009 <- read.csv("PS09Medidas pino 2009 con coordenadas.csv", sep=";", dec=",")
  traits_2009 <- traits_2009[traits_2009$PARCELA == paste0("A",i),]
  # en este año no encuentro alturas que se hayan medido
  
  # solo para árboles testigo y no todas las parcelas (3, 4, 9)
  traits_2014 <- read.csv("Aspurz inventory 2014.csv", sep=";", dec=".")
  traits_2014 <- traits_2014[traits_2014$Plot == paste0("A",i),]
  traits_2014$Height <- traits_2014$Ht_14
  # no se que diferencia hay con la anterior base de datos
  traits_2014b <- read.csv("Alturas dominantes_2014.csv", sep=";", dec=".")
  traits_2014b <- traits_2014b[traits_2014b$Plot == paste0("A",i),]
  traits_2014b$Height <- traits_2014b$Ht_14
  
  traits_2014 <- traits_2014b #selecciono 1 de las bases de datos
  
  traits_2018 <- read.csv("aspurz inventory 2018_alturas.csv", sep=";", dec=".")
  traits_2018 <- traits_2018[traits_2018$PARCELA == i,]
  traits_2018$Height <- traits_2018$H.PROMEDIO..m.
  
  traits_2021 <- read.csv("Aspurz inventario 2021.csv", sep=",", dec=".")
  traits_2021 <- traits_2021[traits_2021$PARCELA == i,]
  traits_2021.height <- read.csv("Aspurz inventario 2021_alturas.csv", sep=",", dec=".")
  traits_2021.height <- traits_2021.height[traits_2021.height$Parcela == i,]
  traits_2021 <- cbind(traits_2021, traits_2021.height[match(traits_2021$ARBOL, traits_2021.height$Arbol),c(3,10)])
  traits_2021$Height <- traits_2021$X..Altura
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Mapas Arboles")
  
  if (i == 3 | i == 4 | i == 5) {
    
    Studyarea <- readOGR(paste0("parcela ",i, "b poli.shp"))
  
  } else {
    
    Studyarea <- readOGR(paste0("parcela ",i, " poli.shp"))
    
  }
  
  Studyarea <- spTransform(Studyarea, "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs")

  # Coerce study area to owin object
  w <- as.owin(Studyarea)
  
  if (i == 1 | i == 7) {
    
    data_1999 <- readOGR(paste0("Pinos_vivos_P",i,"_1999.shp"))
    
  } else {
    
    data_1999 <- readOGR(paste0("Arboles_P",i,"_1999.shp"))
    
  }
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/Mapas Arboles")
  data_2007 <- readOGR(paste0("Pinos_vivos_P",i,"_2007.shp"))
  data_2007@data$id <- as.numeric(data_2007@data$id)
  data_2021 <- data_2007[match(traits_2021$ARBOL, as.numeric(data_2007$id), nomatch = 0),]#readOGR("Pinos_vivos_P7_2020.shp")
  data_LidarT <- data_LidarT.all[data_LidarT.all$PLOT == paste0("ASP-P-",i,".las"),]
  
  #calcular el punto más alto detectado con el terrestre para compararlo con el aéreo
  
  data_LidarT.boxes <- data_LidarT.boxesall[data_LidarT.boxesall$PLOT == paste0("ASP-P-",i,".las"),]
  data_LidarT.boxes.max <- data_LidarT.boxes[data_LidarT.boxes$ID_GROUP > 200,] #crear un spatial polygons vacio
  ID_TREE <- sort(unique(data_LidarT.boxes$ID_GROUP))
  ID_TREE <- ID_TREE[ID_TREE > 0]
  
  for (z in ID_TREE) {
    
    data.group <- data_LidarT.boxes[data_LidarT.boxes$ID_GROUP == ID_TREE[z],]
    data_LidarT.boxes.max <- rbind(data_LidarT.boxes.max, data.group[data.group$Z == max(data.group$Z),][1,])
    
  }
  
  data_LidarT.top = as(gCentroid(data_LidarT.boxes.max,byid=TRUE),"SpatialPointsDataFrame")
  data_LidarT.top@data <- data_LidarT.boxes.max@data
  
  
  
  setwd("~/Documentos/Datos NO publicados/BioIntForest/Data/LiDAR")
  load(paste0("LiDAR_A",i,".RData"))
  #dsm[dsm < 10] <- 10 # clasificar como 10 metros los valores menores de ese umbral
  
  ## select methodology to locate pine centroids from airborne LiDAR
  data_LidarA <- id.trees.clust
  #data_LidarA <- id.trees.tops
  
  if (i==1) {
    
    #ppp_1999 <- with(data_1999, ppp(x = data_1999@coords[,1], y = data_1999@coords[,2], owin(w)))
    ppp_1999 <- with(data_1999, ppp(x = data_1999@coords[,1], y = data_1999@coords[,2], 
                                    marks = data.frame(id = data_1999$etiqueta,
                                                       Size = traits_1999$DM[match(data_1999@data$etiqueta, traits_1999$cod1)], 
                                                        Height = traits_1999$Height[match(data_1999@data$etiqueta, traits_1999$cod1)]), owin(w)))
    
  } else {
    
    ppp_1999 <- with(data_1999, ppp(x = data_1999@coords[,1], y = data_1999@coords[,2], 
                                    marks = data.frame(id = data_1999$Etiqueta,
                                                       Size = traits_1999$DM[match(data_1999@data$Etiqueta, traits_1999$cod1)], 
                                                        Height = traits_1999$Height[match(data_1999@data$Etiqueta, traits_1999$cod1)]), owin(w)))
    
  }
  
  if (i==7) {
    
    ppp_1999 <- with(data_1999, ppp(x = data_1999@coords[,1], y = data_1999@coords[,2], 
                                    marks = data.frame(id = data_1999$ID,
                                                       Size = traits_1999$DM[match(data_1999@data$ID, traits_1999$cod1)],
                                                        Height = traits_1999$Height[match(data_1999@data$ID, traits_1999$cod1)]), owin(w)))
    
  }
  
  ppp_2007 <- with(data_2007, ppp(x = data_2007@coords[,1], y = data_2007@coords[,2], 
                                  marks = data.frame(id = data_2007$id,
                                  Size = traits_2009$Dn_09[match(data_2007@data$id, traits_2009$ARBOL)],
                                  Height = traits_2014$Height[match(data_2007@data$id, traits_2014$Number)], # altura para 2014
                                  Height.b = traits_2018$Height[match(data_2007@data$id, traits_2018$ARBOL)], #altura para 2018
                                  Height.dsm = extract(dsm, data_2007)), owin(w))) 
  
  ppp_2021 <- with(data_2021, ppp(x = data_2021@coords[,1], y = data_2021@coords[,2], 
                                   marks = data.frame(id = data_2021$id,
                                                      Size = traits_2021$DIAM.Avg[match(data_2021$id, traits_2021$ARBOL)], 
                                                      Height = traits_2021$Height[match(data_2021$id, traits_2021$ARBOL)],
                                                      Height.dsm = extract(dsm, data_2021)), owin(w)))
  
  ppp_lidarT <- with(data_LidarT, ppp(x = data_LidarT@coords[,1], y = data_LidarT@coords[,2], 
                     marks = data.frame(id = data_LidarT$ID_TREE, Size = data_LidarT$DBH, 
                                        Height = data_LidarT$Ht, 
                                        Height.dsm = extract(dsm, data_LidarT)), owin(w)))
  
  ppp_lidarT.top <- with(data_LidarT.top, ppp(x = data_LidarT.top@coords[,1], y = data_LidarT.top@coords[,2], 
                                      marks = data.frame(id = data_LidarT.top$ID_GROUP, Size = rowMeans(cbind(data_LidarT.top$D1_WIDTH, data_LidarT.top$D2_HEIGHT)), 
                                                         Height = data_LidarT.top$Z, 
                                                         Height.dsm = extract(dsm, data_LidarT.top)), owin(w)))
  
  
  # las variables hay que cambiarlas dependiendo de como se calcula la posiciónd e los árboles con LiDAR aereo
  ppp_lidarA <- with(data_LidarA, ppp(x = data_LidarA@coords[,1], y = data_LidarA@coords[,2], 
                                      marks = data.frame(id = as.numeric(data_LidarA@data$Tree), 
                                        Size = data_LidarA@data$crownvolume, 
                                        Height.dsm = data_LidarA@data$GPA), owin(w)))
  
  # Función para detectar si el mismo número se encuentra proximo, y si es así a cuanta distancia
  # Esto solo para censos manuales
  

  # Relaciones de proximidad entre metodologías -----------------------------

  nnid <- function (X, Y, k){
    
    #si no lo calculo por separado, las distancias me dan valor infinito
    nn_XtoY.dist <- nncross(X=X,Y=Y, k=1:k, what=c("dist"))
    nn_XtoY.which <- nncross(X=X,Y=Y, k=1:k, what=c("which"))
    nn_XtoY <- cbind(nn_XtoY.dist, nn_XtoY.which)
    
    for (j in 1:k) {
      
      nn_XtoY <- cbind(nn_XtoY, Y[nn_XtoY[[k + j]]]$marks)
      colnames(nn_XtoY)[2*k + j] <- paste("id.", j, sep="")
      
    }
    
    dist.all <- nn_XtoY[,1:k]
    which.all <- nn_XtoY[,(2*k+1):(2*k+k)]
    
    #Para convertir en array con algunos elementos con cero (named integer(0))
    kk <- unlist(lapply(apply(which.all == X$marks, 1, which), function (x) ifelse(length (x) > 0, x, NA)))
    
    kkk <- array(NA, length(kk))
    
    for (i in 1:length(kk)) {
      
      if (is.finite(kk[i])) {
        
        kkk[i] <- dist.all[i, kk[i]]
      }
      
    }

    kkk <- data.frame(which=kk, dist=ifelse(is.finite(kkk), kkk, NA)) 
    
    return(kkk)

  }
  
  nn_07to99 <- data.frame(id07 = ppp_2007$marks$id, nn_07to99 = nnid(X=subset.ppp(ppp_2007, select = "id"), Y=subset.ppp(ppp_1999, select = "id"), k=20))
  
  #hist(nn_07to99)
  #nn_07to99[!is.na(nn_07to99$nn_07to99.dist),]
  
  par(mfrow=c(1,2))
  #par(mfrow=c(2,2), mar=c(2, 2, 3, 3), oma = c(2, 2, 1, 1))
  plot(subset.ppp(ppp_1999, select = "Size"), main = 1999)
  plot(subset.ppp(ppp_2007, select = "Size"), main = 2007)
  
  hist(nn_07to99$nn_07to99.which, main = "", xlab = "Closest label matching")
  hist(nn_07to99$nn_07to99.dist, main = "", xlab = "Dist. label matching")
  mtext(paste0("Plot", i, " // Spatial matching of tree labels - 2007 vs 1999"), side = 3, line = -2, outer = TRUE) 
  
  plot(subset.ppp(ppp_2021, select = "Size"), main = 2021)
  plot(subset.ppp(ppp_lidarT, select = "Size"), main = "LiDAR.T")
  
  nn_lTto21 <- nncross(X=ppp_lidarT,Y=ppp_2021, what=c("dist","which"))
  nn_lTto21$nn_DBH <- ppp_2021$marks$Size[nn_lTto21$which] # el DBH del árbol más cercano del 2021
  
  hist(table(nn_lTto21$which), main = "", xlab = "Closest label matching")
  hist(nn_lTto21$dist, main = "", xlab = "Dist. label matching")
  mtext(paste0("Plot", i, " // Spatial matching of tree labels - 2021 vs Lidar.T"), side = 3, line = -2, outer = TRUE) 

  plot(subset.ppp(ppp_2021, select = "Size"), main = 2021)
  plot(subset.ppp(ppp_lidarA, select = "Size"), main = "LiDAR.A")
  
  nn_lAto21 <- nncross(X=ppp_lidarA,Y=ppp_2021, what=c("dist","which"))
  nn_lAto21$nn_DBH <- ppp_2021$marks$Size[nn_lAto21$which] # el DBH del árbol más cercano del 2021
  
  hist(table(nn_lAto21$which), main = "", xlab = "Closest label matching")
  hist(nn_lAto21$dist, main = "", xlab = "Dist. label matching")
  mtext(paste0("Plot", i, " // Spatial matching of tree labels - 2021 vs Lidar.A"), side = 3, line = -2, outer = TRUE) 
  
  
  # Comparación de rasgos ---------------------------------------------------

  par(mfrow=c(1,2))
  plot(ppp_lidarT$marks$Size, nn_lTto21$nn_DBH, main = "lidar.T", xlab = "dbh", ylab = "nearest.dbh observed")
  plot(ppp_lidarA$marks$Size, nn_lAto21$nn_DBH, main = "lidar.A", xlab = "crown area", ylab = "nearest.dbh observed")

  plot(ppp_2021$marks$Height, ppp_2021$marks$Height.dsm, main = "data 2021", xlab = "height", ylab = "height.dsm") 
  plot(ppp_lidarT$marks$Height, ppp_lidarT$marks$Height.dsm, main = "lidar.T", xlab = "height", ylab = "height.dsm") 
  
  plot(ppp_2007$marks$Height, ppp_2007$marks$Height.dsm, main = "data 2014", xlab = "height", ylab = "height.dsm") 
  plot(ppp_2007$marks$Height.b, ppp_2007$marks$Height.dsm, main = "data 2018", xlab = "height", ylab = "height.dsm")
  
  plot(ppp_lidarT.top$marks$Height, ppp_lidarT.top$marks$Height.dsm, main = "lidar.T", xlab = "height", ylab = "height.dsm") 
  
  
  dens_2021 <- density(ppp_2021, sigma=bw.diggle(ppp_2021))
  
  #https://www.emilyburchfield.org/courses/gsa/point_pattern_lab
  par(mfrow=c(1,3))
  plot(dens_2021, useRaster=FALSE)
  plot(ppp_lidarT, add = TRUE, cols = "white", cex = 0.5, pch = 16)
  plot(rhohat(unmark(ppp_lidarT), dens_2021))
  ks <- cdf.test(unmark(ppp_lidarT), dens_2021)
  plot(ks)
  
  par(mfrow=c(1,3))
  plot(dens_2021, useRaster=FALSE)
  plot(ppp_lidarA, add = TRUE, cols = "white", cex = 0.5, pch = 16)
  plot(rhohat(unmark(ppp_lidarA), dens_2021))
  ks <- cdf.test(unmark(ppp_lidarA), dens_2021)
  plot(ks)
  
  Plot.grid <- as(raster(ncol = 15*3, nrow = 20*3, crs=NULL, ext = extent(Studyarea@bbox)), "SpatialPolygons")
  Plot.grid <- SpatialPolygonsDataFrame(Plot.grid, data.frame(n=1:length(Plot.grid)))
  
  size.vgm <- variogram(marks.Size ~ 1, as.SpatialPointsDataFrame.ppp(ppp_2021)) 
  size.fit <- fit.variogram(size.vgm, model=vgm("Exp", "Sph", "Gau", "Mat")) # fit model
  plot(size.vgm, size.fit)
  size.kriged <- krige(marks.Size ~ 1, as.SpatialPointsDataFrame.ppp(ppp_2021), Plot.grid, model=size.fit)
  
  source("~/Documentos/Datos NO publicados/BioIntForest/Tree_detection/function Raster_extractGrid.R")
  dens_size2021 <- Raster_extractGrid(size.kriged, colnames(size.kriged@data)[grep(".pred", colnames(size.kriged@data))], 0.5, plot = F)
  rasValue.lidarT <- extract(dens_size2021, as.SpatialPointsDataFrame.ppp(ppp_lidarT))
  rasValue.lidarA <- extract(dens_size2021, as.SpatialPointsDataFrame.ppp(ppp_lidarA))
  
  par(mfrow=c(2,2))
  plot(as.SpatialPointsDataFrame.ppp(ppp_lidarT)$marks.Size, rasValue.lidarT, main = "lidar.T", xlab = "DBH ", ylab = "pred. size from obs.")
  plot(as.SpatialPointsDataFrame.ppp(ppp_lidarA)$marks.Size, rasValue.lidarA, main = "lidar.A", xlab = "crown area", ylab = "pred. DBH from obs.")
  plot(as.SpatialPointsDataFrame.ppp(ppp_lidarA)$marks.Height.dsm, rasValue.lidarA, main = "lidar.A", xlab = "height.dsm", ylab = "pred. DBH from obs.")
  
  dens_size2021 <- Smooth.ppp(subset.ppp(ppp_2021, select = "Size"), sigma=2)
  
  par(mfrow=c(1,2))
  ks <- cdf.test(unmark(ppp_lidarT), (dens_size2021)) #as.im.RasterLayer(dens_size2021)
  plot(ks)
  ks <- cdf.test(unmark(ppp_lidarA), (dens_size2021))
  plot(ks)
  
  par(mfrow=c(1,3))
  plot(as.im.RasterLayer(dsm), useRaster=FALSE)
  plot(ppp_2021, add = TRUE, cols = "white", cex = 0.5, pch = 16)
  plot(rhohat(unmark(ppp_2021), as.im.RasterLayer(dsm)))
  ks <- cdf.test(unmark(ppp_2021), as.im.RasterLayer(dsm))
  plot(ks)
  
  par(mfrow=c(1,3))
  plot(as.im.RasterLayer(dsm), useRaster=FALSE)
  plot(ppp_lidarT, add = TRUE, cols = "white", cex = 0.5, pch = 16)
  plot(rhohat(unmark(ppp_lidarT), as.im.RasterLayer(dsm)))
  ks <- cdf.test(unmark(ppp_lidarT), as.im.RasterLayer(dsm))
  plot(ks)
  
  par(mfrow=c(1,3))
  plot(as.im.RasterLayer(dsm), useRaster=FALSE)
  plot(ppp_lidarT.top, add = TRUE, cols = "white", cex = 0.5, pch = 16)
  plot(rhohat(unmark(ppp_lidarT.top), as.im.RasterLayer(dsm)))
  ks <- cdf.test(unmark(ppp_lidarT.top), as.im.RasterLayer(dsm))
  plot(ks)
  
  
  par(mfrow=c(1,3))
  plot(distfun(ppp_lidarA), useRaster=FALSE)
  plot(ppp_2021, add = TRUE, cols = "white", cex = 0.5, pch = 16)
  plot(rhohat(unmark(ppp_2021), distfun(ppp_lidarA)))
  ks <- cdf.test(unmark(ppp_2021), distfun(ppp_lidarA))
  plot(ks)
  
  
  
  # Agregación espacial entre muestreos y metodologias ----------------------

  nsim = 199
  
  L.E99 <- envelope(ppp_1999, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
  g.E99 <- envelope(ppp_1999, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
  #kNN.E99 <- envelope(ppp_1999, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  #F.E99 <- envelope(ppp_1999, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  ppp_1999.size <- subset.ppp(ppp_1999, select = "Size")
  ppp_1999.size <- ppp_1999.size[!is.na(ppp_1999.size$marks)]
  markcor.E99 <- envelope(ppp_1999.size, markcorr, nsim = nsim, envelope = TRUE)
  
  L.E09 <- envelope(ppp_2007, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
  g.E09 <- envelope(ppp_2007, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
  #kNN.E09 <- envelope(ppp_2007, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  #F.E09 <- envelope(ppp_2007, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  ppp_2007.size <- subset.ppp(ppp_2007, select = "Size")
  ppp_2007.size <- ppp_2007.size[!is.na(ppp_2007.size$marks)]
  markcor.E09 <- envelope(ppp_2007.size, markcorr, nsim = nsim, envelope = TRUE)
  
  L.E21 <- envelope(ppp_2021, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
  g.E21 <- envelope(ppp_2021, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
  #kNN.E21 <- envelope(ppp_2021, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  #F.E21 <- envelope(ppp_2021, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  ppp_2021.size <- subset.ppp(ppp_2021, select = "Size")
  ppp_2021.size <- ppp_2021.size[!is.na(ppp_2021.size$marks)]
  markcor.E21 <- envelope(ppp_2021.size, markcorr, nsim = nsim, envelope = TRUE)
  
  L.ElidarT <- envelope(ppp_lidarT, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
  g.ElidarT <- envelope(ppp_lidarT, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
  #kNN.ElidarT <- envelope(ppp_lidarT, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  #F.ElidarT <- envelope(ppp_lidarT, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  ppp_lidarT.size <- subset.ppp(ppp_lidarT, select = "Size")
  markcor.ElidarT <- envelope(ppp_lidarT.size, markcorr, nsim = nsim, envelope = TRUE)
  
  L.ElidarA <- envelope(ppp_lidarA, Lest, r = seq(0,10,0.05), nsim=nsim, fix.n=TRUE, correction="Ripley")
  g.ElidarA <- envelope(ppp_lidarA, pcf, r = seq(0,4,0.02), nsim=nsim, fix.n=TRUE, correction="Ripley")
  #kNN.ElidarA <- envelope(ppp_lidarA, Gest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  #F.ElidarA <- envelope(ppp_lidarA, Fest, r = seq(0,1,0.01), nsim=nsim, fix.n=TRUE, correction="rs")
  ppp_lidarA.size <- subset.ppp(ppp_lidarA, select = "Size")
  markcor.ElidarA <- envelope(ppp_lidarA.size, markcorr, nsim = nsim, envelope = TRUE)
  
  par(mfrow=c(5,3), mar=c(1, 1, 1.25, 1.25), oma = c(4, 4, 2, 2))
  
  plot(L.E99, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
  plot(g.E99, main = NULL, legend = F)
  plot(markcor.E99, main = NULL, legend = F)
  #plot(kNN.E99, main = NULL, legend = F)
  #plot(F.E99, main = NULL, legend = F)
  
  plot(L.E09, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
  plot(g.E09, main = NULL, legend = F)
  plot(markcor.E09, main = NULL, legend = F)
  #plot(kNN.E09, main = NULL, legend = F)
  #plot(F.E09, main = NULL, legend = F)
  
  plot(L.E21, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
  plot(g.E21, main = NULL, legend = F)
  plot(markcor.E21, main = NULL, legend = F)
  #plot(kNN.E21, main = NULL, legend = F)
  #plot(F.E21, main = NULL, legend = F)
  
  plot(L.ElidarT, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
  plot(g.ElidarT, main = NULL, legend = F)
  plot(markcor.ElidarT, main = NULL, legend = F)
  #plot(kNN.ElidarT, main = NULL, legend = F)
  #plot(F.E20, main = NULL, legend = F)
 
  plot(L.ElidarA, . -r ~ r, shade=c("hi", "lo"), legend = F, main = NULL)
  plot(g.ElidarA, main = NULL, legend = F)
  plot(markcor.ElidarA, main = NULL, legend = F)
  #plot(kNN.ElidarA, main = NULL, legend = F)
  #plot(F.E20, main = NULL, legend = F)
  
  mtext(paste0("Plot", i, " // Spatial aggregation in 1999, 2009, 2021, LiDAR.T and LiDAR.A"), side = 3, line = -1, outer = TRUE) 

    #Size <- marks(ppp_2021)$Size
  #dens.size_2021 <- density(ppp_2021, weights=marks(ppp_2021), sigma=bw.ppl)
  #rhohat(ppp_2021, ~ 1)
  #pairs(list(dens_2021, density(ppp_lidarT, sigma=bw.diggle(ppp_lidarT)), density(ppp_lidarA, sigma=bw.diggle(ppp_lidarA))))
  
  values.lidar <- rbind(
    data.frame(Method = "lidar.T", Data = "Observed", P = dens_2021[ppp_lidarT]),
    data.frame(Method = "lidar.T", Data = "Random", P = dens_2021[(rpoint(ppp_lidarT$n, dens_2021))]),
    data.frame(Method = "lidar.A", Data = "Observed", P = c(dens_2021[ppp_lidarA])),
    data.frame(Method = "lidar.A", Data = "Random", P= dens_2021[(rpoint(ppp_lidarA$n, dens_2021))]))
  
  myData <- aggregate(values.lidar$P,
                      by = list(Method = values.lidar$Method, Data = values.lidar$Data),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm=T),
                                          n = length(x)))
  myData <- do.call(data.frame, myData)
  myData$x.se <- myData$x.sd / sqrt(myData$x.n)
  
  library(ggplot2)
  # Default bar plot
  p<- ggplot(myData, aes(x=Method, y=x.mean, fill=Data)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=x.mean-x.se, ymax=x.mean+x.se), width=.2,
                  position=position_dodge(.9)) 
  
  print(p)
  
  if (save.out ==T) dev.off()
  
}
