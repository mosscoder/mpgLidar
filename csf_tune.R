library(sp)
library(rgdal)
library(lidR)
library(tidyverse)
library(future)
library(parallel)
library(ranger)

setwd('~/mpgPostdoc/projects/bareEarth/data/')

#Loading test points, focusing on shrub sensitivity in this analysis
shrubPts <- readOGR(dsn = path.expand('./test_points/shrub_points_utm.shp'), layer = 'shrub_points_utm')@coords %>%
  as.data.frame() %>% mutate(type = 'shrub')

barePts <- readOGR(dsn = path.expand('./test_points/bare_slope_points.shp'), layer = 'bare_slope_points')@coords %>%
  as.data.frame() %>% mutate(type = 'bare')

pts <- rbind(shrubPts, barePts)

colnames(pts)[1:2] <- c('x','y')

#Loading sfm point cloud
ctg <- catalog('./sourceCloud/')

#Setting environment for the cloud
opt_select(ctg) <- 'xyz'
plan('multisession', workers = 20L)

#Extracting cloud data around test points
# clipd <- lasclipCircle(ctg, xcenter = pts$x, ycenter = pts$y, radius = 5) 
# lapply(FUN = function(x){writeLAS(clipd[[x]], paste0('./test_points/subCloud/',x,'.las'))}, X = 1:nrow(pts))
clipd <- lapply(FUN = function(x){readLAS(paste0('./test_points/subCloud/',x,'.las'))}, X = 1:nrow(pts))

#Identifying nearest pixel in cloud to test points
eucFocal <- function(x){
  focalX <- pts[x,1]
  focalY <- pts[x,2]
  cloudX <- clipd[[x]]$X
  cloudY <- clipd[[x]]$Y
  eD <- sqrt((focalX - cloudX)^2 + (focalY - cloudY)^2)
  which.min(eD)
}

closeCldPt <- sapply(FUN = eucFocal, X = 1:nrow(pts))

#Grid along which we will test cloth resolution
res <- seq(from = .3, to = 2.5, by = 0.1)

groundTruth <- function(x){ #Function to extract sensitivity of the cloth res defined above
  fPoint <- closeCldPt[x] #Grab the closest cloud point
  fCloud <- clipd[[x]] #Grab the associated sub-cloud
  type <- pts$type[x] #Grab the type of point (shrub or bare)
  
  grounder <- function(y){ # get the Classifiation of the test point for a given res
    
    pt <- lasground(fCloud,algorithm = csf(sloop_smooth = T, 
                                           class_threshold = 0.1,
                                           cloth_resolution = res[y]
    ))$Classification[fPoint] 
    
    if(type == 'shrub' & pt == 2){ truth <- 0} #if shrub classed as ground, set to zero
    if(type == 'shrub' & pt == 1){ truth <- 1} #if shrub classed as canopy, set to 1
    if(type == 'bare' & pt == 2){ truth <- 1} #if ground classed as ground, set to 1
    if(type == 'bare' & pt == 1){ truth <- 0} #if ground classed as canopy, set to 0
    
    df <- data.frame(#wrap up result into a data.frame
      id = x, 
      type = type,
      res = res[y], 
      truth = truth
    )
    
    df
    
  }
  
  do.call(rbind, lapply(FUN = grounder, X = 1:length(res))) #bind results together across the res grid search
  
}
groundTruth(500) #test one point

#run the funciton for all test points
onePass <- function(x){ do.call(rbind, mclapply(FUN = groundTruth, X = 1:nrow(pts), mc.cores = 25)) }
n <- onePass() #%>% group_by(res, rig, time, type) %>% summarize(sum = sum(truth))

ggplot(n, aes(x = res, y = time, fill = sum/(nrow(pts)/2))) +
  geom_raster() + 
  scale_fill_viridis_c(option = 'D', name = "Sensitivity") +
  facet_wrap(~type) +
  xlab('Cloth Resolution (m)') +
  ylab('Time Step (seconds)')

#Model sensitivity based on z, slope, and tpi with 100m radius using 10m dsm
opt_output_files(ctg) <- './10mDSM/{ID}'
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 30

#grid_metrics(ctg, func = mean(Z), res = 10)

#load the raster after filling holes in gdal
dsm <- raster('./10mDSM/10mDSMFill.tif')

#calc slope
slope <- terrain(dsm, option = 'slope')

#calc tpi
tpiw <- function(x, w=5) {
  m <- matrix(1/(w^2-1), nc=w, nr=w)
  m[ceiling(0.5 * length(m))] <- 0
  f <- focal(x, m)
  x - f
}

tpi100 <- tpiw(dsm, w = 21)

#stack of rasters
st <- stack(dsm, slope, tpi100)
names(st) <- c('z', 'slope', 'tpi100')

#extract topo data at test points
shrubExt <- raster::extract(x = st, y = shrubPts[,1:2])
bareExt <- raster::extract(x = st, y = barePts[,1:2])

envDat <- data.frame(id = 1:nrow(pts), rbind(shrubExt, bareExt))

#join topo data to truth data
n2 <- left_join(n, envDat) 

#model truth as a function of resolution, z, slope, and tpi
mod <- ranger(truth %>% as.factor() ~ ., 
              data = n2 %>% select(truth,res, z, slope, tpi100), 
              importance = 'permutation', 
              num.trees = 100) #low number for speed, gives same good OOB

#turn our topo raster data into a data.frame
rasPts <- rasterToPoints(st[[c('z','slope','tpi100')]]) %>% as.data.frame() %>% na.omit()
colnames(rasPts)[1:2] <- c('x','y')

resPicker <- function(x){ #a function to test a range of cloth resolutions, identify which gives the correct class
  fz <- rasPts$z[x] #test point z
  fslope <- rasPts$slope[x] #test point slope
  ftpi <- rasPts$tpi[x] #test point tpi
  
  res <- numeric(1L) #create empty scalar to receive res in for loop
  
  for(i in rev(seq(from = .3, to = 2.5, by = 0.1))) { #check truth across res grid, break loop when it truth is 1 (correct)
    pred <- predict(mod,
                    data.frame(res = i,
                               z = fz,
                               slope = fslope,
                               tpi100 = ftpi
                    ),
                    num.threads = 1
    )$predictions
    res[1] <- i
    if(pred == 1){break}
    
  }
  
  data.frame(x = rasPts$x[x], #wrap the best res sesult in data.frame
             y = rasPts$y[x],
             res = res)
}

#predict the appropriate cloth res across the landscape
resPicked <- do.call(rbind, mclapply(FUN = resPicker, X = 1:nrow(rasPts), mc.cores = 20))

#load a dummy raster with the predicted res
resRas <- raster('./10mDSM/10mDSM.tif')
nas <- which(is.na(raster::values(resRas)))
values(resRas) <- NA
cells <- cellFromXY(resRas, resPicked)
resPickedClean <- resPicked %>% mutate(cells = cells) %>% filter(! cells %in% nas)

raster::values(resRas)[resPickedClean$cells ] <- resPickedClean$res

writeRaster(resRas, './clothRes/clothRes.tif', format ='GTiff', overwrite = T)
plot(resRas)

ggplot(resPicked %>% mutate(cells = cells) %>% filter(! cells %in% nas), aes(x = x, y = y, fill = res)) +
  geom_raster() +
  scale_fill_viridis_c()


