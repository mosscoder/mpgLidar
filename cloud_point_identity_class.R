library(lidR)
library(future)
library(parallel)
library(mapview)
library(raster)

setwd('~/mpgPostdoc/projects/bareEarth/data/')

#### First identify outliers using TPI style calculation ####
ctg <- catalog('./sourceCloud/') #note this has a .lax to speed up processing

#Set tiling schema
opt_chunk_buffer(ctg) <- 30
opt_chunk_size(ctg) <- 400
#opt_select(ctg)       <- "xyzr"  
opt_output_files(ctg) <- "./outDetect/{ID}"

#Set parallel processing parameters
plan('multisession', workers = 26L)
set_lidr_threads(26L)

#Fun to filter by tpi (also scan intensity)

tpi <- function(las, ngb) {       
  lmean <- grid_metrics(las, mean(Z), ngb) #calc mean in specified neighborhood
  MeanAddLas <- lasmergespatial(las = las, source = lmean, attribute = "ngbMean") #extract neighborhood val from grid
  tpi <- MeanAddLas$Z - MeanAddLas$ngbMean #subtract neighborhood mean from focal point
  MeanAddLas <- lasadddata(MeanAddLas, tpi, 'tpi') #add tpi to cloud
  tpiQuants <- quantile(tpi, probs = c(.025,.975)) #Upper and lower quantiles of tpi
  #iQuant <- quantile(MeanAddLas$Intensity, probs = 0.025) #Lower quantile of scan intensity (2.5 percentile)
  tpiLowCut <- lasfilter(MeanAddLas, tpi > tpiQuants[1]) #cut extreme pits below 2.5 percentile tpi
  tpiHiCut <- lasfilter(tpiLowCut, tpi < tpiQuants[2]) #cut extreme peaks above 97.5 percentil tpi
  #iLowCut <- lasfilter(tpiHiCut, tpiHiCut$Intensity > iQuant) #cut low scan intensities
  tpiHiCut$tpi <- NULL #delete the tpi data from cloud
  return(tpiHiCut)
}

tpi.LAScluster = function(las, ngb) {
  las <- readLAS(las)                          # Read the LAScluster
  if (is.empty(las)) return(NULL)              # Exit early (see documentation)
  
  las <- tpi(las, ngb)      # calc tpi
  las <- lasfilter(las, buffer == 0) # Don't forget to remove the buffer  
  if(is.empty(las)) return(NULL)  
  return(las)                                  # Return the filtered point cloud
}

options <- list(
  need_output_file = TRUE,    # Throw an error if no output template is provided
  need_buffer = TRUE,         # Throw an error if buffer is 0
  automerge = TRUE)           # Automatically merge the output list (here into a LAScatalog)

#Apply noise filter to the catalog considering 1m neighborhood tpi
ctgClean  <- catalog_apply(ctg, tpi.LAScluster, ngb = 1, .options = options)
lidR:::catalog_laxindex(ctgClean) # generate .lax files to speed up process

### Classify points as canopy or ground with cloth simulator ####
opt_chunk_buffer(ctgClean) <- 30
opt_chunk_size(ctgClean) <- 400
#opt_select(ctg)       <- "xyzr"  
opt_output_files(ctgClean) <- "./classed/{ID}"

c <- csf(sloop_smooth = TRUE, class_threshold = 0.1, cloth_resolution = 0.65)

grd <- lasground(ctgClean, algorithm = c)
#grd <- catalog('./classed/')

opt_chunk_buffer(grd) <- 30
opt_chunk_size(grd) <- 400
lidR:::catalog_laxindex(grd)
opt_output_files(grd) <- "./normalized/{ID}"

nrm <- lasnormalize(las = grd, algorithm = tin(), na.rm = T)
lidR:::catalog_laxindex(nrm)

#### Testing tuning parameters of CSF algorithm ####
#Found that no disadvantage to class threshold set low, 
#Too fine of cloth resolution failed to classify dense canopy as canopy, 
#Too coarse rendered pixelated canopy mask
#.65 (meters?) best compromise
# test <- readLAS('./classed/323.las')
# 
# set_lidr_threads(1L)
# clothTester <- function(cr){
#   #cr <- g$cr[x]
#   #ct <- g$ct[x]
#   c <- csf(sloop_smooth = T, class_threshold = 0.1, cloth_resolution = cr,iterations = 1000)
#   grnd <- lasground(test, algorithm = c)
#   writeLAS(grnd, paste0('./clothTest/','cr',cr,'.las'))
# }
# 
# crs <- seq(from = .5, to =1, by = 0.05)
# cts <- seq(from = .1, to =.5, by = 0.1)
# 
# g <- expand.grid(crs, cts)
# colnames(g) <- list('cr','ct')
# 
# mclapply(FUN = clothTester, X =crs, mc.cores = length(crs))

#### Test interp alg ####

# test <-readLAS('./classed/323.las')
# 
# set_lidr_threads(1L)
# 
# # test knn function
# 
# neighbors <- 10:20
# 
# krigTest <- function(n){
# 
#   ras <- grid_terrain(test, res = 0.15, algorithm = kriging(k = n))
#   slp <- terrain(ras)
#   asp <- terrain(ras, opt='aspect')
#   hill <- hillShade(slope = slp, aspect = asp, angle = 30, direction = 215)
#   
# }
# 
# st <- stack(mclapply(FUN = krigTest, X = neighbors, mc.cores = 11))
# 
# install.packages('leafsync')
# 
# lapply(X = st, FUN = function(x){mapview(st[[x]])}) %>% sync()




