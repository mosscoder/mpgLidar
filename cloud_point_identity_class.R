library(lidR)
library(future)
library(parallel)
library(mapview)
library(raster)

ctg <- catalog('~/mpgPostdoc/projects/bareEarth/data/sourceCloud/')

opt_chunk_buffer(ctg) <- 30
opt_chunk_size(ctg) <- 300
opt_output_files(ctg) <- "~/mpgPostdoc/projects/bareEarth/data/classed/{ID}"

plan('multisession', workers = 26L)
set_lidr_threads(26L)

c <- csf(sloop_smooth = TRUE, class_threshold = 0.1, cloth_resolution = 0.65)

grd <- lasground(ctg, algorithm = c)
#grd <- catalog('~/mpgPostdoc/projects/bareEarth/data/classed/')

opt_chunk_buffer(grd) <- 30
opt_chunk_size(grd) <- 300
#lidR:::catalog_laxindex(grd)
opt_output_files(grd) <- "~/mpgPostdoc/projects/bareEarth/data/normalized/{ID}"

nrm <- lasnormalize(las = grd, algorithm = knnidw(), na.rm = T)

opt_chunk_buffer(nrm) <- 30
opt_chunk_size(nrm) <- 300
#lidR:::catalog_laxindex(nrm)
opt_output_files(nrm) <- "~/mpgPostdoc/projects/bareEarth/data/terrain/{ID}"

dtm <- grid_terrain(nrm, res = 0.15, algorithm = knnidw())

plot(grd, mapview = T, map.type = "Esri.WorldImagery")

#### Testing tuning parameters of CSF algorithm ####
#Found that no disadvantage to class threshold set low, 
#Too fine of cloth resolution failed to classify dense canopy as canopy, 
#Too coarse rendered pixelated canopy mask
#.65 (meters?) best compromise
# test <- readLAS('~/mpgPostdoc/projects/bareEarth/data/classed/323.las')
# 
# set_lidr_threads(1L)
# clothTester <- function(cr){
#   #cr <- g$cr[x]
#   #ct <- g$ct[x]
#   c <- csf(sloop_smooth = T, class_threshold = 0.1, cloth_resolution = cr,iterations = 1000)
#   grnd <- lasground(test, algorithm = c)
#   writeLAS(grnd, paste0('~/mpgPostdoc/projects/bareEarth/data/clothTest/','cr',cr,'.las'))
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

test <-readLAS('~/mpgPostdoc/projects/bareEarth/data/classed/323.las')

set_lidr_threads(1L)

# test knn function

neighbors <- 10:20

krigTest <- function(n){

  ras <- grid_terrain(test, res = 0.15, algorithm = kriging(k = n))
  slp <- terrain(ras)
  asp <- terrain(ras, opt='aspect')
  hill <- hillShade(slope = slp, aspect = asp, angle = 30, direction = 215)
  
}

st <- stack(mclapply(FUN = krigTest, X = neighbors, mc.cores = 11))

install.packages('leafsync')

lapply(X = st, FUN = function(x){mapview(st[[x]])}) %>% sync()




