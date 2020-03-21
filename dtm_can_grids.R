library(lidR)
library(future)
library(parallel)
library(mapview)
library(raster)

grd <- catalog('~/mpgPostdoc/projects/bareEarth/data/classed/')

opt_chunk_buffer(grd) <- 30
opt_chunk_size(grd) <- 300
opt_select(grd) <- 'xyzc'
opt_filter(grd) <- '-keep_class 2 9'
opt_output_files(grd) <- "~/mpgPostdoc/projects/bareEarth/data/terrain/{ID}"

plan('multisession', workers = 26L)
set_lidr_threads(26L)

ter <- grid_terrain(grd, res = 0.15, algorithm = knnidw())

nrml <- catalog('~/mpgPostdoc/projects/bareEarth/data/normalized/')

opt_chunk_buffer(nrml) <- 30
opt_chunk_size(nrml) <- 300
opt_select(nrml) <- 'xyzc'
opt_output_files(nrml) <- "~/mpgPostdoc/projects/bareEarth/data/canopy/{ID}"
canopy <- grid_canopy(nrml, 0.15, pitfree(c(0,2,5,10,15), c(0,1), subcircle = 0.2))