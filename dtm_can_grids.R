library(lidR)
library(future)
library(parallel)
library(mapview)
library(raster)

setwd('~/mpgPostdoc/projects/bareEarth/data/')

grd <- catalog('./classed/')

opt_chunk_buffer(grd) <- 30
opt_chunk_size(grd) <- 300 
opt_select(grd) <- 'xyzc'
opt_filter(grd) <- '-keep_class 2 9'
opt_output_files(grd) <- "./terrain/{ID}"

plan('multisession', workers = 20L)
set_lidr_threads(20L)

ter <- grid_terrain(grd, res = 0.15, algorithm = knnidw())

nrml <- catalog('./normalized/')

opt_chunk_buffer(nrml) <- 30
opt_chunk_size(nrml) <- 300
opt_select(nrml) <- 'xyzc'
opt_filter(nrml) <- "-drop_z_below 0"
opt_output_files(nrml) <- "./canopy/{ID}"
canopy <- grid_canopy(nrml, 0.15, pitfree( subcircle = 0.2))