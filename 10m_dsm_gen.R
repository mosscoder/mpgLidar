library(sp)
library(rgdal)
library(lidR)
library(tidyverse)
library(future)
library(parallel)
library(ranger)

setwd('~/mpgPostdoc/projects/bareEarth/data/')

ctg <- catalog('./sourceCloud/')

opt_select(ctg) <- 'xyz'
plan('multisession', workers = 20L)
opt_output_files(ctg) <- './10mDSM/{ID}'
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 300

grid_metrics(ctg, func = mean(Z), res = 10)





