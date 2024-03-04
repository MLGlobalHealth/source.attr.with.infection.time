require(tidyverse)
require(data.table)
require(ggplot2)
require(rstan)
require(knitr)
require(bayesplot)
library(lubridate)
require(dplyr)
require(ggsci)
source('R/functions.R')

analysis <- 'analysis_220713'
results <- 'agegps_sensanalysis_210216_MSM-2010_2022'
indir_data <- '/Users/alexb/Box Sync/Roadmap'
in.dir <- file.path('~/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam',results)

args <- list(
  source_dir ='/Users/alexb/Documents/GitHub/source.attr.with.infection.time.public',
  indir = '/Users/alexb/Box Sync/Roadmap',
  mig_groups=T,
  trsm='MSM',
  sens=T # for sensitivity analysis with patristic distances of 0
)

# load pairs data ----

pairs <- readRDS(file.path(in.dir,paste0(args$trsm,"_pairs.rds")))

anon <- subset(pairs,select=c('PAIR_ID','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','GEN_DIST','TIME_ELAPSED','FROM_AGE','TO_AGE','FROM_AGE_GP','TO_AGE_GP'))

saveRDS(anon,file=file.path(in.dir,paste0(args$trsm,"_anonymised_pairs.rds")))

