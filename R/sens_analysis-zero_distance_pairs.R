

## preamble ----
require(data.table)  # data mangling
require(bayesplot)
require(hexbin)
require(knitr)
require(ggplot2)
require(rstan)  # run Stan from R
require(cmdstanr)
require(ggsci)
require(scales)
require(grid)
require(viridis)
require(loo)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    #outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_TE16_MSM-2010_2022-1665619',
    outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_220815',
    #stanModelFile = 'mm_bgUnif_piGP_221020',
    #stanModelFile = 'mm_sigHierG_bgUnif_piTE_220421',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_230111b',
    #stanModelFile = 'mm_sigHierG_bgUnif_piVanilla_220408b',
    #stanModelFile = 'mm_sigHierG_bgUnif_pi1DGP_230224',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    #job_tag = 'agegps_updated_criteria_MSM-2010_2022',
    #job_tag = 'agegps_TE16_MSM-2010_2022',
    job_tag = 'agegps_sensanalysis_210216_MSM',
    weights = T,
    tpair <- 'tpair_prob_w'
    #tpair <- 'tpair_prob'
  )
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  stopifnot(args_line[[5]]=='-indir')
  stopifnot(args_line[[7]]=='-outdir')
  stopifnot(args_line[[9]]=='-job_tag')
  stopifnot(args_line[[11]]=='-scenario')
  stopifnot(args_line[[13]]=='-rep')
  stopifnot(args_line[[15]]=='-weights')
  
  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['indir']] <- args_line[[6]]
  args[['outdir']] <- args_line[[8]]
  args[['job_tag']] <- args_line[[10]]
  args[['scenario']] <- as.integer(args_line[[12]])
  args[['rep']] <- as.integer(args_line[[14]])
  args[['weights']] <- as.integer(args_line[[16]])
}
args
replicate <- args$rep
cat(paste0("rep ", replicate))
## load functions
source(file.path(args$source_dir, 'R', 'functions_simulation_scenarios_bplace.R'))
source(file.path(args$source_dir, 'R', 'functions_fit_model.R'))

cat(" \n --------------------------------  with arguments -------------------------------- \n")

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('rep_',replicate,'_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

args$tpair <- 'tpair_prob_w'

do <- do[order(PAIR_ID),]

tmp <- paste0(outfile.base,'-rep_',replicate,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

cat(" \n -------------------------------- \n Transmission pair probabilities \n -------------------------------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w','\\[([0-9]+)\\]'),'\\1',as.character(variable)))]

tmp <- subset(do, select = c('PAIR_ID','FROM_AGE','TO_AGE','TO_AGE_GP','FROM_AGE_GP'))

po <- merge(po, tmp, by = 'PAIR_ID')

po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('PAIR_ID','FROM_AGE_GP')
]
po <- dcast.data.table(po, PAIR_ID+FROM_AGE_GP~stat, value.var = 'q')
tmp <- subset(po, PAIR_ID %in% c(2174,286,850,810,2074,621,842))
write.csv(tmp,file=paste0(outfile.base,'-pairprob_pairs_gendist_zero.csv'))

tmp[, PAIR_ID_fct:= as.factor(PAIR_ID)]
ggplot(tmp, aes(x = PAIR_ID_fct, color = FROM_AGE_GP)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  ggsci::scale_color_npg() +
  labs(x = 'pair', y = 'transmission pair probability', colour = 'infector stage') +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-pi_ij_zerodistancepairs_age.png'), p, w = 12, h = 7)

