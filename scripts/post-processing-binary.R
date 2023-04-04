
## preamble ----
require(data.table)
require(bayesplot)
require(hexbin)
require(knitr)
require(ggplot2)
require(rstan)
require(cmdstanr)
require(ggsci)
require(scales)
require(lubridate)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    #outdir = '~/Box\ Sync/Roadmap/source_attribution/mm_sigHierG_bgUnif_piTE_220421-simulations_mm_sigHierG_bgUnif_piTE_220421_networks-scenario_1-1203896',
    outdir = '/Users/alexb/Box Sync/Roadmap/source_attribution/mm_sigHierG_bgUnif_piReg_230118-simulations_500truepairs_perfectpred_exclTE16_subsample0pct-1155684',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    #stanModelFile = 'mm_sigHierG_bgUnif_piTE_220421',
    #stanModelFile = 'mm_sigHierG_bgUnif_piVanilla_220408',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_220325',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_221123',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_230111',
    stanModelFile = 'mm_sigHierG_bgUnif_piReg_230118',
    #stanModelFile = 'mm_bgUnif_piGP_221027_KD',
    scenario = 1,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'simulations_500truepairs_perfectpred_exclTE16_subsample0pct'
  )
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  #stopifnot(args_line[[7]]=='-chain', !is.na(as.integer(args_line[[8]])))
  stopifnot(args_line[[5]]=='-indir')
  stopifnot(args_line[[7]]=='-outdir')
  stopifnot(args_line[[9]]=='-job_tag')

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  #args[['chain']] <- as.integer(args_line[[8]])
  args[['indir']] <- args_line[[6]]
  args[['outdir']] <- args_line[[8]]
  args[['job_tag']] <- args_line[[10]]
}
args
replicate <- args$rep
cat(paste0("rep ", replicate))
## load functions
source(file.path(args$source_dir, 'R', 'functions_simulation_scenarios.R'))

cat(" \n --------------------------------  with arguments -------------------------------- \n")
#str(args_dir)

#do <- data.table(F=list.files(args_dir$out_dir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
#cat(paste("\n", nrow(do),"/",args_dir$numb_chains, "chains are finished \n"))

#outfile.base <- unique( do[, file.path(dirname(dirname(F)), paste0(args$stanModelFile,'-',args$job_tag))] )
#stopifnot(length(outfile.base)==1 )

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

## diagnostics ----
cat(" \n -------------------------------- \n Check convergence and divergences \n -------------------------------- \n")

# Time of execution
time <- model_fit$time()
run_time <- seconds_to_period(time$total)
write.csv(run_time,file=paste0(outfile.base,'-run-time.csv'))

if(grepl('Vanilla',args$stanModelFile)){
  fit.target.pars <- c('logit_y_mix_0','y_mix','log_alpha1_pair[1]','log_phi_pair[1]',"lp__")
}else if(grepl('Reg',args$stanModelFile)){
  fit.target.pars <- c('y_mix[1]','logit_y_mix_0','beta_age',"lp__")
  fit.target.pars <- c('y_mix[1]','logit_y_mix_0','beta_src[1]','beta_src[2]',"lp__")
}else if(grepl('GP',args$stanModelFile)){
  fit.target.pars <- c('logit_y_mix_0','gpscale','lscale[1]','lscale[2]',"lp__")
}
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = fit.target.pars
)
su <- as.data.table(posterior::summarise_draws(po))
tmp <- paste0(outfile.base,"-convergence.csv")
write.csv(su, file = tmp, row.names = TRUE)
su[,min(ess_bulk)]
su[,max(rhat)]

## traces----
po <- model_fit$draws(inc_warmup = TRUE,
                      variables = fit.target.pars
)
tmp <- su$variable[which.min(su$ess_bulk)]
if(grepl('Reg',args$stanModelFile)){
  fit.target.pars <- c('y_mix[1]','logit_y_mix_0','beta_src[1]','beta_src[2]',"lp__")
}
tmp <- unique(fit.target.pars,tmp)
bayesplot:::color_scheme_set("mix-blue-pink")
p <- bayesplot:::mcmc_trace(po,
                            pars = tmp,
                            n_warmup = 5e2,
                            facet_args = list(ncol = 1, labeller = label_parsed)
)
ggsave(file = paste0(outfile.base,'-traces.pdf'), p, w = 12, h = 20)

tmp <- su$variable[which.min(su$ess_bulk)]
po <- model_fit$draws(inc_warmup = TRUE,
                      variables = tmp
)
bayesplot:::color_scheme_set("mix-blue-pink")
min <- quantile(po,0.01)
max <- max(po)
range <- max - min
p <- bayesplot:::mcmc_trace(po,
                            pars = tmp,
                            n_warmup = 5e2,
                            facet_args = list(ncol = 1, labeller = label_parsed)
)+
  bayesplot:::legend_move("bottom")+ggplot2::labs(x="Iteration")+ggplot2::theme(text = element_text(size = 16))+
  xaxis_text(size = 16)+facet_text(size = 16)+legend_text(size = 16)+yaxis_text(size = 16)+xaxis_ticks(size = 14)+
  yaxis_ticks(size = 14)+coord_cartesian(ylim=c(min-range,max+range))
ggsave(file = paste0(outfile.base,'-trace_lwstneff.pdf'), p, w = 10, h = 5)

#
# Pairs plots ----
#
cat("\n ----------- make pairs plots: start ----------- \n")
pd <- model_fit$draws(inc_warmup = FALSE,
                      variables = c(fit.target.pars))
bayesplot:::color_scheme_set("mix-blue-pink")
p <- bayesplot:::mcmc_pairs(pd,
                            pars = c(fit.target.pars),
                            diag_fun = "dens",
                            off_diag_fun = "hex"
)
ggsave(p, file = paste0(outfile.base, "-HMC-pairs_transmission_pars.pdf"), w=length(fit.target.pars)*2, h=length(fit.target.pars)*2)
cat("\n ----------- make pairs plots: end ----------- \n")

## plot posterior probabilities of being a pair----
cat(" \n -------------------------------- \n Plot posterior probabilities of being a pair \n -------------------------------- \n")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('variable')
]
po <- dcast.data.table(po, variable~stat, value.var = 'q')
po[, PAIR_ID := as.integer(gsub('tpair_prob_w\\[([0-9]+)\\]','\\1',as.character(variable)))]
po <- merge(po, sim_scenario, by = 'PAIR_ID')
tmp <- po[, .(PAIR_ID = PAIR_ID, PAIR_ID2 = seq_along(PAIR_ID)), by = 'TRANSMISSION_PAIR']
po <- merge(po, tmp, by = c('TRANSMISSION_PAIR','PAIR_ID'))

po[, SOURCE:= 'Source category 1']
po[BIN_COV=='cat2', SOURCE:= 'Source category 2']

p <- ggplot(po, aes(x = PAIR_ID2, color = SOURCE)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'pair', y = 'transmission pair probability', colour = 'source category of transmitter') +
  facet_grid(~factor(TRANSMISSION_PAIR, levels = c('No','Yes'), labels = c('non-transmission pair', 'true transmission pair')), scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-pi_ij.png'), p, w = 12, h = 7)

p <- ggplot(po, aes(x = TIME_ELAPSED, color = SOURCE)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'time elapsed (years)', y = 'transmission pair probability', colour = 'source category of transmitter') +
  facet_grid(~factor(TRANSMISSION_PAIR, levels = c('No','Yes'), labels = c('non-transmission pair', 'true transmission pair')), scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-pi_ij_timeelapsed.png'), p, w = 12, h = 7)

po[SOURCE_PAIR=='cat1 - No', SOURCE_PAIR:='Source category 1 - No' ]
po[SOURCE_PAIR=='cat2 - Yes', SOURCE_PAIR:='Source category 2 - Yes' ]

cat( '\n plot simulated data coloured by posterior probabilities trsm pair \n ' )
make_plot_simulated_data_colour_prob_tpair(args$rep,po,dps,sim_scenario,outfile.base)

## violin plot posterior probabilities of being a pair----
cat(" \n -------------------------------- \n Violin plot of posterior median of being a pair \n -------------------------------- \n")

p <- ggplot(po, aes(x = TRANSMISSION_PAIR, y = M)) +
  geom_jitter(aes(color = SOURCE), width = 0.3, height = 0, alpha = 0.7) +
  geom_violin(fill = 'transparent') +
  scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1)) +
  ggsci::scale_color_npg() +
  labs(x = 'true transmission pair', y = 'transmission pair probability\n(posterior median)', colour = 'source category of transmitter') +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-pi_ij_violin.png'), p, w = 9, h = 7)

## Number of infections ----
cat(" \n -------------------------------- \n Number of infections \n -------------------------------- \n")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub('tpair_prob_w\\[([0-9]+)\\]','\\1',as.character(variable)))]
tmp <- subset(sim_scenario, select = c('PAIR_ID','BIN_COV','TRANSMISSION_PAIR'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','BIN_COV')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('BIN_COV')
]
po <- dcast.data.table(po, BIN_COV~stat, value.var = 'q')
tmp <- sim_scenario[TRANSMISSION_PAIR == 'Yes',
                    list(value = length(PAIR_ID)),
                    by = 'BIN_COV'
]
tmp <- merge(data.table(BIN_COV = unique(sim_scenario$BIN_COV)), tmp, by = 'BIN_COV', all.x = TRUE)
set(tmp, tmp[, which(is.na(value))], 'value', 0L)
po <- merge(po, subset(tmp, select = c('BIN_COV','value')), by = 'BIN_COV')
saveRDS(po,file=paste0(outfile.base,'-number_pairs','.RDS'))

## PAF ----
cat(" \n -------------------------------- \n PAF using all pairs \n -------------------------------- \n")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'pflows_from'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, TRANS_STAGE_2 := as.factor(gsub('pflows_from\\[([0-9]+)\\]','\\1',as.character(variable)))]
po <- merge(po, src_map, by = 'TRANS_STAGE_2')
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('TRANS_STAGE')
]
po <- dcast.data.table(po, TRANS_STAGE~stat, value.var = 'q')
tmp <- sim_scenario[TRANSMISSION_PAIR == 'Yes',
                    list(value = length(PAIR_ID)),
                    by = 'TRANS_STAGE'
]
tmp <- model_fit$draws(inc_warmup = FALSE,
                       format = 'draws_df',
                       variables = 'true_flows'
)
tmp <- as.data.table(tmp)
setnames(tmp, colnames(tmp), gsub('^\\.','',colnames(tmp)))
tmp <- melt(tmp, id.vars = c('chain','iteration','draw'))
tmp[, TRANS_STAGE_2 := as.factor(gsub('true_flows\\[([0-9]+)\\]','\\1',as.character(variable)))]
tmp <- merge(tmp, src_map, by = 'TRANS_STAGE_2')
tmp <- tmp[,
         list( true_p = quantile(value, probs = c(0.5) ),
               stat = c('M')
         ),
         by = c('TRANS_STAGE')
]
po <- merge(po, subset(tmp, select = c('TRANS_STAGE','true_p')), by = 'TRANS_STAGE')
saveRDS(po,file=paste0(outfile.base,'-PAF_GQ_all_pairs','.RDS'))

po[, BIN_COV:= factor(TRANS_STAGE,levels=c('Undiagnosed','Diagnosed'),labels=c('Source category 1','Source category 2'))]
plot_estimated_attributable_fraction(po,plot='-PAF_GQ_all_pairs',outfile.base)

## MAE ----
cat(" \n -------------------------------- \n PAF using all pairs \n -------------------------------- \n")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'MAE_from'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))

po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         )]
po <- dcast.data.table(po, .~stat, value.var = c('q'))
saveRDS(po,file=paste0(outfile.base,'-MAE_GQ_all_pairs','.RDS'))

## compare with PAF when using fixed thresholds ----

sim_scenario[, STAGE_thrsh:= 0]
sim_scenario[GEN_DIST<0.015, STAGE_thrsh:= 1]

tmp <- sim_scenario[, list(SCENARIO=args$scenario,
                           SRC_1_thrsh=sum(STAGE_thrsh[BIN_COV=='cat1'])/sum(STAGE_thrsh),
                           SRC_2_thrsh=sum(STAGE_thrsh[BIN_COV=='cat2'])/sum(STAGE_thrsh))]
saveRDS(tmp,file=paste0(outfile.base,'-PAF_fixed_threshold','.RDS'))
