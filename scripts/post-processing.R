
## preamble ----
require(data.table)
require(bayesplot)
require(hexbin)
require(knitr)
require(ggplot2)
require(rstan)
require(cmdstanr)
require(ggpubr)
require(ggsci)
require(scales)
require(lubridate)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    outdir = '/Users/alexb/Box Sync/Roadmap/source_attribution/mm_bgUnif_pi1DGP_sim_230224b-simulations_408truepairs_agesrcrec_exclTE16_subsample14pct-261528',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    #stanModelFile = 'mm_sigHierG_bgUnif_piTE_220421',
    #stanModelFile = 'mm_sigHierG_bgUnif_piVanilla_220408',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_220325',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_221123',
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_230111',
    #stanModelFile = 'mm_bgUnif_piGP_221027',
    stanModelFile = 'mm_bgUnif_pi1DGP_sim_230224b',
    scenario = 1,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'simulations_408truepairs_agesrcrec_exclTE16_subsample14pct'
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

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
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

sim_scenario[,FROM_AGE_INT:= round(FROM_AGE)]
sim_scenario[,TO_AGE_INT:= round(TO_AGE)]

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
  fit.target.pars <- c('y_mix[1]','logit_y_mix_0','beta_age_src[1]','beta_age_src[2]','beta_age_src[3]','beta_age_src[4]','beta_age_src[5]',
                       'beta_age_rcp[1]','beta_age_rcp[2]','beta_age_rcp[3]','beta_age_rcp[4]','beta_age_rcp[5]',"lp__")
}else if(grepl('GP',args$stanModelFile)){
  fit.target.pars <- c('logit_y_mix_0','gpscale','lscale[1]','lscale[2]',"lp__")
}
if(grepl('1DGP',args$stanModelFile)){
  if(grepl('_agesrcrec_',args$job_tag)){
    fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
  }else{
    fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
  }
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
cat("\n ----------- make trace plots: start ----------- \n")

po <- model_fit$draws(inc_warmup = TRUE,
                      variables = fit.target.pars
)
tmp <- su$variable[which.min(su$ess_bulk)]
if(grepl('Reg',args$stanModelFile)){
  fit.target.pars <- c('y_mix[1]','logit_y_mix_0','beta_age_src[1]','beta_age_src[2]','beta_age_src[3]','beta_age_src[4]','beta_age_src[5]',
                       'beta_age_rcp[1]','beta_age_rcp[2]','beta_age_rcp[3]','beta_age_rcp[4]','beta_age_rcp[5]',"lp__")
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

cat("\n ----------- make trace plots: end ----------- \n")

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

## plot mixing weights----
cat("\n ----------- plot mixing weights ----------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'y_mix'
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

if(grepl('Vanilla|Reg',args$stanModelFile)){
  p <- ggplot(po, aes(x = variable)) +
    geom_point(aes(y=M)) +
    geom_errorbar( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'time elapsed (years)', y = 'posterior mixing weight') +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-postmixweight.png'), p, w = 7, h = 7)
}else if(grepl('Reg',args$stanModelFile)){
  po[, PAIR_ID := as.integer(gsub('y_mix\\[([0-9]+)\\]','\\1',as.character(variable)))]
  po <- merge(po, sim_scenario, by = 'PAIR_ID')
  p <- ggplot(po, aes(x = variable)) +
    geom_point(aes(y=M)) +
    geom_errorbar( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'time elapsed (years)', y = 'posterior mixing weight') +
    facet_wrap(.~FROM_AGE_GP) +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-postmixweight.png'), p, w = 7, h = 7)
}else{
  po[, PAIR_ID := as.integer(gsub('y_mix\\[([0-9]+)\\]','\\1',as.character(variable)))]
  po <- merge(po, sim_scenario, by = 'PAIR_ID')
  p <- ggplot(po, aes(x = TIME_ELAPSED)) +
    geom_ribbon( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line( aes(y = M )) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'time elapsed (years)', y = 'posterior mixing weight') +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-postmixweight.png'), p, w = 7, h = 7)
}

## plot mixing weights (GP) ----
if(grepl('GP',args$stanModelFile)){
  po <- model_fit$draws(inc_warmup = FALSE,
                        format = 'draws_df',
                        variables = 'y_mix'
  )
  po <- as.data.table(po)
  setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
  po <- melt(po, id.vars = c('chain','iteration','draw'))
  po[, PAIR_ID := as.integer(gsub('y_mix\\[([0-9]+)\\]','\\1',as.character(variable)))]
  po <- merge(po, sim_scenario, by = 'PAIR_ID')
  po <- po[,
           list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                 stat = c('M','CL','IL', 'IU', 'CU')
           ),
           by = c('FROM_AGE_INT')
  ]
  po <- dcast.data.table(po, FROM_AGE_INT~stat, value.var = 'q')

  p1 <- ggplot(po, aes(x = FROM_AGE_INT)) +
    geom_ribbon( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line( aes(y = M )) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'age of source', y = 'posterior mixing weight') +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-postmixweight_FROMage.png'), p1, w = 7, h = 7)

  po <- model_fit$draws(inc_warmup = FALSE,
                        format = 'draws_df',
                        variables = 'y_mix'
  )
  po <- as.data.table(po)
  setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
  po <- melt(po, id.vars = c('chain','iteration','draw'))
  po[, PAIR_ID := as.integer(gsub('y_mix\\[([0-9]+)\\]','\\1',as.character(variable)))]
  po <- merge(po, sim_scenario, by = 'PAIR_ID')
  po <- po[,
           list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                 stat = c('M','CL','IL', 'IU', 'CU')
           ),
           by = c('TO_AGE_INT')
  ]
  po <- dcast.data.table(po, TO_AGE_INT~stat, value.var = 'q')
  p2 <- ggplot(po, aes(x = TO_AGE_INT)) +
    geom_ribbon( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line( aes(y = M )) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'age of recipient', y = 'posterior mixing weight') +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-postmixweight_TOAGE.png'), p2, w = 7, h = 7)

  p <- ggarrange(p1, p2, nrow = 2,align='hv')
  ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'mixingweight_posterior_median_ages.png'), p, w=12, h=7)
  ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'mixingweight_posterior_median_ages.pdf'), p, w=12, h=7)
}

## plot GP function ----

if(grepl('GP',args$stanModelFile)){
  cat("\n ----------- plot GP function ----------- \n")

  if(grepl('1DGP',args$stanModelFile)){
    check <- any(grepl('hsgp_f',model_fit$metadata()$stan_variables))
    check2 <- any(grepl('hsgp_f_src',model_fit$metadata()$stan_variables))
    if(check2){
      fit.target.pars = 'hsgp_f_src'
    }else if(check){
      fit.target.pars = 'hsgp_f'
    }else{
      fit.target.pars = c("f")
    }
    po <- model_fit$draws(inc_warmup = FALSE,
                          format = 'draws_df',
                          variables = fit.target.pars
    )
    po <- as.data.table(po)
    setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
    po <- melt(po, id.vars = c('chain','iteration','draw'))
    

    po[, x_index := as.integer(gsub(paste0(fit.target.pars[1],'\\[([0-9]+)\\]'),'\\1',as.character(variable)))]

    if(grepl('srcage',args$job_tag)){
      po <- merge(po, unique(subset(sim_scenario,select=c('FROM_AGE_INT','IDX_UNIQUE_FROM_AGE'))), by.x = 'x_index',by.y='IDX_UNIQUE_FROM_AGE')
      po <- po[,
               list( q = quantile((value), probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                     stat = c('M','CL','IL', 'IU', 'CU')
               ),
               by = c("FROM_AGE_INT")
      ]
      po <- dcast.data.table(po, FROM_AGE_INT~stat, value.var = 'q')
      p1 <- ggplot(subset(po)) + geom_line( aes(x = FROM_AGE_INT, y = M )) +
        geom_ribbon( aes(x = FROM_AGE_INT,ymin = CL, ymax = CU), alpha = 0.3, fill = "light pink") +
        labs(x = 'Age of Source', y = 'Posterior Median of f',
             title = "Posterior Median of HSGP Random Function by Age of Source") +
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position='bottom',
              text = element_text(size = 10),plot.title = element_text(face="bold",hjust = 0.5,size=11)
        )
      p1
      ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_FROMage.png'), p1, w=12, h=7)
    }else{
      po <- merge(po, unique(subset(sim_scenario,select=c('TO_AGE_INT','IDX_UNIQUE_TO_AGE'))), by.x = 'x_index',by.y='IDX_UNIQUE_TO_AGE')
      po <- po[,
               list( q = quantile((value), probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                     stat = c('M','CL','IL', 'IU', 'CU')
               ),
               by = c("TO_AGE_INT")
      ]
      po <- dcast.data.table(po, TO_AGE_INT~stat, value.var = 'q')
      p1 <- ggplot(subset(po)) + geom_line( aes(x = TO_AGE_INT, y = M )) +
        geom_ribbon( aes(x = TO_AGE_INT,ymin = CL, ymax = CU), alpha = 0.3, fill = "light pink") +
        labs(x = 'Age of Source', y = 'Posterior Median of f',
             title = "Posterior Median of HSGP Random Function by Age of Recipient") +
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position='bottom',
              text = element_text(size = 10),plot.title = element_text(face="bold",hjust = 0.5,size=11)
        )
      p1
      ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_TOage.png'), p1, w=12, h=7)
    }
  }else{
    fit.target.pars = c("f")
    po <- model_fit$draws(inc_warmup = FALSE,
                          format = 'draws_df',
                          variables = fit.target.pars
    )
    po <- as.data.table(po)
    setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
    po <- melt(po, id.vars = c('chain','iteration','draw'))
    
    po[, x_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\1',as.character(variable)))]
    po[, y_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\2',as.character(variable)))]
    po <- merge(po, grid, by = c('x_index','y_index'))
    setnames(po,c('x','y'),c('FROM_AGE_INT','TO_AGE_INT'))
    
    po <- po[,
             list( q = quantile((value), probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                   stat = c('M','CL','IL', 'IU', 'CU')
             ),
             by = c("FROM_AGE_INT")
    ]
    po <- dcast.data.table(po, FROM_AGE_INT~stat, value.var = 'q')
    p1 <- ggplot(subset(po)) + geom_line( aes(x = FROM_AGE_INT, y = M )) +
      geom_ribbon( aes(x = FROM_AGE_INT,ymin = CL, ymax = CU), alpha = 0.3, fill = "light pink") +
      labs(x = 'Age of Source', y = 'Posterior Median of f',
           title = "Posterior Median of HSGP Random Function by Age of Source") +
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position='bottom',
            text = element_text(size = 10),plot.title = element_text(face="bold",hjust = 0.5,size=11)
      )
    p1
    ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_FROMage.png'), p1, w=12, h=7)
    
    po <- model_fit$draws(inc_warmup = FALSE,
                          format = 'draws_df',
                          variables = fit.target.pars
    )
    po <- as.data.table(po)
    setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
    po <- melt(po, id.vars = c('chain','iteration','draw'))
    
    po[, x_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\1',as.character(variable)))]
    po[, y_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\2',as.character(variable)))]
    po <- merge(po, grid, by = c('x_index','y_index'))
    setnames(po,c('x','y'),c('FROM_AGE_INT','TO_AGE_INT'))
    
    po <- po[,
             list( q = quantile((value), probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                   stat = c('M','CL','IL', 'IU', 'CU')
             ),
             by = c("TO_AGE_INT")
    ]
    po <- dcast.data.table(po, TO_AGE_INT~stat, value.var = 'q')
    
    p2 <- ggplot(po) + geom_line( aes(x = TO_AGE_INT, y = M ))+
      geom_ribbon( aes(x = TO_AGE_INT, ymin = CL, ymax = CU), alpha = 0.3, fill = "light pink") +
      labs(x = 'Age of Recipient', y = 'Posterior Median of f',
           title = "Posterior Median of HSGP Random Function by Age of Recipient") +
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position='bottom',
            text = element_text(size = 10),plot.title = element_text(face="bold",hjust = 0.5,size=11)
      )
    p <- ggarrange(p1, p2, nrow = 2,align='hv')
    ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_ages.png'), p, w=12, h=7)
    ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_ages.pdf'), p, w=12, h=7)
    
  }
}

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

p <- ggplot(po, aes(x = PAIR_ID2, color = FROM_AGE_GP)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'pair', y = 'transmission pair probability', colour = 'source category of transmitter') +
  facet_grid(~factor(TRANSMISSION_PAIR, levels = c('No','Yes'), labels = c('non-transmission pair', 'true transmission pair')), scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-pi_ij.png'), p, w = 12, h = 7)

p <- ggplot(po, aes(x = TIME_ELAPSED, color = FROM_AGE_GP)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'time elapsed (years)', y = 'transmission pair probability', colour = 'source category of transmitter') +
  facet_grid(~factor(TRANSMISSION_PAIR, levels = c('No','Yes'), labels = c('non-transmission pair', 'true transmission pair')), scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-pi_ij_timeelapsed.png'), p, w = 12, h = 7)

## read in clock quantiles
dps <- readRDS(file = paste0(outfile.base,'-clock_quantiles.rds'))
po[, d_TSeqT:= round(TIME_ELAPSED,1)]
po <- merge(po,dps,by=c('d_TSeqT'),all.x=T)

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
tmp <- subset(sim_scenario, select = c('PAIR_ID','FROM_AGE_GP','TRANSMISSION_PAIR'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','FROM_AGE_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_AGE_GP')
]
po <- dcast.data.table(po, FROM_AGE_GP~stat, value.var = 'q')
tmp <- sim_scenario[TRANSMISSION_PAIR == 'Yes',
                    list(value = length(PAIR_ID)),
                    by = 'FROM_AGE_GP'
]
tmp <- merge(data.table(FROM_AGE_GP = unique(sim_scenario$FROM_AGE_GP)), tmp, by = 'FROM_AGE_GP', all.x = TRUE)
set(tmp, tmp[, which(is.na(value))], 'value', 0L)
po <- merge(po, subset(tmp, select = c('FROM_AGE_GP','value')), by = 'FROM_AGE_GP')
saveRDS(po,file=paste0(outfile.base,'-number_pairs','.RDS'))

## PAF all pairs ----
cat(" \n -------------------------------- \n PAF using all pairs (GQs) \n -------------------------------- \n")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'pflows_from'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, FROM_AGE_GP_2 := as.factor(gsub('pflows_from\\[([0-9]+)\\]','\\1',as.character(variable)))]
po <- merge(po, src_map, by = 'FROM_AGE_GP_2')
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_AGE_GP')
]
po <- dcast.data.table(po, FROM_AGE_GP~stat, value.var = 'q')
tmp <- sim_scenario[TRANSMISSION_PAIR == 'Yes',
                    list(value = length(PAIR_ID)),
                    by = 'FROM_AGE_GP'
]
tmp <- model_fit$draws(inc_warmup = FALSE,
                       format = 'draws_df',
                       variables = 'true_flows'
)
tmp <- as.data.table(tmp)
setnames(tmp, colnames(tmp), gsub('^\\.','',colnames(tmp)))
tmp <- melt(tmp, id.vars = c('chain','iteration','draw'))
tmp[, FROM_AGE_GP_2 := as.factor(gsub('true_flows\\[([0-9]+)\\]','\\1',as.character(variable)))]
tmp <- merge(tmp, src_map, by = 'FROM_AGE_GP_2')
tmp <- tmp[,
         list( true_p = quantile(value, probs = c(0.5) ),
               stat = c('M')
         ),
         by = c('FROM_AGE_GP')
]
po <- merge(po, subset(tmp, select = c('FROM_AGE_GP','true_p')), by = 'FROM_AGE_GP')
saveRDS(po,file=paste0(outfile.base,'-PAF_GQ_all_pairs','.RDS'))

sim_scenario[, SRC_1_5pct:= 0]
sim_scenario[GEN_DIST<0.015, SRC_1_5pct:= 1]

tmp <- sim_scenario[, list(N_thrsh=as.numeric(sum(SRC_1_5pct))), by=c('FROM_AGE_GP')]
thrsh <- tmp[, list(FROM_AGE_GP=FROM_AGE_GP,
                    pct_thrsh=N_thrsh/sum(N_thrsh))]

plot_PAF_age(po,thrsh,plot='-PAF_GQs_age_all_pairs',outfile.base)

## MAE all pairs ----
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

cat(" \n -------------------------------- \n MAE 5yr ages \n -------------------------------- \n")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub('tpair_prob_w\\[([0-9]+)\\]','\\1',as.character(variable)))]
sim_scenario[, FROM_AGE_FIVE:= cut(FROM_AGE,breaks=seq(15,70,5),include.lowest=T,right=F,
                                   labels=c('15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70'))]
tmp <- subset(sim_scenario, select = c('PAIR_ID','FROM_AGE_FIVE','TRANSMISSION_PAIR'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','FROM_AGE_FIVE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]

tmp <- sim_scenario[TRANSMISSION_PAIR == 'Yes',
                    list(value = length(PAIR_ID)),
                    by = 'FROM_AGE_FIVE'
]
tmp <- merge(data.table(FROM_AGE_FIVE = unique(sim_scenario$FROM_AGE_FIVE)), tmp, by = 'FROM_AGE_FIVE', all.x = TRUE)
set(tmp, tmp[, which(is.na(value))], 'value', 0L)
tmp[, total := sum(value)]
tmp[, true_p := value/total]
po <- merge(po, subset(tmp, select = c('FROM_AGE_FIVE','true_p')), by = 'FROM_AGE_FIVE')

po[, ae := abs(true_p - paf)]
po[, ape := abs(true_p - paf)/true_p]

po <- po[, list( mae= mean(ae)),by=c('draw')]
po <- po[,
         list( q = quantile(mae, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         )]
po <- dcast.data.table(po, .~stat, value.var = c('q'))
saveRDS(po,file=paste0(outfile.base,'-MAE_5yrage_all_pairs','.RDS'))

## compare with PAF when using fixed thresholds ----

sim_scenario[, STAGE_thrsh:= 0]
sim_scenario[GEN_DIST<0.015, STAGE_thrsh:= 1]

tmp <- sim_scenario[, list(SCENARIO=args$scenario,
                           SRC_1_thrsh=sum(STAGE_thrsh[BIN_COV=='cat1'])/sum(STAGE_thrsh),
                           SRC_2_thrsh=sum(STAGE_thrsh[BIN_COV=='cat2'])/sum(STAGE_thrsh))]
saveRDS(tmp,file=paste0(outfile.base,'-PAF_fixed_threshold','.RDS'))
