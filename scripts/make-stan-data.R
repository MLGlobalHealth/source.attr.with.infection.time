cat(" \n ------------- \n \n Running make-stan-data.R \n \n ------------ \n")
## preamble ----
require(data.table)
require(igraph)
require(fitdistrplus)
require(EnvStats)
require(ggplot2)
require(ggsci)
require(boot)
require(abind)
require(lubridate)
require(rstan)
require(cmdstanr)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time.public',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    outdir = '~/Box\ Sync/Roadmap/source_attribution',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_sigHierG_bgUnif_piVanilla_220408', # vanilla model
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_230111', # covariate model
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_230118', # binary predictor
    #stanModelFile = 'mm_bgUnif_pi1DGP_sim_221027', # 2D HSGP model
    #stanModelFile = 'mm_bgUnif_pi1DGP_sim_230224', # 1D HSGP model
    #stanModelFile = 'mm_bgUnif_pi1DGP_sim_230224b', # 2 * 1D HSGP model
    hmc_stepsize = 0.02,
    hmc_num_samples = 15,
    hmc_num_warmup = 10,
    seed = 42,
    chain = 1,
    scenario = 1,
    reps = 1,
    ncases = 500,
    #ncases = 450,    # to simulate Amsterdam data configuration
    p_pairs = 0.5,  # 0.5 - 0.01
    #p_pairs = 0.14,  # to simulate Amsterdam data configuration
    simulate_data = T,
    job_tag = 'simulations_500truepairs_srcage_exclTE16_subsample50pct',
    local=1,
    networks=T,
    sim_binary='fixed', # or after
    bg = 'unif',
    no_corr=0,
    excl_improb=T,
    all_pairs=F,
    late=F,
    network_data_file='V1.2_patch0_Rand10_Run1_PCseed0_0_transmission_events.csv',
    m1 = 24,
    m2 = 24,
    B=1.2,
    src_cat='age' #'binary'/'age'
  )
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  stopifnot(args_line[[5]]=='-seed', !is.na(as.integer(args_line[[6]])))
  stopifnot(args_line[[7]]=='-chain', !is.na(as.integer(args_line[[8]])))
  stopifnot(args_line[[9]]=='-indir')
  stopifnot(args_line[[11]]=='-outdir')
  stopifnot(args_line[[13]]=='-job_tag')
  stopifnot(args_line[[15]]=='-clock_model')
  stopifnot(args_line[[17]]=='-local')
  stopifnot(args_line[[19]]=='-ncases')
  stopifnot(args_line[[21]]=='-p_pairs')
  stopifnot(args_line[[23]]=='-network_data_file')
  stopifnot(args_line[[25]]=='-m1')
  stopifnot(args_line[[27]]=='-m2')
  stopifnot(args_line[[29]]=='-B')
  stopifnot(args_line[[31]]=='-src_cat')

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['seed']] <- as.integer(args_line[[6]])
  args[['chain']] <- as.integer(args_line[[8]])
  args[['indir']] <- args_line[[10]]
  args[['outdir']] <- args_line[[12]]
  args[['job_tag']] <- args_line[[14]]
  args[['clock_model']] <- args_line[[16]]
  args[['local']] <- args_line[[18]]
  args[['ncases']] <- as.integer(args_line[[20]])
  args[['p_pairs']] <- as.numeric(args_line[[22]])
  args[['network_data_file']] <- args_line[[24]]
  args[['m1']] <- args_line[[26]]
  args[['m2']] <- args_line[[28]]
  args[['B']] <- args_line[[30]]
  args[['src_cat']] <- args_line[[32]]
}

cat(" \n ------------- \n setup \n ------------ \n")

args

## load functions
source(file.path(args$source_dir, 'R', 'functions_simulation_scenarios.R'))


## set other args
args$file_stanModel <- file.path(args$source_dir, 'stan_model_files',paste0(args$stanModelFile,'.stan'))
tmp <- Sys.getenv("PBS_JOBID")
args$job_id <- ifelse(tmp != '', tmp, as.character(abs(round(rnorm(1) * 1e6))) )
args$job_dir <- args$outdir
if(args$local==1) args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag,'-',args$job_id))
args$savedata <- TRUE
if (grepl("\\[1\\]", args$job_id))
  args$savedata <- TRUE

out.dir <- args$job_dir
if(args$local==1) dir.create( out.dir )
outfile.base <- file.path(out.dir, paste0(args$stanModelFile,"-",args$job_tag))

## load Belgian transmission chain data ----
cat(" \n ------------- \n Load Belgian transmission chain data \n ------------ \n")
file <- file.path(args$source_dir,'data_other','140921_set7_INFO_TRM.R')
load(file)

# exclude B->A, same as A->B
trm.pol.nA <- subset(trm.pol, withA == FALSE, select = c(d_SeqT, d_TSeqT, BRL, FROM, TO))
setkey(trm.pol.nA, d_TSeqT)
trm.pol.nA[, id := seq_len(nrow(trm.pol.nA))]
trm.pol.nA[, pair_id := as.integer(factor(paste0(FROM,'->',TO)))]



## load estimated molecular clock ----
cat(" \n ------------- \n Load quantiles from fitted molecular clock model \n ------------ \n")

cm <- readRDS(file.path(args$clock_model,'clock_model_gamma_hier_220315-stan_fit.rds'))
pd <- cm$draws(inc_warmup = FALSE)
po <- list()
tmp <- pd[,,which(grepl('log_alpha1$',dimnames(pd)[[3]]))]
po$log_alpha1 <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_alpha1_pair_sd',dimnames(pd)[[3]]))]
po$log_alpha1_pair_sd <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_phi$',dimnames(pd)[[3]]))]
po$log_phi <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_phi_pair_sd',dimnames(pd)[[3]]))]
po$log_phi_pair_sd <- as.vector(tmp)
# get quantiles of mean/median distances at prediction points
dpr <- data.table(d_TSeqT = seq(0.1,18,0.1))
dpr[, id := seq_len(nrow(dpr))]
dpo <- data.table(iter = seq_along(po$log_alpha1),
                  log_alpha1 = po$log_alpha1,
                  log_alpha1_pair_sd = po$log_alpha1_pair_sd,
                  log_phi = po$log_phi,
                  log_phi_pair_sd = po$log_phi_pair_sd,
                  DUMMY = 1
)
dpo[, log_alpha1_pair := rnorm(nrow(dpo), 0, log_alpha1_pair_sd)]
dpo[, log_phi_pair := rnorm(nrow(dpo), 0, log_phi_pair_sd)]
dpo[, er := exp(log_alpha1 + log_alpha1_pair)]
dpo[, beta := exp( -(log_phi + log_phi_pair))]

dpr[, DUMMY := 1]
dpo <- merge(dpo, dpr, by = 'DUMMY', allow.cartesian = TRUE)
set(dpr, NULL, 'DUMMY', NULL)
dpo[, y_pr := rgamma(nrow(dpo), shape = er * d_TSeqT * beta, rate = beta)]

dps <- dpo[,
           list(
             p = quantile(y_pr, prob = c(0.025, seq(0.1, .9, 0.1), 0.975)),
             qlabel = paste0('q',c(0.025, seq(0.1, .9, 0.1), 0.975)*100)
           ),
           by = c('d_TSeqT')
]
dps <- dcast.data.table(dps, d_TSeqT ~ qlabel, value.var = 'p')

p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
p.alpha <- 0.7
p <- ggplot(dps, aes(x = d_TSeqT)) +
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_line(data = dps, aes(y = q50)) +
  geom_point(data = trm.pol.nA, aes(y = BRL), size = 1.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0),labels = scales::label_percent(accuracy = 1L)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,13) ) +
  theme_bw() +
  labs(x = 'Time elapsed\n(years)', y = 'Genetic distance\n(subst/site)') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4),
        panel.grid.minor = element_line(colour = "grey70", size = 0.2)
  )
ggsave(file = paste0(outfile.base,'-pr_clock_metapop.png'), p, w = 6, h = 6)

## predict distances using posterior medians only
dpo <- data.table(iter = seq_along(po$log_alpha1),
                  log_alpha1 = median(po$log_alpha1),
                  log_alpha1_pair_sd = median(po$log_alpha1_pair_sd),
                  log_phi = median(po$log_phi),
                  log_phi_pair_sd = median(po$log_phi_pair_sd),
                  DUMMY = 1
)
dpo[, log_alpha1_pair := rnorm(nrow(dpo), 0, log_alpha1_pair_sd)]
dpo[, log_phi_pair := rnorm(nrow(dpo), 0, log_phi_pair_sd)]
dpo[, er := exp(log_alpha1 + log_alpha1_pair)]
dpo[, beta := exp( -(log_phi + log_phi_pair))]

dpr <- data.table(d_TSeqT = seq(0.01,40,0.01))
dpr[, id := seq_len(nrow(dpr))]
dpr[, DUMMY := 1]
dpo <- merge(dpo, dpr, by = 'DUMMY', allow.cartesian = TRUE)
set(dpr, NULL, 'DUMMY', NULL)
dpo[, y_pr := rgamma(nrow(dpo), shape = er * d_TSeqT * beta, rate = beta)]

dps <- dpo[,
           list(
             p = quantile(y_pr, prob = c(0.025, seq(0.1, .9, 0.1), 0.975, 0.25, 0.75)),
             qlabel = paste0('q',c(0.025, seq(0.1, .9, 0.1), 0.975, 0.25, 0.75)*100)
           ),
           by = c('d_TSeqT')
]
dps <- dcast.data.table(dps, d_TSeqT ~ qlabel, value.var = 'p')
saveRDS(dps,file = paste0(outfile.base,'-clock_quantiles.rds'))

p <- ggplot(dps, aes(x = d_TSeqT)) +
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_line(data = dps, aes(y = q50)) +
  geom_point(data = trm.pol.nA, aes(y = BRL), size = 1.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0),labels = scales::label_percent(accuracy = 1L)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,13) ) +
  theme_bw() +
  labs(x = 'Time elapsed\n(years)', y = 'Genetic distance\n(subst/site)') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4),
        panel.grid.minor = element_line(colour = "grey70", size = 0.2)
  )
ggsave(file = paste0(outfile.base,'-pr_clock_metapop_posterior_median.png'), p, w = 6, h = 6)

# specify median posterior of parameters to simulate with (from convergence.csv for gamma clock model)
log_alpha1 <- median(po$log_alpha1) 
log_alpha1_pair_sd <- median(po$log_alpha1_pair_sd)
log_phi <- median(po$log_phi) 
log_phi_pair_sd <- median(po$log_phi_pair_sd)

## simulate data ----
cat(" \n ------------- \n simulate data \n ------------ \n")
args$sim_binary='fixed'
args$excl_improb=T
args$all_pairs=F
networks <- file.path(args$source_dir,'data_other','transmission_network',args$network_data_file)
infile.inftime <- file.path(args$source_dir,'data_Ams','analysis_220713','roadmap_cd4_v3_est.csv')
simulate_scenarios_networks(networks,log_alpha1, log_alpha1_pair_sd, log_phi, log_phi_pair_sd, ncases=args$ncases, p_pairs=args$p_pairs, all_pairs=args$all_pairs, dps,
                            sim_binary=args$sim_binary, excl_improb=args$excl_improb, infile.inftime, outfile.base)

# load data
sim_scenario <- readRDS(file = paste0(outfile.base,'.rds'))

# define groupings for model ----
sim_scenario[, FROM_AGE_GP:= cut(FROM_AGE,breaks=c(15,30,40,50,60,100),include.lowest=T,right=F,
                                 labels=c('[15-30)','[30-40)','[40-50)','[50-60)','[60+)'))]
sim_scenario[, TO_AGE_GP:= cut(TO_AGE,breaks=c(15,30,40,50,60,100),include.lowest=T,right=F,
                               labels=c('[15-30)','[30-40)','[40-50)','[50-60)','[60+)'))]

sim_scenario[,FROM_AGE_GP_2:= factor(FROM_AGE_GP,labels=c(1,2,3,4,5))]
sim_scenario[,TO_AGE_GP_2:= factor(TO_AGE_GP,labels=c(1,2,3,4,5))]

if(grepl('GP',args$stanModelFile)){
  # modify ages
  # max age 90
  sim_scenario[TO_AGE>90, TO_AGE:=90]
  sim_scenario[FROM_AGE>90, FROM_AGE:=90]

  sim_scenario[,FROM_AGE_INT:= round(FROM_AGE)]
  sim_scenario[,TO_AGE_INT:= round(TO_AGE)]

  # get unique pairs ages
  tmp <- data.table(unique(cbind(sim_scenario$FROM_AGE_INT,sim_scenario$TO_AGE_INT)))
  tmp[, IDX_UNIQUE_PAIR := seq_len(nrow(tmp))]
  setnames(tmp, "V1", "FROM_AGE_INT")
  setnames(tmp, "V2", "TO_AGE_INT")
  sim_scenario <- merge(sim_scenario, tmp, by = c("FROM_AGE_INT","TO_AGE_INT"))

  # re-order by ages
  sim_scenario <- sim_scenario[order(FROM_AGE_INT,TO_AGE_INT),]
  sim_scenario <- data.table(sim_scenario)
  sim_scenario[, PAIR_ID := seq(1,nrow(sim_scenario)) ]

  }

## generate stan data ----
cat(" \n -------------------------------- \n generate stan data \n -------------------------------- \n")

  stan_data <- list()
  stan_data$N = nrow(sim_scenario)
  stan_data$P <- length(unique(sim_scenario$TO_ID))
  if(args$src_cat=='age'){
    stan_data$S <- length(unique(sim_scenario$FROM_AGE_GP_2))
  }
  stan_data$y = sim_scenario$GEN_DIST
  stan_data$x = sim_scenario$TIME_ELAPSED
  tmp <- unique(sim_scenario, by = 'IDX_UNIQUE_TE')[order(IDX_UNIQUE_TE)]
  stan_data$Nux = nrow(tmp)
  stan_data$ux = tmp$TIME_ELAPSED
  stan_data$x_2_ux = sim_scenario$IDX_UNIQUE_TE
  stan_data$idx_true_pairs = as.numeric(factor(sim_scenario$TRANSMISSION_PAIR))-1
  # posterior medians from clock model
  stan_data$log_alpha1 <- log_alpha1
  stan_data$log_alpha1_pair_sd <- log_alpha1_pair_sd
  stan_data$log_phi <- log_phi
  stan_data$log_phi_pair_sd <- log_phi_pair_sd


  stan_data$obs_to_src_idx <- sim_scenario$TRANS_STAGE # for binary predictor
  stan_data$obs_to_age_src_idx <- sim_scenario$FROM_AGE_GP_2
  stan_data$obs_to_age_rcp_idx <- sim_scenario$TO_AGE_GP_2
  stan_data$A <- length(unique(sim_scenario$FROM_AGE_GP)) # coarse age groups
  
  if(grepl('perfectpred',args$job_tag)){
    stan_data$A <- length(unique(sim_scenario$TRANS_STAGE))
  }
  if(grepl('GP',args$stanModelFile)){

    age_vec <- seq(1, 90)
    stan_data$a <- max(age_vec)

    tmp <- sim_scenario[, list(PAIR_ID), by = c("FROM_AGE_INT")]
    tmp2 <- tmp[, list(num = length(PAIR_ID)), by = c("FROM_AGE_INT")]
    setkey(tmp2, "FROM_AGE_INT")
    tmp <- data.table(age = rep(0,90))
    tmp[which(!(age_vec %in% tmp2$FROM_AGE_INT)), num:=0]
    tmp[which(age_vec %in% tmp2$FROM_AGE_INT), num:=tmp2$num]
    stan_data$units_to_obs_length <- tmp$num #Units_to_obs_length
    stan_data$Nk_max <- max(stan_data$units_to_obs_length)

    # age
    tmp <- sim_scenario[, list(PAIR_ID), by = c("FROM_AGE_INT")]
    tmp2 <- tmp[, list(num = length(PAIR_ID)), by = c("FROM_AGE_INT")]
    setkey(tmp2, "FROM_AGE_INT")
    tmp <- data.table(age = seq(1,90,1))
    tmp[which(!(age_vec %in% tmp2$FROM_AGE_INT)), num:=0]
    tmp[which(age_vec %in% tmp2$FROM_AGE_INT), num:=tmp2$num]
    stan_data$units_to_obs_length <- tmp$num #Units_to_obs_length
    stan_data$Nk_max <- max(stan_data$units_to_obs_length)

    # create a matrix of size #ages X max(count age a) and fill in pair IDs for each age of source
    tmp <- matrix(0, nrow = stan_data$a, ncol = stan_data$Nk_max )
    for (a in 1:stan_data$a)
    {tmp3 <- subset(sim_scenario, FROM_AGE_INT == a)
    if(nrow(tmp3) != 0)
    {tmp[a, 1:nrow(tmp3)] <- tmp3[,PAIR_ID]}
    }
    stan_data$age_to_obs_idx <- tmp

    stan_data$A <- length(unique(sim_scenario$FROM_AGE_GP)) # coarse age groups
    stan_data$obs_to_age_idx <- sim_scenario$FROM_AGE_GP

    # create an array with the pair IDs for sources of age a and recipients of age j
    mat <- list()
    for (a in 1:stan_data$a){
      tmp <- subset(sim_scenario, FROM_AGE_INT == a)
      mat[[a]] <- list()
      if(nrow(tmp) != 0){
        for (j in 1:stan_data$a){
          tmp2 <- subset(sim_scenario, FROM_AGE_INT == a & TO_AGE_INT == j)
          if(nrow(tmp2) != 0)
            mat[[a]][[j]] <- tmp2[,PAIR_ID]
          else{
            mat[[a]][[j]] <- 0
          }
        }}
    }

    # for each age of sources, get max # pairs across recipient ages
    max_length = rep(0, a)
    for(a in 1:stan_data$a){
      max_length[a] = max(sapply(seq(1,90), function(x) length(unlist(mat[[a]][x]))))
    }
    row_l <- max(max_length) # max # pairs across source/recipient age matrix
    empty <- array(NA,c(90, row_l, 0))
    ls <- list()
    for(a in 1:stan_data$a){
      tmp <- matrix(0, nrow = row_l, ncol = stan_data$a)
      for(j in 1: stan_data$a){
        ind_length <- length(which(sim_scenario$FROM_AGE_INT == a & sim_scenario$TO_AGE_INT == j)) # number of pairs in entry a,j
        if(ind_length>0){
          tmp[1:ind_length, j] <- unlist(mat[[a]][j])
        }
      }
      ls[[a]] <- tmp
      empty[a,,] <- ls[[a]]
    }
    newarray <- abind( ls, along = 3)

    stan_data$age_recip <- newarray
    stan_data$Nj_max <- row_l

    #Coarse age band index matrix
    mat <- list()
    for (a in 1:stan_data$A){
      tmp <- subset(sim_scenario, FROM_AGE_GP_2 == a)
      mat[[a]] <- list()
      if(nrow(tmp) != 0){
        for (j in 1:stan_data$A){
          tmp2 <- subset(sim_scenario, FROM_AGE_GP_2 == a & TO_AGE_GP_2 == j)
          if(nrow(tmp2) != 0)
            mat[[a]][[j]] <- tmp2[,PAIR_ID]
          else{
            mat[[a]][[j]] <- 0
          }
        }}
    }

    max_length <- rep(0,stan_data$A)
    for(a in 1:stan_data$A){
      max_length[a] = max(sapply(seq(1,stan_data$A), function(x) length(unlist(mat[[a]][x]))))
    }
    row_l <- max(max_length)
    empty <- array(NA,c(stan_data$A, row_l, 0))
    ls <- list()
    for(a in 1:stan_data$A){
      tmp <- matrix(0, nrow = row_l, ncol = stan_data$A)
      for(j in 1: stan_data$A){
        ind_length <- length(which(sim_scenario$FROM_AGE_GP_2 == a & sim_scenario$TO_AGE_GP_2 == j))
        if(ind_length>0){
          tmp[1:ind_length, j] <- unlist(mat[[a]][j])
        }
      }
      ls[[a]] <- tmp
      empty[a,,] <- ls[[a]]
    }
    newarray <- abind( ls, along = 3)

    stan_data$age_recip_coarse <- newarray
    stan_data$Nj_max_coarse <- row_l

    indices <- matrix(NA, args$m1*args$m2, 2)
    mm = 0
    for (m1 in 1:args$m1){
      for (m2 in 1:args$m2){
        mm = mm +1
        indices[mm,] = c(m1, m2)
      }
    }
    stan_data$indices <- indices
    stan_data$D <- 2
    stan_data$L <- rep(args$B*max(sim_scenario$FROM_AGE_INT,sim_scenario$TO_AGE_INT),2) # make grid symmetric for ages of sources/recipients
    stan_data$M <- c(args$m1,args$m2)
    stan_data$M_nD <- stan_data$M[1] * stan_data$M[2]
    tmp <- unique(sim_scenario, by = 'IDX_UNIQUE_TE')[order(FROM_AGE)]
    stan_data$sd1 <- sd(tmp$FROM_AGE)

    tmp <- unique(sim_scenario, by = 'IDX_UNIQUE_TE')[order(TO_AGE)]
    stan_data$sd2 <- sd(tmp$TO_AGE)

    stan_data$Nu_pairs <- max(sim_scenario$IDX_UNIQUE_PAIR)
    stan_data$IDX_UNIQUE_PAIR <- sim_scenario$IDX_UNIQUE_PAIR

    stan_data$n <- 90
    stan_data$m <- 90

    A <- stan_data$n

    age_idx <- seq.int(1,A,1)

    age_idx_std <- (age_idx - mean(age_idx))/sd(age_idx)

    stan_data$ages <- matrix(data=NA,stan_data$n,2)
    stan_data$ages[,1] <- age_idx_std
    stan_data$ages[,2] <- age_idx_std

    ages1_grid = data.table(x = seq.int(1,A,1))
    ages1_grid[, x_index := 1:nrow(ages1_grid)]
    ages2_grid = data.table(y = seq.int(1,A,1))
    ages2_grid[, y_index := 1:nrow(ages2_grid)]

    # find all coordinates
    grid = as.data.table( expand.grid(x_index = ages1_grid$x_index,
                                      y_index = ages2_grid$y_index) )
    grid = merge(grid, ages1_grid, by = 'x_index')
    grid = merge(grid, ages2_grid, by = 'y_index')
    cat('the number of entries on the grid, N, is', nrow(grid))

    sim_scenario <- merge(sim_scenario,grid,by.x=c('FROM_AGE_INT','TO_AGE_INT'),by.y=c('x','y'),all.x=T)

    stan_data$coordinates <- matrix(NA,nrow(sim_scenario),2)
    stan_data$coordinates[,1] <- sim_scenario$x_index
    stan_data$coordinates[,2] <- sim_scenario$y_index

    cat("Number of unique age combinations is : ", stan_data$Nu_pairs)
    stan_data$L <- rep(args$B*max(stan_data$ages[,1],stan_data$ages[,2]),2) # make grid symmetric for ages of sources/recipients


  }
  # create mapping for multiple sources per recipient
  stan_data$pt_idx <- factor(sim_scenario$TO_ID,levels=unique(sim_scenario$TO_ID),labels=seq(1,length(unique(sim_scenario$TO_ID)),1))
  my.matrix <- matrix(0, nrow=length(unique(sim_scenario$TO_ID)), ncol=nrow(sim_scenario))
  for(i in 1:length(unique(sim_scenario$TO_ID))) {my.matrix[i,stan_data$pt_idx==i]=1}
  stan_data$pt_map <- my.matrix

  # create index for source/recipient categories
  sim_scenario[, FROM_BIN_COV:= as.numeric(factor(BIN_COV,levels=c('cat1','cat2')))] # set recipient binary covariate to same as source
  sim_scenario[, TO_BIN_COV:= as.numeric(factor(BIN_COV,levels=c('cat1','cat2')))]
  stan_data$S <- length(unique(c(sim_scenario$FROM_BIN_COV,sim_scenario$TO_BIN_COV)))
  if(args$src_cat=='age'){
    stan_data$S <- length(unique(sim_scenario$FROM_AGE_GP_2))
  }
  if(grepl('perfectpred',args$job_tag)) stan_data$S <- length(unique(sim_scenario$TRANS_STAGE))
  stan_data$idx_src <- matrix(NA,stan_data$S,stan_data$N)
  stan_data$idx_rec <- matrix(NA,stan_data$S,stan_data$N)
  for(a in 1:stan_data$S){
    stan_data$idx_src[a,] <- as.numeric(sim_scenario$FROM_BIN_COV==a)
    stan_data$idx_rec[a,] <- as.numeric(sim_scenario$TO_BIN_COV==a)
  }
  src_map <- unique(subset(sim_scenario,select=c('BIN_COV','FROM_BIN_COV')))

  if(args$src_cat=='age'){
    for(a in 1:stan_data$S){
      stan_data$idx_src[a,] <- as.numeric(sim_scenario$FROM_AGE_GP_2==a)
      stan_data$idx_rec[a,] <- as.numeric(sim_scenario$TO_AGE_GP_2==a)
    }
    src_map <- unique(subset(sim_scenario,select=c('FROM_AGE_GP','FROM_AGE_GP_2')))
    src_map <- src_map[order(FROM_AGE_GP_2),]
  }
  if(grepl('perfectpred',args$job_tag)){
    sim_scenario[, TRANS_STAGE_2:= factor(TRANS_STAGE,levels=c('Diagnosed','Undiagnosed'),labels=c(1,2))]
    for(a in 1:stan_data$S){
      stan_data$idx_src[a,] <- as.numeric(sim_scenario$TRANS_STAGE_2==a)
    }
    src_map <- unique(subset(sim_scenario,select=c('TRANS_STAGE','TRANS_STAGE_2')))
    src_map <- src_map[order(TRANS_STAGE_2),]
  }
  
  if(grepl('GP',args$stanModelFile)){
    sim_scenario[, FROM_AGE_STD:= (FROM_AGE_INT - mean(FROM_AGE_INT))/sd(FROM_AGE_INT)]
    sim_scenario[, TO_AGE_STD:= (TO_AGE_INT - mean(TO_AGE_INT))/sd(TO_AGE_INT)]
    if(grepl('_agesrcrec_',args$job_tag)){
    tmp <- data.table(FROM_AGE_STD = sort(unique(sim_scenario$FROM_AGE_STD)))
    tmp[, IDX_UNIQUE_FROM_AGE := seq_len(nrow(tmp))]
    sim_scenario <- merge(sim_scenario,tmp,by='FROM_AGE_STD',all.x=T)
    sim_scenario <- sim_scenario[order(PAIR_ID)]
    stan_data$ux_src = tmp$FROM_AGE_STD
    stan_data$Nux_src = nrow(tmp)
    stan_data$x_2_ux_src = sim_scenario$IDX_UNIQUE_FROM_AGE
    
    tmp <- data.table(TO_AGE_STD = sort(unique(sim_scenario$TO_AGE_STD)))
    tmp[, IDX_UNIQUE_TO_AGE := seq_len(nrow(tmp))]
    sim_scenario <- merge(sim_scenario,tmp,by='TO_AGE_STD',all.x=T)
    sim_scenario <- sim_scenario[order(PAIR_ID)]
    stan_data$ux_rec = tmp$TO_AGE_STD
    stan_data$Nux_rec = nrow(tmp)
    stan_data$x_2_ux_rec = sim_scenario$IDX_UNIQUE_TO_AGE
    
    stan_data$L <- rep(args$B*max(sim_scenario$FROM_AGE_STD,sim_scenario$FROM_AGE_STD),2)
    stan_data$D <- 2
  }else if(grepl('_agerec_',args$job_tag)){
    tmp <- data.table(TO_AGE_STD = sort(unique(sim_scenario$TO_AGE_STD)))
    tmp[, IDX_UNIQUE_TO_AGE := seq_len(nrow(tmp))]
    sim_scenario <- merge(sim_scenario,tmp,by='TO_AGE_STD',all.x=T)
    sim_scenario <- sim_scenario[order(PAIR_ID)]
    stan_data$ux= tmp$TO_AGE_STD
    stan_data$Nux = nrow(tmp)
    stan_data$x_2_ux = sim_scenario$IDX_UNIQUE_TO_AGE
    stan_data$L <- args$B*max(sim_scenario$TO_AGE_STD)
    stan_data$M <- stan_data$M[2]
    stan_data$D <- 1
    stan_data$c <- args$B
  }else{
    tmp <- data.table(FROM_AGE_STD = sort(unique(sim_scenario$FROM_AGE_STD)))
    tmp[, IDX_UNIQUE_FROM_AGE := seq_len(nrow(tmp))]
    sim_scenario <- merge(sim_scenario,tmp,by='FROM_AGE_STD',all.x=T)
    sim_scenario <- sim_scenario[order(PAIR_ID)]
    stan_data$ux = tmp$FROM_AGE_STD
    stan_data$Nux = nrow(tmp)
    stan_data$x_2_ux = sim_scenario$IDX_UNIQUE_FROM_AGE
    stan_data$L <- args$B*max(sim_scenario$FROM_AGE_STD)
    stan_data$M <- stan_data$M[1]
    stan_data$D <- 1
    stan_data$c <- args$B
  }
}

  ## init values
  stan_init <- list()
  logit_y_mix_0 <- logit(c(0.05,0.1,0.3,0.5))
  log_alpha1_pair <- rep(median(dpo$log_alpha1_pair),4)
  log_phi_pair <- rep(median(dpo$log_alpha1_pair),4)
  lscale <- c(0.5,0.7,0.9,1.1)
  stan_init <- function() list(logit_y_mix_0=logit_y_mix_0,
                               lscale=lscale)

  ## save image before running Stan
  tmp <- names(.GlobalEnv)
  tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
  save(list=tmp, file=paste0(outfile.base,'_stanin.RData') )

  # save stan_data object
  rstan::stan_rdump( names(stan_data), file=paste0(outfile.base, '_cmdstanin.R'), envir=list2env(stan_data))


## fit model ----
cat(" \n -------------------------------- \n fit model if running locally \n -------------------------------- \n")
if(args$local==1){
  options(mc.cores = parallel::detectCores())
  sim_mixture_compiled <- cmdstanr::cmdstan_model(args$file_stanModel,
                                                force_recompile = TRUE,
                                                include_paths = dirname(args$file_stanModel)
  )

  # run Stan using cmdstan
  model_fit <- sim_mixture_compiled$sample(
    data = stan_data,
    iter_warmup = 5e2,
    iter_sampling = 2e3,
    refresh = 100,
    parallel_chains = 4,
    chains = 4,
    adapt_delta = 0.9,
    save_warmup = TRUE,
    init = stan_init
  )

  tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
  cat("\n Save fitted data to file ", tmp , "\n")
  model_fit$save_object(file = tmp)
}else{
 # save stan.data object
  rstan::stan_rdump( names(stan_init), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstaninit.R')), envir=list2env(stan_init))
  rstan::stan_rdump( names(stan_data), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstanin.R')), envir=list2env(stan_data))
}
