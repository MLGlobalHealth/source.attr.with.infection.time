require(tidyverse)
require(data.table)
require(ggplot2)
require(abind)
require(rstan)
require(knitr)
require(bayesplot)
require(ggsci)
require(ggpubr)
require(cmdstanr)
require(EMCluster)

# setup args ----
args <- list(
  #source_dir = '~/Documents/GitHub/source.attr.with.infection.time.fork',
  source_dir = '~/Documents/GitHub/source.attr.with.infection.time.public',
  indir = '~/Box\ Sync/Roadmap/source_attribution',
  outdir = 'out_Amsterdam',
  #pairs_dir = 'agegps_updated_criteria_210216_MSM-2010_2022',
  analysis = 'analysis_220713',
  pairs_dir = 'agegps_sensanalysis_210216_MSM-2010_2022',
  job_tag = 'agegps_sensanalysis_210216_MSM',
  #job_tag = 'sens_bg_2DMM_MSM',
  #job_tag = 'sens_bg_lnorm_MSM',
  trsm = 'MSM',
  clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
  #stanModelFile = 'mm_sigHierG_bgUnif_piVanilla_220408b', # vanilla model
  #stanModelFile = 'mm_sigHierG_bgUnif_piReg_230111b', # covariate model
  stanModelFile = 'mm_bgUnif_piGP_221027b', # 2D HSGP model
  #stanModelFile = 'mm_bglnorm_piGP_240624', # 2D HSGP model
  #stanModelFile = 'mm_bgUnif_pi1DGP_Ams_230224', # 1D HSGP model
  #stanModelFile = 'mm_bgUnif_pi1DGP_Ams_230224b', # 2 * 1D HSGP model
  #stanModelFile = 'mm_bgMVMM_piGP_240710', # 2D HSGP model
  hmc_stepsize = 0.02,
  hmc_num_samples = 15,
  hmc_num_warmup = 10,
  seed = 42,
  chain = 1,
  reps = 1,
  bg = 'unif',
  local = 1,
  time_period='2010-2022',
  m1 = 24,
  m2 = 24,
  B = 1.2
)

in.dir <- file.path(args$outdir,args$pairs_dir)

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

r <- 1

# load pairs data ----

pairs <- fread(file.path('data_Ams',args$analysis,paste0(args$trsm,"_anonymised_pairs.csv")))

tmp <- pairs[, list(N=length(unique(FROM_SEQUENCE_ID))),by='TO_SEQUENCE_ID']

cat(paste0("number of recipients (true pairs) = ", length(unique(pairs$TO_SEQUENCE_ID))))
cat(paste0("mean number of sources per recipient = ", mean(tmp$N)))
cat(paste0("median number of sources per recipient = ", median(tmp$N)))
cat(paste0("number of total pairs = ", nrow(pairs)))

# Prepare data
do <- pairs
do[,LOG_TIME_ELAPSED:= log(TIME_ELAPSED)]

# quantiles of time elapsed
quantile(do$TIME_ELAPSED,probs=c(0.5,0.025,0.975))

# for how many pairs is only one direction possible?
dp <- subset(do,select=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID'))
dp[, FROM_SEQUENCE_ID:= as.numeric(FROM_SEQUENCE_ID)]
dp[, TO_SEQUENCE_ID:= as.numeric(TO_SEQUENCE_ID)]
dp <- dp[, list(pairs= paste0(min(FROM_SEQUENCE_ID,TO_SEQUENCE_ID),'-', max(FROM_SEQUENCE_ID,TO_SEQUENCE_ID))),by=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID')]
cp <- subset(dp,select=c('pairs'))
cat(paste0('proportion of IDs in which only one direction is possible: ', round(nrow(unique(cp))/nrow(dp)*100,0),'%'))

# define groupings for model ----
do[,FROM_AGE_INT:= round(FROM_AGE)]
do[,TO_AGE_INT:= round(TO_AGE)]
do[,FROM_AGE_GP_2:= factor(FROM_AGE_GP,labels=c(1,2,3,4,5))]
do[,TO_AGE_GP_2:= factor(TO_AGE_GP,labels=c(1,2,3,4,5))]

# get unique time elapsed

tmp <- data.table(TIME_ELAPSED = sort(unique(do$TIME_ELAPSED)))
tmp[, IDX_UNIQUE_TE := seq_len(nrow(tmp))]
do <- merge(do, tmp, by = 'TIME_ELAPSED')
do <- do[order(PAIR_ID),]

# get unique pairs ages
tmp <- data.table(unique(cbind(do$FROM_AGE_INT,do$TO_AGE_INT)))
tmp[, IDX_UNIQUE_PAIR := seq_len(nrow(tmp))]
setnames(tmp, "V1", "FROM_AGE_INT")
setnames(tmp, "V2", "TO_AGE_INT")
do <- merge(do, tmp, by = c("FROM_AGE_INT","TO_AGE_INT"))

# re-order by ages
do <- do[order(FROM_AGE_INT,TO_AGE_INT),]
do <- data.table(do)
do[, PAIR_ID := seq(1,nrow(do)) ]

# load signal cone from clock model ----

dps_clock <- readRDS(file = file.path(in.dir,'clock_quantiles.rds'))

do2 <- copy(do)

pal_3 <- pal_npg("nrc")(4)[c(1,3)]
p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
p.alpha <- 0.7

p <- ggplot(data=dps_clock,aes(x=d_TSeqT)) +
  geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=do2,aes(x=TIME_ELAPSED,y=GEN_DIST,colour=FROM_AGE_GP)) +
  theme_bw(base_size=26)+
  scale_colour_npg() +
  labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',colour='Age of probable\ntransmitter')+
  theme(legend.position='bottom',
        legend.text = element_text(size=18))+
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.18,0.02),labels=scales::label_percent(accuracy = 1L)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,18,2))
p
ggsave(file=file.path(out.dir,paste0('global_cluster_in_cone_age_',args$time_period,'.pdf')), p, w=11, h=9)
ggsave(file=file.path(out.dir,paste0('global_cluster_in_cone_age_',args$time_period,'.png')), p, w=11, h=9)


p2 <- ggplot(data=dps_clock,aes(x=d_TSeqT)) +
  geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=do,aes(x=TIME_ELAPSED,y=GEN_DIST,colour=FROM_AGE_GP)) +
  theme_bw(base_size=26)+
  scale_colour_npg() +
  labs(x='\n Time elapsed (in years)',y='',colour='Age of probable\ntransmitter')+
  theme(legend.position='bottom')+
  coord_cartesian(xlim=c(0,10),ylim=c(0,0.08)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.08,0.01),labels=scales::label_percent(accuracy = 1L)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,10,2))
p2

g <- ggarrange(p,p2,ncol=2,align='hv')
ggsave(file=file.path(out.dir,paste0('global_cluster_in_cone_age_panel_',args$time_period,'.pdf')), g, w=25, h=10)

# load estimated molecular clock ----
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

log_alpha1 <- median(po$log_alpha1)
log_alpha1_pair_sd <- median(po$log_alpha1_pair_sd)
log_phi <- median(po$log_phi)
log_phi_pair_sd <- median(po$log_phi_pair_sd)

# plot observed sources ----

tmp <- do[, list(N=length(TO_AGE_GP)),by='FROM_AGE_GP']
tmp <- tmp[, list(FROM_AGE_GP=FROM_AGE_GP,pct=N/sum(N))]
tmp[, TO_AGE_GP:= 'Overall']
g1 <- ggplot(subset(tmp)) + geom_bar(aes(x=TO_AGE_GP,y=pct,fill=FROM_AGE_GP),stat='identity',position=position_dodge(width=0.9)) +
  scale_fill_npg() +
  labs(x='Age of recipient', y='Proportion of attributable\ninfections to age group',fill='Age of likely source') +
  theme_bw(base_size=28) +
  theme(legend.position='bottom') + #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
ggsave(file = paste0(outfile.base,'-obs_sources_age_source.png'), g1, w = 10, h = 8)

tmp <- do[, list(N=length(PAIR_ID)),by=c('FROM_AGE_GP','TO_AGE_GP')]
tmp <- tmp[, list(FROM_AGE_GP=FROM_AGE_GP,pct=N/sum(N)),by='TO_AGE_GP']
g2 <- ggplot(subset(tmp)) + geom_bar(aes(x=TO_AGE_GP,y=pct,fill=FROM_AGE_GP),stat='identity',position=position_dodge(width=0.9)) +
  scale_fill_npg() +
  labs(x='Age of recipient', y='Proportion of attributable\ninfections to age group',fill='Age of likely source') +
  theme_bw(base_size=28) +
  theme(legend.position='bottom') + #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
g2
ggsave(file = paste0(outfile.base,'-obs_sources_age_source_recipient.png'), g2, w = 10, h = 8)

legend_t <- cowplot::get_legend(g2 + theme(legend.position = "bottom"))

g <- ggarrange(g1 + rremove("xlab")+ theme(legend.position='none'),g2+ theme(legend.position='none'),ncol=2,widths=c(0.35,0.65),align='hv')
g <- ggarrange(g, legend_t,ncol=1,heights=c(0.8,0.2))
ggsave(file = paste0(outfile.base,'-obs_sources_age_source_recipient_panel.png'), g, w = 16, h = 8)
ggsave(file = paste0(outfile.base,'-obs_sources_age_source_recipient_panel.pdf'), g, w = 16, h = 8)

saveRDS(dps_clock,file = paste0(outfile.base,'-clock_quantiles.rds'))

# define stan data ----

stan_data <- list()
stan_data$N <- nrow(do)
stan_data$P <- length(unique(do$TO_SEQUENCE_ID))
stan_data$y <- do$GEN_DIST
stan_data$S <- length(unique(do$FROM_AGE_GP))
stan_data$x <- do$TIME_ELAPSED
stan_data$pt_idx <- factor(do$TO_SEQUENCE_ID,levels=unique(do$TO_SEQUENCE_ID),labels=seq(1,length(unique(do$TO_SEQUENCE_ID)),1))
tmp <- unique(do, by = 'IDX_UNIQUE_TE')[order(IDX_UNIQUE_TE)]
stan_data$Nux = nrow(tmp)
stan_data$ux = tmp$TIME_ELAPSED
age_vec <- seq(1, 90)
stan_data$a <- max(age_vec)

# posterior medians from clock model
stan_data$log_alpha1 <- log_alpha1
stan_data$log_alpha1_pair_sd <- log_alpha1_pair_sd
stan_data$log_phi <- log_phi
stan_data$log_phi_pair_sd <- log_phi_pair_sd

stan_data$obs_to_age_src_idx <- do$FROM_AGE_GP_2
stan_data$obs_to_age_rcp_idx <- do$TO_AGE_GP_2

tmp <- do[, list(PAIR_ID), by = c("FROM_AGE_INT")]
tmp2 <- tmp[, list(num = length(PAIR_ID)), by = c("FROM_AGE_INT")]
setkey(tmp2, "FROM_AGE_INT")
tmp <- data.table(age = rep(0,90))
tmp[which(!(age_vec %in% tmp2$FROM_AGE_INT)), num:=0]
tmp[which(age_vec %in% tmp2$FROM_AGE_INT), num:=tmp2$num]
stan_data$units_to_obs_length <- tmp$num #Units_to_obs_length
stan_data$Nk_max <- max(stan_data$units_to_obs_length)

my.matrix <- matrix(0, nrow=length(unique(do$TO_SEQUENCE_ID)), ncol=nrow(do))
for(i in 1:length(unique(do$TO_SEQUENCE_ID))) {my.matrix[i,stan_data$pt_idx==i]=1}
stan_data$pt_map <- my.matrix

# age
tmp <- do[, list(PAIR_ID), by = c("FROM_AGE_INT")]
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
{tmp3 <- subset(do, FROM_AGE_INT == a)
if(nrow(tmp3) != 0)
{tmp[a, 1:nrow(tmp3)] <- tmp3[,PAIR_ID]}
}
stan_data$age_to_obs_idx <- tmp

stan_data$A <- length(unique(do$FROM_AGE_GP)) # coarse age groups
stan_data$obs_to_age_idx <- do$FROM_AGE_GP

# create an array with the pair IDs for sources of age a and recipients of age j
mat <- list()
for (a in 1:stan_data$a){
  tmp <- subset(do, FROM_AGE_INT == a)
  mat[[a]] <- list()
  if(nrow(tmp) != 0){
    for (j in 1:stan_data$a){
      tmp2 <- subset(do, FROM_AGE_INT == a & TO_AGE_INT == j)
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
    ind_length <- length(which(do$FROM_AGE_INT == a & do$TO_AGE_INT == j)) # number of pairs in entry a,j
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
  tmp <- subset(do, FROM_AGE_GP_2 == a)
  mat[[a]] <- list()
  if(nrow(tmp) != 0){
    for (j in 1:stan_data$A){
      tmp2 <- subset(do, FROM_AGE_GP_2 == a & TO_AGE_GP_2 == j)
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
    ind_length <- length(which(do$FROM_AGE_GP_2 == a & do$TO_AGE_GP_2 == j))
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
stan_data$L <- rep(args$B*max(do$FROM_AGE_INT,do$TO_AGE_INT),2) # make grid symmetric for ages of sources/recipients
stan_data$M <- c(args$m1,args$m2) #Maybe increase to 60
stan_data$M_nD <- stan_data$M[1] * stan_data$M[2] #900, maybe increase to 3600
tmp <- unique(do, by = 'IDX_UNIQUE_TE')[order(FROM_AGE)]
stan_data$sd1 <- sd(tmp$FROM_AGE)

tmp <- unique(do, by = 'IDX_UNIQUE_TE')[order(TO_AGE)]
stan_data$sd2 <- sd(tmp$TO_AGE)

stan_data$Nu_pairs <- max(do$IDX_UNIQUE_PAIR)
stan_data$IDX_UNIQUE_PAIR <- do$IDX_UNIQUE_PAIR

stan_data$n <- 90 #75 - 17 + 1
stan_data$m <- 90 #75 - 17 + 1

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

do <- merge(do,grid,by.x=c('FROM_AGE_INT','TO_AGE_INT'),by.y=c('x','y'),all.x=T)

stan_data$coordinates <- matrix(NA,nrow(do),2)
stan_data$coordinates[,1] <- do$x_index
stan_data$coordinates[,2] <- do$y_index

cat("Number of unique age combinations is : ", stan_data$Nu_pairs)
stan_data$L <- rep(args$B*max(stan_data$ages[,1],stan_data$ages[,2]),2) # make grid symmetric for ages of sources/recipients

stan_data$idx_src <- matrix(NA,stan_data$S,stan_data$N)
stan_data$idx_rec <- matrix(NA,stan_data$S,stan_data$N)
for(a in 1:stan_data$S){
  stan_data$idx_src[a,] <- as.numeric(do$FROM_AGE_GP_2==a)
  stan_data$idx_rec[a,] <- as.numeric(do$TO_AGE_GP_2==a)
}

if(grepl('_agesrcrec_',args$job_tag)){
  tmp <- data.table(FROM_AGE_STD = sort(unique(age_idx_std)))
  tmp[, FROM_AGE_INT := seq_len(nrow(tmp))]
  do <- merge(do,tmp,by='FROM_AGE_INT',all.x=T)
  do <- do[order(PAIR_ID)]
  stan_data$ux_src = tmp$FROM_AGE_STD
  stan_data$Nux_src = nrow(tmp)
  stan_data$x_2_ux_src = do$FROM_AGE_INT
  
  tmp <- data.table(TO_AGE_STD = sort(unique(age_idx_std)))
  tmp[, TO_AGE_INT := seq_len(nrow(tmp))]
  do <- merge(do,tmp,by='TO_AGE_INT',all.x=T)
  do <- do[order(PAIR_ID)]
  stan_data$ux_rec = tmp$TO_AGE_STD
  stan_data$Nux_rec = nrow(tmp)
  stan_data$x_2_ux_rec = do$TO_AGE_INT
  
  stan_data$L <- rep(args$B*max(do$FROM_AGE_STD,do$FROM_AGE_STD),2)
  stan_data$D <- 2
}else if(grepl('_agerec_',args$stanModelFile)){
  tmp <- data.table(TO_AGE_INT = sort(unique(age_idx_std)))
  tmp[, TO_AGE_INT := seq_len(nrow(tmp))]
  do <- merge(do,tmp,by='TO_AGE_INT',all.x=T)
  do <- do[order(PAIR_ID)]
  stan_data$ux = tmp$TO_AGE_INT
  stan_data$Nux = nrow(tmp)
  stan_data$x_2_ux = do$TO_AGE_INT
  stan_data$L <- args$B*max(do$TO_AGE_INT)
  stan_data$M <- stan_data$M[1]
  stan_data$D <- 1
  stan_data$c <- args$B
}else{
  tmp <- data.table(FROM_AGE_STD = sort(unique(age_idx_std)))
  tmp[, FROM_AGE_INT := seq_len(nrow(tmp))]
  do <- merge(do,tmp,by='FROM_AGE_INT',all.x=T)
  do <- do[order(PAIR_ID)]
  stan_data$ux = tmp$FROM_AGE_STD
  stan_data$Nux = nrow(tmp)
  stan_data$x_2_ux = do$FROM_AGE_INT
  stan_data$L <- args$B*max(do$FROM_AGE_STD)
  stan_data$M <- stan_data$M[1]
  stan_data$D <- 1
  stan_data$c <- args$B
}

# optional sensitivity analysis for different BG distribution ----

# fit 2D bivariate normal distribution to times and distances
if(grepl('_bgGMM_',args$stanModelFile)){
  
  # fit a mixture of gaussians using EMCluster package
  data <- matrix(c(do$TIME_ELAPSED,do$GEN_DIST),ncol=2)
  n_clusters <- 3  # Number of clusters
  init_params <- init.EM(data, nclass = n_clusters)
  gmm <- emcluster(data,init_params,assign.class=T)

  plot(data, col = gmm$class)
  g <- ggplot(data.table(data)) + geom_point(aes(x=V1,y=V2,col=factor(gmm$class))) + 
    scale_color_discrete() + theme_bw() + labs(x='time elapsed',y='patristic distance',col='cluster membership')
  ggsave(file = paste0(outfile.base,'-BNMM_classifications.png'),g, w = 7, h = 5)
  
  sigma <- LTSigma2variance(gmm$LTSigma)
  
  # store parameters:
  stan_data$K <- n_clusters
  stan_data$mu <- t(gmm$Mu)
  stan_data$gmm_pi <- gmm$pi
  stan_data$Sigma <- array(sigma[, , ], dim = c(ncol(data), ncol(data), n_clusters))
  stan_data$Sigma <- aperm(stan_data$Sigma, c(3, 1, 2))

}

# init values ----
stan_init <- list()
y_mix <- rep(c(0.05,0.1,0.3), 3)
stan_init$y_mix <- y_mix[args$chain]

# save inputs ----
## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=paste0(outfile.base, '_stanin.RData') )

# save stan_data object
rstan::stan_rdump( names(stan_data), file=paste0(outfile.base, '_cmdstanin.R'), envir=list2env(stan_data))

# compile stan model ----
options(mc.cores = parallel::detectCores())
sim_mixture_compiled <- cmdstanr::cmdstan_model(args$file_stanModel,
                                                force_recompile = TRUE,
                                                include_paths = dirname(args$file_stanModel)
)


# run Stan ----
options(mc.cores=parallel::detectCores())

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
  init = function() list(y_mix=y_mix)
)

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Save fitted data to file ", tmp , "\n")
model_fit$save_object(file = tmp)

