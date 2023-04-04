
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
require(grid)
require(ggpubr)
require(ggExtra)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    pairs.dir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time/out_Amsterdam/agegps_updated_criteria_MSM-2010_2022',
    outdir= '/Users/alexb/Documents/GitHub/source.attr.with.infection.time/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_TE16_MSM-2010_2022-1665619',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_bgUnif_piGP_221027_KD',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'agegps_TE16_MSM-2010_2022'
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

do <- do[order(PAIR_ID),]

do[, FROM_AGE_GP:= cut(FROM_AGE,breaks=c(15,30,40,50,100),include.lowest=T,right=F,
                       labels=c('15-29','30-39','40-49','50+'))]

tmp <- do[, list(N=length(unique(FROM_SEQUENCE_ID))),by='TO_SEQUENCE_ID']
mean <- data.table(x=mean(tmp$N),y=-5)
count <- tmp[, list(count=length(TO_SEQUENCE_ID)),by='N']

# load all pairs
pairs <- readRDS(file=file.path('out_Amsterdam',args$pairs_dir, 'all_pairs.rds'))
tmp2 <- pairs[, list(N=length(unique(FROM_SEQUENCE_ID))),by='TO_SEQUENCE_ID']
mean2 <- data.table(x=mean(tmp2$N),y=-5)

# plot cluster sizes
pal <- pal_npg("nrc")(9)
p1 <- ggplot(tmp) +
  geom_histogram(aes(x=N),fill=pal[2],col=pal[2],binwidth=2,alpha=0.5) +
  geom_point(data=mean, aes(y=y, x=x, colour='Average number of plausible sources (following exclusions)'), alpha=1,pch=6, size=2,stroke=1)+
  theme_bw(base_size=16) + labs(x='Number of plausible sources per recipient',y='Number of recipients') +
  scale_x_continuous(breaks=seq(0,120,10),labels=c(seq(0,120,10))) +
  scale_y_continuous(limits=c(-12, 250), expand = c(0, 0)) +
  scale_colour_manual(values=c(pal[1],pal[2]),
                      labels=c('Average number of plausible sources',
                               'Average number of plausible sources (following exclusions)'),name='') +
  guides(fill = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 2)) +
  theme(legend.position='bottom')

# dodge bars
tmp[, data:= 'after_excl']
tmp2[, data:= 'before_excl']
tmp <- rbind(tmp,tmp2)
tmp[N>=50, N:= 50]

tmp <- tmp[, list(y=length(TO_SEQUENCE_ID)),by=c('data','N')]
po <- expand.grid(data=c('before_excl','after_excl'),N=seq(1,50))
tmp <- merge(tmp,po,all=T)

mean[, data:= 'after_excl']
mean2[, data:= 'before_excl']
mean <- rbind(mean,mean2)

tmp[, data:= factor(data,levels=c('before_excl','after_excl'))]
mean[, data:= factor(data,levels=c('before_excl','after_excl'))]

p1 <- ggplot(tmp) +
  geom_bar(aes(x=N,y=y,fill=data),stat = "identity", position = "dodge",alpha=0.75) +
  geom_point(data=mean, aes(y=y, x=x, colour=data), alpha=1,pch=6, size=2,stroke=1)+
  theme_bw(base_size=16) + labs(x='Phylogenetically plausible sources per recipient',y='Number of incident cases') +
  scale_x_continuous(breaks=seq(0,50,10),labels=c(seq(0,40,10),expression("">=50))) +
  scale_y_continuous(limits=c(-10, 140), expand = c(0, 0)) +
  scale_fill_manual(values=c('grey50',pal[2]),
                      labels=c('Plausible sources',
                               'Plausible sources following exclusion criteria'),name='') +
  scale_colour_manual(values=c('grey50',pal[2]),
                      labels=c('Average number of sources',
                               'Average number of sources'),name='') +
  guides(fill = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 2)) +
  theme(legend.position='bottom')
ggsave(file=paste0(outfile.base,'-Amsterdam_cluster_sizes_before_after_exclusions_dodge.png'), p1, w=7, h=5)
ggsave(file=paste0(outfile.base,'-Amsterdam_cluster_sizes_before_after_exclusions_dodge.pdf'), p1, w=7, h=5)


# plot distances

dps_clock <- readRDS(file = file.path(in.dir,'clock_quantiles.rds'))

do2 <- copy(do)

pal_3 <- pal_npg("nrc")(4)[c(1,3)]
p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
p.alpha <- 0.7
q.palette <- RColorBrewer::brewer.pal(6,'Blues')[2:6]

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
  geom_point(data=do2,aes(x=TIME_ELAPSED,y=GEN_DIST,colour=FROM_AGE_GP)) +
  theme_bw(base_size=16)+
  scale_colour_manual(values=q.palette) +
  labs(x='\n Time elapsed (in years)',y='Patristic distance of\nphylogenetically possible pair \n',colour='Age of probable\ntransmitter')+
  theme(legend.position='bottom') +#,
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.16,0.02),labels=scales::label_percent(accuracy = 1L)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,16,2)) +
  coord_cartesian(xlim=c(0,16)) +
  guides(color = guide_legend(override.aes = list(size = 5) ) )
p2
ggsave(file=file.path(out.dir,paste0('distances_signalcone','.pdf')), p2, w=6, h=8)
ggsave(file=file.path(out.dir,paste0('distances_signalcone','.png')), p2, w=6, h=8)


g <- ggarrange(p1,p2 + labs(colour='Age of probable transmitter'),ncol=1,align='hv',
               labels=c('B','C'),font.label=list(size=20))
ggsave(file=file.path(out.dir,paste0('Amsterdam_cluster_sizes_distances','.pdf')), g, w=9, h=11)


#plot ages

p2 <- ggplot(do, aes(x=FROM_AGE_INT, y=TO_AGE_INT)) + geom_point(colour=pal[2],alpha=0.6) +
  labs(x='Age of source at\nestimated date of infection\nof recipient', y= 'Age of recipient at\nestimated date of infection') +
  scale_x_continuous(breaks=seq(20,80,10),labels=seq(20,80,10)) +
  scale_y_continuous(breaks=seq(20,80,10),labels=seq(20,80,10)) + theme_bw(base_size=15)
p2 <- ggExtra::ggMarginal(p2, type = "histogram",fill=pal[2],alpha=0.7,col=pal[2])

p <- ggarrange(p1,p2,nrow=1,labels=c('B','C'),font.label=list(size=20),vjust=1,widths=c(0.5,0.5))
ggsave(file=paste0(outfile.base,'-Amsterdam_cluster_sizes_ages.png'), p, w=10, h=5)
ggsave(file=paste0(outfile.base,'-Amsterdam_cluster_sizes_ages.pdf'), p, w=10, h=5)
