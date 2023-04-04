
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
require(ggpubr)
require(ggExtra)
require(cowplot)
require(Hmisc)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    pairs.dir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/agegps_updated_criteria_MSM-2010_2022',
    #outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_TE16_MSM-2010_2022-1665619',
    outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022'
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

# reclassify age groups ----

do <- do[order(PAIR_ID),]

do[, FROM_AGE_GP:= cut(FROM_AGE,breaks=c(15,30,40,50,100),include.lowest=T,right=F,
                       labels=c('15-29','30-39','40-49','50+'))]
do[, TO_AGE_GP:= cut(TO_AGE,breaks=c(15,30,40,50,100),include.lowest=T,right=F,
                     labels=c('15-29','30-39','40-49','50+'))]

#GP PAF for each age band----
# statify by age gp of recipient
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w','\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_AGE_GP','TRANS_STAGE'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','FROM_AGE_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = c('draw'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_AGE_GP')
]
po <- dcast.data.table(po, FROM_AGE_GP~stat, value.var = 'q')
po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]
saveRDS(po,file=paste0(outfile.base,'-PAF_agegp','.RDS'))
write.csv(po,file=paste0(outfile.base,'-PAF_agegp','.csv'))

pal <- pal_npg("nrc")(4)[c(1,3,4)]
po[, TO_AGE_GP:='Any age']
g1 <- ggplot(subset(po)) + geom_bar(aes(x=TO_AGE_GP,y=M,fill=FROM_AGE_GP),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_AGE_GP,ymin=CL, ymax=CU,fill=FROM_AGE_GP),width=0.5, colour="black",position=position_dodge(width=0.9))	+
  scale_fill_npg() +
  labs(x='', y='Contribution to transmission',
       fill='Age of source at infection date') +
  theme_bw() +
  theme(legend.pos='none',
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.margin=margin(t=-10)) + #,
  coord_cartesian(ylim = c(0,0.5)) +
  scale_y_continuous(expand=c(0,0),labels = scales::label_percent(accuracy = 1L),breaks=seq(0,0.5,0.1))


#GP PAF recipients ----
# statify by age gp of recipient
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w','\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_AGE_GP','TO_AGE_GP','TRANS_STAGE'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','FROM_AGE_GP','TO_AGE_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_AGE_GP')]
po <- merge(po, tmp, by = c('draw','TO_AGE_GP'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('TO_AGE_GP','FROM_AGE_GP')
]
po <- dcast.data.table(po, TO_AGE_GP+FROM_AGE_GP~stat, value.var = 'q')
po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]
saveRDS(po,file=paste0(outfile.base,'-PAF_stratify_agegp','.RDS'))
write.csv(po,file=paste0(outfile.base,'-PAF_stratify_agegp','.csv'))

g2 <- ggplot(subset(po)) + geom_bar(aes(x=TO_AGE_GP,y=M,fill=FROM_AGE_GP),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_AGE_GP,ymin=CL, ymax=CU,fill=FROM_AGE_GP),width=0.5, colour="black",position=position_dodge(width=0.9))	+
  scale_fill_npg() +
  labs(x='Age of new Amsterdam MSM case', y='\nContribution to transmission',fill='Age of source at infection date') +
  theme_bw() +
  theme(legend.pos='bottom',
        legend.margin=margin(t=-10)) + #,
  coord_cartesian(ylim = c(0,0.5)) +
  scale_y_continuous(expand=c(0,0),labels = scales::label_percent(accuracy = 1L),breaks=seq(0,0.5,0.1))
g2

legend_t <- cowplot::get_legend(g2 + theme(legend.position = "bottom"))

g <- ggarrange(g1 + theme(legend.margin=margin(t=-3))  ,g2+ theme(legend.position='none',legend.margin=margin(t=-3)),ncol=2,widths=c(0.35,0.65),align='hv',
               labels='AUTO',font.label = list(size=14),vjust=1)
g_paf <- ggarrange(g, legend_t,ncol=1,heights=c(0.8,0.2))
g_paf

## boxplot of age difference ----

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w','\\[([0-9]+)\\]'),'\\1',as.character(variable)))]

tmp <- subset(do, select = c('PAIR_ID','FROM_AGE','TO_AGE','TO_AGE_GP'))

po <- merge(po, tmp, by = 'PAIR_ID')
po[, age_diff:= FROM_AGE - TO_AGE]
po[, age_diff_int:= round(age_diff,0)]

po <- po[, list(value = sum(value)), by = c('draw','age_diff_int','TO_AGE_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_AGE_GP')]
po <- merge(po, tmp, by = c('draw','TO_AGE_GP'))
po[, paf := value/total]
po <- po[order(draw,TO_AGE_GP,age_diff_int),]

# get cumsum of paf per recipient age group
po <- po[, list(age_diff=age_diff_int,
                 paf = paf,
                cumpaf=cumsum(paf)),by=c('draw','TO_AGE_GP')]
po <- po[, list(q = c(max(age_diff[cumpaf<=0.5]),
                      max(age_diff[cumpaf<=0.25]),
                      max(age_diff[cumpaf<=0.75]),
                      min(age_diff),
                      max(age_diff),
                      max(age_diff[cumpaf<=0.025]),
                      max(age_diff[cumpaf<=0.975])),
                stat = c('M','IL', 'IU', 'min','max','CL','CU')), by=c('draw','TO_AGE_GP')]

po <- po[,
list( q = quantile(q, probs = 0.5)),
         by = c('TO_AGE_GP','stat')
]

q.palette <- RColorBrewer::brewer.pal(6,'Blues')[2:6]

ss <- do[, list(N=length(PAIR_ID)), by=c('TO_AGE_GP')]
po <- merge(po,ss,by='TO_AGE_GP')
po[, prop:= N/sum(N)]

po <- dcast.data.table(po, TO_AGE_GP+N~stat, value.var = 'q')

write.csv(po,file=paste0(outfile.base,'-sources_agegap.csv'))
saveRDS(po,file=paste0(outfile.base,'-sources_agegap.rds'))

p3 <- ggplot(po) +
  geom_hline(yintercept=0,linetype='dashed') +
  geom_boxplot(data=subset(po,TO_AGE_GP=='15-29'), aes(x=TO_AGE_GP, ymin=CL, ymax=CU, lower=IL, upper=IU, middle=M, fill=TO_AGE_GP, weight=sqrt(N)),
               stat = "identity",width=0.3)+
  geom_boxplot(data=subset(po,TO_AGE_GP=='30-39'), aes(x=TO_AGE_GP, ymin=CL, ymax=CU, lower=IL, upper=IU, middle=M, fill=TO_AGE_GP, weight=sqrt(N)),
               stat = "identity",width=0.26)+
  geom_boxplot(data=subset(po,TO_AGE_GP=='40-49'), aes(x=TO_AGE_GP, ymin=CL, ymax=CU, lower=IL, upper=IU, middle=M, fill=TO_AGE_GP, weight=sqrt(N)),
               stat = "identity",width=0.30)+
  geom_boxplot(data=subset(po,TO_AGE_GP=='50+'), aes(x=TO_AGE_GP, ymin=CL, ymax=CU, lower=IL, upper=IU, middle=M, fill=TO_AGE_GP, weight=sqrt(N)),
               stat = "identity",width=0.12)+
  scale_fill_manual(values=q.palette) +
  scale_y_continuous(breaks=seq(-50,50,5),labels=seq(-50,50,5)) +
  labs(x="Age of new Amsterdam MSM case\n", y = "Age difference of\nprobable sources", fill='') +
  theme_bw() +
  theme(legend.position='none')
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-sources_agegap_boxplots.png'), p3, w = 6, h = 6)


g <- ggarrange(g1 + labs(y='Contribution to transmission') + theme_bw(base_size=16) + theme(legend.position='none',legend.margin=margin(t=-3),axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) + guides(fill = guide_legend(title.position = 'top')),
               g2 + labs(y='Contribution to transmission') + theme_bw(base_size=16) + theme(legend.position='none',legend.margin=margin(t=-3),axis.title.x=element_blank()) + guides(fill = guide_legend(title.position = 'top')),
               p3 + labs(y = "Age difference between\nphylogenetically likely source\nand incident case") + theme_bw(base_size=16) + theme(axis.title.x=element_blank(),legend.position='none'),
               ncol=3,widths=c(0.25,0.35,0.35),align='hv',
               labels='AUTO',font.label = list(size=14),vjust=1)
g <- annotate_figure(g,
                     bottom = text_grob("Age of new Amsterdam MSM case",
                                        hjust = 0.5))
g <- ggarrange(g,legend_t,ncol=1,heights=c(0.85,0.15),align='v')
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-panel_sources_boxplot_agediff.pdf'), g, w = 10.5, h = 5)

