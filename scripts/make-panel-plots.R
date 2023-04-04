
require(ggplot2)
require(ggsci)
require(scales)
require(data.table)
require(kableExtra)
require(ggpubr)
require(dplyr)
require(gridExtra)

indir <- '/Users/alexb/Box Sync/Roadmap/source_attribution'
args <- list(model='covariate')
source('R/functions_plotting.R')

## figure 3: 1 competing pair

if(args$model=='vanilla'){
  outfile.base <- '/Users/alexb/Box Sync/Roadmap/source_attribution/mm_sigHierG_bgUnif_piVanilla_220408-simulations_500truepairs_srcbin_exclTE16_subsample50pct-334228/mm_sigHierG_bgUnif_piVanilla_220408-simulations_500truepairs_srcbin_exclTE16_subsample50pct'
}
if(args$model=='covariate'){
  outfile.base <- '/Users/alexb/Box Sync/Roadmap/source_attribution/mm_sigHierG_bgUnif_piReg_230118-simulations_500truepairs_perfectpred_exclTE16_subsample0pct-1155684/mm_sigHierG_bgUnif_piReg_230118-simulations_500truepairs_perfectpred_exclTE16_subsample0pct'
}

sim_scenario <- readRDS(file = paste0(outfile.base,'.rds'))

dps <- readRDS(file = paste0(outfile.base,'-clock_quantiles.rds'))

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
model_fit <- readRDS(file = tmp)

# subfigure A ----
po <- posterior_prob_pair(model_fit,sim_scenario,dps)

g1 <- make_plot_simulated_data_colour_prob_tpair(po,dps,sim_scenario,outfile.base)

# subfigure B ----

# plot false signal

sim_scenario[, d_TSeqT:= round(TIME_ELAPSED,2)]
dps[, d_TSeqT:= round(d_TSeqT,2)]
sim_scenario <- merge(sim_scenario,subset(dps,select=c('d_TSeqT','q2.5','q97.5')), by='d_TSeqT',all.x=T)

# flag the false positives
sim_scenario[TRANSMISSION_PAIR=='No', cat:= ifelse(GEN_DIST>=q2.5 & GEN_DIST<=q97.5,'false_pos',
                                   'true_neg')]
sim_scenario[TRANSMISSION_PAIR=='Yes', cat:= ifelse(GEN_DIST<q2.5 | GEN_DIST>q97.5,'false_neg',
                                                   'true_pos')]
sim_scenario[, cat:=factor(cat,levels=c('true_neg','false_pos','true_pos','false_neg'),
                           labels=c('Unlinked pair with\nno signal',
                                    'Unlinked pair with\nfalse signal',
                                    'Transmission pair\nwith signal',
                                    'Transmission pair\nno signal'))]

tmp <- sim_scenario[, list(p=length(PAIR_ID)/nrow(sim_scenario)), by=c('cat')]
pal <- pal_npg("nrc")(5)

g2 <- ggplot(tmp) +
  geom_bar(aes(x=cat,y=p,fill=cat),stat='identity', position = position_dodge(width=0.9),alpha=0.7) +
  theme_bw(base_size=16) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L)) +
  labs(x='',y='\nPhylogenetically possible pairs') +
  scale_fill_manual(name="",values = c(pal[4],pal[5],pal[2],pal[1]),
                    labels=c('Unlinked pair with\nno signal',
                             'Unlinked pair with\nfalse signal',
                             'Transmission pair\nwith signal',
                             'Transmission pair\nno signal')) +
  theme(legend.pos='none',
        axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.7))

g2_leg <- cowplot::get_legend(g2 + theme(legend.position="bottom",#legend.margin=margin(t=-10,r=-20,b=0,l=-20),
                                           axis.text.x=element_blank(), axis.ticks.x=element_blank()))

# subfigure C ----

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

pal <- pal_npg("nrc")(5)

# plot prob of being a pair
g3 <- make_plot_tpair_violin(po,pal)

g3_leg <- cowplot::get_legend(g3)

# subfigure D ----

# plot PAF
po <- readRDS(file=paste0(outfile.base,'-PAF_GQ_all_pairs','.RDS'))
po[, p_pairs:= 0.5]
po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]

sim_scenario <- readRDS(file = paste0(outfile.base,'.rds')) # use 2 mean sources per case scenario for threshold
sim_scenario[, p_pairs:= 0.5]

sim_scenario[, STAGE_thrsh:= 0]
sim_scenario[GEN_DIST<0.015, STAGE_thrsh:= 1]

sim_scenario[, SOURCE:= factor(BIN_COV,levels=c('cat1','cat2'),labels=c('Category\n1','Category\n2'))]

tot <- sim_scenario[, list(N_thrsh=as.numeric(sum(STAGE_thrsh)),
                           N_truth=as.numeric(sum(TRANSMISSION_PAIR=='Yes'))), by=c('p_pairs')]
sim_scenario <- merge(sim_scenario,tot,by=c('p_pairs'))
tmp <- sim_scenario[, list(p_thrsh=as.numeric(sum(STAGE_thrsh)/unique(N_thrsh)),
                           p_truth=as.numeric(sum(TRANSMISSION_PAIR=='Yes')/unique(N_truth))), by=c('p_pairs','SOURCE')]

if(is.null(po$BIN_COV)) po[, BIN_COV:= factor(TRANS_STAGE,levels=c('Undiagnosed','Diagnosed'),labels=c('cat1','cat2'))]
g4 <- plot_estimated_flows_p(po,tmp)

g4_leg <- cowplot::get_legend(g4 + theme(legend.text=element_text(size=12)))

po[,SOURCE:=factor(BIN_COV,levels=c('cat1','cat2'),labels=c('Source\ncategory 1','Source\ncategory 2'))]
po[,SOURCE:=factor(BIN_COV,levels=c('cat1','cat2'),labels=c('Category\n1','Category\n2'))]
setnames(tmp,'p_pairs','p_pairs_thrsh')
tmp <- merge(po,tmp,by=c('SOURCE'))
write.csv(tmp,file=paste0(outfile.base,'-Prop_sources','.csv'))


## combine ----

g_leg <- ggarrange(g2_leg,g4_leg,widths=c(0.66,0.33))
g <- grid.arrange(arrangeGrob(ggarrange(g1+ theme(legend.title=element_text(size=rel(1))),
                                        ggarrange(g2 + theme(legend.position="none",axis.text.x=element_blank(), axis.ticks.x=element_blank()), g3 + theme(legend.position="none"),
                                                  g4+theme(legend.position="none"),labels=c('B','C','D'),font.label = list(size=18), ncol = 3, align='hv'),
                                        g_leg,
                                        labels=c('A','',''),font.label = list(size=18),ncol=1,heights=c(0.55,0.4,0.05))))

ggsave(g,filename=paste0(outfile.base,'-4panelplot.png'),w=11,h=10)
ggsave(g,filename=paste0(outfile.base,'-4panelplot.pdf'),w=11,h=10)


