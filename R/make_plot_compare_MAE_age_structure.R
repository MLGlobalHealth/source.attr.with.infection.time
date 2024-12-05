
require(ggplot2)
require(ggpubr)
require(gridExtra)
args <- list(
  out.dir = '/Users/alexb/Box Sync/Roadmap/source_attribution'
)

outfile.base <- file.path(args$out.dir,"mm_bgUnif_piGP_221027-simulations_500truepairs_srcage_exclTE16_subsample50pct-1070774/mm_bgUnif_piGP_221027-simulations_500truepairs_srcage_exclTE16_subsample50pct")
tmp1 <- readRDS(file = paste0(outfile.base,'.rds'))
mae1 <- readRDS(file=paste0(outfile.base,'-MAE_GQ_all_pairs','.RDS'))

outfile.base <- file.path(args$out.dir,"mm_bgUnif_piGP_221027-simulations_500truepairs_srcage_lesscorrages_exclTE16_subsample50pct-560476/mm_bgUnif_piGP_221027-simulations_500truepairs_srcage_lesscorrages_exclTE16_subsample50pct")
tmp2 <- readRDS(file = paste0(outfile.base,'.rds'))
mae2 <- readRDS(file=paste0(outfile.base,'-MAE_GQ_all_pairs','.RDS'))

# plot ages ----
tmp1[, TRSM:= 'Non-transmission pair']
tmp1[TRANSMISSION_PAIR=='Yes', TRSM:= 'True transmission pair']

pal <- pal_npg("nrc")(2)
g1 <- ggplot(tmp1) +
  geom_point(aes(x=TO_AGE,y=FROM_AGE,colour=TRANSMISSION_PAIR),alpha=0.7) +
  labs(x='Age of recipient on estimated\ninfection date of recipient',
       y='Age of source on estimated\ninfection date of recipient',
       colour = 'Transmission pair') +
  scale_x_continuous(breaks=seq(15,80,10),labels=seq(15,80,10)) +
  scale_y_continuous(breaks=seq(15,80,10),labels=seq(15,80,10)) +
  scale_colour_manual(values=c(pal[1],pal[2])) +
  theme_bw(base_size=16) +
  theme(strip.background=element_blank(),
        legend.position='bottom',
        legend.margin=margin(t=-10)) +
  guides(colour = guide_legend(nrow = 2, title.position='top'))
g1

tmp2[, TRSM:= 'Non-transmission pair']
tmp2[TRANSMISSION_PAIR=='Yes', TRSM:= 'True transmission pair']

pal <- pal_npg("nrc")(2)
g2 <- ggplot(tmp2) +
  geom_point(aes(x=TO_AGE,y=FROM_AGE,colour=TRANSMISSION_PAIR),alpha=0.7) +
  labs(x='Age of recipient on estimated\ninfection date of recipient',
       y='Age of source on estimated\ninfection date of recipient',
       colour = 'Transmission pair') +
  scale_x_continuous(breaks=seq(15,80,10),labels=seq(15,80,10)) +
  scale_y_continuous(breaks=seq(15,80,10),labels=seq(15,80,10)) +
  scale_colour_manual(values=c(pal[1],pal[2])) +
  theme_bw(base_size=16) +
  theme(strip.background=element_blank(),
        legend.position='bottom',
        legend.margin=margin(t=-10)) +
  guides(colour = guide_legend(nrow = 2, title.position='top'))
g2

mae1[, cat:= 'Distinguished from\nunlinked pairs']
mae2[, cat:= 'More similar to\nunlinked pairs']

mae <- merge(mae1,mae2,by=c('CL','CU','IL','IU','M','cat'),all=T)

pal4 <- pal_npg("nrc")(4)[c(2,4)]

g3 <- ggplot(mae,aes(x=cat),fill=pal3) +
  geom_bar(aes(y=M,fill=cat),stat='identity', position = position_dodge(width=0.9),alpha=0.5,width=0.8) +
  geom_errorbar(aes(y=M,ymin=CL, ymax=CU),color='grey50', width=0.3, position = position_dodge(width=0.9),size=1) +
  scale_fill_manual(name="",values = c(pal4)) +
  theme_bw(base_size=16) + theme(legend.position = "none",
                                 legend.text = element_blank())+
  scale_y_continuous(breaks=c(seq(0,0.1,0.02))) +
  labs(x='\nAge structure of actual transmission pairs', y='\nMean absolute error \n(infections attributed to sources)') +
  coord_cartesian(ylim = c(0,0.1))#scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
  g3

g <- grid.arrange(arrangeGrob(ggarrange(g1,labels='A',font.label = list(size=18)), ggarrange(g2,labels='B',font.label = list(size=18)), ncol = 1),
                              ggarrange(g3,labels='C',font.label = list(size=18)),# Second row with 2 plots in 2 different columns
             ncol = 2)
ggsave(g,file=paste0(outfile.base,'-compare_ages_source_recipients_MAE_bluered','.pdf'),w=12,h=10)
