
posterior_prob_pair <- function(model_fit,sim_scenario,dps){
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

  po[, d_TSeqT:= round(TIME_ELAPSED,1)]
  po <- merge(po,dps,by=c('d_TSeqT'),all.x=T)
  # identify % points in signal cone with trsm prob >50%
  po[, incone:=0]
  po[GEN_DIST>q20 & GEN_DIST<q80, incone:=1]
  po[, prob_pair:=0]
  po[M>0.5, prob_pair:=1]
  return(po)
}

make_plot_simulated_data_colour_prob_tpair <- function(po,dps_clock,sim_scenario,outfile.base){

  po[SOURCE_PAIR=='Source category 1 - No', SOURCE_PAIR:='No' ]
  po[SOURCE_PAIR=='Source category 2 - Yes', SOURCE_PAIR:='Yes' ]

  p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
  p.alpha <- 0.7

  p <- ggplot(data = dps_clock,aes(x=d_TSeqT))+
    geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
    theme_bw(base_size=16)+geom_point(data=po,aes(shape=TRANSMISSION_PAIR,x=TIME_ELAPSED,y=GEN_DIST,colour=M))+
    scale_colour_gradient(name="Transmission pair\nprobability under model\n(posterior median)",
                          low = "grey90",
                          high = muted("blue"),
                          breaks = seq(0,1,0.25),
                          n.breaks=4,
                          labels = scales::label_number(accuracy = 0.01),
                          guide=guide_colourbar(direction='horizontal',barwidth=15)) +
    labs(x='Time elapsed (in years)',y='\nPatristic distance of\nphylogenetically possible pair',shape="True transmission pair")+
    theme(legend.position='bottom',
          legend.title=element_text(size=rel(0.7)),
          legend.text=element_text(size=rel(0.7)))+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0),labels = scales::percent, limits = c(0,1)) +
    coord_cartesian(xlim=c(0,15),ylim=c(0,0.2)) +
    guides(shape = guide_legend(nrow = 2,order=1),by.col=T,colour = guide_colorbar(order=2,barwidth=15))
  p

}

make_plot_tpair_violin <- function(po,pal){

  # plot prob of being a pair
  p <- ggplot(po, aes(x = TRANSMISSION_PAIR, y = M)) +
    geom_jitter(aes(colour=cat), width = 0.3, height = 0, alpha = 0.7) +
    geom_violin(fill = 'transparent') +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggsci::scale_color_npg() +
    labs(x = 'Transmission pair', y = '\nTransmission pair probability\n(posterior median)') +
    scale_colour_manual(name = "", values = c(pal[4],pal[5],pal[2],pal[1]),
                        labels=c('Unlinked pair, no signal',
                          'Unlinked pair, false signal',
                          'Transmission pair, with signal',
                          'Transmission pair, no signal')) +
                        #guide = guide_legend(override.aes = list(color = 'white'))) +
    theme_bw(base_size=16) +
    theme( legend.position = 'bottom',
           #legend.text = element_blank(),
           legend.margin=margin(t=-10,r=-20,b=0,l=-20)) +
    guides(col = guide_legend(nrow = 2),by.col=T)
  return(p)

}

plot_estimated_flows_p <- function(tmp,thrsh){

  tmp[, SOURCE:= factor(BIN_COV,levels=c('cat1','cat2'),labels=c('Category\n1','Category\n2'))]
  tmp[, type:= 'Mixture model']

  thrsh <- melt(thrsh,id.vars=c('p_pairs','SOURCE'))
  setnames(thrsh,'value','M')
  thrsh[, type:= factor(variable,levels=c('p_truth','p_thrsh'),labels=c('Simulated data (truth)','1.5% distance threshold'))]

  tmp <- merge(tmp,thrsh,by=c('SOURCE','p_pairs','type','M'),all=T)
  tmp[, type:= factor(type,levels=c('Simulated data (truth)','Mixture model','1.5% distance threshold'))]
  pal3 <- pal_npg("nrc")(1)

  pal1 <- pal_npg("nrc")(4)[1]
  pal2 <- pal_npg("nrc")(4)[2]
  pal4 <- pal_npg("nrc")(4)[4]

  p <- ggplot(tmp) +
    geom_bar(aes(x=as.factor(SOURCE),y=M,fill=as.factor(type)),stat='identity', position = position_dodge(width=0.9),alpha=0.5, width=0.8) +
    geom_errorbar(aes(x=as.factor(SOURCE),y=M,ymin=CL, ymax=CU,color=as.factor(type)), width=0.3, position = position_dodge(width=0.9), size=1) +
    scale_fill_manual(name="",values = c('grey50',pal2,pal4),labels=c('Simulated data (truth)','Mixture model','1.5% distance threshold')) +
    scale_colour_manual(name="",values = c('grey50','grey50','grey50'),labels=c('Simulated data (truth)','Mixture model','1.5% distance threshold')) +
    scale_shape_manual(name="Estimated proportions\ndistance threshold",values = c(17),labels=c('1.5% threshold')) +
    theme_bw(base_size=16) + theme(legend.position = "bottom", legend.margin=margin(t=-10)) +
                                   #axis.text.x=element_text(angle=35,vjust=0.7)) +
    labs(x='\n', y='\nPopulation-level\nsources of infection') +
    coord_cartesian(ylim = c(0,1))+scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1)) +
    guides(colour = "none",
           fill = guide_legend(ncol = 1),
           by.col=T)

  p

}

plot_MAE_compare_models <- function(tmp,thrsh,ylim=0.4){

  setnames(thrsh,'MAE','value')

  tmp[, type:= 'Model estimates']
  thrsh[, type:= '1.5% distance threshold']
  tmp[, type:= factor(type,levels=c('Model estimates','1.5% distance threshold'))]
  tmp[, competing_pairs := round(1/p_pairs,2)]
  thrsh[, competing_pairs := round(1/p_pairs,2)]

  p2 <- ggplot(tmp,aes(x=competing_pairs),fill=pal3) +
    geom_line(aes(y=M,color=model),size=1) +
    geom_errorbar(aes(y=M,ymin=CL, ymax=CU,color=model), width=0,size=1) +
    geom_point(aes(y=M,color=model), size=2) +
    geom_line(data=thrsh,aes(y=value,color='grey50',linetype = "1.5% distance threshold"), colour = 'grey50', size=1.2) +
    geom_point(data=thrsh,aes(y=value,color='grey50'),color='grey50', size=2) +
    scale_linetype_manual(name = "", values = c(2),
                          guide = guide_legend(override.aes = list(color = c('grey50')),title.position='top',nrow = 1)) +
    scale_colour_npg() +
    theme_bw(base_size=16) + theme(legend.position = "right",
                                   legend.margin=margin(t=-10))+
    labs(x='Average number of phylogenetically\npossible sources per new infection', y='Mean absolute error \n(flows from sources)  \n',
         color='Model') +
    scale_y_continuous(breaks=seq(0,ylim,0.02)) +
    coord_cartesian(ylim = c(0,ylim))+
    guides(line = guide_legend(nrow = 1,order=2),
           colour = guide_legend(nrow = 3,title.position='top',order=1),
           by.col=T)
}
