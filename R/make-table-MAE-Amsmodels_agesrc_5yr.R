## panel MAE ----

require(ggplot2)
require(ggsci)
require(scales)
require(data.table)
require(kableExtra)
require(ggpubr)
require(dplyr)

args <- list(
  indir = '/Users/alexb/Box Sync/Roadmap/source_attribution',
  out.dir = '/Users/alexb/Box Sync/Roadmap/source_attribution/figures',
  #job.name = 'simulations_network_500truepairs_prop_subsample',
  job.name = 'simulations_408truepairs_srcage_agerec_exclTE16_subsample14pct',
  stan.model.vanilla = 'mm_sigHierG_bgUnif_piVanilla_220408',
  stan.model.cov = 'mm_sigHierG_bgUnif_piReg_230111',
  stan.model.hsgp.2d = 'mm_bgUnif_piGP_221027',
  stan.model.hsgp.1da = 'mm_bgUnif_pi1DGP_sim_230224',
  stan.model.hsgp.1db = 'mm_bgUnif_pi1DGP_sim_230224',
  stan.model.hsgp.1dc = 'mm_bgUnif_pi1DGP_sim_230224b'
)

source('R/functions_plotting.R')

tab <- data.table(F=list.dirs('/Users/alexb/Box Sync/Roadmap/source_attribution'))
tab <- tab[ grepl('exclTE16_subsample14pct', F) ,]
tab[, model:= 'vanilla']
tab[grepl('piReg',F), model:= 'covariate']
tab[grepl('piGP',F), model:= 'hsgp2d']
tab[grepl('pi1DGP_sim_230224',F) & grepl('_srcage_',F), model:= 'hsgp1d_src']
tab[grepl('pi1DGP_sim_230224',F) & grepl('_agerec_',F), model:= 'hsgp1d_rec']
tab[grepl('mm_bgUnif_pi1DGP_sim_230224b',F), model:= 'hsgp_srcrec']
tab[, prop:= as.numeric(gsub('^.*_subsample([0-9]+)pct-([0-9]+)','\\1',F))]
tab[prop==0, prop:= 100]
regex <- '^.*attribution/([A-Za-z0-9_]+)-([A-Za-z0-9_]+)-([0-9]+)'
tab[, outfile.base:= file.path(F,paste0(gsub(regex,'\\1',F),'-',gsub(regex,'\\2',F)))]

po <- NULL
src <- NULL
thrsh <- NULL
for(i in 1:nrow(tab)){
  tmp <- readRDS(file=paste0(tab[i,outfile.base],'-MAE_5yrage_all_pairs','.RDS'))
  tmp[, p_pairs:= tab[i,prop]/100]
  tmp[p_pairs==1, p_pairs:= 0.08]

  tmp[, model:= tab[i,model]]

  tmp2 <- readRDS(file = paste0(tab[i,outfile.base],'.rds'))
  tmp3 <- tmp2[, list(N_sources=length(unique(FROM_ID))),by='TO_ID']
  tmp3 <- tmp3[, list(p_pairs=tmp[, p_pairs],
                      N_sources=mean(N_sources))]
  tmp3[, model:= tab[i,model]]

  sim_scenario <- readRDS(file=paste0(tab[i,outfile.base],'.RDS'))
  sim_scenario[, FROM_AGE_FIVE:= cut(FROM_AGE,breaks=seq(15,75,5),include.lowest=T,right=F,
                                     labels=c('15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75'))]

  sim_scenario[, STAGE_thrsh:= 0]
  sim_scenario[GEN_DIST<0.015, STAGE_thrsh:= 1]
  sim_scenario[, p_pairs:= tab[i,prop]/100]
  sim_scenario[p_pairs==1, p_pairs:= 0.08]

  tmp4 <- sim_scenario[, list(N=length(FROM_ID[STAGE_thrsh==1]),
                              N_true=length(FROM_ID[TRANSMISSION_PAIR=='Yes'])), by=c('p_pairs','FROM_AGE_FIVE')]
  tmp4 <- tmp4[, list(FROM_AGE_FIVE=FROM_AGE_FIVE,
                      pct_thrsh=N/nrow(sim_scenario[STAGE_thrsh==1]),
                      pct_true=N_true/nrow(sim_scenario[TRANSMISSION_PAIR=='Yes'])), by=c('p_pairs')]

  tmp4[, thrsh_error:= abs(pct_true - pct_thrsh)]
  tmp4 <- tmp4[, list(MAE=mean(thrsh_error)),by='p_pairs']
  tmp4[, model:= tab[i,model]]

  po <- rbind(po,tmp)
  src <- rbind(src,tmp3)
  thrsh <- rbind(thrsh,tmp4)
}

# save table
po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]
thrsh[, mae_thrsh:= paste0(round(MAE*100,1),'%')]
po <- dcast(po,p_pairs~model,value.var='L')
po <- merge(subset(thrsh,select=c('p_pairs','mae_thrsh'),model=='vanilla'),
            po,
            by=c('p_pairs'))
po <- po[order(-p_pairs)]
po <- melt(po,id.vars='p_pairs')
po[, variable:= factor(variable,levels=c('mae_thrsh','vanilla','covariate','hsgp1d_src','hsgp1d_rec','hsgp_srcrec','hsgp2d'))]
po <- po[order(variable),]
saveRDS(po,file=file.path(args$out.dir,paste0(args$job.name,'-compare_MAE_Amsmodels_ages_source_5yr.RDS')))

