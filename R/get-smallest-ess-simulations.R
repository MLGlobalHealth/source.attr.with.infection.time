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
  job.name = 'simulations_500truepairs_srcage_exclTE16_subsample',
  stan.model.vanilla = 'mm_sigHierG_bgUnif_piVanilla_220408',
  stan.model.cov = 'mm_sigHierG_bgUnif_piReg_230111',
  stan.model.hsgp = 'mm_bgUnif_piGP_221027'
)

source('R/functions_plotting.R')

tab <- data.table(F=list.dirs('/Users/alexb/Box Sync/Roadmap/source_attribution'))
tab <- tab[ grepl('simulations_500truepairs_srcage_exclTE16_subsample', F) ,]
tab[, model:= 'vanilla']
tab[grepl('piReg',F), model:= 'covariate']
tab[grepl('piGP',F), model:= 'hsgp']
tab[, prop:= as.numeric(gsub('^.*_subsample([0-9]+)pct-([0-9]+)','\\1',F))]
tab[prop==0, prop:= 100]
regex <- '^.*attribution/([A-Za-z0-9_]+)-([A-Za-z0-9_]+)-([0-9]+)'
tab[, outfile.base:= file.path(F,paste0(gsub(regex,'\\1',F),'-',gsub(regex,'\\2',F)))]

po <- NULL
for(i in 1:nrow(tab)){
  tmp <- fread(file=paste0(tab[i,outfile.base],'-convergence','.csv'))
  tmp[, model:= tab[i,model]]
  tmp[, prop:= tab[i,prop]]
  po <- rbind(po,tmp)
}

cat(paste0("smallest ESS = ",floor(po[ ess_bulk==min(ess_bulk), ess_bulk])))
