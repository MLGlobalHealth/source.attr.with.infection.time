
require(data.table)


args <- list(
  source_dir= '~/Documents/GitHub/source.attr.with.infection.time',
  analysis = 'analysis_220713',
  indir = '/Users/alexb/Documents/Roadmap/refactor_code',
  data='~/Box\ Sync/Roadmap/RQ1 Estimating introductions'
)

cat("\nRead phylogenetic subgraphs \n")

indir.phsc <- file.path(args$indir, args$analysis, 'subgraphs')
infile.meta <- file.path(args$indir, args$analysis, 'misc', '220713_sequence_labels.rda')
infile.seq <-	file.path(args$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
infile.georeg <- file.path('/Users/alexb/Box Sync/Roadmap','misc','NEWGEO_220713.csv')

cat('\n Extract subgraph taxa\n')
#dsubgraphtaxa <- extract_subgraphs(indir.phsc)
load(indir.phsc)

cat("\nRead patient metadata \n")
load(infile.seq)
load(infile.meta)
setnames(dind, gsub("CENS_D", "RECART_D", names(dind)))
dind <- unique(dind)

# flag which diagnosed individuals have a sequence
dind$SEQ <- dind$PATIENT %in% ds$PATIENT
dind <- merge(dind,unique(subset(dsubgraphtaxa,select=c(ID,ST))),by.x='PATIENT',by.y='ID',all.x=T)

# add meta data from pre-processed persons file
dind <- as.data.table(dind)
setnames(dind, c('PATIENT','BIRTH_Y','BIRTH_CNTRY'), c('ID','BIRTH_YEAR','BIRTH_COUNTRY'))
tmp <- subset(dind, select=c(ID,BIRTH_YEAR,BIRTH_COUNTRY,ORIGIN,LOC_BIRTH,CITY,TRANSM,GENDER,RECART_D,HIV1_POS_D,HIV1_POS_D_lower,HIV1_POS_D_upper))

dsubgraphtaxa <- merge(dsubgraphtaxa,tmp,by='ID')
setnames(dsubgraphtaxa, gsub("CENS_D", "RECART_DATE", names(dsubgraphtaxa)))

tmp <- dsubgraphtaxa
tmp <- subset(tmp,REP=='000' & SELECT!='Ams')

regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z-]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9-]+)_([0-9]+)$'
tmp[, SEQID:= gsub(regex.tip.label,'\\2',TAXA)]
tmp[, SAMPLING_DATE:= as.Date(gsub(regex.tip.label,'\\9',TAXA),format="%Y-%m-%d")]

# define migrant groups
geo <- data.table(read.csv(infile.georeg))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))
tmp <- merge(tmp,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
## msm
tmp[TRANSM=='MSM', mwmb:="Other"]
tmp[TRANSM=='MSM' & ORIGIN %in% c("NL"), mwmb:="NL"]
tmp[TRANSM=='MSM' & WRLD_born %in% c("WEurope","NorthAm","Oceania") & ORIGIN!='NL', mwmb:="G1"]# western countries (non-NL)
tmp[TRANSM=='MSM' & WRLD_born %in% c("EEurope", "CEurope"), mwmb:="G2"]# eastern and central europe
tmp[TRANSM=='MSM' & WRLD_born %in% c("LaAmCar","DutchCarSuriname"), mwmb:="G3"]# caribbean and south america

## hsx
tmp[TRANSM=='HSX', mwmb:="Other"]
tmp[TRANSM=='HSX' & ORIGIN %in% c("NL"), mwmb:="NL"]
tmp[TRANSM=='HSX' & WRLD_born %in% c("Africa"), mwmb:="G4"]# sub-saharan africa
tmp[TRANSM=='HSX' & WRLD_born %in% c("LaAmCar","DutchCarSuriname"), mwmb:="G5"]# caribbean and south america

meta <- subset(tmp, select=c('ID','ST','ST_CLADE','NAME','FULL_NAME','SELECT','CITY','BIRTH_YEAR','GENDER','TRANSM','HIV1_POS_D','RECART_D','SEQID','SAMPLING_DATE','ORIGIN','BIRTH_COUNTRY','mwmb'))
setnames(meta,c('TRANSM','FULL_NAME','RECART_D','CITY','ST','mwmb'),c('RISK_GROUP','CLUSTER_NUMBER','TREATMENT_DATE','LOCATION','SUBTYPE','BPLACE_MG'))
meta[, bplace:= 'migrant']
meta[ORIGIN=='NL', bplace:= 'non-migrant']
saveRDS(meta,file=file.path('data_Ams',args$analysis,'meta_data_mg_country.rds'))

age <- subset(tmp, select=c('ID','BIRTH_YEAR'))
age[, BIRTH_YEAR_DEC:=BIRTH_YEAR]
saveRDS(age,file=file.path('data_Ams',args$analysis,'data_age.rds'))

# make infection date dataset
if(args$analysis=='analysis_220713'){
  dt <- data.table(read.csv(file.path('data_Ams',args$analysis,'roadmap_cd4_v3_est.csv')))
} else{
  dt <- data.table(read.csv(file.path('data_Ams',args$analysis,'roadmap_cd4_vl_est.csv')))
}
dt <- unique(subset(dt,select=c('id','hiv_pos_d','estsctodiagMedian','estsctodiagLL','estsctodiagUL')))
write.csv(dt,file.path('data_Ams',args$analysis,'Infection_date_est_rec.csv'))
