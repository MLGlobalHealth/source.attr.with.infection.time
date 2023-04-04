library(data.table)
library(tidyverse)
library(tidytree)
library(phyloscannerR)
library(ggtree)
library(pals)
library(glue)
library(ggsci)

indir <- "~/Box Sync/Roadmap/analysis_220713/trees/"
outdir <- '/Users/alexb/Box Sync/Roadmap/source_attribution/figures'
trsm <- 'MSM'
rdas <- list.files(indir,pattern = "*AmsMSM__workspace.rda")

# load metadata
pairs.dir <- file.path('~/Documents/GitHub/source.attr.with.infection.time.fork','out_Amsterdam','agegps_updated_criteria_MSM-2010_2022')
all_pairs <- readRDS(file=file.path(pairs.dir, 'all_pairs.rds'))
pairs <- readRDS(file.path(pairs.dir,paste0(trsm,"_pairs.rds")))

for(i in 1:length(rdas)){
  #i <- 5 # bc2
 a.rda <- rdas[i]
  load(file.path(indir,a.rda))

  subtype <- str_match(a.rda, ".*_subtype_([A-Za-z0-9]+)+_wOutgroup_*")[,2]
  print(subtype)
  ptree <- phyloscanner.trees[[1]]

  tree <- ptree$tree

  read.counts <- ptree$read.counts#[1:6487]

  attr(tree, 'BLACKLISTED') <- c(is.na(read.counts), rep(FALSE, tree$Nnode))

  read.counts <- c(read.counts, rep(1, tree$Nnode))
  read.counts[which(is.na(read.counts))] <- 1

  attr(tree, 'READ_COUNT') <- read.counts

  # load metadata
  d1 <- data.table(tip=tree$tip.label)
  regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9]+)-([0-9-]+)-([0-9-]+)_([0-9]+)$'
  d1[, ID:= as.numeric(gsub(regex.tip.label,'\\3',tip))]
  d1[, TRSM:= gsub(regex.tip.label,'\\1',tip)]
  d1 <- merge(d1,unique(subset(pairs,select=c('TO_SEQUENCE_ID','TO_AGE','TO_AGE_GP','TO_EST_INFECTION_DATE'))),by.x='ID',by.y='TO_SEQUENCE_ID',all.x=T)
  d1 <- merge(d1,unique(subset(pairs,select=c('FROM_SEQUENCE_ID','FROM_AGE','FROM_AGE_GP','FROM_EST_INFECTION_DATE'))),by.x='ID',by.y='FROM_SEQUENCE_ID',all.x=T)
  d1 <- merge(d1,unique(subset(all_pairs,select=c('TO_SEQUENCE_ID','TO_CLUSTER_NUMBER'))),by.x='ID',by.y='TO_SEQUENCE_ID',all.x=T)
  d1 <- merge(d1,unique(subset(all_pairs,select=c('FROM_SEQUENCE_ID','FROM_CLUSTER_NUMBER'))),by.x='ID',by.y='FROM_SEQUENCE_ID',all.x=T)

  # keep first obs for sources that have multiple recipients
  setkey(d1,tip)
  d1 <- d1[,.SD[1],by = tip]

  d1[, prs:= 0]
  d1[TRSM=='AmsMSM', prs := 1]
  d1[ID %in% pairs$FROM_SEQUENCE_ID, prs:= 2]
  d1[ID %in% pairs$TO_SEQUENCE_ID, prs:= 3]

  # flag the individuals who are possible sources to the recipients
  tmp <- unique(subset(d1,select=c('TO_CLUSTER_NUMBER'),prs==3))
  tmp[, poss_src:= 1]
  setnames(tmp,'TO_CLUSTER_NUMBER','FROM_CLUSTER_NUMBER')
  d1 <- merge(d1,tmp,by='FROM_CLUSTER_NUMBER',all.x=T)
  d1[is.na(poss_src), poss_src:= 0]

  x <- full_join(as_tibble(tree), d1, by = c("label" = "tip"))
  tree2 <- as.treedata(x)

  new.branch.colours <- attr(tree, "BRANCH_COLOURS")
  new.branch.colours <- as.character(new.branch.colours)
  new.branch.colours[new.branch.colours!='AmsMSM'] <- 'non-MSM'
  new.branch.colours <- factor(new.branch.colours, levels = c("AmsMSM",
                                                              "non-MSM"),
                               labels=c('Amsterdam MSM subgraphs','Amsterdam non-MSM or\noutside Amsterdam'))

  new.individual <- tree2@data$prs
  new.individual <- factor(new.individual, levels=c(0,1,2,3),
                           labels=c('Non-Amsterdam MSM','Amsterdam MSM',
                                    'Amsterdam MSM - potential source' ,'Amsterdam MSM - infected 2010-2021'))
  border <- tree2@data$poss_src
  border <- factor(border, levels=c(0,1),
                           labels=c('Not a recipient or possible source','Recipient or possible source'))

    huh2 <- c(pal_npg("nrc")(1),'grey50')
    pal <- pal_npg("nrc")(9)
    huh2 <- c('grey50',pal[4],pal[2],pal[1])
    pal_branch <- c(pal[3],'grey50')

    tips <- tree$tip.label %>% length()
      tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.3) +
      geom_tippoint(aes(color=new.individual), size=0.75) +

      geom_treescale(x=-0.01, y = length(tree$tip.label)*0.5,  offset=4) +
      scale_colour_manual(values = c("Amsterdam MSM subgraphs" = pal[3],
                                     "Amsterdam MSM - infected 2010-2021"= pal[1],
                                     "Amsterdam MSM - potential source" = pal[2],
                                     "Amsterdam MSM" = pal[4],
                                     "Non-Amsterdam MSM" = 'grey50',
                                     'NA' = 'grey50'),
                          labels=c("Amsterdam MSM subgraphs" = 'Amsterdam MSM subgraphs',
                                   "Amsterdam MSM - infected 2010-2021" = 'Amsterdam MSM infected in 2010-2021',
                                   "Amsterdam MSM - potential source" = 'Plausible source of an Amsterdam MSM infection in 2010-2021',
                                   "Amsterdam MSM" = 'Other Amsterdam MSM',
                                   "Non-Amsterdam MSM" = 'Other individuals from outside Amsterdam or a different Amsterdam risk group',
                                   'NA' = ''),
                          na.value = "black", drop = F, name = "") +
                          theme(legend.text = element_text(size=11),legend.position='right') +
      guides(col=guide_legend(ncol=1)) +
      ylim(-1, length(tree$tip.label) +1)

    ggsave(file.path(outdir,glue("type_{subtype}_tree_newlabels_{trsm}_coloursubgraphs_tips.pdf")), width = 10, height = tips*0.02, limitsize = F)

}
