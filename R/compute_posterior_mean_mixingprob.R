
# pre-amble ----

require(data.table)

# set dir ----
out.dir <- '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam'

args <- list(bground=='unif'
             #bground=='lnorm'
             #bground='gmm'
             )

# load data ----

if(args$bground=='unif'){fit <- readRDS(file.path(out.dir,'mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-fitted_stan_model.rds'))}
if(args$bground=='lnorm'){fit <- readRDS(file.path(out.dir,'mm_bglnorm_piGP_240624-sens_bg_lnorm_MSM-146757/mm_bglnorm_piGP_240624-sens_bg_lnorm_MSM-fitted_stan_model.rds'))}
if(args$bground=='gmm'){fit <- readRDS(file.path(out.dir,'mm_bgGMM_piGP_240711-sens_bg_2DMM_excl_low_res_MSM-315895/mm_bgGMM_piGP_240711-sens_bg_2DMM_excl_low_res_MSM-fitted_stan_model.rds'))}


# summarise omega ----

po <- fit$draws(inc_warmup = FALSE,
                format = 'draws_df',
                variables = 'logit_y_mix_0'
)

po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, omega:= inv.logit(value)]
po <- po[,
         list( q = quantile(omega, probs = c(0.5, 0.025, 0.975) ),
               stat = c('M','CL','CU')
         )
]
po <- dcast.data.table(po, .~stat, value.var = 'q')
po[, L:= paste0(round(M*100,0),'% [',round(CL*100,0),'-',round(CU*100,0),'%]')]

po

