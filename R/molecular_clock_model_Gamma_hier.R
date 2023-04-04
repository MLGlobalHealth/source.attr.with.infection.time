require(data.table)
require(ggplot2)
require(cmdstanr)
require(posterior)
require(bayesplot)

## Preamble----
if (1)
{
  pkg.dir <- '~/Documents/GitHub/source.attr.with.infection.time.fork'
  indir <- file.path(pkg.dir,'data_other')
  outdir <- '~/Documents/GitHub/source.attr.with.infection.time/out_simulated'
  outdir <- '~/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical'
  stan.model <- 'clock_model_gamma_hier_220315'
}
if (0)
{
  pkg.dir <- '~/git/source.attr.with.infection.time'
  indir <- file.path(pkg.dir,'data_other')
  outdir <- '/Users/or105/Library/CloudStorage/Box-Box/OR_Work/2022/2022_source_attribution_with_TSI/molecular_clock'
}

outfile.base <- file.path(outdir, "clock_model_gamma_hier_220315")

## load Belgium data----
infile <- '140921_set7_INFO_TRM.R'
file <- paste(indir, '/',infile, sep = '')
load(file)

# exclude B->A, same as A->B
trm.pol.nA <- subset(trm.pol, withA == FALSE, select = c(d_SeqT, d_TSeqT, BRL, FROM, TO))
setkey(trm.pol.nA, d_TSeqT)
trm.pol.nA[, id := seq_len(nrow(trm.pol.nA))]
trm.pol.nA[, pair_id := as.integer(factor(paste0(FROM,'->',TO)))]

# plot raw data
p <- ggplot(trm.pol, aes(x = d_TSeqT, y = BRL)) +
  geom_point(size = 1.2, alpha = 0.5, data = subset(trm.pol, withA == FALSE & BRL > 0.003)) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,17) ) +
  theme_bw() +
  labs(x = 'time elapsed\n(years)', y = 'genetic distance\n(subst/site)', title = 'epidemiologically confirmed\ntransmission pairs\n') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4), panel.grid.minor = element_line(colour = "grey70", size = 0.2)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x^2)
p
ggsave(file = file.path(outdir,'Genetic_distances_against_time_elapsed_confirmed_pairs.png'), p, w = 6, h = 6)

# plot raw data with individuals
p <- ggplot(trm.pol, aes(x = d_TSeqT, y = BRL)) +
  geom_point(aes(col = paste0(FROM,'->',TO)), size = 1.2, alpha = 0.5, data = subset(trm.pol, withA == FALSE & BRL > 0.003)) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0)) +
  ggsci::scale_color_aaas() +
  coord_trans(ylim = c(0,0.1), xlim = c(0,17) ) +
  theme_bw() +
  labs( x = 'time elapsed\n(years)',
        y = 'genetic distance\n(subst/site)',
        colour = 'pair',
        title = 'epidemiologically confirmed\ntransmission pairs\n'
        ) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4), panel.grid.minor = element_line(colour = "grey70", size = 0.2)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x^2)
p
ggsave(file = file.path(outdir,'Genetic_distances_against_time_elapsed_confirmed_pairs_with_col.png'), p, w = 6, h = 6)


dpr <- data.table(d_TSeqT = seq(0.1,18,0.1))
dpr[, id := seq_len(nrow(dpr))]

cat('---------- fit in stan -----------')

options(mc.cores = parallel::detectCores())
model <- cmdstanr::cmdstan_model(file.path(pkg.dir,'stan_model_files',paste0(stan.model,'.stan')),
                                 force_recompile = TRUE
                                 )

# Fit bayesian model in stan----
stan_data <- list(
  N = nrow(trm.pol.nA),
  x = trm.pol.nA$d_TSeqT,
  y = trm.pol.nA$BRL,
  K = length(unique(trm.pol.nA$pair_id)),
  pair_idx = trm.pol.nA$pair_id,
  N_pr = nrow(dpr),
  x_pr = dpr$d_TSeqT
  )

stan_init <- list(
  log_alpha1 = log(10^(-2.5)),
  log_alpha1_pair = rep(0, stan_data$K),
  log_alpha1_pair_sd = abs((log(10^(-2.5)) - log(10^(-2.1)))/2),
  log_phi = log(5e-3),
  log_phi_pair = rep(0, stan_data$K)
  )

model_fit <- model$sample(
  data = stan_data,
  seed = 42L,
  chains = 4,
  parallel_chains = min(4, parallel::detectCores()),
  refresh = 5e2,
  iter_warmup = 5e2,
  iter_sampling = 2e3,
  save_warmup = TRUE,
  init = list(stan_init, stan_init, stan_init, stan_init)
)

# Save output----
tmp <- paste0(outfile.base,"-stan_fit.rds")
cat("\n Save fitted data to file ", tmp , "\n")
model_fit$save_object(file = tmp)

# Check convergence and divergences----
# save summary of parameters
su <- as.data.table(posterior::summarise_draws(
  model_fit$draws(
    variables = c('log_alpha1','log_alpha1_pair','log_alpha1_pair_sd','log_phi','log_phi_pair','log_phi_pair_sd','lp__'),
    inc_warmup = FALSE)
))
tmp <- paste0(outfile.base,"-convergence.csv")
write.csv(su, file = tmp, row.names = TRUE)
cat("\nmin ess = ", su[,min(ess_bulk)], "max rhat = ", su[,max(rhat)])

cat("\nmarginal summaries:")
print(su)


# save summary of diagnostics
diagnostics <- posterior::summarise_draws(posterior::as_draws_df(model_fit$sampler_diagnostics()), ~quantile(.x, probs = c(0.01, 0.5, 0.99)))
diagnostics <- rbind(diagnostics, c('min_ess_bulk', su[,min(ess_bulk)],0,0))
diagnostics <- rbind(diagnostics, c('max_rhat', su[,max(rhat)],0,0))
tmp <- paste0(outfile.base,"-diagnostics.csv")
write.csv(diagnostics, file = tmp, row.names = TRUE)

#
# trace plots
pd <- model_fit$draws(inc_warmup = TRUE)
fit.target.pars <- c('log_alpha0','log_alpha1','phi')
bayesplot:::color_scheme_set("mix-blue-pink")
p <- bayesplot:::mcmc_trace(pd,
                            pars = c(fit.target.pars,'lp__'),
                            n_warmup = 5e2,
                            facet_args = list(ncol = 1, labeller = label_parsed)
                            )
ggsave(file = paste0(outfile.base,'-traces.pdf'), p, w = 12, h = 16)

# pair plots
p <- bayesplot:::mcmc_pairs(pd,
                            pars = c(fit.target.pars,'logphi','beta','alpha[1]','lp__'),
                            diag_fun = "dens",
                            off_diag_fun = "hex"
                            )
ggsave(file = paste0(outdir,'-pairs.pdf'), p, w = 12, h = 12)


# reshape HMC samples----
pd <- model_fit$draws(inc_warmup = FALSE)
po <- list()
tmp <- pd[,,which(grepl('mu\\[',dimnames(pd)[[3]]))]
po$mu <- unname(apply(tmp[,,], 3, rbind))
tmp <- pd[,,which(grepl('alpha\\[',dimnames(pd)[[3]]))]
po$alpha <- unname(apply(tmp[,,], 3, rbind))
tmp <- pd[,,which(grepl('log_alpha1$',dimnames(pd)[[3]]))]
po$log_alpha1 <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_alpha1_pair\\[',dimnames(pd)[[3]]))]
po$log_alpha1_pair <- unname(apply(tmp[,,], 3, rbind))
tmp <- pd[,,which(grepl('log_alpha1_pair_sd',dimnames(pd)[[3]]))]
po$log_alpha1_pair_sd <- as.vector(tmp)
tmp <- pd[,,which(grepl('beta\\[',dimnames(pd)[[3]]))]
po$beta <- unname(apply(tmp[,,], 3, rbind))
tmp <- pd[,,which(grepl('log_phi$',dimnames(pd)[[3]]))]
po$log_phi <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_phi_pair\\[',dimnames(pd)[[3]]))]
po$log_phi_pair <- unname(apply(tmp[,,], 3, rbind))
tmp <- pd[,,which(grepl('log_phi_pair_sd',dimnames(pd)[[3]]))]
po$log_phi_pair_sd <- as.vector(tmp)
tmp <- pd[,,which(grepl('mu_pr',dimnames(pd)[[3]]))]
po$mu_pr <- unname(apply(tmp[,,], 3, rbind))
tmp <- pd[,,which(grepl('y_pr',dimnames(pd)[[3]]))]
po$y_pr <- unname(apply(tmp[,,], 3, rbind))

# mean evolutionary rates----

cat('\nmean evolutionary rate =')
quantile(exp(po$log_alpha1), prob = c(0.025,0.5,0.975) )

dpo <- data.table(reshape2::melt(po$log_alpha1_pair))
setnames(dpo,1:3,c('iter','pair_id','log_alpha1_pair'))
tmp <- data.table(iter = seq_along(po$log_alpha1), log_alpha1 = po$log_alpha1)
dpo <- merge(dpo, tmp, by = 'iter')
dpo[, er := exp(log_alpha1_pair + log_alpha1)]
dps <- dpo[,
           list(p = quantile(er,prob = c(0.025,0.5,0.975)),
                qlabel = c('CL','M','CU')
                ),
           by = 'pair_id']
dps <- dcast.data.table(dps, pair_id ~ qlabel, value.var = 'p')

# posterior predictive check----
dpo <- data.table(reshape2::melt(po$alpha))
setnames(dpo,1:3,c('iter','id','alpha'))
tmp <- data.table(reshape2::melt(po$beta))
setnames(tmp,1:3,c('iter','id','beta'))
dpo <- merge(dpo, tmp, by = c('iter','id'))
dpo[, y_ppr := rgamma(nrow(dpo), shape = alpha, rate = beta)]
dps <- dpo[,
          list( p = quantile(y_ppr,prob = c(0.025,0.5,0.975)),
                qlabel = c('CL','M','CU')
                ),
          by = c('id')
          ]
dps <- dcast.data.table(dps, id ~ qlabel, value.var = 'p')
dps <- merge(trm.pol.nA, dps, by = 'id')
cat('\nposterior pred check, points outside 95%ppr = ', dps[, mean(BRL < CL | CU < BRL)])

# plot
p.palette <- RColorBrewer::brewer.pal(5,'Oranges')

tmp <- dps[, list(CL = median(CL), M = median(M), CU = median(CU)), by = 'd_TSeqT']
p2 <- ggplot(trm.pol, aes(x = d_TSeqT)) +
  geom_errorbar(data = tmp, aes(ymin = CL, ymax = CU), colour = p.palette[3], width = 0.5,linewidth=0.9) +
  geom_point(aes(y = BRL), colour='grey50', size = 1.2, alpha = 0.5, data = subset(trm.pol, withA == FALSE & BRL > 0.003)) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0)) +
  #coord_trans(ylim = c(0,0.1), xlim = c(0,13) ) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,17) ) +
  theme_bw() +
  labs(x = 'Time elapsed\n(years)', y = 'Genetic distance\n(subst/site)') +
  theme(panel.grid.major = element_line(colour = "grey70", linewidth = 0.4),
        panel.grid.minor = element_line(colour = "grey70", linewidth = 0.2)
        )
p2
ggsave(file = paste0(outfile.base,'-ppr_check.pdf'), p2, w = 6, h = 6)


# get quantiles of distances at prediction points----
#   we first characterise the posterior predicted distances around the
#   average evolutionary rate in the meta-analysis
dpo <- data.table(reshape2::melt(po$y_pr))
setnames(dpo,1:3,c('iter','id','value'))
dpo[, variable := 'y_pr']
dps <- dpo[,
           list(
             p = quantile(value, prob = c(0.025, seq(0.1, .9, 0.1), 0.975)),
             qlabel = paste0('q',c(0.025, seq(0.1, .9, 0.1), 0.975)*100)
             ),
           by = c('variable','id')
           ]
dps <- dcast.data.table(dps, variable + id ~ qlabel, value.var = 'p')
dps <- merge(dpr,dps, by = 'id')


p.alpha <- 0.7
tmp <- subset(dps, variable == 'y_pr')
p <- ggplot(tmp, aes(x = d_TSeqT)) +
  geom_ribbon(data = tmp, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = tmp, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_line(data = tmp, aes(y = q50)) +
  geom_point(data = trm.pol.nA, aes(y = BRL), size = 1.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0),labels = scales::label_percent(accuracy = 1L)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,13) ) +
  theme_bw() +
  labs(x = 'Time elapsed\n(years)', y = 'Genetic distance\n(subst/site)') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4),
        panel.grid.minor = element_line(colour = "grey70", size = 0.2)
        )
ggsave(file = paste0(outfile.base,'-pr_clock.png'), p, w = 6, h = 6)

#   we second characterise the posterior predicted distances around
#   hypothetical evolutionary rates in the inferred meta-population
dpo <- data.table(iter = seq_along(po$log_alpha1),
                  log_alpha1 = po$log_alpha1,
                  log_alpha1_pair_sd = po$log_alpha1_pair_sd,
                  log_phi = po$log_phi,
                  log_phi_pair_sd = po$log_phi_pair_sd,
                  DUMMY = 1
                  )
dpo[, log_alpha1_pair := rnorm(nrow(dpo), 0, log_alpha1_pair_sd)]
dpo[, log_phi_pair := rnorm(nrow(dpo), 0, log_phi_pair_sd)]
dpo[, er := exp(log_alpha1 + log_alpha1_pair)]
dpo[, beta := exp( -(log_phi + log_phi_pair))]

dpr[, DUMMY := 1]
dpo <- merge(dpo, dpr, by = 'DUMMY', allow.cartesian = TRUE)
set(dpr, NULL, 'DUMMY', NULL)
dpo[, y_pr := rgamma(nrow(dpo), shape = er * d_TSeqT * beta, rate = beta)]

dps <- dpo[,
           list(
             p = quantile(y_pr, prob = c(0.025, seq(0.1, .9, 0.1), 0.975)),
             qlabel = paste0('q',c(0.025, seq(0.1, .9, 0.1), 0.975)*100)
           ),
           by = c('d_TSeqT')
           ]
dps <- dcast.data.table(dps, d_TSeqT ~ qlabel, value.var = 'p')
cm_dps <- dps
p <- ggplot(dps, aes(x = d_TSeqT)) +
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_line(data = dps, aes(y = q50)) +
  geom_point(data = trm.pol.nA, aes(y = BRL), size = 1.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0),labels = scales::label_percent(accuracy = 1L)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,13) ) +
  theme_bw() +
  labs(x = 'Time elapsed\n(years)', y = 'Genetic distance\n(subst/site)') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4),
        panel.grid.minor = element_line(colour = "grey70", size = 0.2)
  )
p
ggsave(file = paste0(outfile.base,'-pr_clock_metapop.png'), p, w = 6, h = 6)

# plot raw data with individuals over clock signal
p1 <- ggplot(trm.pol, aes(x = d_TSeqT)) +
  #geom_point(data = subset(trm.pol, withA == FALSE & BRL > 0.003),aes(y = BRL,col = paste0(FROM,'->',TO)), size = 1.2, alpha = 0.5) +
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data = subset(trm.pol.nA, BRL > 0.003),aes(y = BRL,col = paste0(FROM,'->',TO)), size = 1.2, alpha = 0.5) +
  #geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol, withA==FALSE & BRL>0.003)) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0)) +
  ggsci::scale_color_aaas() +
  coord_trans(ylim = c(0,0.1), xlim = c(0,17) ) +
  theme_bw() +
  labs( x = 'Time elapsed\n(years)',
        y = 'Patristic distance\n(subst/site)',
        colour = 'Pair'
        #title = 'epidemiologically confirmed\ntransmission pairs\n'
  ) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4), panel.grid.minor = element_line(colour = "grey70", size = 0.2)) +
  #geom_smooth(aes(y = BRL),method = "lm", se = FALSE, color = "black", formula = y ~ x^2)
geom_line(data = dps,aes(y = q50),color = "black")
p1
ggsave(file = file.path(outdir,'Genetic_distances_against_time_elapsed_confirmed_pairs_with_col_signalcone.png'), p1, w = 6, h = 6)

legend_p1 <- cowplot::get_legend(p1 + theme_bw(base_size=16) + theme(legend.position = "bottom"))

g <- ggarrange(
        ggarrange(p1+labs(x='Time elapsed (years)')+theme_bw(base_size=16)+theme(legend.position='none'),
                  p2 +theme_bw(base_size=16) + labs(x='Time elapsed (years)',y=''),
                  nrow=1,labels='AUTO',align='hv',font.label=list(size=20),hjust=-0.1),
        ggarrange(legend_p1,NULL,nrow=1),nrow=2,heights=c(0.9,0.1))
ggsave(g,filename=file.path(outdir,'distances_signal_ppc_panel.png'),w=11,h=5)
ggsave(g,filename=file.path(outdir,'distances_signal_ppc_panel.pdf'),w=11,h=5)

# summarise # obs outside crIs
dps <- dpo[,
           list( p = quantile(y_pr,prob = c(0.025,0.5,0.975)),
                 qlabel = c('CL','M','CU')
           ),
           by = c('id')
]
dps <- dcast.data.table(dps, id ~ qlabel, value.var = 'p')
dps <- merge(trm.pol.nA, dps, by = 'id')
cat('\nposterior pred check, points outside 95%ppr = ', dps[, mean(BRL < CL | CU < BRL)])

## predict distances using posterior medians only ----
dpo <- data.table(iter = seq_along(po$log_alpha1),
                  log_alpha1 = median(po$log_alpha1),
                  log_alpha1_pair_sd = median(po$log_alpha1_pair_sd),
                  log_phi = median(po$log_phi),
                  log_phi_pair_sd = median(po$log_phi_pair_sd),
                  DUMMY = 1
)
dpo[, log_alpha1_pair := rnorm(nrow(dpo), 0, log_alpha1_pair_sd)]
dpo[, log_phi_pair := rnorm(nrow(dpo), 0, log_phi_pair_sd)]
dpo[, er := exp(log_alpha1 + log_alpha1_pair)]
dpo[, beta := exp( -(log_phi + log_phi_pair))]

dpr[, DUMMY := 1]
dpo <- merge(dpo, dpr, by = 'DUMMY', allow.cartesian = TRUE)
set(dpr, NULL, 'DUMMY', NULL)
dpo[, y_pr := rgamma(nrow(dpo), shape = er * d_TSeqT * beta, rate = beta)]

dps <- dpo[,
           list(
             p = quantile(y_pr, prob = c(0.025, seq(0.1, .9, 0.1), 0.975)),
             qlabel = paste0('q',c(0.025, seq(0.1, .9, 0.1), 0.975)*100)
           ),
           by = c('d_TSeqT')
]
dps <- dcast.data.table(dps, d_TSeqT ~ qlabel, value.var = 'p')
dps_new <- dps

p <- ggplot(dps_new, aes(x = d_TSeqT)) +
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_line(data = dps, aes(y = q50)) +
  geom_point(data = trm.pol.nA, aes(y = BRL), size = 1.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0),labels = scales::label_percent(accuracy = 1L)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,13) ) +
  theme_bw() +
  labs(x = 'Time elapsed\n(years)', y = 'Genetic distance\n(subst/site)') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4),
        panel.grid.minor = element_line(colour = "grey70", size = 0.2)
  )
ggsave(file = paste0(outfile.base,'-pr_clock_metapop_posterior_median.png'), p, w = 6, h = 6)

colnames(cm_dps)[2:length(colnames(cm_dps))] <- paste0('cm_',colnames(cm_dps)[2:length(colnames(cm_dps))])
dps <- merge(dps_new,cm_dps,by=c('d_TSeqT'),all=T)

tmp2 <- melt(dps,id.vars=c('d_TSeqT'))
tmp2[,model:='Median']
tmp2[grepl('cm_',variable),model:='Samples']
tmp2[,quantile:=gsub('cm_q','',variable)]
tmp2[,quantile:=gsub('*q','',quantile)]
tmp2 <- dcast(tmp2,d_TSeqT+quantile~model,value.var='value')

p <- ggplot(tmp2, aes(x = d_TSeqT)) +
  geom_point(data = tmp2, aes(x = Samples, y = Median), col = p.palette[4], alpha = p.alpha) +
  geom_line(data = tmp2, aes(x = Samples, y = Samples), col = 'black', alpha = p.alpha) +
  facet_wrap(.~quantile) +
  #geom_line(data = tmp, aes(y = q50)) +
  #geom_point(data = trm.pol.nA, aes(y = BRL), size = 1.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0),labels = scales::label_percent(accuracy = 1L)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0),labels = scales::label_percent(accuracy = 1L)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,0.1) ) +
  theme_bw() +
  labs(x = 'Quantiles Genetic distance\n(subst/site)\n clock model samples', y = 'Quantiles Genetic distance\n(subst/site)\n using posterior median') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4),
        panel.grid.minor = element_line(colour = "grey70", size = 0.2)
  )
ggsave(file = paste0(outfile.base,'-QQ_clock_facet_line.png'), p, w = 9, h = 9)


