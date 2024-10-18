args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  decisionModel <- args[1]
  learningModel <- args[2]
  trendType <- args[3]
  trendPar <- args[4]
  nTrendPar <- args[5]
  task <- args[6]
  savePlots <- TRUE
} else {
  rm(list=ls())
  decisionModel <- 'RDM'
  learningModel <- c('NULL', 'zSM', 'bV', 'zSMuAHbV', 'zSMaAHbV')[4]
  trendType <- c('NULL', 'DCT', 'poly', 'exptrend') [2]
  trendPar <- c('B', 'v', 'u')[1]
  nTrendPar <- n_cosines <- 10
  task <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2', 'mileticvanmaanen2019exp2block1')[1]
  cores_for_chains <- 3
  cores_per_chain <- 7
  max_trys <- 20
  savePlots <- FALSE
}
if(trendType == 'NULL') {
  trendPar <- 'NULL'
  nTrendPar <- n_cosines <- 0
}

library(EMC2)
library(emcAdapt)
library(Hmisc)
library(stringr)
library(forecast)
source('./utils_plotting.R')


# Load all ----------------------------------------------------------------
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMs')
samples_dir = file.path(root_dir,  'samples')
figures_dir = file.path(root_dir, 'figures')

print(load(file.path(root_dir, 'datasets', paste0(task, '.RData'))))
if(trendType == 'NULL') {
  save_fn_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-', learningModel, '-DCT-NULL.RData'))
} else {
  save_fn_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-', learningModel, '-trend-', trendType, '-', trendPar, '-', nTrendPar, '.RData'))
}
#if(trendModel != 'NULL') save_fn_samples <- gsub('-DCT-NULL', paste0('-exptrend-', trendModel), save_fn_samples)
#if(n_cosines > 0) save_fn_samples <- gsub('-DCT-', paste0('-DCT-',n_cosines), save_fn_samples)

save_fn_pp <- sub('.RData', '_pp.RData', save_fn_samples)
if(!file.exists(save_fn_samples)) stop('Samples not found')
load(save_fn_samples); nchains <- chain_n(samplers); filter=colnames(nchains)[nchains[1,]>0][sum(nchains[1,]>0)];
if(!savePlots) plot_chains(samplers, filter=filter, selection='mu', layout=c(3,3))
if(!savePlots) plot_chains(samplers, filter=filter, selection='alpha', subject=names(samplers[[1]]$data)[6], layout=c(3,3))
if(!savePlots) plot_chains(samplers, filter=filter, selection='LL', layout=c(3,3))

# pairs_posterior(samplers, selection='alpha', subject='6')

dat$accuracy <- as.numeric(dat$S) == as.numeric(dat$R)      ## ensure logical
dat$choice <- dat$R == levels(dat$R)[1]
if(!file.exists(save_fn_pp)) {
  # undebug(EMC2:::update_pars)
  # undebug(EMC2:::make_data)
  # undebug(EMC2:::post_predict)
  pp <- post_predict(samplers, n_cores = 10, filter=filter) #, subfilter=700) ## subfilter=nchains[1,4]-1000)
  pp <- pp[order(pp$subjects,pp$postn,pp$trials),]
  pp$accuracy <- as.numeric(pp$S) == as.numeric(pp$R)       ## ensure logical
  pp$choice <- pp$R==dat[dat$choice==1,'R'][1]
  save(pp, file=save_fn_pp)
} else {
  load(save_fn_pp)
}

# Plots -------------------------------------------------------------------
figures_dir = file.path(figures_dir, paste0('dataset-', task))
if(!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
save_fn_figures <- file.path(figures_dir, paste0('model-', decisionModel, '-', learningModel, '-trend-', trendType, '-', trendPar, '-', nTrendPar, '_XXXX.pdf'))
makeDefaultPlots(dat=dat, pp=pp, samplers=samplers,
                 save_fn_template=save_fn_figures,
                 include_posterror=grepl('AH', learningModel))  # accuracy history effect

# trial_durations <- list('wagenmakers2004_CS'=1.265,
#                          'forstmann2008'=2.850,
#                          'mileticvanmaanen2019exp2block2'=1.627,
#                          'wagenmakers2008exp2'=0.837)
# par(mfrow=c(1,1))
# plotSpectrum(dat, pp, trial_duration=trial_durations[[task]], plot.log=TRUE, ylim=c(-5,-1))
# par(mfrow=c(2,3))
# for(subjects in 1:6) {
#   plotSpectrum(droplevels(dat[dat$subjects==subjects,]), droplevels(pp[pp$subjects==subjects,]), trial_duration=trial_durations[[task]], plot.log=TRUE)
# }
# plot_fit(dat, pp, factors = c('S'))

# Diagnostics -------------------------------------------------------------
gds = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) round(EMC2:::gd_pmwg(samplers,selection=x,print_summary = FALSE),2))
iats = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) EMC2:::iat_pmwg(samplers,selection=x, print_summary=FALSE))
min_es = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) EMC2:::es_pmwg(samplers,selection=x,summary_alpha=min, print_summary=FALSE))
names(gds) <- names(iats) <- names(min_es) <- c('mu', 'alpha', 'correlation', 'variance')
all <- list('gelman diags'=gds, 'iats'=iats, 'minimum effective sizes'=min_es)

sink(gsub('_XXXX.pdf', '_diagnostics.txt', save_fn_figures))
print(all)
sink()
