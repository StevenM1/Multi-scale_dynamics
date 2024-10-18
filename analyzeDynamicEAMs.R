args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  decisionModel <- args[1]
  learningModel <- args[2]
  dctModel <- args[3]
  task <- args[4]
  n_cosines <- as.numeric(args[5])
  trendModel <- args[6]
  savePlots <- TRUE
} else {
  # manual run
  rm(list=ls())
  decisionModel <- 'RDM'
  learningModel <- c('NULL', 'zSM', 'zSMuAH', 'zSMaAH', 'zSMuAHuFH', 'zSMuAHaFH', 'zSMuAHuFH', 'zSMaAHaFH', 'zSMuAHaRR', 'bV', 'zSMuAHbV', 'zSMuAHbV_with_baseline')[11]
  dctModel <- 'NULL' #'NULL'
  trendModel <- 'NULL' #'B'
  n_cosines = ''
  task <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2', 'mileticvanmaanen2019exp2block1')[1]
  savePlots <- FALSE
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
save_fn_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-', learningModel, '-DCT-', dctModel, '.RData'))
if(trendModel != 'NULL') save_fn_samples <- gsub('-DCT-NULL', paste0('-exptrend-', trendModel), save_fn_samples)
if(n_cosines > 0) save_fn_samples <- gsub('-DCT-', paste0('-DCT-',n_cosines), save_fn_samples)

save_fn_pp <- sub('.RData', '_pp.RData', save_fn_samples)
if(!file.exists(save_fn_samples)) stop('Samples not found')
load(save_fn_samples); nchains <- chain_n(samplers); filter=colnames(nchains)[nchains[1,]>0][sum(nchains[1,]>0)];
if(!savePlots) plot_chains(samplers, filter=filter, selection='mu', layout=c(3,3))
if(!savePlots) plot_chains(samplers, filter=filter, selection='alpha', subject=names(samplers[[1]]$data)[6], layout=c(3,3))
if(!savePlots) plot_chains(samplers, filter=filter, selection='LL', layout=c(3,3))


#plot_density(samplers, selection = 'mu', layout=c(4,4))
#plot_density(samplers, selection = 'alpha', layout=c(4,4))


dat$accuracy <- as.numeric(dat$S) == as.numeric(dat$R)      ## ensure logical
dat$choice <- dat$R == levels(dat$R)[1]
if(!file.exists(save_fn_pp)) {
  pp <- post_predict(samplers, n_cores = 20, filter=filter) #, subfilter=700) ## subfilter=nchains[1,4]-1000)
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
fn_ <- file.path(figures_dir, paste0('model-', decisionModel, '-', learningModel, '_dct-', dctModel, '_XXXX.pdf'))
if(trendModel != 'NULL') fn_ <- gsub('_dct-NULL', paste0('_exptrend-', trendModel), fn_)
if(n_cosines > 0) fn_ <- gsub('_dct-', paste0('_dct-',n_cosines), fn_)

fn_ <- file.path(figures_dir, paste0('model-', decisionModel, '-', learningModel, '_dct-', dctModel, '_XXXX.pdf'))
if(trendModel!='NULL') fn_ <- gsub('_dct-', '_exptrend-', fn_)
makeDefaultPlots(dat=dat, pp=pp, samplers=samplers,
                 save_fn_template=fn_,
                 include_posterror=grepl('AH', learningModel))  # accuracy history effect

# Diagnostics -------------------------------------------------------------
gds = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) round(EMC2:::gd_pmwg(samplers,selection=x,print_summary = FALSE),2))
iats = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) EMC2:::iat_pmwg(samplers,selection=x, print_summary=FALSE))
min_es = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) EMC2:::es_pmwg(samplers,selection=x,summary_alpha=min, print_summary=FALSE))
names(gds) <- names(iats) <- names(min_es) <- c('mu', 'alpha', 'correlation', 'variance')
all <- list('gelman diags'=gds, 'iats'=iats, 'minimum effective sizes'=min_es)

sink(gsub('_XXXX.pdf', '_diagnostics.txt', fn_))
print(all)
sink()
