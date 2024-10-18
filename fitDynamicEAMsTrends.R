#### Dynamic
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  #### Run from command line
  decisionModel <- args[1]
  learningModel <- args[2]
  task <- args[3]
  trendType <- args[4]
  trendPar <- args[5]
  nTrendPar <- args[6]
  cores_for_chains = args[7]
  cores_per_chain = args[8]
  if(length(args) == 9) max_trys <- as.integer(args[9]) else max_trys <- 15
#  if(length(args) == 9) n_cosines <- as.integer(args[9]) else n_cosines <- 0
  debugMode <- FALSE
} else {
  #### Run interactively
  rm(list=ls())
  decisionModel <- 'RDM'
  learningModel <- c('NULL', 'zSM', 'bV', 'zSMuAHbV', 'zSMaAHbV')[1]
  trendType <- c('NULL', 'DCT', 'poly', 'exptrend') [2]
  trendPar <- c('B', 'v', 'u')[3]
  nTrendPar <- n_cosines <- 10
  task <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2', 'mileticvanmaanen2019exp2block1')[1]
  cores_for_chains <- 3
  cores_per_chain <- 7
  max_trys <- 20
  debugMode <- FALSE
}
if(trendType == 'NULL') {
  trendPar <- 'NULL'
  nTrendPar <- n_cosines <- 0
}

init_chains <- function(samplers, start_mu = NULL, start_var = NULL,
                        verbose = FALSE, particles = 1000,
                        cores_per_chain=1,cores_for_chains = length(samplers))
{
  # save attributes
  model_list <- attr(samplers, 'model_list')
  design_list <- attr(samplers, 'design_list')
  data_list <- attr(samplers, 'data_list')

  samplers <- parallel::mclapply(samplers, EMC2:::init, start_mu = start_mu, start_var = start_var,
                                 verbose = verbose, particles = particles,
                                 n_cores = cores_per_chain, mc.cores=cores_for_chains)

  attr(samplers, 'model_list') <- model_list
  attr(samplers, 'design_list') <- design_list
  attr(samplers, 'data_list') <- data_list
  samplers
}


# Set-up directories
print(paste('Running:', decisionModel, learningModel, task, trendType, trendPar, nTrendPar))
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMs')
samples_dir = file.path(root_dir, 'samples')
figures_dir = file.path(root_dir, 'figures')

## Libraries
library(EMC2)
library(emcAdapt)
library(abind)

# File names
save_fn_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-', learningModel, '-trend-', trendType, '-', trendPar, '-', nTrendPar, '.RData'))

# save_fn_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-', learningModel, '-DCT-', dctModel, '.RData'))
# if(trendModel != 'NULL') save_fn_samples <- gsub('-DCT-NULL', paste0('-exptrend-', trendModel), save_fn_samples)
# if(n_cosines > 0) save_fn_samples <- gsub('-DCT-', paste0('-DCT-',n_cosines), save_fn_samples)

save_fn_h5 <- gsub('\\.RData', '.h5', save_fn_samples)
doSample <- TRUE

# if(task %in% c('mileticvanmaanen2019exp2block2', 'forstmann2008') & (dctModel != 'NULL')) {
#   # hard models on hard datasets, search a bit wider and slower
#   p_accept = 0.7
#   particle_factor = 60
# } else {
# Sampling parameters -- these are the default values
p_accept = 0.8
particle_factor = 40
# }
print(paste0('Using p_accept=', p_accept, ' and particle_factor=', particle_factor))
preburn = 150
iter <- 2000

# tux18 has an issue with h5 files
save_fn_h5 <- ifelse(Sys.info()[['nodename']] == 'tux18psy', save_fn_samples, save_fn_h5)


# Check for existing samples
if(file.exists(save_fn_samples)) {
#  stop('Samples already exit.')
  load(save_fn_samples)
  niter <- chain_n(samplers)[1,4]
  if(niter > 0) {
    if(niter >= 1500) stop('Done - not converged, but stopping anyway')
    if(niter < 1500) {print(niter); max_trys = round((1500-niter)/100); print(paste0('new max trys = ',max_trys))}
    # progress <- EMC2:::check_progress(samplers, stage='sample', iter=iter, max_gd=1.1, mean_gd=1.05, min_es=0, min_unique=600, max_trys=max_trys, step_size=100, n_cores=cores_per_chain, verbose=FALSE)
    # progress$iters_total <- niter
    # progress$trys <- niter/progress$step_size
    # progress <- EMC2:::check_progress(samplers, stage='sample', iter=iter, max_gd=1.1, mean_gd=1.05, min_es=0, min_unique=600, max_trys=max_trys, step_size=100, n_cores=cores_per_chain, verbose=TRUE, progress=progress)
    # doSample <- !progress$done
    #
    # ##
    # if(progress$done & progress$trys >= max_trys) {
    #   gd <- EMC2:::check_gd(samplers, stage='sample', max_gd=1.1, mean_gd=1.05, trys=max_trys, verbose=F)
    #   if(!gd$gd_done) {
    #     ## maximum iterations reached, but not converged: try increasing particle_factor and decreasing p_star?
    #     p_accept = 0.8
    #     particle_factor = 40
    #     max_trys = 10
    #     doSample = TRUE
    #     print(paste0('Maximum iterations reached with default settings;
    #                  trying ', max_trys, ' more iterations with more particles
    #                  (particle_factor=', particle_factor, ')
    #                  and lower acceptance rate (p_accept=', p_accept, ')'))
    #   }
    # }
    # rm(progress)
  }
  rm(samplers)  # free memory!
} else {
  source(file.path(root_dir, 'makeSamplersTrends.R'))
  makeSamplers(decisionModel=decisionModel, learningModel=learningModel,
               trendType=trendType, trendPar=trendPar, nTrendPar=nTrendPar, task=task, samples_dir=samples_dir, n_chains=3, trend_nuisance=TRUE)
}

if(debugMode) {
  debug(EMC2:::update_pars)
  debug(EMC2:::get_pars)

  #debug(run_samplers)
  run_emc(samplers=save_fn_samples,
          fileName=save_fn_h5, iter = 100,
          cores_for_chains=1,
          cores_per_chain=1, verbose = TRUE, verboseProgress = TRUE,
          max_trys=max_trys)
  doSample <- FALSE
}



if(doSample) {
  # Check for init
  if(learningModel != 'NULL') {
    samplers <- EMC2:::loadRData(save_fn_samples)
    if(!samplers[[1]]$init & trendType=='NULL') {  # NB: nuisance is an issue at the moment
      # load corresponding null model
      save_fn_null_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-NULL-DCT-NULL.RData'))
      if(!file.exists(save_fn_null_samples)) {
        warning('Could not find any samples corresponding to a null model -- proceeding with start points initialized at the default values. May be slower to converge.')
      } else {
        samplers_null <- EMC2:::loadRData(save_fn_null_samples)

        # find mu_mean and mu_var
        idx <- samplers_null[[1]]$samples$stage == 'sample'
        if((sum(idx) <= 1)) {
          warning('Null model doesnt have samples at the sample stage yet -- proceeding with start points initialized at the default values. May be slower to converge.')
        } else {
          mu_samples <- do.call(cbind, lapply(samplers_null, function(x) x$samples$theta_mu[,idx]))
          var_samples <- abind(lapply(samplers_null, function(x) x$samples$theta_var[,,idx]), along=3)
          mu_mean <- apply(mu_samples, 1, mean)
          mu_var <- apply(var_samples, 1:2, mean)

          # which parameters are in the new sampler, that are also in the old sampler?
          par_names_new <- samplers[[1]]$par_names # %in% samplers_null[[1]]$par_names
          par_names_null <- samplers_null[[1]]$par_names
          select_pars <- par_names_new[par_names_new%in%par_names_null]
          new_pars <- par_names_new[!par_names_new%in%par_names_null]

          # fill in 0 or 1 for new parameters
          start_mu <- setNames(rep(0, length(par_names_new)), par_names_new)
          start_mu[select_pars] <- mu_mean[select_pars]

          # start with identity, fill in parameters that are known
          start_var <- diag(length(par_names_new))
          dimnames(start_var) <-  list(par_names_new, par_names_new)
          start_var[select_pars, select_pars] <- mu_var[select_pars, select_pars]

          # init
          samplers <- init_chains(samplers, start_mu=start_mu, start_var=start_var,
                                  particles = 1000,
                                  cores_for_chains=cores_for_chains,
                                  cores_per_chain=cores_per_chain)

          save(samplers, file=save_fn_samples)
        }
      }
    }
  }

  samplers <- run_emc(samplers=save_fn_samples, fileName=save_fn_h5,
                      iter = iter,
                      cores_for_chains=cores_for_chains,
                      cores_per_chain=cores_per_chain,
                      preburn=preburn,
                      verbose = TRUE, verboseProgress = TRUE,
                      p_accept = p_accept,
                      particle_factor = particle_factor,
                      max_trys=max_trys)

  # reload samplers
  samplers <- EMC2:::loadRData(save_fn_samples)

  ## Post predicts
  pp_fn = sub('.RData', '_pp.RData', save_fn_samples)

  load(file.path(root_dir, 'datasets', paste0(task, '.RData')))
  dat$stim <- as.numeric(dat$S)
  dat$stim[dat$stim==2] <- -1
  dat$accuracy <- as.numeric(dat$S) == as.numeric(dat$R)   ## ensure logical
  dat$choice <- dat$R == levels(dat$R)[1]

  if(file.exists(pp_fn) & !doSample) {
    load(pp_fn)
  } else {
    nchains <- chain_n(samplers)
    pp <- post_predict(samplers, subfilter=nchains[1,4]-1000, n_cores = 20)
    pp <- pp[order(pp$subjects,pp$postn,pp$trials),]
    pp$accuracy <- as.numeric(pp$S) == as.numeric(pp$R)      ## ensure logical
    pp$choice <- pp$R==dat[dat$choice==1,'R'][1]

    save(pp, file=pp_fn)
  }

  ## make plots ##
  source(file.path(root_dir, 'utils_plotting.R'))

  figure_dir = file.path(figures_dir, paste0('dataset-', task))
  if(!dir.exists(figure_dir)) dir.create(figure_dir, recursive=TRUE)

  save_fn_figures <- file.path(figure_dir, paste0('model-', decisionModel, '-', learningModel, '-trend-', trendType, '-', trendPar, '-', nTrendPar, '_XXXX.pdf'))
  makeDefaultPlots(dat=dat, pp=pp, samplers=samplers, save_fn_template=save_fn_figures)

  ## run diagnostics ##
  gds = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) round(EMC2:::gd_pmwg(samplers,selection=x,print_summary = FALSE), 2))
  iats = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) EMC2:::iat_pmwg(samplers,selection=x, print_summary=FALSE))
  min_es = lapply(c('mu', 'alpha', 'correlation', 'variance'), function(x) EMC2:::es_pmwg(samplers,selection=x,summary_alpha=min, print_summary=FALSE))
  names(gds) <- names(iats) <- names(min_es) <- c('mu', 'alpha', 'correlation', 'variance')
  all <- list('gelman diags'=gds, 'iats'=iats, 'minimum effective sizes'=min_es)
  sink(gsub("XXXX.pdf", 'diagnostics.txt', save_fn_figures))
  print(all)
  sink()
}



