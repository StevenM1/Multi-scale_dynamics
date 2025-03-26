# get start points from pervious sampler
get_start_mu_var <- function(samplers_null, par_names_new, prior_var) {
  library(abind)
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
    # par_names_new <- samplers[[1]]$par_names # %in% samplers_null[[1]]$par_names
    par_names_null <- samplers_null[[1]]$par_names
    select_pars <- par_names_new[par_names_new%in%par_names_null]
    new_pars <- par_names_new[!par_names_new%in%par_names_null]

    # fill in 0 or 1 for new parameters
    start_mu <- setNames(rep(0, length(par_names_new)), par_names_new)
    start_mu[select_pars] <- mu_mean[select_pars]
    start_mu[grepl('alpha', names(start_mu))] <- qnorm(.5)  # start at 0.5 for learning rates

    # These are just some guestimates based on earlier experience.
    # Note that these are *NOT* priors and have no influence on the likelihood or posterior
    # but in some cases, they can help burn in a bit
    start_mu[grepl('weight1', names(start_mu))] <- -.1
    start_mu[grepl('weight2', names(start_mu))] <- -.5
    start_mu[grepl('weight3', names(start_mu))] <- -.2

    # start with identity, fill in parameters that are known
    start_var <- diag(length(par_names_new))   #*prior_var[par_names_new]
    dimnames(start_var) <-  list(par_names_new, par_names_new)
    start_var[select_pars, select_pars] <- mu_var[select_pars, select_pars]
    for(pname in new_pars) {
      if(pname %in% names(prior_var)) start_var[pname, pname] <- prior_var[pname] # *.1  # start close to the original value
    }
    start_var['B','B'] <- 1  # Not for threshold, this most likely has to increase quite a bit compared to the null model, so generate some particles further away as well
    return(list(start_mu=start_mu, start_var=start_var))
  }
}




generate_filenames <- function(dataset, learningModel, trendModel='NULL', trendPar='NULL', nTrendPars=0, samples_dir='./samples') {
  if(trendModel == 'NULL') {
    samplers_fn <- paste0(dataset, '_model-', learningModel, '_trend-NULL.RData')
  } else {
    samplers_fn <- paste0(dataset, '_model-', learningModel, '_trend-', trendModel, trendPar, nTrendPars, '.RData')
  }
  samplers_fn <- file.path(samples_dir, samplers_fn)
  return(samplers_fn)
}

loadDataSamplersPP <- function(task, learningModel, trendModel='NULL', pp_conditional=FALSE,
                               trendPar='NULL', nTrendPars=0, samples_dir=c('samples', 'samples_mc123', 'samples_trends'), old_location=FALSE,
                               return_fns=FALSE) {

  # Load data, posterior predictives, samples
  load(paste0('./datasets/', task, '.RData'))
  dat <- EMC2:::add_trials(dat)

  foundSamplers <- FALSE
  if(length(samples_dir) > 1) {
    for(samples_dir_ in samples_dir) {
      this_samplers_fn <- generate_filenames(dataset=task, learningModel=learningModel,
                                        trendModel=trendModel, trendPar=trendPar,
                                        nTrendPars=nTrendPars,
                                        samples_dir = samples_dir_)
      if(file.exists(this_samplers_fn)) {
        if(foundSamplers) {
          print(samplers_fn); print(this_samplers_fn); warning('Multiple samplers found in multiple directories! Taking the first I found.')
        } else {
          foundSamplers <- TRUE
          samplers_fn <- this_samplers_fn
        }
      }
    }
  } else {
    samplers_fn <- generate_filenames(dataset=task, learningModel=learningModel,
                                      trendModel=trendModel, trendPar=trendPar,
                                      nTrendPars=nTrendPars,
                                      samples_dir = samples_dir)
  }
  # samplers_fn <- generate_filenames(dataset=task, learningModel=learningModel,
  #                                   trendModel=trendModel, trendPar=trendPar,
  #                                   nTrendPars=nTrendPars,
  #                                   samples_dir = samples_dir)

  if(!old_location) pp_loc <- gsub('samples', 'posteriorpredictives', samplers_fn) else pp_loc <- samplers_fn
  if(!pp_conditional) {
    pp_fn <- gsub('.RData', '_pp-unconditional.RData', pp_loc)
  } else {
    pp_fn <- gsub('.RData', '_pp-conditional.RData', pp_loc)
  }
  if(return_fns) return(list(samplers_fn=samplers_fn, pp_fn=pp_fn))

  samplers <- EMC2:::loadRData(samplers_fn)
  pp <- EMC2:::loadRData(pp_fn)

  # ensure pp and data is correctly ordered
  pp <- pp[order(pp$subjects,pp$postn,pp$trials),]
  # ensure the exact same factors
  dat$choice <- dat$R == levels(dat$R)[1]
  pp$choice <- pp$R==dat[dat$choice==1,'R'][1]

  return(list(dat=dat, pp=pp, samplers=samplers))
}

