#### Dynamic
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  #### Run from command line
  learningModel <- args[1]
  dataset <- args[2]

  # sampling settings
  cores_for_chains = as.integer(args[3])
  cores_per_chain = as.integer(args[4])
  max_tries = as.integer(args[5])

  # trend model settings
  if(length(args) >= 6) trendModel <- args[6] else trendModel <- 'NULL'
  if(length(args) >= 7) trendPar <- args[7] else trendPar <- 'NULL'
  if(length(args) >= 8) nTrendPars <- as.integer(args[8]) else nTrendPars <- 3
  # print(args)
  debugMode <- FALSE
} else {

  #### Run interactively
  rm(list=ls())

  # model and dataset
  learningModel <- 'vSMbAHbV'
  ds_num <- 1
  dataset <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')[ds_num]

  # Sampling settings
  cores_per_chain=6
  cores_for_chains=3
  max_tries <- 20

  # DCT settings
  trendModel='NULL'   #'NULL'  # or 'DCT'
  trendPar='NULL'       # 'NULL'    # or 'v', 'u', 'B'
  nTrendPars=3       # probably always 3

  debugMode <- FALSE
}


# Specifics for dataset 4: allow bias manipulation to affect thresholds?
ds4_with_baseline <- FALSE


# Source code -------------------------------------------------------------
library(EMC2)
library(emcAdapt)
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMsNewEMC')
source(file.path(root_dir, 'extra_EMC2_functions/adaptive.R'))
source(file.path(root_dir, 'extra_EMC2_functions/model_RDMdynamic.R'))
source(file.path(root_dir, 'extra_EMC2_functions/utils.R'))


# Directories, file names -------------------------------------------------
samples_dir = file.path(root_dir, 'samples')
if(trendModel == 'DCT') samples_dir <- gsub('samples', 'samples_trends', samples_dir)
figures_dir = file.path(root_dir, 'figures')
ds_path <- file.path(root_dir, 'datasets', paste0(dataset, '.RData'))
samplers_fn <- generate_filenames(dataset=dataset, learningModel=learningModel, trendModel=trendModel, trendPar=trendPar, nTrendPars=nTrendPars, samples_dir = samples_dir)
print(samplers_fn)

if(!file.exists(samplers_fn)) {
  # Make samplers
  ## Load data
  load(ds_path)
  if(!'stim' %in% colnames(dat)) {
    dat$stim <- as.numeric(dat$S)
    dat$stim[dat$stim==2] <- -1
    save(dat, file=ds_path)
  }
  if('proportion_word' %in% colnames(dat)) dat$proportionWord <- dat$proportion_word  # new EMC doesnt like underscores in variable names

  # Default Flist (formula)
  Flist=list(v~lM,B~1,A~1,t0~1,s~lM,q01~1,q02~1,q03~1,
             alpha1~1,weight1~1,alpha2~1,weight2~1,alpha3~1,weight3~1)

  # Adapt component
  adapt <- list(useC=TRUE, useSMsApproach=TRUE, simulateByTrial=TRUE, includeDCT=TRUE,
                dynamic=list(
                  learningModel=learningModel,
                  columns_to_include=c('stim', 'accuracy',  'rt'),
                  output_name=c('B','v'),
                  adapt_fun_name=ifelse(grepl('bV', learningModel), 'deltadeltathreshold', 'delta')
                ))

  # Default constants
  constants <- c(s=log(1), A=log(0),
                 q01=0,
                 q02=1,
                 q03=3)  # NB: approximately the mean drift rate of correct accumulator across datasets (NULL models)
  if(grepl('uV', learningModel)) constants[['q03']] <- 0.7  # in case of uV, q03 is the RT estimate (not drift rate estimate). Take a rough estimate of RTs

  # which learning rates/weights to estimate?
  if(!grepl('SM', learningModel)) constants <- c(constants, alpha1=qnorm(0), weight1=0)
  if(!grepl('AH', learningModel)) constants <- c(constants, alpha2=qnorm(0), weight2=0)
  if(!(grepl('uV', learningModel) | grepl('bV', learningModel))) constants <- c(constants, alpha3=qnorm(0), weight3=0)

  ## add DCT/cosines
  if(trendModel == 'DCT') {
    Flist <- append(Flist, lapply(paste0(paste0(trendPar, 'cos', 1:nTrendPars),'~1'), as.formula))
  }

  # Generate model
  p_types_default <-c("v"=1, "B"=1, "A"=0, "t0"=0.1, "s"=1,
                      "alpha1"=-Inf, "alpha2"=-Inf, "alpha3"=-Inf,
                      "q01"=0, "q02"=1, "q03"=3,
                      "weight1"=0, "weight2"=0, "weight3"=0,
                      setNames(rep(0,3), paste0(rep('ucos', 3), 1:3)),     # effect on URGENCY
                      setNames(rep(0,3), paste0(rep('Bcos', 3), 1:3)),     # effect on THRESHOLD
                      setNames(rep(0,3), paste0(rep('vcos', 3), 1:3)),     # effect on DRIFT
                      setNames(rep(0,2), paste0(rep('vtheta', 2), 1:2)),   # exponential trend on drift
                      setNames(rep(0,2), paste0(rep('utheta', 2), 1:2)),   # exponential trend on urgency
                      setNames(rep(0,2), paste0(rep('Btheta', 2), 1:2)))
  p_types <- p_types_default[sapply(Flist, function(x) all.vars(x)[attr(terms(x), "response")])]

  alpha2min <- ifelse(grepl('AH', learningModel), 0.05, 0)  # only apply minimum LR to models that include the mechanism
  alpha3min <- ifelse(grepl('V', learningModel), 0.05, 0)
  model <- generateRDMDynamic(
    p_types=p_types,
    transform=list(func=c(alpha1='pnorm', alpha2='pnorm', alpha3='pnorm',
                          v='exp', A='exp', t0='exp', s='exp'),
                   lower=c(alpha2=alpha2min, alpha3=alpha3min)
    ),
    bound=list(minmax=cbind(v=c(1e-3,20),B=c(0.001,20),A=c(1e-4,Inf),t0=c(0.05,Inf),s=c(0,Inf),
                            alpha1=c(pnorm(-4), pnorm(4)),
                            alpha2=c(alpha2min+pnorm(-4), pnorm(4)),  # Lower bound on second learning rate
                            alpha3=c(alpha3min+pnorm(-4), pnorm(4))),
               exception=c(A=0, v=0, alpha1=0, alpha2=0, alpha3=0))
  )


  # Generate design ------------------------------------------------------------------
  ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  # average+difference
  if(dataset == 'wagenmakers2008exp2') {
    Flist[[1]] = v~lM*W
    if(!grepl('SM', learningModel) | ds4_with_baseline) {
      ## the null model to compare the z and v models against is one where
      ## we tell the model that thresholds should differ between conditions
      Flist[[2]] = B~proportionWord*lR
    }
    Ffactors=list(subjects=levels(dat$subjects),
                  S=c('nonword', 'word'),
                  W=levels(dat$W),
                  proportionWord=levels(dat$proportionWord))
    Ffunctions=list(errorAccumulator=function(d) as.character(d$lM)=="FALSE")
    Flist[[5]] = s~errorAccumulator

    constants=c(constants, v_Wvlf=0, v_Whf=0, v_Wnonword=0)
    Clist=NULL
  } else if(dataset == 'forstmann2008') {
    # model speed-accuracy tradeoff
    Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
    dimnames(Emat) <- list(NULL,c("a-n","a-s"))
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S),
                  E=levels(dat$E))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat)
    Flist[[2]] = B~E+lR
  } else if(dataset == 'degee2014') {
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat)
    Flist[[1]] = v~lM*S  # drift rates for S=left and S=right are *very* different!
  } else if(dataset == 'mileticvanmaanen2019exp2block1') {
    Flist[[1]] = v~lM*difficulty
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S),
                  difficulty=levels(dat$difficulty))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat)
    constants=c(constants, v_difficultyD2=0, v_difficultyD3=0, v_difficultyD4=0, v_difficultyD5=0,
                v_lMd=0)
  } else {
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat)
  }

  # Design
  design <- design(
    factors=Ffactors,
    Rlevels=levels(dat$R),
    matchfun=function(d) d$S==d$lR,
    contrasts=Clist,
    formula=Flist,
    functions=Ffunctions,
    constants=constants,
    model=model)

  # Prior -------------------------------------------------------------------
  prior_sd <- c(alpha1=.1,  alpha2=.1, alpha3=.1, weight1=.1, weight2=.1, weight3=.1,
                ucos1=1e-4,ucos2=1e-4,ucos3=1e-4,
                Bcos1=1e-4,Bcos2=1e-4,Bcos3=1e-4,
                vcos1=1e-4,vcos2=1e-4,vcos3=1e-4)
  prior <- prior(design,  psd = sqrt(prior_sd), type="standard")
  # plot(prior)

  # Make emc object ----------------------------------------------------------------
  emc <- make_emc(dat, design, compress=FALSE, prior=prior, rt_resolution=0.001)
  emc <- augment_emc(emc=emc, adapt=adapt)   ## add adapt and augment to dadms

  # Initialize with reasonable start points ---------------------------------------
  # try to find null model
  null_fn <- file.path(root_dir, '..', 'dynamicEAMs',  'samples', paste0(dataset, '_model-RDM-NULL-DCT-NULL.RData'))
  if(file.exists(null_fn)) {
    out = get_start_mu_var(EMC2:::loadRData(null_fn), par_names_new=emc[[1]]$par_names, prior_var=sqrt(prior_sd))
    emc <- parallel::mclapply(emc, FUN=EMC2:::init, start_mu=out[['start_mu']], start_var = out[['start_var']], n_cores=cores_per_chain, mc.cores=cores_for_chains)
    class(emc) <- 'emc'
  }
  save(emc, file=samplers_fn)
} else {
  emc <- EMC2:::loadRData(samplers_fn)
  # check if already done
  stop_criteria=list(max_gd=1.1, omit_mpsrf=TRUE, selection=c('alpha', 'mu'), iter=1000)
  done = EMC2:::check_progress(emc, stage='sample', iter=1000, stop_criteria = stop_criteria,
                               max_tries=50, step_size=100, n_cores=cores_for_chains*cores_per_chain, verbose=FALSE)$done
  if(done) stop('Already finished!')
}


# Ready for sampling
# Fit ---------------------------------------------------------------------
stop_criteria <- vector("list", length = 4)
stages_names <- c('preburn', 'burn', 'adapt', 'sample')
names(stop_criteria) <- stages_names
stop_criteria <- mapply(EMC2:::get_stop_criteria, stages_names,
                        stop_criteria, MoreArgs = list(type = emc[[1]]$type))
stop_criteria$burn$max_gd = 1.75   ## burn well!

emc <- fit(emc, cores_per_chain = cores_per_chain, cores_for_chains = cores_for_chains, stop_criteria=stop_criteria,
           fileName=samplers_fn, verbose=TRUE, verboseProgress=TRUE, max_tries = max_tries)


## Posterior predictives
source(file.path(root_dir, 'extra_EMC2_functions/make_data.R'))
pp_unconditional <- predict(emc, n_post = 100, n_cores=cores_for_chains*cores_per_chain, conditionalOnData=FALSE)
save(pp_unconditional, file=gsub('.RData', '_pp-unconditional.RData', gsub('/samples', '/posteriorpredictives', samplers_fn)))
pp_conditional <- predict(emc, n_post = 100, n_cores=cores_for_chains*cores_per_chain, conditionalOnData=TRUE)
save(pp_conditional, file=gsub('.RData', '_pp-conditional.RData', gsub('/samples', '/posteriorpredictives', samplers_fn)))


## done
