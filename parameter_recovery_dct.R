# parameter recovery
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  #### Run from command line
  i <- as.integer(args[1])+1
  # sampling settings
  cores_for_chains = as.integer(args[2])
  cores_per_chain = as.integer(args[3])
}

library(EMC2)
library(emcAdapt)
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMsNewEMC')
source(file.path(root_dir, 'extra_EMC2_functions/adaptive.R'))
source(file.path(root_dir, 'extra_EMC2_functions/model_RDMdynamic.R'))
source(file.path(root_dir, 'extra_EMC2_functions/utils.R'))
source(file.path(root_dir, 'extra_EMC2_functions/make_data.R'))

samples_dir=file.path(root_dir, 'samples_trends')
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
learningModels <- c('zSMuAHbV', 'zSMuAHbV', 'zSMuAHbV','zSMuAHbV')
trendModels <- rep('DCT', 4)
trendPars <- c('B', 'v', 'v', 'B')
nTrendParss <- rep(3, 4)
cores_for_chains = 3
cores_per_chain = 15


dataset <- tasks[i]
learningModel <- learningModels[i]
trendModel <- trendModels[i]
trendPar <- trendPars[i]
nTrendPars <- nTrendParss[i]

print(dataset)

simData_fn <- file.path(root_dir, 'parameter_recovery', 'datasets', paste0(dataset, '_', learningModel, '_', trendModel, nTrendPars, trendPar, '_simulated_data.RData'))
samplers_fn <- generate_filenames(dataset=dataset, learningModel=learningModel, trendModel=trendModel, trendPar=trendPar, nTrendPars=nTrendPars, samples_dir = samples_dir)
samplers <- EMC2:::loadRData(samplers_fn)

if(!file.exists(simData_fn)) {
  simulated_data <- predict(samplers, n_post=1, stat = 'mean', conditionalOnData=FALSE)
  if('errorAccumulator' %in% colnames(simulated_data)) {
    pars <- attr(simulated_data, 'pars')
    npars <- attr(simulated_data, 'npars')
    simulated_data <- simulated_data[,!grepl('errorAccumulator', colnames(simulated_data))]
    attr(simulated_data, 'pars') <- pars
    attr(simulated_data, 'npars') <- npars
  }
  save(simulated_data, file=simData_fn)
} else {
  load(simData_fn)
}

samplers_fn <- gsub('samples_trends', 'parameter_recovery/samples', samplers_fn)
## Fitting routine copied from fitDynamicEAMs.R
if(!file.exists(samplers_fn)) {
  # Make samplers
  ## Load data
  dat <- EMC2:::loadRData(simData_fn)
  if(!'stim' %in% colnames(dat)) {
    dat$stim <- as.numeric(dat$S)
    dat$stim[dat$stim==2] <- -1
    save(dat, file=simData_fn)
  }
  if('proportion_word' %in% colnames(dat)) dat$proportionWord <- dat$proportion_word  # new EMC doesnt like underscores in variable names
  if('errorAccumulator' %in% colnames(dat)) dat <- dat[,!grepl('errorAccumulator', colnames(dat))]

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
    if(!grepl('SM', learningModel)) {  # | ds4_with_baseline) {
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
  save(emc, file=samplers_fn)
  do_fit <- TRUE
} else {
  emc <- EMC2:::loadRData(samplers_fn)
  # check if already done
  stop_criteria=list(max_gd=1.1, omit_mpsrf=TRUE, selection=c('alpha', 'mu'), iter=1000)
  done = EMC2:::check_progress(emc, stage='sample', iter=1000, stop_criteria = stop_criteria,
                               max_tries=50, step_size=100, n_cores=cores_for_chains*cores_per_chain, verbose=FALSE)$done
  if(done) do_fit <- FALSE else do_fit <- TRUE
}

if(do_fit) emc <- fit(emc, cores_per_chain = cores_per_chain, cores_for_chains = cores_for_chains,
                      fileName=samplers_fn, verbose=TRUE, verboseProgress=TRUE, max_tries = 20)


## plot
trendModels <- rep('DCT', 4)
trendPars <- c('B', 'v', 'v', 'B')
nTrendParss <- rep(3, 4)
learningModels <- c('zSMuAHbV', 'zSMuAHbV', 'zSMuAHbV','zSMuAHbV')
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
for(i in 1:length(tasks)) {
  dataset <- task <- tasks[i]
  learningModel <- learningModels[i]
  trendModel <- trendModels[i]
  trendPar <- trendPars[i]
  nTrendPars <- nTrendParss[i]

  simData_fn <- file.path(root_dir, 'parameter_recovery', 'datasets', paste0(dataset, '_', learningModel, '_', trendModel, nTrendPars, trendPar, '_simulated_data.RData'))
  simulatedData <- EMC2:::loadRData(simData_fn)
  samplers_fn <- generate_filenames(dataset=dataset, learningModel=learningModel, trendModel=trendModel, trendPar=trendPar, nTrendPars=nTrendPars, samples_dir = './parameter_recovery/samples')
  samples <- EMC2:::loadRData(samplers_fn)

  #subfilter=0
  idx <- which(samples[[1]]$samples$stage=='sample')
  # idx <- idx[(length(idx)-subfilter):length(idx)]
  alpha <- abind(lapply(samples, function(x) x$samples$alpha[,,idx]), along=3)
  mu <- abind(lapply(samples, function(x) x$samples$theta_mu[,idx]), along=2)
  fitParsMean <- apply(alpha,1:2,mean)
  fitParsQps <- apply(alpha,1:2,quantile, c(0.025, .5, 0.975))
  fitMuQps <- apply(mu, 1, quantile, c(0.025, 0.5, 0.975))

  truePars <- attr(simulatedData, 'pars')

  pname_mapping = list('v'='$u$',
                       'v_lMd'='$dv$',
                       'B'='b',
                       'B_Ea-n'='$b_{a-n}$',
                       'B_Ea-s'='$b_{a-s}$',
                       'B_lRd'='$b_{l-r}$',
                       't0'='$t_0$',
                       's_lMd'='$ds$',
                       'alpha3'='$\\alpha_{FM}$',
                       'weight3'='$w_{FM}$',
                       'alpha2'='$\\alpha_{AM}$',
                       'weight2'='$w_{AM}$',
                       'alpha1'='$\\alpha_{SM}$',
                       'weight1'='$w_{SM}$',
                       's_errorAccumulatorTRUE'='$ds$',
                       'v_lMTRUE'='$dv$',
                       'v_lMTRUE:Wvlf'='$dv:W_{vlf}$',
                       'v_lMTRUE:Whf'='$dv:W_{hf}$',
                       'v_lMTRUE:Wnonword'='$dv:W_{hf}$')

  for(ii in 1:10) {
    pname_mapping[paste0('vcos', ii)] = paste0('$\\beta_{',ii,'}$')
    pname_mapping[paste0('Bcos', ii)] = paste0('$\\beta_{',ii,'}$')
  }


  nrows <- 3
  if(task %in% c('forstmann2008', 'wagenmakers2008exp2')) {
    if(trendModel!='NULL') nrows <- nrows+1
  }

  for(ftype in c('jpeg', 'pdf')) {
    fn = paste0('./figures/parameter_recovery_', task, '_dct.pdf')
    if(ftype == 'pdf') pdf(file=fn, width=8, height=nrows*2)
    if(ftype == 'jpeg') jpeg(gsub('.pdf','.jpeg',fn), width=8, height=nrows*2, units='in', quality=100, res=500)

    par(mfrow=c(nrows,5), mar=c(3,3,2,.5), bty='l', mgp=c(2,.5,0), tck=-0.05)
    for(pname in colnames(truePars)) {
      pname_label = pname_mapping[[pname]]
      x<-truePars[,pname]
      y<-fitParsQps[2,pname,]
      xlim <- ylim <- range(c(x,fitParsQps[,pname,]))

      CI_contains_x <- fitParsQps[3,pname,] > x & fitParsQps[1,pname,] < x

      xmean <- mean(x)

      xlab <- ifelse(which(pname == colnames(truePars))>(nrows-1)*5, 'Data generating', '')
      ylab <- ifelse(which(pname == colnames(truePars)) %% 5 == 1, 'Posterior (95% CI)', '')
      plot(x,y, main=TeX(pname_label), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
      arrows(x0=x, y0=fitParsQps[1,pname,], y1=fitParsQps[3,pname,], code=3, angle=90, length=.025, col=adjustcolor(1, alpha.f=.5))
      legend('topleft', paste0('r=',round(cor(x,y),2)), bty='n', cex=.8)
      legend('bottomright', paste0(round(mean(CI_contains_x)*100,1), '%'), bty='n', cex=.8)
      abline(a=0,b=1)

      # group-level recovery plot
      # if(pname %in% colnames(fitMuQps)) {
      #   usr <- par('usr')
      #   arrows(x0=usr[1],
      #          y0=fitMuQps[1,pname],
      #          y1=fitMuQps[3,pname],
      #          code=3,
      #          angle=90,
      #          length=0.025, col=2, lwd=2,
      #          xpd=TRUE
      #   )
      #   segments(x0=xmean, x1=xmean, y0=usr[3], y1=xmean, col=2, lwd=2)
      #   segments(x0=xmean, x1=usr[1], y0=xmean, y1=xmean, col=2, lwd=2)
      # }
    }
    dev.off()
  }
}
