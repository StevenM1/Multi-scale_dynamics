#### Make samplers for dynamic models
library(EMC2)
library(emcAdapt)


# makeDCT <- function(frame_times, high_pass=1/128) {
#   n_frames = length(frame_times)
#   n_times <- 1:n_frames
#
#   dt = (frame_times[length(frame_times)]-frame_times[1]) / (n_frames-1)  # should be 1 with trials
#   order = pmin(n_frames-1, floor(2*n_frames*high_pass*dt))
#   cosine_drift = matrix(0, nrow=n_frames, ncol=order+1)
#   normalizer = sqrt(2/n_frames)
#
#   for(k in seq(1, order+1)) {
#     cosine_drift[,k] = normalizer*cos((pi/n_frames)*(n_times+0.5)*k)
#   }
#   cosine_drift = cosine_drift/max(cosine_drift)
#   return(cosine_drift)
# }

#decisionModel <- 'RDM'
#learningModel <- 'zSMuAHbV'
#dctModel <- 'NULL'
#task <- 'wagenmakers2008exp2'
makeSamplers <- function(decisionModel, learningModel, task, trendType, trendPar, nTrendPar, samples_dir='./samples',
                         n_chains=3, dat=NULL, do_return=FALSE, dct_highpass=1/400, overwrite=FALSE, trend_nuisance=FALSE) {

  save_fn_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-', learningModel, '-trend-', trendType, '-', trendPar, '-', nTrendPar, '.RData'))

  if(file.exists(save_fn_samples) & !overwrite) {
    print('samplers file already exists, returning')
    return(0)
  }
  if(!dir.exists(samples_dir)) dir.create(samples_dir, recursive=TRUE)

  if(is.null(dat)) {
    print(load(file.path(samples_dir, '..', 'datasets', paste0(task, '.RData'))));
    dat$stim <- as.numeric(dat$S)
    dat$stim[dat$stim==2] <- -1
    dat$resp <- as.numeric(dat$R)
    dat$resp[dat$resp==2] <- -1

    # response frequency?
    dat$rf <- 1/dat$rt
  }

  ## general set-up
  constants1=c(s=log(1), A=log(0))
  constants <- c(alpha1min=0, alpha2min=0, alpha3min=0)
  probit_scale <- c('alpha1', 'alpha2', 'alpha3', 'q01', 'q02',
                    'Btheta2', 'utheta2', 'vtheta2')  # exponents
  logit_scale <- c('v', "A", "t0", "s", 'q03')     # don't add thresholds here
  ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
  # average+difference
  adapt <- list(useC=TRUE, useSMsApproach=TRUE, simulateByTrial=FALSE, includeDCT=TRUE,
                dynamic=list(
                  columns_to_include='accuracy',
                  adapt_fun_name='delta')
  )

  ## Combine learning effects & DCT in Flist
  Flist=list(v~lM,B~1,A~1,t0~1,s~lM) # standard RDM parameters

  # Set-up LEARNING model ---------------------------------------------------
  if(learningModel == 'zSM') {
    # z ~ stim memory ------------------------------------------------------------
    adapt$dynamic=list(
      columns_to_include='stim',
      adapt_fun_name='delta',
      output_name=c('B'),
      output_fun=function(npars, Qvalues, da) {
        Qvalues <- Qvalues[,1]
        Qvalues[as.numeric(da$lR)==2] <- -Qvalues[as.numeric(da$lR)==2]
        B = npars[,'B'] + npars[,'weight1'] * Qvalues
        B = pmin(B, 10)  # hacky hacky
        B = pmax(B, 0.01)
        return(B)
      }
    )
    constants <- c(constants,
                   q01=qnorm(0), #alpha1=qnorm(0), weight1=0,
                   q02=qnorm(0), alpha2=qnorm(0), weight2=0,
                   q03=qnorm(0), alpha3=qnorm(0), weight3=0)
    Flist <- append(Flist, list(alpha1~1, weight1~1))
  } else if(learningModel == 'vSM') {
    # v ~ stim memory ------------------------------------------------------------
    adapt$dynamic=list(
      columns_to_include='stim',
      adapt_fun_name='delta',
      output_name=c('v'),
      output_fun=function(npars, Qvalues, da) {
        Qvalues <- Qvalues[,1]
        Qvalues[as.numeric(da$lR)==2] <- -Qvalues[as.numeric(da$lR)==2]
        v = npars[,'v'] + npars[,'weight1'] * Qvalues
        v = pmax(v, 0.1)  ## hacky hacky
        v = pmin(v, 10)   ## hacky hacky
        return(v)
      }
    )
    constants <- c(constants,
                   q01=qnorm(0), #alpha1=qnorm(0), weight1=0,
                   q02=qnorm(0), alpha2=qnorm(0), weight2=0,
                   q03=qnorm(0), alpha3=qnorm(0), weight3=0)
    Flist <- append(Flist, list(alpha1~1, weight1~1))
  } else if(learningModel == 'zSMuAH')  {
    # z ~ SM, u~AH ------------------------------------------------------------
    adapt$simulateByTrial = TRUE
    adapt$dynamic=list(
      columns_to_include=c('stim', 'accuracy'),
      adapt_fun_name='delta',
      output_name=c('B', 'v'),
      output_fun=function(npars, Qvalues, da) {
        # start point effects
        Qvalues[as.numeric(da$lR)==2,1] <- -Qvalues[as.numeric(da$lR)==2,1]
        B = npars[,'B'] + npars[,'weight1'] * Qvalues[,1]
        B = pmin(B, 10)    # hacky hacky
        B = pmax(B, 0.01)

        ## urgency ~ AH and FH
        ## assuming ARD, then v_1 = v0 + (Q1-Q2)
        ## then               v_2 = v0 + (Q2-Q1)
        ## To get v0, take mean(v_1,v_2) -- (Q1-Q2 + Q2-Q1 = 0, so that cancels out)
        v <- npars[,'v']
        mean_vs <- apply(cbind(v[seq(1,length(v),2)], v[seq(2,length(v),2)]),1,mean)
        urgency <- rep(mean_vs, each=2)
        quality <- v - urgency  # v = v0 + quality

        # only the urgency component can be reduced/increased
        urgency <- urgency + npars[,'weight2']*round(1-Qvalues[,2],2) #+ npars[,'weight3']*(1/Qvalues[,3])

        # but must remain non-negative
        urgency[urgency < -.25] <- 20  # smoothen likelihood landscape: anything that is 'a bit' negative will be set to 0, but anything very negative to 20
        urgency <- pmax(urgency, 0.001)

        v <- urgency + quality
        v = pmin(v, 20)   # hacky hacky to make sure no crazy drift values occur
        v = pmax(v, 0.001)
        return(cbind(B,v))
      }
    )

    constants <- c(constants,
                   q01=qnorm(0),
                   q02=qnorm(1),
                   q03=qnorm(0.5), alpha3=qnorm(0), weight3=0)
    constants[['alpha2min']] <- 0.15   # For accuracy learning, set minimum learning rate of 0.15
    Flist <- append(Flist, list(alpha1~1, weight1~1))  # SM-related parameters
    Flist <- append(Flist, list(alpha2~1, weight2~1))  # AH-related weight, accuracy learning rate
    #Flist <- append(Flist, list(alpha3~1, weight3~1))  # RT learning rate
  } else if(learningModel == 'zSMaAH')  {
    # z ~ SM, a~AH ------------------------------------------------------------
    adapt$simulateByTrial = TRUE
    adapt$dynamic=list(
      columns_to_include=c('stim', 'accuracy'),
      adapt_fun_name='delta',
      output_name=c('B'),
      output_fun=function(npars, Qvalues, da) {
        # start point effects
        Qvalues[as.numeric(da$lR)==2,1] <- -Qvalues[as.numeric(da$lR)==2,1]
        B = npars[,'B'] + npars[,'weight1'] * Qvalues[,1] + npars[,'weight2']*round(1-Qvalues[,2],2)
        B = pmin(B, 10)    # hacky hacky
        B = pmax(B, 0.01)

        return(B)
      }
    )

    constants <- c(constants,
                   q01=qnorm(0),
                   q02=qnorm(1),
                   q03=qnorm(0.5), alpha3=qnorm(0), weight3=0)
    constants[['alpha2min']] <- 0.15   # For accuracy learning, set minimum learning rate of 0.15
    Flist <- append(Flist, list(alpha1~1, weight1~1))  # SM-related parameters
    Flist <- append(Flist, list(alpha2~1, weight2~1))  # AH-related weight, accuracy learning rate
    #Flist <- append(Flist, list(alpha3~1, weight3~1))  # RT learning rate
  } else if(learningModel == 'zSMvAH')  {
    # z ~ SM, v~AH ------------------------------------------------------------
    adapt$simulateByTrial = TRUE
    adapt$dynamic=list(
      columns_to_include=c('stim', 'accuracy'),
      adapt_fun_name='delta',
      output_name=c('B', 'v'),
      output_fun=function(npars, Qvalues, da) {
        # start point effects
        Qvalues[as.numeric(da$lR)==2,1] <- -Qvalues[as.numeric(da$lR)==2,1]
        B = npars[,'B'] + npars[,'weight1'] * Qvalues[,1] #+ npars[,'weight2']*round(1-Qvalues[,2],2)
        B = pmin(B, 10)    # hacky hacky
        B = pmax(B, 0.01)

        ## drift ~ AH
        ## assuming ARD, then v_1 = v0 + (Q1-Q2)
        ## then               v_2 = v0 + (Q2-Q1)
        ## To get v0, take mean(v_1,v_2) -- (Q1-Q2 + Q2-Q1 = 0, so that cancels out)
        v <- npars[,'v']
        mean_vs <- apply(cbind(v[seq(1,length(v),2)], v[seq(2,length(v),2)]),1,mean)
        urgency <- rep(mean_vs, each=2)
        quality <- v - urgency  # v = v0 + quality

        # only the urgency component can be reduced/increased
        quality <- quality + npars[,'weight2']*round(1-Qvalues[,2],2)*sign(quality)

        v <- urgency + quality
        v = pmin(v, 20)   # hacky hacky to make sure no crazy drift values occur
        v = pmax(v, 0.001)

        return(cbind(B, v))
      }
    )

    constants <- c(constants,
                   q01=qnorm(0),
                   q02=qnorm(1),
                   q03=qnorm(0.5), alpha3=qnorm(0), weight3=0)
    constants[['alpha2min']] <- 0.15   # For accuracy learning, set minimum learning rate of 0.15
    Flist <- append(Flist, list(alpha1~1, weight1~1))  # SM-related parameters
    Flist <- append(Flist, list(alpha2~1, weight2~1))  # AH-related weight, accuracy learning rate
  }  else if(learningModel == 'bV')  {
    # b~drift estimate ------------------------------------------------------------
    adapt$simulateByTrial = TRUE
    adapt$dynamic=list(
      columns_to_include=c('stim', 'accuracy',  'rt'),
      adapt_fun_name='threshold',
      output_name=c('B'),
      output_fun=function(npars, Qvalues, da) {
        # threshold effects
        B = npars[,'B'] + npars[,'weight3'] * Qvalues[,3]
        B = pmin(B, 10)    # hacky hacky
        B = pmax(B, 0.01)
        return(B)
      }
    )

    constants <- c(constants,
                   q01=qnorm(0), alpha1=qnorm(0), weight1=0,
                   q02=qnorm(1), alpha2=qnorm(0), weight2=0,
                   q03=log(3))   # NB: approximately the mean drift rate of correct accumulator across datasets (NULL models)
    constants[['alpha2min']] <- 0.15   # For accuracy learning, set minimum learning rate of 0.15
    #    Flist <- append(Flist, list(alpha1~1, weight1~1))  # SM-related parameters
    #    Flist <- append(Flist, list(alpha2~1, weight2~1))  # AH-related weight, accuracy learning rate
    Flist <- append(Flist, list(alpha3~1, weight3~1))  # RT learning rate
  } else if(learningModel == 'zSMuAHbV')  {
    # z ~ SM, u~AH, b~drift estimate ------------------------------------------------------------
    adapt$simulateByTrial = TRUE
    adapt$dynamic=list(
      columns_to_include=c('stim', 'accuracy',  'rt'),
      adapt_fun_name='deltadeltathreshold',
      output_name=c('B','v'),
      output_fun=function(npars, Qvalues, da) {
        # start point effects
        Qvalues[as.numeric(da$lR)==2,1] <- -Qvalues[as.numeric(da$lR)==2,1]
        B = npars[,'B'] + npars[,'weight1'] * Qvalues[,1] + npars[,'weight3']*Qvalues[,3]
        B = pmin(B, 10)     # hacky hacky
        B = pmax(B, 0.01)

        ## urgency ~ AH
        ## assuming ARD, then v_1 = v0 + (Q1-Q2)
        ## then               v_2 = v0 + (Q2-Q1)
        ## To get v0, take mean(v_1,v_2) -- (Q1-Q2 + Q2-Q1 = 0, so that cancels out)
        v <- npars[,'v']
        mean_vs <- apply(cbind(v[seq(1,length(v),2)], v[seq(2,length(v),2)]),1,mean)
        urgency <- rep(mean_vs, each=2)
        quality <- v - urgency  # v = v0 + quality

        # only the urgency component can be reduced/increased
        urgency <- urgency + npars[,'weight2']*round(1-Qvalues[,2],2) #+ npars[,'weight3']*(1/Qvalues[,3])

        # but must remain non-negative
        urgency[urgency < -.25] <- 20  # smoothen likelihood landscape: anything that is 'a bit' negative will be set to 0, but anything very negative to 20
        urgency <- pmax(urgency, 0.001)

        v <- urgency + quality
        v = pmin(v, 20)   # hacky hacky to make sure no crazy drift values occur
        v = pmax(v, 0.001)
        return(cbind(B,v))
      }
    )

    constants <- c(constants,
                   q01=qnorm(0), #alpha1=qnorm(0), weight1=0,
                   q02=qnorm(1), #alpha2=qnorm(0), weight2=0,
                   q03=log(3))  # NB: approximately the mean drift rate of correct accumulator across datasets (NULL models)
    constants[['alpha2min']] <- 0.15   # For accuracy learning, set minimum learning rate of 0.15
    Flist <- append(Flist, list(alpha1~1, weight1~1))  # SM-related parameters
    Flist <- append(Flist, list(alpha2~1, weight2~1))  # AH-related weight, accuracy learning rate
    Flist <- append(Flist, list(alpha3~1, weight3~1))  # RT/difficulty learning rate
  } else if(learningModel == 'zSMaAHbV')  {
    # z ~ SM, a~AH, b~drift estimate ------------------------------------------------------------
    adapt$simulateByTrial = TRUE
    adapt$dynamic=list(
      columns_to_include=c('stim', 'accuracy',  'rt'),
      adapt_fun_name='deltadeltathreshold',
      output_name=c('B'),
      output_fun=function(npars, Qvalues, da) {
        # start point effects
        Qvalues[as.numeric(da$lR)==2,1] <- -Qvalues[as.numeric(da$lR)==2,1]
        B = npars[,'B'] + npars[,'weight1'] * Qvalues[,1] + npars[,'weight3']*Qvalues[,3] + npars[,'weight2']*round(1-Qvalues[,2],2)
        B = pmin(B, 10)     # hacky hacky
        B = pmax(B, 0.01)

        return(B)
      }
    )

    constants <- c(constants,
                   q01=qnorm(0), #alpha1=qnorm(0), weight1=0,
                   q02=qnorm(1), #alpha2=qnorm(0), weight2=0,
                   q03=log(3))  # NB: approximately the mean drift rate of correct accumulator across datasets (NULL models)
    constants[['alpha2min']] <- 0.15   # For accuracy learning, set minimum learning rate of 0.15
    Flist <- append(Flist, list(alpha1~1, weight1~1))  # SM-related parameters
    Flist <- append(Flist, list(alpha2~1, weight2~1))  # AH-related weight, accuracy learning rate
    Flist <- append(Flist, list(alpha3~1, weight3~1))  # RT/difficulty learning rate
  }  else {
    # no learning
    adapt$dynamic=list(
      columns_to_include='stim',
      adapt_fun_name='delta'
    )
    constants <- c(constants,
                   q01=qnorm(0), alpha1=qnorm(0), weight1=0,
                   q02=qnorm(0), alpha2=qnorm(0), weight2=0,
                   q03=qnorm(0), alpha3=qnorm(0), weight3=0)
  }

  # DCT model ---------------------------------------------------------------
  if(trendType != 'DCT') {
    # no dct, all constants
    constants <- c(constants, setNames(rep(0, 10), paste0('vcos', 1:10)))
    constants <- c(constants, setNames(rep(0, 10), paste0('Bcos', 1:10)))
    constants <- c(constants, setNames(rep(0, 10), paste0('ucos', 1:10)))
  } else {
    if(nTrendPar == 0) {
      averageNTrials <- mean(aggregate(rt~subjects,dat,length)$rt)

      if(dct_highpass==1/250) {
        if(task == 'wagenmakers2008exp2') {
          nTrendPar = 10
        } else {
          nTrendPar = ncol(EMC2:::makeDCT(1:averageNTrials, high_pass=1/250))
        }
      } else {
        nTrendPar = ncol(EMC2:::makeDCT(1:averageNTrials, high_pass=dct_highpass))
        if(task == 'wagenmakers2008exp2') nTrendPar = 7
      }
    }

    if(trendPar == 'B') {
      ## drift only affects thresholds, so set effect on drift rates to 0
      constants <- c(constants, setNames(rep(0, nTrendPar), paste0('vcos', 1:nTrendPar)))
      constants <- c(constants, setNames(rep(0, nTrendPar), paste0('ucos', 1:nTrendPar)))
      Flist2 <- lapply(paste0(paste0('Bcos', 1:nTrendPar),'~1'), as.formula) # cosines for DCT
      Flist <- append(Flist, Flist2)
    } else if(trendPar == 'v') {
      ## drift only affects rates, so set effect on thresholds to 0
      constants <- c(constants, setNames(rep(0, nTrendPar), paste0('Bcos', 1:nTrendPar)))
      constants <- c(constants, setNames(rep(0, nTrendPar), paste0('ucos', 1:nTrendPar)))
      Flist2 <- lapply(paste0(paste0('vcos', 1:nTrendPar),'~1'), as.formula) # cosines for DCT
      Flist <- append(Flist, Flist2)
    } else if(trendPar == 'u') {
      constants <- c(constants, setNames(rep(0, nTrendPar), paste0('Bcos', 1:nTrendPar)))
      constants <- c(constants, setNames(rep(0, nTrendPar), paste0('vcos', 1:nTrendPar)))
      Flist2 <- lapply(paste0(paste0('ucos', 1:nTrendPar),'~1'), as.formula) # cosines for DCT
      Flist <- append(Flist, Flist2)
    }
  }

  # Exponential trend model -------------------------------------------------
  if(trendType != 'exp') {
    # constants <- c(constants, setNames(rep(0, 2), paste0('vtheta', 1:2)))
    # constants <- c(constants, setNames(rep(0, 2), paste0('Btheta', 1:2)))
    # constants <- c(constants, setNames(rep(0, 2), paste0('utheta', 1:2)))
    constants <- c(constants,
                   'vtheta1'=0, 'Btheta1'=0, 'utheta1'=0,
                   'vtheta2'=qnorm(0),'Btheta2'=qnorm(0),'utheta2'=qnorm(0))
  } else {
    if(trendPar == 'B') {
      constants <- c(constants, setNames(rep(0, 2), paste0('vtheta', 1:2)))
      constants <- c(constants, setNames(rep(0, 2), paste0('utheta', 1:2)))
      Flist2 <- lapply(paste0(paste0('Btheta', 1:2),'~1'), as.formula) # cosines for DCT
    } else if(trendPar == 'v') {
      constants <- c(constants, setNames(rep(0, 2), paste0('Btheta', 1:2)))
      constants <- c(constants, setNames(rep(0, 2), paste0('utheta', 1:2)))
      Flist2 <- lapply(paste0(paste0('vtheta', 1:2),'~1'), as.formula) # cosines for DCT
    } else if(trendPar == 'u') {
      constants <- c(constants, setNames(rep(0, 2), paste0('vtheta', 1:2)))
      constants <- c(constants, setNames(rep(0, 2), paste0('Btheta', 1:2)))
      Flist2 <- lapply(paste0(paste0('utheta', 1:2),'~1'), as.formula) # cosines for DCT
    }
    Flist <- append(Flist, Flist2)
  }

  if(trendType != 'poly') {
    constants <- c(constants, setNames(rep(0, 10), paste0('vpoly', 1:10)))
    constants <- c(constants, setNames(rep(0, 10), paste0('Bpoly', 1:10)))
    constants <- c(constants, setNames(rep(0, 10), paste0('upoly', 1:10)))
  } else {
    if(trendPar == 'B') {
      constants <- c(constants, setNames(rep(0, 10), paste0('vpoly', 1:10)))
      constants <- c(constants, setNames(rep(0, 10), paste0('upoly', 1:10)))
      Flist2 <- lapply(paste0(paste0('Bpoly', 1:nTrendPar),'~1'), as.formula) # polynomials
    } else if(trendPar == 'v') {
      constants <- c(constants, setNames(rep(0, 10), paste0('Bpoly', 1:10)))
      constants <- c(constants, setNames(rep(0, 10), paste0('upoly', 1:10)))
      Flist2 <- lapply(paste0(paste0('vpoly', 1:nTrendPar),'~1'), as.formula) # polynomials
    } else if(trendPar == 'u') {
      constants <- c(constants, setNames(rep(0, 10), paste0('vpoly', 1:10)))
      constants <- c(constants, setNames(rep(0, 10), paste0('Bpoly', 1:10)))
      Flist2 <- lapply(paste0(paste0('upoly', 1:nTrendPar),'~1'), as.formula) # polynomials
    }
    Flist <- append(Flist, Flist2)

  }

  p_types <- sapply(Flist, function(x) all.vars(x)[1])
  model <- generateRDMdynamic(p_types=p_types, constantsFixed=constants, probit_scale = probit_scale, logit_scale = logit_scale)

  # Design ------------------------------------------------------------------
  if(task == 'wagenmakers2008exp2') {
    Flist[[1]] = v~lM*W
    if(!grepl('SM', learningModel)) {
      ## the null model to compare the z and v models against is one where
      ## we tell the model that thresholds should differ between conditions
      Flist[[2]] = B~proportion_word*lR
    }
    Ffactors=list(subjects=levels(dat$subjects),
                  S=c('nonword', 'word'),
                  W=levels(dat$W),
                  proportion_word=levels(dat$proportion_word))
    Ffunctions=list(errorAccumulator=function(d) as.character(d$lM)=="FALSE")
    Flist[[5]] = s~errorAccumulator

    constants1=c(constants1, v_Wvlf=0, v_Whf=0, v_Wnonword=0)
    Clist=NULL
  } else if(task == 'forstmann2008') {
    # model speed-accuracy tradeoff
    Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
    dimnames(Emat) <- list(NULL,c("a-n","a-s"))
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S),
                  E=levels(dat$E))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat)
    Flist[[2]] = B~E+lR
  } else if(task == 'degee2014') {
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat)
    Flist[[1]] = v~lM*S  # drift rates for S=left and S=right are *very* different!
  } else if(task == 'mileticvanmaanen2019exp2block1') {
    Flist[[1]] = v~lM*difficulty
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S),
                  difficulty=levels(dat$difficulty))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat)
    constants1=c(constants1, v_difficultyD2=0, v_difficultyD3=0, v_difficultyD4=0, v_difficultyD5=0,
                 v_lMd=0)
  } else {
    Ffunctions=NULL
    Ffactors=list(subjects=levels(dat$subjects),
                  S=levels(dat$S))
    Clist=list(lM=ADmat,lR=ADmat,S=ADmat)
  }

  ### THIS PRESERVES MEMORY
  # all closures (ie functions & formulas) have their own enclosing environment -- set that to an empty environment for much better memory performance
  if('output_fun' %in% names(adapt$dynamic)) {
    environment(adapt$dynamic$output_fun) <- new.env(parent=globalenv())
  }
  matchfun = function(d) d$S==d$lR
  environment(matchfun) <- new.env(parent=globalenv())
  Flist <- lapply(Flist, function(x) {environment(x) <- new.env(parent=globalenv()); return(x)})
  ###

  ## Design
  design <- make_design(
    Ffactors=Ffactors,
    Rlevels=levels(dat$R),
    matchfun=matchfun,
    Clist=Clist,
    Flist=Flist,
    Ffunctions=Ffunctions,
    constants=constants1,
    adapt=adapt,   # Special adapt component
    model=model)
  hyper <- names(attr(design, 'p_vector'))

  if(trend_nuisance) {
     nuisance_non_hyper <- which(grepl('(Bcos|ucos|vcos|poly)', hyper)) # | grepl('vcos', names(attr(design, 'p_vector'))))
     if(length(nuisance_non_hyper) > 0) {
       hyper <- hyper[-nuisance_non_hyper]
     }
  }
  theta_mu_mean <- setNames(rep(0, length(hyper)), hyper)
  theta_mu_var <- setNames(rep(1, length(hyper)), hyper)

  ## priors on weights and learning rates mildly informative, centered at 0
  ## helps stabilize sampling for some models/datasets
  theta_mu_var[grepl('weight1', hyper)] <- 0.1
  theta_mu_var[grepl('alpha1', hyper)] <- 0.1
  theta_mu_var[grepl('weight2', hyper)] <- 0.1
  theta_mu_var[grepl('alpha2', hyper)] <- 0.1
  theta_mu_var[grepl('weight3', hyper)] <- 0.1
  theta_mu_var[grepl('alpha3', hyper)] <- 0.1

  ## Trends
  ## In case there's any exponential trends, set prior var to 0.01
  theta_mu_var[grepl('theta1', hyper)] <- 0.01  # natural scale
  theta_mu_var[grepl('theta2', hyper)] <- 0.1  # NB: this parameter is on the probit scale
  # cosines, polynomials    # NB: I changed this so that the basis sets now max to 0.01 instead of 1
  # maybe easier for the sampler because the parameters will be in a slightly more 'similar' range as the other pars?
  # theta_mu_var[grepl('(cos|poly)', hyper)] <- 0.01

  if(!trend_nuisance) {
    prior <- list(
      theta_mu_mean=theta_mu_mean,
      theta_mu_var=diag(theta_mu_var)
    )
  } else {
    if(length(nuisance_non_hyper) > 0) {
      if(length(nuisance_non_hyper) > 1) {
        prior <- list(
          theta_mu_mean=theta_mu_mean,
          theta_mu_var=diag(theta_mu_var),
          prior_nuis=list(
            theta_mu_mean = rep(0, length(nuisance_non_hyper)),
            theta_mu_var = diag(rep(0.1, length(nuisance_non_hyper)))
          )
        )
      } else {
        prior <- list(
          theta_mu_mean=theta_mu_mean,
          theta_mu_var=diag(theta_mu_var),
          prior_nuis=list(
            theta_mu_mean = 0,
            theta_mu_var = 0.1
          )
        )
      }
    } else {
      prior <- list(
        theta_mu_mean=theta_mu_mean,
        theta_mu_var=diag(theta_mu_var)
      )
    }
  }
  if(trend_nuisance) {
    samplers <- make_samplers(dat, design, n_chains = n_chains,
                             nuisance_non_hyper = nuisance_non_hyper,
                              type='standard', prior=prior, rt_resolution = 0.001)
  } else {
    samplers <- make_samplers(dat, design, n_chains = n_chains,
                              # nuisance_non_hyper = nuisance_non_hyper,
                              type='standard', prior=prior, rt_resolution = 0.001)
  }
  if(do_return) return(samplers)
  save(samplers, file=save_fn_samples)
  print(paste0("I should have saved to", save_fn_samples))
}
