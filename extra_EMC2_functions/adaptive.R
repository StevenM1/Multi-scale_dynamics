## ALL FUNCTIONS OVERWRITTEN BY SM
augment_emc <- function(emc, adapt) {
  for(chain in 1:length(emc)) {
    for(subject in names(emc[[1]]$data)) {
      to_add <- list(augment(subject, emc[[chain]]$data[[subject]], design=list(adapt=adapt)),
                     design=adapt)
      names(to_add)[1] <- subject
      attr(emc[[chain]]$data[[subject]], 'adapt') <- to_add
    }
  }
  return(emc)
}

apply_mechanisms <- function(npars, Qvalues, da, learningModel) {
  # SM: start point effects
  Qvalues[as.numeric(da$lR)==2,1] <- -Qvalues[as.numeric(da$lR)==2,1]
  B = npars[,'B']
  if(grepl('zSM', learningModel)) B <- B + npars[,'weight1'] * Qvalues[,1]
  if(grepl('bAH', learningModel)) B <- B + npars[,'weight2'] * round(1-Qvalues[,2],2)
  if(grepl('bV', learningModel))  B <- B + npars[,'weight3'] * Qvalues[,3]
  if(grepl('(zSM|bV|bAH)', learningModel)) {
    B = pmin(B, 10)     # hacky hacky
    B = pmax(B, 0.01)
  }

  ## urgency ~ AH
  ## assuming ARD, then v_1 = v0 + (Q1-Q2)
  ## then               v_2 = v0 + (Q2-Q1)
  ## To get v0, take mean(v_1,v_2) -- (Q1-Q2 + Q2-Q1 = 0, so that cancels out)
  v <- npars[,'v']
  mean_vs <- apply(cbind(v[seq(1,length(v),2)], v[seq(2,length(v),2)]),1,mean)
  urgency <- rep(mean_vs, each=2)
  quality <- v - urgency  # v = v0 + quality

  if(grepl('vAH', learningModel)) quality <- quality + npars[,'weight2'] * round(1-Qvalues[,2],2)*sign(quality)
  if(grepl('vSM', learningModel)) quality <- quality + npars[,'weight1'] * Qvalues[,1]

  # only the urgency component can be reduced/increased
  if(grepl('uAH', learningModel)) urgency <- urgency + npars[,'weight2']*round(1-Qvalues[,2],2)
  mean_bs <- apply(cbind(B[seq(1,length(B),2)], B[seq(2,length(B),2)]),1,mean)
  if(grepl('uV',  learningModel)) urgency <- urgency + npars[,'weight3']*(mean_bs/Qvalues[,3])

  # but must remain non-negative
  if(grepl('(uAH|uV)', learningModel)) {
    urgency[urgency < -.25] <- 20  # smoothen likelihood landscape: anything that is 'a bit' negative will be set to 0, but anything very negative to 20
    urgency <- pmax(urgency, 0.001)
  }

  v <- urgency + quality
  if(grepl('(uAH|uV|vAH|vSM)', learningModel)) {
    v = pmin(v, 20)   # hacky hacky to make sure no crazy drift values occur
    v = pmax(v, 0.01)
  }
  return(cbind(B,v))
}

makeDCT <- function(frame_times, high_pass=1/128) {
  n_frames = length(frame_times)
  n_times <- 1:n_frames

  dt = (frame_times[length(frame_times)]-frame_times[1]) / (n_frames-1)  # should be 1 with trials
  order = pmin(n_frames-1, floor(2*n_frames*high_pass*dt))
  cosine_drift = matrix(0, nrow=n_frames, ncol=order+1)
  normalizer = sqrt(2/n_frames)

  for(k in seq(1, order+1)) {
    cosine_drift[,k] = normalizer*cos((pi/n_frames)*(n_times+0.5)*k)
  }
  cosine_drift = cosine_drift/max(cosine_drift)#/100    # THIS IS NEW -- instead of priors on 0.01
  return(cosine_drift)
}

augment = function(s,da,design)
  # Adds attributes to augmented data
  # learn: empty array for Q values with dim = choice alternative (low,high) x
  #   stimulus x trials (max across stimuli)
  # index: look up (row number) for stimuli in da, a matrix dim  = max trials x
  #   stimulus matrix (rows for each choice alternative contiguous)
{
  getIndex <- function(typei,cname,da,maxn) {
    out <- which(da[,cname]==typei)
    c(out,rep(NA,maxn-length(out)))
  }

  getIndexThis <- function(da, outcomes) {
    ## gets index of Q-value corresponding to the accumulator
    # NB: outcomes is sorted by trial order, da is *not*
    # so we need to return an index that takes this into account
    #da$lS <- da[cbind(1:nrow(da), match(as.character(da$lR), colnames(da)))]
    da$colNumbers <- match(as.character(da$lS), colnames(outcomes))

    return(cbind(da$trials, da$colNumbers))
  }

  getIndexOther <- function(da, outcomes) {
    ## ONLY WORKS FOR 2AFC!!
    ## gets index of Q-value corresponding to the *OTHER* accumulator (in an advantage framework)
    lROptions <- unique(as.character(da$lR))
    da$lRother <- ifelse(da$lR==lROptions[1], lROptions[2], lROptions[1])
    da$lSOther <- da[cbind(1:nrow(da), match(paste0('s_', as.character(da$lRother)), colnames(da)))]
    da$colNumbers <- match(as.character(da$lSOther), colnames(outcomes))
    return(cbind(da$trials, da$colNumbers))
    # return(cbind(da$trials, match(da[cbind(da$trials, match(da$lRother, colnames(da)))], colnames(outcomes))))
  }

  makepArray <- function(x) {
    stim <- design$adapt$stimulus$targets
    x <- x[x$R==x$lR,]          # NB: this ONLY allows the winning accumulator/stimulus to receive feedback!!
    x <- x[order(x$trials),]    # *must* be ordered by trial n for this matrix

    pArray <- matrix(NA, nrow=max(x$trials), ncol=length(stim))
    colnames(pArray) <- stim
    for(i in 1:nrow(x)) {
      trial <- x[i,'trials']
      pArray[trial,as.character(x[i,'s_left'])] <- x[i,'p_left']
      pArray[trial,as.character(x[i,'s_right'])] <- x[i,'p_right']
    }
    return(pArray)
  }

  makeOutcomes <- function(x) {
    stim <- design$adapt$stimulus$targets
    x <- x[x$R==x$lR,]          # NB: this ONLY allows the winning accumulator/stimulus to receive feedback!!
    # x <- x[x$winner,]
    x <- x[order(x$trials),]    # *must* be ordered by trial n for this matrix

    outcomes <- data.frame(matrix(NA, nrow=max(x$trials), ncol=length(stim)))
    colnames(outcomes) <- stim
    for(trial in unique(x$trials)) {
      outcomes[trial,as.character(x[trial,'lS'])] <- x[trial,'reward']
    }
    return(outcomes)
  }

  if (!is.null(design$adapt$stimulus)) {
    ## FOR RL ONLY -- not for dynamic
    targets <- design$adapt$stimulus$targets
    par <- design$adapt$stimulus$output_name
    #    maxn <- max(sapply(dimnames(targets)[[1]], function(x){table(da[da$subjects==s,x])}))

    # da index x stimulus
    # out <- sapply(targets[1,],getIndex,cname=dimnames(targets)[[1]][1],
    #               da=da[da$subjects==s,],maxn=maxn)
    stimulus <- list() #index=out)
    # accumulator x stimulus x trials
    ## SM: why not stimulus x trials, and match to accumulators at a later point (via getIndex)? would make using c much easier
    # stimulus$learn <- array(NA,dim=c(dim(targets),maxn/dim(targets)[1]),
    #                          dimnames=list(rownames(targets),targets[1,],NULL))
    stimulus$targets <- targets
    stimulus$par <- par

    if('trials' %in% colnames(da)) {
      outcomes <- makeOutcomes(da[da$subjects==s,])
      pArray <- makepArray(da[da$subjects==s,])
      return(list(stimulus=stimulus, outcomes=outcomes, pArray=pArray,
                  index=getIndexThis(da[da$subjects==s,], outcomes),
                  indexOther=getIndexOther(da[da$subjects==s,], outcomes)))
    }
  } else if(!is.null(design$adapt$dynamic)) {
    # For dynamic, NOT RL
    columns_to_include <- design$adapt$dynamic$columns_to_include
    if('trials' %in% colnames(da)) {
      tmp <- da[da$subjects==s,]
      tmp <- tmp[order(tmp$trials),]

      # unique trials only!
      tmp <- tmp[tmp$lM==TRUE,]

      outcomes <- tmp[, columns_to_include]
      if(length(columns_to_include) == 1) outcomes <- matrix(outcomes, ncol=1)
      index <- match(da[da$subjects==s,'trials'], tmp$trials)

      if(is.null(design$adapt$includeDCT)) {
        return(list(outcomes=outcomes, index=index))
      } else {
        dct <- makeDCT(high_pass=1/50, frame_times=1:nrow(tmp))
        polyX <- poly(1:nrow(tmp), order=10)
        polyX <- apply(polyX, 2, function(x) x/max(x)) # normalize to max 1 for each column
        return(list(outcomes=outcomes, index=index, dct=dct, poly=polyX))
      }
    }
  } else {
    return(list(stimulus=stimulus))
  }
}

## very funny Niek
# adapt.c.emc <- function(...){
#   return(NA)
# }

update_pars = function(s,npars,da,rfun=NULL,mapped_p=FALSE) #, #return_learning=FALSE,mapped_p=FALSE,
                       #return_all=FALSE)
  # for subject s
  # Either return da filling in responses and RT if da has no responses (in
  # in which case rfun must be supplied), or return npars filling in adapted
  # parameters or if return_learning returns learning (e.g., Q values)
{
  adapt <- attr(da,"adapt")$design

  if(adapt$useC & !mapped_p & !any(is.na(da$R))) {
    outcomes <- attr(da, 'adapt')[[s]]$outcomes
    index <- attr(da, 'adapt')[[s]]$index
    npars <- npars[da$subjects==s,]
    da <- da[da$subjects==s,]

    # update
    if(!is.null(adapt$stimulus)) {
      ###### RL ONLY -- IGNORE FOR DYNAMIC
      if(adapt$stimulus$adapt_fun_name=='delta') {
        learningRates <- matrix(npars[1,'alpha'], nrow=nrow(outcomes), ncol=ncol(outcomes))  # TODO make this more flexible
        startValues <- rep(npars[1,'q0'], ncol(outcomes))  # TODO make this more flexible
        updated <- adapt.c.emc(feedback=as.matrix(outcomes),
                               arguments=list(startValues = startValues,
                                              learningRates = learningRates),
                               learningRule='delta')
        colnames(updated$adaptedValues) <- colnames(outcomes)
        allQs <- updated$adaptedValues
      } else if(adapt$stimulus$adapt_fun_name=='vkf') {
        volatilityLearningRates <- matrix(npars[1,'alpha'], nrow=nrow(outcomes), ncol=ncol(outcomes))  # TODO make this more flexible
        predictionsStartValues <- rep(npars[1,'q0'], ncol(outcomes))  # TODO make this more flexible
        volatilitiesStartValues <- rep(npars[1,'volatility0'], ncol(outcomes))  # TODO make this more flexible
        uncertaintiesStartValues <- rep(npars[1,'w0'], ncol(outcomes))
        updated <- adapt.c.emc(feedback=as.matrix(outcomes),
                               arguments=list(volatilityLearningRates = volatilityLearningRates,
                                              predictionsStartValues = predictionsStartValues,
                                              volatilitiesStartValues = volatilitiesStartValues,
                                              uncertaintiesStartValues = uncertaintiesStartValues),
                               learningRule='vkf')
        allQs <- updated$adaptedPredictions
      } else if(adapt$stimulus$adapt_fun_name=='vkfbinary') {
        volatilityLearningRates <- matrix(npars[1,'alpha'], nrow=nrow(outcomes), ncol=ncol(outcomes))  # TODO make this more flexible
        predictionsStartValues <- rep(npars[1,'q0'], ncol(outcomes))  # TODO make this more flexible
        volatilitiesStartValues <- rep(npars[1,'volatility0'], ncol(outcomes))  # TODO make this more flexible
        uncertaintiesStartValues <- rep(npars[1,'w0'], ncol(outcomes))
        updated <- adapt.c.emc(feedback=as.matrix(outcomes),
                               arguments=list(volatilityLearningRates = volatilityLearningRates,
                                              predictionsStartValues = predictionsStartValues,
                                              volatilitiesStartValues = volatilitiesStartValues,
                                              uncertaintiesStartValues = uncertaintiesStartValues),
                               learningRule='vkf')
        allQs <- updated$adaptedPredictions
      }
      if(return_learning) {
        return(updated)
      }

      Q <- allQs[index]

      # Advantage framework (2AFC ONLY!!)
      if(('ws' %in% adapt$stimulus$output_par_names) & ('wd' %in% adapt$stimulus$output_par_names)) {
        indexOther <- attr(da, 'adapt')[[s]]$indexOther
        Q <- cbind(Q, allQs[indexOther])
      }

      ## function
      npars[,adapt$stimulus$output_name] <- adapt$stimulus$output_fun(npars[,adapt$stimulus$output_par_names], Q)
      ## hacky way of preventing errors due to extreme values
      npars[,'v'] <- pmin(npars[,'v'], 1e3)
      npars[,'v'] <- pmax(npars[,'v'], 0)

      ## add prediction errors and other latent learning variables!
      attr(npars, 'learn') <- updated
      return(npars)

    } else if(!is.null(adapt$dynamic)) {
      ###### DYNAMIC ONLY
      if(adapt$dynamic$adapt_fun_name=='delta') {
        # Delta rule for updating --------------------------------------------------------------
        da$lM <- as.logical(da$lM)
        learningRates <- as.matrix(cbind(npars[da$lM,'alpha1'], npars[da$lM,'alpha2'], npars[da$lM,'alpha3']))
        startValues <- npars[1, c('q01', 'q02', 'q03')]

        # if any of learning rates <= pnorm(-4) or >= pnorm(4), do not accept -- return parameters that are impossible so we get NAs
        alpha1 <- npars[1,c('alpha1', 'alpha2', 'alpha3')]
        nm <- colnames(npars)
#        alphamin <- npars[1,c('alpha1min', 'alpha2min', 'alpha3min')] # find corresponding minimum alpha values
        alphadiff <- alpha1#-alphamin
        if(any(alphadiff != 0)) {    # NB: this filters out learning rates that were set to a constant e.g., (0, 0.15)
          idx <- which(alphadiff != 0)
          if(any((alphadiff[idx] <= pnorm(-4)) | (alpha1[idx] >= pnorm(4)))) {
            npars[,'t0'] <- npars[,'t0']+3
            return(npars)
          }
        }

        # remove not-updated columns
        alphaIsZero <- round(learningRates[1,],6)==0
        learningRates <- learningRates[,!alphaIsZero]
        startValues <- startValues[!alphaIsZero]

        outcomes <- as.matrix(outcomes)
        if(ncol(outcomes)>1) outcomes <- as.matrix(outcomes[,!alphaIsZero[1:ncol(outcomes)]])

        if(!all(alphaIsZero)) {
          updated <- adapt.c.emc(feedback=outcomes,
                                 arguments=list(startValues = startValues,
                                                learningRates = learningRates),
                                 learningRule='delta')
        } else {
          #### Nothing to update, return npars!
          return(npars)
        }

        # construct updated Q-value matrix
        allQs <- matrix(NA, nrow=nrow(outcomes), ncol=3)
        if(sum(alphaIsZero)>0) allQs[,alphaIsZero] <- npars[da$lM, c('q01', 'q02', 'q03')][,alphaIsZero]    # re-add not-updated columns
        if(!all(alphaIsZero)) allQs[,!alphaIsZero] <- updated$adaptedValues

      } else if(adapt$dynamic$adapt_fun_name=='threshold') {
        # Threshold learning. ONLY uses Q-value 3 --------------------------------------------------------------
        da$lM <- as.logical(da$lM)
        learningRates <- as.matrix(cbind(npars[da$lM,'alpha1'], npars[da$lM,'alpha2'], npars[da$lM,'alpha3'])) # TODO make this more flexible
        startValues <- npars[1, c('q01', 'q02', 'q03')]

        # if any of learning rates <= pnorm(-4) or >= pnorm(4), do not accept -- return parameters that are impossible so we get NAs
        alpha1 <- npars[1,c('alpha1', 'alpha2', 'alpha3')]
        nm <- colnames(npars)
#        alphamin <- npars[1,c('alpha1min', 'alpha2min', 'alpha3min')] # grepl('alpha.*min', nm)]   # find corresponding minimum alpha values
        alphadiff <- alpha1#-alphamin
        if(any(alphadiff != 0)) {
          idx <- which(alphadiff != 0)
          if(any((alphadiff[idx] <= pnorm(-4)) | (alpha1[idx] >= pnorm(4)))) {
            npars[,'t0'] <- npars[,'t0']+3
            return(npars)
          }
        }

        # find 'base' threshold, order to trial order.
        b <- npars[order(da$trials),'B']

        # get mean threshold of both accumulators
        b <- apply(cbind(b[seq(1,length(b),2)], b[seq(2,length(b),2)]),1,mean)


        updated <- adaptb.c.emc(b0=b,
                                q0=startValues[3],
                                weight=npars[1,'weight3'],
                                learningRates=matrix(learningRates[,3], nrow=nrow(learningRates)),
                                rts=outcomes[,3])

        # reconstruct Q-value matrix
        allQs <- matrix(0, nrow=nrow(outcomes), ncol=3)
        allQs[,3] <- updated$adaptedValues[,2]           # Q3 from threshold updating
      } else if(adapt$dynamic$adapt_fun_name=='deltadeltathreshold') {
        # Combines 2 delta rules with threshold learning. --------------------------------------------------------------
        da$lM <- as.logical(da$lM)
        learningRates <- as.matrix(cbind(npars[da$lM,'alpha1'], npars[da$lM,'alpha2'], npars[da$lM,'alpha3']))
        startValues <- npars[1, c('q01', 'q02', 'q03')]

        # if any of learning rates <= pnorm(-4) or >= pnorm(4), do not accept -- return parameters that are impossible so we get NAs
        alpha1 <- npars[1,c('alpha1', 'alpha2', 'alpha3')]
        nm <- colnames(npars)
#        alphamin <- npars[1,c('alpha1min', 'alpha2min', 'alpha3min')] # grepl('alpha.*min', nm)]   # find corresponding minimum alpha values
        alphadiff <- alpha1#-alphamin
        if(any(alphadiff != 0)) {
          idx <- which(alphadiff != 0)
          if(any((alphadiff[idx] <= pnorm(-4)) | (alpha1[idx] >= pnorm(4)))) {
            return(npars)
          }
        }

        # update threshold FIRST
        # find 'base' threshold, order to trial order.
        b <- npars[order(da$trials),'B']
        # get mean threshold of both accumulators
        b <- apply(cbind(b[seq(1,length(b),2)], b[seq(2,length(b),2)]),1,mean)

        outcomes <- as.matrix(outcomes)

        ## new
        out = emcAdapt::adaptdeltab.c.emc(feedback=outcomes[,1:2],
                                          startValues=startValues[1:2],
                                          learningRates=learningRates[,1:2],

                                          b0=b,
                                          q0=startValues[3],
                                          weight=npars[1,'weight3'],
                                          learningRatesB=matrix(learningRates[,3], nrow=nrow(learningRates)),
                                          rts=outcomes[,3])
        allQs <- matrix(NA, nrow=nrow(outcomes), ncol=3)
        allQs[,1:2] <- out$adaptedValues[,1:2]  # Q1, Q2 from delta rule
        allQs[,3] <- out$adaptedValuesB[,2]     # Q3 from threshold updating
      }

      # apply output function ('control')
      if('collapse_across_SM' %in% names(adapt)) {
          if('stim_identifier' %in% colnames(da)) {
          # Ugly way to identify dataset 4, for which we need to collapse Q_SM within blocks
          daQsm <- cbind(da[da$lR=='word',], Qsm=allQs[,1])
          QsmByBlock <- aggregate(Qsm~block, daQsm, mean)
          for(block in unique(QsmByBlock$block)) allQs[daQsm$block == block,1] <- QsmByBlock[QsmByBlock$block == block,'Qsm']
        } else {
          allQs[,1] <- mean(allQs[,1])
        }
      }
      if('collapse_across_AM' %in% names(adapt)) allQs[,2] <- mean(allQs[,2])
      if('collapse_across_FM' %in% names(adapt)) allQs[,3] <- mean(allQs[,3])


      npars[,adapt$dynamic$output_name] <- apply_mechanisms(npars, allQs[index,], da, learningModel=adapt$dynamic$learningModel)
      npars[,c('q01', 'q02', 'q03')] <- allQs[index,]
      npars[is.na(npars)] <- 0  # remove NA values -- these shouldn't be here, but sometimes are and cause errors in likelihood calculations
      return(npars)
    }
  } else if(adapt$useSMsApproach) {
    ## This is for generating Posterior Predictives!
    if(!is.null(adapt$dynamic)) {
      outcomes <- attr(da, 'adapt')[[s]]$outcomes
      index <- attr(da, 'adapt')[[s]]$index  # index goes from trial-order to da-order
      npars <- npars[da$subjects==s,]
      da <- da[da$subjects==s,]

      ## Ugh ugly!! Do DCT transform here (if applicable)
      npars <- DCTtransform(npars=npars, da=da, s)
      ## potentially, polynomial basis functions
      npars <- polytransform(npars=npars, da=da, s)
      ## potentially, exponential trends
      npars <- trendTransform(npars=npars, da=da, s)


      ## add some columns that could be useful for updating
      da$R <- factor(da$R, levels=levels(da$lR))

      ##
      startValues <- npars[1, c('q01', 'q02', 'q03')]

      ## add response
      add_response <- any(is.na(da$R))
      if(is.null(rfun)) rfun <- attr(da, 'model')()$rfun

      if(adapt$dynamic$adapt_fun_name != 'delta2LR') {
        learningRates <- as.matrix(cbind(npars[,'alpha1'], npars[,'alpha2'], npars[,'alpha3']))  # TODO make this more flexible
      } else {
        learningRatesPos <- as.matrix(cbind(npars[,'alpha1Pos'], npars[,'alpha2Pos'], npars[,'alpha3Pos']))
        learningRatesNeg <- as.matrix(cbind(npars[,'alpha1Neg'], npars[,'alpha2Neg'], npars[,'alpha3Neg']))
        learningRatesNeg[learningRatesNeg==-1] <- learningRatesPos[learningRatesNeg==-1]
      }

      for(trial in 1:nrow(outcomes)) {
        Ri <- da$trials==trial    # Ri goes from da-order to trial-order

        # what is Q?
        if(trial == 1) {
          allQs <- startValues
        }

        # update npars
        npars[Ri,adapt$dynamic$output_name] <- apply_mechanisms(npars[Ri,],
                                                                matrix(allQs, nrow=nrow(npars[Ri,]), ncol=3, byrow=TRUE),
                                                                da[Ri,],
                                                                learningModel = attr(da, "adapt")$design$dynamic$learningModel)
        npars[Ri,c('q01', 'q02', 'q03')] <- matrix(allQs, nrow=nrow(npars[Ri,]), ncol=3, byrow=TRUE)

        # simulate trial
        Rrt <- rfun(da[Ri,'lR'], npars[Ri,])
        da[Ri, 'rt'] <- Rrt[,'rt']
        da[Ri, 'R'] <- Rrt[,'R']

        da[Ri, 'choice'] <- Rrt[,'R']==levels(da$R)[[1]]
        da[Ri, 'accuracy'] <- as.character(da[Ri, 'S']) == as.character(da[Ri,'R'])

        # update Q-values, npars
        if(trial < nrow(outcomes)) {
          feedback <- as.numeric(da[Ri,adapt$dynamic$columns_to_include][1,])
          if(adapt$dynamic$adapt_fun_name %in% c('threshold', 'deltadeltathreshold')) {
            feedback[3] <- mean(npars[Ri,'B'])/da[Ri,'rt'][1]   # estimated difficulty/drift
          }

          # we're assuming delta here
          for(ii in 1:length(allQs)) {
            PE <- feedback[ii]-allQs[ii]
            if(is.na(PE)) next  # no update, continue
            if(adapt$dynamic$adapt_fun_name%in%c('delta', 'threshold', 'deltadeltathreshold')) {
              allQs[ii] <- allQs[ii] + learningRates[trial,ii]*PE
            } else {
              if(PE < 0) {
                allQs[ii] <- allQs[ii] + learningRatesNeg[trial,ii]*PE
              } else {
                allQs[ii] <- allQs[ii] + learningRatesPos[trial,ii]*PE
              }
            }
          }
        }
      }
      attr(da, 'npars') <- npars
      return(da)
    }
  }
}

