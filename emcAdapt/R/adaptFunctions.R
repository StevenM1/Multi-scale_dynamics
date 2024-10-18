adapt.c.emc <- function(feedback, arguments, learningRule='delta') {
  nTrials <- nrow(feedback)
  nAdapt <- ncol(feedback)

  if(learningRule == 'delta') {
    # Simple Delta rule
    startValues <- arguments$startValues
    learningRates <- arguments$learningRates

    ## empty output arrays
    adaptedValues <- predictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)
    out = .C('adaptDelta',
             nTrials=nTrials,
             nChoices=nAdapt,
             values=as.double(startValues),
             adaptedValues=as.double(adaptedValues),
             predictionErrors=as.double(predictionErrors),
             outcomes=as.double(feedback),
             learningRates=as.double(learningRates),
             NAOK=TRUE)

    adaptedValues <- matrix(out$adaptedValues, nrow=nTrials, ncol=nAdapt)
    predictionErrors <- matrix(out$predictionErrors, nrow=nTrials, ncol=nAdapt)
    return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors))
  } else if(learningRule == 'delta2LR') {
    # Simple Delta rule
    startValues <- arguments$startValues
    learningRatesPos <- arguments$learningRatesPos
    learningRatesNeg <- arguments$learningRatesNeg

    ## empty output arrays
    adaptedValues <- predictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)
    out = .C('adaptDelta2LR',
             nTrials=nTrials,
             nChoices=nAdapt,
             values=as.double(startValues),
             adaptedValues=as.double(adaptedValues),
             predictionErrors=as.double(predictionErrors),
             outcomes=as.double(feedback),
             learningRatesPos=as.double(learningRatesPos),
             learningRatesNeg=as.double(learningRatesNeg),
             NAOK=TRUE)

    adaptedValues <- matrix(out$adaptedValues, nrow=nTrials, ncol=nAdapt)
    predictionErrors <- matrix(out$predictionErrors, nrow=nTrials, ncol=nAdapt)
    return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors))
  } else if(learningRule %in% c('vkf', 'vkf2', 'vkfbinary', 'vkfbinary2')) {
    # Variational Kalman filter
    predictionsStartValues <- arguments$predictionsStartValues
    volatilitiesStartValues <- arguments$volatilitiesStartValues
    volatilityLearningRates <- arguments$volatilityLearningRates
    uncertaintiesStartValues <- arguments$uncertaintiesStartValues

    if(learningRule == 'vkfbinary') {
      ccall <- 'adaptVKFbinary'
    } else if(learningRule == 'vkfbinary2') {
      ccall <- 'adaptVKFbinary2'
    } else if(learningRule == 'vkf') {
      ccall <- 'adaptVKF'
    } else if(learningRule == 'vkf2') {
      ccall <- 'adaptVKF2'
    } else {
      stop('learningRule not understood')
    }

    ## empty output arrays
    predictions <- predictionErrors <- learningRates <- volatilities <- volatilityPredictionErrors <- uncertainties <- matrix(nrow=nTrials, ncol=nAdapt)

    out = .C(ccall,
             nTrials=nTrials,
             nChoices=nAdapt,
             # Predictions
             predictions=as.double(predictionsStartValues),
             adaptedPredictions=as.double(predictions),
             predictionErrors=as.double(predictionErrors),
             learningRates=as.double(learningRates),
             # volatility
             volatilities=as.double(volatilitiesStartValues),
             adaptedVolatilities=as.double(volatilities),
             volatilityPredictionErrors=as.double(volatilityPredictionErrors),
             volatilityLearningRates=as.double(volatilityLearningRates),
             # uncertainty
             uncertainties=as.double(uncertaintiesStartValues),
             adaptedUncertainties=as.double(uncertainties),
             # feedback
             outcomes=as.double(feedback),
             NAOK=TRUE)

    adaptedPredictions <- matrix(out$adaptedPredictions, nrow=nTrials, ncol=nAdapt)
    predictionErrors <- matrix(out$predictionErrors, nrow=nTrials, ncol=nAdapt)

    adaptedVolatilities <- matrix(out$adaptedVolatilities, nrow=nTrials, ncol=nAdapt)
    volatilityPredictionErrors <- matrix(out$volatilityPredictionErrors, nrow=nTrials, ncol=nAdapt)

    learningRates <- matrix(out$learningRates, nrow=nTrials, ncol=nAdapt)
    adaptedUncertainties <- matrix(out$adaptedUncertainties, nrow=nTrials, ncol=nAdapt)

    return(list(adaptedPredictions=adaptedPredictions, predictionErrors=predictionErrors, learningRates=learningRates,
                adaptedVolatilities=adaptedVolatilities, volatilityPredictionErrors=volatilityPredictionErrors, adaptedUncertainties=adaptedUncertainties))
  }
}

adapt.r.test <- function(startValues, learningRates, feedback, learningRule='delta',
                         learningRatesNeg=NULL, riskStartValues=NULL, riskLearningRates=NULL) {
  nTrials <- nrow(feedback)
  nAdapt <- ncol(feedback)

  # declare output array
  adaptedValues <- predictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)

  if(learningRule == 'SARSARisk') {
    adaptedRiskValues <- riskPredictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)
    riskValues <- riskStartValues
  }

  values <- startValues
  # Update all values (can we vectorize this? Hopefully)
  for(trial in 1:nTrials) {
    c_ = which(!is.na(feedback[trial,])) #choice[trial] # c_ = choice
    o_ = feedback[trial,c_] # o_ = outcome

    # Keep track of current value, used later for drift rates
    adaptedValues[trial, c_] = values[c_]
    adaptedValues[trial, -c_] = values[-c_]

    if(learningRule == 'delta' | learningRule == 'SARSAvarLR' | learningRule == 'SARSARisk') {
      predicted <- values[c_]
      if(learningRule == 'SARSARisk') {
        adaptedRiskValues[trial, c_] = riskValues[c_]
        adaptedRiskValues[trial, -c_] = riskValues[-c_]
      }
    } else if(learningRule == 'Qlearning')  {
      predicted <- max(values)
    }
    # calculate PE
    dv = o_-predicted  # prediction error = outcome (reward) - predicted value
    predictionErrors[trial,c_] <- dv  # keep track of this

    if(!is.null(learningRatesNeg)) {
      LR = ifelse(dv>0, learningRates[trial,c_], learningRatesNeg[trial,c_])
    } else if(learningRule == 'SARSAvarLR') {
      LR = learningRates[trial,c_] * (predicted*(1-predicted))
    } else if (learningRule == 'delta') {
      LR = learningRates[trial,c_]
    } else if(learningRule == 'SARSARisk') {
      LR = learningRates[trial,c_]
      riskLR = riskLearningRates[trial, c_]
      # update risk
      riskPE <- dv^2 - riskValues[c_]
      riskValues[c_] = riskValues[c_]+riskLR*riskPE
      riskPredictionErrors[trial,c_] = riskPE
      dv <- dv /sqrt(riskValues[c_])  # scale PE
    }
    values[c_] = values[c_] + LR*dv

    # # update values
    # if(!is.null(learningRatesNeg)) {
    #   values[c_] = values[c_] + ifelse(dv>0, learningRates[trial,c_], learningRatesNeg[trial,c_])*dv
    # } else {
    #   if(learningRule == 'SARSAvarLR') {
    #     LR
    #   }
    #   values[c_] = values[c_] + learningRates[trial,c_]*dv
    # }
  }
  if(learningRule == 'SARSARisk') {
    return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors,
           adaptedRiskValues=adaptedRiskValues, riskPredictionErrors=riskPredictionErrors))
  } else {
    return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors))
  }
}

vkf.r <- function(x, data, sample = FALSE, updateOnly=TRUE, implementationNumber=1) {
  #       data: column-vector of outcomes
  #       0<lambda<1, volatility learning rate
  #       v0>0, initial volatility
  #       sigma2>0, outcome noise
  #       sigmar = simulation response noise (sd)
  #       mu0 = initial value of mu (assumed to be 0 in code below)
  # also see https://github.com/payampiray/VKF/blob/master/vkf.m

  if (is.null(names(x))) {
    names(x) <- c("lambda", "v0", "sigma2", "sigmar", "mu0")
  }

  # extract parameters
  lambda <- x[["lambda"]]
  v0 <- x[["v0"]]
  sigma2 <- x[["sigma2"]]

  # number of trials
  nt <- nrow(data)

  # initial values
  m <- 300 * x[["mu0"]]
  w0 <- x[["sigma2"]]
  w <- w0
  v <- x[["v0"]]

  predictions <- numeric(nt)
  learning_rate <- numeric(nt)
  volatility <- numeric(nt)
  prediction_error <- numeric(nt)
  volatility_error <- numeric(nt)
  uncertainty <- numeric(nt)

  for (t in seq_len(nt)) {

    o <- data$o[t]
    predictions[t] <- m
    volatility[t] <- v
    uncertainty[t] <- w

    mpre <- m
    wpre <- w

    delta_m <- o - m
    k       <- (w + v) / (w + v + sigma2)                         # Eq 9
    m       <- m + k * delta_m                                    # Eq 10
    w       <- (1 - k) * (w + v)                                  # Eq 11

    if(implementationNumber == 1) {
      wcov    <-  (1 - k) * wpre                                    # Eq 12
      delta_v <-  (m - mpre)^2 + w + wpre - 2 * wcov - v
    } else {
      delta_v = k^2 * delta_m^2 + k*wpre - k*v
    }
    v       <-  v + lambda * delta_v                              # Eq 13

    learning_rate[t] <- k
    prediction_error[t] <- delta_m
    volatility_error[t] <- delta_v
  }

  if(updateOnly) {
    return(list(predictions = predictions,
                volatility = volatility,
                learning_rate = learning_rate,
                prediction_error = prediction_error,
                volatility_error = volatility_error,
                uncertainty = uncertainty))
  }
  # response model
  if (sample) {

    # generate responses
    data$r <- rnorm(nt, mean = predictions, sd = x[["sigmar"]])
    attr(data, "latent_state_pars") <- list(predictions = predictions,
                                            volatility = volatility,
                                            learning_rate = learning_rate,
                                            prediction_error = prediction_error,
                                            volatility_error = volatility_error,
                                            uncertainty = uncertainty)
    return(data)


  } else {

    out <- sum(dnorm(data$r,
                     mean = predictions,
                     sd = x[["sigmar"]],
                     log = TRUE))
    attr(data, "latent_state_pars") <- list(predictions = predictions,
                                            volatility = volatility,
                                            learning_rate = learning_rate,
                                            prediction_error = prediction_error,
                                            volatility_error = volatility_error,
                                            uncertainty = uncertainty)
    return(out)

  }
}

vkf_binary.r <- function(x, data, sample = FALSE, updateOnly=TRUE, implementationNumber=1) {
  #       data: column-vector of outcomes
  #       0<lambda<1, volatility learning rate
  #       v0>0, initial volatility
  #       sigma2>0, outcome noise
  #       sigmar = simulation response noise (sd)
  #       mu0 = initial value of mu (assumed to be 0 in code below)
  # also see https://github.com/payampiray/VKF/blob/master/vkf.m

  if (is.null(names(x))) {
    names(x) <- c("lambda", "v0", "omega", "sigmar", "mu0")
  }

  # extract parameters
  lambda <- x[["lambda"]]
  v0 <- x[["v0"]]
  omega <- x[["omega"]]

  # number of trials
  nt <- nrow(data)

  # initial values
  m <- 300 * x[["mu0"]]
  w0 <- x[["omega"]]
  w <- w0
  v <- x[["v0"]]

  predictions <- numeric(nt)
  learning_rate <- numeric(nt)
  volatility <- numeric(nt)
  prediction_error <- numeric(nt)
  volatility_error <- numeric(nt)
  uncertainty <- numeric(nt)

  for (t in seq_len(nt)) {

    o <- data$o[t]
    predictions[t] <- m
    volatility[t] <- v
    uncertainty[t] <- w

    mpre <- m
    wpre <- w

    k       <- (w + v) / (w + v + omega)                          # Eq 14
    alpha   <- sqrt(w + v)                                        # Eq 15
    delta_m <- o - (1/(1+exp(-m)))                                # Eq 16 (sigmoid)
    m       <- m + alpha * delta_m                                # Eq 16
    w       <- (1 - k) * (w + v)                                  # Eq 17

    if(implementationNumber == 1) {
      wcov    <-  (1 - k) * wpre                                    # Eq 18
      delta_v <-  (m - mpre)^2 + w + wpre - 2 * wcov - v            # Eq 19.1
    } else {
      delta_v = delta_m^2 + k*wpre - k*v
      }
    v       <-  v + lambda * delta_v                              # Eq 19.2

    learning_rate[t] <- alpha                                     # learning rate is now alpha!!
    prediction_error[t] <- delta_m
    volatility_error[t] <- delta_v
  }

  if(updateOnly) {
    return(list(predictions = predictions,
                volatility = volatility,
                learning_rate = learning_rate,
                prediction_error = prediction_error,
                volatility_error = volatility_error,
                uncertainty = uncertainty))
  }
  # response model
  if (sample) {

    # generate responses
    data$r <- rnorm(nt, mean = predictions, sd = x[["sigmar"]])
    attr(data, "latent_state_pars") <- list(predictions = predictions,
                                            volatility = volatility,
                                            learning_rate = learning_rate,
                                            prediction_error = prediction_error,
                                            volatility_error = volatility_error,
                                            uncertainty = uncertainty)
    return(data)


  } else {

    out <- sum(dnorm(data$r,
                     mean = predictions,
                     sd = x[["sigmar"]],
                     log = TRUE))
    attr(data, "latent_state_pars") <- list(predictions = predictions,
                                            volatility = volatility,
                                            learning_rate = learning_rate,
                                            prediction_error = prediction_error,
                                            volatility_error = volatility_error,
                                            uncertainty = uncertainty)
    return(out)

  }
}

vkf.r.multicolumns <- function(x, feedback, sample = FALSE, updateOnly=TRUE, implementationNumber=1) {
  #       data: column-vector of outcomes
  #       0<lambda<1, volatility learning rate
  #       v0>0, initial volatility
  #       sigma2>0, outcome noise
  #       sigmar = simulation response noise (sd)
  #       mu0 = initial value of mu (assumed to be 0 in code below)
  # also see https://github.com/payampiray/VKF/blob/master/vkf.m

  if (is.null(names(x))) {
    names(x) <- c("lambda", "v0", "sigma2", "sigmar", "mu0")
  }

  # extract parameters
  lambda <- x[["lambda"]]
  v0 <- x[["v0"]]
  sigma2 <- x[["sigma2"]]

  # number of trials
  nt <- nrow(feedback)
  nc <- ncol(feedback)

  # initial values
  m <- 300 * x[["mu0"]]
  w0 <- x[["sigma2"]]
  w <- w0
  v <- x[["v0"]]

  predictions <- matrix(0, nrow=nt, ncol=nc)
  learning_rate <- matrix(NA, nrow=nt, ncol=nc)
  volatility <- matrix(v0, nrow=nt, ncol=nc)
  prediction_error <- matrix(NA, nrow=nt, ncol=nc)
  volatility_error <- matrix(NA, nrow=nt, ncol=nc)
  uncertainty <- matrix(w0, nrow=nt, ncol=nc)

  for(column in seq_len(nc)) {
    m <- 300 * x[["mu0"]]
    w0 <- x[["sigma2"]]
    w <- w0
    v <- x[["v0"]]
    for (t in seq_len(nt)) {

      o <- feedback[t,column]
      predictions[t,column] <- m
      volatility[t,column] <- v
      uncertainty[t,column] <- w

      mpre <- m
      wpre <- w

      delta_m <- o - m
      k       <- (w + v) / (w + v + sigma2)                         # Eq 9
      m       <- m + k * delta_m                                    # Eq 10
      w       <- (1 - k) * (w + v)                                  # Eq 11

      if(implementationNumber == 1) {
        wcov    <-  (1 - k) * wpre                                    # Eq 12
        delta_v <-  (m - mpre)^2 + w + wpre - 2 * wcov - v
      } else {
        delta_v = (k^2) * (delta_m^2) + k*wpre - k*v
      }
      v       <-  v + lambda * delta_v                              # Eq 13

      learning_rate[t,column] <- k
      prediction_error[t,column] <- delta_m
      volatility_error[t,column] <- delta_v

    }
  }

  if(updateOnly) {
    return(list(predictions = predictions,
                volatility = volatility,
                learning_rate = learning_rate,
                prediction_error = prediction_error,
                volatility_error = volatility_error,
                uncertainty = uncertainty))
  }
}

vkf_binary.r.multicolumns <- function(x, feedback, sample = FALSE, updateOnly=TRUE, implementationNumber=1) {
  #       data: column-vector of outcomes
  #       0<lambda<1, volatility learning rate
  #       v0>0, initial volatility
  #       sigma2>0, outcome noise
  #       sigmar = simulation response noise (sd)
  #       mu0 = initial value of mu (assumed to be 0 in code below)
  # also see https://github.com/payampiray/VKF/blob/master/vkf.m

  if (is.null(names(x))) {
    names(x) <- c("lambda", "v0", "omega", "sigmar", "mu0")
  }

  # extract parameters
  lambda <- x[["lambda"]]
  v0 <- x[["v0"]]
  omega <- x[["omega"]]

  # number of trials
  nt <- nrow(feedback)
  nc <- ncol(feedback)

  # initial values
  m <- 300 * x[["mu0"]]
  w0 <- x[["omega"]]
  w <- w0
  v <- x[["v0"]]

  predictions <- matrix(0, nrow=nt, ncol=nc)
  learning_rate <- matrix(NA, nrow=nt, ncol=nc)
  volatility <- matrix(v0, nrow=nt, ncol=nc)
  prediction_error <- matrix(NA, nrow=nt, ncol=nc)
  volatility_error <- matrix(NA, nrow=nt, ncol=nc)
  uncertainty <- matrix(w0, nrow=nt, ncol=nc)

  for(column in seq_len(nc)) {
    m <- 300 * x[["mu0"]]
    w0 <- x[["omega"]]
    w <- w0
    v <- x[["v0"]]

    for (t in seq_len(nt)) {

      o <- feedback[t,column]
      predictions[t,column] <- m
      volatility[t,column] <- v
      uncertainty[t,column] <- w

      mpre <- m
      wpre <- w

      k       <- (w + v) / (w + v + omega)                          # Eq 14
      alpha   <- sqrt(w + v)                                        # Eq 15
      delta_m <- o - (1/(1+exp(-m)))                                # Eq 16 (sigmoid)
      m       <- m + alpha * delta_m                                # Eq 16
      w       <- (1 - k) * (w + v)                                  # Eq 17

      if(implementationNumber == 1) {
        wcov    <-  (1 - k) * wpre                                    # Eq 18
        delta_v <-  (m - mpre)^2 + w + wpre - 2 * wcov - v            # Eq 19.1
      } else {
        delta_v = (alpha*delta_m)^2 + k*wpre - k*v
      }
      v       <-  v + lambda * delta_v                              # Eq 19.2

      learning_rate[t,column] <- alpha                              # learning rate is now alpha!!
      prediction_error[t,column] <- delta_m
      volatility_error[t,column] <- delta_v
    }
  }
  if(updateOnly) {
    return(list(predictions = predictions,
                volatility = volatility,
                learning_rate = learning_rate,
                prediction_error = prediction_error,
                volatility_error = volatility_error,
                uncertainty = uncertainty))
  }
}

#


adaptb.r.test <- function(b0, q0, rts, weight, learningRates) {
  nTrials <- length(rts)
  nAdapt <- 2

  # declare output array
  values <- c(b0[1], q0)

  adaptedValues <- predictionErrors <- matrix(nrow=nTrials, ncol=nAdapt)
  for(trial in 1:nTrials) {
    # if(trial > 1) {
    #   # store threshold and Q-value of last trial
    # }

    rt = rts[trial]

    # values[2] = Q-value3 = running average of drift rates
    predicted <- values[2]
    b <- b0[trial] + weight*predicted     # update threshold
    dv <- (b/rt) - predicted   # was the drift on last trial lower or higher than running average of drifts? PE
    values[1] <- b

    # update
    predictionErrors[trial,2] <- dv

    # store
    adaptedValues[trial,] <- values  # nb: values[1] has been updated here, values[2] not yet -- because we're looking at *last* trial

    # update Q-value for new trial
    values[2] <- values[2] + learningRates[trial,2]*dv

  }
  return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors, b0=b0))
}

adaptb.c.emc <- function(b0, q0, weight, learningRates, rts) {
  nTrials <- length(rts)

  startValues <- c(b0[1], q0)
  #learningRates <- arguments$learningRates

  ## empty output arrays
  adaptedValues <- predictionErrors <- matrix(nrow=nTrials, ncol=2)
  out = .C('adaptThresholds',
           nTrials=nTrials,
           b0=b0,
           weight=weight,
           values=as.double(startValues),
           adaptedValues=as.double(adaptedValues),
           predictionErrors=as.double(predictionErrors),
           rts=as.double(rts),
           learningRates=as.double(learningRates),
           NAOK=TRUE)

  adaptedValues <- matrix(out$adaptedValues, nrow=nTrials, ncol=2)
  predictionErrors <- matrix(out$predictionErrors, nrow=nTrials, ncol=2)
  return(list(adaptedValues=adaptedValues, predictionErrors=predictionErrors))
}


#### FOR TESTING ONLY
# simulate_dynEAM <- function(#v, v_diff, b, t0, A,
#   v1, v2, b, t0, A,
#   weight1, weight2, weight3,
#   alpha1, alpha2, alpha3,
#   q01, q02, q03,
#   s1=1, s2=1,
#   stimuli, nTrials) {
#
#   allPars <- cbind(v1, v2, b, t0, A, weight1, weight2, weight3, alpha1, alpha2, alpha3, q01, q02, q03, s1, s2, PE=NA)
#   if(nrow(allPars) == 1) {
#     allPars <- matrix(rep(allPars, each=nTrials), nTrials, byrow=FALSE)
#   }
#   colnames(allPars) <- c('v1', 'v2', 'B', 't0', 'A', 'weight1', 'weight2', 'weight3', 'alpha1', 'alpha2', 'alpha3', 'q1', 'q2', 'q3', 's1', 's2', 'PE')
#
#   # data array (output)
#   rtdata <- array(NA, dim=c(nTrials, 2), dimnames=list(1:nTrials, c('R', 'rt')))
#
#   # loop over trials
#   for(trial in 1:nTrials) {
#     if(trial > 1) {
#       allPars[trial,'q1'] <- allPars[trial-1,'q1'] + allPars[trial-1,'alpha1']*(stimuli[trial-1] - allPars[trial-1,'q1'])
#       allPars[trial,'q2'] <- allPars[trial-1,'q2'] + allPars[trial-1,'alpha2']*((Rrt$R==1) - allPars[trial-1,'q2'])
#       # dv <- allPars[trial-1,'B']/Rrt$rt - allPars[trial-1,'q3']
#       allPars[trial,'q3'] <- allPars[trial-1,'q3'] + allPars[trial-1,'alpha3']*(dv)
#     }
#
#     ## Get parameters of trial N of all iters
#     #pars <- matrix(aperm(allPars[,,trial,], c(2,1,3)), nrow=nIter*2, ncol=length(parNames), byrow=TRUE)
#     #colnames(pars) <- dimnames(allPars)[[2]]
#     pars <- rbind(allPars[trial,c('v1', 'B', 't0', 'A', 's1', 'weight1', 'weight2', 'weight3', 'q1', 'q2', 'q3')],
#                   allPars[trial,c('v2', 'B', 't0', 'A', 's2', 'weight1', 'weight2', 'weight3', 'q1', 'q2', 'q3')])
#     colnames(pars) <- c('v', 'B', 't0', 'A', 's', 'weight1', 'weight2', 'weight3', 'q1', 'q2', 'q3')
#
#     ## dynamic component
#     pars[,'B'] <- pars[,'B'] + pars[,'weight3']*(pars[,'q3'])
#     pars[,'v'] <- pars[,'v'] + pars[,'weight2']*(1-pars[,'q2'])
#     #pars[,'B'] <- pars[,'B'] + pars[,'weight1']*rep(c(1,-1), nTrials)*pars[,'q1']
#     pars[,'B'] <- pars[,'B'] + pars[,'weight1']*c(1,-1)*pars[,'q1']
#     pars[,'B'] <- pmax(pars[,'B'], 0.01)
#
#     #print(pars[,'B'])
#     ## simulate
#     Rrt <- EMC2:::rRDM(factor(levels=c(1,2)), pars)
#     Rrt$R <- as.numeric(Rrt$R)
#
#     # store data
#     rtdata[trial,] <- as.matrix(Rrt)
#     allPars[trial,'B'] <- mean(pars[,'B'])
#     allPars[trial,'v1'] <- pars[1,'v']
#     allPars[trial,'v2'] <- pars[2,'v']
#     dv <- allPars[trial,'B']/Rrt$rt - allPars[trial,'q3']
#     allPars[trial,'PE'] <- dv
#
#   }
#
#   return(list(rtdata, allPars))
# }
#
# v <- 1
# v_diff <- .5
# b <- 3
# t0 <- 0.15
# weight1 <- 0.0
# alpha1 <- 0.3
# weight2 <- 0   ## accuracy history weight
# alpha2 <- 0.3    ## learning rate for accuracy history
# weight3 <- -.7   ## rate history weight
# alpha3 <- 0.5    ## learning rate for rate history
#
# q01 <- 0         ## initial Q-value (stimulus)
# q02 <- 1         ## initial Q-value (accuracy)
# q03 <- 3/2       ## initial Q-value (rate)
#
# nTrials <- 1e2
# stimuli <- sample(c(0,1), size=nTrials, replace = TRUE)
#
# b <- c(rep(2, nTrials/2), rep(3, nTrials/2))
# #b[sample(1:nTrials, nTrials/2)] <- 1
# #alpha3=0
#
# #b <- c(2, 2, rep(3, 200))
# #debugonce(simulate_dynEAM)
# out <- simulate_dynEAM(v1=v+v_diff, v2=v-v_diff, A=0, b=b, t0=t0, weight1=weight1, weight2=weight2, weight3=weight3,
#                        alpha1=alpha1, alpha2=alpha2, alpha3=alpha3, q01=q01, q02=q02, q03=q03, s1=1, s2=1,
#                        nTrials = nTrials, stimuli = stimuli)
# df <- data.frame(cbind(out[[1]], out[[2]][,c('v1', 'v2', 'B', 'q3', 'PE')]))
# df$estimated_v <- df$B/df$rt
#
# out1 = adaptb.c.emc(b0=b, q0=3/2, weight=weight3, learningRates=matrix(alpha3, nrow=nTrials, ncol=1), rts=df$rt)
# #debugonce(adaptb.r.test)
# out2 = adaptb.r.test(b0=b, q0=3/2, weight=weight3, learningRates = matrix(alpha3, nrow=nTrials, ncol=2), rts=df$rt)
#
# ## combine all useful data
# df = cbind(df, data.frame(PE.c=out1[[2]][,1],
#                            b.c=out1[[1]][,1],
#                            q3.c=out1[[1]][,2],
#                            PE.r=out2[[2]][,2],
#                            b.r=out2[[1]][,1],
#                            q3.r=out2[[1]][,2]))
# pairs(df[,c('PE.c', 'PE.r', 'PE')])
# pairs(df[,c('b.c', 'b.r', 'B')])
# pairs(df[,c('q3.c', 'q3.r', 'q3')])
# tail(df)



## testing 2LR
# nTrials <- 10000
# notChosen <- sample(c(1, 2), size=nTrials, replace=TRUE)
# eta1 <- .3  # positive LR
# eta2 <- 1  # negative LR
# startValues <- c(.5, .5)
# feedback <- matrix(sample(c(0,1), size=nTrials*2, replace=TRUE, prob = c(.1,.9)), nrow=nTrials)
#
# out = adapt.c.emc(feedback, arguments=list(startValues=c(0,0),
#                                            learningRatesPos=matrix(eta1, nrow=nTrials, ncol=2),
#                                            learningRatesNeg=matrix(eta2, nrow=nTrials, ncol=2)),
#                   learningRule='delta2LR')
# out
#
# update_function_for_mb <- function() {
#   nTrials <- 10000
#   notChosen <- sample(c(1, 2), size=nTrials, replace=TRUE)
#   eta1 <- .3  # positive LR
#   eta2 <- 1  # negative LR
#   startValues <- c(.5, .5)
#   feedback <- matrix(sample(c(0,1), size=nTrials*2, replace=TRUE, prob = c(.1,.9)), nrow=nTrials)
#
#   out = adapt.c.emc(feedback, arguments=list(startValues=c(0,0),
#                                              learningRatesPos=matrix(eta1, nrow=nTrials, ncol=2),
#                                              learningRatesNeg=matrix(eta2, nrow=nTrials, ncol=2)),
#                     learningRule='delta2LR')
#   return(out)
# }
#
# library(microbenchmark)
# microbenchmark(update_function_for_mb(), times = 5e3)
#


# feedback[cbind(1:nTrials, as.numeric(notChosen))] <- NA
# choice <- ifelse(notChosen==1, 2, 1)


# Testing VKF -------------------------------------------------------------
# feedback = matrix(rnorm(40), ncol=2)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# out1.1 = adapt.c.emc(feedback=feedback, learningRule='vkf',
#                   arguments=list(predictionsStartValues=predictionsStartValues,
#                                  volatilitiesStartValues=volatilitiesStartValues,
#                                  uncertaintiesStartValues=uncertaintiesStartValues,
#                                  volatilityLearningRates=volatilityLearningRates))
# out1.2 = adapt.c.emc(feedback=feedback, learningRule='vkf2',
#                    arguments=list(predictionsStartValues=predictionsStartValues,
#                                   volatilitiesStartValues=volatilitiesStartValues,
#                                   uncertaintiesStartValues=uncertaintiesStartValues,
#                                   volatilityLearningRates=volatilityLearningRates))
#
#
# # ## compare with Q's implementation
# column <- 2
# x <- c("lambda"=volatilityLearningRates[1,column],    # volatility learning rate
#        "v0"=volatilitiesStartValues[column],        # volatility start value
#        "sigma2"=uncertaintiesStartValues[column],    # uncertainty start value
#        "sigmar"=uncertaintiesStartValues[column],    # simulation noise, unused for updating
#        "mu0"=predictionsStartValues[column]      # prediction start value
#        )
#
# out2.1 <- vkf.r(x, data=data.frame(o=feedback[,column]), updateOnly=TRUE, implementationNumber = 1)
# out2.2 <- vkf.r(x, data=data.frame(o=feedback[,column]), updateOnly=TRUE, implementationNumber = 2)

#
# all(out2.1$predictions==out1.1$adaptedPredictions[,column])
# all(out2.1$volatility==out1.1$adaptedVolatilities[,column])
# all(out2.1$uncertainty==out1.1$adaptedUncertainties[,column])

# all(out2.2$predictions==out1.2$adaptedPredictions[,column])
# all(out2.2$volatility==out1.2$adaptedVolatilities[,column])
# all(out2.2$uncertainty==out1.2$adaptedUncertainties[,column])

#

## speed comparison
# feedback = matrix(rnorm(2000), ncol=1)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# library(microbenchmark)
# microbenchmark(adapt.c.emc(feedback=feedback, learningRule='vkf',
#                            arguments=list(predictionsStartValues=predictionsStartValues,
#                                           volatilitiesStartValues=volatilitiesStartValues,
#                                           uncertaintiesStartValues=uncertaintiesStartValues,
#                                           volatilityLearningRates=volatilityLearningRates)),
#                vkf.r(x, data=data.frame(o=feedback), updateOnly=TRUE))
# Unit: microseconds
#       min       lq      mean   median        uq      max neval
#   67.978   83.517  107.8521   92.824  103.9555 1274.075   100
# 1183.014 1215.014 1349.7270 1252.796 1352.3645 4342.310   100


# Testing BINARY VKF -------------------------------------------------------------
# feedback = matrix(rbinom(20, size=1, prob=.5), ncol=2)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# out1.1 = adapt.c.emc(feedback=feedback, learningRule='vkfbinary',
#                   arguments=list(predictionsStartValues=predictionsStartValues,
#                                  volatilitiesStartValues=volatilitiesStartValues,
#                                  uncertaintiesStartValues=uncertaintiesStartValues,
#                                  volatilityLearningRates=volatilityLearningRates))
#
# out1.2 = adapt.c.emc(feedback=feedback, learningRule='vkfbinary2',
#                   arguments=list(predictionsStartValues=predictionsStartValues,
#                                  volatilitiesStartValues=volatilitiesStartValues,
#                                  uncertaintiesStartValues=uncertaintiesStartValues,
#                                  volatilityLearningRates=volatilityLearningRates))
#
#
# ## compare with Q's implementation
# column <- 2
# x <- c("lambda"=volatilityLearningRates[1,column],    # volatility learning rate
#        "v0"=volatilitiesStartValues[column],        # volatility start value
#        "omega"=uncertaintiesStartValues[column],    # uncertainty start value
#        "sigmar"=uncertaintiesStartValues[column],    # simulation noise, unused for updating
#        "mu0"=predictionsStartValues[column]      # prediction start value
#        )
# out2.1 <- vkf_binary.r(x, data=data.frame(o=feedback[,column]), updateOnly=TRUE, implementationNumber = 1)
# out2.2 <- vkf_binary.r(x, data=data.frame(o=feedback[,column]), updateOnly=TRUE, implementationNumber = 2)
#
#
# all(out2.1$predictions==out1.1$adaptedPredictions[,column])
# all(out2.2$volatility==out1.2$adaptedVolatilities[,column])
# #
# all(out1.2$adaptedPredictions[,column]==out1.1$adaptedPredictions[,column])

#
# # speed comparison
# feedback = matrix(rnorm(2000), ncol=1)
# predictionsStartValues = rep(0, ncol(feedback))
# volatilitiesStartValues = rep(.1, ncol(feedback))
# uncertaintiesStartValues = rep(1, ncol(feedback))
# volatilityLearningRates = matrix(.1, nrow=nrow(feedback), ncol=ncol(feedback))
#
# library(microbenchmark)
# microbenchmark(adapt.c.emc(feedback=feedback, learningRule='vkfbinary',
#                            arguments=list(predictionsStartValues=predictionsStartValues,
#                                           volatilitiesStartValues=volatilitiesStartValues,
#                                           uncertaintiesStartValues=uncertaintiesStartValues,
#                                           volatilityLearningRates=volatilityLearningRates)),
#                vkf_binary.r(x, data=data.frame(o=feedback), updateOnly=TRUE))
# Unit: microseconds
# min        lq       mean   median        uq      max neval
# 75.973   89.7285   98.45125   97.088  103.3405  168.428   100
# 1346.276 1365.3820 1481.72811 1382.028 1410.5435 3725.588   100



# Testing delta rule ------------------------------------------------------
#
# nTrials <- 100
# notChosen <- sample(c(1, 2), size=nTrials, replace=TRUE)
# eta1 <- .3
# eta2 <- .4
# startValues <- c(.5, .5)
# feedback <- matrix(rnorm(nTrials*2, 5, 3), nrow=nTrials)
# feedback[cbind(1:nTrials, as.numeric(notChosen))] <- NA
# choice <- ifelse(notChosen==1, 2, 1)
# #
# tmpSarsaC = adapt.c.emc(feedback=feedback, arguments=list(startValues=startValues, learningRates=matrix(eta1, nrow=nTrials, ncol=2)), learningRule='delta')
# tmpSarsaR = adapt.r.test(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2), learningRule='delta')
#
# tmpSarsaC$adaptedValues == tmpSarsaR$adaptedValues
#
# #
# ## negative learning rates?
# posLR <- matrix(eta1, nrow=nTrials, ncol=2)
# negLR <- matrix(eta2, nrow=nTrials, ncol=2)
# tmpSarsaC = adapt.c.emc(feedback=feedback, startValues=startValues, posLR, learningRule='SARSA', negLR)
# tmpSarsaR = adapt.r.test(feedback=feedback, startValues=startValues, posLR, learningRule='SARSA', negLR)
#
#
# ## variable learning rates?
# tmpSarsaC = adapt.c.emc(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2), learningRule='SARSAvarLR')
# tmpSarsaR = adapt.r.test(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2), learningRule='SARSAvarLR')
#
# tmpSarsaC$adaptedValues == tmpSarsaR$adaptedValues
#
#
# ## Risk learning?
# riskStartValues <- c(1,1)
# tmpSarsaRiskC = adapt.c.emc(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2),
#                             learningRule='SARSARisk', riskStartValues=riskStartValues, riskLearningRates = matrix(eta1, nrow=nTrials, ncol=2))
# tmpSarsaRiskR = adapt.r.test(feedback=feedback, startValues=startValues, matrix(eta1, nrow=nTrials, ncol=2),
#                              learningRule='SARSARisk', riskStartValues=riskStartValues, riskLearningRates = matrix(eta1, nrow=nTrials, ncol=2))
# tmpSarsaRiskC$adaptedValues == tmpSarsaRiskR$adaptedValues
# tmpSarsaRiskC$adaptedRiskValues == tmpSarsaRiskR$adaptedRiskValues
# tmpSarsaRiskC$riskPredictionErrors == tmpSarsaRiskR$riskPredictionErrors
