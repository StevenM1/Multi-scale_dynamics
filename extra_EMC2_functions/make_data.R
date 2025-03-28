predict.emc <-function (object, hyper = FALSE, n_post = 50, n_cores = 1,
                        stat = c("random", "mean", "median")[1], ...)
{
  emc <- object
  dots <- list(...)
  data <- get_data(emc)
  design <- get_design(emc)

  ## SM added
  design[[1]]$adapt <- lapply(emc[[1]]$data, function(x) attr(x,'adapt'))
  if('conditionalOnData' %in% names(dots)) {
    if(dots$conditionalOnData) {
      for(subj in 1:length(design[[1]]$adapt)) design[[1]]$adapt[[subj]]$design$simulateByTrial <- FALSE
    } else {
      for(subj in 1:length(design[[1]]$adapt)) design[[1]]$adapt[[subj]]$design$simulateByTrial <- TRUE
    }
  }

  if('collapse_across_SM' %in% names(dots)) for(subj in 1:length(design[[1]]$adapt)) design[[1]]$adapt[[subj]]$design$collapse_across_SM <- TRUE
  if('collapse_across_AM' %in% names(dots)) for(subj in 1:length(design[[1]]$adapt)) design[[1]]$adapt[[subj]]$design$collapse_across_AM <- TRUE
  if('collapse_across_FM' %in% names(dots)) for(subj in 1:length(design[[1]]$adapt)) design[[1]]$adapt[[subj]]$design$collapse_across_FM <- TRUE
  ## END

  if (is.null(data$subjects)) {
    jointModel <- TRUE
    all_samples <- emc
  }
  else {
    jointModel <- FALSE
    data <- list(data)
  }
  post_out <- vector("list", length = length(data))
  for (j in 1:length(data)) {
    if (jointModel)
      emc <- single_out_joint(all_samples, j)
    subjects <- levels(data[[j]]$subjects)
    if (hyper) {
      pars <- vector(mode = "list", length = n_post)
    }
    else {
      dots$selection <- "alpha"
      dots$merge_chains <- TRUE
      dots$by_subject <- TRUE
      samps <- do.call(get_pars, c(list(emc), EMC2:::fix_dots(dots,
                                                       get_pars)))
      if (stat != "random") {
        p <- do.call(rbind, lapply(samps, function(x) apply(x[[1]],
                                                            2, stat)))
      }
      pars <- vector(mode = "list", length = n_post)
      for (i in 1:n_post) {
        if (stat != "random")
          pars[[i]] <- p
        else {
          pars[[i]] <- do.call(rbind, lapply(samps, function(x) {
            x[[1]][sample(1:nrow(x[[1]]), 1), ]
          }))
        }
      }
    }
    simDat <- parallel::mclapply(1:n_post, function(i) {
      do.call(EMC2:::make_data, c(list(pars[[i]], design = design[[j]],
                                data = data[[j]]), EMC2:::fix_dots(dots, EMC2:::make_data)))
    }, mc.cores = n_cores)
    out <- cbind(postn = rep(1:n_post, times = unlist(lapply(simDat,
                                                             function(x) dim(x)[1]))), do.call(rbind, simDat))
    ## SM: find updated parameters
    if(!is.null(attr(simDat[[1]], 'npars'))) {
      parsList <- lapply(simDat, attr, 'npars')
      parsList <- lapply(1:length(parsList), function(x) cbind(postn=x, parsList[[x]]))
      outPars <- do.call(rbind, parsList)
      attr(out, 'npars') <- outPars
    }

    if (n_post == 1)
      pars <- pars[[1]]
    attr(out, "pars") <- pars
    post_out[[j]] <- out
  }
  if (!jointModel)
    post_out <- post_out[[1]]
  return(post_out)
}



make_data <- function (parameters, design = NULL, n_trials = NULL, data = NULL,
          expand = 1, mapped_p = FALSE, hyper = FALSE, ...)
{
  LT <- 0
  UT <- Inf
  LC <- 0
  UC <- Inf
  LCresponse <- TRUE
  UCresponse <- TRUE
  LCdirection <- TRUE
  UCdirection <- TRUE
  force_direction <- FALSE
  force_response <- FALSE
  rtContaminantNA <- FALSE
  return_Ffunctions <- FALSE
  Fcovariates = NULL
  optionals <- list(...)
  for (name in names(optionals)) {
    assign(name, optionals[[name]])
  }
  if (is(parameters, "emc")) {
    if (is.null(design))
      design <- get_design(parameters)
    if (is.null(data))
      data <- get_data(parameters)
    if (!hyper) {
      parameters <- do.call(rbind, credint(parameters,
                                           probs = 0.5, selection = "alpha", by_subject = TRUE))
    }
    else {
      mu <- get_pars(parameters, selection = "mu", merge_chains = T,
                     return_mcmc = F)
      Sigma <- get_pars(parameters, selection = "Sigma",
                        merge_chains = T, return_mcmc = F)
      mu <- rowMeans(mu)
      Sigma <- apply(Sigma, 1:2, mean)
      parameters <- make_random_effects(design, group_means = mu,
                                        covariances = Sigma)
    }
  }
  sampled_p_names <- names(sampled_pars(design))
  if (is.null(dim(parameters))) {
    if (is.null(names(parameters)))
      names(parameters) <- sampled_p_names
  }
  else {
    if (length(rownames(parameters)) != length(design$Ffactors$subjects)) {
      stop("input parameter matrix must have number of rows equal to number of subjects specified in design")
    }
    if (is.null(colnames(parameters)))
      colnames(parameters) <- sampled_p_names
    rownames(parameters) <- design$Ffactors$subjects
  }
  if (!is.null(attr(design, "custom_ll"))) {
    data <- list()
    for (i in 1:nrow(parameters)) {
      data[[i]] <- design$model()$rfun(parameters[i, ],
                                       n_trials = n_trials, subject = i)
    }
    return(do.call(rbind, data))
  }
  model <- design$model
  if (is.data.frame(parameters))
    parameters <- as.matrix(parameters)
  if (!is.matrix(parameters))
    parameters <- make_pmat(parameters, design)
  if (is.null(data)) {
    design$Ffactors$subjects <- rownames(parameters)
    if (mapped_p)
      n_trials <- 1
    if (is.null(n_trials))
      stop("If data is not provided need to specify number of trials")
    Ffactors = c(design$Ffactors, list(trials = 1:n_trials))
    data <- as.data.frame.table(array(dim = unlist(lapply(Ffactors,
                                                          length)), dimnames = Ffactors))
    for (i in names(design$Ffactors)) data[[i]] <- factor(data[[i]],
                                                          levels = design$Ffactors[[i]])
    names(data)[dim(data)[2]] <- "R"
    data$R <- factor(data$R, levels = design$Rlevels)
    data$trials <- as.numeric(as.character(data$trials))
    if (!is.null(design$Fcovariates)) {
      if (!is.null(Fcovariates) & !all(unlist(lapply(Fcovariates,
                                                     is.null)))) {
        if (!(all(names(Fcovariates) %in% names(design$Fcovariates))))
          stop("All Fcovariates must be named in design$Fcovariates")
        if (!is.data.frame(Fcovariates)) {
          if (!all(unlist(lapply(Fcovariates, is.function))))
            stop("Fcovariates must be either a data frame or list of functions")
          nams <- names(Fcovariates)
          Fcovariates <- do.call(cbind.data.frame, lapply(Fcovariates,
                                                          function(x) {
                                                            x(data)
                                                          }))
          names(Fcovariates) <- nams
        }
        n <- dim(Fcovariates)[1]
        if (n != nrow(data))
          Fcovariates <- Fcovariates[sample(1:n, nrow(data),
                                            replace = TRUE), , drop = F]
        data <- cbind.data.frame(data, Fcovariates)
      }
      empty_covariates <- names(design$Fcovariates)[!(names(design$Fcovariates) %in%
                                                        names(data))]
      if (length(empty_covariates) > 0)
        data[, empty_covariates] <- 0
    }
  }
  else {
    LT <- attr(data, "LT")
    if (is.null(LT))
      LT <- 0
    UT <- attr(data, "UT")
    if (is.null(UT))
      UT <- Inf
    LC <- attr(data, "LC")
    if (is.null(LC))
      LC <- 0
    UC <- attr(data, "UC")
    if (is.null(UC))
      UC <- Inf
    if (!force_direction) {
      ok <- data$rt == -Inf
      ok[is.na(ok)] <- FALSE
      LCdirection <- any(ok)
      ok <- data$rt == Inf
      ok[is.na(ok)] <- FALSE
      UCdirection = any(ok)
    }
    if (!force_response) {
      if (!any(is.infinite(data$rt)) & any(is.na(data$R))) {
        LCresponse <- UCresponse <- FALSE
      }
      else {
        ok <- data$rt == -Inf
        bad <- is.na(ok)
        LCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
        ok <- data$rt == Inf
        bad <- is.na(ok)
        UCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
      }
    }
    data <- EMC2:::add_trials(data[order(data$subjects), ])
  }
  if (!is.factor(data$subjects))
    data$subjects <- factor(data$subjects)
  if (!is.null(model)) {
    if (!is.function(model))
      stop("model argument must  be a function")
    if (is.null(model()$p_types))
      stop("model()$p_types must be specified")
    if (is.null(model()$Ttransform))
      stop("model()$Ttransform must be specified")
  }

  ## SM
  ## if we want to simulate adaptive models conditional on the empirical RTs, we need the RTs to be present in the 'outcomes' variable
  simulate=TRUE
  if(!is.null(design$adapt)) {
    adapt <- design$adapt[[1]]$design
    if(!is.null(adapt$simulateByTrial)) {
      if(adapt$simulateByTrial==FALSE) {
        simulate=FALSE
      }
    }
  } else {
    simulate=TRUE
  }

  ## OLD (note that simulate is a variable now)
  data <- EMC2:::design_model(EMC2:::add_accumulators(data, design$matchfun,
                                                simulate = simulate, type = model()$type, Fcovariates = design$Fcovariates),
                               design, model, add_acc = FALSE, compress = FALSE, verbose = FALSE,
                               rt_check = FALSE)
  ## NEW
  ## and now we can remove RTs again
  if(!simulate) data$rt <- NA
  if(!is.null(design$adapt)) {
    # if there's an adapt component, there's two ways of simulating: either trial by trial or in one go.
    # trial by trial is preferred in RL designs because it simulates feedback on every trial;
    # in one go (simply uses npars) is faster

    ## Add adapt component to new data
    if (is.null(attr(data,"adapt"))) {
      attr(data,"adapt") <- setNames(lapply(levels(data$subjects), function(x) design$adapt[[x]][[x]]),
                                     levels(data$subjects))
      #setNames(
      #lapply(levels(data$subjects),augment,da=data,design=design),
      #levels(data$subjects))
      attr(data,"adapt")$design <- design$adapt[[1]]$design
    }

    if(!is.null(adapt$simulateByTrial)) {
      if(adapt$simulateByTrial) {
        ## assumes this is FALSE. No Ttransform here!
        pars <- t(apply(parameters, 1, EMC2:::do_pre_transform, model()$pre_transform))
        pars <- EMC2:::map_p(EMC2:::add_constants(pars, design$constants), data,
                             model())
        pars <- EMC2:::do_transform(pars, model()$transform)

        if (expand>1) {
          expand <- 1
          warning("Expand does not work with this type of model")
        }

        ##add response
        data$R <- NA
        daList <- setNames(
          lapply(levels(data$subjects),update_pars,npars=pars,da=data,#return_all=FALSE,
                 rfun=model()$rfun,#return_learning=FALSE,
                 mapped_p=mapped_p),
          levels(data$subjects))

        # Extract parameters, combine with data
        npars <- do.call(rbind, lapply(daList, attr, 'npars'))
        npars <- cbind(data, npars)

        for (i in names(daList)) adapt[[i]] <- attr(daList[[i]],"adapt")[[i]]
        data <- do.call(rbind,daList)
        attr(data,"adapt") <- adapt

        #data <- adapt_data(data,design,model,pars,mapped_p=mapped_p,add_response = TRUE)
        if (mapped_p) return(data)
        data <- data[data$lR==levels(data$lR)[1],!(names(data) %in% c("lR","lM"))]

        if('Qvalues' %in% names(attributes(pars))) attr(data, 'Qvalues') <- attr(pars, 'Qvalues')
        if('predictionErrors' %in% names(attributes(pars))) attr(data, 'predictionErrors') <- attr(pars, 'predictionErrors')
        attr(data,"adapt") <- adapt
        attr(data, 'npars') <- npars
        return(data)
      }
    }
  }
  ## END NEW

  pars <- t(apply(parameters, 1, EMC2:::do_pre_transform, model()$pre_transform))
  pars <- EMC2:::map_p(EMC2:::add_constants(pars, design$constants), data,
                model())
  pars <- EMC2:::do_transform(pars, model()$transform)
  ## SM ADDED A COPY OF NPARS TO SAVE
  pars <- npars <- model()$Ttransform(pars, data)
  npars <- cbind(data, npars)


  pars <- EMC2:::add_bound(pars, model()$bound)
  if (any(dimnames(pars)[[2]] == "pContaminant") && any(pars[,
                                                             "pContaminant"] > 0))
    pc <- pars[data$lR == levels(data$lR)[1], "pContaminant"]
  else pc <- NULL
  if (mapped_p)
    return(cbind(data[, !(names(data) %in% c("R", "rt"))],
                 pars))
  if (expand > 1) {
    data <- cbind(rep = rep(1:expand, each = dim(data)[1]),
                  data.frame(lapply(data, rep, times = expand)))
    lR <- rep(data$lR, expand)
    pars <- apply(pars, 2, rep, times = expand)
  }
  else lR <- data$lR
  if (any(names(data) == "RACE")) {
    Rrt <- matrix(ncol = 2, nrow = dim(data)[1]/length(levels(data$lR)),
                  dimnames = list(NULL, c("R", "rt")))
    RACE <- data[data$lR == levels(data$lR)[1], "RACE"]
    ok <- as.numeric(data$lR) <= as.numeric(as.character(data$RACE))
    for (i in levels(RACE)) {
      pick <- data$RACE == i
      lRi <- factor(data$lR[pick & ok])
      tmp <- pars[pick & ok, ]
      attr(tmp, "ok") <- rep(T, nrow(tmp))
      Rrti <- model()$rfun(lRi, tmp)
      Rrti$R <- as.numeric(Rrti$R)
      Rrt[RACE == i, ] <- as.matrix(Rrti)
    }
    Rrt <- data.frame(Rrt)
    Rrt$R <- factor(Rrt$R, labels = levels(lR), levels = 1:length(levels(lR)))
  }
  else Rrt <- model()$rfun(lR, pars)
  dropNames <- c("lR", "lM", "lSmagnitude")
  if (!return_Ffunctions && !is.null(design$Ffunctions))
    dropNames <- c(dropNames, names(design$Ffunctions))
  data <- data[data$lR == levels(data$lR)[1], !(names(data) %in%
                                                  dropNames)]
  for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[, i]
  data <- EMC2:::make_missing(data[, names(data) != "winner"], LT,
                       UT, LC, UC, LCresponse, UCresponse, LCdirection, UCdirection)
  if (!is.null(pc)) {
    if (!any(is.infinite(data$rt)) & any(is.na(data$R)))
      stop("Cannot have contamination and censoring with no direction and response")
    contam <- runif(length(pc)) < pc
    data[contam, "R"] <- NA
    if (LC != 0 | is.finite(UC)) {
      if ((LCdirection & UCdirection) & !rtContaminantNA)
        stop("Cannot have contamination with a mixture of censor directions")
      if (rtContaminantNA & ((is.finite(LC) & !LCresponse &
                              !LCdirection) | (is.finite(UC) & !UCresponse &
                                               !UCdirection)))
        stop("Cannot have contamination and censoring with no direction and response")
      if (rtContaminantNA | (!LCdirection & !UCdirection))
        data[contam, "rt"] <- NA
      else if (LCdirection)
        data[contam, "rt"] <- -Inf
      else data[contam, "rt"] <- Inf
    }
    else data[contam, "rt"] <- NA
  }
  ## SM ADDED
  attr(data, 'npars') <- npars
  ## END

  attr(data, "p_vector") <- parameters
  data
}

assignInNamespace("make_data", make_data, ns="EMC2")
assignInNamespace("predict.emc", predict.emc, ns="EMC2")
# debug(EMC2:::make_data)
# predict(emc, stage='adapt')
