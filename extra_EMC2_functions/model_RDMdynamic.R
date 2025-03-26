trendTransform = function(npars, da, i) {
  index <- attr(da, 'adapt')[[i]]$index
  if(is.null(index)) return(npars)

  utheta <- npars[1,grepl('utheta', colnames(npars))]
  vtheta <- npars[1,grepl('vtheta', colnames(npars))]
  Btheta <- npars[1,grepl('Btheta', colnames(npars))]

  if(any(utheta != 0)) {
    trend <- utheta[1]*exp(utheta[2]*-(1:(nrow(da)/2)))
    ## urgency: add values to rates
    npars[,'v'] <- npars[,'v']+trend[index]
  }
  if(any(vtheta != 0)) {
    trend <- vtheta[1]*exp(vtheta[2]*-(1:(nrow(da)/2)))
    ## quality: add/subtract values to rates depending on match
    lM <- as.logical(as.character(da[,'lM']))
    direction <- ifelse(lM, 1, -1)
    npars[,'v'] <- npars[,'v']+trend[index]*direction
  }
  if(any(Btheta != 0)) {
    trend <- Btheta[1]*exp(Btheta[2]*-(1:(nrow(da)/2)))
    ## threshold: add values to threshold
    npars[,'B'] <- npars[,'B']+trend[index]
  }
  return(npars)
}

DCTtransform =function(npars, da, i) {
  dct <- attr(da, 'adapt')[[i]]$dct
  index <- attr(da, 'adapt')[[i]]$index
  if(is.null(index)) return(npars)

  ucos <- npars[1,grepl('ucos', colnames(npars))]
  vcos <- npars[1,grepl('vcos', colnames(npars))]
  Bcos <- npars[1,grepl('Bcos', colnames(npars))]

  if(all(c(vcos, Bcos, ucos) == 0)) {
    return(npars)
  }

  if(any(ucos != 0)) {
    ## drifts in trial order
    # ensure correct sizes
    nucos <- length(ucos)
    ndct <- ncol(dct)
    if(nucos < ndct) dct <- dct[,1:nucos]     # less cosines estimated than provided
    if(nucos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nucos-ndct))  # more cosines estimated than provided
    if(nucos == 1) dct <- as.matrix(dct)

    udrift <- dct %*% ucos
    npars[,'v'] <- npars[,'v'] + udrift[index]
  }
  if(any(vcos != 0)) {
    ## drifts in trial order
    # ensure correct sizes
    nvcos <- length(vcos)
    ndct <- ncol(dct)
    if(nvcos < ndct) dct <- dct[,1:nvcos]     # less cosines estimated than provided
    if(nvcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nvcos-ndct))  # more cosines estimated than provided
    if(nvcos == 1) dct <- as.matrix(dct)

    vdrift <- dct %*% vcos
    lM <- as.logical(as.character(da[,'lM']))
    direction <- ifelse(lM, 1, -1)
    npars[,'v'] <- npars[,'v'] + vdrift[index]*direction    # direction = [1, -1] depending on match
  }
  if(any(Bcos != 0)) {
    nBcos <- length(Bcos)
    ndct <- ncol(dct)
    if(nBcos < ndct) dct <- dct[,1:nBcos]     # less cosines estimated than provided
    if(nBcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nBcos-ndct))  # more cosines estimated than provided
    if(nBcos == 1) dct <- as.matrix(dct)

    Bdrift <- dct %*% Bcos
    npars[,'B'] <- npars[,'B'] + Bdrift[index]
  }

  return(npars)
}

polytransform =function(npars, da, i) {
  polyX <- attr(da, 'adapt')[[i]]$poly
  index <- attr(da, 'adapt')[[i]]$index
  if(is.null(index)) return(npars)

  upoly <- npars[1,grepl('upoly', colnames(npars))]
  vpoly <- npars[1,grepl('vpoly', colnames(npars))]
  Bpoly <- npars[1,grepl('Bpoly', colnames(npars))]

  if(all(c(upoly, vpoly, Bpoly) == 0)) {
    return(npars)
  }

  if(any(upoly != 0)) {
    ## drifts in trial order
    # ensure correct sizes
    nucos <- length(upoly)
    ndct <- ncol(polyX)
    if(nucos < ndct) polyX <- polyX[,1:nucos]     # less cosines estimated than provided
    if(nucos > ndct) polyX <- cbind(polyX, matrix(0, nrow=nrow(polyX), ncol=nucos-ndct))  # more cosines estimated than provided

    if(nucos == 1) polyX <- as.matrix(polyX)
    udrift <- polyX %*% upoly
    npars[,'v'] <- npars[,'v'] + udrift[index]
  }
  if(any(vpoly != 0)) {
    ## drifts in trial order
    # ensure correct sizes
    nvcos <- length(vpoly)
    ndct <- ncol(polyX)
    if(nvcos < ndct) polyX <- polyX[,1:nvcos]     # less cosines estimated than provided
    if(nvcos > ndct) polyX <- cbind(polyX, matrix(0, nrow=nrow(polyX), ncol=nvcos-ndct))  # more cosines estimated than provided

    if(nvcos == 1) polyX <- as.matrix(polyX)
    vdrift <- polyX %*% vpoly
    lM <- as.logical(as.character(da[,'lM']))
    direction <- ifelse(lM, 1, -1)
    npars[,'v'] <- npars[,'v'] + vdrift[index]*direction    # direction = [1, -1] depending on match
  }
  if(any(Bpoly != 0)) {
    nBcos <- length(Bpoly)
    ndct <- ncol(polyX)
    if(nBcos < ndct) polyX <- polyX[,1:nBcos]     # less cosines estimated than provided
    if(nBcos > ndct) polyX <- cbind(polyX, matrix(0, nrow=nrow(polyX), ncol=nBcos-ndct))  # more cosines estimated than provided
    if(nBcos == 1) polyX <- as.matrix(polyX)

    Bdrift <- polyX %*% Bpoly
    npars[,'B'] <- npars[,'B'] + Bdrift[index]
  }

  return(npars)
}

rdmDynamic <- function() {
  list(
    type="RACE",
    transform=list(func=c(alpha1='pnorm', alpha2='pnorm', alpha3='pnorm', q01='pnorm', q02='pnorm', Btheta2='pnorm', utheta2='pnorm', vtheta2='pnorm',
                          A='exp', t0='exp', s='exp', q03='exp')),
    bound=list(minmax=cbind(v=c(1e-3,Inf),
                            B=c(0,Inf),
                            A=c(1e-4,Inf),
                            t0=c(0.05,Inf),
                            s=c(0,Inf)),
               exception=c(A=0, v=0)),
    p_types=c("v",
              "B",
              "A",
              "t0",
              "s",
              "alpha1", "alpha2", "alpha3",
              # "alpha1min", "alpha2min", "alpha3min",
              "q01", "q02", "q03",
              "weight1", "weight2", "weight3",
              paste0(rep('ucos', 3), 1:3),     # effect on URGENCY
              paste0(rep('Bcos', 3), 1:3),     # effect on THRESHOLD
              paste0(rep('vcos', 3), 1:3),     # effect on DRIFT
              paste0(rep('vtheta', 2), 1:2),  # exponential trend on drift
              paste0(rep('utheta', 2), 1:2),  # exponential trend on urgency
              paste0(rep('Btheta', 2), 1:2)),   # exponential trend on threshold

    # Trial dependent parameter transform
    Ttransform = function(pars, dadm) {
      parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                           levels(dadm$subjects))
      for (i in levels(dadm$subjects)) {
        if(is.null(attr(dadm, 'adapt'))) {
          parsList[[i]] <- pars[dadm$subjects==i,]
        } else {
          ## transform cosines first
          npars <- DCTtransform(npars=pars[dadm$subjects==i,], da=dadm[dadm$subjects==i,], i)

          ## potentially, polynomial basis functions
          npars <- polytransform(npars=npars, da=dadm[dadm$subjects==i,], i)

          ## potentially, exponential trends
          npars <- trendTransform(npars=npars, da=dadm[dadm$subjects==i,], i)

          ## and to updating next
          pars[dadm$subjects==i,] <- npars
          parsList[[i]] <- update_pars(i,pars,dadm)
        }
      }

      pars <- do.call(rbind,parsList)
      # pars[,'v'] <- pmax(pars[,'v'],0)
      # pars[,'v'] <- pmin(pars[,'v'],30)
      # pars[,'B'] <- pmax(pars[,'B'],0.1)
      # pars[,'B'] <- pmin(pars[,'B'],30)
      # attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) # & (pars[,'v'] > 1e-6) & (pars[,'v'] < 1e3) & (pars[,'B'] > .1) & (pars[,'alpha'] < 1)
      pars
    },
    # p_vector transform
    transform = function(x) x,
    # Random function for racing accumulators
    rfun=function(lR,pars) EMC2:::rRDM(lR,pars),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) EMC2:::dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) EMC2:::pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      EMC2:::log_likelihood_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}


generateRDMDynamic <- function(
    transform=list(func=c(alpha1='pnorm', alpha2='pnorm', alpha3='pnorm', q01='pnorm', q02='pnorm', A='exp', t0='exp', s='exp', q03='exp')),
    bound=list(minmax=cbind(v=c(1e-3,Inf),B=c(0,Inf),A=c(1e-4,Inf),t0=c(0.05,Inf),s=c(0,Inf)),
               exception=c(A=0, v=0)),
    p_types=c("v"=1,
              "B"=1,
               "A"=0,
               "t0"=0.1,
               "s"=1,
               "alpha1"=-Inf,
               "alpha2"=-Inf,
               "alpha3"=-Inf,
               # "alpha1min", "alpha2min", "alpha3min",
               "q01"=qnorm(0.5),
               "q02"=qnorm(1),
               "q03"=log(3),
               "weight1"=0,
               "weight2"=0,
               "weight3"=0,
#               'alpha1min'=0,
#               'alpha2min'=0,
#               'alpha3min'=0,
               setNames(rep(0,3), paste0(rep('ucos', 3), 1:3)),     # effect on URGENCY
               setNames(rep(0,3), paste0(rep('Bcos', 3), 1:3)),     # effect on THRESHOLD
               setNames(rep(0,3), paste0(rep('vcos', 3), 1:3)),     # effect on DRIFT
               setNames(rep(0,2), paste0(rep('vtheta', 2), 1:2)),  # exponential trend on drift
               setNames(rep(0,2), paste0(rep('utheta', 2), 1:2)),  # exponential trend on urgency
               setNames(rep(0,2), paste0(rep('Btheta', 2), 1:2)))) {
    model <- rdmDynamic()
    model$p_types <- p_types
    model$transform <- transform
    model$bound <- bound

    f <- function() return(model)
    environment(f) <- list2env(list(transform=transform, model=model, p_types=p_types, bound=bound, transform=transform), parent = globalenv())   # empty environment
    return(f)
}
