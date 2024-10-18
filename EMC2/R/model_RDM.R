

# n-choice uniformly varying start point (0-A) Wald
#    race, with t0, v, A, b (boundary) parameterixaiton

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with:
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review. Comments and changes added by Andrew Heathcote. Trish's
# code is for k = threshold, a = half width of uniform threshold variability,
# l = rate of accumulation. Note that Wald mean = k/l and shape = k^2.
#
# I use a different parameterization in terms of v=l (rate),
# uniform start point variability from 0-A (A>=0), threshold b (>0) and hence
# B=b-A (>=0) as a threshold gap. Hence k = b-A/2 = B + A/2 and a=A
#
# Added ability to set s = diffusive standard deviation, by scaling A, B and v
# i.e., given s, same result for A/s, B/s, v/s with s=1

#### distribution functions

# # Moved to C++ model_RDM.cpp
#
# digt0 <- function(t,k=1,l=1) {
#       # pdf of inverse gaussian at t with no k variability
#       # much faster than statmod's dinvgauss funciton
#
#       lambda <- k^2
#       l0 <- l==0
#       e <- numeric(length(t))
#       if ( any(!l0) ) {
#         mu <- k[!l0]/l[!l0]
#         e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
#       }
#       if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
#       x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
#       x[t<=0] <- 0
#       x
# }
#
#
# digt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
#     # pdf of inverse gaussian at t with k +/- a/2 uniform variability
#     # returns digt0 if a<1e-10
#
#
#     options(warn=-1)
#     if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
#     if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
#     if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
#
#     tpos <- t<=0
#
#     atiny <- a<=tiny & !tpos
#     a[atiny] <- 0
#
#     ltiny <- (l<=tiny) & !atiny & !tpos
#     notltiny <- (l>tiny) & !atiny & !tpos
#     l[l<=tiny] <- 0
#
#     x <- numeric(length(t))
#
#     # No threshold variability
#     if ( any(atiny) )
#       x[atiny] <- digt0(t=t[atiny],k=k[atiny],l=l[atiny])
#
#     # Threshold variability
#     if ( any(!atiny) ) {
#
#       if ( any(notltiny) ) { # rate non-zero
#
#         sqr.t <- sqrt(t[notltiny])
#
#         term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
#         term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
#         term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])
#
#         term.2a <- log(.5)+log(l[notltiny])
#         term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
#         term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.2d <- term.2b+term.2c
#         term.2 <- exp(term.2a)*term.2d
#
#         term.3 <- term.1+term.2
#         term.4 <- log(term.3)-log(2)-log(a[notltiny])
#         x[notltiny] <- exp(term.4)
#       }
#
#       if ( any(ltiny) ) {  # rate zero
#         log.t <- log(t[ltiny])
#         term.1 <- -.5*(log(2)+log(pi)+log.t)
#         term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
#         term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
#         term.4 <- (exp(-term.2)-exp(-term.3))
#         term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
#         x[ltiny] <- exp(term.5)
#       }
#
#     }
#
#     x[x<0 | is.nan(x) ] <- 0
#     x
# }
#
#
# pigt0 <- function(t,k=1,l=1) {
#       # cdf of inverse gaussian at t with no k variability
#       # much faster than statmod's pinvgauss funciton
#
#       mu <- k/l
#       lambda <- k^2
#
#       e <- exp(log(2*lambda) - log(mu))
#       add <- sqrt(lambda/t) * (1 + t/mu)
#       sub <- sqrt(lambda/t) * (1 - t/mu)
#
#       p.1 <- 1 - pnorm(add)
#       p.2 <- 1 - pnorm(sub)
#       x <- exp(e + log(p.1)) + p.2
#
#       x[t<0] <- 0
#       x
# }
#
#
# pigt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
#     # cdf of inverse gaussian at t with k +/- a/2 uniform variability
#     # returns pigt0 if a<=0
#
#     options(warn=-1)
#     if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
#     if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
#     if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
#
#     tpos <- t<=0
#
#     atiny <- a<=tiny & !tpos
#     a[atiny] <- 0
#
#     ltiny <- (l<=tiny) & !atiny & !tpos
#     notltiny <- (l>tiny) & !atiny & !tpos
#     l[l<=tiny] <- 0
#
#     x <- numeric(length(t))
#
#     # No threshold variability
#     if ( any(atiny) )
#       x[atiny] <- pigt0(t[atiny],k[atiny],l[atiny])
#
#     # Threshold variability
#     if ( any(!atiny) ) {
#
#       if ( any(notltiny) ) { # rate non-zero
#
#         log.t <- log(t[notltiny])
#         sqr.t <- sqrt(t[notltiny])
#
#         term.1a <- .5*log.t-.5*log(2*pi)
#         term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
#         term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
#         term.1 <- exp(term.1a)*(term.1b-term.1c)
#
#         term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) +
#                          log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
#         term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) +
#                          log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
#         term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
#
#         term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
#         term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
#         term.4 <- term.4c*term.4a + term.4d*term.4b
#
#         x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
#       }
#
#       if ( any(ltiny) ) {  # rate zero
#         sqr.t <- sqrt(t[ltiny])
#         log.t <- log(t[ltiny])
#         term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
#         term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
#         term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
#
#         term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
#         term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
#         term.6 <- 1 + exp(term.6b) - exp(term.6a)
#
#         x[ltiny] <- term.5 + term.6
#       }
#
#     }
#
#     x[x<0 | is.nan(x) ] <- 0
#     x
# }
#
#
dRDM <- function(rt,pars)
  # density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[,"t0"] &
    !pars[,"v"] < 0  # code handles rate zero case
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[ok,c("A","B","v")] <- pars[ok,c("A","B","v")]/pars[ok,"s"]
  out[ok] <- dWald(rt[ok],v=pars[ok,"v"],B=pars[ok,"B"],A=pars[ok,"A"],t0=pars[ok,"t0"])
  out
}

pRDM <- function(rt,pars)
  # cumulative density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[,"t0"] &
    !pars[,"v"] < 0  # code handles rate zero case
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[ok,c("A","B","v")] <- pars[ok,c("A","B","v")]/pars[ok,"s"]
  out[ok] <- pWald(rt[ok],v=pars[ok,"v"],B=pars[ok,"B"],A=pars[ok,"A"],t0=pars[ok,"t0"])
  out
}

#### random

rWald <- function(n,B,v,A)
  # random function for single accumulator
{

  rwaldt <- function(n,k,l,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss

    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }

    flag <- l>tiny
    x <- rep(NA,times=n)

    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- k^2

    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]

    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)

    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }

  # Kluge to return Inf for negative rates
  out <- numeric(n)
  ok <- !v<0
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok])
  out[!ok] <- Inf
  out
}


rRDM <- function(lR,pars,p_types=c("v","B","A","t0"))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. "s" parameter will be used but can be omotted
  #
  # test
  # pars=cbind(B=c(1,2),v=c(1,1),A=c(0,0),t0=c(.2,.2)); lR=factor(c(1,2))
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
  pars[,"B"][pars[,"B"]<0] <- 0 # Protection for negatives
  pars[,"A"][pars[,"A"]<0] <- 0
  dt <- matrix(rWald(dim(pars)[1],B=pars[,"B"],v=pars[,"v"],A=pars[,"A"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(pars[,"t0"],nrow=length(levels(lR)))[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}


#' The Racing Diffusion Model (RDM)
#'
#' The Racing Diffusion Model, also known as the Racing Wald Model, proposes that for each choice alternative, noisy accumulators race towards a common bound.
#' The first accumulator to reach the bound determines the choice made. The time taken to reach the threshold determines the response times. For details see `Tillman, Van Zandt & Logan, 2020`
#'
#' The core parameters of the RDM are the drift rate `v`, the response threshold `B`,
#' within trial variation in drift rate `s`, between trial variation in startpoint of the drift rate `A`, and non-decision time `t0`.
#' Frequently `s` is fixed to 1 to satisfy scaling constraints.
#'
#' Here we use the b = B + A parameterization, which ensures that the response threshold is always higher than the between trial variation in start point of the drift rate.
#'
#' @return
#' @export
#'
#' @examples
rdmB <- function(){
  list(
    type="RACE",
    c_name = "rdmB",
    p_types=c("v","B","A","t0","s"),
    # Transform to natural scale
    Ntransform=function(x) {
      # transform parameters back to real line
      exp(x)
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR,pars) rRDM(lR,pars),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

#' The Racing Diffusion Model (RDM)
#'
#' The Racing Diffusion Model, also known as the Racing Wald Model, proposes that for each choice alternative, noisy accumulators race towards a common bound.
#' The first accumulator to reach the bound determines the choice made. The time taken to reach the threshold determines the response times. For details see `Tillman, Van Zandt & Logan, 2020`
#'
#' The core parameters of the RDM are the drift rate `v`, the response threshold `B`,
#' within trial variation in drift rate `s`, between trial variation in startpoint of the drift rate `A`, and non-decision time `t0`.
#' Frequently `s` is fixed to 1 to satisfy scaling constraints.
#'
#' Here we use the b = B + A parameterization, which ensures that the response threshold is always higher than the between trial variation in start point of the drift rate.
#'
#' @return
#' @export
#'
#' @examples
rdmBnoc <- function(){
  list(
    type="RACE",
#    c_name = "rdmB",
    p_types=c("v","B","A","t0","s"),
    # Transform to natural scale
    Ntransform=function(x) {
      # transform parameters back to real line
      exp(x)
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR,pars) rRDM(lR,pars),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}


# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' Title
#'
#' @return
#' @export
#'
#' @examples
rdmBt0natural <- function(){
  list(
    type="RACE",
    p_types=c("v","B","A","t0","s"),
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]  != "t0"] <- exp(x[,dimnames(x)[[2]]  != "t0"])
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR,pars) rRDM(lR,pars),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}


# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Added v0 parameter for initial rates in RL version.
# rdmRL <- list(
#   type="RACE",
#   p_types=c("v0","v","B","A","t0","alpha","w","q0","s"),
#
#   Ntransform=function(x) {
#     # Transform to natural scale
#     exp(x)
#   },
#   # Trial dependent parameter transform
#   Ttransform = function(pars,dadm) {
#     parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
#                          levels(dadm$subjects))
#     for (i in levels(dadm$subjects))
#       parsList[[i]] <- update_pars(i,pars,dadm)
#     pars <- do.call(rbind,parsList)
#     attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
#     pars
#   },
#   # p_vector transform
#   transform = function(x) x,
#   # Random function for racing accumulators
#   rfun=function(lR,pars) rRDM(lR,pars),
#   # Density function (PDF) for single accumulator
#   dfun=function(rt,pars) dRDM(rt,pars),
#   # Probability function (CDF) for single accumulator
#   pfun=function(rt,pars) pRDM(rt,pars),
#   # Race likelihood combining pfun and dfun
#   log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
#     log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
# )


## SM: RL and ARD combined
# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Added q0 parameter for initial rates in RL version.
#' RL-ARD with RDM as model
#'
#' @return
#' @export
#'
#' @examples
rdmRLARD <- function(){list(
  type="RACE",
  p_types=c("v0","v","B","A","t0","alpha","q0","s", "wd", "ws"),

  Ntransform=function(x) {
    # Transform to natural scale
    probit_scale <- c('alpha', 'q0')
    x[,!dimnames(x)[[2]] %in% probit_scale] <- exp(x[,!dimnames(x)[[2]] %in% probit_scale])
    x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    parsList <- Qlist <- PElist <- setNames(vector(mode="list",
                                                   length=length(levels(dadm$subjects))),
                                            levels(dadm$subjects))

    for (i in levels(dadm$subjects)) {
      npars <- update_pars(i,pars,dadm)
      parsList[[i]] <- npars
      if('learn' %in% names(attributes(npars))) {
        Qlist[[i]] <- attr(npars, 'learn')$adaptedValues
        PElist[[i]] <- apply(attr(npars, 'learn')$predictionErrors, 1, sum, na.rm=TRUE)
      }
    }
    pars <- do.call(rbind,parsList)
    if('learn' %in% names(attributes(npars))) {
      attr(pars, 'predictionErrors') <- do.call(rbind, PElist)
      attr(pars, 'Qvalues') <- do.call(rbind, Qlist)
    }
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) & (pars[,'v0']>0.01)
    pars
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rRDM(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dRDM(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pRDM(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)
}


## SM: RL and RD combined
# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' RL-RD with RDM as model
#'
#' @return
#' @export
#'
#' @examples
rdmRL <- function() {list(
  type="RACE",
  p_types=c("v0","v","B","A","t0","alpha","w","q0","s"),

  Ntransform=function(x) {
    # Transform to natural scale
    probit_scale <- c('alpha', 'q0')
    x[,!dimnames(x)[[2]] %in% probit_scale] <- exp(x[,!dimnames(x)[[2]] %in% probit_scale])
    x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                         levels(dadm$subjects))
    for (i in levels(dadm$subjects))
      parsList[[i]] <- update_pars(i,pars,dadm)
    pars <- do.call(rbind,parsList)
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) # & (pars[,'v'] > 1e-6) & (pars[,'v'] < 1e3) & (pars[,'B'] > .1) & (pars[,'alpha'] < 1)
    pars
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rRDM(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dRDM(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pRDM(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)
}

## SM: VKF and RD combined
# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' RL-RD with VKF and RDM as model
#'
#' @return
#' @export
#'
#' @examples
rdmRLvkf <- function() {list(
  type="RACE",
  p_types=c("v0","v","B","A","t0","w","s",
            "alpha", "q0", "w0", "volatility0"),

  Ntransform=function(x) {
    # Transform to natural scale
    probit_scale <- c('alpha', 'q0', 'volatility0', 'w0')
    x[,!dimnames(x)[[2]] %in% probit_scale] <- exp(x[,!dimnames(x)[[2]] %in% probit_scale])
    x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                         levels(dadm$subjects))
    for (i in levels(dadm$subjects))
      parsList[[i]] <- update_pars(i,pars,dadm)
    pars <- do.call(rbind,parsList)
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
    pars
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rRDM(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dRDM(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pRDM(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)
}

## SM: VKF and ARD combined
# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' RL-ARD with VKF and RDM as model
#'
#' @return
#' @export
#'
#' @examples
rdmRLARDvkf <- function(){list(
  type="RACE",
  p_types=c("v0","v","B","A","t0","s", "wd", "ws",
            "alpha", "q0", "w0", "volatility0"),

  Ntransform=function(x) {
    # Transform to natural scale
    probit_scale <- c('alpha', 'q0', 'volatility0', 'w0')
    # actually, q0 can take more flexible norms outside of the binomial reward function?
    x[,!dimnames(x)[[2]] %in% probit_scale] <- exp(x[,!dimnames(x)[[2]] %in% probit_scale])
    x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                         levels(dadm$subjects))
    for (i in levels(dadm$subjects))
      parsList[[i]] <- update_pars(i,pars,dadm)
    pars <- do.call(rbind,parsList)
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
    pars
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rRDM(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dRDM(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pRDM(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)
}

#'
#' ### SM: RDM dynamic with DCT for slow waves
#' # RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' # Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' #' RL-RD with RDM as model
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' rdmDynamicDCT3 <- function() {list(
#'   type="RACE",
#'   p_types=c("v","B","A","t0","s",
#'             "alpha1", "alpha2", "alpha3",
#'             "alpha1min", "alpha2min", "alpha3min",
#'             "q01", "q02", "q03",
#'             "weight1", "weight2", "weight3",
#'             paste0(rep('vcos', 10), 1:10),    # effect on DRIFT
#'             paste0(rep('Bcos', 10), 1:10)     # effect on THRESHOLD
#'   ),  # maximally 10 cosines, should really be enough
#'
#'   Ntransform=function(x) {
#'     # Transform to natural scale
#'     probit_scale <- c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     logit_scale = c("v", "B", "A", "t0", "s")
#'
#'     pnames <- dimnames(x)[[2]]
#'     ## weights on normal scale?
#'     x[,pnames %in% logit_scale] <- exp(x[,pnames %in% logit_scale])
#'     if(any(pnames %in% probit_scale)) {
#'       x[,pnames %in% probit_scale] <- pnorm(x[,pnames %in% probit_scale])
#'     }
#'
#'     # scale alphas to set range [alphamin, 1]
#'     if('alpha1min' %in% pnames) x[,pnames == 'alpha1'] <- x[,pnames == 'alpha1']*(1-x[,pnames == 'alpha1min']) + x[,pnames == 'alpha1min']
#'     if('alpha2min' %in% pnames) x[,pnames == 'alpha2'] <- x[,pnames == 'alpha2']*(1-x[,pnames == 'alpha2min']) + x[,pnames == 'alpha2min']
#'     if('alpha3min' %in% pnames) x[,pnames == 'alpha3'] <- x[,pnames == 'alpha3']*(1-x[,pnames == 'alpha3min']) + x[,pnames == 'alpha3min']
#'
#'     x
#'
#'   },
#'   DCTtransform =function(npars, da, i) {
#'     dct <- attr(da, 'adapt')[[i]]$dct
#'     index <- attr(da, 'adapt')[[i]]$index
#'
#'     vcos <- npars[1,grepl('vcos', colnames(npars))]
#'     Bcos <- npars[1,grepl('Bcos', colnames(npars))]
#'
#'     if(all(c(vcos, Bcos) == 0)) {
#'       return(npars)
#'     }
#'
#'     if(any(vcos != 0)) {
#'       ## drifts in trial order
#'       # ensure correct sizes
#'       nvcos <- length(vcos)
#'       ndct <- ncol(dct)
#'       if(nvcos < ndct) dct <- dct[,1:nvcos]     # less cosines estimated than provided
#'       if(nvcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nvcos-ndct))  # more cosines estimated than provided
#'
#'       vdrift <- dct %*% vcos
#'       npars[,'v'] <- npars[,'v'] + vdrift[index]  # to "EMC-order"
#'     }
#'     if(any(Bcos != 0)) {
#'       nBcos <- length(Bcos)
#'       ndct <- ncol(dct)
#'       if(nBcos < ndct) dct <- dct[,1:nBcos]     # less cosines estimated than provided
#'       if(nBcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nBcos-ndct))  # more cosines estimated than provided
#'
#'       Bdrift <- dct %*% Bcos
#'       npars[,'B'] <- npars[,'B'] + Bdrift[index]
#'     }
#'
#'     return(npars)
#'   },
#'   # Trial dependent parameter transform
#'   Ttransform = function(pars,dadm) {
#'     parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
#'                          levels(dadm$subjects))
#'     for (i in levels(dadm$subjects)) {
#'       ## transform cosines first
#'       npars <- rdmDynamicDCT3()$DCTtransform(npars=pars[dadm$subjects==i,], da=dadm[dadm$subjects==i,], i)
#'
#'       ## and to updating next
#'       pars[dadm$subjects==i,] <- npars
#'       parsList[[i]] <- update_pars(i,pars,dadm)
#'     }
#'
#'     pars <- do.call(rbind,parsList)
#'     pars[,'v'] <- pmax(pars[,'v'],0)
#'     pars[,'v'] <- pmin(pars[,'v'],30)
#'     pars[,'B'] <- pmax(pars[,'B'],0.1)
#'     pars[,'B'] <- pmin(pars[,'B'],30)
#'     attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) # & (pars[,'v'] > 1e-6) & (pars[,'v'] < 1e3) & (pars[,'B'] > .1) & (pars[,'alpha'] < 1)
#'     pars
#'   },
#'   # p_vector transform
#'   transform = function(x) x,
#'   # Random function for racing accumulators
#'   rfun=function(lR,pars) rRDM(lR,pars),
#'   # Density function (PDF) for single accumulator
#'   dfun=function(rt,pars) dRDM(rt,pars),
#'   # Probability function (CDF) for single accumulator
#'   pfun=function(rt,pars) pRDM(rt,pars),
#'   # Race likelihood combining pfun and dfun
#'   log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
#'     log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
#' )
#' }
#'
#' ### SM: RDM dynamic with DCT for slow waves
#' # RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' # Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' #' RL-RD with RDM as model
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' rdmDynamicDCT4 <- function() {list(
#'   type="RACE",
#'   p_types=c("v","B","A","t0","s",
#'             "alpha1", "alpha2", "alpha3",
#'             "alpha1min", "alpha2min", "alpha3min",
#'             "q01", "q02", "q03",
#'             "weight1", "weight2", "weight3",
#'             paste0(rep('vcos', 10), 1:10),    # effect on DRIFT
#'             paste0(rep('Bcos', 10), 1:10)     # effect on THRESHOLD
#'   ),  # maximally 10 cosines, should really be enough
#'
#'   Ntransform=function(x) {
#'     # Transform to natural scale
#'     probit_scale <- c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     logit_scale = c("v", "B", "A", "t0", "s")
#'
#'     ## set urgency signal v to minimally exp(-3) [~0.05] - the likelihood landscape is extremely flat below this value, causing sampling issues
#'     x[x[,'v'] < -3,'v'] <- 10  # set to extremely high value so likelihood always becomes 0
#'
#'     pnames <- dimnames(x)[[2]]
#'     ## weights on normal scale
#'     x[,pnames %in% logit_scale] <- exp(x[,pnames %in% logit_scale])
#'     if(any(pnames %in% probit_scale)) {
#'       x[,pnames %in% probit_scale] <- pnorm(x[,pnames %in% probit_scale])
#'     }
#'
#'     # scale alphas to set range [alphamin, 1]
#'     if('alpha1min' %in% pnames) x[,pnames == 'alpha1'] <- x[,pnames == 'alpha1']*(1-x[,pnames == 'alpha1min']) + x[,pnames == 'alpha1min']
#'     if('alpha2min' %in% pnames) x[,pnames == 'alpha2'] <- x[,pnames == 'alpha2']*(1-x[,pnames == 'alpha2min']) + x[,pnames == 'alpha2min']
#'     if('alpha3min' %in% pnames) x[,pnames == 'alpha3'] <- x[,pnames == 'alpha3']*(1-x[,pnames == 'alpha3min']) + x[,pnames == 'alpha3min']
#'
#'     x
#'
#'   },
#'   DCTtransform =function(npars, da, i) {
#'     dct <- attr(da, 'adapt')[[i]]$dct
#'     index <- attr(da, 'adapt')[[i]]$index
#'
#'     vcos <- npars[1,grepl('vcos', colnames(npars))]
#'     Bcos <- npars[1,grepl('Bcos', colnames(npars))]
#'
#'     if(all(c(vcos, Bcos) == 0)) {
#'       return(npars)
#'     }
#'
#'     if(any(vcos != 0)) {
#'       ## drifts in trial order
#'       # ensure correct sizes
#'       nvcos <- length(vcos)
#'       ndct <- ncol(dct)
#'       if(nvcos < ndct) dct <- dct[,1:nvcos]     # less cosines estimated than provided
#'       if(nvcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nvcos-ndct))  # more cosines estimated than provided
#'
#'       vdrift <- dct %*% vcos
#'       lM <- as.logical(as.character(da[,'lM']))
#'       direction <- ifelse(lM, 1, -1)
#'       npars[,'v'] <- npars[,'v'] + vdrift[index]*direction    # direction = [1, -1] depending on match
#'     }
#'     if(any(Bcos != 0)) {
#'       nBcos <- length(Bcos)
#'       ndct <- ncol(dct)
#'       if(nBcos < ndct) dct <- dct[,1:nBcos]     # less cosines estimated than provided
#'       if(nBcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nBcos-ndct))  # more cosines estimated than provided
#'
#'       Bdrift <- dct %*% Bcos
#'       npars[,'B'] <- npars[,'B'] + Bdrift[index]
#'     }
#'
#'     return(npars)
#'   },
#'   # Trial dependent parameter transform
#'   Ttransform = function(pars,dadm) {
#'     parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
#'                          levels(dadm$subjects))
#'     for (i in levels(dadm$subjects)) {
#'       ## transform cosines first
#'       npars <- rdmDynamicDCT4()$DCTtransform(npars=pars[dadm$subjects==i,], da=dadm[dadm$subjects==i,], i)
#'
#'       ## and to updating next
#'       pars[dadm$subjects==i,] <- npars
#'       parsList[[i]] <- update_pars(i,pars,dadm)
#'     }
#'
#'     pars <- do.call(rbind,parsList)
#'     pars[,'v'] <- pmax(pars[,'v'],0)
#'     pars[,'v'] <- pmin(pars[,'v'],30)
#'     pars[,'B'] <- pmax(pars[,'B'],0.1)
#'     pars[,'B'] <- pmin(pars[,'B'],30)
#'     attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) # & (pars[,'v'] > 1e-6) & (pars[,'v'] < 1e3) & (pars[,'B'] > .1) & (pars[,'alpha'] < 1)
#'     pars
#'   },
#'   # p_vector transform
#'   transform = function(x) x,
#'   # Random function for racing accumulators
#'   rfun=function(lR,pars) rRDM(lR,pars),
#'   # Density function (PDF) for single accumulator
#'   dfun=function(rt,pars) dRDM(rt,pars),
#'   # Probability function (CDF) for single accumulator
#'   pfun=function(rt,pars) pRDM(rt,pars),
#'   # Race likelihood combining pfun and dfun
#'   log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
#'     log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
#' )
#' }
#'
#' ### SM: Generate RDM dynamic with DCT for slow waves
#' # RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' # Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' #' RL-RD with RDM as model
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' generateRDMdynamic <- function(p_types, constantsFixed=NULL,
#'                                probit_scale=c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03"),
#'                                logit_scale=c("v", "B", "A", "t0", "s")) {
#'   model <- rdmDynamicDCT3()
#'   model$p_types <- p_types
#'   model <- model[names(model)!='Ntransform']
#'   model$Ntransform=function(x,
#'                             constants=constantsFixed,
#'                             probit_scale_ = probit_scale,
#'                             logit_scale_ = logit_scale) {
#'     if(!is.null(constants)) {
#'       for(nm in names(constants)) {
#'         x <- cbind(x, nm=constants[[nm]])
#'         colnames(x)[ncol(x)] <- nm
#'       }
#'     }
#'     # Transform to natural scale
#'     # probit_scale=c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     # logit_scale = c("v", "B", "A", "t0", "s")
#'
#'     pnames <- dimnames(x)[[2]]
#'     ## weights on normal scale?
#'     x[,pnames %in% logit_scale_] <- exp(x[,pnames %in% logit_scale_])
#'     if(any(pnames %in% probit_scale_)) {
#'       x[,pnames %in% probit_scale_] <- pnorm(x[,pnames %in% probit_scale_])
#'     }
#'
#'     # scale alphas to set range [alphamin, 1]
#'     if('alpha1min' %in% pnames) x[,pnames == 'alpha1'] <- x[,pnames == 'alpha1']*(1-x[,pnames == 'alpha1min']) + x[,pnames == 'alpha1min']
#'     if('alpha2min' %in% pnames) x[,pnames == 'alpha2'] <- x[,pnames == 'alpha2']*(1-x[,pnames == 'alpha2min']) + x[,pnames == 'alpha2min']
#'     if('alpha3min' %in% pnames) x[,pnames == 'alpha3'] <- x[,pnames == 'alpha3']*(1-x[,pnames == 'alpha3min']) + x[,pnames == 'alpha3min']
#'
#'     #
#'     if('alpha1Posmin' %in% pnames) x[,pnames == 'alpha1Pos'] <- x[,pnames == 'alpha1Pos']*(1-x[,pnames == 'alpha1Posmin']) + x[,pnames == 'alpha1Posmin']
#'     if('alpha2Posmin' %in% pnames) x[,pnames == 'alpha2Pos'] <- x[,pnames == 'alpha2Pos']*(1-x[,pnames == 'alpha2Posmin']) + x[,pnames == 'alpha2Posmin']
#'     if('alpha3Posmin' %in% pnames) x[,pnames == 'alpha3Pos'] <- x[,pnames == 'alpha3Pos']*(1-x[,pnames == 'alpha3Posmin']) + x[,pnames == 'alpha3Posmin']
#'     if('alpha1Negmin' %in% pnames) x[,pnames == 'alpha1Neg'] <- x[,pnames == 'alpha1Neg']*(1-x[,pnames == 'alpha1Negmin']) + x[,pnames == 'alpha1Negmin']
#'     if('alpha2Negmin' %in% pnames) x[,pnames == 'alpha2Neg'] <- x[,pnames == 'alpha2Neg']*(1-x[,pnames == 'alpha2Negmin']) + x[,pnames == 'alpha2Negmin']
#'     if('alpha3Negmin' %in% pnames) x[,pnames == 'alpha3Neg'] <- x[,pnames == 'alpha3Neg']*(1-x[,pnames == 'alpha3Negmin']) + x[,pnames == 'alpha3Negmin']
#'
#'     x
#'   }
#'
#'   f <- function() return(model)
#'   environment(f) <- list2env(list(constantsFixed=constantsFixed, model=model,probit_scale=probit_scale, logit_scale=logit_scale), parent = globalenv())   # empty environment
#'   return(f)
#' }
#'
#' ### SM: Generate RDM dynamic with DCT for slow waves
#' # RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' # Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' #' RL-RD with RDM as model
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' generateRDMdynamic2 <- function(p_types, constantsFixed=NULL,
#'                                probit_scale=c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03"),
#'                                logit_scale=c("v", "B", "A", "t0", "s")) {
#'   model <- rdmDynamicDCT4()
#'   model$p_types <- p_types
#'   model <- model[names(model)!='Ntransform']
#'   model$Ntransform=function(x,
#'                             constants=constantsFixed,
#'                             probit_scale_ = probit_scale,
#'                             logit_scale_ = logit_scale) {
#'     if(!is.null(constants)) {
#'       for(nm in names(constants)) {
#'         x <- cbind(x, nm=constants[[nm]])
#'         colnames(x)[ncol(x)] <- nm
#'       }
#'     }
#'     # Transform to natural scale
#'     # probit_scale=c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     # logit_scale = c("v", "B", "A", "t0", "s")
#'
#'     pnames <- dimnames(x)[[2]]
#'     ## weights on normal scale?
#'     x[,pnames %in% logit_scale_] <- exp(x[,pnames %in% logit_scale_])
#'     if(any(pnames %in% probit_scale_)) {
#'       x[,pnames %in% probit_scale_] <- pnorm(x[,pnames %in% probit_scale_])
#'     }
#'
#'     # scale alphas to set range [alphamin, 1]
#'     if('alpha1min' %in% pnames) x[,pnames == 'alpha1'] <- x[,pnames == 'alpha1']*(1-x[,pnames == 'alpha1min']) + x[,pnames == 'alpha1min']
#'     if('alpha2min' %in% pnames) x[,pnames == 'alpha2'] <- x[,pnames == 'alpha2']*(1-x[,pnames == 'alpha2min']) + x[,pnames == 'alpha2min']
#'     if('alpha3min' %in% pnames) x[,pnames == 'alpha3'] <- x[,pnames == 'alpha3']*(1-x[,pnames == 'alpha3min']) + x[,pnames == 'alpha3min']
#'
#'     #
#'     if('alpha1Posmin' %in% pnames) x[,pnames == 'alpha1Pos'] <- x[,pnames == 'alpha1Pos']*(1-x[,pnames == 'alpha1Posmin']) + x[,pnames == 'alpha1Posmin']
#'     if('alpha2Posmin' %in% pnames) x[,pnames == 'alpha2Pos'] <- x[,pnames == 'alpha2Pos']*(1-x[,pnames == 'alpha2Posmin']) + x[,pnames == 'alpha2Posmin']
#'     if('alpha3Posmin' %in% pnames) x[,pnames == 'alpha3Pos'] <- x[,pnames == 'alpha3Pos']*(1-x[,pnames == 'alpha3Posmin']) + x[,pnames == 'alpha3Posmin']
#'     if('alpha1Negmin' %in% pnames) x[,pnames == 'alpha1Neg'] <- x[,pnames == 'alpha1Neg']*(1-x[,pnames == 'alpha1Negmin']) + x[,pnames == 'alpha1Negmin']
#'     if('alpha2Negmin' %in% pnames) x[,pnames == 'alpha2Neg'] <- x[,pnames == 'alpha2Neg']*(1-x[,pnames == 'alpha2Negmin']) + x[,pnames == 'alpha2Negmin']
#'     if('alpha3Negmin' %in% pnames) x[,pnames == 'alpha3Neg'] <- x[,pnames == 'alpha3Neg']*(1-x[,pnames == 'alpha3Negmin']) + x[,pnames == 'alpha3Negmin']
#'
#'     x
#'   }
#'
#'   f <- function() return(model)
#'   environment(f) <- list2env(list(constantsFixed=constantsFixed, model=model,probit_scale=probit_scale, logit_scale=logit_scale), parent = globalenv())   # empty environment
#'   return(f)
#' }
#'
#'


### SM: RDM dynamic
# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' RL-RD with RDM as model
#'
#' @return
#' @export
#'
#' @examples
rdmDynamic <- function() {list(
  type="RACE",
  p_types=c("v","B","A","t0","s",
            "alpha1", "alpha2", "alpha3",
            "alpha1min", "alpha2min", "alpha3min",
            "q01", "q02", "q03",
            "weight1", "weight2", "weight3",
            paste0(rep('ucos', 10), 1:10),     # effect on URGENCY
            paste0(rep('Bcos', 10), 1:10),     # effect on THRESHOLD
            paste0(rep('vcos', 10), 1:10),     # effect on DRIFT
            paste0(rep('vtheta', 2), 1:2),     # exponential trend on drift
            paste0(rep('utheta', 2), 1:2),     # exponential trend on urgency
            paste0(rep('Btheta', 2), 1:2)      # exponential trend on threshold
  ),  # maximally 10 cosines

  Ntransform=function(x) {
    # Transform to natural scale
    probit_scale <- c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
    logit_scale = c("v", "B", "A", "t0", "s")

    if('v' %in% logit_scale) {
      ## set urgency signal v to minimally exp(-3) [~0.05] - the likelihood landscape is extremely flat below this value, causing sampling issues
      x[x[,'v'] < -3,'v'] <- 10  # set to extremely high value so likelihood always becomes 0
    }

    pnames <- dimnames(x)[[2]]
    ## weights on normal scale
    x[,pnames %in% logit_scale] <- exp(x[,pnames %in% logit_scale])
    if(any(pnames %in% probit_scale)) {
      x[,pnames %in% probit_scale] <- pnorm(x[,pnames %in% probit_scale])
    }

    # scale alphas to set range [alphamin, 1]
    if('alpha1min' %in% pnames) x[,pnames == 'alpha1'] <- x[,pnames == 'alpha1']*(1-x[,pnames == 'alpha1min']) + x[,pnames == 'alpha1min']
    if('alpha2min' %in% pnames) x[,pnames == 'alpha2'] <- x[,pnames == 'alpha2']*(1-x[,pnames == 'alpha2min']) + x[,pnames == 'alpha2min']
    if('alpha3min' %in% pnames) x[,pnames == 'alpha3'] <- x[,pnames == 'alpha3']*(1-x[,pnames == 'alpha3min']) + x[,pnames == 'alpha3min']

    x

  },
  trendTransform = function(npars, da, i) {
    index <- attr(da, 'adapt')[[i]]$index

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
  },
  DCTtransform =function(npars, da, i) {
    dct <- attr(da, 'adapt')[[i]]$dct
    index <- attr(da, 'adapt')[[i]]$index

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
  },
  polytransform =function(npars, da, i) {
    polyX <- attr(da, 'adapt')[[i]]$poly
    index <- attr(da, 'adapt')[[i]]$index

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
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                         levels(dadm$subjects))
    for (i in levels(dadm$subjects)) {
      ## transform cosines first
      npars <- rdmDynamic()$DCTtransform(npars=pars[dadm$subjects==i,], da=dadm[dadm$subjects==i,], i)

      ## potentially, polynomial basis functions
      npars <- rdmDynamic()$polytransform(npars=npars, da=dadm[dadm$subjects==i,], i)

      ## potentially, exponential trends
      npars <- rdmDynamic()$trendTransform(npars=npars, da=dadm[dadm$subjects==i,], i)

      ## and to updating next
      pars[dadm$subjects==i,] <- npars
      parsList[[i]] <- update_pars(i,pars,dadm)
    }

    pars <- do.call(rbind,parsList)
    pars[,'v'] <- pmax(pars[,'v'],0)
    pars[,'v'] <- pmin(pars[,'v'],30)
    pars[,'B'] <- pmax(pars[,'B'],0.1)
    pars[,'B'] <- pmin(pars[,'B'],30)
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) # & (pars[,'v'] > 1e-6) & (pars[,'v'] < 1e3) & (pars[,'B'] > .1) & (pars[,'alpha'] < 1)
    pars
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rRDM(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dRDM(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pRDM(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)
}


### SM: Generate RDM dynamic with DCT for slow waves
# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' RL-RD with RDM as model
#'
#' @return
#' @export
#'
#' @examples
generateRDMdynamic <- function(p_types, constantsFixed=NULL,
                               probit_scale=c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03"),
                               logit_scale=c("v", "B", "A", "t0", "s")) {
  model <- rdmDynamic()
  model$p_types <- p_types
  ## overwrite Ntransform
  model <- model[names(model)!='Ntransform']
  model$Ntransform=function(x,
                            constants=constantsFixed,
                            probit_scale_ = probit_scale,
                            logit_scale_ = logit_scale) {
    if(!is.null(constants)) {
      tmp <- matrix(rep(constants, each=nrow(x)), ncol=length(constants), dimnames=list(row.names(x), names(constants)))
      x <- cbind(x, tmp)
      # for(nm in names(constants)) {
      #   x <- cbind(x, nm=constants[[nm]])
      #   colnames(x)[ncol(x)] <- nm
      # }
    }


    if('v' %in% logit_scale_) {
      ## set urgency signal v to minimally exp(-4) [~0.019] -
      # the likelihood landscape is extremely flat below this value, causing sampling issues
      x[x[,'v'] < -4,'v'] <- 10  # set to extremely high value so likelihood always becomes 0
    }

    ## Transformations
    pnames <- dimnames(x)[[2]]
    x[,pnames %in% logit_scale_] <- exp(x[,pnames %in% logit_scale_])
    if(any(pnames %in% probit_scale_)) {
      x[,pnames %in% probit_scale_] <- pnorm(x[,pnames %in% probit_scale_])
    }

    # scale alphas to set range [alphamin, 1]
    if('alpha1min' %in% pnames) x[,pnames == 'alpha1'] <- x[,pnames == 'alpha1']*(1-x[,pnames == 'alpha1min']) + x[,pnames == 'alpha1min']
    if('alpha2min' %in% pnames) x[,pnames == 'alpha2'] <- x[,pnames == 'alpha2']*(1-x[,pnames == 'alpha2min']) + x[,pnames == 'alpha2min']
    if('alpha3min' %in% pnames) x[,pnames == 'alpha3'] <- x[,pnames == 'alpha3']*(1-x[,pnames == 'alpha3min']) + x[,pnames == 'alpha3min']

    #
    if('alpha1Posmin' %in% pnames) x[,pnames == 'alpha1Pos'] <- x[,pnames == 'alpha1Pos']*(1-x[,pnames == 'alpha1Posmin']) + x[,pnames == 'alpha1Posmin']
    if('alpha2Posmin' %in% pnames) x[,pnames == 'alpha2Pos'] <- x[,pnames == 'alpha2Pos']*(1-x[,pnames == 'alpha2Posmin']) + x[,pnames == 'alpha2Posmin']
    if('alpha3Posmin' %in% pnames) x[,pnames == 'alpha3Pos'] <- x[,pnames == 'alpha3Pos']*(1-x[,pnames == 'alpha3Posmin']) + x[,pnames == 'alpha3Posmin']
    if('alpha1Negmin' %in% pnames) x[,pnames == 'alpha1Neg'] <- x[,pnames == 'alpha1Neg']*(1-x[,pnames == 'alpha1Negmin']) + x[,pnames == 'alpha1Negmin']
    if('alpha2Negmin' %in% pnames) x[,pnames == 'alpha2Neg'] <- x[,pnames == 'alpha2Neg']*(1-x[,pnames == 'alpha2Negmin']) + x[,pnames == 'alpha2Negmin']
    if('alpha3Negmin' %in% pnames) x[,pnames == 'alpha3Neg'] <- x[,pnames == 'alpha3Neg']*(1-x[,pnames == 'alpha3Negmin']) + x[,pnames == 'alpha3Negmin']

    if('Btheta2' %in% pnames) x[,pnames == 'Btheta2'] <- x[,pnames == 'Btheta2']/4  # [0, 0.25]
    if('utheta2' %in% pnames) x[,pnames == 'utheta2'] <- x[,pnames == 'utheta2']/4  # [0, 0.25]
    if('vtheta2' %in% pnames) x[,pnames == 'vtheta2'] <- x[,pnames == 'vtheta2']/4  # [0, 0.25]

    x
  }

  f <- function() return(model)
  environment(f) <- list2env(list(constantsFixed=constantsFixed, model=model,probit_scale=probit_scale, logit_scale=logit_scale), parent = globalenv())   # empty environment
  return(f)
}

