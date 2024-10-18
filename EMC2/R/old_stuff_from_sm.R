#'
#' ### SM: RDM dynamic
#' # RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' # Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' #' RL-RD with RDM as model
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' rdmDynamic <- function() {list(
#'   type="RACE",
#'   p_types=c("v","B","A","t0","s",
#'             "alpha1", "alpha2", "alpha3",
#'             "q01", "q02", "q03",
#'             "weight1", "weight2", "weight3"),
#'
#'   Ntransform=function(x) {
#'     # Transform to natural scale
#'     probit_scale <- c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     logit_scale = c("v", "B", "A", "t0", "s")
#'     ## weights on normal scale?
#'     x[,dimnames(x)[[2]] %in% logit_scale] <- exp(x[,dimnames(x)[[2]] %in% logit_scale])
#'     x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
#'     x
#'   },
#'   # Trial dependent parameter transform
#'   Ttransform = function(pars,dadm) {
#'     parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
#'                          levels(dadm$subjects))
#'     for (i in levels(dadm$subjects))
#'       parsList[[i]] <- update_pars(i,pars,dadm)
#'     pars <- do.call(rbind,parsList)
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
#'
#' ### SM: RDM dynamic2
#' # RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
#' # Added q0 parameter for initial rates in RL version. Note that this is the *SAME* code as RLARD except for the p_types :-/ that should be easier!
#' #' RL-RD with RDM as model
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' rdmDynamic2 <- function() {list(
#'   type="RACE",
#'   p_types=c("v","B","A","t0","s", "alphamin",
#'             "alpha1", "alpha2", "alpha3",
#'             "q01", "q02", "q03",
#'             "weight1", "weight2", "weight3"),
#'
#'   Ntransform=function(x) {
#'     # Transform to natural scale
#'     probit_scale <- c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     logit_scale = c("v", "B", "A", "t0", "s")
#'     ## weights on normal scale?
#'     x[,dimnames(x)[[2]] %in% logit_scale] <- exp(x[,dimnames(x)[[2]] %in% logit_scale])
#'     x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
#'
#'     # scale alpha1 to set range [alphamin, 1]
#'     x[,dimnames(x)[[2]] == 'alpha1'] <- x[,dimnames(x)[[2]] == 'alpha1']*(1-x[,dimnames(x)[[2]] == 'alphamin']) + x[,dimnames(x)[[2]] == 'alphamin']
#'     x
#'   },
#'   # Trial dependent parameter transform
#'   Ttransform = function(pars,dadm) {
#'     parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
#'                          levels(dadm$subjects))
#'     for (i in levels(dadm$subjects))
#'       parsList[[i]] <- update_pars(i,pars,dadm)
#'     pars <- do.call(rbind,parsList)
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
#'
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
#' rdmDynamicDCT <- function() {list(
#'   type="RACE",
#'   p_types=c("v","B","A","t0","s",
#'             "alpha1", "alpha2", "alpha3",
#'             "q01", "q02", "q03",
#'             "weight1", "weight2", "weight3",
#'             paste0(rep('wcos', 10), 1:10)),  # maximally 10 cosines, should really be enough
#'
#'   Ntransform=function(x) {
#'     # Transform to natural scale
#'     probit_scale <- c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     logit_scale = c("v", "B", "A", "t0", "s")
#'     ## weights on normal scale?
#'     x[,dimnames(x)[[2]] %in% logit_scale] <- exp(x[,dimnames(x)[[2]] %in% logit_scale])
#'     x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
#'     x
#'   },
#'   DCTtransform =function(npars, da, i) {
#'     dct <- attr(da, 'adapt')[[i]]$dct
#'     index <- attr(da, 'adapt')[[i]]$index
#'
#'     wcos <- npars[1,grepl('wcos', colnames(npars))]
#'     nwcos <- length(wcos) # sum(grepl('wcos', colnames(npars)))
#'     ndct <- ncol(dct)
#'
#'     # ensure correct sizes
#'     if(nwcos < ndct) dct <- dct[,1:nwcos]     # less cosines estimated than provided
#'     if(nwcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nwcos-ndct))  # more cosines estimated than provided
#'
#'     ## drift in trial order
#'     drift <- dct %*% wcos
#'
#'     ## to "EMC-order"
#'     drift <- drift[index]
#'     npars[,'v'] <- npars[,'v'] + drift
#'     return(npars)
#'   },
#'   # Trial dependent parameter transform
#'   Ttransform = function(pars,dadm) {
#'     parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
#'                          levels(dadm$subjects))
#'     for (i in levels(dadm$subjects)) {
#'       ## OLD cosines here
#'       # npars <- pars[dadm$subjects==i,]
#'       # da <- dadm[dadm$subjects==i,]
#'
#'       # dct <- attr(da, 'adapt')[[i]]$dct
#'       # index <- attr(da, 'adapt')[[i]]$index
#'       #
#'       # nwcos <- sum(grepl('wcos', colnames(npars)))
#'       # ncdt <- ncol(dct)
#'       # dct <- dct[,1:nwcos]
#'       #
#'       # ## drift in trial order
#'       # drift <- dct %*% npars[1,grepl('wcos', colnames(npars))]
#'       #
#'       # ## to "EMC-order"
#'       # drift <- drift[index]
#'       # npars[,'v'] <- npars[,'v'] + drift
#'
#'       ## transform cosines first
#'       npars <- rdmDynamicDCT()$DCTtransform(npars=pars[dadm$subjects==i,], da=dadm[dadm$subjects==i,], i)
#'
#'       ## and to updating next
#'       pars[dadm$subjects==i,] <- npars
#'       parsList[[i]] <- update_pars(i,pars,dadm)
#'     }
#'
#'     pars <- do.call(rbind,parsList)
#'     pars[,'v'] <- pmax(pars[,'v'],0)
#'     pars[,'v'] <- pmin(pars[,'v'],30)
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
#' rdmDynamicDCT2 <- function() {list(
#'   type="RACE",
#'   p_types=c("v","B","A","t0","s",
#'             "alpha1", "alpha2", "alpha3",
#'             "q01", "q02", "q03",
#'             "weight1", "weight2", "weight3",
#'             paste0(rep('vcos', 10), 1:10),    # effect on DRIFT
#'             paste0(rep('Bcos', 10), 1:10)     # effect on threshold
#'   ),  # maximally 10 cosines, should really be enough
#'
#'   Ntransform=function(x) {
#'     # Transform to natural scale
#'     probit_scale <- c('alpha1', "alpha2", "alpha3", "q01", "q02", "q03")
#'     logit_scale = c("v", "B", "A", "t0", "s")
#'     ## weights on normal scale?
#'     x[,dimnames(x)[[2]] %in% logit_scale] <- exp(x[,dimnames(x)[[2]] %in% logit_scale])
#'     x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
#'     x
#'   },
#'   DCTtransform =function(npars, da, i) {
#'     dct <- attr(da, 'adapt')[[i]]$dct
#'     index <- attr(da, 'adapt')[[i]]$index
#'
#'     vcos <- npars[1,grepl('vcos', colnames(npars))]
#'     Bcos <- npars[1,grepl('Bcos', colnames(npars))]
#'     nvcos <- length(vcos)
#'     ndct <- ncol(dct)
#'
#'     # ensure correct sizes
#'     if(nvcos < ndct) dct <- dct[,1:nvcos]     # less cosines estimated than provided
#'     if(nvcos > ndct) dct <- cbind(dct, matrix(0, nrow=nrow(dct), ncol=nvcos-ndct))  # more cosines estimated than provided
#'
#'     ## drifts in trial order
#'     vdrift <- dct %*% vcos
#'     Bdrift <- dct %*% Bcos
#'
#'     ## to "EMC-order"
#'     npars[,'v'] <- npars[,'v'] + vdrift[index]
#'     npars[,'B'] <- npars[,'B'] + Bdrift[index]
#'     return(npars)
#'   },
#'   # Trial dependent parameter transform
#'   Ttransform = function(pars,dadm) {
#'     parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
#'                          levels(dadm$subjects))
#'     for (i in levels(dadm$subjects)) {
#'       ## transform cosines first
#'       npars <- rdmDynamicDCT2()$DCTtransform(npars=pars[dadm$subjects==i,], da=dadm[dadm$subjects==i,], i)
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
