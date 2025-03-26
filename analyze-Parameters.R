## parameter values
library(EMC2)
library(emcAdapt)
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMsNewEMC')
source(file.path(root_dir, 'extra_EMC2_functions/adaptive.R'))
source(file.path(root_dir, 'extra_EMC2_functions/model_RDMdynamic.R'))
source(file.path(root_dir, 'extra_EMC2_functions/utils.R'))


get_theta_mu_values <- function(samplers) {
  samples <- EMC2:::merge_chains(samplers)
  idx <- samples$samples$stage=='sample'

  transform <- samplers[[1]]$model()$transform
  theta_mu_ <- t(samples$samples$theta_mu[,idx])
  for(cname in colnames(theta_mu_)) {
    if(!cname %in% names(transform$func)) {
      transform$func[cname] <- 'identity'
      transform$upper[cname] <- Inf
      transform$lower[cname] <- -Inf
    }
  }
  theta_mu_transformed <- EMC2:::do_transform(theta_mu_, transform)
  theta_mu_transformed <- theta_mu_transformed[,c('alpha1', 'weight1', 'alpha2', 'weight2', 'alpha3', 'weight3')]
  colnames(theta_mu_transformed) <- c('alpha_SM', 'weight_SM', 'alpha_AM', 'weight_AM', 'alpha_FM', 'weight_FM')
  apply(theta_mu_transformed, 2, quantile, c(.025, .5, .975))
}


get_alpha_values <- function(samplers, return_range=TRUE) {
  transform <- samplers[[1]]$model()$transform
  alphas <- EMC2:::parameters(samplers, selection='alpha')
  for(cname in colnames(alphas)) {
    if(!cname %in% names(transform$func)) {
      transform$func[cname] <- 'identity'
      transform$upper[cname] <- Inf
      transform$lower[cname] <- -Inf
    }
  }

  subjs <- alphas[,1]
  alphas <- data.frame(EMC2:::do_transform(as.matrix(alphas[,grepl('(alpha|weight)', colnames(alphas))]), transform))
  alphas$subjects <- subjs

  if(return_range) {
    agg <- aggregate(.~subjects, alphas, median)[2:ncol(alphas)]
    return(apply(agg,2,range))
  } else {
    return(alphas)
  }
}

datasets <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
samplers <- lapply(datasets, loadDataSamplersPP, 'zSMuAHbV', pp_conditional=FALSE)

out = lapply(samplers, function(x) get_theta_mu_values(x$samplers))
names(out) <- c('dataset 1', 'dataset 2', 'dataset 3', 'dataset 4')
parameter_values <- t(do.call(cbind, lapply(out, data.frame)))
# write.csv(parameter_values, file='./tables/parameter_tables/parameter_values.csv')
# write.xlsx(round(parameter_values,7), file='./tables/parameter_tables/parameter_values.xlsx', colNames = TRUE, rowNames = TRUE)

# get range of medians per subject
out2 = lapply(samplers, function(x) get_alpha_values(x$samplers))
names(out2) <- c('dataset 1', 'dataset 2', 'dataset 3', 'dataset 4')
alpha_values <- t(do.call(cbind, lapply(out2, data.frame)))


# Make table for paper ----------------------------------------------------
tabl <- data.frame(matrix(NA, nrow=6, ncol=4))
for(i in 1:6) {
  for(ds in 1:4) {
    tabl[i,ds] <- paste0(round(out[[ds]][2,i], 2),  # median hyper
                         ' (', round(out[[ds]][1,i],2), ' - ', round(out[[ds]][3,i],2), ') [', # CI hyper
                         round(out2[[ds]][1,i],2),  ' - ', round(out2[[ds]][2,i],2), ']') # range median subject-level
  }
}
row.names(tabl) <- c('SM, alpha',
                     'SM, weight',
                     'AM, alpha',
                     'AM, weight',
                     'FM, alpha',
                     'FM, weight')
colnames(tabl) <- c('Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4')
write.csv(tabl, file='./tables/parameter_values_incl_range.csv')
write.xlsx(tabl, file='./tables/parameter_values_incl_range.xlsx', colNames = TRUE, rowNames = TRUE)



# Parameters: Dynamic vs static -------------------------------------------
process_npars <- function(npars) {
  npars$db <-  (npars[npars$lR==levels(npars$lR)[1],'B'] - npars[npars$lR==levels(npars$lR)[2],'B']) #npars$q01*npars$weight1
  npars$u <- (npars[npars$lR==levels(npars$lR)[1],'v'] + npars[npars$lR==levels(npars$lR)[2],'v'])/2
  npars$dr <- npars[npars$lR==levels(npars$lR)[1],'v'] - npars[npars$lR==levels(npars$lR)[2],'v']
  npars$dv <- npars[npars$lM==TRUE,'v'] - npars[npars$lM==FALSE,'v']
  npars$meanb <- (npars[npars$lR==levels(npars$lR)[1],'B'] + npars[npars$lR==levels(npars$lR)[2],'B'])/2
  npars$ds <- npars[npars$lM==TRUE,'s'] - npars[npars$lM==FALSE,'s']
  return(npars)
}
## null models
samplers0 <- lapply(datasets, loadDataSamplersPP, 'NULL', pp_conditional=FALSE, samples_dir='./samples_mc123')
npars0 <- lapply(samplers0, function(x) process_npars(attr(x$pp, 'npars'))) #process_npars(npars)
nullagg0 <- lapply(npars0, function(x) aggregate(cbind(u,dv,meanb,db)~subjects*trials, x, quantile, c(0.025, .5, .975)))
nullagg1 <- lapply(npars0, function(x) aggregate(cbind(t0=t0, ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, x, mean))  # across-trial mean
nullagg2 <- lapply(nullagg1, function(x) aggregate(.~subjects, x, quantile, c(0.025, 0.5, .975))) # across-trial mean


## dynamic models
nparsms3 <- lapply(samplers, function(x) process_npars(attr(x$pp, 'npars'))) #process_npars(npars)
ms3agg0 <- lapply(nparsms3, function(x) aggregate(cbind(u,dv,meanb,db)~subjects*trials, x, quantile, c(0.025, .5, .975)))
ms3agg1 <- lapply(nparsms3, function(x) aggregate(cbind(t0=t0, ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, x, mean))  # across-trial mean
ms3agg2 <- lapply(ms3agg1, function(x) aggregate(.~subjects, x, quantile, c(0.025, 0.5, .975))) # across-trial mean

pretty_labels <- list('dv'='dv', 'meanb'='b', 'u'='u', 't0'='t0', 'ds'='ds', 'db'='db',
                      'alpha1'=TeX('$\\alpha_{SM}$'),
                      'alpha2'=TeX('$\\alpha_{AM}$'),
                      'alpha3'=TeX('$\\alpha_{FM}$'),
                      'weight1'=TeX('$w_{SM}$'),
                      'weight2'=TeX('$w_{AM}$'),
                      'weight3'=TeX('$w_{FM}$'))
for(ftype in c('jpeg', 'pdf')) {
  fn = './figures/parameters_static-vs-dynamic.pdf'
  if(ftype == 'pdf') pdf(file='./figures/parameters_static-vs-dynamic.pdf', width=8, height=5.5)
  if(ftype == 'jpeg') jpeg(gsub('.pdf','.jpeg', fn), width=8, height=5.5, units='in', quality=100, res=500)

  par(mfrow=c(4,6), mar=c(2,2,1,.5), oma=c(1,2,1,0),bty='l', mgp=c(2,.6,0)) # loop over datasets then trials
  for(dataset in 1:4) {
    null <- nullagg2[[dataset]] #all_info[[dataset]]$agg0b1
    ms3 <- ms3agg2[[dataset]] #all_info[[dataset]]$agg1b1
    for(cn in c('dv', 'u', 'meanb', 'db', 't0', 'ds')) { #s_lMd', 'db')) {  #}, 'db')) {
      # xlab <- ifelse(dataset==4, 'IID-RDM', '')
      # ylab <- ifelse(cn=='dv', 'MS3-RDM', '')
      main <- ifelse(dataset == 1, pretty_labels[[cn]], '')
      xlim <- range(null[,cn])
      ylim <- range(ms3[,cn]) #, agg11[,cn])
      plot(null[,cn][,2], ms3[,cn][,2], main='', xlab='', ylab='', xlim=xlim, ylim=ylim)
      if(cn=='dv') {
        mtext(paste("Dataset", dataset), side=2, cex=.66*1.2, line=3, font=2)
        mtext('MS3-RDM', side=2, line=2, cex=.66)
      }
      if(dataset == 1) mtext(main, side=3, line=.5, cex=.66*1.2, font=2)
      if(dataset == 4) mtext('IID-RDM', side=1, line=1.5, cex=.66)
      if(!(cn=='db' & dataset %in% c(1,3))) abline(a=0,b=1)
      arrows(x0=null[,cn][,2], y0=ms3[,cn][,1], y1=ms3[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
      arrows(x0=null[,cn][,1], y0=ms3[,cn][,2], x1=null[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
    }
  }
  dev.off()
}



# Parameters: Trends vs no trends -----------------------------------------
## trend models
samplersT <- list(
  loadDataSamplersPP('wagenmakers2004_CS', 'zSMuAHbV', 'DCT', trendPar='B', nTrendPars=3, pp_conditional=FALSE, samples_dir='./samples_trends'),
  loadDataSamplersPP('forstmann2008', 'zSMuAHbV', 'DCT', trendPar='v', nTrendPars=3, pp_conditional=FALSE, samples_dir='./samples_trends'),
  loadDataSamplersPP('mileticvanmaanen2019exp2block2', 'zSMuAHbV', 'DCT', trendPar='v', nTrendPars=3, pp_conditional=FALSE, samples_dir='./samples_trends'),
  loadDataSamplersPP('wagenmakers2008exp2', 'zSMuAHbV', 'DCT', trendPar='B', nTrendPars=3, pp_conditional=FALSE, samples_dir='./samples_trends')
)
nparsT <- lapply(samplersT, function(x) process_npars(attr(x$pp, 'npars')))
trendagg1 <- lapply(nparsT, function(x) aggregate(cbind(t0, ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db),
                                                        alpha1, alpha2, alpha3, weight1, weight2, weight3)~subjects*postn, x, mean))  # across-trial mean
trendagg2 <- lapply(trendagg1, function(x) aggregate(.~subjects, x, quantile, c(0.025, 0.5, .975))) # across-trial mean

# add alphas
nparsms3 <- lapply(samplers, function(x) process_npars(attr(x$pp, 'npars')))
ms3agg1 <- lapply(nparsms3, function(x) aggregate(cbind(t0, ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db),
                                                        alpha1, alpha2, alpha3, weight1, weight2, weight3)~subjects*postn, x, mean))  # across-trial mean
ms3agg2 <- lapply(ms3agg1, function(x) aggregate(.~subjects, x, quantile, c(0.025, 0.5, .975))) # across-trial mean


for(ftype in c('jpeg', 'pdf')) {
  fn = './figures/parameters_trend_vs_notrend-1.pdf'
  if(ftype == 'pdf') pdf(file=fn, width=8, height=10)
  if(ftype == 'jpeg') jpeg(gsub('.pdf','.jpeg', fn), width=8, height=10, units='in', quality=100, res=500)
  par(mfrow=c(8,6), mar=c(2,2,1,.5), oma=c(1,2,1,0),bty='l', mgp=c(2,.6,0)) # loop over datasets then trials
  for(dataset in 1:4) {
    notrend <- ms3agg2[[dataset]]
    trend <- trendagg2[[dataset]]
    for(cn in c('dv', 'u', 'meanb', 'db', 't0', 'ds')) { #s_lMd', 'db')) {  #}, 'db')) {
      # xlab <- ifelse(dataset==4, 'No trends', '')
      # ylab <- ifelse(cn=='dv', 'Trends', '')
      main <- ifelse(dataset == 1, pretty_labels[[cn]], '')
      xlim <- range(notrend[,cn])
      ylim <- range(trend[,cn]) #, agg11[,cn])
      plot(notrend[,cn][,2], trend[,cn][,2], main='', xlab='', ylab='', xlim=xlim, ylim=ylim)
      if(cn=='dv') {
        mtext(paste("Dataset", dataset), side=2, cex=.66*1.2, line=3, font=2)
        mtext('Trends', side=2, line=2, cex=.66)
      }
      if(dataset == 1) mtext(main, side=3, line=.5, cex=.66*1.2, font=2)
      # if(dataset == 4) mtext(xlab, side=1, line=1.5, cex=.66)
      #if(!(cn=='db' & dataset %in% c(1,3)))
      abline(a=0,b=1)
      arrows(x0=notrend[,cn][,2], y0=trend[,cn][,1], y1=trend[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
      arrows(x0=notrend[,cn][,1], y0=trend[,cn][,2], x1=notrend[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
    }
  }

  for(dataset in 1:4) {
    notrend <- ms3agg2[[dataset]]
    trend <- trendagg2[[dataset]]
    for(cn in c('alpha1', 'alpha2', 'alpha3', 'weight1', 'weight2', 'weight3')) { #s_lMd', 'db')) {  #}, 'db')) {
      xlab <- ifelse(dataset==4, 'No trends', '')
      ylab <- ifelse(cn=='dv', 'Trends', '')
      main <- ifelse(dataset == 1, pretty_labels[[cn]], '')
      xlim <- range(notrend[,cn])
      ylim <- range(trend[,cn]) #, agg11[,cn])
      plot(notrend[,cn][,2], trend[,cn][,2], main='', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
      if(cn=='alpha1') {
        mtext(paste("Dataset", dataset), side=2, cex=.66*1.2, line=3, font=2)
        mtext('Trends', side=2, line=2, cex=.66)
      }
      if(dataset == 1) mtext(main, side=3, line=.5, cex=.66*1.2, font=2)
      if(dataset == 4) mtext(xlab, side=1, line=1.5, cex=.66)
      if(!(cn=='db' & dataset %in% c(1,3))) abline(a=0,b=1)
      arrows(x0=notrend[,cn][,2], y0=trend[,cn][,1], y1=trend[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
      arrows(x0=notrend[,cn][,1], y0=trend[,cn][,2], x1=notrend[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
    }
  }
  dev.off()
}



# Correlations between parameters -----------------------------------------
datasets <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
samplers <- lapply(datasets, loadDataSamplersPP, 'zSMuAHbV', pp_conditional=FALSE)


pname_mapping = list('v'='u',
                     'v_lMd'='dv',
                     'B'='b',
                     'B_Ea-n'='$b[a-n]',
                     'B_Ea-s'='$b[a-s]',
                     'B_lRd'='$b[l-r]',
                     't0'='$t[0]',
                     's_lMd'='ds',
                     'alpha3'='$alpha[FM]',
                     'weight3'='$w[FM]',
                     'alpha2'='$alpha[AM]',
                     'weight2'='$w[AM]',
                     'alpha1'='$alpha[SM]',
                     'weight1'='$w[SM]',
                     's_errorAccumulatorTRUE'='ds',
                     'v_lMTRUE'='dv',
                     'v_lMTRUE:Wvlf'='$dv:W[vlf]',
                     'v_lMTRUE:Whf'='$dv:W[hf]',
                     'v_lMTRUE:Wnonword'='$dv:W[hf]')
# pname_mapping <- lapply(pname_mapping, TeX)

# Plot correlation matrix
pdf(file='./figures/correlation_matrices.pdf', width=6, height=6)
par(mfrow=c(2,2))
EMC2:::plot_relations(samplers[[1]]$samplers, only_cred = T, nice_names = unlist(pname_mapping[samplers[[1]]$samplers[[1]]$par_names]), plot_cred = F, plot_means=F)
mtext('Dataset 1', side=3, at=par("usr")[1]+0.2*diff(par("usr")[1:2]))
EMC2:::plot_relations(samplers[[2]]$samplers, only_cred = T, nice_names = unlist(pname_mapping[samplers[[2]]$samplers[[1]]$par_names]), plot_cred = F, plot_means=F)
mtext('Dataset 2', side=3, at=par("usr")[1]+0.2*diff(par("usr")[1:2]))
EMC2:::plot_relations(samplers[[3]]$samplers, only_cred = T, nice_names = unlist(pname_mapping[samplers[[3]]$samplers[[1]]$par_names]), plot_cred = F, plot_means=F)
mtext('Dataset 3', side=3, at=par("usr")[1]+0.2*diff(par("usr")[1:2]))
EMC2:::plot_relations(samplers[[4]]$samplers, only_cred = T, nice_names = unlist(pname_mapping[samplers[[4]]$samplers[[1]]$par_names]), plot_cred = F, plot_means=F)
mtext('Dataset 4', side=3, at=par("usr")[1]+0.2*diff(par("usr")[1:2]))
dev.off()

## uuhhr okay, what to do wiht this..?





# OLD below ---------------------------------------------------------------
# Unconditioned on data ---------------------------------------------------
## unfortunately it's rather annoying to have post_predict also return the *updated* parameter values
# so instead we use the posterior predictive data here to reconstruct the parameter values by trial
# samplers <- samplers1
# pp1 <- get_pp('wagenmakers2004_CS', learningModel, dctModel = dctModel)$pp
# samplers0 <- get_pp('wagenmakers2004_CS', 'NULL', dctModel = dctModel)$samplers
# pp0 <- get_pp('wagenmakers2004_CS', 'NULL', dctModel = dctModel)$pp

# get_trialwise_parameters <- function(samplers, pp, subject, condition_on_data=FALSE) {
#   dadm <- samplers[[1]]$data[[subject]]
#   adapt <- attr(dadm, 'adapt')
#   pps <- pp[pp$subjects==subject,]
#
#   for(postn in unique(pps$postn)) {
#     p_vector <- t(attr(pp, 'pars')[[postn]][subject,])
#
#     dimnames(p_vector)[[1]] <- levels(dadm$subjects)
#     npars <- attr(dadm,"model")()$Ntransform(
#       EMC2:::map_p(
#         attr(dadm,"model")()$transform(EMC2:::add_constants(p_vector,attr(dadm,"constants"))),
#         dadm)
#     )
#     npars <- npars[order(dadm$trials),]   # order!
#
#
#     npars_match <- npars[seq(1,nrow(npars),2),] #dadm$lM==TRUE,]
#     pp_ <- pps[pps$postn==postn,]
#     if(condition_on_data) {
#       # enforce data RT and accuracy
#       dadm_ <- dadm[order(dadm$trials),]
#       dadm_ <- dadm_[dadm$winner,]
#       pp_$rt <- dadm_$rt
#       pp_$accuracy <- dadm_$accuracy
#     }
#
#     updated1 <- adaptb.c.emc(b0=npars_match[,'B'],
#                              q0=npars_match[1,'q03'],
#                              weight=npars_match[1,'weight3'],
#                              learningRates=npars_match[,'alpha3'],
#                              rts=pp_$rt)
#     q03 = updated1$adaptedValues[,2]
#     startValues <- npars_match[,c('q01', 'q02')]
#     learningRates <- npars_match[,c('alpha1', 'alpha2')]
#     outcomes <- as.matrix(pp_[,c('stim', 'accuracy')])
#
#     updated2 <- adapt.c.emc(feedback=outcomes,
#                             arguments=list(startValues = startValues,
#                                            learningRates = learningRates),
#                             learningRule='delta')
#
#     allQs <- matrix(NA, nrow=nrow(outcomes), ncol=3)
#     allQs[,3] <- updated1$adaptedValues[,2]   # Q3 from threshold updating
#     allQs[,1:2] <- updated2$adaptedValues[,1:2] # Q1, Q2 from delta rule
#
#     ## apply adaptive mapping
#     idx <- pps$subjects==subject&pps$postn==postn
#     pp_da <- EMC2:::add_accumulators(pp_)
#     if(!is.null(adapt$design$dynamic$output_fun)) {
#       b_v <- adapt$design$dynamic$output_fun(npars, allQs[rep(1:nrow(allQs), each=2),], pp_da)
#
#       # add to pps
#       pps[idx, c('q01', 'q02', 'q03')] <- allQs
#       pps[idx, c('B1', 'v1')] <- b_v[seq(1, nrow(b_v), 2),]
#       pps[idx, c('B2', 'v2')] <- b_v[seq(2, nrow(b_v), 2),]
#     } else {
#       pps[idx, c('B1', 'v1')] <- npars[seq(1, nrow(npars), 2), c('B', 'v')]
#       pps[idx, c('B2', 'v2')] <- npars[seq(2, nrow(npars), 2), c('B', 'v')]
#     }
#
#     pps[idx, 's1'] <- npars[seq(1, nrow(npars), 2), 's']
#     pps[idx, 's2'] <- npars[seq(2, nrow(npars), 2), 's']
#     pps[idx, dimnames(p_vector)[[2]]] <- matrix(p_vector[1,], nrow=sum(idx), ncol=length(p_vector), byrow=TRUE)
#
#   }
#   pps[,'u'] <- (pps$v1 + pps$v2)/2      # apply(cbind(pps$v1, pps$v2), 1, mean)
#   pps[,'dv'] <- pps$v1-pps$v2           # apply(pp2_with_parameters, 1, function(x) mean(c(x$v1, x$v2))) #$v1-pp2_with_parameters$v1
#   pps[,'meanb'] <- (pps$B1 + pps$B2)/2  # apply(cbind(pps$B1, pps$B2),1,mean)
#   pps[,'db'] <- pps$B1-pps$B2
#   pps[,'ds'] <- pps$s1-pps$s2
#   return(pps)
# }
#
#
#
#
# # For visualization of posterior parameters across trials ----------
# # condition on data
# get_pp2 <- function(task, learningModel, dctModel=NULL, get_pp=FALSE, decisionModel='RDM') {
#   ## NB: This *needs* to be conditioned on the data, otherwise it's not possible to know the subjective difficulty in trial t-1 in the *data*
#   # Load all ----------------------------------------------------------------
#   root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMs')
#   samples_dir = file.path(root_dir,  'samples')
#   figures_dir = file.path(root_dir, 'figures')
#
#   print(load(file.path(root_dir, 'datasets', paste0(task, '.RData'))))
#   save_fn_samples <- file.path(samples_dir, paste0(task, '_model-', decisionModel, '-', learningModel, '-DCT-', dctModel, '.RData'))
#   save_fn_pp <- sub('.RData', '_pp.RData', save_fn_samples)
#   if(!file.exists(save_fn_samples)) stop('Samples not found')
#   #if(!file.exists(save_fn_pp)) stop('pp not found')
#   load(save_fn_samples); nchains <- chain_n(samplers); filter=colnames(nchains)[nchains[1,]>0][sum(nchains[1,]>0)];
#   #  plot_chains(samplers, filter=filter, selection='mu', layout=c(3,3))
#
#   if(get_pp) {
#     load(save_fn_pp)
#     return(list(samplers=samplers,pp=pp))
#   }
#
#   #
#   for(samplerN in 1:length(samplers)) {
#     for(subject in length(samplers[[samplerN]]$data)) {
#       attr(samplers[[samplerN]]$data[[subject]], 'adapt')$design$simulateByTrial <- FALSE
#     }
#   }
#   attr(samplers, 'design_list')[[1]]$adapt$simulateByTrial <- FALSE
#
#   # undebug(EMC2:::update_pars)
#   # undebug(EMC2:::make_data)
#   pp2 <- post_predict(samplers, n_cores = 20, filter=filter) ## subfilter=nchains[1,4]-1000)
#   pp2 <- pp2[order(pp2$subjects, pp2$postn, pp2$trials),]
#   pp2$accuracy <- as.numeric(pp2$S) == as.numeric(pp2$R)       ## ensure logical
#   pp2$choice <- pp2$R==dat[dat$choice==1,'R'][1]
#
#   return(list(samplers=samplers, pp2=pp2))
# }
#
# plot_posterior_parameter <- function(cis, col=2, do_plot=FALSE, ...) {
#   if(do_plot) plot(cis[,2],type='l',col=col, ...) else lines(cis[,2], col=col)
#   polygon(c(1:nrow(cis), nrow(cis):1), c(cis[,3], rev(cis[,1])),col=adjustcolor(col,alpha.f=.2), border=NA)
# }
#
#
# get_fns <- function(task, learningModel, trendModel='DCT', trendPar='B', nTrendPars=3, samples_dir='./samples/') {
#   if(trendModel=='NULL') {
#     pp_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '-DCT-NULL_pp.RData')
#     samplers_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '-DCT-NULL.RData')
#   } else {
#     pp_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '-trend-', trendModel, '-', trendPar, '-', nTrendPars, '_pp.RData')
#     samplers_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '-trend-', trendModel, '-', trendPar, '-', nTrendPars, '.RData')
#   }
#   return(list(pp_fn=pp_fn, samplers_fn=samplers_fn))
# }
#
# get_model_comparison <- function(fns, subfilter=0) {
#   allSamplers <- list()
#   allGds <- c()
#   for(fn in fns) {
#     load(fn)
#     modelName <- strsplit(strsplit(fn, 'model-')[[1]][2], '.RData')[[1]][1]
#     nchains <- chain_n(samplers)
#     if(nchains[1,4]>=100) {
#       allSamplers[[modelName]] <- samplers
#       allGds <- c(allGds, max(EMC2:::gd_pmwg(samplers, subfilter=subfilter, print_summary = FALSE)))
#     }
#   }
#   if(length(allSamplers)>0) {
#     mc = EMC2:::compare_IC(allSamplers, print_summary=FALSE, subfilter=subfilter)
#     mc <- cbind(mc, max.gd=allGds)
#     print(round(mc,3))
#     mc
#   }
# }
#
# get_all_parameters <- function(task, dctModel='NULL') {
#   samplers0 <- get_pp(task, 'NULL', dctModel = dctModel)$samplers
#   samplers1 <- get_pp(task, 'zSMuAHbV', dctModel = dctModel)$samplers
#   pp0 <- get_pp(task, 'NULL', dctModel = dctModel)$pp
#   pp1 <- get_pp(task, 'zSMuAHbV', dctModel = dctModel)$pp
#   pp2 <- get_pp2(task, 'zSMuAHbV', dctModel = dctModel)$pp2
#
#   # get trend model
#   all_fns <- c(get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar='B', nTrendPars=3, samples_dir='./samples/')$samplers_fn,
#                get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar='v', nTrendPars=3, samples_dir='./samples/')$samplers_fn,
#                get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar='u', nTrendPars=3, samples_dir='./samples/')$samplers_fn)
#   mc <- get_model_comparison(all_fns)
#   winner <- row.names(mc)[which.max(mc$wBPIC)]
#   winningPar <- strsplit(winner, '-')[[1]][5]
#   fns <- get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar=winningPar, nTrendPars=3, samples_dir='./samples/')
#   pp3 <- EMC2:::loadRData(fns[[1]])
#   samplers3 <- EMC2:::loadRData(fns[[2]])
#
#   #augment with parameter values
#   pp01 <- do.call(rbind, mclapply(names(samplers0[[1]]$data), get_trialwise_parameters, samplers=samplers0, pp=pp0, mc.cores=length(samplers0[[1]]$data)))
#   pp11 <- do.call(rbind, mclapply(names(samplers0[[1]]$data), get_trialwise_parameters, samplers=samplers1, pp=pp1, mc.cores=length(samplers0[[1]]$data)))
#   pp21 <- do.call(rbind, mclapply(names(samplers0[[1]]$data), get_trialwise_parameters, samplers=samplers1, pp=pp2, condition_on_data=TRUE, mc.cores=length(samplers0[[1]]$data)))
#   pp31 <- do.call(rbind, mclapply(names(samplers0[[1]]$data), get_trialwise_parameters, samplers=samplers3, pp=pp3, mc.cores=length(samplers0[[1]]$data)))
#
#   # Aggregate across posterior predictives
#   agg0a <- aggregate(cbind(u,dv,meanb,db)~subjects*trials, pp01, quantile,c(.025, .5, .975))
#   agg1a <- aggregate(cbind(u,dv,meanb,db)~subjects*trials, pp11, quantile,c(.025, .5, .975))
#   agg2a <- aggregate(cbind(u,dv,meanb,db)~subjects*trials, pp21, quantile,c(.025, .5, .975))
#   agg3a <- aggregate(cbind(u,dv,meanb,db)~subjects*trials, pp31, quantile,c(.025, .5, .975))
#
#   # Aggregate across trials
#   agg0b = aggregate(cbind(t0=exp(t0), ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, pp01, mean)  # across-trial mean
#   agg1b = aggregate(cbind(t0=exp(t0), ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, pp11, mean)  # across-trial mean
#   agg2b = aggregate(cbind(t0=exp(t0), ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, pp21, mean)  # across-trial mean
#   agg3b = aggregate(cbind(t0=exp(t0), ds=abs(ds), u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, pp31, mean)  # across-trial mean
#
#   agg0b1 <- aggregate(.~subjects, agg0b, quantile, c(0.025, 0.5, .975))
#   agg1b1 <- aggregate(.~subjects, agg1b, quantile, c(0.025, 0.5, .975))
#   agg2b1 <- aggregate(.~subjects, agg2b, quantile, c(0.025, 0.5, .975))
#   agg3b1 <- aggregate(.~subjects, agg3b, quantile, c(0.025, 0.5, .975))
#
#   return(list(samplers0=samplers0,
#               pp01=pp01, pp11=pp11, pp21=pp21, pp31=pp31,
#               agg0a=agg0a, agg1a=agg1a, agg2a=agg2a, agg3a=agg3a,
#               agg0b=agg0b, agg1b=agg1b, agg2b=agg2b, agg3b=agg3b,
#               agg0b1=agg0b1, agg1b1=agg1b1, agg2b1=agg2b1, agg3b1=agg3b1
#   ))
# }
#
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2') #, 'mileticvanmaanen2019exp2block1')
# all_info <- mclapply(tasks, get_all_parameters, mc.cores=4)
#
# ## plots
# pretty_labels <- list('dv'='dv', 'meanb'='b', 'u'='u', 't0'='t0', 'ds'='ds', 'db'='db',
#                       'alpha1'='alpha (SM)', 'alpha2'='alpha (AM)', 'alpha3'='alpha (FM)',
#                       'weight1'='weight (SM)', 'weight2'='weight (AM)', 'weight3'='weight (FM)')
# for(ftype in c('jpeg', 'pdf')) {
#   fn = './figures/parameters_static-vs-dynamic.pdf'
#   if(ftype == 'pdf') pdf(file='./figures/parameters_static-vs-dynamic.pdf', width=8, height=5.5)
#   if(ftype == 'jpeg') jpeg(gsub('.pdf','.jpeg', fn), width=8, height=5.5, units='in', quality=100, res=500)
#
#
#   par(mfrow=c(4,6), mar=c(2,2,1,.5), oma=c(1,2,1,0),bty='l', mgp=c(2,.6,0)) # loop over datasets then trials
#   for(dataset in 1:4) {
#     agg01 <- all_info[[dataset]]$agg0b1
#     agg11 <- all_info[[dataset]]$agg1b1
#     for(cn in c('dv', 'u', 'meanb', 'db', 't0', 'ds')) { #s_lMd', 'db')) {  #}, 'db')) {
#       # xlab <- ifelse(dataset==4, 'IID-RDM', '')
#       # ylab <- ifelse(cn=='dv', 'MS3-RDM', '')
#       main <- ifelse(dataset == 1, pretty_labels[[cn]], '')
#       xlim <- range(agg01[,cn])
#       ylim <- range(agg11[,cn]) #, agg11[,cn])
#       plot(agg01[,cn][,2], agg11[,cn][,2], main='', xlab='', ylab='', xlim=xlim, ylim=ylim)
#       if(cn=='dv') {
#         mtext(paste("Dataset", dataset), side=2, cex=.66*1.2, line=3, font=2)
#         mtext('MS3-RDM', side=2, line=2, cex=.66)
#       }
#       if(dataset == 1) mtext(main, side=3, line=.5, cex=.66*1.2, font=2)
#       if(dataset == 4) mtext('IID-RDM', side=1, line=1.5, cex=.66)
#       if(!(cn=='db' & dataset %in% c(1,3))) abline(a=0,b=1)
#       arrows(x0=agg01[,cn][,2], y0=agg11[,cn][,1], y1=agg11[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#       arrows(x0=agg01[,cn][,1], y0=agg11[,cn][,2], x1=agg01[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#     }
#   }
#   dev.off()
# }
#
#
# # add alphas and weights
# samplers1DCT <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMuAHbV-trend-DCT-B-3.RData')
# samplers2DCT <- EMC2:::loadRData('./samples/forstmann2008_model-RDM-zSMuAHbV-trend-DCT-v-3.RData') #get_pp('forstmann2008', learningModel, dctModel = dctModel)$samplers
# samplers3DCT <- EMC2:::loadRData('./samples/mileticvanmaanen2019exp2block2_model-RDM-zSMuAHbV-trend-DCT-u-3.RData') #get_pp('mileticvanmaanen2019exp2block2', learningModel, dctModel = dctModel)$samplers
# samplers4DCT <- EMC2:::loadRData('./samples/wagenmakers2008exp2_model-RDM-zSMuAHbV-trend-DCT-B-3.RData') #get_pp('wagenmakers2008exp2', learningModel, dctModel = dctModel)$samplers
#
# out_notrend = lapply(list(samplers1, samplers2, samplers3, samplers4), get_alpha_values, return_range=FALSE)
# out_trend = lapply(list(samplers1DCT, samplers2DCT, samplers3DCT, samplers4DCT), get_alpha_values, return_range=FALSE)
# names(out_trend) <- names(out_notrend) <- c('dataset 1', 'dataset 2', 'dataset 3', 'dataset 4')
#
# out_trend <- lapply(out_trend, function(x) aggregate(.~subjects, x, quantile, c(0.025, .5, .975)))
# out_notrend <- lapply(out_notrend, function(x) aggregate(.~subjects, x, quantile, c(0.025, .5, .975)))
#
# # Trend vs no trend -------------------------------------------------------
# for(ftype in c('jpeg', 'pdf')) {
#   fn = './figures/parameters_trend_vs_notrend-1.pdf'
#   if(ftype == 'pdf') pdf(file=fn, width=8, height=10)
#   if(ftype == 'jpeg') jpeg(gsub('.pdf','.jpeg', fn), width=8, height=10, units='in', quality=100, res=500)
#   par(mfrow=c(8,6), mar=c(2,2,1,.5), oma=c(1,2,1,0),bty='l', mgp=c(2,.6,0)) # loop over datasets then trials
#   for(dataset in 1:4) {
#     agg01 <- all_info[[dataset]]$agg1b1
#     agg11 <- all_info[[dataset]]$agg3b1
#     for(cn in c('dv', 'u', 'meanb', 'db', 't0', 'ds')) { #s_lMd', 'db')) {  #}, 'db')) {
#       # xlab <- ifelse(dataset==4, 'No trends', '')
#       # ylab <- ifelse(cn=='dv', 'Trends', '')
#       main <- ifelse(dataset == 1, pretty_labels[[cn]], '')
#       xlim <- range(agg01[,cn])
#       ylim <- range(agg11[,cn]) #, agg11[,cn])
#       plot(agg01[,cn][,2], agg11[,cn][,2], main='', xlab='', ylab='', xlim=xlim, ylim=ylim)
#       if(cn=='dv') {
#         mtext(paste("Dataset", dataset), side=2, cex=.66*1.2, line=3, font=2)
#         mtext('Trends', side=2, line=2, cex=.66)
#       }
#       if(dataset == 1) mtext(main, side=3, line=.5, cex=.66*1.2, font=2)
#       # if(dataset == 4) mtext(xlab, side=1, line=1.5, cex=.66)
#       #if(!(cn=='db' & dataset %in% c(1,3)))
#       abline(a=0,b=1)
#       arrows(x0=agg01[,cn][,2], y0=agg11[,cn][,1], y1=agg11[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#       arrows(x0=agg01[,cn][,1], y0=agg11[,cn][,2], x1=agg01[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#     }
#   }
#   # dev.off()
#   # }
#
#   #
#   # for(ftype in c('jpeg', 'pdf')) {
#   #   fn = './figures/parameters_trend_vs_notrend-2.pdf'
#   #   if(ftype == 'pdf') pdf(file=fn, width=8, height=5.5)
#   #   if(ftype == 'jpeg') jpeg(gsub('.pdf','.jpeg', fn), width=8, height=5.5, units='in', quality=100, res=500)
#   # par(mfrow=c(4,6), mar=c(2,2,1,.5), oma=c(1,2,1,0),bty='l', mgp=c(2,.6,0)) # loop over datasets then trials
#   for(dataset in 1:4) {
#     agg01 <- out_notrend[[dataset]] # X axis: dynamic without trends
#     agg11 <- out_trend[[dataset]]  # y axis: dynamic with trends
#     for(cn in c('alpha1', 'alpha2', 'alpha3', 'weight1', 'weight2', 'weight3')) { #s_lMd', 'db')) {  #}, 'db')) {
#       xlab <- ifelse(dataset==4, 'No trends', '')
#       ylab <- ifelse(cn=='dv', 'Trends', '')
#       main <- ifelse(dataset == 1, pretty_labels[[cn]], '')
#       xlim <- range(agg01[,cn])
#       ylim <- range(agg11[,cn]) #, agg11[,cn])
#       plot(agg01[,cn][,2], agg11[,cn][,2], main='', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
#       if(cn=='alpha1') {
#         mtext(paste("Dataset", dataset), side=2, cex=.66*1.2, line=3, font=2)
#         mtext('Trends', side=2, line=2, cex=.66)
#       }
#       if(dataset == 1) mtext(main, side=3, line=.5, cex=.66*1.2, font=2)
#       if(dataset == 4) mtext(xlab, side=1, line=1.5, cex=.66)
#       if(!(cn=='db' & dataset %in% c(1,3))) abline(a=0,b=1)
#       arrows(x0=agg01[,cn][,2], y0=agg11[,cn][,1], y1=agg11[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#       arrows(x0=agg01[,cn][,1], y0=agg11[,cn][,2], x1=agg01[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#     }
#   }
#   dev.off()
# }
#
#
# plot_trialwise_subject <- function(dat, ci0, ci1, xlim=c(100,150), xlab='Trials') {
#   #  dat <- dat[dat$trials %in% xlim[1]:xlim[2],]
#   #  idx1 <- ci0$trials %in% xlim[1]:xlim[2]
#   #  idx2 <- ci1$trials %in% xlim[1]:xlim[2]
#   idx <- dat$trials %in% xlim[1]:xlim[2]
#   for(param in c('db', 'u', 'meanb')) {
#     ylim <- range(rbind(ci0[,param],ci1[,param]))
#     plot_posterior_parameter(ci0[,param], do_plot=TRUE, ylim=ylim, col=1, xlab='', ylab=pretty_labels[[param]], xlim=xlim, xaxt='n', main='')
#     mtext(xlab, side=1, line=1, cex=par()$cex*par()$cex.lab)
#     plot_posterior_parameter(ci1[,param], do_plot=FALSE)
#     axis(side=1, at=pretty(seq(0, 1000)), labels=NULL)
#
#     if(param == 'db') points(dat[idx, 'trials'], dat[idx,]$stim*c(ylim[1]*.9), pch=4)
#     if(param == 'u') abline(v=dat[idx&dat$accuracy==0,'trials'],lty=2,col=1)
#   }
# }
#
#
# ##
# dataset <- 1
# samplers0 <- all_info[[dataset]]$samplers0
# agg0 <- all_info[[dataset]]$agg0a
# agg1 <- all_info[[dataset]]$agg2a
#
# subject = 1
# dat <- samplers0[[1]]$data[[subject]]; dat<-dat[order(dat$trials),]; dat<-dat[dat$winner,]
# par(mfrow=c(2,3), mar=c(3,3,2,1), bty='l', mgp=c(2,1,0))
# plot_trialwise_subject(dat=dat, ci0=agg0[agg0$subjects==subject,],ci1=agg1[agg1$subjects==subject,])
#
# ## single subject for paper
# for(colors_ in c('black', 'white')) {
#   if(colors_ == 'black') {
#     pdf(file=paste0('./figures/parameters_by_trial_by_subject_ds-', dataset, '_subject-', subject, '_black.pdf'), width=6, height=2.5)
#     if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white')
#     cols <- palette()
#     cols[1] <- 'white'
#     palette(cols)
#   } else {
#     cols <- palette()
#     cols[1] <- 'black'
#     palette(cols)
#     pdf(file=paste0('./figures/parameters_by_trial_by_subject_ds-', dataset, '_subject-', subject, '.pdf'), width=9, height=2.5)
#   }
#   par(mfrow=c(1,3), mar=c(5,4,4,1), bty='l', mgp=c(2,1,0))
#   plot_trialwise_subject(dat=dat, ci0=agg0[agg0$subjects==subject,],ci1=agg1[agg1$subjects==subject,])
#   dev.off()
# }
#
#
# ## trials
# pdf(file=paste0('./figures/dataset-',task,'_parameters_by_trial_by_subject.pdf'), width=12, height=6)
# par(mfrow=c(2,3))
# xlim <- c(100, 150)
# for(subject in unique(agg0$subjects)) { #[1:10]) {
#   dat <- samplers0[[1]]$data[[subject]]; dat<-dat[order(dat$trials),]; dat<-dat[dat$winner,]
#
#   ylim <- range(rbind(agg0[agg0$subjects==subject,'db'],agg2[agg2$subjects==subject,'db']))
#   plot_posterior_parameter(agg0[agg0$subjects==subject,'db'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='db')
#   plot_posterior_parameter(agg2[agg2$subjects==subject,'db'], do_plot=FALSE)
#   points(dat$trials, dat$stim*ylim[1], pch=4)
#
#   ylim <- range(rbind(agg0[agg0$subjects==subject,'meanb'],agg2[agg2$subjects==subject,'meanb']))
#   plot_posterior_parameter(agg0[agg0$subjects==subject,'meanb'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='meanb')
#   plot_posterior_parameter(agg2[agg2$subjects==subject,'meanb'], do_plot=FALSE)
#
#   ylim <- range(rbind(agg0[agg0$subjects==subject,'u'],agg2[agg2$subjects==subject,'u']))
#   plot_posterior_parameter(agg0[agg0$subjects==subject,'u'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='u')
#   plot_posterior_parameter(agg2[agg2$subjects==subject,'u'], do_plot=FALSE)
#   abline(v=dat[dat$accuracy==0,'trials'],lty=2,col='grey')
# }
# dev.off()
#
#
#
#
# #### Combine everything
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2', 'mileticvanmaanen2019exp2block1')
# #task <- tasks[3]
# for(task in tasks) {
#   samplers0 <- get_pp(task, 'NULL', dctModel = dctModel)$samplers
#   samplers1 <- get_pp(task, 'zSMuAHbV', dctModel = dctModel)$samplers
#   pp0 <- get_pp(task, 'NULL', dctModel = dctModel)$pp
#   pp1 <- get_pp(task, 'zSMuAHbV', dctModel = dctModel)$pp
#   pp2 <- get_pp2(task, 'zSMuAHbV', dctModel = dctModel)$pp2
#
#   #augment with parameter values
#   pp01 <- do.call(rbind, mclapply(names(samplers0[[1]]$data), get_trialwise_parameters, samplers=samplers0, pp=pp0, mc.cores=length(samplers0[[1]]$data)))
#   pp11 <- do.call(rbind, mclapply(names(samplers0[[1]]$data), get_trialwise_parameters, samplers=samplers1, pp=pp1, mc.cores=length(samplers0[[1]]$data)))
#   pp21 <- do.call(rbind, mclapply(names(samplers0[[1]]$data), get_trialwise_parameters, samplers=samplers1, pp=pp2, condition_on_data=TRUE, mc.cores=length(samplers0[[1]]$data)))
#
#   agg0 <- aggregate(cbind(u,dv,meanb,db)~subjects*trials, pp01, quantile,c(.025, .5, .975))
#   agg1 <- aggregate(cbind(u,dv,meanb,db)~subjects*trials, pp11, quantile,c(.025, .5, .975))
#   agg2 <- aggregate(cbind(u,dv,meanb,db)~subjects*trials, pp21, quantile,c(.025, .5, .975))
#
#   #  agg2 <- aggregate(cbind(u,dv,meanb,db,q01,q02,q03)~subjects*trials, pp21, quantile,c(.025, .5, .975))
#
#   #
#   # subject <- 2; xlim=c(250,300)
#   # dat[dat$subjects==subject,][250:300,]
#
#   pdf(file=paste0('./figures/dataset-',task,'_parameters_by_trial_by_subject.pdf'), width=12, height=6)
#   par(mfrow=c(2,3))
#   xlim <- c(100, 150)
#   for(subject in unique(agg0$subjects)) { #[1:10]) {
#     dat <- samplers0[[1]]$data[[subject]]; dat<-dat[order(dat$trials),]; dat<-dat[dat$winner,]
#
#     ylim <- range(rbind(agg0[agg0$subjects==subject,'db'],agg2[agg2$subjects==subject,'db']))
#     plot_posterior_parameter(agg0[agg0$subjects==subject,'db'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='db')
#     plot_posterior_parameter(agg2[agg2$subjects==subject,'db'], do_plot=FALSE)
#     points(dat$trials, dat$stim*ylim[1], pch=4)
#
#     ylim <- range(rbind(agg0[agg0$subjects==subject,'meanb'],agg2[agg2$subjects==subject,'meanb']))
#     plot_posterior_parameter(agg0[agg0$subjects==subject,'meanb'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='meanb')
#     plot_posterior_parameter(agg2[agg2$subjects==subject,'meanb'], do_plot=FALSE)
#
#     ylim <- range(rbind(agg0[agg0$subjects==subject,'u'],agg2[agg2$subjects==subject,'u']))
#     plot_posterior_parameter(agg0[agg0$subjects==subject,'u'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='u')
#     plot_posterior_parameter(agg2[agg2$subjects==subject,'u'], do_plot=FALSE)
#     abline(v=dat[dat$accuracy==0,'trials'],lty=2,col='grey')
#
#     #ylim <- range(rbind(agg0[agg0$subjects==subject,'dv'],agg2[agg2$subjects==subject,'dv']))
#     #plot_posterior_parameter(agg0[agg0$subjects==subject,'dv'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='dv')
#     #plot_posterior_parameter(agg2[agg2$subjects==subject,'dv'], do_plot=FALSE)
#
#
#     # ylim <- range(agg2[agg2$subjects==subject,'q02'])
#     # plot_posterior_parameter(agg0[agg0$subjects==subject,'q02'], do_plot=TRUE, xlim=xlim, ylim=ylim, col=1, xlab='Trials', ylab='u')
#     # plot_posterior_parameter(agg2[agg2$subjects==subject,'q02'], do_plot=FALSE)
#     # abline(v=dat[dat$accuracy==0,'trials'],lty=2,col='grey')
#   }
#   dev.off()
#
#   ## Compare null parameters with dynamic parameters
#   agg0 = aggregate(cbind(t0=exp(t0), ds, u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, pp01, mean)  # across-trial mean
#   agg1 = aggregate(cbind(t0=exp(t0), ds, u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, pp11, mean)  # across-trial mean
#   agg2 = aggregate(cbind(t0=exp(t0), ds, u, dv=abs(dv), meanb, db=abs(db))~subjects*postn, pp21, mean)  # across-trial mean
#
#   agg01 <- aggregate(.~subjects, agg0, quantile, c(0.025, 0.5, .975))
#   agg11 <- aggregate(.~subjects, agg1, quantile, c(0.025, 0.5, .975))
#   agg21 <- aggregate(.~subjects, agg2, quantile, c(0.025, 0.5, .975))
#
#   pdf(file=paste0('./figures/dataset-', task, '_parameter_comparison.pdf'), width=6, height=6)
#   par(mfrow=c(2,3), mar=c(4,3,3,1), bty='l', mgp=c(2,1,0))
#   for(cn in c('dv', 'u', 'meanb', 't0', 'db', 'ds')) { #s_lMd', 'db')) {  #}, 'db')) {
#     plot(agg01[,cn][,2], agg11[,cn][,2], main=cn, xlab='Null', ylab='Dynamic')
#     abline(a=0,b=1)
#     arrows(x0=agg01[,cn][,2], y0=agg11[,cn][,1], y1=agg11[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#     arrows(x0=agg01[,cn][,1], y0=agg11[,cn][,2], x1=agg01[,cn][,3], length=.025, angle=90, code=3, col=adjustcolor(1, alpha.f=.2))
#   }
#   dev.off()
# }
#
