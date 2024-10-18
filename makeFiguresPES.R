## Figures for paper
rm(list=ls())
#source('./utils_plotting.R')
figures_dir = './figures'

add_trial_nrs <- function(dat) {
  dat$trial_nr <- NA
  for(s in unique(dat$subjects)) {
    dat[dat$subjects==s, 'trials'] <- 1:(sum(dat$subjects==s))
  }
  return(dat)
}

loadDataAndPP <- function(task, learningModel, dctModel, samples_dir='samples') {
  # Load data, posterior predictives, samples
  load(paste0('./datasets/', task, '.RData'))
  load(file.path(samples_dir, paste0(task, '_model-RDM-', learningModel, '-DCT-', dctModel, '_pp.RData')))
  load(file.path(samples_dir, paste0(task, '_model-RDM-', learningModel, '-DCT-', dctModel, '.RData')))

  # ensure pp and data is correctly ordered
  pp <- pp[order(pp$subjects,pp$postn,pp$trials),]
  pp$accuracy <- as.numeric(pp$S) == as.numeric(pp$R)       ## ensure logical
  dat$accuracy <- as.numeric(dat$S) == as.numeric(dat$R)    ## ensure logical
  dat$choice <- dat$R == levels(dat$R)[1]
  pp$choice <- pp$R==dat[dat$choice==1,'R'][1]

  return(list(dat=dat, pp=pp, samplers=samplers))
}

getPostPreDifference <- function(dat, rtcolname='rt', acccolname='accuracy') {
  errors <- dat$accuracy == 0
  ## remove first trial (per participant) and last trial (per participant)
  if(!'trials' %in% colnames(dat)) dat <- add_trial_nrs(dat)
  trialRangeBySub <- aggregate(trials~subjects, dat, range)

  dat$exclude_trial <- FALSE
  for(subject in unique(trialRangeBySub$subjects)) {
    ## don't include errors on the first and last trial; can't calculate post-error slowing because of lack of a pre / post error trial
    dat[dat$subjects==subject,'exclude_trial'] = dat[dat$subjects==subject,'trials'] %in% trialRangeBySub[trialRangeBySub$subjects==subject,'trials']
  }
  errors <- errors & !dat$exclude_trial
  errors <- which(errors)
  preError <- errors-1

  ## moving accuracy
  filterOrder <- 10
  movingaverage <- stats::filter(dat[,'accuracy'], rep(1,filterOrder), sides = 1)/filterOrder #meanAccuracy <- mean(dat[,acccolname])
  movingaverage[1] <- dat[1,'accuracy']
  movingaveragert <- stats::filter(dat[,'rt'], rep(1,filterOrder), sides = 1)/filterOrder #meanAccuracy <- mean(dat[,acccolname])
  movingaveragert[1] <- dat[1,'rt']
  for(i in 2:filterOrder) {
    movingaverage[i] <- (stats::filter(dat[,'accuracy'], rep(1,i), sides = 1)/i)[i]
    movingaveragert[i] <- (stats::filter(dat[,'rt'], rep(1,i), sides = 1)/i)[i]
  }


  f <- mean #f <- function(x, ...) quantile(x, probs=c(.1, .5, .9), ...)
  rtdf <- expand.grid(trialNposterror=1:5, probs=c(.1, .5, .9), measure='rt', value=NA)
  rtabsdf <- expand.grid(trialNposterror=seq(-3,7), probs='mean', measure='rtabs', value=NA) #c(.1, .5, .9), measure='rtabs', value=NA)
  accabsdf <- expand.grid(trialNposterror=seq(-3,7), probs='mean', measure='accabs', value=NA)
  accdf <- expand.grid(trialNposterror=1:5, probs='mean', measure='accuracy', value=NA)
  for(trialNposterror in 1:5) {
    rtdf[rtdf$trialNposterror==trialNposterror,'value'] <- f(dat[errors+trialNposterror,rtcolname] - movingaveragert[preError], na.rm=TRUE)
    accdf[accdf$trialNposterror==trialNposterror,'value'] <- mean(dat[errors+trialNposterror,acccolname] - movingaverage[preError], na.rm=TRUE)
  }
  for(trialNposterror in seq(-3,7)) {
    trialIdx <- errors+trialNposterror
    trialIdx <- trialIdx[trialIdx>0]  ## remove
    rtabsdf[rtabsdf$trialNposterror==trialNposterror,'value'] <- f(dat[trialIdx,rtcolname], na.rm=TRUE)
    accabsdf[accabsdf$trialNposterror==trialNposterror,'value'] <- mean(dat[trialIdx,acccolname], na.rm=TRUE)
  }
  #
  return(rbind(rtdf, rtabsdf, accdf, accabsdf))
}

getPESDistributions <- function(dat, pp, pp2, averageOnly, plotpp) {
  if(!averageOnly) {
    postpreData = lapply(unique(dat$subjects), function(x) {
      out <- getPostPreDifference(dat[dat$subjects==x,])
      out$subjects <- x
      out})
    postpreData <- do.call(rbind, postpreData)

    if(plotpp) {
      postprePP <-  data.frame(expand.grid(postn=unique(pp$postn), subjects=unique(pp$subjects)))
      diffs = parallel::mclapply(1:nrow(postprePP),
                                 function(x) {
                                   subject <- postprePP[x,'subjects']
                                   postn <- postprePP[x,'postn']
                                   out <- getPostPreDifference(pp[pp$subjects==subject&pp$postn==postn,])
                                   out$subjects <- subject
                                   out$postn <- postn
                                   out
                                 }, mc.cores=15)
      diffs <- do.call(rbind, diffs)
      postprePPall <- diffs
      postprePPsummary <- aggregate(value~subjects*measure*probs*trialNposterror,
                                    diffs,
                                    quantile, c(.025, .5, .975))

      if(!is.null(pp2)) {
        postprePP <-  data.frame(expand.grid(postn=unique(pp2$postn), subjects=unique(pp2$subjects)))
        diffs = parallel::mclapply(1:nrow(postprePP),
                                   function(x) {
                                     subject <- postprePP[x,'subjects']
                                     postn <- postprePP[x,'postn']
                                     out <- getPostPreDifference(pp2[pp2$subjects==subject&pp2$postn==postn,])
                                     out$subjects <- subject
                                     out$postn <- postn
                                     out
                                   }, mc.cores=15)
        diffs <- do.call(rbind, diffs)
        postprePPall2 <- diffs
        postprePPsummary2 <- aggregate(value~subjects*measure*probs*trialNposterror,
                                       diffs,
                                       quantile, c(.025, .5, .975))
      }
    }
  }

  ## across subs
  #  if(averageOnly) {
  postpreDataAverage <- getPostPreDifference(dat)
  postpreDataAverage$subjects <- 'Average'

  diffs <- lapply(unique(pp$postn), function(x) {
    out <- getPostPreDifference(pp[pp$postn==x,])
    out$postn <- x
    out
  })
  diffs <- do.call(rbind, diffs)
  diffs$subjects <- 'Average'
  postprePPAverage <- aggregate(value~subjects*measure*probs*trialNposterror,
                                diffs,
                                quantile, c(.025, .5, .975))

  if(!is.null(pp2)) {
    diffs <- lapply(unique(pp2$postn), function(x) {
      out <- getPostPreDifference(pp2[pp2$postn==x,])
      out$postn <- x
      out
    })
    diffs <- do.call(rbind, diffs)
    diffs$subjects <- 'Average'
    postprePPAverage2 <- aggregate(value~subjects*measure*probs*trialNposterror,
                                   diffs,
                                   quantile, c(.025, .5, .975))
  }
  #  }

  if(is.null(pp2)) {
    postprePPAverage2 <- NULL
    postprePPsummary2 <- NULL
    postprePPall2 <- NULL
  }
  if(averageOnly) {
    postpreData <- NULL
    postprePPsummary <- NULL
    postprePPsummary2 <- NULL
    postprePPall <- postprePPall2 <- NULL
  }
  # else {
  #   postpreDataAverage <- NULL
  #   postprePPAverage <- NULL
  #   postprePPAverage2 <- NULL
  # }
  return(list(postpreData=postpreData,
              postpreDataAverage=postpreDataAverage,
              postprePPAverage=postprePPAverage,
              postprePPAverage2=postprePPAverage2,
              postprePPsummary=postprePPsummary,
              postprePPsummary2=postprePPsummary2,
              postprePPall=postprePPall,
              postprePPall2=postprePPall2))
}

plotPES <- function(distributions=NULL, dat=NULL, pp=NULL, pp2=NULL,
                    orderSubjects=TRUE, firstTrialOnly=FALSE, use_mean=TRUE,
                    averageOnly=TRUE, meanOnly=FALSE,
                    xlab1='Trial relative to error', xlab2='Trial relative to error',
                    rtmeasure='rt', accmeasure='accuracy', # or rtabs / accabs
                    main1='', main2='', ylim1=NULL, ylim2=NULL, setPar=TRUE,
                    plotAxisTickLabels=TRUE,
                    ylab1='', ylab2='', plotAccuracy=TRUE, add.legend=FALSE, plotpp=TRUE) {
  dat <- add_trial_nrs(dat)

  if(is.null(distributions)) {
    distributions <- getDistributionsToPlot(dat, pp, pp2, averageOnly, plotpp)
  }
  postprePPAverage2 <- distributions$postprePPAverage2
  postprePPAverage <- distributions$postprePPAverage
  postpreDataAverage <- distributions$postpreDataAverage
  postpreData <- distributions$postpreData
  postprePPsummary2 <- distributions$postprePPsummary2
  postprePPsummary <- distributions$postprePPsummary

  if(meanOnly) {
    postprePPAverage <- postprePPAverage[postprePPAverage$probs %in% c('mean', '0.5'),]
    postprePPAverage2 <- postprePPAverage2[postprePPAverage2$probs %in% c('mean', '0.5'),]
    postpreData <- postpreData[postpreData$probs %in% c('mean', '0.5'),]
    postprePPsummary <- postprePPsummary[postprePPsummary$probs %in% c('mean', '0.5'),]
    postprePPsummary2 <- postprePPsummary2[postprePPsummary2$probs %in% c('mean', '0.5'),]
  }

  ## Plot
  if(averageOnly) {
    if(setPar){
      mar4 <- ifelse(main1=='Dataset 4', .3, .1)
      par(mar=c(2,3,2,mar4), bty='l')
    }
    if(is.null(ylim1)) ylim1 <- range(c(postprePPAverage[postprePPAverage$measure==rtmeasure,'value'],
                                        postpreDataAverage[postpreDataAverage$measure==rtmeasure,'value']))
    xs <- unique(postprePPAverage[postprePPAverage$measure==rtmeasure, 'trialNposterror'])

    plot(xs, xs, type='n', xlab=xlab1, ylab=ylab1,  xaxt='n', ylim=ylim1, main=main1, xlim=range(xs)+c(-.25, .25))
    grid()
    abline(h=0, lty=2, col='grey')
    axis(side=1, at=xs)
    ## data
    points(postpreDataAverage[postpreDataAverage$measure==rtmeasure,'trialNposterror'],
           postpreDataAverage[postpreDataAverage$measure==rtmeasure,'value'], pch=4, lwd=2)

    ## posterior predictives
    arrows(x0=postprePPAverage[postprePPAverage$measure==rtmeasure,'trialNposterror']-.15*(!is.null(pp2)),
           y0=postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,1],
           y1=postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,3],
           angle=90, code=3, length=.02, lwd=2, col=2)
    if(!is.null(pp2)) {
      arrows(x0=postprePPAverage2[postprePPAverage2$measure==rtmeasure,'trialNposterror']+.15,
             y0=postprePPAverage2[postprePPAverage2$measure==rtmeasure,'value'][,1],
             y1=postprePPAverage2[postprePPAverage2$measure==rtmeasure,'value'][,3],
             angle=90, code=3, length=.02, lwd=2, col=4)
      if(add.legend) legend('topright', c('Urgency', 'Threshold', 'Data'), lwd=c(1.5,1.5,1.5),
                            lty=c(1,1,NA), pch=c(NA,NA,4), col=c(2,4,1), bty='n', cex=.75)
    }

    ## Accuracy
    if(plotAccuracy) {
      if(setPar){
        par(mar=c(4,3,1,mar4))
      }
      if(is.null(ylim2)) {
        tmp1 <- postprePPAverage[postprePPAverage$measure==accmeasure,'value']
        tmp2 <- postpreDataAverage[postpreDataAverage$measure==accmeasure,'value']
        ylim2 <- range(c(tmp1[tmp1>0], tmp2[tmp2>0]))
      }
      xs <- unique(postprePPAverage[postprePPAverage$measure==accmeasure, 'trialNposterror'])
      plot(xs, xs, type='n', xlab=xlab2, ylab=ylab2, xaxt='n',xlim=range(xs)+c(-.25, .25), ylim=ylim2)
      abline(h=0, lty=2, col='grey')
      axis(side=1, at=xs)

      ## data
      points(postpreDataAverage[postpreDataAverage$measure==accmeasure,'trialNposterror'],
             postpreDataAverage[postpreDataAverage$measure==accmeasure,'value'], pch=4, lwd=2)

      ## PP
      arrows(x0=postprePPAverage[postprePPAverage$measure==accmeasure,'trialNposterror']-.15*(!is.null(pp2)),
             y0=postprePPAverage[postprePPAverage$measure==accmeasure,'value'][,1],
             y1=postprePPAverage[postprePPAverage$measure==accmeasure,'value'][,3],
             angle=90, code=3, length=.02, lwd=2, col=2)

      ##
      if(!is.null(pp2)) {
        arrows(x0=postprePPAverage2[postprePPAverage2$measure==accmeasure,'trialNposterror']+.15,
               y0=postprePPAverage2[postprePPAverage2$measure==accmeasure,'value'][,1],
               y1=postprePPAverage2[postprePPAverage2$measure==accmeasure,'value'][,3],
               angle=90, code=3, length=.02, lwd=2, col=4
        )
      }
    }
  } else {
    ## plot all subjects
    idx <- 1:nrow(postpreData)
    if(setPar){
      par(mar=c(2,3,2,1), bty='l')
    }
    if(is.null(ylim1)) ylim1 <- range(c(postprePPsummary[postprePPsummary$measure=='rt','value'],
                                        postpreData[postpreData$measure=='rt','value']))

    if(orderSubjects) {
      postpreData$subjects <- as.numeric(postpreData$subjects)
      postprePPsummary$subjects <- as.numeric(postprePPsummary$subjects)
      tmp <- postpreData[postpreData$probs=='0.5'&postpreData$measure=='rt'&postpreData$trialNposterror==1,]
      mapping <- data.frame(plotOrder=1:nrow(tmp), subjects=as.numeric(tmp[order(tmp$value),'subjects']))
      postpreData <- merge(postpreData, mapping)
      postprePPsummary <- merge(postprePPsummary, mapping)
      postpreData$subjects <- postpreData$plotOrder
      postprePPsummary$subjects <- postprePPsummary$plotOrder
      if(!is.null(pp2)) {
        postprePPsummary2$subjects <- as.numeric(postprePPsummary2$subjects)
        postprePPsummary2 <- merge(postprePPsummary2, mapping)
        postprePPsummary2$subjects <- postprePPsummary2$plotOrder
      }
    }
    xs <- as.numeric(unique(postpreData$subjects))
    plot(xs, xs, type='n', xlab=xlab1, ylab=ylab1,  xaxt='n', ylim=ylim1, main=main1, xlim=range(xs)+c(-.25, +.25))
    abline(h=0, lty=2, col='grey')
    if(plotAxisTickLabels) {
      axis(side=1, at=xs, labels=unique(postpreData$subjects))
    } #else {
      #axis(side=1, at=xs, labels=rep('', length(xs)))
    #}

    ## data
    points(as.numeric(postpreData[postpreData$measure=='rt'&postpreData$trialNposterror==1,'subjects']),
           postpreData[postpreData$measure=='rt'&postpreData$trialNposterror==1,'value'], pch=4, lwd=2)

    ## PP
    arrows(x0=as.numeric(postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'subjects'])-.15*(!is.null(pp2)),
           y0=postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'value'][,1],
           y1=postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'value'][,3],
           angle=90, code=3, length=.02, lwd=2, col=2)
    xs <- as.numeric(postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'subjects'])
    ys <- postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'value'][,2]
    lines(xs[order(xs)], ys[order(xs)], lty=1,col=2)

    if(!is.null(pp2)) {
      arrows(x0=as.numeric(postprePPsummary2[postprePPsummary2$measure=='rt'&postprePPsummary2$trialNposterror==1,'subjects'])+.15,
             y0=postprePPsummary2[postprePPsummary2$measure=='rt'&postprePPsummary2$trialNposterror==1,'value'][,1],
             y1=postprePPsummary2[postprePPsummary2$measure=='rt'&postprePPsummary2$trialNposterror==1,'value'][,3],
             angle=90, code=3, length=.02, lwd=2, col=4)
      if(add.legend) legend('topright', c('Urgency', 'Threshold', 'Data'), lwd=c(1.5,1.5,1.5),
                            lty=c(1,1,NA), pch=c(NA,NA,4), col=c(2,4,1), bty='n', cex=.75)
    }
  }
}


plotPES2 <- function(distributions=NULL, dat=NULL, pp=NULL, pp2=NULL,
                    orderSubjects=TRUE, firstTrialOnly=FALSE, use_mean=TRUE,
                    meanOnly=FALSE,
                    xlab1='Trial relative to error', xlab2='Trial relative to error',
                    rtmeasure='rt', accmeasure='accuracy', # or rtabs / accabs
                    main1='', main2='', ylim1=NULL, ylim2=NULL, setPar=TRUE,
                    plotAxisTickLabels=TRUE,
                    ylab1='', ylab2='', plotAccuracy=TRUE, add.legend=FALSE, plotpp=TRUE) {
  dat <- add_trial_nrs(dat)

#  meanrt <-

  if(is.null(distributions)) {
    distributions <- getDistributionsToPlot(dat, pp, pp2, averageOnly, plotpp)
  }
  postprePPAverage <- distributions$postprePPAverage
  postpreDataAverage <- distributions$postpreDataAverage
  postpreData <- distributions$postpreData
  postprePPsummary <- distributions$postprePPsummary

  if(meanOnly) {
    postprePPAverage <- postprePPAverage[postprePPAverage$probs %in% c('mean', '0.5'),]
    postpreData <- postpreData[postpreData$probs %in% c('mean', '0.5'),]
    postpreDataAverage <- postpreDataAverage[postpreDataAverage$probs %in% c('mean', '0.5'),]
    postprePPsummary <- postprePPsummary[postprePPsummary$probs %in% c('mean', '0.5'),]
  }

  ## Plot
  if(setPar){
    mar4 <- ifelse(main1=='Dataset 4', .3, .1)
    par(mar=c(2,3,2,mar4), bty='l')
  }
  if(is.null(ylim1)) ylim1 <- range(c(postprePPAverage[postprePPAverage$measure==rtmeasure,'value'],
                                      postpreDataAverage[postpreDataAverage$measure==rtmeasure,'value']))
  ylim1 <- ylim1*c(.95, 1.05)
  xs <- unique(postprePPAverage[postprePPAverage$measure==rtmeasure, 'trialNposterror'])

  plot(xs, xs, type='n', xlab=xlab1, ylab=ylab1,  xaxt='n', ylim=ylim1, main=main1, xlim=range(xs)+c(-.25, .25))
#  grid()
#  abline(h=0, lty=2, col='grey')
  axis(side=1, at=xs)
  abline(h=mean(aggregate(rt~subjects,dat,mean)[,2]), col=par()$col, lty=2)
  abline(v=0, lty=2, col=par()$col)

  ## data
  points(postpreDataAverage[postpreDataAverage$measure==rtmeasure,'trialNposterror'],
         postpreDataAverage[postpreDataAverage$measure==rtmeasure,'value'], pch=4, lwd=2)

  ## posterior predictives
  arrows(x0=postprePPAverage[postprePPAverage$measure==rtmeasure,'trialNposterror']-.15*(!is.null(pp2)),
         y0=postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,1],
         y1=postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,3],
         angle=90, code=3, length=.02, lwd=2, col=2)
  xs <- postprePPAverage[postprePPAverage$measure==rtmeasure,'trialNposterror']
  lines(xs, postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,2], col=2)
}

getPostPreDifferenceBySATCondition <- function(dat, rtcolname='rt', acccolname='accuracy', conditionColumn='E') {
  errors <- dat$accuracy == 0
  ## remove first trial (per participant) and last trial (per participant)
  if(!'trials' %in% colnames(dat)) dat <- add_trial_nrs(dat)
  trialRangeBySub <- aggregate(trials~subjects, dat, range)

  dat$exclude_trial <- FALSE
  for(subject in unique(trialRangeBySub$subjects)) {
    ## don't include errors on the first and last trial; can't calculate post-error slowing because of lack of a pre / post error trial
    dat[dat$subjects==subject,'exclude_trial'] = dat[dat$subjects==subject,'trials'] %in% trialRangeBySub[trialRangeBySub$subjects==subject,'trials']
  }
  errors <- errors & !dat$exclude_trial

  ## Here we split
  conditions <- unique(dat[,conditionColumn])
  # is_neutral = dat[,conditionColumn] == 'neutral'
  # is_speed = dat[,conditionColumn] == 'speed'
  # is_accuracy = dat[,conditionColumn] == 'accuracy'
  #
  # errors_neutral = which(errors & is_neutral)
  # errors_speed = which(errors & is_speed)
  # errors_accuracy = which(errors & is_accuracy)

  #  errors <- which(errors)
  #  preError <- errors-1

  f <- mean #function(x, ...) quantile(x, probs=c(.1, .5, .9), ...)
  rtabsdf <- expand.grid(condition=conditions, trialNposterror=seq(-3,7), probs='mean', measure='rtabs', value=NA)
  for(condition in conditions) {
    errors_ = which(errors & dat[,conditionColumn] == condition)
    for(trialNposterror in seq(-3,7)) {
      trialIdx <- errors_+trialNposterror
      trialIdx <- trialIdx[trialIdx>0]  ## remove
      rtabsdf[rtabsdf$condition==condition&rtabsdf$trialNposterror==trialNposterror,'value'] <- f(dat[trialIdx,rtcolname], na.rm=TRUE)
    }
  }
  #
  return(rtabsdf)
}

getPESDistributionsBySAT <- function(dat, pp, pp2, averageOnly, plotpp) {
  if(!averageOnly) {
    postpreData = lapply(unique(dat$subjects), function(x) {
      out <- getPostPreDifferenceBySATCondition(dat[dat$subjects==x,])
      out$subjects <- x
      out})
    postpreData <- do.call(rbind, postpreData)

    if(plotpp) {
      postprePP <-  data.frame(expand.grid(postn=unique(pp$postn), subjects=unique(pp$subjects)))
      diffs = parallel::mclapply(1:nrow(postprePP),
                                 function(x) {
                                   subject <- postprePP[x,'subjects']
                                   postn <- postprePP[x,'postn']
                                   out <- getPostPreDifferenceBySATCondition(pp[pp$subjects==subject&pp$postn==postn,])
                                   out$subjects <- subject
                                   out$postn <- postn
                                   out
                                 }, mc.cores=15)
      diffs <- do.call(rbind, diffs)
      postprePPall <- diffs
      postprePPsummary <- aggregate(value~subjects*measure*probs*trialNposterror,
                                    diffs,
                                    quantile, c(.025, .5, .975))

      if(!is.null(pp2)) {
        postprePP <-  data.frame(expand.grid(postn=unique(pp2$postn), subjects=unique(pp2$subjects)))
        diffs = parallel::mclapply(1:nrow(postprePP),
                                   function(x) {
                                     subject <- postprePP[x,'subjects']
                                     postn <- postprePP[x,'postn']
                                     out <- getPostPreDifferenceBySATCondition(pp2[pp2$subjects==subject&pp2$postn==postn,])
                                     out$subjects <- subject
                                     out$postn <- postn
                                     out
                                   }, mc.cores=15)
        diffs <- do.call(rbind, diffs)
        postprePPall2 <- diffs
        postprePPsummary2 <- aggregate(value~subjects*measure*probs*trialNposterror,
                                       diffs,
                                       quantile, c(.025, .5, .975))
      }
    }
  }

  ## across subs
  #  if(averageOnly) {
  postpreDataAverage <- getPostPreDifferenceBySATCondition(dat)
  postpreDataAverage$subjects <- 'Average'

  diffs <- lapply(unique(pp$postn), function(x) {
    out <- getPostPreDifferenceBySATCondition(pp[pp$postn==x,])
    out$postn <- x
    out
  })
  diffs <- do.call(rbind, diffs)
  diffs$subjects <- 'Average'
  postprePPAverage <- aggregate(value~subjects*measure*probs*trialNposterror,
                                diffs,
                                quantile, c(.025, .5, .975))

  if(!is.null(pp2)) {
    diffs <- lapply(unique(pp2$postn), function(x) {
      out <- getPostPreDifferenceBySATCondition(pp2[pp2$postn==x,])
      out$postn <- x
      out
    })
    diffs <- do.call(rbind, diffs)
    diffs$subjects <- 'Average'
    postprePPAverage2 <- aggregate(value~subjects*measure*probs*trialNposterror,
                                   diffs,
                                   quantile, c(.025, .5, .975))
  }
  #  }

  if(is.null(pp2)) {
    postprePPAverage2 <- NULL
    postprePPsummary2 <- NULL
    postprePPall2 <- NULL
  }
  if(averageOnly) {
    postpreData <- NULL
    postprePPsummary <- NULL
    postprePPsummary2 <- NULL
    postprePPall <- postprePPall2 <- NULL
  }
  # else {
  #   postpreDataAverage <- NULL
  #   postprePPAverage <- NULL
  #   postprePPAverage2 <- NULL
  # }
  return(list(postpreData=postpreData,
              postpreDataAverage=postpreDataAverage,
              postprePPAverage=postprePPAverage,
              postprePPAverage2=postprePPAverage2,
              postprePPsummary=postprePPsummary,
              postprePPsummary2=postprePPsummary2,
              postprePPall=postprePPall,
              postprePPall2=postprePPall2))
}



# MS3-RDM -----------------------------------------------------------------
get_PES_measures <- function(model) {
  wagenmakers2004_CS <- loadDataAndPP('wagenmakers2004_CS', model, 'NULL')
  forstmann2008 <- loadDataAndPP('forstmann2008', model, 'NULL')
  mileticvanmaanen2019exp2block2 <- loadDataAndPP('mileticvanmaanen2019exp2block2', model, 'NULL')
  # filter these to get rid of subjects with very few errors
  dat <- mileticvanmaanen2019exp2block2$dat
  dat$error <- !dat$accuracy
  nerrors <- aggregate(error~subjects,dat,sum)
  goodsubs <- nerrors[nerrors[,2]>40,1]
  dat <- dat[dat$subjects %in% goodsubs,]
  mileticvanmaanen2019exp2block2$dat <- mileticvanmaanen2019exp2block2$dat[mileticvanmaanen2019exp2block2$dat$subjects %in% goodsubs,]
  mileticvanmaanen2019exp2block2$pp <- mileticvanmaanen2019exp2block2$pp[mileticvanmaanen2019exp2block2$pp$subjects %in% goodsubs,]
  wagenmakers2008exp2 <- loadDataAndPP('wagenmakers2008exp2', model, 'NULL')

  # Estimate PES distributions
  distributions1 <- getPESDistributions(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, pp2=NULL, averageOnly = FALSE, plotpp=TRUE)
  distributions2 <- getPESDistributions(dat=forstmann2008$dat, pp=forstmann2008$pp, pp2=NULL, averageOnly = FALSE, plotpp=TRUE)
  distributions3 <- getPESDistributions(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, pp2=NULL, averageOnly = FALSE, plotpp=TRUE)
  distributions4 <- getPESDistributions(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, pp2=NULL, averageOnly = FALSE, plotpp=TRUE)

  # Aggregate
  distributions1$postpreDataAverage <- aggregate(value~measure*probs*trialNposterror, distributions1$postpreData, mean)
  distributions1$postprePPAverage <- aggregate(value~measure*probs*trialNposterror, aggregate(value~measure*probs*trialNposterror*postn, distributions1$postprePPall, mean), quantile, c(0.025, 0.5, 0.975))

  distributions2$postpreDataAverage <- aggregate(value~measure*probs*trialNposterror, distributions2$postpreData, mean)
  distributions2$postprePPAverage <- aggregate(value~measure*probs*trialNposterror, aggregate(value~measure*probs*trialNposterror*postn, distributions2$postprePPall, mean), quantile, c(0.025, 0.5, 0.975))

  distributions3$postpreDataAverage <- aggregate(value~measure*probs*trialNposterror, distributions3$postpreData, mean)
  distributions3$postprePPAverage <- aggregate(value~measure*probs*trialNposterror, aggregate(value~measure*probs*trialNposterror*postn, distributions3$postprePPall, mean), quantile, c(0.025, 0.5, 0.975))

  distributions4$postpreDataAverage <- aggregate(value~measure*probs*trialNposterror, distributions4$postpreData, mean)
  distributions4$postprePPAverage <- aggregate(value~measure*probs*trialNposterror, aggregate(value~measure*probs*trialNposterror*postn, distributions4$postprePPall, mean), quantile, c(0.025, 0.5, 0.975))


  # For dataset 2, get aggregates of SAT separate
  # collapse across accuracy/neutral
  forstmann2008$dat[forstmann2008$dat$E=='neutral','E'] <- 'accuracy'  # collapse ACC and NEU
  forstmann2008$dat <- droplevels(forstmann2008$dat)
  forstmann2008$pp[forstmann2008$pp$E=='neutral','E'] <- 'accuracy'
  forstmann2008$pp <- droplevels(forstmann2008$pp)

  distributions2b <- getPESDistributionsBySAT(dat=forstmann2008$dat, pp=forstmann2008$pp, pp2=NULL, averageOnly = FALSE, plotpp=TRUE)
  distributions2b$postpreDataAverage <- aggregate(value~condition*measure*probs*trialNposterror, distributions2b$postpreData, mean)
  distributions2b$postprePPAverage <- aggregate(value~condition*measure*probs*trialNposterror, aggregate(value~condition*measure*probs*trialNposterror*postn, distributions2b$postprePPall, mean), quantile, c(0.025, 0.5, 0.975))

  # split for plotting
  distributions2bSPD <- distributions2bACC <- distributions2b
  distributions2bSPD$postpreDataAverage <- distributions2bSPD$postpreDataAverage[distributions2bSPD$postpreDataAverage$condition=='speed',]
  distributions2bACC$postpreDataAverage <- distributions2bACC$postpreDataAverage[distributions2bACC$postpreDataAverage$condition=='accuracy',]
  distributions2bSPD$postprePPAverage <- distributions2bSPD$postprePPAverage[distributions2bSPD$postprePPAverage$condition=='speed',]
  distributions2bACC$postprePPAverage <- distributions2bACC$postprePPAverage[distributions2bACC$postprePPAverage$condition=='accuracy',]
  return(list(wagenmakers2004_CS=wagenmakers2004_CS,
              forstmann2008=forstmann2008,
              mileticvanmaanen2019exp2block2=mileticvanmaanen2019exp2block2,
              wagenmakers2008exp2=wagenmakers2008exp2,

              distributions1=distributions1,
              distributions2=distributions2,
              distributions3=distributions3,
              distributions4=distributions4,

              distributions2b=distributions2b,
              distributions2bSPD=distributions2bSPD,
              distributions2bACC=distributions2bACC
              ))
}
model1 <- 'zSMuAHbV'
MS3 <- get_PES_measures(model1)
model2 <- 'zSMuAH'
MS2 <- get_PES_measures(model2)


# Plot
# Make plot combining group-level and individual-level effects
## New plot: Top row = group, bottom row = individual
m <- matrix(nrow=100, ncol=100)
row1 <- 1:50
row2 <- row1+50
column1Upper <- 1:20
column2Upper <- 21:40
column3Upper <- 41:60
column4Upper <- 61:80
column5Upper <- 81:100
## Upper row
m[row1,column1Upper] <- 1
m[row1,column2Upper] <- 2
m[row1,column3Upper] <- 3
m[row1,column4Upper] <- 4
m[row1,column5Upper] <- 5
## Lower row
column1Lower <- 1:14
column2Lower <- 15:36
column3Lower <- 37:82
column4Lower <- 83:100
m[row2,column1Lower] <- 6
m[row2,column2Lower] <- 7
m[row2,column3Lower] <- 8
m[row2,column4Lower] <- 9

for(model in c(model1, model2)) {
  if(model == model1) list2env(MS3,globalenv())
  if(model == model2) list2env(MS2,globalenv())

  for(ftype in c('pdf','jpeg')) {
    if(ftype == 'pdf') pdf(file=file.path(figures_dir, paste0('accuracy_history_group_and_individual_', model, '.pdf')), width=8.5, height=5)
    if(ftype == 'jpeg') jpeg(file=file.path(figures_dir, paste0('accuracy_history_group_and_individual_', model, '.jpeg')), width=8.5, height=5, units='in', res=500)
      l <- layout(m)
      #layout.show(l)
      par(mgp=c(2,1,0), bty='l', mar=c(3,2,2,.5), oma=c(0,1,0,0))
      plotPES2(distributions=distributions1, dat=wagenmakers2004_CS$dat, pp2=NULL, main1='Dataset 1', rtmeasure = 'rtabs', accmeasure='accabs', ylab1='RT (s)', ylab='Accuracy', plotAccuracy = FALSE, meanOnly=TRUE, setPar = FALSE)
      mtext('A. RT (s)', side=2, cex=.66, line=2)
      plotPES2(distributions=distributions2bSPD, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (spd)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
      plotPES2(distributions=distributions2bACC, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (acc/neu)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
      plotPES2(distributions=distributions3, dat=mileticvanmaanen2019exp2block2$dat, pp2=NULL, main1='Dataset 3', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
      plotPES2(distributions=distributions4, dat=wagenmakers2008exp2$dat, pp2=NULL, main1='Dataset 4', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)

      #layout(matrix(1:4, ncol=4), widths=c(.15, .25, .5, .2))
      par(mar=c(2,2,3,.5), mgp=c(2,.75,0), bty='l') # , oma=c(0,1,0,0)
      plotPES(distributions=distributions1, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='', main1='Dataset 1')
      mtext('Participant', side=1, line=1, cex=.66)
      abline(h=0, col=adjustcolor(1, alpha.f=.5))
      mtext('B. RT difference (s)', side=2, cex=.66, line=2)
      plotPES(distributions=distributions2, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 2')
      mtext('Participant', side=1, line=1, cex=.66)
      abline(h=0, col=adjustcolor(1, alpha.f=.5))
      plotPES(distributions=distributions3, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 3')
      mtext('Participant', side=1, line=1, cex=.66)
      abline(h=0, col=adjustcolor(1, alpha.f=.5))
      plotPES(distributions=distributions4, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 4')
      mtext('Participant', side=1, line=1, cex=.66)
      abline(h=0, col=adjustcolor(1, alpha.f=.5))
      dev.off()
  }


## Group-level only
for(colors_ in c('black', 'white')) {
  if(colors_ == 'black') {
    pdf(file=file.path(figures_dir, paste0('accuracy_history_group_', model, '_black.pdf')), width=8.5, height=2.75)
    if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white')
  } else {
    pdf(file=file.path(figures_dir, paste0('accuracy_history_group_', model, '.pdf')), width=8.5, height=2.75)
  }
  par(mgp=c(2,1,0), mfrow=c(1,5), bty='l', mar=c(3,3,2,.1))
  plotPES2(distributions=distributions1, dat=wagenmakers2004_CS$dat, pp2=NULL, main1='Dataset 1', rtmeasure = 'rtabs', accmeasure='accabs', ylab1='RT (s)', ylab='Accuracy', plotAccuracy = FALSE, meanOnly=TRUE, setPar = FALSE)
  plotPES2(distributions=distributions2bSPD, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (spd)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions2bACC, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (acc/neu)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions3, dat=mileticvanmaanen2019exp2block2$dat, pp2=NULL, main1='Dataset 3', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  plotPES2(distributions=distributions4, dat=wagenmakers2008exp2$dat, pp2=NULL, main1='Dataset 4', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  dev.off()
}

# Individual only
for(colors_ in c('black', 'white')) {
  if(colors_ == 'black') {
    pdf(file=file.path(figures_dir, paste0('accuracy_history_individual_', model, '_black.pdf')), width=8, height=2.75)
    if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white')
  } else {
    pdf(file=file.path(figures_dir, paste0('accuracy_history_individual_', model, '.pdf')), width=8, height=2.5)
  }
  layout(matrix(1:4, ncol=4), widths=c(.15, .25, .5, .2))
  par(mar=c(2,2,3,.5), mgp=c(2,.75,0), bty='l', oma=c(0,1,0,0))
  plotPES(distributions=distributions1, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='', main1='Dataset 1')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  mtext('RT difference (s)', side=2, cex=.66, line=2)
  plotPES(distributions=distributions2, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 2')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions3, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 3')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions4, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 4')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  dev.off()
}
}


## Combine MS2-RDM group, individual with MS3-RDM group
m <- matrix(nrow=100, ncol=100)
row1 <- 1:33
row2 <- 34:67
row3 <- 67:100
column1Upper <- 1:20
column2Upper <- 21:40
column3Upper <- 41:60
column4Upper <- 61:80
column5Upper <- 81:100
## Upper row
m[row1,column1Upper] <- 1
m[row1,column2Upper] <- 2
m[row1,column3Upper] <- 3
m[row1,column4Upper] <- 4
m[row1,column5Upper] <- 5
## Middle row
column1Lower <- 1:14
column2Lower <- 15:36
column3Lower <- 37:82
column4Lower <- 83:100
m[row2,column1Lower] <- 6
m[row2,column2Lower] <- 7
m[row2,column3Lower] <- 8
m[row2,column4Lower] <- 9
# Bottom row
m[row3,column1Upper] <- 10
m[row3,column2Upper] <- 11
m[row3,column3Upper] <- 12
m[row3,column4Upper] <- 13
m[row3,column5Upper] <- 14


## Group 1
for(ftype in c('jpeg', 'pdf')) {
  if(ftype == 'pdf') pdf(file=file.path(figures_dir, paste0('accuracy_history_group_and_individual_zSMuAH_and_zSMuAHbV.pdf')), width=8.5, height=7)
  if(ftype == 'jpeg') jpeg(file=file.path(figures_dir, paste0('accuracy_history_group_and_individual_zSMuAH_and_zSMuAHbV.jpeg')), width=8.5, height=7, units='in', res=500, quality=100)
  list2env(MS2,globalenv())
  l <- layout(m)
  par(mgp=c(2,1,0), bty='l', mar=c(3,2,2,.5), oma=c(0,2,0,0))
  plotPES2(distributions=distributions1, dat=wagenmakers2004_CS$dat, pp2=NULL, main1='Dataset 1', rtmeasure = 'rtabs', accmeasure='accabs', ylab1='RT (s)', ylab='Accuracy', plotAccuracy = FALSE, meanOnly=TRUE, setPar = FALSE)
  mtext('RT (s)', side=2, cex=.66, line=2)
  mtext('A. MS2-RDM', side=2, cex=par()$cex*par()$cex.main, line=3, font=2)
  plotPES2(distributions=distributions2bSPD, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (spd)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions2bACC, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (acc/neu)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions3, dat=mileticvanmaanen2019exp2block2$dat, pp2=NULL, main1='Dataset 3', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  plotPES2(distributions=distributions4, dat=wagenmakers2008exp2$dat, pp2=NULL, main1='Dataset 4', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)

  ## individual 1
  par(mar=c(2,2,3,.5), mgp=c(2,.75,0), bty='l') # , oma=c(0,1,0,0)
  plotPES(distributions=distributions1, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='', main1='Dataset 1')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  mtext('RT difference (s)', side=2, cex=.66, line=2)
  mtext('B. MS2-RDM', side=2, cex=par()$cex*par()$cex.main, line=3, font=2)
  plotPES(distributions=distributions2, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 2')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions3, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 3')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions4, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 4')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))

  ## group 2
  list2env(MS3,globalenv())
  par(mgp=c(2,1,0), bty='l', mar=c(3,2,2,.5))
  plotPES2(distributions=distributions1, dat=wagenmakers2004_CS$dat, pp2=NULL, main1='Dataset 1', rtmeasure = 'rtabs', accmeasure='accabs', ylab1='RT (s)', ylab='Accuracy', plotAccuracy = FALSE, meanOnly=TRUE, setPar = FALSE)
  mtext('RT (s)', side=2, cex=.66, line=2)
  mtext('C. MS3-RDM', side=2, cex=par()$cex*par()$cex.main, line=3, font=2)
  plotPES2(distributions=distributions2bSPD, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (spd)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions2bACC, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (acc/neu)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions3, dat=mileticvanmaanen2019exp2block2$dat, pp2=NULL, main1='Dataset 3', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  plotPES2(distributions=distributions4, dat=wagenmakers2008exp2$dat, pp2=NULL, main1='Dataset 4', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  dev.off()
}



# group-individual tiwce --------------------------------------------------
m <- matrix(nrow=100, ncol=100)
row1 <- 1:25
row2 <- row1+25 #$34:67
row3 <- row2+25 # 67:100
row4 <- row3+25
column1Upper <- 1:20
column2Upper <- 21:40
column3Upper <- 41:60
column4Upper <- 61:80
column5Upper <- 81:100
## Upper row
m[row1,column1Upper] <- 1
m[row1,column2Upper] <- 2
m[row1,column3Upper] <- 3
m[row1,column4Upper] <- 4
m[row1,column5Upper] <- 5
## Middle row
column1Lower <- 1:14
column2Lower <- 15:36
column3Lower <- 37:82
column4Lower <- 83:100
m[row2,column1Lower] <- 6
m[row2,column2Lower] <- 7
m[row2,column3Lower] <- 8
m[row2,column4Lower] <- 9
# Bottom row
m[row3,column1Upper] <- 10
m[row3,column2Upper] <- 11
m[row3,column3Upper] <- 12
m[row3,column4Upper] <- 13
m[row3,column5Upper] <- 14
# Bottom bototm row
m[row4,column1Lower] <- 15
m[row4,column2Lower] <- 16
m[row4,column3Lower] <- 17
m[row4,column4Lower] <- 18


for(ftype in c('jpeg', 'pdf')) {
  if(ftype == 'pdf') pdf(file=file.path(figures_dir, paste0('accuracy_history_group_and_individual_zSMuAH_and_zSMuAHbV_full.pdf')), width=8.5, height=8.5)
  if(ftype == 'jpeg') jpeg(file=file.path(figures_dir, paste0('accuracy_history_group_and_individual_zSMuAH_and_zSMuAHbV_full.jpeg')), width=8.5, height=8.5, units='in', res=500, quality=100)

  #pdf(file=file.path(figures_dir, paste0('accuracy_history_group_and_individual_zSMuAH_and_zSMuAHbV_full.pdf')), width=8.5, height=8.5)
  list2env(MS2,globalenv())
  l <- layout(m)
  par(mgp=c(2,1,0), bty='l', mar=c(3,2,2,.5), oma=c(0,2,0,0))
  plotPES2(distributions=distributions1, dat=wagenmakers2004_CS$dat, pp2=NULL, main1='Dataset 1', rtmeasure = 'rtabs', accmeasure='accabs', ylab1='RT (s)', ylab='Accuracy', plotAccuracy = FALSE, meanOnly=TRUE, setPar = FALSE)
  mtext('RT (s)', side=2, cex=.66, line=2)
  mtext('A. MS2-RDM', side=2, cex=par()$cex*par()$cex.main, line=3, font=2)
  plotPES2(distributions=distributions2bSPD, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (spd)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions2bACC, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (acc/neu)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions3, dat=mileticvanmaanen2019exp2block2$dat, pp2=NULL, main1='Dataset 3', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  plotPES2(distributions=distributions4, dat=wagenmakers2008exp2$dat, pp2=NULL, main1='Dataset 4', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)

  ## individual 1
  par(mar=c(2,2,3,.5), mgp=c(2,.75,0), bty='l') # , oma=c(0,1,0,0)
  plotPES(distributions=distributions1, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='', main1='Dataset 1')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  mtext('RT difference (s)', side=2, cex=.66, line=2)
  mtext('B. MS2-RDM', side=2, cex=par()$cex*par()$cex.main, line=3, font=2)
  plotPES(distributions=distributions2, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 2')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions3, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 3')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions4, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 4')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))

  ## group 2
  list2env(MS3,globalenv())
  par(mgp=c(2,1,0), bty='l', mar=c(3,2,2,.5))
  plotPES2(distributions=distributions1, dat=wagenmakers2004_CS$dat, pp2=NULL, main1='Dataset 1', rtmeasure = 'rtabs', accmeasure='accabs', ylab1='RT (s)', ylab='Accuracy', plotAccuracy = FALSE, meanOnly=TRUE, setPar = FALSE)
  mtext('RT (s)', side=2, cex=.66, line=2)
  mtext('C. MS3-RDM', side=2, cex=par()$cex*par()$cex.main, line=3, font=2)
  plotPES2(distributions=distributions2bSPD, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (spd)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions2bACC, dat=forstmann2008$dat, pp2=NULL, main1='Dataset 2 (acc/neu)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
  plotPES2(distributions=distributions3, dat=mileticvanmaanen2019exp2block2$dat, pp2=NULL, main1='Dataset 3', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  plotPES2(distributions=distributions4, dat=wagenmakers2008exp2$dat, pp2=NULL, main1='Dataset 4', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)

  # individual 2
  par(mar=c(2,2,3,.5), mgp=c(2,.75,0), bty='l') # , oma=c(0,1,0,0)
  plotPES(distributions=distributions1, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='', main1='Dataset 1')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  mtext('RT difference (s)', side=2, cex=.66, line=2)
  mtext('D. MS3-RDM', side=2, cex=par()$cex*par()$cex.main, line=3, font=2)
  plotPES(distributions=distributions2, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 2')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions3, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 3')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  plotPES(distributions=distributions4, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='Dataset 4')
  mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  dev.off()
}

## Dataset 5 -----
mileticvanmaanen2019exp2block1 <- loadDataAndPP('mileticvanmaanen2019exp2block1', 'zSMuAHbV', 'NULL')
# filter these to get rid of subjects with very few errors
mileticvanmaanen2019exp2block1$pp$rt <- pmin(mileticvanmaanen2019exp2block1$pp$rt, 3) # Censor, like data
dat <- mileticvanmaanen2019exp2block1$dat
dat$error <- !dat$accuracy
nerrors <- aggregate(error~subjects,dat,sum)
goodsubs <- nerrors[nerrors[,2]>40,1]
dat <- dat[dat$subjects %in% goodsubs,]
mileticvanmaanen2019exp2block1$dat <- mileticvanmaanen2019exp2block1$dat[mileticvanmaanen2019exp2block1$dat$subjects %in% goodsubs,]
mileticvanmaanen2019exp2block1$pp <- mileticvanmaanen2019exp2block1$pp[mileticvanmaanen2019exp2block1$pp$subjects %in% goodsubs,]

distributions3b <- getPESDistributions(dat=mileticvanmaanen2019exp2block1$dat, pp=mileticvanmaanen2019exp2block1$pp, pp2=NULL, averageOnly = FALSE, plotpp=TRUE)
distributions3b$postpreDataAverage <- aggregate(value~measure*probs*trialNposterror, distributions3b$postpreData, mean)
distributions3b$postprePPAverage <- aggregate(value~measure*probs*trialNposterror, aggregate(value~measure*probs*trialNposterror*postn, distributions3b$postprePPall, mean), quantile, c(0.025, 0.5, 0.975))


for(ftype in c('jpeg','pdf')) {
  fn = file.path(figures_dir, 'accuracy_history_group_and_individual_dataset5_zSMuAHbV.pdf')
  if(ftype == 'pdf') pdf(file=fn, width=5, height=3)
  if(ftype == 'jpeg') jpeg(file=gsub('.pdf', '.jpeg', fn), width=5, height=3, units='in', res=500, quality=100)
  par(mfrow=c(1,2), bty='l', cex=.66, mar=c(4,3,2,1))
  plotPES2(distributions=distributions3b, dat=mileticvanmaanen2019exp2block1$dat, pp2=NULL, main1='A. Dataset 5 (group)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
  mtext('RT (s)', side=2, cex=.66, line=2)
  plotPES(distributions=distributions3b, pp2=NULL, averageOnly=FALSE, plotAccuracy = FALSE, setPar=FALSE, plotAxisTickLabels = FALSE, orderSubjects = TRUE, meanOnly=TRUE, xlab1='Participant', main1='B. Dataset 5 (individual)')
  #mtext('Participant', side=1, line=1, cex=.66)
  abline(h=0, col=adjustcolor(1, alpha.f=.5))
  mtext('RT difference (s)', side=2, cex=.66, line=2)
  dev.off()
}


# Post-error accuracy -----------------------------------------------------
get_errorrate1 <- function(x) {
  tmp <- rle(x$accuracy)
  n_consecutive_errors <- tmp$lengths[tmp$values==0]
  return(n_consecutive_errors)
}

get_rate_proportions <- function(rates, n, max_errors=99) {
  num <- sum(rates>n & rates<max_errors)
  den <- sum(rates>(n-1) & rates<max_errors)
  if(den < 5) return(NA)  ## calculate a proportion only if there's at least 10 occurances of the denominator, otherwise it hardly makes sense
  x <- num/den #sum(rates>n & rates<max_errors)/sum(rates>(n-1) & rates<max_errors)
  if(is.infinite(x) | is.na(x)) {
    return(NA)
  } else{
    return(x)
  }
}


get_stats <- function(dat) {
  rates <- get_errorrate1(dat)
  p_error <- mean(dat$accuracy==0)
  p_error_given_1error <- get_rate_proportions(rates, 1, max_errors=5)
  p_error_given_2error <- get_rate_proportions(rates, 2, max_errors=5)
  p_error_given_3error <- get_rate_proportions(rates, 3, max_errors=5)

  return(data.frame(p_error=p_error,
                    p_error_given_1error=p_error_given_1error,
                    p_error_given_2error=p_error_given_2error,
                    p_error_given_3error=p_error_given_3error))
}


tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
#pdf(file='./figures/post_error_accuracy.pdf', width=8, height=2.5)
par(mfrow=c(1,4))
for(task in tasks) {
  dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
  pp <- EMC2:::loadRData(paste0('../dynamicEAMs_tmp/samples/', task, '_model-RDM-zSMuAHbV-DCT-NULL_pp.RData'))

  lst <- lapply(unique(dat$subjects), function(x) get_stats(dat[dat$subjects==x,]))
  out <- do.call(rbind, lst)
  #out[is.na(out)] <- 0
  idx <- !is.na(out[,3])
  #out[,2:4] <- out[,2:4]-out[,1]  # get subjectwise difference first so we get the correct variance
  mns <- apply(out[idx,], 2, mean)
  #mns[2:4] <- mns[2:4]+mns[1]   # add back to means
  ses <- apply(out[idx,], 2, function(x) sd(x)/length(x))

  # aggregation is a pain here.....
  allmns <- mclapply(unique(pp$postn), function(postn) {
    tmp <- lapply(unique(dat$subjects), function(x) get_stats(pp[pp$subjects==x&pp$postn==postn,]))
    out2 <- do.call(rbind, tmp)
    #idx2 <- !is.na(out[,3])
    #out2[,2:4] <- out2[,2:4]-out2[,1]
    out2[is.na(out2)] <- 0  # NB: NA here means that there were 0 occurrences
    mns2 <- apply(out2[idx,], 2, mean)  # same set of subjects!
    mns2
    #out2
  }, mc.cores=20)
  allmns <- do.call(rbind, allmns)
  #allmns[,2:4] <- allmns[,2:4] + allmns[,1] # add back p(error)
  mns2 <- apply(allmns, 2, mean)
  qps <- apply(allmns, 2, quantile, c(.025, .5, .975), na.rm=TRUE)

  ## add back p(error)
  # mns2[2:4] <- mns2[2:4]+mns2[1]
  # qps[,2:3] <- qps[,2:3]+qps[,1]

  ylim = range(qps, mns, mns2, na.rm=TRUE)#*c(.9, 1.1)
  ylim = ylim*c(1-0.1*sign(ylim)[1], 1.1)

  mns <- mns[1:3]
  qps <- qps[,1:3]
  plot(1:length(mns), mns, ylim=ylim, xaxt='n', bty='l', xlab='', ylab='p', main=task)
  axis(side=1, at=1:4, labels=c('p(error)', 'p(error|1error)', 'p(error|2errors)', 'p(error|3errors)'))
  abline(h=0, lty=2)
  arrows(1:ncol(qps), y0=qps[1,], y1=qps[3,], code=3, angle=90, length=.05)
}
#dev.off()

