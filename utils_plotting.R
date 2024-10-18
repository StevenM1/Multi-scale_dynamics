##plotting functions
library(stats)
library(abind)
library(stringr)
library(Hmisc)
library(forecast)
library(parallel)

add_trial_nrs <- function(dat) {
  dat$trial_nr <- NA
  for(s in unique(dat$subjects)) {
    dat[dat$subjects==s, 'trials'] <- 1:(sum(dat$subjects==s))
  }
  return(dat)
}


# Sequential effects ------------------------------------------------------
plotSequentialEffects3 <- function(dat=dat, pp=pp, stimulusCoded = TRUE, nHistory=4, dataColumn='S', sub='', main1='', main2='') {
  dat <- add_trial_nrs(dat)
  dat$Sminus5 <- str_sub(Lag(dat[,dataColumn], 5),1,1)
  dat$Sminus4 <- str_sub(Lag(dat[,dataColumn], 4),1,1)
  dat$Sminus3 <- str_sub(Lag(dat[,dataColumn], 3),1,1)
  dat$Sminus2 <- str_sub(Lag(dat[,dataColumn], 2),1,1)
  dat$Sminus1 <- str_sub(Lag(dat[,dataColumn], 1),1,1)
  if(nHistory == 4) dat$pattern <- paste(dat$Sminus4,dat$Sminus3,dat$Sminus2,dat$Sminus1, sep='')
  if(nHistory == 3) dat$pattern <- paste(dat$Sminus3,dat$Sminus2,dat$Sminus1, sep='')
  if(nHistory == 2) dat$pattern <- paste(dat$Sminus2,dat$Sminus1, sep='')
  dat <- dat[!dat$trials<6,]

  pp <- pp[order(pp$subjects,pp$postn,pp$trials),]
  pp$Sminus5 <- str_sub(Lag(pp[,dataColumn], 5),1,1)
  pp$Sminus4 <- str_sub(Lag(pp[,dataColumn], 4),1,1)
  pp$Sminus3 <- str_sub(Lag(pp[,dataColumn], 3),1,1)
  pp$Sminus2 <- str_sub(Lag(pp[,dataColumn], 2),1,1)
  pp$Sminus1 <- str_sub(Lag(pp[,dataColumn], 1),1,1)
  if(nHistory == 4) pp$pattern <- paste(pp$Sminus4,pp$Sminus3,pp$Sminus2,pp$Sminus1, sep='')
  if(nHistory == 3) pp$pattern <- paste(pp$Sminus3,pp$Sminus2,pp$Sminus1, sep='')
  if(nHistory == 2) pp$pattern <- paste(pp$Sminus2,pp$Sminus1, sep='')
  pp <- pp[!pp$trials<6,]

  # define stim
  dat$stim <- dat$S == dat[dat$choice==1,'R'][1]
  pp$stim <- pp$S == dat[dat$choice==1,'R'][1]

  # define possible patterns
  allPatterns <- expand.grid(lapply(1:nHistory, function(x) return(unique(dat$Sminus1))), stringsAsFactors = FALSE)# expand.grid(unique(dat$Sminus1), unique(dat$Sminus1), unique(dat$Sminus1), unique(dat$Sminus1), stringsAsFactors = FALSE)
  # allPatterns <- allPatterns[,ncol(allPatterns):1]
  # allPatterns <- allPatterns[nrow(allPatterns):1,]
  allPatterns <- factor(apply(allPatterns, 1, paste0, collapse=''))# paste0(allPatterns$Var1, allPatterns$Var2, allPatterns$Var3 , allPatterns$Var4))

  # aggregate
  dataBySub <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*subjects, dat, mean)
  dataMeans <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, dataBySub, mean)
  dataSEs <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, dataBySub, function(x) sd(x)/sqrt(length(x)))
  dataMeans$pattern <- factor(dataMeans$pattern, levels=allPatterns)
  dataSEs$pattern <- factor(dataSEs$pattern, levels=allPatterns)

  dataBySubByChoice <- aggregate(cbind(rt,accuracy)~pattern*subjects*choice, dat, mean)
  dataMeansByChoice <- aggregate(cbind(rt,accuracy)~pattern*choice, dataBySubByChoice, mean)
  dataMeansByChoice$pattern <- factor(dataMeansByChoice$pattern, levels=allPatterns)

  ppByPostN <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*postn, pp, mean)
  ppQuantiles <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, ppByPostN, quantile, c(.025, .5, .975))
  ppQuantiles$pattern <- factor(ppQuantiles$pattern, levels=allPatterns)

  ppByChoiceByPostN <- aggregate(cbind(rt,accuracy,stim)~pattern*postn*choice, pp, mean)
  ppQuantilesByChoice <- aggregate(cbind(rt,accuracy,stim)~pattern*choice, ppByChoiceByPostN, quantile, c(.025, .5, .975))
  ppQuantilesByChoice$pattern <- factor(ppQuantilesByChoice$pattern, levels=allPatterns)

  # plot RTs
  par(mar=c(4,4,4,1))
  shift <- .1
  ylim <- range(c(dataMeansByChoice$rt+dataSEs$rt, dataMeansByChoice$rt-dataSEs$rt, ppQuantilesByChoice$rt[,1], ppQuantilesByChoice$rt[,3]))
  idx1 <- dataMeansByChoice$choice==1
  plot(x=as.numeric(dataMeansByChoice$pattern[idx1])+shift, y=dataMeansByChoice$rt[idx1], pch=4, lwd=2, ylim=ylim, xaxt='n', ylab='RT', xlab='Pattern', main=main1)
  axis(side=1, at=1:length(dataMeans$pattern), labels=levels(dataMeans$pattern), las=2)
  idx2 <- dataMeansByChoice$choice==0
  points(x=as.numeric(dataMeansByChoice$pattern[idx2])-shift, y=dataMeansByChoice$rt[idx2], pch=4, lwd=2,col=2)

  idx1 <- ppQuantilesByChoice$choice==1
  arrows(x0=as.numeric(ppQuantilesByChoice$pattern[idx1])+shift, y0=ppQuantilesByChoice$rt[idx1,1], y1=ppQuantilesByChoice$rt[idx1,3], code=3, angle=90, length=.05, lwd=2)
  idx2 <- ppQuantilesByChoice$choice==0
  arrows(x0=as.numeric(ppQuantilesByChoice$pattern[idx2])-shift, y0=ppQuantilesByChoice$rt[idx2,1], y1=ppQuantilesByChoice$rt[idx2,3], code=3, angle=90, length=.05, lwd=2,col=2)
  mtext(side = 3, line = 0.25, text=sub, cex=.66)
  legend('bottomright', c('Post. pred. (95CI)', 'Data',
                          paste0('choice=', as.character(dat[dat$choice==1,'R'][1])),
                          paste0('choice=', as.character(dat[dat$choice==0,'R'][1]))), lwd=c(2,2,2,2), lty=c(1, NA, 1, 1), pch=c(NA, 4, NA, NA), col=c(1,1,1,2), bty='n', cex=.75)
  # legend('topleft', c('repeat/alternate'), lty=1, col=2, bty='n')

  par(mar=c(4,4,4,1))
  shift = 0.1
  ylim = range(c(dataMeans$choice-dataSEs$choice, dataMeans$choice+dataSEs$choice, dataMeans$stim-dataSEs$stim, dataMeans$stim+dataSEs$stim,
                 ppQuantiles$choice))
  plot(x=as.numeric(dataMeans$pattern), y=dataMeans$choice, lwd=2, ylim=ylim, xaxt='n', ylab='Proportion', xlab='Pattern', pch=4, main=main2)
  axis(side=1, at=1:length(dataMeans$pattern), labels=levels(dataMeans$pattern), las=2)
  arrows(x0=as.numeric(ppQuantiles$pattern),
         y0=ppQuantiles$choice[,1],
         y1=ppQuantiles$choice[,3], angle=90, code=3, length=.05, lwd=2)
  #arrows(x0=as.numeric(dataMeans$pattern)-shift, y0=dataMeans$choice-dataSEs$choice, y1=dataMeans$choice+dataSEs$choice, angle=90, code=3, length=.1, lwd=2)
  #points(x=as.numeric(dataMeans$pattern)+shift, y=dataMeans$stim, lwd=2)
  #arrows(x0=as.numeric(dataMeans$pattern)+shift, y0=dataMeans$stim-dataSEs$stim, y1=dataMeans$stim+dataSEs$stim, angle=90, code=3, length=.1, lty=2, lwd=2)
  abline(h=0.5, lty=2, col='grey')
  legend('topleft', c(paste0('p(choice=', as.character(dat[dat$choice==1,'R'][1]), ')')), lwd=c(2), lty=c(1), col=1, bty='n')
  par(mar=c(5, 4, 4, 2) + 0.1)
}


# Waves / cosines -------------------------------------------------------------------
getWaves <- function(samplers, subject, filter='sample', pname='Bcos', nwaves=200) {
  idxStage<- samplers[[1]]$samples$stage==filter
  pnames <- samplers[[1]]$par_names
  idxPars <- grepl(pname, pnames)

  cossamps <- do.call(cbind, lapply(samplers, function(x) x$samples$alpha[idxPars,subject,idxStage]))

  dcts <- attr(samplers[[1]]$data[[subject]], 'adapt')[[1]]$dct
  if(nrow(cossamps) == 0) return(matrix(0, nrow=nrow(dcts), ncol=2))
  dcts <- dcts[,1:nrow(cossamps)]

  ## get 200 random waves
  if(ncol(cossamps) > nwaves) {
    cossamps <- cossamps[,sample(1:ncol(cossamps), size=nwaves)]
  }

  waves <- NULL
  for(i in 1:ncol(cossamps)) waves <- cbind(waves, dcts%*%cossamps[,i])
  return(waves)
}
plotWaves <- function(subject, samplers, pp, filter='sample') {

  par(mar=c(0,4,3,1)+.1)
  waves <- getWaves(samplers, subject, filter, pname='Bcos')
  plot(waves[,1], type='l', col=adjustcolor('black', alpha.f=.2), main=paste0('subject ', subject), ylim=range(waves), ylab='"Threshold drift"', xaxt='n')
  for(i in 2:ncol(waves)) lines(waves[,i], col=adjustcolor('black', alpha.f=.1))

  par(mar=c(3,4,0,1)+.1)
  waves <- getWaves(samplers, subject, filter, pname='vcos')
  plot(waves[,1], type='l', col=adjustcolor('black', alpha.f=.2), main='', ylab='"Rates drift/attention"', ylim=range(waves))
  for(i in 2:ncol(waves)) lines(waves[,i], col=adjustcolor('black', alpha.f=.1))

  par(mar=c(0,4,0,1)+.1)
  plot(dat[dat$subjects==subject,'rt'], ylab='RT (s)', xaxt='n')
  ppByTrial <- aggregate(rt~trials, pp[pp$subjects==subject,], quantile, c(0.025, 0.5, 0.975))
  polygon(c(ppByTrial$trials, rev(ppByTrial$trials)),
          c(ppByTrial$rt[,1], rev(ppByTrial$rt[,3])), col=adjustcolor(2, alpha.f=.5), border=adjustcolor(2, alpha.f=.5))
  legend('topleft', c('Data', 'Posterior predictions'), pch=c(1, NA), lty=c(NA, 1), col=c(1,22), bty='n')


  ppByTrial <- aggregate(accuracy~postn, pp[pp$subjects==subject,], ma, order=30)$accuracy
  quants <- apply(ppByTrial, 2, quantile, c(.025, .5, .975), na.rm=TRUE)
  trials <- which(!is.na(quants[1,]))

  par(mar=c(5,4,0,1)+.1)
  plot(ma(dat$accuracy[dat$subjects==subject], order=30), lwd=2, ylab='Accuracy (moving average)', xlab='Trial')
  lines(trials, quants[2,trials], col='red', lwd=2)
  polygon(c(trials, rev(trials)),
          c(quants[1,trials], rev(quants[3,trials])), col=adjustcolor(2, alpha.f=.5), border=adjustcolor(2, alpha.f=.5))
  #  ppByTrial <- aggregate(accuracy~trials, pp[pp$subjects==subject,], quantile, c(0.025, 0.5, 0.975))
}


# Post-error slowing ------------------------------------------------------
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


  f <- function(x, ...) quantile(x, probs=c(.1, .5, .9), ...)
  rtdf <- expand.grid(trialNposterror=1:5, probs=c(.1, .5, .9), measure='rt', value=NA)
  rtabsdf <- expand.grid(trialNposterror=seq(-3,7), probs=c(.1, .5, .9), measure='rtabs', value=NA)
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

getDistributionsToPlot <- function(dat, pp, pp2, averageOnly, plotpp) {
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

plotPostError8 <- function(distributions=NULL,
                           dat=NULL, pp=NULL, pp2=NULL,
                           orderSubjects=TRUE,
                           firstTrialOnly=FALSE,
                           use_mean=TRUE,
                           averageOnly=TRUE,
                           xlab1='Trial after error', xlab2='Trial after error',
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
    } else {
      axis(side=1, at=xs, labels=rep('', length(xs)))
    }

    ## data
    points(as.numeric(postpreData[postpreData$measure=='rt'&postpreData$trialNposterror==1,'subjects']),
           postpreData[postpreData$measure=='rt'&postpreData$trialNposterror==1,'value'], pch=4, lwd=2)

    ## PP
    arrows(x0=as.numeric(postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'subjects'])-.15*(!is.null(pp2)),
           y0=postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'value'][,1],
           y1=postprePPsummary[postprePPsummary$measure=='rt'&postprePPsummary$trialNposterror==1,'value'][,3],
           angle=90, code=3, length=.02, lwd=2, col=2)

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







# Power spectrum ----------------------------------------------------------
# Functions for FFT on data and models ------------------------------------
interpolatePowers <- function(freqBySubByPostN, powerBySubByPostN) {
  # make approxfuncs for every postN
  funcs <- mapply(function(x,y) approxfun(x=x, y=y), freqBySubByPostN$rt, powerBySubByPostN$rt)

  # find frequencies of subject/pp with fewest trials (and frequencies)
  frequencyRange <- range(unlist(freqBySubByPostN$rt))
  nFrequencies <- min(sapply(freqBySubByPostN$rt, length))
  frequencies <- freqBySubByPostN$rt[[which.min(sapply(freqBySubByPostN$rt, length))]]

  # interpolate at these frequencies
  interpolatedPowers <- lapply(funcs, function(x) x(frequencies))
  return(list(frequencies=frequencies, interpolatedPowers=interpolatedPowers))
}

# getPowerSpectra <- function(data, mean.pp=FALSE, by.postn=FALSE,
#                             spans=c(3,5), detrend=TRUE, demean=FALSE, get.log=FALSE) {
#   if(mean.pp) {
#     pp <- data
#     pp <- pp[order(pp$subjects, pp$postn, pp$trials),]
#     powerBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$spec))
#     freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
#
#     if(is.list(powerBySubByPostN$rt)) {
#       powerBySubByPostN$rt <- do.call(rbind, powerBySubByPostN$rt)
#     }
#     powerBySub <- lapply(unique(powerBySubByPostN$subjects), function(x) apply(powerBySubByPostN[powerBySubByPostN$subjects==x,'rt'],2,mean))
#     meanPower <- apply(do.call(rbind, powerBySub), 2, mean)
#
#     freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
#     if(is.list(freqBySubByPostN$rt)) {
#       freqBySubByPostN$rt <- do.call(rbind, freqBySubByPostN$rt)
#     }
#     freqBySub <- lapply(unique(freqBySubByPostN$subjects), function(x) apply(freqBySubByPostN[freqBySubByPostN$subjects==x,'rt'],2,mean))
#     meanFreqs <- apply(do.call(rbind, freqBySub), 2, mean)
#     return(data.frame(freq=meanFreqs, power=meanPower))
#   } else if(by.postn) {
#     pp <- data   # the 'data' that were passed are actually posterior predictives
#     pp <- pp[order(pp$subjects, pp$postn, pp$trials),]  # ensure correct ordering
#
#     # get spectra by subject by postN
#     spectraByPostN <- lapply(unique(pp$subjects), function(subject) lapply(unique(pp[pp$subjects==subject,'postn']), function(postn) {
#       spectrum(pp[pp$subjects==subject&pp$postn==postn,'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)
#     }))
#
#     freqBySubByPostN <- lapply(unique(pp$subjects), function(subject) spectraByPostN[[subject]][[1]]$freq)  # these should all be the same within each subject
#     powerBySubByPostN <- lapply(spectraByPostN, function(x) lapply(x, function(y) y$spec))
#
#
#
#     #freqBySubByPostN <- lapply(unique(pp$subjects), function(subject) spectraByPostN[[subject]][[1]]$spec)
#     #powerBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$spec))
#     #freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
#     if(is.list(powerBySubByPostN$rt)) {
#       ## varying numbers of trials across subjects -> varying numbers of estimated frequencies
#       ## need to interpolate to get everyone on the same frequency range
#       tmp <- interpolatePowers(freqBySubByPostN, powerBySubByPostN)
#       powerBySubByPostN$rt <- tmp[[2]]
#       meanFreqs <- tmp[[1]]
#       powerByPostN <- lapply(unique(powerBySubByPostN$postn), function(x) apply(do.call(rbind, powerBySubByPostN[powerBySubByPostN$postn==x,'rt']),2,mean))
#       } else {
#         powerByPostN <- lapply(unique(powerBySubByPostN$postn), function(x) apply(powerBySubByPostN[powerBySubByPostN$postn==x,'rt'],2,mean))
#         freqBySub <- lapply(unique(freqBySubByPostN$subjects), function(x) apply(freqBySubByPostN[freqBySubByPostN$subjects==x,'rt'],2,mean))
#         meanFreqs <- freqBySub[[1]] ## SHOULD BE THE SAME FOR ALL SUBS!
#       }
#     #    if(is.list(powerBySubByPostN$rt)) {
#     #      powerBySubByPostN$rt <- do.call(rbind, powerBySubByPostN$rt)
#     #    }
#     #    meanPower <- apply(do.call(rbind, powerByPostN), 2, mean)
#
#     #    freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, log=get.log)$freq))
#     #    if(is.list(freqBySubByPostN$rt)) {
#     #      freqBySubByPostN$rt <- do.call(rbind, freqBySubByPostN$rt)
#     #    }
#     #    freqByPostN <- lapply(unique(freqBySubByPostN$postn), function(x) apply(freqBySubByPostN[freqBySubByPostN$postn==x,'rt'],2,mean))
#     #    meanFreqs <- apply(do.call(rbind, freqByPostN), 2, mean)
#     return(list(meanFreqs, powerByPostN))
#   } else {
#     # estimate power spectrum by subject
#     spectra <- lapply(unique(data$subjects), function(x) spectrum(data[data$subjects==x, 'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log))
#     freqsBySub <- lapply(spectra, function(x) x$freq)
#     powerBySub <- lapply(spectra, function(x) x$spec)
#
#     # not all subjects have the same trial numbers.
#     # As a result, the actual frequencies at which we have a power estimate differ per subject, so we need to interpolate
#     approximators <- mapply(function(x,y) approxfun(x=x,y=y), freqsBySub, powerBySub)
#     frequencies <- freqsBySub[[which.min(sapply(freqsBySub, length))]]    # get the frequencies corresponding to the shortest range(?)
#     interpolatedPowers <- do.call(rbind, lapply(approximators, function(x) x(frequencies))) # this has shape [nSubjects x nFrequencies]
#
#     meanPower <- apply(interpolatedPowers, 2, mean)
#     return(data.frame(freq=frequencies, power=meanPower))
#   }
# }

getPowerSpectra <- function(data, mean.pp=FALSE, by.postn=FALSE,
                            spans=c(3,5), detrend=TRUE, demean=FALSE, get.log=FALSE) {
  if(mean.pp) {
    pp <- data
    pp <- pp[order(pp$subjects, pp$postn, pp$trials),]
    powerBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$spec))
    freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))

    if(is.list(powerBySubByPostN$rt)) {
      powerBySubByPostN$rt <- do.call(rbind, powerBySubByPostN$rt)
    }
    powerBySub <- lapply(unique(powerBySubByPostN$subjects), function(x) apply(powerBySubByPostN[powerBySubByPostN$subjects==x,'rt'],2,mean))
    meanPower <- apply(do.call(rbind, powerBySub), 2, mean)

    freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
    if(is.list(freqBySubByPostN$rt)) {
      freqBySubByPostN$rt <- do.call(rbind, freqBySubByPostN$rt)
    }
    freqBySub <- lapply(unique(freqBySubByPostN$subjects), function(x) apply(freqBySubByPostN[freqBySubByPostN$subjects==x,'rt'],2,mean))
    meanFreqs <- apply(do.call(rbind, freqBySub), 2, mean)
    return(data.frame(freq=meanFreqs, power=meanPower))
  } else if(by.postn) {
    pp <- data   # the 'data' that were passed are actually posterior predictives
    pp <- pp[order(pp$subjects, pp$postn, pp$trials),]  # ensure correct ordering

    # get spectra by subject by postN
    spectraByPostN <- mclapply(unique(pp$subjects), function(subject) lapply(unique(pp[pp$subjects==subject,'postn']), function(postn) {
      spectrum(pp[pp$subjects==subject&pp$postn==postn,'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)
    }), mc.cores=10)

    freqBySubByPostN <- lapply(unique(pp$subjects), function(subject) spectraByPostN[[subject]][[1]]$freq)  # these should all be the same within each subject
    powerBySubByPostN <- lapply(spectraByPostN, function(x) lapply(x, function(y) y$spec))
    npostn <- length(powerBySubByPostN[[1]])
    nsubs <- length(freqBySubByPostN)

    # Are are frequencies the same?
    nUniqueLengths <- length(unique(sapply(freqBySubByPostN, length)))
    if(nUniqueLengths > 1) {
      # We need to interpolate because varying subjects have varying trial numbers. Find the subject with fewest trials
      nFrequencies <- min(sapply(freqBySubByPostN, length))
      frequencies <- freqBySubByPostN[[which.min(sapply(freqBySubByPostN, length))]]

      # interpolation: ugly for loop, sorry
      interpolatedPowers <- vector(mode='list', length=nsubs)
      for(subject in 1:nsubs) {
        interpolatedPowers[[subject]] <- matrix(NA, nrow=npostn, ncol=nFrequencies)
        for(postn in 1:npostn) {
          f <- approxfun(x=freqBySubByPostN[[subject]], y=powerBySubByPostN[[subject]][[postn]])
          interpolatedPowers[[subject]][postn,] <- f(frequencies)
        }
      }

      # now get across-subject mean per posterior predictive
      meanpowerByPostN <- do.call(rbind, lapply(1:npostn, function(postn) apply(do.call(rbind, lapply(interpolatedPowers, function(x) x[postn,])), 2, mean)))
    } else {
      frequencies <- freqBySubByPostN[[1]]
      tmp <- lapply(powerBySubByPostN, function(x) do.call(rbind, x))
      meanpowerByPostN <- do.call(rbind, lapply(1:npostn, function(postn) apply(do.call(rbind, lapply(tmp, function(x) x[postn,])), 2, mean)))
    }
    return(list(freq=frequencies, power=meanpowerByPostN))
  } else {
    # estimate power spectrum by subject
    spectra <- lapply(unique(data$subjects), function(x) spectrum(data[data$subjects==x, 'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log))
    freqsBySub <- lapply(spectra, function(x) x$freq)
    powerBySub <- lapply(spectra, function(x) x$spec)

    # not all subjects have the same trial numbers.
    # As a result, the actual frequencies at which we have a power estimate differ per subject, so we need to interpolate
    approximators <- mapply(function(x,y) approxfun(x=x,y=y), freqsBySub, powerBySub)
    frequencies <- freqsBySub[[which.min(sapply(freqsBySub, length))]]    # get the frequencies corresponding to the shortest range(?)
    interpolatedPowers <- do.call(rbind, lapply(approximators, function(x) x(frequencies))) # this has shape [nSubjects x nFrequencies]

    meanPower <- apply(interpolatedPowers, 2, mean)
    return(data.frame(freq=frequencies, power=meanPower))
  }
}


makeChainPlots <- function(samplers, save_fn=NULL) {
  nchains <- chain_n(samplers)
  filter=colnames(nchains)[nchains[1,]>0][sum(nchains[1,]>0)]

  if(!is.null(save_fn)) pdf(file=save_fn, width=10, height=10)
  plot_chains(samplers, filter=filter, selection='mu', layout=c(3,2))
  plot_chains(samplers, filter=filter, selection='alpha', layout=c(4,4))
  plot_chains(samplers, filter=filter, selection='LL', layout=c(3,3))
  if(!is.null(save_fn)) dev.off()
}


makeSequentialPlots2 <- function(dat, pp, save_fn=NULL, include_posterror=TRUE) {
  # 2. sequential
  if(!is.null(save_fn)) pdf(file=save_fn, width=10, height=10)
  def.par <- par(no.readonly = TRUE)
  layout(matrix(c(1,2,3,4,5,5), byrow=TRUE, ncol=2)) #layout(rbind(matrix(c(1,2), ncol=2),3,4))
  plotSequentialEffects3(dat, pp, nHistory=3, main1=paste0('RT ~ stim hist'), main2='Accuracy ~ stim hist')
  if(include_posterror) {
    distributions <- getDistributionsToPlot(dat, pp=pp, pp2=NULL, averageOnly = FALSE, plotpp = TRUE)
    ## Plot 1: absolute RTs pre, error, post-error
    #nf <- layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
    plotPostError8(distributions, pp2=NULL, averageOnly = TRUE, main1='', ylim1=NULL, ylim2=NULL,
                   ylab1='RT (s)',ylab2='Accuracy',
                   xlab1='Trial after error', xlab2='Trial after error',
                   accmeasure='accabs', rtmeasure='rtabs', plotAccuracy = TRUE, setPar=FALSE)
    plotPostError8(distributions, pp2=NULL, averageOnly = FALSE, main1='', xlab1='Subject', ylim1=NULL, ylim2=NULL, ylab1='RT difference (s)', ylab2='Accuracy difference', setPar=FALSE)
  }
  par(def.par)
  if(!is.null(save_fn)) dev.off()
}



makeWavePlots <- function(dat, pp, samplers, save_fn=NULL) {
  if(!is.null(save_fn)) pdf(file=save_fn, width=10, height=10)
  def.par <- par(no.readonly = TRUE)
  layout(mat=matrix(c(1,2,3,4,5,6,7,8), nrow=4, ncol=2), heights=c(2,2,3,3))
  for(subject in names(samplers[[1]]$data)) {
    plotWaves(subject=subject, samplers=samplers, pp=pp)
  }
  par(def.par)
  if(!is.null(save_fn)) dev.off()
}


# plotSpectrum <- function(dat, pp, pp2=NULL, xlab='', ylab='', main='', add.legend=FALSE, ylim=NULL, plot.log=TRUE, xlim=NULL, detrend=FALSE) {
#   powerSpectraData <- getPowerSpectra(dat, detrend=detrend, demean=TRUE)
#   powerSpectraData <- powerSpectraData[order(powerSpectraData$freq),]
#
#   if(plot.log) { f <- function(x) log(x) } else { f <- function(x) x }
#   plot(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim)
#
#   if(!is.null(pp)) {
#     powerSpectraModel <- getPowerSpectra(pp, by.postn=TRUE, detrend=detrend, demean=TRUE)
#     for(i in 1:100) {
#       lines(x=f(powerSpectraModel[[1]]),
#             y=f(powerSpectraModel[[2]][i,]),
#             col=adjustcolor(3, alpha.f=.1))
#     }
#
#     if(!is.null(pp2)) {
#       powerSpectraModel2 <- getPowerSpectra(pp2, by.postn=TRUE, detrend=FALSE, demean=TRUE)
#       for(i in 1:100) {
#         lines(x=f(powerSpectraModel2[[1]]),
#               y=f(powerSpectraModel2[[2]][i,]),
#               col=adjustcolor(2, alpha.f=.1))
#       }
#     }
#
#     ## get mean of pp
#     freqs <- powerSpectraModel[[1]]
#     powers <- apply(powerSpectraModel[[2]],2,mean)
#     lines(x=f(freqs), y=f(powers), col='dark green')
#
#     if(!is.null(pp2)) {
#       freqs <- powerSpectraModel2[[1]]
#       powers <- apply(powerSpectraModel2[[2]],2,mean)
#       lines(x=f(freqs), y=f(powers), col='dark red')
#     }
#   }
#   lines(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), col=1)
#   if(!is.null(pp)) lines(x=f(freqs), y=f(powers), col='dark green')
# }
plotSpectrum <- function(dat, pp, pp2=NULL, xlab='', ylab='', main='', add.legend=FALSE, ylim=NULL, plot.log=TRUE, xlim=NULL, detrend=FALSE,
                         trial_duration=NULL) {
  powerSpectraData <- getPowerSpectra(dat, detrend=detrend, demean=TRUE)
  powerSpectraData <- powerSpectraData[order(powerSpectraData$freq),]

  if(plot.log) { f <- function(x) log(x) } else { f <- function(x) x }
  if(is.null(trial_duration)) {
    plot(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim)
  } else {
    # frequencies in Hz
    x_axis_ticks_seconds <- c(3600, 30*60, 15*60, 5*60, 2*60, 60, 30, 5, 1)
    x_axis_ticks_hz <- 1/x_axis_ticks_seconds
    x_axis_ticks_1_over_trial <- x_axis_ticks_hz*trial_duration    ## NB
    plot(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim, xaxt='n')
    axis(side=1, at=f(x_axis_ticks_1_over_trial), labels=c('1 hr', '30 m', '15 m', '5 m', '2 m', '1 m', '30 s', '5 s', '1 s'), las=2)
  }


  if(!is.null(pp)) {
    powerSpectraModel <- getPowerSpectra(pp, by.postn=TRUE, detrend=detrend, demean=TRUE)
    for(i in 1:100) {
      lines(x=f(powerSpectraModel[[1]]),
            y=f(powerSpectraModel[[2]][i,]),
            col=adjustcolor(3, alpha.f=.1))
    }

    if(!is.null(pp2)) {
      powerSpectraModel2 <- getPowerSpectra(pp2, by.postn=TRUE, detrend=FALSE, demean=TRUE)
      for(i in 1:100) {
        lines(x=f(powerSpectraModel2[[1]]),
              y=f(powerSpectraModel2[[2]][i,]),
              col=adjustcolor(2, alpha.f=.1))
      }
    }

    ## get mean of pp
    freqs <- powerSpectraModel[[1]]
    powers <- apply(powerSpectraModel[[2]],2,mean)
    lines(x=f(freqs), y=f(powers), col='dark green')

    if(!is.null(pp2)) {
      freqs <- powerSpectraModel2[[1]]
      powers <- apply(powerSpectraModel2[[2]],2,mean)
      lines(x=f(freqs), y=f(powers), col='dark red')
    }
  }
  lines(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), col=par()$col)

  if(!is.null(pp)) lines(x=f(freqs), y=f(powers), col='dark green')
}

# makeSpectrumPlots <- function(dat, pp, save_fn=NULL) {
#   if(!is.null(save_fn)) pdf(file=save_fn, width=10, height=10)
#   powerSpectraData <- getPowerSpectra(dat, detrend=FALSE, demean=TRUE)
#   powerSpectraData <- powerSpectraData[order(powerSpectraData$freq),]
#   powerSpectraModel <- getPowerSpectra(pp, by.postn=TRUE, detrend=FALSE, demean=TRUE)
#
#   plot(x=log(powerSpectraData$freq), y=log(powerSpectraData$power), xlab='Log frequency', ylab='log power', type='n', main='')
#   for(i in 1:100) {
#     lines(x=log(powerSpectraModel[[1]]),
#           y=log(powerSpectraModel[[2]][[i]]),
#           col=adjustcolor(2, alpha.f=.2))
#   }
#   lines(x=log(powerSpectraData$freq), y=log(powerSpectraData$power), col=1)
#
#   ## get mean of pp
#   freqs <- powerSpectraModel[[1]]
#   powers <- apply(do.call(rbind, powerSpectraModel[[2]]),2,mean)[order(freqs)]
#   freqs <- freqs[order(freqs)]
#   lines(x=log(freqs), y=log(powers), col='dark red')
#   legend('topright', c('Data', 'PP'), lwd=c(1,1), col=c(1,2), bty='n', cex=1)
#   if(!is.null(save_fn)) dev.off()
# }

makeDefaultPlots <- function(dat, pp, samplers, save_fn_template=paste0('./figures/dataset-', task, '/model-', decisionModel, '-', learningModel, '_dct-', dctModel, '_XXXX.pdf'),
                             include_posterror=TRUE) {

  # Plot all ------------------------------------------------------------------
  if(!is.null(save_fn_template)) {
    if(!dir.exists(dirname(save_fn_template))) dir.create(dirname(save_fn_template), recursive = TRUE)
    save_fn_chains = gsub('_XXXX.pdf', '_chains.pdf', save_fn_template)
    save_fn_fits = gsub('_XXXX.pdf', '_fit.pdf', save_fn_template)
    save_fn_sequential = gsub('_XXXX.pdf', '_sequential.pdf', save_fn_template)
    save_fn_waves = gsub('_XXXX.pdf', '_waves.pdf', save_fn_template)
    save_fn_spectrum = gsub('_XXXX.pdf', '_spectrum.pdf', save_fn_template)
  } else {
    save_fn_chains = save_fn_fits = save_fn_sequential = save_fn_waves = save_fn_spectrum = NULL
  }

  # Chains
  makeChainPlots(samplers, save_fn=save_fn_chains)

  ## effects of interest
  if(!is.null(save_fn_fits)) pdf(file=save_fn_fits, width=10, height=10)
  # 1. overall fit
  if('W' %in% colnames(samplers[[1]]$data[[1]])) {
    par(mfrow=c(2,2))
    plot_fit(dat, pp, factors = c('proportion_word', 'W', 'S'))
  } else if('E' %in% colnames(samplers[[1]]$data[[1]])) {
    par(mfrow=c(2,2))
    plot_fit(dat, pp, factors = c('E', 'S'))
  } else {
    par(mfrow=c(2,1))
    plot_fit(dat, pp, factors = c('S'))
  }
  if(!is.null(save_fn_fits)) dev.off()

  # 2. sequential
  makeSequentialPlots2(dat=dat, pp=pp, save_fn=save_fn_sequential, include_posterror=include_posterror)

  # 4. Slow waves
  makeWavePlots(dat=dat, pp=pp, samplers=samplers, save_fn=save_fn_waves)

  # Power spectra
  if(!is.null(save_fn_spectrum)) pdf(file=save_fn_spectrum, width=10, height=10)
  plotSpectrum(dat=dat, pp=pp, xlab='Log frequency', ylab='Log power')
  dev.off()
}

save_model_comparisons <- function(task, overwrite=FALSE) {
  # Model comparison --------------------------------------------------------
  save_fn <- paste0('./figures/dataset-', task, '/modelcomparison.csv')
  if(!file.exists(save_fn) | overwrite) {
    allSamplers <- list()
    all_fns <- Sys.glob(paste0('./samples/', task, '*.RData'))
    all_fns <- all_fns[!sapply(all_fns, function(x) grepl('_pp.RData', x))]
    for(fn in all_fns) {
      load(fn)
      modelName <- strsplit(strsplit(fn, 'model-')[[1]][2], '.RData')[[1]][1]
      nchains <- chain_n(samplers)
      if(colnames(nchains)[nchains[1,]>0][sum(nchains[1,]>0)]=='sample') {
         allSamplers[[modelName]] <- samplers
      }
    }
    mc <- EMC2:::compare_IC(allSamplers)
    write.csv(x=round(mc,3), file=save_fn)
    save(mc, file=gsub('.csv', '.RData', save_fn))
  }
}
