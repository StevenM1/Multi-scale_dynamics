## Figures for paper
rm(list=ls())
library(EMC2)
library(emcAdapt)
library(parallel)
source('./extra_EMC2_functions/adaptive.R')
source('./extra_EMC2_functions/make_data.R')
source('./extra_EMC2_functions/model_RDMdynamic.R')
source('./extra_EMC2_functions/utils.R')
figures_dir = './figures'


## We unfortunately need all plotting code here...
load_save_spectrum <- function(dat, dat_fn, by.postn=FALSE, detrend, demean, overwrite=FALSE) {
  powerSpectraData <- NULL
  if(!is.null(dat_fn) & !overwrite) {
    fourier_fn <- gsub('.RData', '_fourier.RData', dat_fn)
    if(file.exists(fourier_fn)) {
      powerSpectraData <- EMC2:::loadRData(fourier_fn)
    }
  }

  if(is.null(powerSpectraData)) {
    powerSpectraData <- getPowerSpectra(dat, by.postn = by.postn, detrend=detrend, demean=demean)
    if(!by.postn) powerSpectraData <- powerSpectraData[order(powerSpectraData$freq),]  # only order data(?)
    if(!is.null(dat_fn)) {
      fourier_fn <- gsub('.RData', '_fourier.RData', dat_fn)
      save(powerSpectraData, file=fourier_fn)
    }
  }
  return(powerSpectraData)
}

## Fourier plots
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

plotSpectrum <- function(dat, pp, pp2=NULL,
                         dat_fn=NULL, pp_fn=NULL, pp2_fn=NULL,
                         xlab='', ylab='', main='', add.legend=FALSE, ylim=NULL, plot.log=TRUE, xlim=NULL, detrend=FALSE,
                         demean=TRUE,
                         trial_duration=NULL, full_x=FALSE, plot_xticklabels=TRUE) {
  powerSpectraData <- load_save_spectrum(dat, dat_fn, detrend=detrend, demean=demean)

  if(plot.log) { f <- function(x) log(x) } else { f <- function(x) x }
  if(is.null(trial_duration)) {
    plot(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim)
  } else {
    # frequencies in Hz
    x_axis_ticks_seconds <- c(3600, 30*60, 15*60, 5*60, 2*60, 60, 30, 5, 1)
    x_axis_ticks_hz <- 1/x_axis_ticks_seconds
    x_axis_ticks_1_over_trial <- x_axis_ticks_hz*trial_duration    ## NB
    # ylim <- range(f(powerSpectraData$power))

    # which()
    powers <- f(powerSpectraData$power)
    freqs <- f(powerSpectraData$freq)

    if(!full_x) {
      idx <- freqs>=(log(x_axis_ticks_1_over_trial)[4]-0.2)
    } else {
      idx <- rep(TRUE, length(powerSpectraData$freq))
    }
    #    xlim <- c(log(x_axis_ticks_1_over_trial)[4], max(freqs))
    #    ylim <- range(powers[freqs>=xlim[1]], na.rm=TRUE)

    plot(x=f(powerSpectraData$freq)[idx], y=f(powerSpectraData$power)[idx], xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim, xaxt='n')
    if(plot_xticklabels) {
      axis(side=1, at=log(x_axis_ticks_1_over_trial), labels=c('1 hr', '30 m', '15 m', '5 m', '2 m', '1 m', '30 s', '5 s', '1 s'), las=2)
    } else {
      axis(side=1, at=log(x_axis_ticks_1_over_trial), labels=rep('', length(x_axis_ticks_1_over_trial)), las=2)
    }
  }


  if(!is.null(pp)) {
    powerSpectraModel <- load_save_spectrum(pp, pp_fn, by.post=TRUE, detrend=detrend, demean=demean)
    #    powerSpectraModel <- powerSpectraModel(pp, by.postn=TRUE, detrend=detrend, demean=demean)
    for(i in 1:100) {
      lines(x=f(powerSpectraModel[[1]]),
            y=f(powerSpectraModel[[2]][i,]),
            col=adjustcolor(3, alpha.f=.1))
    }

    if(!is.null(pp2)) {
      powerSpectraModel2 <- load_save_spectrum(pp2, pp2_fn, by.post=TRUE, detrend=detrend, demean=demean)
      #      powerSpectraModel2 <- getPowerSpectra(pp2, by.postn=TRUE, detrend=detrend, demean=demean)
      for(i in 1:100) {
        lines(x=f(powerSpectraModel2[[1]]),
              y=f(powerSpectraModel2[[2]][i,]),
              col=adjustcolor(3, alpha.f=.1))
      }
    }

    ## get mean of pp
    freqs <- powerSpectraModel[[1]]
    powers <- apply(powerSpectraModel[[2]],2,mean)
    lines(x=f(freqs), y=f(powers), col='dark red')

    if(!is.null(pp2)) {
      freqs2 <- powerSpectraModel2[[1]]
      powers2 <- apply(powerSpectraModel2[[2]],2,mean)
      lines(x=f(freqs2), y=f(powers2), col='dark red')
    }
  }
  lines(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), col=par()$col)

  if(!is.null(pp)) lines(x=f(freqs), y=f(powers), col='dark green')
  if(!is.null(pp2)) lines(x=f(freqs2), y=f(powers2), col='dark red')
}

plot_fit <- function(data,pp,subject=NULL,factors=NULL,functions=NULL,
                     stat=NULL,stat_name="",adjust=1,
                     quants=c(.025,.5,.975),do_plot=TRUE,
                     xlim=NULL,ylim=NULL,
                     layout=NULL,mfcol=FALSE,
                     probs=c(1:99)/100,
                     data_lwd=2,fit_lwd=1,
                     q_points=c(.1,.3,.5,.7,.9),
                     qp_cex=1,pqp_cex=.5,lpos="topleft", main = "")
{
  if (!is.null(stat) & is.null(factors)) factors <- NA
  if (!is.null(subject)) {
    snams <- levels(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject(s) not present\n")
    dat <- droplevels(data[data$subjects %in% subject,])
    pp <- droplevels(pp[pp$subjects %in% subject,])
    if (length(subject>1))
      fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))] else
        fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt"))]
  } else {
    dat <- data
    fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))]
  }

  okd <- !is.na(dat$R) & is.finite(dat$rt)
  okpp <- !is.na(pp$R) & is.finite(pp$rt)

  if (!is.null(functions)) for (i in 1:length(functions)) {
    dat <- cbind.data.frame(functions[[i]](dat),dat)
    names(dat)[1] <- names(functions)[i]
    pp <- cbind.data.frame(functions[[i]](pp),pp)
    names(pp)[1] <- names(functions)[i]
    fnams <- c(names(functions)[i],fnams)
  }

  if (!is.null(factors)) {
    if (any(is.na(factors))) fnams <- NA else {
      if (!all(factors %in% fnams))
        stop("factors must name factors in data")
      fnams <- factors
    }
  }
  if (!any(is.na(layout))) if (!is.null(layout))
    if (mfcol) par(mfcol=layout) else par(mfrow=layout)

  if (all(is.na(data$rt))) stop("Use plot_fit_choice if no rt data")

  if (!any(is.na(fnams))) {
    cells <- dat[,fnams,drop=FALSE]
    for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
    cells <- apply(cells,1,paste,collapse=" ")
    pp_cells <- pp[,fnams,drop=FALSE]
    for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
    pp_cells <- apply(pp_cells,1,paste,collapse=" ")
  }
  if (!is.null(stat)) { # statistic
    postn <- unique(pp$postn)
    if (any(is.na(fnams))) ucells <- "" else ucells <- sort(unique(cells))
    tab <- matrix(nrow=length(ucells),ncol=4,
                  dimnames=list(ucells,c("Observed",names(quantile(1:5,quants)))))
    for (i in ucells) {
      if (i=="") {
        dati <- dat
        ppi <- pp
        obs <- stat(dati)
        pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
        tab[1,] <- c(obs,quantile(pred,quants))
      } else {
        dati <- dat[cells==i,]
        ppi <- pp[pp_cells==i,]
        obs <- stat(dati)
        pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
        tab[i,] <- c(obs,quantile(pred,quants))
      }
      if (do_plot) {
        dens <- density(pred,adjust=adjust)
        if (!is.null(xlim)) xlimi <- xlim else
          xlimi <- c(pmin(obs,min(dens$x)),pmax(obs,max(dens$x)))
        plot(dens,main=paste(main,i),xlab=stat_name,xlim=xlimi)
        abline(v=obs)
      }
    }
    invisible(tab)
  } else { # cdf
    if (any(is.na(fnams))) cells <- ""
    pok <- probs %in% q_points
    R <- levels(dat$R)
    if (is.null(ylim)) ylim <- c(0,1)
    # set common xlim
    if (is.null(xlim)) {
      xlim <- c(Inf,-Inf)
      for (i in sort(unique(cells))) {
        if (i=="") {
          dati <- dat
          ppi <- pp
        } else {
          dati <- dat[cells==i & okd,]
          ppi <- pp[pp_cells==i & okpp,]
        }
        pqs <- pq <- qs <- setNames(vector(mode="list",length=length(R)),R)
        for (j in R) if (length(dati$rt[dati$R==j])>=length(q_points)) {
          qs[[j]] <- quantile(dati$rt[dati$R==j],probs=probs)
          pq[[j]] <- quantile(ppi$rt[ppi$R==j],probs=probs)
          pqs[[j]] <- tapply(ppi$rt[ppi$R==j],ppi$postn[ppi$R==j],
                             quantile,probs=probs[pok])
        } else qs[[j]] <- pq[[j]] <- pqs[[j]] <- NA
        rx <- cbind(do.call(rbind,lapply(qs,function(x){x[c(1,length(probs))]})),
                    do.call(rbind,lapply(pq,function(x){x[c(1,length(probs))]})))
        xlimi <- c(min(rx,na.rm=TRUE),max(rx,na.rm=TRUE))
        if (!any(is.na(xlimi))) {
          xlim[1] <- pmin(xlim[1],xlimi[1])
          xlim[2] <- pmax(xlim[2],xlimi[2])
        }
      }
    }
    for (i in sort(unique(cells))) {
      if (i=="") {
        dati <- dat
        ppi <- pp
        okdi <- okd
        okppi <- okpp
      } else {
        dati <- dat[cells==i,]
        ppi <- pp[pp_cells==i,]
        okdi <- okd[cells==i]
        okppi <- okpp[pp_cells==i]
      }
      pR <- tapply(okdi,dati$R,sum)/dim(dati)[1]
      ppR <- tapply(okppi,ppi$R,sum)/dim(ppi)[1]
      dati <- dati[okdi,]
      ppi <- ppi[okppi,]
      pqs <- pq <- qs <- setNames(vector(mode="list",length=length(R)),R)
      for (j in R) if (length(dati$rt[dati$R==j])>=length(q_points)) {
        isj <- ppi$R==j
        qs[[j]] <- quantile(dati$rt[dati$R==j],probs=probs)
        pq[[j]] <- quantile(ppi$rt[isj],probs=probs)
        pqs[[j]] <- tapply(ppi$rt[isj],ppi$postn[isj],quantile,probs=probs[pok])
      } else qs[[j]] <- pq[[j]] <- pqs[[j]] <- NA
      if ( !any(is.na(pq[[1]])) ) {
        plot(pq[[1]],probs*ppR[1],xlim=xlim,ylim=ylim,main=main,#paste(main,i),
             xlab="RT",type="l",
             lwd=fit_lwd,ylab="",lty=1, col='lightgreen') #"p(R)",lty=1)
        tmp=lapply(pqs[[1]],function(x){
          points(x,probs[pok]*ppR[1],col='darkgreen',bg=adjustcolor('lightgreen', alpha.f=.8),pch=21,cex=pqp_cex)})
        points(pq[[1]][pok],probs[pok]*ppR[1],cex=pqp_cex*3,pch=21,bg="lightgreen", col="darkgreen")
        #   points(x,probs[pok]*ppR[1],col=adjustcolor('lightgreen', alpha.f=.8),pch=16,cex=pqp_cex)})
        # points(pq[[1]][pok],probs[pok]*ppR[1],cex=pqp_cex*3,pch=16,col="lightgreen")
        lines(qs[[1]],probs*pR[1],lwd=data_lwd,lty=1)
        points(qs[[1]][pok],probs[pok]*pR[1],cex=qp_cex,pch=16)
        do_plot=FALSE
      } else do_plot=TRUE
      if (length(qs)>1) {
        for (j in 2:length(qs)) if (!any(is.na(pq[[j]]))) {
          if (do_plot) {
            plot(pq[[j]],probs*ppR[j],xlim=xlim,ylim=ylim,main=main, #paste(main,i),
                 xlab="RT",type="l",
                 lwd=fit_lwd,ylab="", lty=j, col='lightgreen') #p(R)",lty=j)
            do_plot <- FALSE
          } else lines(pq[[j]],probs*ppR[j],lwd=fit_lwd,lty=j)
          tmp=lapply(pqs[[j]],function(x){
            #   points(x,probs[pok]*ppR[j],col=adjustcolor('lightgreen', alpha.f=.8),pch=16,cex=pqp_cex)})
            # points(pq[[j]][pok],probs[pok]*ppR[j],cex=pqp_cex*3,pch=16,col='lightgreen')
            points(x,probs[pok]*ppR[j],col='darkgreen', bg=adjustcolor('lightgreen', alpha.f=.8),pch=21,cex=pqp_cex)})
          points(pq[[j]][pok],probs[pok]*ppR[j],cex=pqp_cex*3,pch=21,bg='lightgreen', col='darkgreen')
          #   points(x,probs[pok]*ppR[1],col='darkgreen',bg=adjustcolor('lightgreen', alpha.f=.8),pch=21,cex=pqp_cex)})
          # points(pq[[1]][pok],probs[pok]*ppR[1],cex=pqp_cex*3,pch=21,bg="lightgreen", col="darkgreen")
          lines(qs[[j]],probs*pR[j],lwd=data_lwd,lty=j)
          points(qs[[j]][pok],probs[pok]*pR[j],cex=qp_cex,pch=16)
        }
      }
      if(!is.null(lpos)) legend(lpos,as.character(R),lty=1:length(R),bty="n")#,title="Response")
    }
  }
}


library(Hmisc)
library(stringr)

# Sequential effects plots
addPreviousTrialInfo <- function(dat, dataColumn='S', nHistory=3, repeatAlternate=FALSE) {
  if(!'trials' %in% colnames(dat)) {
    dat <- EMC2:::add_trials(dat)
  }
  if('postn' %in% colnames(dat)) {
    # data passed is posterior predictives
    # triple check to make sure trial order is correct
    dat <- dat[order(dat$subjects,dat$postn,dat$trials),]
  }

  # Previous stimuli
  for(i in (nHistory+1):0) {
    dat[,paste0('Sminus', i)] <- Hmisc:::Lag(dat[,dataColumn], i)
  }

  if(repeatAlternate) {
    # recode into repeat / alternate
    tmp <- dat
    for(i in (nHistory):0) {
      tmp[,paste0('Sminus',i)] <- ifelse(dat[,paste0('Sminus',i+1)]==dat[,paste0('Sminus',i)], 'R', 'A')
    }
    # code a response
    # does the response repeat last stimulus, or alternate?
    tmp[,'R_repeat_alternate'] <- ifelse(tmp$R == dat$Sminus1, 'R', 'A')

    dat[,c(paste0('Sminus', (nHistory+1):0), 'R_repeat_alternate')] <- tmp[,c(paste0('Sminus', (nHistory+1):0), 'R_repeat_alternate')]
  }

  for(i in (nHistory+1):0) {
    dat[,paste0("Sminus",i)] <- str_sub(dat[,paste0("Sminus",i)], 1, 1)
  }

  dat$pattern <- apply(dat[,paste0('Sminus', nHistory:1)], 1, paste, collapse='')

  # remove trials without previous history. NB: for repeatAlternate, one *more* trial has to be removed; nHistory = 3 essentially looks 4 trials back
  dat <- dat[!dat$trials<(nHistory+1+repeatAlternate),]
  return(dat)
}

getSequentialEffectsBySub <- function(dat, pp, dataColumn='S', nHistory=3, repeatAlternate=FALSE, use_median=FALSE) {
  dat <- addPreviousTrialInfo(dat, dataColumn=dataColumn, nHistory=nHistory, repeatAlternate=repeatAlternate)
  pp <- addPreviousTrialInfo(pp, dataColumn=dataColumn, nHistory=nHistory, repeatAlternate=repeatAlternate)

  # define stim, choice in terms of repeat / alternate
  if(repeatAlternate) {
    # 1 = repeat, 0 = alternate
    dat$stim <- as.numeric(dat$Sminus0 == 'R')
    dat$choice <- as.numeric(dat$R_repeat_alternate == 'R')  # the response REPEATS the previous STIMULUS
    pp$stim <- as.numeric(pp$Sminus0 == 'R')
    pp$choice <- as.numeric(pp$R_repeat_alternate == 'R')
  } else {
    # stim is the same as the choice of R level 1
    dat$stim <- dat$S == dat[dat$choice==1,'R'][1]
    pp$stim <- pp$S == dat[dat$choice==1,'R'][1]
  }

  # define possible patterns
  uniqueStim <- sort(unique(dat$Sminus1))
  #  if('w' %in% uniqueStim) uniqueStim <- uniqueStim[c(2,1)]  # flip so order of plotting is consistent
  allPatterns <- expand.grid(lapply(1:nHistory, function(x) return(uniqueStim)), stringsAsFactors = FALSE)
  allPatterns <- factor(apply(allPatterns, 1, paste0, collapse=''))
  if('nnn' %in% allPatterns) allPatterns <- factor(c('nnn', 'wnn', 'nwn', 'nnw', 'wwn', 'wnw', 'nww', 'www'), levels=c('nnn', 'wnn', 'nwn', 'nnw', 'wwn', 'wnw', 'nww', 'www'))

  # Lots of aggregations. First, the data
  dataBySub <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*subjects, dat, mean)
  if(use_median) dataBySub$rt <- aggregate(rt~pattern*subjects,dat,median)$rt   # overwrite RT with median instead of mean
  dataMeans <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, dataBySub, mean)
  dataSEs <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, dataBySub, function(x) sd(x)/sqrt(length(x)))
  dataBySub$pattern <- factor(dataBySub$pattern, levels=allPatterns)
  dataMeans$pattern <- factor(dataMeans$pattern, levels=allPatterns)
  dataSEs$pattern <- factor(dataSEs$pattern, levels=allPatterns)

  # Split data also by choice
  dataBySubByChoice <- aggregate(cbind(rt,accuracy)~pattern*subjects*choice, dat, mean)
  if(use_median) dataBySubByChoice$rt <- aggregate(rt~pattern*subjects*choice,dat,median)$rt  # overwrite RT with median instead of mean
  dataBySubByChoice$pattern <- factor(dataBySubByChoice$pattern, levels=allPatterns)
  dataMeansByChoice <- aggregate(cbind(rt,accuracy)~pattern*choice, dataBySubByChoice, mean)
  dataMeansByChoice$pattern <- factor(dataMeansByChoice$pattern, levels=allPatterns)

  # Next, the posterior predictives
  # Group-level, across choices
  ppByPostN <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*subjects*postn, pp, mean)
  if(use_median) ppByPostN$rt <- aggregate(rt~pattern*subjects*postn, pp, median)$rt    # overwrite with median
  ppByPostN <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*postn, ppByPostN, mean)
  ppQuantiles <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, ppByPostN, quantile, c(.025, .5, .975))
  ppQuantiles$pattern <- factor(ppQuantiles$pattern, levels=allPatterns)

  # Group-level, by choice
  ppByChoiceByPostN <- aggregate(cbind(rt,accuracy,stim)~pattern*postn*choice, pp, mean)
  if(use_median) ppByChoiceByPostN$rt <- aggregate(rt~pattern*postn*choice, pp, median)$rt
  ppQuantilesByChoice <- aggregate(cbind(rt,accuracy,stim)~pattern*choice, ppByChoiceByPostN, quantile, c(.025, .5, .975))
  ppQuantilesByChoice$pattern <- factor(ppQuantilesByChoice$pattern, levels=allPatterns)

  # Subject-level, by choice
  ppBySubByChoiceByPostN <- aggregate(cbind(rt,accuracy,stim)~pattern*subjects*postn*choice, pp, mean)
  if(use_median) ppBySubByChoiceByPostN$rt <- aggregate(rt~pattern*subjects*postn*choice, pp, median)$rt

  #  ppBySubByChoiceQuantiles <- aggregate(cbind(rt,accuracy,stim)~pattern*postn*choice, ppBySubByChoiceByPostN, quantile, c(.025, .5, .975))
  #  ppBySubByChoiceQuantiles$pattern <- factor(ppBySubByChoiceQuantiles$pattern, levels=allPatterns)

  return(list(
    dat=dat,

    dataMeans=dataMeans,
    dataSEs=dataSEs,
    dataBySub=dataBySub,

    # Data by choice, group and subjects
    dataMeansByChoice=dataMeansByChoice,
    dataBySubByChoice=dataBySubByChoice,

    # Posterior predictives group
    ppQuantiles=ppQuantiles,
    # Group by choice
    ppQuantilesByChoice=ppQuantilesByChoice,

    # Subject by choice
    ppBySubByChoiceByPostN=ppBySubByChoiceByPostN,
    pp=pp  # needed for choice by subject later...
  ))
}

plotSequentialEffectsGroup <- function(dat=dat, pp=pp, seqEffectsSummary=NULL, nHistory=3, plotStimulusPattern=FALSE, repeatAlternate=FALSE,
                                       dataColumn='S', sub='', main1='', main2='', ylab1='', ylab2='', full.legend=TRUE, use_median=FALSE) {
  if(is.null(seqEffectsSummary)) {
    seqEffectsSummary <- getSequentialEffectsBySub(dat=dat, pp=pp, nHistory=3, dataColumn=dataColumn, repeatAlternate=repeatAlternate)
  }
  # unpack
  dat <- seqEffectsSummary$dat
  dataMeans = seqEffectsSummary$dataMeans
  dataSEs=seqEffectsSummary$dataSEs
  dataMeansByChoice = seqEffectsSummary$dataMeansByChoice
  ppQuantiles=seqEffectsSummary$ppQuantiles
  ppQuantilesByChoice=seqEffectsSummary$ppQuantilesByChoice
  ppBySubByChoiceQuantiles=seqEffectsSummary$ppBySubByChoiceQuantiles

  response_column_for_legend <- ifelse(repeatAlternate, 'R_repeat_alternate', 'R')

  if(use_median) meanRTsbyChoice <- aggregate(rt~choice, aggregate(rt~subjects*choice,dat,median),mean) else  meanRTsbyChoice <- aggregate(rt~choice, aggregate(rt~subjects*choice,dat,mean),mean)

  # plot RTs
  par(mar=c(1,2,3,1))
  shift <- .1
  ylim <- range(c(dataMeansByChoice$rt+dataSEs$rt, dataMeansByChoice$rt-dataSEs$rt, ppQuantilesByChoice$rt[,1], ppQuantilesByChoice$rt[,3])) #*c(.985, 1)
  if(ylim[1] > 0.7) ylim[1] <- 0.74  # dataset 3 needs some tweaking to prevent overlap of legend label and data
  if(ylim[2] < 0.54) ylim[1] <- 0.445  # dataset 2 also needs some tweaking to prevent overlap of legend label and data


  # data
  idx1 <- dataMeansByChoice$choice==1
  plot(x=as.numeric(dataMeansByChoice$pattern[idx1])+shift, y=dataMeansByChoice$rt[idx1], pch=4, lwd=2, ylim=ylim, xaxt='n', ylab=ylab1, xlab='', main=main1, xlim=range(as.numeric(dataMeansByChoice$pattern[idx1]))+c(-.2, .2))
  mtext(ylab1, side=2, line=2.5, cex=.66)
  axis(side=1, at=1:length(dataMeans$pattern), labels=levels(dataMeans$pattern), las=2)
  abline(h=meanRTsbyChoice[meanRTsbyChoice$choice==1,'rt'], col=par()$col, lwd=1, lty=2)
  idx2 <- dataMeansByChoice$choice==0
  points(x=as.numeric(dataMeansByChoice$pattern[idx2])-shift, y=dataMeansByChoice$rt[idx2], pch=4, lwd=2,col=2)
  abline(h=meanRTsbyChoice[meanRTsbyChoice$choice==0,'rt'], col=2, lwd=1, lty=2)


  # pp
  idx1 <- ppQuantilesByChoice$choice==1
  xs <- as.numeric(ppQuantilesByChoice$pattern[idx1])+shift
  arrows(x0=xs, y0=ppQuantilesByChoice$rt[idx1,1], y1=ppQuantilesByChoice$rt[idx1,3], code=3, angle=90, length=.05, lwd=2)
  lines(xs[order(xs)], ppQuantilesByChoice$rt[idx1,2][order(xs)], lwd=1)

  idx2 <- ppQuantilesByChoice$choice==0
  xs <- as.numeric(ppQuantilesByChoice$pattern[idx2])-shift
  arrows(x0=as.numeric(ppQuantilesByChoice$pattern[idx2])-shift, y0=ppQuantilesByChoice$rt[idx2,1], y1=ppQuantilesByChoice$rt[idx2,3], code=3, angle=90, length=.05, lwd=2,col=2)
  lines(xs[order(xs)], ppQuantilesByChoice$rt[idx2,2][order(xs)], lwd=1, col=2)
  mtext(side = 3, line = 0.25, text=sub, cex=.66)
  if(full.legend) {
    legend('bottomleft', c('Post. pred. (95% CI)', 'Data',
                           paste0('choice=', as.character(dat[dat$choice==1,response_column_for_legend][2])),
                           paste0('choice=', as.character(dat[dat$choice==0,response_column_for_legend][2]))), lwd=c(2,2,2,2), lty=c(1, NA, 1, 1), pch=c(NA, 4, NA, NA), col=c(par()$col,par()$col,par()$col,2), bty='n', cex=.75)
  } else {
    legend('bottomleft', c(paste0('choice=', as.character(dat[dat$choice==1,response_column_for_legend][2])),
                           paste0('choice=', as.character(dat[dat$choice==0,response_column_for_legend][2]))), lwd=c(2,2), lty=c(1, 1), pch=c(NA, NA), col=c(par()$col,2), bty='n', cex=.75)
  }

  par(mar=c(5,2,2,1))
  shift = 0.1
  ylim = range(c(dataMeans$choice-dataSEs$choice, dataMeans$choice+dataSEs$choice, dataMeans$stim-dataSEs$stim, dataMeans$stim+dataSEs$stim,
                 ppQuantiles$choice))
  if(plotStimulusPattern) {
    plot(x=as.numeric(dataMeans$pattern), y=dataMeans$stim, lwd=2, ylim=ylim, xaxt='n', ylab=ylab2, xlab='', pch=4, main=main2)
    warning('Plotting *stim* as a function of previous stimuli')
  } else {
    plot(x=as.numeric(dataMeans$pattern), y=dataMeans$choice, lwd=2, ylim=ylim, xaxt='n', ylab=ylab2, xlab='', pch=4, main=main2)
  }
  mtext(ylab2, side=2, line=2.5, cex=.66)
  mtext('Pattern', side=1, line=2.5, cex=.66)
  axis(side=1, at=1:length(dataMeans$pattern), labels=levels(dataMeans$pattern), las=2)
  xs <- as.numeric(ppQuantiles$pattern)
  arrows(x0=as.numeric(ppQuantiles$pattern),
         y0=ppQuantiles$choice[,1],
         y1=ppQuantiles$choice[,3], angle=90, code=3, length=.05, lwd=2)
  lines(xs[order(xs)], ppQuantiles$choice[,2][order(xs)])
  abline(h=0.5, lty=2, col=1)
  legend('topright', c(paste0('p(choice=', as.character(dat[dat$choice==1,response_column_for_legend][2]), ')')), lwd=c(2), lty=c(1), col=par()$col, bty='n', cex=.75)
  # par(mar=c(4, 4, 4, 2) + 0.1)
}

plotSequentialEffectsIndividual <- function(seqEffectsSummaryRA, xlab1='', ylab1='', xlab2='', ylab2='', main='', orderBy='rtchoice') {
  # map back 'response' to choice, easier to read
  dataBySubByChoice <- seqEffectsSummaryRA$dataBySubByChoice
  # get differences
  dataBySubByChoice <- dataBySubByChoice[dataBySubByChoice$pattern=='RRR',]
  dataBySubByChoiceW <- reshape(dataBySubByChoice, direction='wide', idvar=c('subjects', 'pattern'), timevar='choice', v.names=c('rt', 'accuracy'))
  dataBySubByChoiceW$rtDifference <- dataBySubByChoiceW$rt.0 - dataBySubByChoiceW$rt.1  # alternative minus repeat. This should be *positive* overall, as repeats are faster
  #dataBySubByChoiceW$rt.1 - dataBySubByChoiceW$rt.0  # repeat minus alternate. This should be *negative* overall, as alternates are slower

  # now choice
  dataBySub <- seqEffectsSummaryRA$dataBySub
  dataBySub <- dataBySub[dataBySub$pattern=='RRR',]  # choice is all we need, no more reshaping

  # pp
  ppBySubByChoiceByPostN <- seqEffectsSummaryRA$ppBySubByChoiceByPostN
  ppBySubByChoiceByPostN <- ppBySubByChoiceByPostN[ppBySubByChoiceByPostN$pattern=='RRR',]
  ppBySubByChoiceByPostNW <- reshape(ppBySubByChoiceByPostN, direction='wide', idvar=c('subjects', 'pattern', 'postn'), timevar='choice', v.names=c('rt', 'accuracy', 'stim'))
  # calculate rt differences
  ppBySubByChoiceByPostNW$rtDifference <- ppBySubByChoiceByPostNW$rt.0 - ppBySubByChoiceByPostNW$rt.1  # alternative minus repeat. This should be *positive* overall, as repeats are faster
  #ppBySubByChoiceByPostNW$rt.1 - ppBySubByChoiceByPostNW$rt.0  # repeat minus alternate. This should be *negative* overall, as alternates are slower
  ppBySubByChoiceQuantiles <- aggregate(rtDifference~subjects, ppBySubByChoiceByPostNW, quantile, c(0.025, 0.5, 0.975))

  # now choice
  pp <- seqEffectsSummaryRA$pp
  ppBySubjectByPostN <- aggregate(choice~pattern*subjects*postn,pp,mean)
  ppBySubjectByPostN <- ppBySubjectByPostN[ppBySubjectByPostN$pattern=='RRR',]
  ppBySubjectByPostNQuantiles <- aggregate(choice~subjects, ppBySubjectByPostN, quantile, c(0.025, 0.5, 0.975))

  # plot rt difference
  if(grepl('rt', orderBy)) {
    ordering1 <- order(dataBySubByChoiceW$rtDifference)
  } else {
    ordering1 <- order(dataBySub$choice)
  }
  ppBySubByChoiceQuantiles <- ppBySubByChoiceQuantiles[ordering1,]
  dataBySubByChoiceW <- dataBySubByChoiceW[ordering1,]

  def.mar <- par()$mar
  new.mar <- def.mar
  new.mar[1] <- 1
  par(mar=new.mar)
  ylims <- range(ppBySubByChoiceQuantiles$rtDifference)
  if(ylims[1] > 0) ylims[1] <- 0
  plot(1:nrow(dataBySubByChoiceW), dataBySubByChoiceW$rtDifference, ylim=ylims, #xlab='Subject', ylab='RT advantage (repeat - alternate) after 3 repeats')
       xlab=xlab1, ylab=ylab1, main=main, pch=4, xaxt='n', lwd=2)
  mtext(ylab1, side=2, line=2.5, cex=.66)
  abline(h=0, lty=2, col=1) #'grey')
  arrows(1:nrow(dataBySubByChoiceW),
         y0=ppBySubByChoiceQuantiles$rtDifference[,1], y1=ppBySubByChoiceQuantiles$rtDifference[,3], angle=90, code=3, length=0.025,
         col=2, lwd=2)
  lines(1:nrow(dataBySubByChoiceW), ppBySubByChoiceQuantiles$rtDifference[,2], col=2)


  # plot choices
  if(orderBy == 'rt') {
    ordering2 <- ordering1
  } else {
    ordering2 <- order(dataBySub$choice)
  }
  ppBySubjectByPostNQuantiles <- ppBySubjectByPostNQuantiles[ordering2,]
  dataBySub <- dataBySub[ordering2,]

  new.mar[1] <- def.mar[1]
  new.mar[3] <- 1
  par(mar=new.mar)
  ylims <- range(ppBySubjectByPostNQuantiles$choice)
  plot(1:nrow(dataBySub), dataBySub$choice, pch=4, ylim=ylims, # xlab='Subject', ylab='Proportion repeat choices after 3 repeats')
       xlab=xlab2, ylab=ylab2, xaxt='n', lwd=2)
  mtext(ylab2, side=2, line=2.5, cex=.66)
  abline(h=0.5, lty=2, col=1) #'grey')
  if(ylims[1] > 0.7) abline(h=0.75, lty=2, col=1)# 'grey')
  arrows(1:nrow(dataBySub),
         y0=ppBySubjectByPostNQuantiles$choice[,1], y1=ppBySubjectByPostNQuantiles$choice[,3], angle=90, code=3, length=0.025, lwd=2, col=2)
  lines(1:nrow(dataBySub), ppBySubjectByPostNQuantiles$choice[,2], col=2)
}


getPostPreDifference <- function(dat, rtcolname='rt', acccolname='accuracy', trialRange=c(-3,7)) {
  errors <- dat$accuracy == 0
  ## remove first trial (per participant) and last trial (per participant)
  trialRangeBySub <- aggregate(trials~subjects, dat, range)

  dat$exclude_trial <- FALSE
  for(subject in unique(trialRangeBySub$subjects)) {
    ## don't include errors on the first and last trial; can't calculate post-error slowing because of lack of a pre / post error trial
    dat[dat$subjects==subject,'exclude_trial'] = dat[dat$subjects==subject,'trials'] %in% trialRangeBySub[trialRangeBySub$subjects==subject,'trials']
  }
  errors <- errors & !dat$exclude_trial
  errors <- which(errors)

  # RT difference post-pre
  rtdf <- expand.grid(trialNposterror=1, probs='mean', measure='rt', value=NA)
  rtdf[rtdf$trialNposterror==1,'value'] <- mean(dat[errors+1,rtcolname] - dat[errors-1,rtcolname])

  # Absolute RTs
  rtabsdf <- expand.grid(trialNposterror=seq(trialRange[1],trialRange[2]), probs='mean', measure='rtabs', value=NA)
  for(trialNposterror in seq(trialRange[1],trialRange[2])) {
    trialIdx <- errors+trialNposterror
    trialIdx <- trialIdx[trialIdx>0]  ## remove potential trial numbers below 0
    if(trialNposterror != 0) trialIdx <- trialIdx[!trialIdx %in% errors]  ## only inspect accurate responses, dont include errors in pre/post error trials
    rtabsdf[rtabsdf$trialNposterror==trialNposterror,'value'] <- mean(dat[trialIdx,rtcolname], na.rm=TRUE)
  }

  return(rbind(rtabsdf, rtdf))
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
  # dat <- EMC2:::add_trials(dat)

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
      tmp <- postpreData[postpreData$probs=='mean'&postpreData$measure=='rt'&postpreData$trialNposterror==1,]
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


plotPES2 <- function(distributions=NULL, dat=NULL, pp=NULL, pp2=NULL, do.plot=TRUE, shift=FALSE,
                     orderSubjects=TRUE, firstTrialOnly=FALSE, use_mean=TRUE,
                     meanOnly=FALSE,
                     xlab1='Trial relative to error', xlab2='Trial relative to error',
                     rtmeasure='rt', accmeasure='accuracy', # or rtabs / accabs
                     main1='', main2='', ylim1=NULL, ylim2=NULL, setPar=TRUE,
                     plotAxisTickLabels=TRUE,
                     ylab1='', ylab2='', plotAccuracy=TRUE, add.legend=FALSE, plotpp=TRUE) {
  # dat <- EMC2:::add_trials(dat)

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

  if(do.plot) plot(xs, xs, type='n', xlab=xlab1, ylab=ylab1,  xaxt='n', ylim=ylim1, main=main1, xlim=range(xs)+c(-.25, .25))
  #  grid()
  #  abline(h=0, lty=2, col='grey')
  axis(side=1, at=xs)
  abline(h=mean(aggregate(rt~subjects,dat,mean)[,2]), col=par()$col, lty=2)
  abline(v=0, lty=2, col=par()$col)

  ## data
  points(postpreDataAverage[postpreDataAverage$measure==rtmeasure,'trialNposterror'],
         postpreDataAverage[postpreDataAverage$measure==rtmeasure,'value'], pch=4, lwd=2)

  ## posterior predictives
  arrows(x0=postprePPAverage[postprePPAverage$measure==rtmeasure,'trialNposterror']-.15*(!do.plot)+shift*0.15,
         y0=postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,1],
         y1=postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,3],
         angle=90, code=3, length=.02, lwd=2, col=2+(!do.plot))
  xs <- postprePPAverage[postprePPAverage$measure==rtmeasure,'trialNposterror']
  lines(xs, postprePPAverage[postprePPAverage$measure==rtmeasure,'value'][,2], col=2+(!do.plot))
}

getPostPreDifferenceBySATCondition <- function(dat, rtcolname='rt', acccolname='accuracy', conditionColumn='E') {
  errors <- dat$accuracy == 0
  ## remove first trial (per participant) and last trial (per participant)
  if(!'trials' %in% colnames(dat)) dat <- EMC2:::add_trials(dat)
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
get_PES_measures <- function(model, ...) {
  opts <- list(...)
  if('trendModel' %in% names(opts)) {
    wagenmakers2004_CS <- loadDataSamplersPP('wagenmakers2004_CS', model, trendModel='DCT', trendPar = 'B', nTrendPars = 3)
    forstmann2008 <- loadDataSamplersPP('forstmann2008', model, trendModel='DCT', trendPar = 'v', nTrendPars = 3)
    mileticvanmaanen2019exp2block2 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', model, trendModel='DCT', trendPar = 'v', nTrendPars = 3)
    wagenmakers2008exp2 <- loadDataSamplersPP('wagenmakers2008exp2', model, trendModel='DCT', trendPar = 'B', nTrendPars = 3)
  } else {
    wagenmakers2004_CS <- loadDataSamplersPP('wagenmakers2004_CS', model)
    forstmann2008 <- loadDataSamplersPP('forstmann2008', model)
    mileticvanmaanen2019exp2block2 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', model)
    wagenmakers2008exp2 <- loadDataSamplersPP('wagenmakers2008exp2', model)
  }
  # filter these to get rid of subjects with very few errors
  dat <- mileticvanmaanen2019exp2block2$dat
  dat$error <- !dat$accuracy
  nerrors <- aggregate(error~subjects,dat,sum)
  goodsubs <- nerrors[nerrors[,2]>40,1]
  dat <- dat[dat$subjects %in% goodsubs,]
  mileticvanmaanen2019exp2block2$dat <- mileticvanmaanen2019exp2block2$dat[mileticvanmaanen2019exp2block2$dat$subjects %in% goodsubs,]
  mileticvanmaanen2019exp2block2$pp <- mileticvanmaanen2019exp2block2$pp[mileticvanmaanen2019exp2block2$pp$subjects %in% goodsubs,]

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


plot_4_tasks <- function(tasks, learningModel, trendModel='NULL', nTrendPars=3, trial_durations,
                         debug=FALSE, full_x=FALSE, margin_text='', margin_text2='', plot_xticklabels=TRUE, ...) {
  opts <- list(...)

  for(taskn in 1:length(tasks)) {
    task <- tasks[taskn]
    dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))

    if(trendModel=='winning') {
      if(task %in% c('wagenmakers2004_CS', 'wagenmakers2008exp2')) {
        trendPar <- 'B'
      } else {
        trendPar <- 'v'
      }
      samplers_fn <- loadDataSamplersPP(task, learningModel, trendModel='DCT', trendPar=trendPar, nTrendPars = 3, return_fns = TRUE)
    } else {
      samplers_fn <- loadDataSamplersPP(task, learningModel, trendModel='NULL', trendPar=trendPar, nTrendPars = 3, return_fns = TRUE)
    }
    pp_fn <- samplers_fn$pp_fn
    samplers_fn <- samplers_fn$samplers_fn
    # pp_fn <- gsub('.RData', '_pp-unconditional.RData', gsub('/samples', '/posteriorpredictives', samplers_fn))
    # pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
    if(!file.exists(pp_fn)) { plot.new(); next}

    pp <- EMC2:::loadRData(pp_fn)
    # samplers <- EMC2:::loadRData(samplers_fn)
    if(debug) pp <- NULL

    plotSpectrum(dat, pp=pp, pp2=NULL,
                 dat_fn=paste0('./datasets/', task, '.RData'),
                 pp_fn=pp_fn, pp2_fn=NULL,
                 xlab='', ylab='', main='', plot.log=TRUE, detrend=FALSE, trial_duration = trial_durations[taskn], full_x = full_x,
                 plot_xticklabels=plot_xticklabels)
    if(taskn == 1) mtext('Log power', side=2, line=2, cex=.66, las=0)
    if(taskn == 1) mtext(text=margin_text, side=2, line=4, las=0, cex=par()$cex*par()$cex.main, font=2)
    # if(taskn == 1) mtext(text=margin_text2, side=2, line=5, las=0, cex=par()$cex*par()$cex.main, font=2)
    if('xlab' %in% names(opts)) mtext(opts$xlab, side=1, line=3, cex=.66, las=0)
    if('mains' %in% names(opts)) mtext(opts$mains[taskn], side=3, cex=par()$cex*par()$cex.main, font=2)
  }
}



# Hard-easy ---------------------------------------------------------------

applyQuintiles <- function(x) {
  cut(x, breaks=c(quantile(x, probs = seq(0, 1, by = 0.20), na.rm=TRUE)),
      labels=c("5","4","3","2","1"), include.lowest=TRUE)
}

get_aggregates_single_subject <- function(subject, samplers, pp, npost=100) {
  dadm <- samplers[[1]]$data[[subject]]

  out = do.call(rbind, lapply(1:npost, function(x) {
    npars <- attr(pp, 'npars')
    npars <- npars[npars$subjects==subject&npars$postn==x,]
    q <- npars[,'q03']

    ## merge data and Q-value of this postn
    tmp <- cbind(dadm, q)
    tmp$qbin <- as.numeric(applyQuintiles(q))  # bin Q-values
    tmp$qbin <- 6-tmp$qbin  # flip to have increasing difficulty instead of decreasing


    ## merge data and pp
    data_pp <- merge(x=tmp, y=pp[pp$postn==x&pp$subjects==subject,],
                     by.x=c('subjects', 'S', 'trials', 'stim'),
                     by.y=c('subjects', 'S', 'trials', 'stim'),
                     suffixes = c('_data', '_pp'))
    data_pp <- data_pp[order(data_pp$trials),]
    data_pp <- data_pp[data_pp$winner==TRUE,]
    data_pp$qbin_prev <- Lag(data_pp$qbin, 1)

    ## aggregate data and pp
    agg <- aggregate(cbind(rt_data, rt_pp)~qbin_prev, data_pp, mean)
    agg$subject <- subject
    agg$postn <- x
    return(agg)
  }))
}

get_da_with_Q3bins <- function(samplers, pp) {
  aggregates <- mclapply(names(samplers[[1]]$data), #1:length(samplers[[1]]$data),
                         function(subject) get_aggregates_single_subject(subject=subject, samplers=samplers, pp=pp),
                         mc.cores = 20)

  return(do.call(rbind, aggregates))
}

get_slope <- function(rt_data, fip_difficulty=TRUE) {
  if(flip_difficulty) rt_data$D <- max(rt_data$D)-rt_data$D+1
  lm(rt~as.numeric(D), dat=rt_data)$coef[2]
}

get_stats <- function(pp2, da, use_median=FALSE, flip=FALSE) {
  ## Calculate quantiles for (rt, accuracy) x (data, pp) to get credible intervals
  ## assumes column name is 'D'
  if(!'D' %in% colnames(da)) da$D <- da$q3p5prevbin
  combined_data <- merge(pp2, da, by=c('trials', 'subjects', 'S'))
  if(use_median) {
    rt_pp <- aggregate(rt.x~D*postn*subjects, combined_data, quantile, c(.1, .5, .9)) # RT quantiles per trial
    rt_pp <- aggregate(rt.x~D*subjects, rt_pp, quantile, c(.025, .5, .975))
    rt_pp <- aggregate(cbind(`10%`, `50%`, `90%`)~D, rt_pp, mean)

    ## data
    rt_data <- aggregate(rt~D*subjects, da, quantile, c(0.1, .5, .9))
    rt_data <- aggregate(rt~D, rt_data, mean)
  } else {
    ## data
    rt_data <- aggregate(rt~D*subjects,da, mean)

    # get slope by subject
    slopes_data <- sapply(unique(rt_data$subjects), function(subject) lm(rt~as.numeric(D), dat=rt_data[rt_data$subjects==subject,])$coef[2])

    rt_data <- aggregate(rt~D, rt_data, mean)

    if(!flip) {
      rt_pp <- aggregate(rt.x~D*postn*subjects, combined_data, mean)  # mean RT per trial
      rt_pp <- aggregate(rt.x~D*subjects, rt_pp, quantile, c(.025, .5, .975))
      rt_pp <- aggregate(rt.x~D, rt_pp, mean)
    } else {
      rt_pp <- aggregate(rt.x~D*postn*subjects, combined_data, mean)  # mean RT per trial
      # slopes x subject
      slopes_pp <- sapply(unique(rt_pp$subjects), function(subject) {
        sapply(unique(rt_pp$postn), function(postn) {
          lm(rt.x~as.numeric(D), dat=rt_pp[rt_pp$subjects==subject&rt_pp$postn==postn,])$coef[2]
        })
      })

      rt_pp <- aggregate(rt.x~D*postn, rt_pp, mean)    # mean across subjects first
      rt_pp <- aggregate(rt.x~D, rt_pp, quantile, c(.025, .5, .975))   # then take the credible interval across iterations
      rt_pp[,2:4] <- rt_pp$rt.x
    }
  }


  return(list(rt_pp=rt_pp, #acc_pp=acc_pp,
              rt_data=rt_data,
              slopes_data=slopes_data,
              slopes_pp=slopes_pp,
              slopes_ci=apply(slopes_pp,2,quantile,c(0.025, 0.975)))) #, acc_data=acc_data))
}


plot_q3_effect <- function(stats, mrt=0, main1='',ylab1=ylab, xlab1='', xlab2='', plotAccuracy=FALSE) {
  rt_pp <- stats[['rt_pp']]
  acc_pp <- stats[['acc_pp']]
  rt_data <- stats[['rt_data']]
  acc_data <- stats[['acc_data']]

  nc <- ncol(rt_pp)
  if(nc == 4) {
    rt_pp <- cbind(rt_pp, rt_pp[,2:4])
    rt_data <- cbind(rt_data, rt_data[,2])
  }
  ylim1 <- range(rt_pp[,5:7])*c(.95, 1.05)
  ylim2 <- c(.55, .95) #range(range(accstats[,2:ncol(accstats)]))

  xs <- 1:nrow(rt_pp)
  plot(xs, rt_data[,3], xlab=xlab1, ylab=ylab1, pch=4, ylim=ylim1, main=main1, lwd=2, xaxt='n', type='n')
  points(xs, rt_data[,3], pch=4, lwd=2)
  abline(h=mrt, col=par()$col, lty=2)
  arrows(xs, y0=rt_pp[,5], y1=rt_pp[,7], code=3, angle=90, length=.05, col=2, lwd=2)
  lines(xs, y=rt_pp[,6], col=2, lwd=1)

  if(plotAccuracy) {
    plot(xs, acc_data[,2], xlab=xlab2, ylab='Accuracy current trial', pch=4,
         ylim=ylim2,
         lwd=2, xaxt='n')
    axis(side=1, at=xs, labels=rt_pp$D)

    arrows(xs, y0=acc_pp[,2], y1=acc_pp[,4], code=3, angle=90, length=.05, col=2, lwd=2)
    lines(xs, y=acc_pp[,3], col=2, lwd=1)
  }
}

get_difficulty_effects <- function(dat, pp) {
  # aggregate posterior predictives: Get mean RT per previous difficulty per posterior predictive
  aggBySubjectByPostN <- aggregate(rt~difficulty_prev*subjects*postn,pp,mean)
  aggByPostN <- aggregate(rt~difficulty_prev*postn, aggBySubjectByPostN, mean)

  # and 2.5th and 97.5th quantile across posterior predictives to get uncertainty, and 0.5 for median
  aggPp <- aggregate(rt~difficulty_prev, aggByPostN, quantile, c(0.025, 0.5, 0.975))

  # get uncertainty per subject
  aggBySubjectPp <- aggregate(rt~difficulty_prev*subjects, aggBySubjectByPostN, quantile, c(0.025, 0.5, 0.975))

  # For data, do the same, just get mean
  aggBySubjectData <- aggregate(rt~difficulty_prev*subjects,dat,mean)
  aggData <- aggregate(rt~difficulty_prev, aggBySubjectData, mean)

  # Slopes
  dat$difficulty_prev <- as.numeric(dat$difficulty_prev)
  coefficientsBySubject <- sapply(sort(unique(dat$subjects)), function(subject) lm(rt~difficulty_prev, dat[dat$subjects==subject,])$coef)
  coefficientsBySubjectData <- reshape2::melt(coefficientsBySubject, value.name='coefficient', varnames=c('parameter', 'subjects'))

  # slopes of pp
  pp$difficulty_prev <- as.numeric(pp$difficulty_prev)
  coefficientsBySubjectByPostN <- parallel::mclapply(1:100, function(postn) {
    coeffs_ <- sapply(sort(unique(dat$subjects)), function(subject) lm(rt~difficulty_prev, pp[pp$subjects==subject&pp$postn==postn,])$coef)
    coeffs_ <- reshape2::melt(coeffs_, value.name='coefficient', varnames=c('parameter', 'subjects'))
    coeffs_$postn <- postn
    return(coeffs_)
  }, mc.cores=20)
  coefficientsBySubjectByPostN <- do.call(rbind, coefficientsBySubjectByPostN)

  ## aggregate slopes as well
  coefficientsBySubjectPp <- aggregate(coefficient~subjects*parameter, coefficientsBySubjectByPostN, quantile, c(0.025, 0.5, 0.975))

  ## and mean slope?

  return(list(## Data aggregates: Mean across subjects, mean per subject
    aggData=aggData,
    aggBySubjectData=aggBySubjectData,

    ## Pp aggregates: median+CI across group
    aggPp=aggPp,
    aggBySubjectPp=aggBySubjectPp,  # median+CI per subject

    # slopes
    coefficientsBySubjectData=coefficientsBySubjectData,
    coefficientsBySubjectByPostN=coefficientsBySubjectByPostN,
    coefficientsBySubjectPp=coefficientsBySubjectPp  # CI
  ))
}

prep_data <- function(dat) {
  if(!'trials' %in% colnames(dat)) dat <- EMC2:::add_trials(dat)
  if('difficulty' %in% colnames(dat)) {
    dat$difficulty <- as.numeric(dat$difficulty)
    dat$difficulty <- 6-dat$difficulty ## flip difficulty for plotting
    dat$difficulty_prev <- Lag(dat$difficulty, 1)
  } else if('W' %in% colnames(dat)) {
    dat$difficulty <- as.numeric(factor(dat$W, levels=c('hf', 'lf', 'vlf', 'nonword')))
    dat$difficulty_prev <- Lag(dat$difficulty, 1)
    dat <- dat[dat$difficulty_prev < 4,]  # exclude non-words; difficulty is undetermined here
  }
  dat <- dat[dat$trials>1,]  # remove first trial, no previous difficulty

  dat
}


# Plot 1. Group level rt ~ difficulty on previous trial
plot_group_effect <- function(stats, mrt=NULL, ...) {
  ylim <- range(c(stats$aggData$rt, stats$aggPp$rt))
  plot(as.integer(stats$aggData$difficulty_prev), stats$aggData$rt, pch=4, lwd=2, xlab='Difficulty', ylab='RT (s)', ylim=ylim, ...)
  if(!is.null(mrt)) abline(h=mrt, lty=2, col=1)

  arrows(x0=stats$aggPp$difficulty_prev, y0=stats$aggPp$rt[,1], y1=stats$aggPp$rt[,3], code=3, angle=90, length=0.025, col=2, lwd=2)
  lines(x=stats$aggPp$difficulty_prev, stats$aggPp$rt[,2], lwd=2, col=2)
}

# Plot 2. subject level slopes
plot_subject_slopes <- function(stats, ...) {
  coefficients <- stats$coefficientsBySubjectData[stats$coefficientsBySubjectData$parameter=='difficulty_prev',]
  ordering <- order(coefficients$coefficient)
  coefficientsPp <- stats$coefficientsBySubjectPp[stats$coefficientsBySubjectPp$parameter=='difficulty_prev',]

  ylim <- range(c(coefficients$coefficient, coefficientsPp$coefficient))
  plot(1:nrow(coefficients), coefficients[ordering, 'coefficient'], pch=4, ylab='Slope', xlab='Participant', lwd=2, ...)
  abline(h=0)

  arrows(x0=1:nrow(coefficients),
         y0=coefficientsPp$coefficient[ordering,1],
         y1=coefficientsPp$coefficient[ordering,3],
         code=3, angle=90, length=0.025, col=2, lwd=2
  )
  lines(1:nrow(coefficients), coefficientsPp$coefficient[ordering,2], lwd=2,col=2)
}

plot_subjective_group_effect <- function(aggregates, xlab='Difficulty', ylab='RT (s)', mrt=NULL, ...) {
  group <- aggregate(rt_data~qbin_prev*subject, aggregates, mean)
  group <- aggregate(rt_data~qbin_prev, group, mean)
  ci <- aggregate(rt_pp~qbin_prev*postn, aggregates, mean)
  ci <- aggregate(rt_pp~qbin_prev, ci, quantile, c(0.025, 0.5, 0.975))
  ylim <- range(c(group$rt_data, ci$rt_pp))
  plot(as.numeric(group$qbin_prev), group$rt_data, pch=4, lwd=2, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  if(!is.null(mrt)) abline(h=mrt, lty=2, col=1)
  arrows(x0=as.numeric(ci$qbin_prev), y0=ci$rt_pp[,1],y1=ci$rt_pp[,2], angle=90, code=3, length=0.025, col=2,lwd=2)
}



# First get all plotting data (except Fourier) ---------------------------------------------
learningModel <- 'uAH'

for(learningModel in c('zSM', 'uAH','bV', 'NULL')) {
  # if(learningModel == '')
  if(learningModel == 'NULL') {
    PESmeasures <- get_PES_measures(learningModel, trendModel='DCT')
    wagenmakers2004_CS <- loadDataSamplersPP('wagenmakers2004_CS', learningModel, trendModel='DCT', trendPar='B', nTrendPars = 3)
    forstmann2008 <- loadDataSamplersPP('forstmann2008', learningModel, trendModel='DCT', trendPar='v', nTrendPars = 3)
    mileticvanmaanen2019exp2block2 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', learningModel, trendModel='DCT', trendPar='v', nTrendPars = 3)
    wagenmakers2008exp2 <- loadDataSamplersPP('wagenmakers2008exp2', learningModel, trendModel='DCT', trendPar='B', nTrendPars = 3)
  } else {
    PESmeasures <- get_PES_measures(learningModel)
    wagenmakers2004_CS <- loadDataSamplersPP('wagenmakers2004_CS', learningModel)
    forstmann2008 <- loadDataSamplersPP('forstmann2008', learningModel)
    mileticvanmaanen2019exp2block2 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', learningModel)
    wagenmakers2008exp2 <- loadDataSamplersPP('wagenmakers2008exp2', learningModel)
  }

  # Objective difficulties: hard-easy
  dat4 <- prep_data(wagenmakers2008exp2$dat)
  pp4 <- prep_data(wagenmakers2008exp2$pp)
  stats4 <- get_difficulty_effects(dat4, pp4)


  # By choice
  seqEffectsSummary1 <- getSequentialEffectsBySub(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, repeatAlternate = FALSE)
  seqEffectsSummary2 <- getSequentialEffectsBySub(dat=forstmann2008$dat, pp=forstmann2008$pp, repeatAlternate = FALSE)
  seqEffectsSummary3 <- getSequentialEffectsBySub(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, repeatAlternate = FALSE)
  seqEffectsSummary4 <- getSequentialEffectsBySub(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, repeatAlternate = FALSE)

  # error related

  ## big-ass layout
  m <- matrix(nrow=100, ncol=100)
  row1 <- 1:20
  row2 <- 21:45
  row3 <- 46:70
  row4 <- 71:100

  # Fourier, stimulus history
  column1Upper <- 1:25
  column2Upper <- 26:50
  column3Upper <- 51:75
  column4Upper <- 76:100
  m[row1,column1Upper] <- 1
  m[row1,column2Upper] <- 2
  m[row1,column3Upper] <- 3
  m[row1,column4Upper] <- 4

  # stimulus-history: fill column-wise
  m[row2,column1Upper] <- 5
  m[row3,column1Upper] <- 6
  m[row2,column2Upper] <- 7
  m[row3,column2Upper] <- 8
  m[row2,column3Upper] <- 9
  m[row3,column3Upper] <- 10
  m[row2,column4Upper] <- 11
  m[row3,column4Upper] <- 12

  # PES
  column1lower <- 1:16
  column2lower <- 17:32
  column3lower <- 33:49
  column4lower <- 50:66
  column5lower <- 67:83
  column6lower <- 84:100

  ## bottom row
  m[row4,column1lower] <- 13
  m[row4,column2lower] <- 14
  m[row4,column3lower] <- 15
  m[row4,column4lower] <- 16
  m[row4,column5lower] <- 17
  m[row4,column6lower] <- 18

  # layout.show(layout(m))

  #  Fourier ----------------------------------------------------------------
  tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
  trial_durations <- c(1.265, 2.850, 1.627, 0.837)
  debug = FALSE

  for(ftype in c('pdf','jpeg')) {
    if(ftype == 'pdf') pdf(file=paste0('./figures/spectra_SH_ER_', ifelse(learningModel=='NULL', 'DCT', learningModel), 'only.pdf'), width=8, height=7)
    if(ftype == 'jpeg') jpeg(file=paste0('./figures/spectra_SH_ER_', ifelse(learningModel=='NULL', 'DCT', learningModel), 'only.jpeg'), width=8, height=7, units='in', quality=100, res=500)
    layout(m)
    par(bty='l', oma=c(0,2,0,0), mgp=c(2,1,0))
    par(mar=c(1,2,2,1))
    plot_4_tasks(tasks, learningModel = learningModel, trendModel=ifelse(learningModel=='NULL', 'winning', 'NULL'), trendPar='NULL', trial_durations = trial_durations, debug=debug, full_x=TRUE,
                 mains=c('Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4'), xlab='Period (s)', margin_text='', plot_xticklabels=TRUE, margin_text2='')


    # stimulus history --------------------------------------------------------
    # Plot
    plotSequentialEffectsGroup(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, seqEffectsSummary=seqEffectsSummary1,
                               nHistory = 3, main1='', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
    plotSequentialEffectsGroup(dat=forstmann2008$dat, pp=forstmann2008$pp, seqEffectsSummary=seqEffectsSummary2,
                               nHistory = 3, main1='', full.legend=FALSE)
    plotSequentialEffectsGroup(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, seqEffectsSummary=seqEffectsSummary3,
                               nHistory = 3, main1='', full.legend=FALSE)
    plotSequentialEffectsGroup(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, seqEffectsSummary=seqEffectsSummary4,
                               nHistory = 3, main1='', full.legend=FALSE)


    # Error-related -----------------------------------------------------------
    ## Plot
    par(mar=c(3,2,1,1))
    plotPES2(distributions=PESmeasures$distributions1, dat=PESmeasures$wagenmakers2004_CS$dat, pp2=NULL, main1='Ds 1', rtmeasure = 'rtabs', accmeasure='accabs', ylab1='RT (s)', ylab='Accuracy', plotAccuracy = FALSE, meanOnly=TRUE, setPar = FALSE)
    mtext('RT (s)', side=2, cex=.66, line=2)
    plotPES2(distributions=PESmeasures$distributions2bSPD, dat=PESmeasures$forstmann2008$dat, pp2=NULL, main1='Ds 2 (spd)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
    plotPES2(distributions=PESmeasures$distributions2bACC, dat=PESmeasures$forstmann2008$dat, pp2=NULL, main1='Ds 2 (acc/neu)', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE, ylim1=c(0.35, 0.54))
    plotPES2(distributions=PESmeasures$distributions3, dat=PESmeasures$mileticvanmaanen2019exp2block2$dat, pp2=NULL, main1='Ds 3', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)
    plotPES2(distributions=PESmeasures$distributions4, dat=PESmeasures$wagenmakers2008exp2$dat, pp2=NULL, main1='Ds 4', rtmeasure = 'rtabs', accmeasure='accabs', plotAccuracy = FALSE,meanOnly=TRUE, setPar = FALSE)


    # Hard-easy ---------------------------------------------------------------
    plot_group_effect(stats4, mrt=mean(wagenmakers2008exp2$dat$rt), main='Ds 4 Hard-Easy', xaxt='n')
    axis(side=1, at=1:3, labels=c('Easy', '', "Hard"))
    dev.off()
  }
}
