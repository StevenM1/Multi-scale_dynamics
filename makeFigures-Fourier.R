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


get_fns <- function(task, learningModel, trendModel='DCT', trendPar='B', nTrendPars=3, samples_dir='./samples') {
  samplers_fn <- generate_filenames(dataset=task, learningModel=learningModel,
                                    trendModel=trendModel, trendPar=trendPar,
                                    nTrendPars=nTrendPars,
                                    samples_dir = samples_dir)
  pp_fn <- gsub('.RData', '_pp-unconditional.RData', samplers_fn)
  # if(trendModel=='NULL') {
  #   pp_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '_trend-NULL_pp-unconditional.RData')
  #   samplers_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '_trend-NULL.RData')
  # } else {
  #   pp_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '-trend-', trendModel, '-', trendPar, '-', nTrendPars, '_pp.RData')
  #   samplers_fn <- paste0(samples_dir, task, '_model-RDM-', learningModel, '-trend-', trendModel, '-', trendPar, '-', nTrendPars, '.RData')
  # }
  return(list(pp_fn=pp_fn, samplers_fn=samplers_fn))
}

get_model_comparison <- function(fns, subfilter=0) {
  allSamplers <- list()
  allGds <- c()
  for(fn in fns) {
    load(fn)
    modelName <- strsplit(strsplit(fn, 'model-')[[1]][2], '.RData')[[1]][1]
    nchains <- chain_n(samplers)
    if(nchains[1,4]>=100) {
      allSamplers[[modelName]] <- samplers
      allGds <- c(allGds, max(EMC2:::gd_pmwg(samplers, subfilter=subfilter, print_summary = FALSE)))
    }
  }
  if(length(allSamplers)>0) {
    mc = EMC2:::compare_IC(allSamplers, print_summary=FALSE, subfilter=subfilter)
    mc <- cbind(mc, max.gd=allGds)
    print(round(mc,3))
    mc
  }
}
# getFns <- function(task, samples_dir='./samples/') {
#   all_fns <- Sys.glob(paste0(samples_dir, task, '*.RData'))
#   fns = all_fns[!grepl('_pp', all_fns)]
#   return(fns)
# }

get_delta_IC <- function(mc) {
  mc$delta_DIC <- mc$DIC-min(mc$DIC)
  mc$delta_BPIC <- mc$BPIC-min(mc$BPIC)
  return(mc)
}



## Plotting wrappers
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
    samplers_fn <- samplers_fn$samplers_fn    # pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
    if(!file.exists(pp_fn)) { plot.new(); next}

    pp <- EMC2:::loadRData(pp_fn)
    # samplers <- EMC2:::loadRData(samplers_fn)
    if(debug) pp <- NULL

    plotSpectrum(dat, pp=pp, pp2=NULL,
                 dat_fn=paste0('./datasets/', task, '.RData'),
                 pp_fn=pp_fn, pp2_fn=NULL,
                 xlab='', ylab='', main='', plot.log=TRUE, detrend=FALSE, trial_duration = trial_durations[taskn], full_x = full_x,
                 plot_xticklabels=plot_xticklabels)
    if(taskn == 1) mtext('Log power', side=2, line=3, cex=.66, las=0)
    if(taskn == 1) mtext(text=margin_text, side=2, line=4, las=0, cex=par()$cex*par()$cex.main, font=2)
    # if(taskn == 1) mtext(text=margin_text2, side=2, line=5, las=0, cex=par()$cex*par()$cex.main, font=2)
    if('xlab' %in% names(opts)) mtext(opts$xlab, side=1, line=3, cex=.66, las=0)
    if('mains' %in% names(opts)) mtext(opts$mains[taskn], side=3, cex=par()$cex*par()$cex.main, font=2)
  }
}

plot_4_tasks_cdfs <- function(tasks, margin_text='') {
  for(taskn in 1:length(tasks)) {
    task <- tasks[taskn]
    dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
    fns <- loadDataSamplersPP(task, learningModel='NULL', trendModel='NULL', return_fns=TRUE)
    pp_fn <- fns$pp_fn; samplers_fn <- fns$samplers_fn
    pp <- EMC2:::loadRData(pp_fn)
    samplers <- EMC2:::loadRData(samplers_fn)

    dat$R <- factor(as.logical(dat$accuracy), levels=c(TRUE, FALSE), labels=c('Correct', 'Error'))
    pp$R <- factor(as.logical(dat$accuracy), levels=c(TRUE, FALSE), labels=c('Correct', 'Error'))
    dat$S <- factor('', levels=c(''))
    pp$S <- factor('', levels=c(''))

    #par(mfrow=c(2,1),bty='l')
    main_ <- paste0('Dataset ', taskn)
    if(taskn==1) lpos = 'right' else lpos=NULL
    # par(mar=c(2.5,3,2,1))
    plot_fit(data=dat, pp=pp, factors='S', layout=NULL, main = main_, lpos=lpos, qp_cex=2, pqp_cex=1)
    if(taskn == 1) mtext('p(R)', side=2, line=2.5, cex=.66,las=0)
    if(taskn == 1) mtext('A. IID-RDM', side=2, line=4, cex=par()$cex*par()$cex.main, las=0, font=2)
    # if(taskn == 1) mtext(margin_text, side=2, line=4, cex=par()$cex*par()$cex.main, las=0, font=2)
  }
}


# New Figure 1: Only CDFs and Fourier Spectra of IID-RDM -----------------------------
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
trial_durations <- c(1.265, 2.850, 1.627, 0.837)
debug = FALSE

ftype <- 'jpeg'
if(ftype == 'pdf') pdf(file='./figures/spectra_IIDRDM.pdf', width=8, height=4)
if(ftype == 'jpeg') jpeg(file='./figures/spectra_IIDRDM.jpeg', width=8, height=4, units='in', quality=100, res=200)
#layout(matrix(1:24, byrow=TRUE), heights=c(1,))
par(mfrow=c(2, length(tasks)), mar=c(3,3,2,0.5), oma=c(1,2,0,0), bty='l', mgp=c(2,1,0), las=1)

# Row 1: CDFs
plot_4_tasks_cdfs(tasks)

par(mar=c(3,3,0,0.5))
# Row 2: Null
plot_4_tasks(tasks, learningModel = 'NULL', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('', '', '', ''), xlab='Period (s)', margin_text='B. IID-RDM', plot_xticklabels=TRUE, margin_text2='')
dev.off()


# New Figure 3: Spectra of MS1-RDM, MS2-RDM, MS3-RDM ----------------------
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
trial_durations <- c(1.265, 2.850, 1.627, 0.837)
debug = FALSE

ftype <- 'jpeg'
if(ftype == 'pdf') pdf(file='./figures/spectra_MS1-MS2-MS3-RDMs.pdf', width=8, height=4)
if(ftype == 'jpeg') jpeg(file='./figures/spectra_MS1-MS2-MS3-RDMs.jpeg', width=8, height=4, units='in', quality=100, res=200)
#layout(matrix(1:24, byrow=TRUE), heights=c(1,))
par(mfrow=c(3, length(tasks)), mar=c(3,3,2,0.5), oma=c(1,2,0,0), bty='l', mgp=c(2,1,0), las=1)
par(mar=c(3,3,1,0.5))
# Row 1: zSM
plot_4_tasks(tasks, learningModel = 'zSM', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4'), xlab='', margin_text='A. MS1-RDM', margin_text2='A',
             plot_xticklabels=FALSE)

# Row 4: zSMuAH
par(mar=c(3,3,0,0.5))
plot_4_tasks(tasks, learningModel = 'zSMuAH', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('', '', '', ''), xlab='', margin_text='B. MS2-RDM', plot_xticklabels=FALSE)

# Row 5: zSMuAHbV
plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('', '', '', ''), xlab='Period (s)', margin_text='C. MS3-RDM')

dev.off()


## MS3-RDM without and with trends, full range
ftype <- 'jpeg'
if(ftype == 'pdf') pdf(file='./figures/spectra_MS3RDM_notrends_vs_trend_fullrange.pdf', width=8, height=3)
if(ftype == 'jpeg') jpeg(file='./figures/spectra_MS3RDM_notrends_vs_trend_fullrange.jpeg', width=8, height=3, units='in', quality=100, res=200)

par(mfrow=c(2, length(tasks)), mar=c(3,2,2,1), oma=c(2,3,0,0), bty='l', mgp=c(2,1,0), las=1)
# Row 5: zSMuAHbV
plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4'), xlab='Period', full_x=TRUE, margin_text='A. MS3-RDM')

# Row 6: zSMuAHbV + trends
plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='winning', trendPar='winning', trial_durations = trial_durations, debug=debug,
             mains=c('', '', '', ''), xlab='Period', full_x=TRUE, margin_text='B. MS3-RDM + Trends')
dev.off()



# What about Dataset 5? (MVM exp 2 block 1)
pp <- EMC2:::loadRData('./posteriorpredictives_trends/mileticvanmaanen2019exp2block1_model-zSMuAHbV_trend-NULL_pp-unconditional.RData')
data <- EMC2:::loadRData('./datasets/mileticvanmaanen2019exp2block1.RData')

ftype <- 'jpeg'
if(ftype == 'pdf') pdf(file='./figures/spectra_MS3RDM_dataset5.pdf', width=6, height=3)
if(ftype == 'jpeg') jpeg(file='./figures/spectra_MS3RDM_dataset5.jpeg', width=6, height=3, units='in', quality=100, res=200)
par(mfrow=c(1, 2), mar=c(4,3,2,0.5), oma=c(1,2,0,0), bty='l', mgp=c(3,1,0), las=1)
data$R <- factor(as.logical(data$accuracy), levels=c(TRUE, FALSE), labels=c('Correct', 'Error'))
pp$R <- factor(as.logical(data$accuracy), levels=c(TRUE, FALSE), labels=c('Correct', 'Error'))
data$S <- factor('', levels=c(''))
pp$S <- factor('', levels=c(''))

#par(mfrow=c(2,1),bty='l')
main_ <- paste0('Dataset 5')
#if(taskn==1) lpos = 'right' else lpos=NULL
# par(mar=c(2.5,3,2,1))
plot_fit(data=data, pp=pp, factors='S', layout=NULL, main = 'A. Defective CDF', lpos='right')
mtext('p(R)', side=2, line=2.5, cex=par()$cex*par()$cex.lab,las=0)
par(mgp=c(3,1,0))
mar=c(4,3,2,0.5)
plotSpectrum(dat=data,pp=pp,
             dat_fn='./datasets/mileticvanmaanen2019exp2block1.RData',
             pp_fn='./posteriorpredictives_trends/mileticvanmaanen2019exp2block1_model-zSMuAHbV_trend-NULL_pp-unconditional.RData', trial_duration = mean(data$rt)+0.8,xlab = 'Period',ylab = 'Log power', main='B. Spectrum')
mtext('Log power', side=2, cex=par()$cex*par()$cex.lab, font=1, las=0, line=2.5)
#mtext(text='Spectrum', side=3, cex=par()$cex*par()$main, font=2)
dev.off()



# Spectra full ------------------------------------------------------------
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# trial_durations <- c(1.265, 2.850, 1.627, 0.837)
# debug = FALSE
#
#
# ftype <- 'jpeg'
# if(ftype == 'pdf') pdf(file='./figures/spectra_all.pdf', width=8, height=8)
# if(ftype == 'jpeg') jpeg(file='./figures/spectra_all.jpeg', width=8, height=8, units='in', quality=100, res=80)
# #layout(matrix(1:24, byrow=TRUE), heights=c(1,))
# par(mfrow=c(6, length(tasks)), mar=c(3,3,2,0.5), oma=c(1,2,0,0), bty='l', mgp=c(2,1,0), las=1)
#
# # Row 1: CDFs
# plot_4_tasks_cdfs(tasks)
#
# par(mar=c(3,3,0,0.5))
# # Row 2: Null
# plot_4_tasks(tasks, learningModel = 'NULL', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='', margin_text='B. IID-RDM', plot_xticklabels=FALSE, margin_text2='')
#
# # Row 3: zSM
# plot_4_tasks(tasks, learningModel = 'zSM', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='', margin_text='C. MS1-RDM', margin_text2='C',
#              plot_xticklabels=FALSE)
#
# # Row 4: zSMuAH
# plot_4_tasks(tasks, learningModel = 'zSMuAH', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='', margin_text='D. MS2-RDM', plot_xticklabels=FALSE)
#
# # Row 5: zSMuAHbV
# plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='', margin_text='E. MS3-RDM')
#
# # Row 6: zSMuAHbV + trends
# plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='winning', trendPar='winning', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='Period', full_x=TRUE, margin_text='F. MS3-RDM + Trends')
# dev.off()
#


# Remove last row
# ftype <- 'pdf'
# if(ftype == 'pdf') pdf(file='./figures/spectra_all_notrends.pdf', width=8, height=6.67)
# if(ftype == 'jpeg') jpeg(file='./figures/spectra_all_notrends.jpeg', width=8, height=6.67, units='in', quality=100, res=100)
# #layout(matrix(1:24, byrow=TRUE), heights=c(1,))
# par(mfrow=c(5, length(tasks)), mar=c(3,2,2,1), oma=c(3,3,0,0), bty='l', mgp=c(2,1,0), las=1)
#
# # Row 1: CDFs
# plot_4_tasks_cdfs(tasks)
# #mtext('A', side=2, font=2)
#
# par(mar=c(1,2,0,1))
# # Row 2: Null
# plot_4_tasks(tasks, learningModel = 'NULL', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='', margin_text='B. IID-RDM', plot_xticklabels=FALSE)
#
# # Row 3: zSM
# plot_4_tasks(tasks, learningModel = 'zSM', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='', margin_text='C. MS1-RDM', plot_xticklabels=FALSE)
#
# # Row 4: zSMuAH
# plot_4_tasks(tasks, learningModel = 'zSMuAH', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='', margin_text='D. MS2-RDM', plot_xticklabels=FALSE)
#
# # Row 5: zSMuAHbV
# plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
#              mains=c('', '', '', ''), xlab='Period', margin_text='E. MS3-RDM')
#
# # Row 6: zSMuAHbV + trends
# # plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='winning', trendPar='winning', trial_durations = trial_durations, debug=debug,
# #              mains=c('', '', '', ''), xlab='Period', full_x=TRUE, margin_text='MS3-RDM + Trends')
# dev.off()










# ACFs --------------------------------------------------------------------
get_acfs <- function(dat, lag.max=20) {
  allAcfs <- lapply(unique(dat$subjects), function(x) acf(dat[dat$subjects==x,'rt'], plot = FALSE, lag.max=lag.max)$acf)
  sharedACFs <- min(sapply(allAcfs, function(x) min(dim(x)[[1]])))
  acfsBySub <- do.call(rbind, lapply(allAcfs, function(x) x[1:sharedACFs,,]))
  #  acfsBySub <- do.call(rbind, lapply(unique(dat$subjects), function(x) acf(dat[dat$subjects==x,'rt'], plot = FALSE)$acf))
  apply(acfsBySub, 2, mean)
}

plot_4_tasks_acfs <- function(tasks, data_only=FALSE, learningModel='NULL', trendModel='NULL', trendPar='NULL', nTrendPars=3,
                              margin_text='', ...) {
  opts <- list(...)
  for(taskn in 1:length(tasks)) {
    task <- tasks[taskn]
    dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
    acfs <- get_acfs(dat)
    if(!data_only) {
      if(trendModel=='winning') {
        all_fns <- c(get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar='B', nTrendPars=3, samples_dir='./samples/')$samplers_fn,
                     get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar='v', nTrendPars=3, samples_dir='./samples/')$samplers_fn,
                     get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar='u', nTrendPars=3, samples_dir='./samples/')$samplers_fn)
        mc <- get_model_comparison(all_fns)
        winner <- row.names(mc)[which.max(mc$wBPIC)]
        winningPar <- strsplit(winner, '-')[[1]][5]
        fns <- get_fns(task, 'zSMuAHbV', trendModel='DCT', trendPar=winningPar, nTrendPars=3, samples_dir='./samples/')
      } else {
        fns <- get_fns(task, learningModel, trendModel, nTrendPars=nTrendPars, samples_dir='./samples/')
      }
      pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
      pp <- EMC2:::loadRData(pp_fn)

      acfsByPostn <- do.call(rbind, lapply(unique(pp$postn), function(x) get_acfs(pp[pp$postn==x,])))
      acfsCI <- apply(acfsByPostn, 2, quantile, c(0.025, 0.975))
    }

    plot(1:length(acfs)-1, acfs, pch=4, xlab='Lag', ylab='ACF', ylim=c(-0.05, 1)); abline(h=0)
    if(!data_only) {
      arrows(x0=2:length(acfs)-1,
             y0=acfsCI[1,2:(length(acfs))], y1=acfsCI[2,2:(length(acfs))], length=0.035, angle=90, code=3)
    }
    if(taskn==1) mtext(text=margin_text, side=2, line=3, las=0, cex=par()$cex)
    if('mains' %in% names(opts)) mtext(opts$mains[taskn], side=3, cex=par()$cex*par()$cex.main, font=2)
  }
}


pdf(file='./figures/spectra_all_acfs.pdf', width=8, height=7)

par(mfrow=c(5, length(tasks)), mar=c(3,3,2,1), oma=c(1,2,0,0), bty='l', mgp=c(2,1,0), las=1)
# Row 1: data
plot_4_tasks_acfs(tasks, data_only = TRUE, learningModel = 'NULL', trendModel='NULL', trendPar='NULL',
                  margin_text = 'Data only', mains=paste0('Dataset ', 1:4))

# Row 2: zSM
plot_4_tasks_acfs(tasks, data_only = FALSE, learningModel = 'zSM', trendModel='NULL', trendPar='NULL',
                  margin_text='MS1-RDM (SM)')

# Row 3: zSMuAH
plot_4_tasks_acfs(tasks, data_only = FALSE, learningModel = 'zSMuAH', trendModel='NULL', trendPar='NULL',
                  margin_text='MS2-RDM (SM, AH)')

# Row 4: zSMuAHbV
plot_4_tasks_acfs(tasks, data_only = FALSE, learningModel = 'zSMuAHbV', trendModel='NULL', trendPar='NULL',
                  margin_text='MS3-RDM')

# Row 5: zSMuAHbV + trends
plot_4_tasks_acfs(tasks, learningModel = 'zSMuAHbV', trendModel='winning', trendPar='winning',
                  margin_text='MS3-RDM + Trends')
dev.off()




## row 1: only null models
pdf(file='./figures/spectra_null_models_durations.pdf', width=8, height=2.5)
par(mfcol=c(1, length(tasks)), mar=c(4,3,2,1), oma=c(0,1,0,0), bty='l', mgp=c(2,1,0), las=1)
for(taskn in 1:length(tasks)) {
  task <- tasks[taskn]
  dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
  fns <- get_fns(task, 'NULL', 'NULL', samples_dir='./samples/')
  pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
  if(!file.exists(pp_fn)) {
    plot.new()
    next
  }
  pp <- EMC2:::loadRData(pp_fn)
  samplers <- EMC2:::loadRData(samplers_fn)

  main_ <- paste0('Dataset ', taskn)
  if(is.null(trial_durations[taskn])) {
    xlab_ <- 'Log frequency'
  } else {
    xlab_ <- 'Period'
  }
  ylab_ <- ''
  # if(taskn == 1) {
  #   ylim <- c(-5, -3)
  # } else if(taskn == 2) {
  #   ylim <- c()
  # }
  plotSpectrum(dat, pp=pp, pp2=NULL, xlab='', ylab=ylab_, main=main_, plot.log=TRUE, detrend=FALSE,
               trial_duration = trial_durations[taskn])
  if(taskn == 1) mtext('Log power', side=2, line=3, cex=.66,las=0)
  mtext(xlab_, side=1, line=3, cex=.66, las=0)
}
dev.off()


# IID models: CDFs with spectra ---------------------------------------------------
#pdf(file='./figures/spectra_iid_models_cdf_and_spectra.pdf', width=8, height=4)
jpeg(file=paste0('./figures/spectra_iid_models_cdf_and_spectra.jpeg'), width=8, height=4, units='in', quality=100, res=250)
#jpeg(file=)
par(mfrow=c(2, length(tasks)), mar=c(3,2,2,1), oma=c(3,3,0,0), bty='l', mgp=c(2,1,0), las=1)
par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')

# Row 1: CDFs
plot_4_tasks_cdfs(tasks)

# Row 2: spectra
plot_4_tasks(tasks, learningModel = 'NULL', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('', '', '', ''), xlab='Period (s)', margin_text='B. IID-RDM', plot_xticklabels=TRUE)

dev.off()



# IID versus winning ------------------------------------------------------
jpeg(file=paste0('./figures/spectra_iid_vs_ms3_models_spectra.jpeg'), width=8, height=6, units='in', quality=100, res=250)
par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
par(mfrow=c(3, length(tasks)), mar=c(3,2,2,1), oma=c(2,3,0,0), bty='l', mgp=c(2,1,0), las=1)
# Row 1: Spectra NULL
plot_4_tasks(tasks, learningModel = 'NULL', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4'), xlab='Period (s)', margin_text='A. IID-RDM', full_x=TRUE)

# Row 5: zSMuAHbV
plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='NULL', trendPar='NULL', trial_durations = trial_durations, debug=debug,
             mains=c('', '', '', ''), xlab='Period', full_x=TRUE, margin_text='B. MS3-RDM')

# Row 6: zSMuAHbV + trends
plot_4_tasks(tasks, learningModel = 'zSMuAHbV', trendModel='winning', trendPar='winning', trial_durations = trial_durations, debug=debug,
             mains=c('', '', '', ''), xlab='Period', full_x=TRUE, margin_text='C. MS3-RDM + Trends', plot_xticklabels=TRUE)

dev.off()





# ## Winning
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# trial_durations <- c(1.265, 2.850, 1.627, 0.837)
# pdf(file='./figures/spectra_winning_models_durations.pdf', width=8, height=2.5)
# par(mfcol=c(1, length(tasks)), mar=c(4,3,2,1), oma=c(0,1,0,0), bty='l', mgp=c(2,1,0), las=1)
# for(taskn in 1:length(tasks)) {
#   task <- tasks[taskn]
#   dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#   fns <- get_fns(task, 'zSMuAHbV', 'NULL', samples_dir='./samples/')
#   pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#   if(!file.exists(pp_fn)) {
#     plot.new()
#     next
#   }
#   pp <- EMC2:::loadRData(pp_fn)
#   samplers <- EMC2:::loadRData(samplers_fn)
#
#   main_ <- paste0('Dataset ', taskn)
#   if(is.null(trial_durations[taskn])) {
#     xlab_ <- 'Log frequency'
#   } else {
#     xlab_ <- 'Period'
#   }
#   ylab_ <- ''
#   plotSpectrum(dat, pp=pp, pp2=NULL, xlab='', ylab=ylab_, main=main_, plot.log=TRUE,
#                detrend=FALSE, trial_duration = trial_durations[taskn])
#   if(taskn == 1) mtext('Log power', side=2, line=3, cex=.66,las=0)
#   mtext(xlab_, side=1, line=3, cex=.66, las=0)
# }
# dev.off()



## only zSMuAH
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# trial_durations <- c(1.265, 2.850, 1.627, 0.837)
# pdf(file='./figures/spectra_MS2-RDM_durations-with-null.pdf', width=8, height=2.5)
# par(mfcol=c(1, length(tasks)), mar=c(4,3,2,1), oma=c(0,1,0,0), bty='l', mgp=c(2,1,0), las=1)
# for(taskn in 1:length(tasks)) {
#   task <- tasks[taskn]
#   dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#   fns <- get_fns(task, 'zSMuAH', 'NULL', samples_dir='./samples/')
#   pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#   if(!file.exists(pp_fn)) {
#     plot.new()
#     next
#   }
#   pp <- EMC2:::loadRData(pp_fn)
#   samplers <- EMC2:::loadRData(samplers_fn)
#
#   pp2 <- EMC2:::loadRData(get_fns(task, 'NULL', 'NULL', samples_dir='./samples/')$pp_fn)
#   main_ <- paste0('Dataset ', taskn)
#   if(is.null(trial_durations[taskn])) {
#     xlab_ <- 'Log frequency'
#   } else {
#     xlab_ <- 'Period'
#   }
#   ylab_ <- ''
#   plotSpectrum(dat, pp=pp, pp2=pp2, xlab='', ylab=ylab_, main=main_, plot.log=TRUE,
#                detrend=FALSE, trial_duration = trial_durations[taskn])
#   if(taskn == 1) mtext('Log power', side=2, line=3, cex=.66,las=0)
#   mtext(xlab_, side=1, line=3, cex=.66, las=0)
# }
# dev.off()

# zSM
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# trial_durations <- c(1.265, 2.850, 1.627, 0.837)
# pdf(file='./figures/spectra_SM-RDM_durations.pdf', width=8, height=2.5)
# par(mfcol=c(1, length(tasks)), mar=c(4,3,2,1), oma=c(0,1,0,0), bty='l', mgp=c(2,1,0), las=1)
# for(taskn in 1:length(tasks)) {
#   task <- tasks[taskn]
#   dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#   fns <- get_fns(task, 'zSM', 'NULL', samples_dir='./samples/')
#   pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#   if(!file.exists(pp_fn)) {
#     plot.new()
#     next
#   }
#   pp <- EMC2:::loadRData(pp_fn)
#   samplers <- EMC2:::loadRData(samplers_fn)
#
# #  pp2 <- EMC2:::loadRData(get_fns(task, 'NULL', 'NULL', samples_dir='./samples/')$pp_fn)
#   main_ <- paste0('Dataset ', taskn)
#   if(is.null(trial_durations[taskn])) {
#     xlab_ <- 'Log frequency'
#   } else {
#     xlab_ <- 'Period'
#   }
#   ylab_ <- ''
#   plotSpectrum(dat, pp=pp, pp2=NULL, xlab='', ylab=ylab_, main=main_, plot.log=TRUE,
#                detrend=FALSE, trial_duration = trial_durations[taskn])
#   if(taskn == 1) mtext('Log power', side=2, line=3, cex=.66,las=0)
#   mtext(xlab_, side=1, line=3, cex=.66, las=0)
# }
# dev.off()


#debug(plotSpectrum)
#axis(side=1, at=c())


# # Dynamic models
# trial_durations <- c(1.265, 2.850, 1.627, 0.837)
# for(ftype in c('pdf', 'png')) {
#   if(ftype == 'pdf') pdf(file='./figures/spectra_combined.pdf', width=8, height=4)
#   if(ftype == 'png') png(file='./figures/spectra_combined.png', width=8, height=4, units='in', res=100)
#   tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
#   par(mfcol=c(2, length(tasks)), mar=c(2,3,2,1), oma=c(1,1,0,0), bty='l', mgp=c(2,1,0), las=1)
#   #par(mfcol=c(2, length(tasks)), mar=c(2,3,2,1), oma=c(1,1,0,0), bty='l')
#
#   for(taskn in 1:length(tasks)) {
#     task <- tasks[taskn]
#
#     all_fns <- getFns(task)
#     mc1 <- get_model_comparison(c(all_fns[!grepl('DCT-NULL', all_fns) & grepl('-trend-DCT', all_fns) & grepl('-3', all_fns)]))
#
#
#     dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#     dctModel <- c('NULL', 'DCT')
#     learningModel <- 'zSMuAHbV' # c('zSMuAHbV', 'zSMuAHbV')
#     samples_dir = './samples/'
#     for(i in 1:2) {
#       if(i == 2) par(mar=c(3,3,1,1)) else par(mar=c(2,3,2,1))
#       if(dctModel[i] == 'DCT') {
# #        samples_dir = './samples_prior_hp1_400/'
#         if(task %in% c('wagenmakers2004_CS', 'mileticvanmaanen2019exp2block2')) dct <- 'u'
#         if(task == 'forstmann2008') dct <- 'v'
#         if(task== 'wagenmakers2008exp2') dct <- 'B'
#       } else {
#         dct <- 'NULL'
#       }
# #      lm <- learningModel[1] # ifelse(task=='mileticvanmaanen2019exp2block2' & learningModel[i] != 'NULL', 'zSMuAHbV', learningModel[i])
#       if(i==2) {
#         samplers_fn <- paste0('./samples/', task, '_model-', row.names(mc1[which.min(mc1$Dmean),]), '.RData')
#         pp_fn <- gsub('.RData', '_pp.RData', samplers_fn)
#       } else {
#         fns <- get_fns(task, learningModel, dct, samples_dir=samples_dir)
#         pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#       }
#       if(!file.exists(pp_fn)) {
#         plot.new()
#         next
#       }
#       pp <- EMC2:::loadRData(pp_fn)
#       samplers <- EMC2:::loadRData(samplers_fn)
#
#       if(i == 1) main_ <- paste0('Dataset ', taskn) else main_ <- ''
#       if(i == 2) {
#         if(is.null(trial_durations[taskn])) {
#           xlab_ <- 'Log frequency'
#         } else {
#           xlab_ <- 'Period'
#         }
#       } else {
#         xlab_ <- ''
#       }
#       plotSpectrum(dat, pp, pp2=NULL, xlab='', ylab='', main=main_, plot.log=TRUE, detrend=FALSE, demean=TRUE,
#                    trial_duration = trial_durations[taskn], full_x = TRUE)
#       if(taskn == 1) mtext('Log power', side=2, line=3, cex=par()$cex, las=0)#  ylab_ <- 'Log power' else ylab_ <- ''
#       mtext(xlab_, side=1, line=3, cex=par()$cex, las=1)#  ylab_ <- 'Log power' else ylab_ <- ''
#     }
#   }
#   dev.off()
# }
#
#
#
# # Combine CDFs with spectra -----------------------------------------------
# # Slightly adjusted plot_fit function
# plot_fit <- function(data,pp,subject=NULL,factors=NULL,functions=NULL,
#                      stat=NULL,stat_name="",adjust=1,
#                      quants=c(.025,.5,.975),do_plot=TRUE,
#                      xlim=NULL,ylim=NULL,
#                      layout=NULL,mfcol=FALSE,
#                      probs=c(1:99)/100,
#                      data_lwd=2,fit_lwd=1,
#                      q_points=c(.1,.3,.5,.7,.9),
#                      qp_cex=1,pqp_cex=.5,lpos="topleft", main = "")
# {
#   if (!is.null(stat) & is.null(factors)) factors <- NA
#   if (!is.null(subject)) {
#     snams <- levels(data$subjects)
#     if (is.numeric(subject)) subject <- snams[subject]
#     if (!all(subject %in% snams)) stop("Subject(s) not present\n")
#     dat <- droplevels(data[data$subjects %in% subject,])
#     pp <- droplevels(pp[pp$subjects %in% subject,])
#     if (length(subject>1))
#       fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))] else
#         fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt"))]
#   } else {
#     dat <- data
#     fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))]
#   }
#
#   okd <- !is.na(dat$R) & is.finite(dat$rt)
#   okpp <- !is.na(pp$R) & is.finite(pp$rt)
#
#   if (!is.null(functions)) for (i in 1:length(functions)) {
#     dat <- cbind.data.frame(functions[[i]](dat),dat)
#     names(dat)[1] <- names(functions)[i]
#     pp <- cbind.data.frame(functions[[i]](pp),pp)
#     names(pp)[1] <- names(functions)[i]
#     fnams <- c(names(functions)[i],fnams)
#   }
#
#   if (!is.null(factors)) {
#     if (any(is.na(factors))) fnams <- NA else {
#       if (!all(factors %in% fnams))
#         stop("factors must name factors in data")
#       fnams <- factors
#     }
#   }
#   if (!any(is.na(layout))) if (!is.null(layout))
#     if (mfcol) par(mfcol=layout) else par(mfrow=layout)
#
#   if (all(is.na(data$rt))) stop("Use plot_fit_choice if no rt data")
#
#   if (!any(is.na(fnams))) {
#     cells <- dat[,fnams,drop=FALSE]
#     for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
#     cells <- apply(cells,1,paste,collapse=" ")
#     pp_cells <- pp[,fnams,drop=FALSE]
#     for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
#     pp_cells <- apply(pp_cells,1,paste,collapse=" ")
#   }
#   if (!is.null(stat)) { # statistic
#     postn <- unique(pp$postn)
#     if (any(is.na(fnams))) ucells <- "" else ucells <- sort(unique(cells))
#     tab <- matrix(nrow=length(ucells),ncol=4,
#                   dimnames=list(ucells,c("Observed",names(quantile(1:5,quants)))))
#     for (i in ucells) {
#       if (i=="") {
#         dati <- dat
#         ppi <- pp
#         obs <- stat(dati)
#         pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
#         tab[1,] <- c(obs,quantile(pred,quants))
#       } else {
#         dati <- dat[cells==i,]
#         ppi <- pp[pp_cells==i,]
#         obs <- stat(dati)
#         pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
#         tab[i,] <- c(obs,quantile(pred,quants))
#       }
#       if (do_plot) {
#         dens <- density(pred,adjust=adjust)
#         if (!is.null(xlim)) xlimi <- xlim else
#           xlimi <- c(pmin(obs,min(dens$x)),pmax(obs,max(dens$x)))
#         plot(dens,main=paste(main,i),xlab=stat_name,xlim=xlimi)
#         abline(v=obs)
#       }
#     }
#     invisible(tab)
#   } else { # cdf
#     if (any(is.na(fnams))) cells <- ""
#     pok <- probs %in% q_points
#     R <- levels(dat$R)
#     if (is.null(ylim)) ylim <- c(0,1)
#     # set common xlim
#     if (is.null(xlim)) {
#       xlim <- c(Inf,-Inf)
#       for (i in sort(unique(cells))) {
#         if (i=="") {
#           dati <- dat
#           ppi <- pp
#         } else {
#           dati <- dat[cells==i & okd,]
#           ppi <- pp[pp_cells==i & okpp,]
#         }
#         pqs <- pq <- qs <- setNames(vector(mode="list",length=length(R)),R)
#         for (j in R) if (length(dati$rt[dati$R==j])>=length(q_points)) {
#           qs[[j]] <- quantile(dati$rt[dati$R==j],probs=probs)
#           pq[[j]] <- quantile(ppi$rt[ppi$R==j],probs=probs)
#           pqs[[j]] <- tapply(ppi$rt[ppi$R==j],ppi$postn[ppi$R==j],
#                              quantile,probs=probs[pok])
#         } else qs[[j]] <- pq[[j]] <- pqs[[j]] <- NA
#         rx <- cbind(do.call(rbind,lapply(qs,function(x){x[c(1,length(probs))]})),
#                     do.call(rbind,lapply(pq,function(x){x[c(1,length(probs))]})))
#         xlimi <- c(min(rx,na.rm=TRUE),max(rx,na.rm=TRUE))
#         if (!any(is.na(xlimi))) {
#           xlim[1] <- pmin(xlim[1],xlimi[1])
#           xlim[2] <- pmax(xlim[2],xlimi[2])
#         }
#       }
#     }
#     for (i in sort(unique(cells))) {
#       if (i=="") {
#         dati <- dat
#         ppi <- pp
#         okdi <- okd
#         okppi <- okpp
#       } else {
#         dati <- dat[cells==i,]
#         ppi <- pp[pp_cells==i,]
#         okdi <- okd[cells==i]
#         okppi <- okpp[pp_cells==i]
#       }
#       pR <- tapply(okdi,dati$R,sum)/dim(dati)[1]
#       ppR <- tapply(okppi,ppi$R,sum)/dim(ppi)[1]
#       dati <- dati[okdi,]
#       ppi <- ppi[okppi,]
#       pqs <- pq <- qs <- setNames(vector(mode="list",length=length(R)),R)
#       for (j in R) if (length(dati$rt[dati$R==j])>=length(q_points)) {
#         isj <- ppi$R==j
#         qs[[j]] <- quantile(dati$rt[dati$R==j],probs=probs)
#         pq[[j]] <- quantile(ppi$rt[isj],probs=probs)
#         pqs[[j]] <- tapply(ppi$rt[isj],ppi$postn[isj],quantile,probs=probs[pok])
#       } else qs[[j]] <- pq[[j]] <- pqs[[j]] <- NA
#       if ( !any(is.na(pq[[1]])) ) {
#         plot(pq[[1]],probs*ppR[1],xlim=xlim,ylim=ylim,main=main,#paste(main,i),
#              xlab="RT",type="l",
#              lwd=fit_lwd,ylab="",lty=1, col='dark green') #"p(R)",lty=1)
#         tmp=lapply(pqs[[1]],function(x){
#           points(x,probs[pok]*ppR[1],col=adjustcolor('darkgreen', alpha.f=.8),pch=16,cex=pqp_cex)})
#         points(pq[[1]][pok],probs[pok]*ppR[1],cex=pqp_cex*3,pch=16,col="darkgreen")
#         lines(qs[[1]],probs*pR[1],lwd=data_lwd,lty=1)
#         points(qs[[1]][pok],probs[pok]*pR[1],cex=qp_cex,pch=16)
#         do_plot=FALSE
#       } else do_plot=TRUE
#       if (length(qs)>1) {
#         for (j in 2:length(qs)) if (!any(is.na(pq[[j]]))) {
#           if (do_plot) {
#             plot(pq[[j]],probs*ppR[j],xlim=xlim,ylim=ylim,main=main, #paste(main,i),
#                  xlab="RT",type="l",
#                  lwd=fit_lwd,ylab="", lty=j, col='darkgreen') #p(R)",lty=j)
#             do_plot <- FALSE
#           } else lines(pq[[j]],probs*ppR[j],lwd=fit_lwd,lty=j)
#           tmp=lapply(pqs[[j]],function(x){
#             points(x,probs[pok]*ppR[j],col=adjustcolor('darkgreen', alpha.f=.8),pch=16,cex=pqp_cex)})
#           points(pq[[j]][pok],probs[pok]*ppR[j],cex=pqp_cex*3,pch=16,col='darkgreen')
#           lines(qs[[j]],probs*pR[j],lwd=data_lwd,lty=j)
#           points(qs[[j]][pok],probs[pok]*pR[j],cex=qp_cex,pch=16)
#         }
#       }
#       if(!is.null(lpos)) legend(lpos,as.character(R),lty=1:length(R),bty="n")#,title="Response")
#     }
#   }
# }
#
#
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# for(colors_ in c('white')) {#}, 'white')) {
#   for(ftype in c('pdf')) {#}, 'pdf')) {
#     if(colors_ == 'black') {
#       if(ftype == 'pdf') {
#         pdf(file=paste0('./figures/spectra_null_models_with_cdf_black.pdf'), width=8, height=4)
#       } else {
#         jpeg(file=paste0('./figures/spectra_null_models_with_cdf_black.jpeg'), width=8, height=4, units='in', quality=100, res=250)
#       }
#       #pdf(file=file.path(figures_dir, 'stimulus_memory_effects_group_stimuluscoded_black.pdf'), width=8, height=5)
#       if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
#     } else {
#       if(ftype == 'pdf') {
#         pdf(file=paste0('./figures/spectra_null_models_with_cdf.pdf'), width=8, height=4)
#       } else {
#         jpeg(file=paste0('./figures/spectra_null_models_with_cdf.jpeg'), width=8, height=4, units='in', quality=100, res=250)
#       }
# #      pdf(file=paste0('./figures/spectra_null_models_with_cdf.pdf'), width=8, height=4)
#     }
#
#   par(mfcol=c(2, length(tasks)), mar=c(3,3,2,1), oma=c(1,1,0,0), bty='l', mgp=c(1.5,0.6,0), las=1)
#   for(taskn in 1:4) {
#     task <- tasks[taskn]
#     dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#     fns <- get_fns(task, 'NULL', 'NULL', samples_dir='./samples/')
#     pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#     pp <- EMC2:::loadRData(pp_fn)
#     samplers <- EMC2:::loadRData(samplers_fn)
#
#     dat$R <- factor(as.logical(dat$accuracy), levels=c(TRUE, FALSE), labels=c('Correct', 'Error'))
#     pp$R <- factor(as.logical(dat$accuracy), levels=c(TRUE, FALSE), labels=c('Correct', 'Error'))
#     dat$S <- factor('', levels=c(''))
#     pp$S <- factor('', levels=c(''))
#
#     #par(mfrow=c(2,1),bty='l')
#     main_ <- paste0('Dataset ', taskn)
#     if(taskn==1) lpos = 'topleft' else lpos=NULL
#     par(mar=c(2.5,3,2,1))
#     plot_fit(data=dat, pp=pp, factors='S', layout=NULL, main = main_, lpos=lpos)
#     if(taskn == 1) mtext('p(R)', side=2, line=2.5, cex=.66,las=0)
#     if(is.null(trial_durations[taskn])) {
#       xlab_ <- 'Log frequency'
#     } else {
#       xlab_ <- 'Period'
#     }
#     ylab_ <- ''
#     par(mar=c(3,3,0.5,1))
#     plotSpectrum(dat, pp=pp, #NULL, #pp,
#                  pp2=NULL, xlab='', ylab=ylab_, main='', plot.log=TRUE, detrend=FALSE,trial_duration = trial_durations[taskn])
#     if(taskn == 1) mtext('Log power', side=2, line=2.5, cex=.66,las=0)
#     mtext(xlab_, side=1, line=3, cex=par()$cex, las=1)#  ylab_ <- 'Log power' else ylab_ <- ''
# #    break
#   }
#   dev.off()
#   }
# }
# out <- plot_fit(data=dat, pp=pp, factors='S', layout=NULL,do_plot = FALSE, main = 'Dataset 1')
#
#
#
# # Null vs winning models --------------------------------------------------
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# for(colors_ in c('white')) {#} 'white')) {
#   for(ftype in c('pdf')) {
#     if(colors_ == 'black') {
#       if(ftype == 'pdf') {
#         pdf(file=paste0('./figures/spectra_null_models_with_winning_black.pdf'), width=8, height=4)
#       } else {
#         jpeg(file=paste0('./figures/spectra_null_models_with_winning_black.jpeg'), width=8, height=4, units='in', quality=100, res=250)
#       }
#       #pdf(file=file.path(figures_dir, 'stimulus_memory_effects_group_stimuluscoded_black.pdf'), width=8, height=5)
#       if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
#     } else {
#       if(ftype == 'pdf') {
#         pdf(file=paste0('./figures/spectra_null_models_with_winning.pdf'), width=8, height=4)
#       } else {
#         jpeg(file=paste0('./figures/spectra_null_models_with_winning.jpeg'), width=8, height=4, units='in', quality=100, res=250)
#       }
#       #      pdf(file=paste0('./figures/spectra_null_models_with_cdf.pdf'), width=8, height=4)
#     }
#
# #    par(mfcol=c(2, length(tasks)), mar=c(3,3,2,1), oma=c(0,1,0,0), bty='l', mgp=c(1.5,0.6,0), las=1)
#
#     par(mfcol=c(2, length(tasks)), mar=c(2,3,2,1), oma=c(1,1,0,0), bty='l', mgp=c(2,1,0), las=1)
#     #par(mfcol=c(2, length(tasks)), mar=c(2,3,2,1), oma=c(1,1,0,0), bty='l')
#     for(taskn in 1:length(tasks)) {
#       task <- tasks[taskn]
#       dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#       dctModel <- c('NULL', 'DCT')
#       learningModel <- c('NULL', 'zSMuAHbV')
#       samples_dir='./samples/'
#       dct <- 'NULL'
#
#       #   fns <- get_fns(task, 'NULL', 'NULL', samples_dir='./samples/')
#       #   pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#       #   pp <- EMC2:::loadRData(pp_fn)
#       #   samplers <- EMC2:::loadRData(samplers_fn)
#
#       for(i in 1:2) {
#         if(i == 2) par(mar=c(3,3,1,1)) else par(mar=c(2,3,2,1))
#         lm <- learningModel[i] # ifelse(task=='mileticvanmaanen2019exp2block2' & learningModel[i] != 'NULL', 'zSMuAHbV', learningModel[i])
#         fns <- get_fns(task, lm, dct, samples_dir=samples_dir)
#         pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#         if(!file.exists(pp_fn)) {
#           plot.new()
#           next
#         }
#         pp <- EMC2:::loadRData(pp_fn)
#         samplers <- EMC2:::loadRData(samplers_fn)
#
#         if(i == 1) main_ <- paste0('Dataset ', taskn) else main_ <- ''
#         if(i == 2) {
#           if(is.null(trial_durations[taskn])) {
#             xlab_ <- 'Log frequency'
#           } else {
#             xlab_ <- 'Period'
#           }
#         } else {
#           xlab_ <- ''
#         }
#  #       plotSpectrum(dat, pp, pp2=NULL, xlab='', ylab='', main=main_, plot.log=TRUE, detrend=FALSE,trial_duration = trial_durations[taskn])
# #        if(i == 2) xlab_ <- 'Log frequency' else xlab_ <- ''
#         plotSpectrum(dat, pp, pp2=NULL, xlab='', ylab='', main=main_, plot.log=TRUE, detrend=FALSE,trial_duration = trial_durations[taskn])
#         if(taskn == 1) mtext('Log power', side=2, line=3, cex=par()$cex, las=0)#  ylab_ <- 'Log power' else ylab_ <- ''
#         mtext(xlab_, side=1, line=3, cex=.66, las=0)
#       }
#     }
#     dev.off()
#   }
# }
#
#
#
#
#
# # Compare trend models ----------------------------------------------------
# ## Plot trends themselves
# getWaves <- function(samplers, subject, filter='sample', pname='Bcos', nwaves=100, return_long=TRUE) {
#   idxStage<- samplers$samples$stage==filter
#   pnames <- samplers$par_names
#   idxPars <- grepl(pname, pnames)
#
#   pars <- samplers$samples$alpha[idxPars,subject,idxStage]
#   if(sum(idxPars) == 1) pars <- t(pars)
#   # pars should be [nPars, nSamples]
#
#   if(grepl('cos', pname)) {
#     X <- attr(samplers$data[[subject]], 'adapt')[[1]]$dct
#   } else if(grepl('poly', pname)) {
#     X <- attr(samplers$data[[subject]], 'adapt')[[1]]$poly
#   }
# #  dcts <- attr(samplers[[1]]$data[[subject]], 'adapt')[[1]]$dct
#   #if(nrow(pars) == 0) return(matrix(0, nrow=nrow(X), ncol=2))
#   X <- X[,1:nrow(pars),drop = F]
#
#   ## get 100 random waves
#   if(ncol(pars) > nwaves) {
#     pars <- pars[,sample(1:ncol(pars), size=nwaves),drop = F]
#   }
#
#   waves <- NULL
#   for(i in 1:ncol(pars)) waves <- cbind(waves, X%*%pars[,i])
#   if(return_long) return(reshape2:::melt(waves, varnames=c('trial', 'postn')))
#   return(waves)
# }
#
# #samplers <- EMC2:::loadRData(gsub('_pp', '', fn))
# #samples <- EMC2::merge_samples(samplers)
# #undebug(getWaves)
#
# plotWaves <- function(samplers, pname, ylab='Threshold', ...) {
#   samples <- merge_samples(samplers)
#   allWaves <- do.call(rbind, lapply(names(samplers[[1]]$data), function(x) {
#     y <- getWaves(samples, x, pname = pname); y$subjects<-x;y}))
#   meanWaveB <- aggregate(value~trial,aggregate(value~trial*postn,allWaves,mean), quantile, c(.025, .5, .975))
#   plot(meanWaveB[,1], meanWaveB[,2][,2], type='l', ylim=range(meanWaveB[2]), xlab='Trials', ylab=ylab, ...)
#   abline(h=0, col='lightgray', lty=2)
#   polygon(c(meanWaveB[,1], rev(meanWaveB[,1])),
#           c(meanWaveB[,2][,1], rev(meanWaveB[,2][,3])), col=adjustcolor(2, alpha.f=0.2), border=FALSE)
# }
#
#
#
#
#
#
# #task <- 'wagenmakers2004_CS'
# pdf(file='./figures/all_trend_shapes-57.pdf', width=12, height=8)
# # taskn <- 1
# for(taskn in 1:4) {
#   task <- tasks[taskn]
#   all_fns <- getFns(task)
#   dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#   mc1 <- get_model_comparison(c(all_fns[!grepl('DCT-NULL', all_fns) & grepl('-trend-', all_fns)], paste0('./samples/', task, '_model-RDM-zSMuAHbV-DCT-NULL.RData')))
#   mc2 <- get_model_comparison(c(all_fns[!grepl('DCT-NULL', all_fns) & grepl('-trend-', all_fns)]))
#
#   # titles:
#   layout(rbind(c(1,1,2,2,3,3),matrix(data=1:18+3, nrow=3),c(1,1,2,2,3,3)), heights=c(0.025, .33,.33,.33,0.01))
#   par(mar=c(2,3,2,1), oma=c(1,2,1,0), bty='l', mgp=c(2,1,0))
#   plot.new(); mtext('Threshold', font=2, cex=par()$cex*par()$cex.main);
#   plot.new(); mtext('Urgency', font=2, cex=par()$cex*par()$cex.main); mtext(paste0("Dataset ", taskn), side=3,font=2, cex=par()$cex*par()$cex.main,line=1.5);
#   plot.new(); mtext('Evidence quality', font=2, cex=par()$cex*par()$cex.main)
#   for(trendPar in c('B', 'u', 'v')) {
#     for(trendModel in c('DCT', 'poly')) {
#       for(nPars in 5:7) { #1:3) {
#         fn <- paste0('./samples/', task, '_model-RDM-zSMuAHbV-trend-', trendModel, '-', trendPar, '-', nPars, '.RData')
#         if(!file.exists(fn)) {plot.new(); next}
#         samplers <- EMC2:::loadRData(fn)
#         if(chain_n(samplers)[1,4] == 0) {plot.new(); next}
#         wbpic = mc1[paste0('RDM-zSMuAHbV-trend-', trendModel, '-',trendPar, '-', nPars),'wBPIC']
#         maxgd = mc1[paste0('RDM-zSMuAHbV-trend-', trendModel, '-',trendPar, '-', nPars),'max.gd']
#         wbpic2 = mc2[paste0('RDM-zSMuAHbV-trend-', trendModel, '-',trendPar, '-', nPars),'wBPIC']
#
#         trendModelName <- ifelse(trendModel=='DCT','cos','poly')
#         pname <- paste0(trendPar, trendModelName)
#         main_ <- paste0(trendPar,'~',trendModelName)
#
#         plotWaves(samplers, ylab='', pname=pname, main=main_)
#         legend('topleft', c(paste0('wBPIC=', round(wbpic,3)),
#                             paste0('wBPIC2=', round(wbpic2,3)),
#                             '=(excl null model)'),
#                lty=c(NA,NA,NA), lwd=c(NA,NA,NA),
#                col=c(1,1,1),
#                text.col=c(ifelse(max(mc1['wBPIC'])==wbpic, 4,1),
#                           ifelse(max(mc2['wBPIC'])==wbpic2, 4,1),
#                           1),
#                bty='n', x.intersp = -2)
#         legend('bottomright', c(paste0('max gelmans=', round(maxgd,2))), lty=NA, lwd=NA, col=1, bty='n', x.intersp = 22,
#                text.col=ifelse(maxgd<1.1, 1, 2))
#         if(trendPar=='B'&trendModel=='DCT') mtext(text=paste0(nPars, ' parameter(s)'), side=2, cex=par()$cex*par()$cex.axis, las=0, line=2.5, font=2)
#       }
#     }
#   }
# }
# dev.off()
#
#
# ## Same thing but now plot spectra
# # DCT x Polynomials
# # B, v, u
# # 2-5 parameters
# ## -> 6 columns, 4 rows
#
# trial_durations <- c(1.265, 2.850, 1.627, 0.837)
# # #tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# #
# # par(mfcol=c(3, 6), mar=c(2,3,2,1), oma=c(1,1,0,0), bty='l', mgp=c(2,1,0), las=1)
# # taskn <- 2
# # task <- tasks[taskn]
# # dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
# # for(trendModel in c('DCT', 'poly')) {
# #   for(trendPar in c('B', 'u', 'v')) {
# #     for(nPars in 1:3) {
# #       #if(trendModel == 'DCT') {
# #       #  fn <- paste0('./samples/', task, '_model-RDM-zSMuAHbV-', trendModel, '-', nPars, trendPar, '_pp.RData')
# #       #} else {
# #       fn <- paste0('./samples/', task, '_model-RDM-zSMuAHbV-trend-', trendModel, '-', trendPar, '-', nPars, '_pp.RData')
# #       #}
# #       if(!file.exists(fn)) {plot.new(); next}
# #       pp <- EMC2:::loadRData(fn)
# #       xlab_ <- ifelse(nPars==5, 'Period', '')
# #       main_ <- paste0(trendModel, '-', nPars, trendPar)
# #       plotSpectrum(dat, pp, pp2=NULL, xlab=xlab_, ylab='', main=main_, plot.log=TRUE, detrend=FALSE, trial_duration = trial_durations[taskn])
# #     }
# #   }
# # }
#
#
# pdf(file='./figures/all_trend_spectra-57.pdf', width=12, height=8)
# # taskn <- 1
# for(taskn in 1:4) {
#   task <- tasks[taskn]
#   all_fns <- getFns(task)
#   dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#   mc1 <- get_model_comparison(c(all_fns[!grepl('DCT-NULL', all_fns) & grepl('-trend-', all_fns)], paste0('./samples/', task, '_model-RDM-zSMuAHbV-DCT-NULL.RData')))
#   mc2 <- get_model_comparison(c(all_fns[!grepl('DCT-NULL', all_fns) & grepl('-trend-', all_fns)]))
#
#   # titles:
#   layout(rbind(c(1,1,2,2,3,3),matrix(data=1:18+3, nrow=3),c(1,1,2,2,3,3)), heights=c(0.025, .33,.33,.33,0.01))
#   par(mar=c(2,3,2,1), oma=c(1,2,1,0), bty='l', mgp=c(2,1,0))
#   plot.new(); mtext('Threshold', font=2, cex=par()$cex*par()$cex.main);
#   plot.new(); mtext('Urgency', font=2, cex=par()$cex*par()$cex.main); mtext(paste0("Dataset ", taskn), side=3,font=2, cex=par()$cex*par()$cex.main,line=1.5);
#   plot.new(); mtext('Evidence quality', font=2, cex=par()$cex*par()$cex.main)
#   for(trendPar in c('B', 'u', 'v')) {
#     for(trendModel in c('DCT', 'poly')) {
#       for(nPars in 5:7) { #1:3) {
# #        fn <- paste0('./samples/', task, '_model-RDM-zSMuAHbV-trend-', trendModel, '-', trendPar, '-', nPars, '.RData')
# #        if(!file.exists(fn)) {plot.new(); next}
#         fn <- paste0('./samples/', task, '_model-RDM-zSMuAHbV-trend-', trendModel, '-', trendPar, '-', nPars, '_pp.RData')
#         if(!file.exists(fn)) {plot.new(); next}
#         pp <- EMC2:::loadRData(fn)
#
#         wbpic = mc1[paste0('RDM-zSMuAHbV-trend-', trendModel, '-',trendPar, '-', nPars),'wBPIC']
#         maxgd = mc1[paste0('RDM-zSMuAHbV-trend-', trendModel, '-',trendPar, '-', nPars),'max.gd']
#         wbpic2 = mc2[paste0('RDM-zSMuAHbV-trend-', trendModel, '-',trendPar, '-', nPars),'wBPIC']
#
#         trendModelName <- ifelse(trendModel=='DCT','cos','poly')
#         pname <- paste0(trendPar, trendModelName)
#         main_ <- paste0(trendPar,'~',trendModelName)
#
#         xlab_ <- ifelse(nPars==3, 'Period', '')
#         #main_ <- paste0(trendModel, '-', nPars, trendPar)
#         plotSpectrum(dat, pp, pp2=NULL, xlab=xlab_, ylab='', main=main_, plot.log=TRUE, detrend=FALSE, trial_duration = trial_durations[taskn], full_x = TRUE)
#
#         # plotWaves(samplers, ylab='', pname=pname, main=main_)
#         legend('topright', c(paste0('wBPIC=', round(wbpic,3)),
#                             paste0('wBPIC2=', round(wbpic2,3)),
#                             '=(excl null model)'),
#                lty=c(NA,NA,NA), lwd=c(NA,NA,NA),
#                col=c(1,1,1),
#                text.col=c(ifelse(max(mc1['wBPIC'])==wbpic, 4,1),
#                           ifelse(max(mc2['wBPIC'])==wbpic2, 4,1),
#                           1),
#                bty='n')
#         legend('bottomleft', c(paste0('max gelmans=', round(maxgd,2))), lty=NA, lwd=NA, col=1, bty='n', x.intersp = 22,
#                text.col=ifelse(maxgd<1.1, 1, 2))
#         if(trendPar=='B'&trendModel=='DCT') mtext(text=paste0(nPars, ' parameter(s)'), side=2, cex=par()$cex*par()$cex.axis, las=0, line=2.5, font=2)
#       }
#     }
#   }
# }
# dev.off()
#
#
# # ## Wagenmakers 2008 exp 2: Baseline difference, zSM, baseline+zSM -------
# tasks <- c('wagenmakers2008exp2', 'wagenmakers2008exp2', 'wagenmakers2008exp2')
# trial_durations <- c(0.837, 0.837, 0.837)
# learningModels <- c('NULL', 'zSMuAHbV', 'zSMuAHbV_with_baseline')
# mains <- c('Bias only', 'MS3-RDM', 'MS3-RDM + bias')
#
# pdf(file='./figures/spectra_wagenmakers2008exp2_comparison.pdf', width=6, height=2.5)
# par(mfcol=c(1, length(tasks)), mar=c(4,3,2,1), oma=c(0,1,0,0), bty='l', mgp=c(2,1,0), las=1)
# for(taskn in 1:length(tasks)) {
#   task <- tasks[taskn]
#   dat <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
#   fns <- get_fns(task, learningModels[taskn], 'NULL', samples_dir='./samples/')
#   pp_fn <- fns[[1]]; samplers_fn <- fns[[2]]
#   if(!file.exists(pp_fn)) { plot.new(); next }
#   pp <- EMC2:::loadRData(pp_fn)
#   samplers <- EMC2:::loadRData(samplers_fn)
#
#   if(is.null(trial_durations[taskn])) {
#     xlab_ <- 'Log frequency'
#   } else {
#     xlab_ <- 'Period'
#   }
#   ylab_ <- ''
#   # if(taskn == 1) {
#   #   ylim <- c(-5, -3)
#   # } else if(taskn == 2) {
#   #   ylim <- c()
#   # }
#   plotSpectrum(dat, pp=pp, pp2=NULL, xlab='', ylab=ylab_, main=mains[taskn], plot.log=TRUE, detrend=FALSE,
#                trial_duration = trial_durations[taskn])
#   if(taskn == 1) mtext('Log power', side=2, line=3, cex=.66,las=0)
#   mtext(xlab_, side=1, line=3, cex=.66, las=0)
# }
# dev.off()
#
#
