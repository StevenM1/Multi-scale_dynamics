## Figures for paper
rm(list=ls())
library(EMC2)
library(emcAdapt)
source('./extra_EMC2_functions/adaptive.R')
source('./extra_EMC2_functions/make_data.R')
source('./extra_EMC2_functions/model_RDMdynamic.R')
source('./extra_EMC2_functions/utils.R')
figures_dir = './figures'

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
  par(mar=c(1,2,2,1))
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

  par(mar=c(4,2,2,1))
  shift = 0.1
  ylim = range(c(dataMeans$choice-dataSEs$choice, dataMeans$choice+dataSEs$choice, dataMeans$stim-dataSEs$stim, dataMeans$stim+dataSEs$stim,
                 ppQuantiles$choice))
  if(plotStimulusPattern) {
    plot(x=as.numeric(dataMeans$pattern), y=dataMeans$stim, lwd=2, ylim=ylim, xaxt='n', ylab=ylab2, xlab='Pattern', pch=4, main=main2)
    warning('Plotting *stim* as a function of previous stimuli')
  } else {
    plot(x=as.numeric(dataMeans$pattern), y=dataMeans$choice, lwd=2, ylim=ylim, xaxt='n', ylab=ylab2, xlab='Pattern', pch=4, main=main2)
  }
  mtext(ylab2, side=2, line=2.5, cex=.66)
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



# Sequential effects MS3-RDM ------------------------------------------------------
wagenmakers2004_CS <- loadDataSamplersPP('wagenmakers2004_CS', 'zSMuAHbV') #, 'NULL')
forstmann2008 <- loadDataSamplersPP('forstmann2008', 'zSMuAHbV')
mileticvanmaanen2019exp2block2 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', 'zSMuAHbV')
wagenmakers2008exp2 <- loadDataSamplersPP('wagenmakers2008exp2', 'zSMuAHbV')

# By choice
seqEffectsSummary1 <- getSequentialEffectsBySub(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, repeatAlternate = FALSE)
seqEffectsSummary1RA <- getSequentialEffectsBySub(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, repeatAlternate = TRUE)
seqEffectsSummary2 <- getSequentialEffectsBySub(dat=forstmann2008$dat, pp=forstmann2008$pp, repeatAlternate = FALSE)
seqEffectsSummary2RA <- getSequentialEffectsBySub(dat=forstmann2008$dat, pp=forstmann2008$pp, repeatAlternate = TRUE)
seqEffectsSummary3 <- getSequentialEffectsBySub(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, repeatAlternate = FALSE)
seqEffectsSummary3RA <- getSequentialEffectsBySub(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, repeatAlternate = TRUE)
seqEffectsSummary4 <- getSequentialEffectsBySub(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, repeatAlternate = FALSE)
seqEffectsSummary4RA <- getSequentialEffectsBySub(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, repeatAlternate = TRUE)

# Plots: Group, by stimulus
for(colors_ in c('black', 'white')) {
  if(colors_ == 'black') {
    pdf(file=file.path(figures_dir, 'SH_MS3RDM_group_stimuluscoded_black.pdf'), width=8, height=5)
    if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
  } else {
    pdf(file=file.path(figures_dir, 'SH_MS3RDM_group_stimuluscoded.pdf'), width=8, height=5)
  }
  par(mfcol=c(2,4), bty='l', oma=c(0,2,0,0))
  plotSequentialEffectsGroup(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, seqEffectsSummary=seqEffectsSummary1,
                             nHistory = 3, main1='Dataset 1', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=forstmann2008$dat, pp=forstmann2008$pp, seqEffectsSummary=seqEffectsSummary2,
                             nHistory = 3, main1='Dataset 2', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, seqEffectsSummary=seqEffectsSummary3,
                             nHistory = 3, main1='Dataset 3', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, seqEffectsSummary=seqEffectsSummary4,
                             nHistory = 3, main1='Dataset 4', full.legend=FALSE)
  dev.off()
}

# Group, by RA
pdf(file=file.path(figures_dir, 'SH_MS3RDM_group_RAcoded.pdf'), width=8, height=5)
par(mfcol=c(2,4), bty='l')
plotSequentialEffectsGroup(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, seqEffectsSummary=seqEffectsSummary1RA,
                           nHistory = 3, main1='Dataset 1', ylab1='RT (s)', ylab2='Choice proportion', full.legend=TRUE, repeatAlternate = TRUE)
plotSequentialEffectsGroup(dat=forstmann2008$dat, pp=forstmann2008$pp, seqEffectsSummary=seqEffectsSummary2RA,
                           nHistory = 3, main1='Dataset 2', full.legend=FALSE, repeatAlternate = TRUE)
plotSequentialEffectsGroup(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, seqEffectsSummary=seqEffectsSummary3RA,
                           nHistory = 3, main1='Dataset 3', full.legend=FALSE, repeatAlternate = TRUE)
plotSequentialEffectsGroup(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, seqEffectsSummary=seqEffectsSummary4RA,
                           nHistory = 3, main1='Dataset 4', full.legend=FALSE, repeatAlternate = TRUE)
dev.off()

# Individual
for(colors_ in c('black', 'white')) {
  if(colors_ == 'black') {
    pdf(file=file.path(figures_dir, 'SH_MS3RDM_individual_RAcoded_ordered_black.pdf'), width=8, height=5)
    if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
  } else {
    pdf(file=file.path(figures_dir, 'SH_MS3RDM_individual_RAcoded_ordered.pdf'), width=8, height=5)
  }
  layout(matrix(1:8, nrow=2), widths = c(.15,.25,.5,.2))
  #par(mar=c(2,3,2,1), bty='l', mgp=c(2,1,0), oma=c(0,0,0,0))
  par(mar=c(2,3,2,.5), mgp=c(2,.75,0), bty='l', oma=c(0,0,0,0))
  plotSequentialEffectsIndividual(seqEffectsSummary1RA, main='Dataset 1', xlab2='', ylab2='P(R|RRR)', ylab1='RT(A|RRR) - RT(R|RRR)', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
  plotSequentialEffectsIndividual(seqEffectsSummary2RA, main='Dataset 2', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
  plotSequentialEffectsIndividual(seqEffectsSummary3RA, main='Dataset 3', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
  plotSequentialEffectsIndividual(seqEffectsSummary4RA, main='Dataset 4', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  dev.off()
}


# Group level of only SM and SM+AH -----------
wagenmakers2004_CS <- loadDataSamplersPP('wagenmakers2004_CS', 'zSM')
forstmann2008 <- loadDataSamplersPP('forstmann2008', 'zSM')
mileticvanmaanen2019exp2block2 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', 'zSM')
wagenmakers2008exp2 <- loadDataSamplersPP('wagenmakers2008exp2', 'zSM')

# By choice
seqEffectsSummary1 <- getSequentialEffectsBySub(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, repeatAlternate = FALSE)
seqEffectsSummary2 <- getSequentialEffectsBySub(dat=forstmann2008$dat, pp=forstmann2008$pp, repeatAlternate = FALSE)
seqEffectsSummary3 <- getSequentialEffectsBySub(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, repeatAlternate = FALSE)
seqEffectsSummary4 <- getSequentialEffectsBySub(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, repeatAlternate = FALSE)

# By repeat/alternate
seqEffectsSummary1RA <- getSequentialEffectsBySub(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, repeatAlternate = TRUE)
seqEffectsSummary2RA <- getSequentialEffectsBySub(dat=forstmann2008$dat, pp=forstmann2008$pp, repeatAlternate = TRUE)
seqEffectsSummary3RA <- getSequentialEffectsBySub(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, repeatAlternate = TRUE)
seqEffectsSummary4RA <- getSequentialEffectsBySub(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, repeatAlternate = TRUE)


# Plots: Group, by stimulus
for(colors_ in c('black', 'white')) {
  if(colors_ == 'black') {
    pdf(file=file.path(figures_dir, 'SH_MS1RDM_group_stimuluscoded_black.pdf'), width=8, height=5)
    if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
  } else {
    pdf(file=file.path(figures_dir, 'SH_MS1RDM_group_stimuluscoded.pdf'), width=8, height=5)
  }
  par(mfcol=c(2,4), bty='l', oma=c(0,2,0,0))
  plotSequentialEffectsGroup(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, seqEffectsSummary=seqEffectsSummary1,
                             nHistory = 3, main1='Dataset 1', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=forstmann2008$dat, pp=forstmann2008$pp, seqEffectsSummary=seqEffectsSummary2,
                             nHistory = 3, main1='Dataset 2', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, seqEffectsSummary=seqEffectsSummary3,
                             nHistory = 3, main1='Dataset 3', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, seqEffectsSummary=seqEffectsSummary4,
                             nHistory = 3, main1='Dataset 4', full.legend=FALSE)
  dev.off()
}

# new figure for Andrew: 4 rows (RT group, choice group, RT individual, choice individual)
layout(rbind(matrix(1:8, nrow=2), matrix(1:8, nrow=2)+8))
par(bty='l', oma=c(0,2,0,0))
plotSequentialEffectsGroup(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, seqEffectsSummary=seqEffectsSummary1,
                           nHistory = 3, main1='Dataset 1', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
plotSequentialEffectsGroup(dat=forstmann2008$dat, pp=forstmann2008$pp, seqEffectsSummary=seqEffectsSummary2,
                           nHistory = 3, main1='Dataset 2', full.legend=FALSE)
plotSequentialEffectsGroup(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, seqEffectsSummary=seqEffectsSummary3,
                           nHistory = 3, main1='Dataset 3', full.legend=FALSE)
plotSequentialEffectsGroup(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, seqEffectsSummary=seqEffectsSummary4,
                           nHistory = 3, main1='Dataset 4', full.legend=FALSE)
layout(matrix(1:8, nrow=2), widths = c(.15,.25,.5,.2))

#par(mar=c(2,3,2,1), bty='l', mgp=c(2,1,0), oma=c(0,0,0,0))
par(mar=c(2,3,2,.5), mgp=c(2,.75,0), bty='l', oma=c(0,0,0,0))
plotSequentialEffectsIndividual(seqEffectsSummary1RA, main='Dataset 1', xlab2='', ylab2='P(R|RRR)', ylab1='RT(A|RRR) - RT(R|RRR)', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
plotSequentialEffectsIndividual(seqEffectsSummary2RA, main='Dataset 2', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
plotSequentialEffectsIndividual(seqEffectsSummary3RA, main='Dataset 3', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
plotSequentialEffectsIndividual(seqEffectsSummary4RA, main='Dataset 4', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)


m <- matrix(nrow=100, ncol=100)
row1 <- 1:25
row2 <- row1+25
row3 <- row2+25
row4 <- row3+25
column1Upper <- row1
column2Upper <- row2
column3Upper <- row3
column4Upper <- row4
## Upper two rows
m[row1,column1Upper] <- 1
m[row2,column1Upper] <- 2
m[row1,column2Upper] <- 3
m[row2,column2Upper] <- 4
m[row1,column3Upper] <- 5
m[row2,column3Upper] <- 6
m[row1,column4Upper] <- 7
m[row2,column4Upper] <- 8
## Lower two rows
column1Lower <- 1:14
column2Lower <- 15:36
column3Lower <- 37:82
column4Lower <- 83:100
m[row3,column1Lower] <- 9
m[row4,column1Lower] <- 10
m[row3,column2Lower] <- 11
m[row4,column2Lower] <- 12
m[row3,column3Lower] <- 13
m[row4,column3Lower] <- 14
m[row3,column4Lower] <- 15
m[row4,column4Lower] <- 16

for(fn in c('pdf', 'jpeg')) {
  if(fn == 'pdf') pdf(file=file.path(figures_dir, 'SH_MS1RDM_group_and_individual.pdf'), width=8, height=8)
  if(fn == 'jpeg') jpeg(file=file.path(figures_dir, 'SH_MS1RDM_group_and_individual.jpeg'), width=8, height=8, units='in', quality=100, res=500)
  l <- layout(m)
  #layout.show(l)

  par(bty='l', oma=c(0,3,0,0), mgp=c(2,.75,0))
  plotSequentialEffectsGroup(dat=wagenmakers2004_CS$dat, pp=wagenmakers2004_CS$pp, seqEffectsSummary=seqEffectsSummary1,
                             nHistory = 3, main1='Dataset 1', ylab1='A. RT (s)', ylab2='B. Choice proportion', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=forstmann2008$dat, pp=forstmann2008$pp, seqEffectsSummary=seqEffectsSummary2,
                             nHistory = 3, main1='Dataset 2', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=mileticvanmaanen2019exp2block2$dat, pp=mileticvanmaanen2019exp2block2$pp, seqEffectsSummary=seqEffectsSummary3,
                             nHistory = 3, main1='Dataset 3', full.legend=FALSE)
  plotSequentialEffectsGroup(dat=wagenmakers2008exp2$dat, pp=wagenmakers2008exp2$pp, seqEffectsSummary=seqEffectsSummary4,
                             nHistory = 3, main1='Dataset 4', full.legend=FALSE)
  #layout(matrix(1:8, nrow=2), widths = c(.15,.25,.5,.2))

  #par(mar=c(2,3,2,1), bty='l', mgp=c(2,1,0), oma=c(0,0,0,0))
  par(mar=c(2,2,2,1), bty='l', mgp=c(2,.75,0))
  plotSequentialEffectsIndividual(seqEffectsSummary1RA, main='Dataset 1', xlab2='', ylab2='D. P(R|RRR)', ylab1='C. RT(A|RRR) - RT(R|RRR)', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
  plotSequentialEffectsIndividual(seqEffectsSummary2RA, main='Dataset 2', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
  plotSequentialEffectsIndividual(seqEffectsSummary3RA, main='Dataset 3', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  par(mar=c(2,2,2,0.5), bty='l', mgp=c(2,1,0))
  plotSequentialEffectsIndividual(seqEffectsSummary4RA, main='Dataset 4', xlab2='', orderBy='rt'); mtext('Participant', side=1, line=1, cex=.66)
  dev.off()
}



#
# # Group level of NULL (bias), MS3-RDM, MS3-RDM + bias  -----------
# wagenmakers2008exp2null <- loadDataAndPP('wagenmakers2008exp2', 'NULL', 'NULL')
# wagenmakers2008exp2MS3 <- loadDataAndPP('wagenmakers2008exp2', 'zSMuAHbV', 'NULL')
# wagenmakers2008exp2MS3bias <- loadDataAndPP('wagenmakers2008exp2', 'zSMuAHbV_with_baseline', 'NULL')
#
# # By choice
# seqEffectsSummary4null <- getSequentialEffectsBySub(dat=wagenmakers2008exp2null$dat, pp=wagenmakers2008exp2null$pp, repeatAlternate = FALSE)
# seqEffectsSummary4MS3 <- getSequentialEffectsBySub(dat=wagenmakers2008exp2MS3$dat, pp=wagenmakers2008exp2MS3$pp, repeatAlternate = FALSE)
# seqEffectsSummary4MS3bias <- getSequentialEffectsBySub(dat=wagenmakers2008exp2MS3bias$dat, pp=wagenmakers2008exp2MS3bias$pp, repeatAlternate = FALSE)
#
#
# # Plots: Group, by stimulus
# for(colors_ in c('black', 'white')) {
#   if(colors_ == 'black') {
#     pdf(file=file.path(figures_dir, 'stimulus_memory_effects_group_stimuluscoded_wagenmakers2008exp2_comparison_black.pdf'), width=6, height=5)
#     if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
#   } else {
#     pdf(file=file.path(figures_dir, 'stimulus_memory_effects_group_stimuluscoded_wagenmakers2008exp2_comparison.pdf'), width=6, height=5)
#   }
#   par(mfcol=c(2,3), bty='l', oma=c(0,2,0,0))
#   plotSequentialEffectsGroup(dat=wagenmakers2008exp2null$dat, pp=wagenmakers2008exp2null$pp, seqEffectsSummary=seqEffectsSummary4null,
#                              nHistory = 3, main1='Bias only', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
#   plotSequentialEffectsGroup(dat=wagenmakers2008exp2MS3$dat, pp=wagenmakers2008exp2MS3$pp, seqEffectsSummary=seqEffectsSummary4MS3,
#                              nHistory = 3, main1='MS3-RDM', full.legend=FALSE)
#   plotSequentialEffectsGroup(dat=seqEffectsSummary4MS3bias$dat, pp=seqEffectsSummary4MS3bias$pp, seqEffectsSummary=seqEffectsSummary4MS3bias,
#                              nHistory = 3, main1='MS3-RDM + bias', full.legend=FALSE)
#   dev.off()
# }


## Dataset 5, full model
mileticvanmaanen2019exp2block1 <- loadDataSamplersPP('mileticvanmaanen2019exp2block1', 'zSMuAHbV')
# mileticvanmaanen2019exp2block1$pp$rt <- pmin(mileticvanmaanen2019exp2block1$pp$rt, 3) # Censor, like data
seqEffectsSummary3b <- getSequentialEffectsBySub(dat=mileticvanmaanen2019exp2block1$dat, pp=mileticvanmaanen2019exp2block1$pp, repeatAlternate = FALSE)
seqEffectsSummary3bRA <- getSequentialEffectsBySub(dat=mileticvanmaanen2019exp2block1$dat, pp=mileticvanmaanen2019exp2block1$pp, repeatAlternate = TRUE)

for(ftype in c('jpeg','pdf')) {
  fn = file.path(figures_dir, 'SH_MS3RDM_group_and_individual_dataset5.pdf')
  if(ftype == 'pdf') pdf(file=fn, width=6, height=4)
  if(ftype == 'jpeg') jpeg(file=gsub('.pdf', '.jpeg', fn), width=6, height=4, units='in', res=500, quality=100)
  layout(matrix(1:4, nrow=2), widths = c(.4, .6))
  par(bty='l', oma=c(0,2,0,0), cex=0.66, cex.axis=1.0, cex.lab=1)
  plotSequentialEffectsGroup(dat=mileticvanmaanen2019exp2block1$dat, pp=mileticvanmaanen2019exp2block1$pp, seqEffectsSummary=seqEffectsSummary3b,
                             nHistory = 3, main1='A. Dataset 5 (group)', full.legend=FALSE, ylab1='RT (s)', ylab2='Choice proportion')
  par(mar=c(4,3,2,1))
  plotSequentialEffectsIndividual(seqEffectsSummary3bRA, main='B. Dataset 5 (individual)', xlab2='', orderBy='rt', ylab2='P(R|RRR)', ylab1='RT(A|RRR) - RT(R|RRR)'); mtext('Participant', side=1, line=1, cex=par()$cex*par()$cex.axis)
  dev.off()
}
