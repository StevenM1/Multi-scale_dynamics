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


# Objective difficulties
dataset4 <- loadDataSamplersPP('wagenmakers2008exp2', 'zSMuAHbV', pp_conditional = TRUE)
#dataset4 <- get_pp2('wagenmakers2008exp2', 'zSMuAHbV', 'NULL')
dat4 <- prep_data(dataset4$dat)
pp4 <- prep_data(dataset4$pp)
stats4 <- get_difficulty_effects(dat4, pp4)

dataset5 <- loadDataSamplersPP('mileticvanmaanen2019exp2block1', 'zSMuAHbV', pp_conditional = TRUE, samples_dir = './samples_trends')
dat5 <- prep_data(dataset5$dat)
pp5 <- prep_data(dataset5$pp)
stats5 <- get_difficulty_effects(dat5, pp5)

# Subjective difficulties: Use model to determine Q-value of last trial (binned into 5 bins)
wagenmakers2004_CS <- loadDataSamplersPP('wagenmakers2004_CS', 'zSMuAHbV', pp_conditional = TRUE)
forstmann2008 <- loadDataSamplersPP('forstmann2008', 'zSMuAHbV', pp_conditional = TRUE)
mileticvanmaanen2019exp2block2 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', 'zSMuAHbV', pp_conditional = TRUE)
wagenmakers2008exp2 <- loadDataSamplersPP('wagenmakers2008exp2', 'zSMuAHbV', pp_conditional = TRUE)

aggregates1 <- with(wagenmakers2004_CS, get_da_with_Q3bins(samplers, pp))
aggregates2 <- with(forstmann2008, get_da_with_Q3bins(samplers, pp))
aggregates3 <- with(mileticvanmaanen2019exp2block2, get_da_with_Q3bins(samplers, pp))
aggregates4 <- with(wagenmakers2008exp2, get_da_with_Q3bins(samplers, pp))
aggregates5 <- with(dataset5, get_da_with_Q3bins(samplers, pp))

par(mfrow=c(2,4))
plot_subjective_group_effect(aggregates1, main='Dataset 1', mrt=mean(wagenmakers2004_CS$dat$rt))
plot_subjective_group_effect(aggregates2, main='Dataset 2', mrt=mean(forstmann2008$dat$rt))
plot_subjective_group_effect(aggregates3, main='Dataset 3', mrt=mean(mileticvanmaanen2019exp2block2$dat$rt))
plot_subjective_group_effect(aggregates4, main='Dataset 4', mrt=mean(wagenmakers2008exp2$dat$rt))


## For presentation
pdf(file='./figures/FM_MS3RDM_objective.pdf', width=6, height=4)
par(mfrow=c(1,2), bty='l')
palette(c('white', palette()[2:length(palette())]))
par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
plot_group_effect(stats4, mrt=mean(wagenmakers2008exp2$dat$rt), main='Dataset 4', xaxt='n')
axis(side=1, at=1:3, labels=c('High frequency', 'Low frequency', "Very low frequency"))
#mtext('A. Objective',cex=par()$cex*par()$cex.main, font=2, side=2, line = 3)
plot_group_effect(stats5, mrt=mean(dataset5$dat$rt), main='Dataset 5')
dev.off()

pdf(file='./figures/FM_MS3RDM_subjective.pdf', width=8, height=4)
par(mfrow=c(1,4), bty='l')
palette(c('white', palette()[2:length(palette())]))
par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
plot_subjective_group_effect(aggregates1, main='Dataset 1', mrt=mean(wagenmakers2004_CS$dat$rt))
plot_subjective_group_effect(aggregates2, main='Dataset 2', mrt=mean(forstmann2008$dat$rt))
plot_subjective_group_effect(aggregates3, main='Dataset 3', mrt=mean(mileticvanmaanen2019exp2block2$dat$rt))
plot_subjective_group_effect(aggregates4, main='Dataset 4', mrt=mean(wagenmakers2008exp2$dat$rt))
dev.off()





## combine objective (top row) with subjective (bottom row)
m <- matrix(nrow=100, ncol=100)
row1 <- 2:49
row2 <- 50:99
column1Upper <- 1:20
column2Upper <- 20:40
#column5Upper <- 40:42 # some whitspace
column3Upper <- 40:60
column4Upper <- 60:100

## Upper two rows
m[row1,column1Upper] <- 1
m[row1,column2Upper] <- 2
m[row1,column3Upper] <- 3
m[row1,column4Upper] <- 4
# m[row1,column5Upper] <- 13
#m[1,column5Upper] <- 5
#m[100,column5Upper] <- 6

## Lower two row
column1Lower <- 1:20
column2Lower <- 20:40
column3Lower <- 40:60
column4Lower <- 60:80
column5Lower <- 80:100
m[row2,column1Lower] <- 5
m[row2,column2Lower] <- 6
m[row2,column3Lower] <- 7
m[row2,column4Lower] <- 8
m[row2,column5Lower] <- 9

m <- m+2
m[1,c(column1Upper, column2Upper)] <- 1
m[100,c(column1Upper, column2Upper)] <- 1
m[1,c(column3Upper, column4Upper)] <- 2
m[100,c(column3Upper, column4Upper)] <- 2

for(ftype in c('jpeg','pdf')) {
  fn = './figures/FM_MS3RDM_objective_and_subjective.pdf'
  if(ftype == 'pdf') pdf(file=fn, width=8, height=4.5)
  if(ftype == 'jpeg') jpeg(file=gsub('.pdf', '.jpeg', fn), width=8, height=4.5, units='in', res=500, quality=100)
  #pdf(file='./figures/bV_effect_objective_subjective_combined.pdf', width=8, height=4.5)
  l <- layout(m)
  par(bty='l', mar=c(3,4,2,1), mgp=c(2,1,0), oma=c(0,0,1,0))
  #layout.show(l)
  plot.new(); mtext('Group-averaged', cex=par()$cex*par()$cex.main, font=2, line = 2)
  plot.new(); mtext('Individual', cex=par()$cex*par()$cex.main, font=2, line = 2)

  #par(mfrow=c(2,4))
  plot_group_effect(stats4, mrt=mean(wagenmakers2008exp2$dat$rt), main='Dataset 4')
  mtext('A. Objective',cex=par()$cex*par()$cex.main, font=2, side=2, line = 3)
  par(mar=c(3,4,2,2))
  plot_group_effect(stats5, mrt=mean(dataset5$dat$rt), main='Dataset 5')
  par(mar=c(3,5,2,1))
  plot_subject_slopes(stats4, main='Dataset 4')
  par(mar=c(3,4,2,1))
  plot_subject_slopes(stats5, main='Dataset 5')

  plot_subjective_group_effect(aggregates1, main='Dataset 1', mrt=mean(wagenmakers2004_CS$dat$rt))
  mtext('B. Subjective',cex=par()$cex*par()$cex.main, font=2, side=2, line = 3)
  plot_subjective_group_effect(aggregates2, main='Dataset 2', mrt=mean(forstmann2008$dat$rt))
  plot_subjective_group_effect(aggregates3, main='Dataset 3', mrt=mean(mileticvanmaanen2019exp2block2$dat$rt))
  plot_subjective_group_effect(aggregates4, main='Dataset 4', mrt=mean(wagenmakers2008exp2$dat$rt))
  plot_subjective_group_effect(aggregates5, main='Dataset 5', mrt=mean(dataset5$dat$rt))
  dev.off()
}





# # append subjective difficulty
# da1 <- with(wagenmakers2004_CS, get_da_with_Q3bins(samplers, pp2))
# da2 <- with(forstmann2008, get_da_with_Q3bins(samplers, pp2))
# da3 <- with(mileticvanmaanen2019exp2block2, get_da_with_Q3bins(samplers, pp2))
# da4 <- with(wagenmakers2008exp2, get_da_with_Q3bins(samplers, pp2))
#
#
# par(mfrow=c(2,4))
# plot_group_effect(stats5)
# plot_group_effect(stats4)
# plot_subject_slopes(stats5)
# plot_subject_slopes(stats4)
#
#
#
#
#
#
# # Subjective difficulty ---------------------------------------------------
# wagenmakers2004_CS <- get_pp2('wagenmakers2004_CS', 'zSMuAHbV', 'NULL')
# forstmann2008 <- get_pp2('forstmann2008', 'zSMuAHbV', 'NULL')
# mileticvanmaanen2019exp2block2 <- get_pp2('mileticvanmaanen2019exp2block2', 'zSMuAHbV', 'NULL')
# wagenmakers2008exp2 <- get_pp2('wagenmakers2008exp2', 'zSMuAHbV', 'NULL')
#
# # append subjective difficulty
# da1 <- with(wagenmakers2004_CS, get_da_with_Q3bins(samplers, pp2))
# da2 <- with(forstmann2008, get_da_with_Q3bins(samplers, pp2))
# da3 <- with(mileticvanmaanen2019exp2block2, get_da_with_Q3bins(samplers, pp2))
# da4 <- with(wagenmakers2008exp2, get_da_with_Q3bins(samplers, pp2))
#
# # Flip q3p5prevbin to match difficulty
# da1$q3p5prevbin <- 6-da1$q3p5prevbin
# da2$q3p5prevbin <- 6-da2$q3p5prevbin
# da3$q3p5prevbin <- 6-da3$q3p5prevbin
# da4$q3p5prevbin <- 6-da4$q3p5prevbin
#
# allstats1b <- get_stats(wagenmakers2004_CS[['pp2']], da1, flip = TRUE)
# allstats2b <- get_stats(forstmann2008[['pp2']], da2, flip=TRUE)
# allstats3b <- get_stats(mileticvanmaanen2019exp2block2[['pp2']], da3, flip=TRUE)
# allstats4b <- get_stats(wagenmakers2008exp2[['pp2']], da4, flip=TRUE)
#
#
# mrt1 <- mean(aggregate(rt~subjects,da1,mean)[,2])
# mrt2 <- mean(aggregate(rt~subjects,da2,mean)[,2])
# mrt3 <- mean(aggregate(rt~subjects,da3,mean)[,2])
# mrt4 <- mean(aggregate(rt~subjects,da4,mean)[,2])
#
#
# #pdf(file='./figures/bV_effect_rtonly.pdf', width=8, height=2.5)
# par(mfrow=c(2,4), bty='l', mar=c(3,4,2,1), mgp=c(2,1,0), oma=c(0,1,0,0))
# plot_q3_effect(allstats1b, mrt=mrt1, main1='Dataset 1', xlab1='Subjective difficulty', ylab1='')
# mtext(side=2, text='RT (s)', line=2.5, cex=par()$cex, las=0)
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
# plot_q3_effect(allstats2b, mrt=mrt2, main1='Dataset 2', xlab1='Subjective difficulty', ylab1='')
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
# plot_q3_effect(allstats3b, mrt=mrt3, main1='Dataset 3', xlab1='Subjective difficulty', ylab1='')
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
# plot_q3_effect(allstats4b, mrt=mrt4, main1='Dataset 4', xlab1='Subjective difficulty', ylab1='')
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
# #dev.off()
#
#
#
# plotslopes <- function(stats, ...) {
#   ordering <- order(stats$slopes_data)
#   plot(stats$slopes_data[ordering], ylim=range(c(stats$slopes_ci, stats$slopes_data)), pch=4, lwd=2, xaxt='n', ...)
#   abline(h=0, lty=2)
#   arrows(x0=1:ncol(stats$slopes_ci), y0=stats$slopes_ci[1,ordering], y1=stats$slopes_ci[2,ordering], code=3, angle=90, length=0.05, col=2, lwd=2)
# }
#
# plotslopes(allstats1b, xlab='Participant', ylab='Slope')
# plotslopes(allstats2b, xlab='Participant', ylab='')
# plotslopes(allstats3b, xlab='Participant', ylab='')
# plotslopes(allstats4b, xlab='Participant', ylab='')
#
#
#
#
#
#
#
#
#
# # Objective difficulty ---------------------------------------------------
# #wagenmakers2008exp2 <- get_pp2('wagenmakers2008exp2', 'zSMuAHbV', 'NULL')
# #da4 <- with(wagenmakers2008exp2, get_da_with_Q3bins(samplers, pp2))
# da4_objective <- da4
# da4_objective$D <- Lag(da4_objective$W, 1)
# da4_objective <- da4_objective[!is.na(da4_objective$D),]
# da4_objective$D <- factor(da4_objective$D, levels=c('hf', 'lf', 'vlf', 'nonword')) #levels=c('vlf', 'lf', 'hf', 'nonword'))
# # remove nonwords -- difficulty is undetermined here
# da4_objective <- droplevels(da4_objective[da4_objective$D!='nonword',])
# allstats4_objective <- get_stats(wagenmakers2008exp2[['pp2']], da4_objective, flip=TRUE)
# # get rid of nonwords
# # allstats4_objective <- lapply(allstats4_objective, function(x) {
# #   if('D' %in% colnames(x)) return(x[1:3,]) else return(x) })
#
# # Dataset 5
# # mvmblock1 <- get_pp2('mileticvanmaanen2019exp2block1', 'zSMuAHbV', 'NULL')
# # mvmblock1dat <- mvmblock1$dat
# # mvmblock1dat$difficulty_prev <- Lag(mvmblock1dat$difficulty, 1)
# # aggData <- aggregate(rt~difficulty_prev, aggregate(rt~difficulty_prev*subjects,mvmblock1dat,mean), mean)
# #
# # ## pp
# # mvmblock1pp2 <- mvmblock1$pp2
# # mvmblock1pp2$difficulty_prev <- Lag(mvmblock1pp2$difficulty, 1)
# # aggByPostN <- aggregate(rt~difficulty_prev*postn, aggregate(rt~difficulty_prev*subjects*postn,mvmblock1pp2,mean), mean)
# # aggPP2 <- aggregate(rt~difficulty_prev, aggByPostN, quantile, c(0.025, 0.5, .975))
# #
# # mvmblock1pp <- mvmblock1$pp
# # mvmblock1pp$difficulty_prev <- Lag(mvmblock1pp$difficulty, 1)
# # aggByPostN <- aggregate(rt~difficulty_prev*postn, aggregate(rt~difficulty_prev*subjects*postn,mvmblock1pp,mean), mean)
# # aggPP <- aggregate(rt~difficulty_prev, aggByPostN, quantile, c(0.025, 0.5, .975))
#
#
#
# # group effects on subjective value
#
#
#
# ## what about dataset 4?
# dataset4 <- get_pp2('wagenmakers2008exp2', 'zSMuAHbV', 'NULL')
# dat4 <- prep_data(dataset4$dat)
# pp4 <- prep_data(dataset4$pp2)
# stats4 <- get_difficulty_effects(dat4, pp4)
#
# stats <- stats4
# #xs <- 1:3
# ylim <- range(c(stats$aggData$rt, stats$aggPp$rt))
# plot(stats$aggData$difficulty_prev, stats$aggData$rt, pch=4, lwd=2, xlab='Difficulty', ylab='RT (s)', ylim=ylim)
# arrows(x0=stats$aggPp$difficulty_prev, y0=stats$aggPp$rt[,1], y1=stats$aggPp$rt[,3], code=3, angle=90, length=0.025, col=2, lwd=2)
# lines(x=stats$aggPp$difficulty_prev, stats$aggPp$rt[,2], lwd=2, col=2)
#
# # Plot 2: individual differences (objective)
# coefficients <- stats$coefficientsBySubjectData[stats$coefficientsBySubjectData$parameter=='difficulty_prev',]
# ordering <- order(coefficients$coefficient)
# plot(1:nrow(coefficients), coefficients[ordering, 'coefficient'], pch=4, ylab='Slope', xlab='Subject', lwd=2)
# abline(h=0)
#
# coefficientsPp <- stats$coefficientsBySubjectPp[stats$coefficientsBySubjectPp$parameter=='difficulty_prev',]
# arrows(x0=1:nrow(coefficients),
#        y0=coefficientsPp$coefficient[ordering,1],
#        y1=coefficientsPp$coefficient[ordering,3],
#        code=3, angle=90, length=0.025, col=2, lwd=2
#        )
# lines(1:nrow(coefficients), coefficientsPp$coefficient[ordering,2], lwd=2,col=2)
#
#
#
#
# # Group-level
# plot(as.numeric(aggData$difficulty_prev)[5:1], aggData$rt, ylim=range(aggData$rt, aggPP$rt, aggPP2$rt), pch=4)
# arrows(x0=as.numeric(aggPP2$difficulty_prev)[5:1],
#        y0=aggPP2$rt[,1], y1=aggPP2$rt[,3], angle=90, code=3, length=0.05)
#
# arrows(x0=as.numeric(aggPP$difficulty_prev)[5:1],
#        y0=aggPP$rt[,1], y1=aggPP$rt[,3], angle=90, code=3, length=0.05, col=2, lwd=2)
# legend('bottomright', c('Conditional pp', 'Unconditional pp'), col=1:2, lty=c(1,1),bty='n')
#
#
# ## what about per subject?
# mvmblock1dat$subjects <- as.numeric(mvmblock1dat$subjects)
# mvmblock1pp$subjects <- as.numeric(mvmblock1pp$subjects)
# mvmblock1pp2$subjects <- as.numeric(mvmblock1pp2$subjects)
# mvmblock1dat$difficulty_prev <- as.numeric(mvmblock1dat$difficulty_prev)
# mvmblock1pp$difficulty_prev <- as.numeric(mvmblock1pp$difficulty_prev)
# mvmblock1pp2$difficulty_prev <- as.numeric(mvmblock1pp2$difficulty_prev)
#
# aggBySubData <- aggregate(rt~difficulty_prev*subjects,mvmblock1dat,mean)
# # get slopes by subject
# coeffs <- sapply(sort(unique(mvmblock1dat$subjects)), function(subject) lm(rt~as.numeric(difficulty_prev), mvmblock1dat[mvmblock1dat$subjects==subject,])$coef)
#
# aggBySubByPostN <- aggregate(rt~difficulty_prev*subjects*postn,mvmblock1pp,mean)
# coeffsPP <- parallel::mclapply(1:100, function(postn) {
#   sapply(sort(unique(mvmblock1dat$subjects)), function(subject) lm(rt~as.numeric(difficulty_prev), mvmblock1pp[mvmblock1pp$subjects==subject&mvmblock1pp$postn==postn,])$coef)
# }, mc.cores=20)
# aggBySub <- aggregate(rt~difficulty_prev*subjects, aggBySubByPostN, quantile, c(0.025, 0.5, .975))
#
# par(mfrow=c(2,3))
# for(subject in 11:16) { #length(unique(aggBySubData$subject))) {
# #  subjectName <- unique(aggBySubData$subject)[[subject]]
#   idx <- aggBySubData$subjects==subject
#   idx2 <- aggBySub$subjects==subject
#   plot(as.numeric(aggBySubData[idx,'difficulty_prev']),aggBySubData[idx,'rt'], ylim=range(c(aggBySubData[idx,'rt'], aggBySub[idx2,]$rt)), main=subject)
#   arrows(x0=as.numeric(aggBySub[idx2,'difficulty_prev']),
#          y0=aggBySub[idx2,]$rt[,1], y1=aggBySub[idx2,]$rt[,3], angle=90, code=3, length=0.05)
#
#   ## add slopes
#   abline(a=coeffs[1,subject], b=coeffs[2,subject])
#
#   ## pp
#   for(postn in 1:100) abline(a=coeffsPP[[postn]][1,subject],b=coeffsPP[[postn]][2,subject], col=adjustcolor(2, alpha.f=.3))
#   ## mean slope?
# }
#
# # plot slopes against slopes
# ordering <- order(coeffs[2,])
# plot(1:ncol(coeffs), coeffs[2,ordering], pch=4); abline(h=0)
# # for(postn in 1:100) points(1:ncol(coeffs), coeffsPP[[postn]][2,ordering], col=adjustcolor(2,alpha.f=.1))
#
# all_slopes <- sapply(1:100, function(x) coeffsPP[[x]][2,])
# all_slopes_CI <- apply(all_slopes, 1, quantile, c(0.025, 0.5, .975))
# arrows(x0=1:ncol(coeffs), y0=all_slopes_CI[1,ordering], y1=all_slopes_CI[3,ordering], code=3, angle=90, length=0.05, col=2)
# lines(1:ncol(coeffs), all_slopes_CI[2,ordering], col=2)
#
#
# #da5_objective <- with(mvmblock1, get_da_with_Q3bins(samplers, pp2))
# da5_objective$D <- Lag(da5_objective$difficulty, 1)
# da5_objective$D <- factor(da5_objective$D, levels=c('D5', 'D4', 'D3', 'D2', 'D1'))
# da5_objective <- da5_objective[!is.na(da5_objective$D),]
# allstats5 <- get_stats(mvmblock1[['pp2']], da5_objective, flip=TRUE)
#
# dat2 <- do.call(rbind, mvmblock1$samplers[[1]]$data)
# dat2 <- dat2[dat2$lM==TRUE,]
# mrt5 <- mean(aggregate(rt~subjects,dat2,mean)[,2])
#
#
# #dev.off()
# #pdf(file='./figures/bV_effect_objectivediff_rtonly.pdf', width=4.5, height=2.5)
# par(mfcol=c(1,2), bty='l', mar=c(3,4,2,1), mgp=c(2,1,0))
# plot_q3_effect(allstats4_objective, mrt=mrt4, main1='Dataset 4', xlab1='Word frequency', ylab1='')
# mtext(side=2, text='RT (s)', line=3, cex=par()$cex, las=0)
# axis(side=1, at=3:1, labels=c('Very low', 'Low', 'High'))
# par(mar=c(3,3,2,1))
# plot_q3_effect(allstats5, mrt=mrt5, main1='Dataset 5', xlab1='Difficulty', ylab='')
# axis(side=1, at=5:1, labels=c('Hard', '', '', '', 'Easy'))
# #dev.off()
#
#
# ## Combine in 6 panels (4 times subjective, 2 times objective)
# m <- matrix(nrow=100, ncol=100)
# row1 <- 1:50
# row2 <- row1+50
# column1Upper <- 1:round(100/6)
# column2Upper <- (1+round(100/6*1)):round(100/6*2)
# column3Upper <- (1+round(100/6*2)):round(100/6*3)
# column4Upper <- (1+round(100/6*3)):round(100/6*4)
# column5Upper <- (1+round(100/6*4)):round(100/6*5)
# column6Upper <- (1+round(100/6*5)):round(100/6*6)
# ## Upper two rows
# m[row1,column1Upper] <- 1
# m[row1,column2Upper] <- 2
# m[row1,column3Upper] <- 3
# m[row1,column4Upper] <- 4
# m[row1,column5Upper] <- 5
# m[row1,column6Upper] <- 6
# ## Lower two rows
# column1Lower <- 1:14
# column2Lower <- 15:36
# column3Lower <- 37:82
# column4Lower <- 83:100
# m[row2,column1Lower] <- 7
# m[row2,column2Lower] <- 8
# m[row2,column3Lower] <- 9
# m[row2,column4Lower] <- 10
#
#
# pdf(file='./figures/bV_effect_group_and_individual.pdf', width=8, height=5)
# l <- layout(m)
# #layout.show(l)
#
# par(bty='l', mar=c(3,2,2,.5), mgp=c(2,1,0), oma=c(0,1,0,1))
# plot_q3_effect(allstats1b, mrt=mrt1, main1='Dataset 1', xlab1='Subjective difficulty', ylab1='')
# mtext(side=2, text='RT (s)', line=2, cex=par()$cex, las=0)
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
# plot_q3_effect(allstats2b, mrt=mrt2, main1='Dataset 2', xlab1='Subjective difficulty', ylab1='')
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
# plot_q3_effect(allstats3b, mrt=mrt3, main1='Dataset 3', xlab1='Subjective difficulty', ylab1='')
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
# plot_q3_effect(allstats4b, mrt=mrt4, main1='Dataset 4', xlab1='Subjective difficulty', ylab1='')
# axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
#
# plot_q3_effect(allstats4_objective, mrt=mrt4, main1='Dataset 4', xlab1='Word frequency', ylab1='')
# #mtext(side=2, text='RT (s)', line=3, cex=par()$cex, las=0)
# axis(side=1, at=3:1, labels=c('Very low', 'Low', 'High'))
# #par(mar=c(3,3,2,1))
# plot_q3_effect(allstats5, mrt=mrt5, main1='Dataset 5', xlab1='Difficulty', ylab='')
# axis(side=1, at=5:1, labels=c('Hard', '', '', '', 'Easy'))
#
# ## Next, by subject
# par(mgp=c(1,1,0))
# plotslopes(allstats1b, xlab='Participant', ylab='', main='Dataset 1')
# mtext(side=2, text='Slope', line=2, cex=par()$cex, las=0)
# plotslopes(allstats2b, xlab='Participant', ylab='', main='Dataset 2')
# plotslopes(allstats3b, xlab='Participant', ylab='', main='Dataset 3')
# plotslopes(allstats4b, xlab='Participant', ylab='', main='Dataset 4')
# dev.off()
#
#
# plotslopes(allstats4_objective, xlab='Participant', ylab='')
# plotslopes(allstats5, xlab='Participant', ylab='')
#
#
#
# # Combine in 5-panel plot -------------------------------------------------
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
# for(colors_ in c('black', 'white')) {
#   if(colors_ == 'black') {
#     pdf(file='./figures/bV_effect_rtonly_5panels_black.pdf', width=8, height=2.5)
#     #pdf(file=file.path(figures_dir, 'stimulus_memory_effects_group_stimuluscoded_black.pdf'), width=8, height=5)
#     if(colors_ == 'black') par(bg = 'black', fg='white', col='white', col.axis='white', col.main='white', col.lab='white', col.sub='white')
#   } else {
#     pdf(file='./figures/bV_effect_rtonly_5panels.pdf', width=8, height=2.5)
#   }
#
#   par(mfcol=c(1,5), bty='l', mar=c(3,3,2,.5), mgp=c(2,1,0), oma=c(0,1,0,0))
#   plot_q3_effect(allstats1, mrt=mrt1, main1='Dataset 1', xlab1='Subjective difficulty', ylab1='')
#   mtext(side=2, text='RT (s)', line=2.5, cex=par()$cex, las=0)
#   axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
#   plot_q3_effect(allstats2, mrt=mrt2, main1='Dataset 2', xlab1='Subjective difficulty', ylab1='')
#   axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
#   plot_q3_effect(allstats3, mrt=mrt3, main1='Dataset 3', xlab1='Subjective difficulty', ylab1='')
#   axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
#   plot_q3_effect(allstats4, mrt=mrt4, main1='Dataset 4', xlab1='Subjective difficulty', ylab1='')
#   axis(1, 1:5, labels=c('Easy', '', '', '', 'Hard'))
#   par(mar=c(3,3,2,1))
#   plot_q3_effect(allstats4_objective, mrt=mrt4, main1='Dataset 4', xlab1='Objective difficulty', ylab1='')
#   axis(side=1, at=1:3, labels=c('Easy', '', 'Hard'))
#   #axis(side=1, at=3:1, labels=c('Very low', 'Low', 'High'))
#   dev.off()
# }
#
#
#
#
# # For supplement: Dataset 5 -----------------------------------------------
# pdf(file='./figures/supplementary_figures/RS_dataset5.pdf', width=5, height=2.5)
# par(mfcol=c(1,2), bty='l', mar=c(3,4,2,1), mgp=c(2,1,0), las=1)
# plot_q3_effect(allstats4_objective, mrt=mrt4, main1='Dataset 4', xlab1='Difficulty', ylab1='')
# mtext(side=2, text='RT (s)', line=3, cex=par()$cex, las=0)
# axis(side=1, at=3:1, labels=c('Hard', '', 'Easy'))
# par(mar=c(3,3,2,1))
# plot_q3_effect(allstats5, mrt=mrt5, main1='Dataset 5', xlab1='Difficulty', ylab='')
# axis(side=1, at=5:1, labels=c('Hard', '', '', '', 'Easy'))
# dev.off()

