rm(list=ls())
library(EMC2)
library(emcAdapt)
library(parallel)
source('./extra_EMC2_functions/adaptive.R')
source('./extra_EMC2_functions/make_data.R')
source('./extra_EMC2_functions/model_RDMdynamic.R')
source('./extra_EMC2_functions/utils.R')
figures_dir = './figures'

# Load all ----------------------------------------------------------------
merge_dat_null_ms3 <- function(task, pp_ms3=NULL, ...) {
  ds_null <- loadDataSamplersPP(task, 'NULL', pp_conditional = TRUE)
  pp_null <- ds_null$pp
  dat_null <- ds_null$dat
  if(is.null(pp_ms3)) {
    ds_ms3 <- loadDataSamplersPP(task, learningModel='zSMuAHbV', pp_conditional=TRUE, ...)
    pp_ms3 <- ds_ms3$pp
  }

  ## trialwise mean RT prediction
  agg_iid <- aggregate(rt~trials*subjects, pp_null, mean)
  agg_ms3 <- aggregate(rt~trials*subjects, pp_ms3, mean)

  agg_combined <- merge(agg_iid, agg_ms3,by=c('subjects','trials'), suffixes=c('_iid', '_ms3'))
  merged <- merge(dat_null, agg_combined, by=c('subjects','trials'))
  merged
}

ds1 <- merge_dat_null_ms3('wagenmakers2004_CS')
ds2 <- merge_dat_null_ms3('forstmann2008')
ds3 <- merge_dat_null_ms3('mileticvanmaanen2019exp2block2')
ds4 <- merge_dat_null_ms3('wagenmakers2008exp2')

get_var_explained_ms3 <- function(x) cor(x[,'rt_ms3'], x[,'rt'])^2
r2_ms3 <- sapply(list('ds1'=ds1, 'ds2'=ds2, 'ds3'=ds3, 'ds4'=ds4), function(ds) mean(sapply(unique(ds$subjects), function(x) get_var_explained_ms3(ds[ds$subjects==x,]))))
r2_ms3_se <- sapply(list('ds1'=ds1, 'ds2'=ds2, 'ds3'=ds3, 'ds4'=ds4), function(ds) sd(sapply(unique(ds$subjects), function(x) get_var_explained_ms3(ds[ds$subjects==x,])))/sqrt(length(unique(ds$subjects))))

get_var_explained_iid <- function(x) cor(x[,'rt_iid'], x[,'rt'])^2
r2_iid <- sapply(list('ds1'=ds1, 'ds2'=ds2, 'ds3'=ds3, 'ds4'=ds4), function(ds) mean(sapply(unique(ds$subjects), function(x) get_var_explained_iid(ds[ds$subjects==x,]))))
r2_iid_se <- sapply(list('ds1'=ds1, 'ds2'=ds2, 'ds3'=ds3, 'ds4'=ds4), function(ds) sd(sapply(unique(ds$subjects), function(x) get_var_explained_iid(ds[ds$subjects==x,])))/sqrt(length(unique(ds$subjects))))


## Trend models
ds1T <- merge_dat_null_ms3('wagenmakers2004_CS', trendModel='DCT', trendPar='B', nTrendPars=3)
ds2T <- merge_dat_null_ms3('forstmann2008',trendModel='DCT', trendPar='v', nTrendPars=3)
ds3T <- merge_dat_null_ms3('mileticvanmaanen2019exp2block2',trendModel='DCT', trendPar='v', nTrendPars=3)
ds4T <- merge_dat_null_ms3('wagenmakers2008exp2',trendModel='DCT', trendPar='B', nTrendPars=3)

r2_ms3T <- sapply(list('ds1'=ds1T, 'ds2'=ds2T, 'ds3'=ds3T, 'ds4'=ds4T), function(ds) mean(sapply(unique(ds$subjects), function(x) get_var_explained_ms3(ds[ds$subjects==x,]))))
r2_ms3T_se <- sapply(list('ds1'=ds1T, 'ds2'=ds2T, 'ds3'=ds3T, 'ds4'=ds4T), function(ds) sd(sapply(unique(ds$subjects), function(x) get_var_explained_ms3(ds[ds$subjects==x,])))/sqrt(length(unique(ds$subjects))))



## Simulate RT data with two mechanisms turned off
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMsNewEMC')
source(file.path(root_dir, 'extra_EMC2_functions/adaptive.R'))
source(file.path(root_dir, 'extra_EMC2_functions/model_RDMdynamic.R'))
source(file.path(root_dir, 'extra_EMC2_functions/utils.R'))
source(file.path(root_dir, 'extra_EMC2_functions/make_data.R'))

ds1 <- loadDataSamplersPP('wagenmakers2004_CS', 'zSMuAHbV')
ds2 <- loadDataSamplersPP('forstmann2008', 'zSMuAHbV')
ds3 <- loadDataSamplersPP('mileticvanmaanen2019exp2block2', 'zSMuAHbV')
ds4 <- loadDataSamplersPP('wagenmakers2008exp2', 'zSMuAHbV')

# get_pps <- function(ds, n_post=100, n_cores=1) {
#   pp_SMonly <- predict(ds$samplers, n_post = n_post, n_cores=n_cores, collapse_across_AM=TRUE, collapse_across_FM=TRUE, conditionalOnData=TRUE)
#   pp_AMonly <- predict(ds$samplers, n_post = n_post, n_cores=n_cores, collapse_across_SM=TRUE, collapse_across_FM=TRUE, conditionalOnData=TRUE)
#   pp_FMonly <- predict(ds$samplers, n_post = n_post, n_cores=n_cores, collapse_across_SM=TRUE, collapse_across_AM=TRUE, conditionalOnData=TRUE)
#   return(list(SMonly=pp_SMonly, AMonly=pp_AMonly, FMonly=pp_FMonly))
# }
# # this takes a while
# all_pps <- lapply(list(ds1, ds2, ds3, ds4), get_pps, n_post=100, n_cores=10)
# # save for future reference
# for(i in 1:4) {
#   dsname <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')[i]
#   this_pp <- all_pps[[i]]
#   smonly <- this_pp[[1]]
#   amonly <- this_pp[[2]]
#   fmonly <- this_pp[[3]]
#   save(smonly, file=paste0('./pp_singlemechanism/', dsname, '-model-SMonly_pp-conditional.RData'))
#   save(amonly, file=paste0('./pp_singlemechanism/', dsname, '-model-AMonly_pp-conditional.RData'))
#   save(fmonly, file=paste0('./pp_singlemechanism/', dsname, '-model-FMonly_pp-conditional.RData'))
# }

## load
all_pps <- vector(length=4, mode='list')
for(i in 1:4) {
  dsname <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')[i]
  all_pps[[i]]$SMonly <- EMC2:::loadRData(paste0('./pp_singlemechanism/', dsname, '-model-SMonly_pp-conditional.RData'))
  all_pps[[i]]$AMonly <- EMC2:::loadRData(paste0('./pp_singlemechanism/', dsname, '-model-AMonly_pp-conditional.RData'))
  all_pps[[i]]$FMonly <- EMC2:::loadRData(paste0('./pp_singlemechanism/', dsname, '-model-FMonly_pp-conditional.RData'))
}

get_var_explained_ms3 <- function(x) cor(x[,'rt_ms3'], x[,'rt'])^2

# Merge with data, find R2 again of single mechanism
null_preds_merged_ds1 <- lapply(all_pps[[1]], function(x) merge_dat_null_ms3('wagenmakers2004_CS', x))
out1 <- lapply(null_preds_merged_ds1, function(y)
  mean(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))
)
out1_se <- lapply(null_preds_merged_ds1, function(y)
  sd(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))/sqrt(length(unique(y$subjects)))
)


null_preds_merged_ds2 <- lapply(all_pps[[2]], function(x) merge_dat_null_ms3('forstmann2008', x))
out2 <- lapply(null_preds_merged_ds2, function(y)
  mean(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))
)
out2_se <- lapply(null_preds_merged_ds2, function(y)
  sd(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))/sqrt(length(unique(y$subjects)))
)

null_preds_merged_ds3 <- lapply(all_pps[[3]], function(x) merge_dat_null_ms3('mileticvanmaanen2019exp2block2', x))
out3 <- lapply(null_preds_merged_ds3, function(y)
  mean(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))
)
out3_se <- lapply(null_preds_merged_ds3, function(y)
  sd(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))/sqrt(length(unique(y$subjects)))
)

null_preds_merged_ds4 <- lapply(all_pps[[4]], function(x) merge_dat_null_ms3('wagenmakers2008exp2', x))
out4 <- lapply(null_preds_merged_ds4, function(y)
  mean(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))
)
out4_se <- lapply(null_preds_merged_ds4, function(y)
  sd(sapply(unique(y$subjects), function(x) get_var_explained_ms3(y[y$subjects==x,])))/sqrt(length(unique(y$subjects)))
)

## make a nice table
tabl <- data.frame(rbind(unlist(out1), unlist(out2), unlist(out3), unlist(out4)))
tabl$iid <- r2_iid
tabl$ms3 <- r2_ms3
tabl$ms3T <- r2_ms3T
tabl <- tabl[,c('iid','SMonly','AMonly', 'FMonly', 'ms3', 'ms3T')]

## make SEs
tabl_se <- data.frame(rbind(unlist(out1_se), unlist(out2_se), unlist(out3_se), unlist(out4_se)))
tabl_se$iid <- r2_iid_se
tabl_se$ms3 <- r2_ms3_se
tabl_se$ms3T <- r2_ms3T_se
tabl_se <- tabl_se[,c('iid','SMonly','AMonly', 'FMonly', 'ms3', 'ms3T')]

# combine
tabl_final <- tabl
for(row in 1:nrow(tabl)) {
  for(col in 1:ncol(tabl)) {
    tabl_final[row, col] <- paste0(round(tabl[row,col],3), ' (', round(tabl_se[row,col], 3), ')')
  }
}
colnames(tabl_final) <- c('IID-RDM', 'SM only', 'AM only', 'FM only', 'MS3-RDM', 'MS3-RDM + Trends')
row.names(tabl_final) <- paste0('Dataset ', 1:4)
write.xlsx(tabl_final, file='./tables/R2_values.xlsx', colNames = TRUE, rowNames = TRUE)




# Latent: Coefficient of variation for u, B, z ------------------------------------------------------------------
get_sds_cvs <- function(task, learningModel='zSMuAHbV') {
  ds1_ms3 <- loadDataSamplersPP(task, learningModel=learningModel, pp_conditional=TRUE)
  npars <- data.frame(attr(ds1_ms3$pp, 'npars'))
  npars$dB <- 2*npars$q01*npars$weight1  # 2x because this is added to one accumulator's threshold and subtracted from the other's
  npars$u <- (1-npars$q02)*npars$weight2
  npars$B_ <- npars$q03*npars$weight3

  is_lR1 = npars$lR==levels(npars$lR)[1]
  if(learningModel == 'NULL') {
    # get dB from variability
    npars$dB <- abs(npars$B[is_lR1]-npars$B[!is_lR1])
  }

  npars$vmean <- rep((npars$v[is_lR1]+npars$v[!is_lR1])/2, each=2)
  npars$Bmean <- rep((npars$B[is_lR1]+npars$B[!is_lR1])/2, each=2)

  sds_per_postn <- aggregate(cbind(dB,u,B_,Bsd=Bmean)~postn*subjects,npars[is_lR1,],sd)
  means_per_postn <- aggregate(cbind(vmean,Bmean=Bmean)~postn*subjects,npars[is_lR1,],mean)
  sds_means <- cbind(sds_per_postn, means_per_postn[,3:4])
  sds_means$dB_cv <- sds_means$dB/sds_means$vmean
  sds_means$u_cv <- sds_means$u/sds_means$vmean
  sds_means$B_cv <- sds_means$B_/sds_means$Bmean

  mean_sd_cv_by_sub <- aggregate(.~subjects, sds_means, mean)
  return(mean_sd_cv_by_sub)
}

ds1 <- get_sds_cvs('wagenmakers2004_CS')
ds2 <- get_sds_cvs('forstmann2008')
ds3 <- get_sds_cvs('mileticvanmaanen2019exp2block2')
ds4 <- get_sds_cvs('wagenmakers2008exp2')

apply(ds1[,2:ncol(ds1)], 2, mean)
apply(ds2[,2:ncol(ds2)], 2, mean)
apply(ds3[,2:ncol(ds3)], 2, mean)
apply(ds4[,2:ncol(ds4)], 2, mean)


## Null models
ds2iid <- get_sds_cvs('forstmann2008', learningModel='NULL')
ds4iid <- get_sds_cvs('wagenmakers2008exp2', learningModel='NULL')
apply(ds2iid[,2:ncol(ds2iid)], 2, mean)
apply(ds4iid[,2:ncol(ds4iid)], 2, mean)



# What is the accuracy of participants in DS4 if we turn off SM? ----------
ds4 <- loadDataSamplersPP('wagenmakers2008exp2', 'zSMuAHbV', pp_conditional = TRUE)
samplers4 <- ds4$samplers
for(samplern in c(1)) {
  for(subjectn in 1:length(samplers4[[1]]$data)) {
    # remove stim_identifier, just a trick to generate data with the Q-SM values
    # meaned across the entire experiment instead of block-wise
    colnames(samplers4[[samplern]]$data[[subjectn]])[4] <- 'idf'
  }
}
# debug(update_pars)
pp_noSM <- predict(samplers4, n_post = 100, n_cores=20,
                   collapse_across_SM=TRUE, conditionalOnData=TRUE)
pp_noSM$accuracy <- pp_noSM$R==pp_noSM$S

## pp of NULL model without any threshold bias
samplers_noB <- EMC2:::loadRData('./samples-ds4-nobias.RData')
pp_noB <- predict(samplers_noB, n_cores=20)

## true pp
pp_withSM <- EMC2:::loadRData('./posteriorpredictives/wagenmakers2008exp2_model-zSMuAHbV_trend-NULL_pp-conditional.RData')
mean(pp_withSM$accuracy)

## pp of IID-RDM
samplers_iidrdm <- loadDataSamplersPP('wagenmakers2008exp2', 'NULL', pp_conditional = TRUE)
pp_iidrdm <- samplers_iidrdm$pp


# Plot these
dat <- EMC2:::loadRData('./datasets/wagenmakers2008exp2.RData')
dat$choice <- dat$R == levels(dat$R)[1]
pp_withSM$choice <- pp_withSM$R==dat[dat$choice==1,'R'][1]
pp_noB$choice <- pp_noB$R==dat[dat$choice==1,'R'][1]
pp_noSM$choice <- pp_noSM$R==dat[dat$choice==1,'R'][1]
pp_iidrdm$choice <- pp_iidrdm$R==dat[dat$choice==1,'R'][1]

get_accuracy_rt_by_condition <- function(dat=dat,pp=pp) {
  dat$pattern <- dat$condition
  pp$pattern <- pp$condition
  dataBySub <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*subjects, dat, mean)
  dataMeans <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, dataBySub, mean)
  dataSEs <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, dataBySub, function(x) sd(x)/sqrt(length(x)))
  dataBySub$pattern <- factor(dataBySub$pattern, levels=levels(dat$pattern))
  dataMeans$pattern <- factor(dataMeans$pattern, levels=levels(dat$pattern))
  dataSEs$pattern <- factor(dataSEs$pattern, levels=levels(dat$pattern))

  # Split data also by choice
  dataBySubByChoice <- aggregate(cbind(rt,accuracy)~pattern*subjects*choice, dat, mean)
  dataBySubByChoice$pattern <- factor(dataBySubByChoice$pattern, levels=levels(dat$pattern))
  dataMeansByChoice <- aggregate(cbind(rt,accuracy)~pattern*choice, dataBySubByChoice, mean)
  dataMeansByChoice$pattern <- factor(dataMeansByChoice$pattern, levels=levels(dat$pattern))

  # Next, the posterior predictives
  # Group-level, across choices
  ppByPostN <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*subjects*postn, pp, mean)
  ppByPostN <- aggregate(cbind(rt,accuracy,choice,stim)~pattern*postn, ppByPostN, mean)
  ppQuantiles <- aggregate(cbind(rt,accuracy,choice,stim)~pattern, ppByPostN, quantile, c(.025, .5, .975))
  ppQuantiles$pattern <- factor(ppQuantiles$pattern, levels=levels(dat$pattern))

  # Group-level, by choice
  ppByChoiceByPostN <- aggregate(cbind(rt,accuracy,stim)~pattern*postn*choice, pp, mean)
  ppQuantilesByChoice <- aggregate(cbind(rt,accuracy,stim)~pattern*choice, ppByChoiceByPostN, quantile, c(.025, .5, .975))
  ppQuantilesByChoice$pattern <- factor(ppQuantilesByChoice$pattern, levels=levels(dat$pattern))

  # Subject-level, by choice
  ppBySubByChoiceByPostN <- aggregate(cbind(rt,accuracy,stim)~pattern*subjects*postn*choice, pp, mean)

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
  par(mar=c(1,2,2,2))
  shift <- .1
  ylim <- range(c(dataMeansByChoice$rt+dataSEs$rt, dataMeansByChoice$rt-dataSEs$rt, ppQuantilesByChoice$rt[,1], ppQuantilesByChoice$rt[,3])) #*c(.985, 1)
  if(ylim[1] > 0.7) ylim[1] <- 0.74  # dataset 3 needs some tweaking to prevent overlap of legend label and data
  if(ylim[2] < 0.54) ylim[1] <- 0.445  # dataset 2 also needs some tweaking to prevent overlap of legend label and data


  # data
  idx1 <- dataMeansByChoice$choice==1
  plot(x=as.numeric(dataMeansByChoice$pattern[idx1])+shift, y=dataMeansByChoice$rt[idx1], pch=4, lwd=2, ylim=ylim, xaxt='n', ylab=ylab1, xlab='', main=main1, xlim=range(as.numeric(dataMeansByChoice$pattern[idx1]))+c(-.2, .2))
  mtext(ylab1, side=2, line=2.5, cex=.66)
  axis(side=1, at=1:length(dataMeans$pattern), labels=levels(dataMeans$pattern), las=1)
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

  par(mar=c(4,2,2,2))
  shift = 0.1
  # ylim = range(c(dataMeans$choice-dataSEs$choice, dataMeans$choice+dataSEs$choice, dataMeans$stim-dataSEs$stim, dataMeans$stim+dataSEs$stim,
  #                ppQuantiles$choice))
  ylim <- c(0.1, .9)
  if(plotStimulusPattern) {
    plot(x=as.numeric(dataMeans$pattern), y=dataMeans$stim, lwd=2, ylim=ylim, xaxt='n', ylab=ylab2, xlab='Block', pch=4, main=main2)
    warning('Plotting *stim* as a function of previous stimuli')
  } else {
    plot(x=as.numeric(dataMeans$pattern), y=dataMeans$choice, lwd=2, ylim=ylim, xaxt='n', ylab=ylab2, xlab='Block', pch=4, main=main2)
  }
  mtext(ylab2, side=2, line=2.5, cex=.66)
  axis(side=1, at=1:length(dataMeans$pattern), labels=levels(dataMeans$pattern), las=1)
  xs <- as.numeric(ppQuantiles$pattern)
  arrows(x0=as.numeric(ppQuantiles$pattern),
         y0=ppQuantiles$choice[,1],
         y1=ppQuantiles$choice[,3], angle=90, code=3, length=.05, lwd=2)
  lines(xs[order(xs)], ppQuantiles$choice[,2][order(xs)])
  # abline(h=0.5, lty=2, col=1)
  legend('topright', c(paste0('p(choice=', as.character(dat[dat$choice==1,response_column_for_legend][2]), ')')), lwd=c(2), lty=c(1), col=par()$col, bty='n', cex=.75)
  # par(mar=c(4, 4, 4, 2) + 0.1)
}


dat[dat$proportion_word=='75','condition'] <- '75% Words'
dat[dat$proportion_word=='25','condition'] <- '25% Words'
pp_withSM[pp_withSM$proportion_word=='75','condition'] <- '75% Words'
pp_withSM[pp_withSM$proportion_word=='25','condition'] <- '25% Words'
pp_iidrdm[pp_iidrdm$proportion_word=='75','condition'] <- '75% Words'
pp_iidrdm[pp_iidrdm$proportion_word=='25','condition'] <- '25% Words'
pp_noB[pp_noB$proportion_word=='75','condition'] <- '75% Words'
pp_noB[pp_noB$proportion_word=='25','condition'] <- '25% Words'
dat$condition <- factor(dat$condition, levels=c('75% Words', '25% Words'))
pp_withSM$condition <- factor(pp_withSM$condition, levels=c('75% Words', '25% Words'))
pp_noB$condition <- factor(pp_noB$condition, levels=c('75% Words', '25% Words'))
pp_iidrdm$condition <- factor(pp_iidrdm$condition, levels=c('75% Words', '25% Words'))

summaryByPword1 <- get_accuracy_rt_by_condition(dat=dat, pp=pp_withSM)
summaryByPword2 <- get_accuracy_rt_by_condition(dat=dat, pp=pp_noB)
summaryByPword3 <- get_accuracy_rt_by_condition(dat=dat, pp=pp_iidrdm)


pdf(file='./figures/SM_effect_ds4.pdf', width=7, height=5)
par(mfcol=c(2,3), bty='l', oma=c(0,2,0,.5))
plotSequentialEffectsGroup(dat=dat, pp=pp_withSM, seqEffectsSummary = summaryByPword1,
                           main1='Dataset 4 (MS1-RDM)', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
abline(h=c(.25,.75), lty=2)
plotSequentialEffectsGroup(dat=dat, pp=pp_iidrdm, seqEffectsSummary = summaryByPword3,
                           main1='Dataset 4 (IID-RDM)', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
abline(h=c(.25,.75), lty=2)
plotSequentialEffectsGroup(dat=dat, pp=pp_noB, seqEffectsSummary = summaryByPword2,
                           nHistory = 3, main1='Dataset 4 (no adaptation)', full.legend=FALSE)
abline(h=c(.25,.75), lty=2)
dev.off()

# seqEffectsSummary1 <- getSequentialEffectsBySub(dat=dat, pp=pp_withSM, repeatAlternate = FALSE)
# seqEffectsSummary2 <- getSequentialEffectsBySub(dat=dat, pp=pp_noB, repeatAlternate = FALSE)
# seqEffectsSummary3 <- getSequentialEffectsBySub(dat=dat, pp=pp_noSM, repeatAlternate = FALSE)
#
# pdf(file='./figures/SM_effect_ds4.pdf', width=5, height=5)
# par(mfcol=c(2,2), bty='l', oma=c(0,2,0,0))
# plotSequentialEffectsGroup(dat=dat, pp=pp_withSM, seqEffectsSummary=seqEffectsSummary1,
#                            nHistory = 3, main1='Dataset 4 (with SM)', ylab1='RT (s)', ylab2='Choice proportion', full.legend=FALSE)
# plotSequentialEffectsGroup(dat=dat, pp=pp_noB, seqEffectsSummary=seqEffectsSummary2,
#                            nHistory = 3, main1='Dataset 4 (no adaptation)', full.legend=FALSE)
# dev.off()
# pp_withSM$R2 <- pp_withSM$R=='word'
# agg1acc <- aggregate(R2~subjects*proportionWord*S, pp_withSM, mean)
# aggregate(accuracy~proportionWord, agg1acc, mean)  ## just use means
# aggregate(accuracy~proportionWord, agg1acc, function(x) sd(x)/sqrt(length(x)))  ## just use SE
#
# pp_noB$R2 <- pp_noB$R=='word'
# agg2acc <- aggregate(accuracy~subjects*proportionWord, pp_noB, mean)
# aggregate(accuracy~proportionWord, agg2acc, mean)
# aggregate(accuracy~proportionWord, agg2acc, function(x) sd(x)/sqrt(length(x)))
#
#
#
# aggregate(rt~proportionWord*S, aggregate(rt~subjects*proportionWord*S, pp_withSM, mean), mean)
# aggregate(rt~proportionWord*S, aggregate(rt~subjects*proportionWord*S, pp_noB, mean), mean)

#
#
# mean(aggregate(accuracy~subjects, aggregate(accuracy~postn*subjects, pp_noSM, mean), mean)$accuracy)
# mean(aggregate(rt~subjects, aggregate(rt~postn*subjects, pp_noSM, mean), mean)$rt)



# ##
# dat <- EMC2:::loadRData('./datasets/wagenmakers2008exp2.RData')
# dat$choice <- dat$R == levels(dat$R)[1]
# pp_noSM$choice <- pp_noSM$R==dat[dat$choice==1,'R'][1]



# pp_ds4
## can we observe SH effects?
seqEffectsSummary4 <- getSequentialEffectsBySub(dat=dat, pp=pp_noSM, repeatAlternate = FALSE)
plotSequentialEffectsGroup(dat=dat, pp=pp_noSM, seqEffectsSummary=seqEffectsSummary4,
                           nHistory = 3, main1='Dataset 4', full.legend=FALSE)

## Ok, what about just a null model with no threshold bias?
# load('./datasets/wagenmakers2008exp2.RData')
# if('proportion_word' %in% colnames(dat)) dat$proportionWord <- dat$proportion_word  # new EMC doesnt like underscores in variable names
#
# # Default Flist (formula)
# Flist=list(v~lM,B~1,A~1,t0~1,s~lM)
# model <- RDM
# constants <- c(s=log(1), A=log(0))
# # Generate design ------------------------------------------------------------------
# ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  # average+difference
# Flist[[1]] = v~lM*W
# Ffactors=list(subjects=levels(dat$subjects),
#               S=c('nonword', 'word'),
#               W=levels(dat$W),
#               proportionWord=levels(dat$proportionWord))
# Ffunctions=list(errorAccumulator=function(d) as.character(d$lM)=="FALSE")
# Flist[[5]] = s~errorAccumulator
# constants=c(constants, v_Wvlf=0, v_Whf=0, v_Wnonword=0)
# Clist=NULL
#
# # Design
# design <- design(
#   factors=Ffactors,
#   Rlevels=levels(dat$R),
#   matchfun=function(d) d$S==d$lR,
#   contrasts=Clist,
#   formula=Flist,
#   functions=Ffunctions,
#   constants=constants,
#   model=model)
# emc <- make_emc(dat, design)
# emc <- fit(emc, cores_for_chains=3, cores_per_chain=19, verbose=TRUE, verboseProgress = TRUE, fileName = './samples-ds4-nobias.RData')
# ppNewNull <- predict(emc, n_cores=20)

ppNewNull$accuracy <- ppNewNull$R==ppNewNull$S
mean(aggregate(accuracy~subjects, aggregate(accuracy~postn*subjects, ppNewNull, mean), mean)$accuracy)


#mean(aggregate(accuracy~subjects, aggregate(accuracy~postn*subjects, all_pps[[4]]$FMonly, mean), mean)$accuracy)

## and if we simply turn off the threshold bias in the NULL model?
ds4_null <- loadDataSamplersPP('wagenmakers2008exp2', 'NULL', pp_conditional = TRUE)
ds4_null <- ds4_null$samplers
ds4_null[[1]]$samples$alpha['B',,] <- c(ds4_null[[1]]$samples$alpha['B',,] + ds4_null[[1]]$samples$alpha['B_lRword',,] + ds4_null[[1]]$samples$alpha['B_proportionWord75',,] + ds4_null[[1]]$samples$alpha['B_proportionWord75:lRword',,])/4
ds4_null[[2]]$samples$alpha['B',,] <- c(ds4_null[[2]]$samples$alpha['B',,] + ds4_null[[2]]$samples$alpha['B_lRword',,] + ds4_null[[2]]$samples$alpha['B_proportionWord75',,] + ds4_null[[2]]$samples$alpha['B_proportionWord75:lRword',,])/4
ds4_null[[3]]$samples$alpha['B',,] <- c(ds4_null[[3]]$samples$alpha['B',,] + ds4_null[[3]]$samples$alpha['B_lRword',,] + ds4_null[[3]]$samples$alpha['B_proportionWord75',,] + ds4_null[[3]]$samples$alpha['B_proportionWord75:lRword',,])/4

pp_nullNoB <- predict(ds4_null, n_post = 100, n_cores=20,
                      conditionalOnData=TRUE)
pp_nullNoB$accuracy <- pp_nullNoB$R==pp_nullNoB$S
mean(aggregate(accuracy~subjects, aggregate(accuracy~postn*subjects, pp_nullNoB, mean), mean)$accuracy)




# old stuff below ---------------------------------------------------------





# ds1.SMonly <- merge_dat_null_ms3('wagenmakers2004_CS', pp_ms3=pp_SMonly)
# ds1.AMonly <- merge_dat_null_ms3('wagenmakers2004_CS', pp_ms3=pp_AMonly)
# ds1.FMonly <- merge_dat_null_ms3('wagenmakers2004_CS', pp_ms3=pp_FMonly)
# mean(sapply(unique(ds1.SMonly$subjects), function(x) get_var_explained_ms3(ds1.SMonly[ds1.SMonly$subjects==x,])))
# mean(sapply(unique(ds1.AMonly$subjects), function(x) get_var_explained_ms3(ds1.AMonly[ds1.AMonly$subjects==x,])))
# mean(sapply(unique(ds1.FMonly$subjects), function(x) get_var_explained_ms3(ds1.FMonly[ds1.FMonly$subjects==x,])))



# plot_cdf(ds1$dat, post_predict=pp_SMonly, functions=list(acc=function(x) x$S==x$R), defective_factor='acc')
# plot_cdf(ds1$dat, post_predict=pp_AMonly, functions=list(acc=function(x) x$S==x$R), defective_factor='acc')
# plot_cdf(ds1$dat, post_predict=pp_FMonly, functions=list(acc=function(x) x$S==x$R), defective_factor='acc')

#aggregate(.~subjects, aggregate(cbind(dv,u,B_)~postn*subjects,npars,sd), mean)

#mean(npars$B_/npars$B)



#
#
# ## dataset 2: SA  get_var_explained <- function(x) cor(x[,'rt_ms3'], x[,'rt'])^2
# #   mean(sapply(unique(merged$subjects), function(x) get_var_explained(merged[merged$subjects==x,])))
# # T manipulation
# ds2_null <- get_pp2(task='forstmann2008', learningModel='NULL')
# dat2_null <- EMC2:::add_trials(ds2_null$dat)
# pp2_null <- ds2_null$pp2
#
# ds2_ms3 <- get_pp2(task='forstmann2008', learningModel='zSMuAHbV')
# pp2_ms3 <- ds2_ms3$pp2
#
# ## Merge mean predictions with data
# agg_iid <- aggregate(cbind(accuracy,rt)~trials*subjects, pp2_null, mean)
# agg_ms3 <-aggregate(cbind(accuracy,rt)~trials*subjects, pp2_ms3, mean)
#
# agg_combined <- merge(agg_iid, agg_ms3,by=c('subjects','trials'), suffixes=c('_iid', '_ms3'))
# merged <- merge(dat2_null, agg_combined, by=c('subjects','trials'))
#
# get_var_explained_ms3 <- function(x, colname='rt') cor(x[,paste0(colname, '_ms3')], x[,colname])^2
# get_var_explained_iid <- function(x, colname='rt') cor(x[,paste0(colname, '_iid')], x[,colname])^2
# mean(sapply(unique(merged$subjects), function(x) get_var_explained_ms3(merged[merged$subjects==x,])))
# mean(sapply(unique(merged$subjects), function(x) get_var_explained_iid(merged[merged$subjects==x,])))
#
# mean(sapply(unique(dat2_null$subjects), function(x) summary(lm(rt~E+R, dat2_null[dat2_null$subjects==x,]))$r.squared))
#
#
# ## Dataset 3
# ds3_null <- get_pp2(task='mileticvanmaanen2019exp2block2', learningModel='NULL')
# dat3_null <- EMC2:::add_trials(ds3_null$dat)
# pp3_null <- ds3_null$pp2
#
# ds3_ms3 <- get_pp2(task='mileticvanmaanen2019exp2block2', learningModel='zSMuAHbV')
# pp3_ms3 <- ds3_ms3$pp2
#
#
# ## dataset 4: base rate, difficulty manipulation
# ds4_null <- get_pp2(task='wagenmakers2008exp2', learningModel='NULL')
# dat4_null <- EMC2:::add_trials(ds4_null$dat)
# pp4_null <- ds4_null$pp2
#
# ds4_ms3 <- get_pp2(task='wagenmakers2008exp2', learningModel='zSMuAHbV')
# pp4_ms3 <- ds4_ms3$pp2
#
# ## Merge mean predictions with data
# agg_iid <- aggregate(cbind(accuracy,rt)~trials*subjects, pp4_null, mean)
# agg_ms3 <-aggregate(cbind(accuracy,rt)~trials*subjects, pp4_ms3, mean)
#
# agg_combined <- merge(agg_iid, agg_ms3,by=c('subjects','trials'), suffixes=c('_iid', '_ms3'))
# merged <- merge(dat4_null, agg_combined, by=c('subjects','trials'))
#
# get_var_explained_ms3 <- function(x, colname='rt') cor(x[,paste0(colname, '_ms3')], x[,colname])^2
# get_var_explained_iid <- function(x, colname='rt') cor(x[,paste0(colname, '_iid')], x[,colname])^2
# mean(sapply(unique(merged$subjects), function(x) get_var_explained_ms3(merged[merged$subjects==x,])))
# mean(sapply(unique(merged$subjects), function(x) get_var_explained_iid(merged[merged$subjects==x,])))
#
# #summary(lm(rt~W, dat4_null[dat4_null$subjects==1,]))$R_squared
#
# mean(sapply(unique(dat4_null$subjects), function(x) summary(lm(rt~W, dat4_null[dat4_null$subjects==x,]))$r.squared))
# mean(sapply(unique(dat4_null$subjects), function(x) summary(lm(rt~proportion_word, dat4_null[dat4_null$subjects==x,]))$r.squared))
# mean(sapply(unique(dat4_null$subjects), function(x) summary(lm(rt~proportion_word+W, dat4_null[dat4_null$subjects==x,]))$r.squared))
#


#
# ## Dataset 5, difficulty manipulation
# ds5_null <- get_pp2(task='mileticvanmaanen2019exp2block1', learningModel='NULL')
# dat5_null <- EMC2:::add_trials(ds5_null$dat)
# pp5_null <- ds5_null$pp2
#
# ds5_ms3 <- get_pp2(task='mileticvanmaanen2019exp2block1', learningModel='zSMuAHbV')
# dat5_null <- EMC2:::add_trials(ds5_ms3$dat)
# pp5_ms3 <- ds5_ms3$pp2
#
# ## Merge mean predictions with data
# #agg_iid <- aggregate(cbind(accuracy,rt)~trials*subjects, pp4_null, mean)
# agg_ms3 <-aggregate(cbind(accuracy,rt)~trials*subjects, pp4_ms3, mean)
# colnames(agg_ms3)<-c('trials','subjects','accuracy_ms3', 'rt_ms3')
# #agg_combined <- merge(agg_iid, agg_ms3,by=c('subjects','trials'), suffixes=c('_iid', '_ms3'))
# merged <- merge(dat4_null, agg_ms3, by=c('subjects','trials'))
#
# get_var_explained_ms3 <- function(x, colname='rt') cor(x[,paste0(colname, '_ms3')], x[,colname])^2
# #get_var_explained_iid <- function(x, colname='rt') cor(x[,paste0(colname, '_iid')], x[,colname])^2
# mean(sapply(unique(merged$subjects), function(x) get_var_explained_ms3(merged[merged$subjects==x,])))
# #mean(sapply(unique(merged$subjects), function(x) get_var_explained_iid(merged[merged$subjects==x,])))
#
# get_var_explained2 <- function(x, colname='rt') sum((x[,paste0(colname, '_ms3')]-x[,colname])^2)/sum((x[,colname]-mean(x[,colname]))^2)
# 1-mean(sapply(unique(merged$subjects), function(x) get_var_explained2(merged[merged$subjects==x,])))
#
# #get_var_explained3 <- function(x, colname='rt') sum((x[,paste0(colname, '_iid')]-x[,colname])^2)/sum((x[,colname]-mean(x[,colname]))^2)
# #1-mean(sapply(unique(merged$subjects), function(x) get_var_explained3(merged[merged$subjects==x,])))
#
#
# mean(sapply(unique(merged$subjects), function(x) summary(lm(rt~rt_ms3, merged[merged$subjects==x,]))$r.squared))
# #mean(sapply(unique(merged$subjects), function(x) summary(lm(rt~1, merged[merged$subjects==x,]))$r.squared))
#
# mean(sapply(unique(dat5_null$subjects), function(x) summary(lm(rt~difficulty, dat5_null[dat5_null$subjects==x,]))$r.squared))
#
#
#
# ## ToDo: Partition RT variance explained by mechanism.
#
#
#
# ## ToDo: Express parameter variability as SDs / coefficient of variation
# ## all_infofrom makeTablesParameters.R
#
# aggregate(cbind(meanb,u,`V3`)~subjects, aggregate(cbind(meanb,u,abs(db))~subjects*trials, all_info[[1]]$pp21,mean), function(x) sd(x)/mean(x))
#
# aggregate(cbind(meanb,u,db)~subjects, aggregate(cbind(meanb,u,db)~subjects*trials, all_info[[2]]$pp21,mean), function(x) sd(x)/mean(x))
# aggregate(cbind(meanb,u,db)~subjects, aggregate(cbind(meanb,u,db)~subjects*trials, all_info[[2]]$pp01,mean), function(x) sd(x)/mean(x))
#
# aggregate(cbind(meanb,u,`V3`)~subjects, aggregate(cbind(meanb,u,abs(db))~subjects*trials, all_info[[3]]$pp21,mean), function(x) sd(x)/mean(x))
# #aggregate(cbind(meanb,u,db)~subjects, aggregate(cbind(meanb,u,db)~subjects*trials, all_info[[3]]$pp21,mean), function(x) sd(x)/mean(x))
#
# aggregate(cbind(meanb,u,db)~subjects, aggregate(cbind(meanb,u,db)~subjects*trials, all_info[[3]]$pp21,mean), function(x) sd(x)/mean(x))
#
#
#
# ## Average SD on a trialwise level? Should be lower for MS3?
# for(pps in list(list(pp1_null,pp1_ms3),
#                 list(pp2_null,pp2_ms3),
#                 list(pp3_null,pp3_ms3),
#                 list(pp4_null,pp4_ms3))) {
#   pp_null <- pps[[1]]
#   pp_ms3 <- pps[[2]]
#   agg_iid <- aggregate(cbind(accuracy,rt)~trials*subjects, pp_null, sd)
#   agg_ms3 <- aggregate(cbind(accuracy,rt)~trials*subjects, pp_ms3, sd)
#
#   agg_combined <- merge(agg_iid, agg_ms3,by=c('subjects','trials'), suffixes=c('_iid', '_ms3'))
#   agg_combined <- aggregate((rt_iid/rt_ms3)~subjects,agg_combined,mean)
#   print(mean(agg_combined[,2]))
# }
