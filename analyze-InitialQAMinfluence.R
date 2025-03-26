rm(list=ls())
library(EMC2)
library(emcAdapt)
library(parallel)
# install.packages('latex2exp')
library(latex2exp)
source('./extra_EMC2_functions/adaptive.R')
source('./extra_EMC2_functions/make_data.R')
source('./extra_EMC2_functions/model_RDMdynamic.R')
source('./extra_EMC2_functions/utils.R')
figures_dir = './figures'


## Load PPs of model with Q_{AM,0} = 1
dss <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
xlims <- list(c(1,20), c(1,20), c(1,20), c(1,20))
ylabs <- c(latex2exp::TeX('Difference in $Q_{AM}$-values'), '', '', '')
pdf(file='./figures/initialQAMvalue.pdf', width=7, height=3.5)
par(mfrow=c(1,4), bty='l', mar=c(4,4,3,0))
for(dsnum in 1:length(dss)) {
  ds <- dss[dsnum]
  ds1_qAM_1 <- loadDataSamplersPP(ds, learningModel='zSMuAHbV', samples_dir='./samples', pp_conditional = TRUE)
  ds1_qAM_p5 <- loadDataSamplersPP(ds, learningModel='zSMuAHbV', samples_dir = './samples_initialQAM', pp_conditional = TRUE)

  npars1 <- attr(ds1_qAM_1$pp, 'npars')
  nparsp5 <- data.frame(attr(ds1_qAM_p5$pp, 'npars'))
  meanq02_1 <- aggregate(q02~subjects*trials, npars1, mean)
  meanq02_p5 <- aggregate(q02~subjects*trials, nparsp5, mean)

  # difference
  meanq02 <- cbind(meanq02_1, q02p5=meanq02_p5[,'q02'])
  meanq02$difference <- meanq02$q02-meanq02$q02p5
  plot(0,0,xlim=xlims[[dsnum]], ylim=c(-.2,.5), xlab='Trial', ylab=ylabs[dsnum],type='n', main=paste0('Dataset ', dsnum))
  abline(h=0, lty=2, col='grey')
  for(subject in unique(meanq02_1$subjects)) {
    lines(meanq02$trials[meanq02$subjects==subject],meanq02$difference[meanq02$subjects==subject])
  }
}
dev.off()
