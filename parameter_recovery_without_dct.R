# parameter recovery
library(EMC2)
library(emcAdapt)
library(abind)
library(latex2exp)

source('./makeSamplers.R')

tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
learningModels <- c('zSMuAHbV', 'zSMuAHbV', 'zSMuAHbV','zSMuAHbV')
dctModels <- c('NULL', 'NULL', 'NULL', 'NULL')

for(i in c(3)) {
  task <- tasks[i]
  learningModel <- learningModels[i]
  dctModel <- dctModels[i]
  print(task)

  simData_fn <- paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData')
  print(load(paste0('./samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel,'.RData')))
  if(!file.exists(simData_fn)) {
    simulated_data <- post_predict(samples=samplers, n_post=1, use_par = 'mean')
    simulated_data <- simulated_data[order(simulated_data$subjects, simulated_data$trials),]

    attr_ <- attr(simulated_data, 'pars')

    # load corresponding 'real' dataset
    tmp <- EMC2:::loadRData(paste0('./datasets/', task, '.RData'))
    simulated_data <- simulated_data[,c(colnames(tmp), 'stim', 'resp', 'trials')]
    attr(simulated_data, 'pars') <- attr_
    save(simulated_data, file=paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData'))
  } else {
    load(simData_fn)
  }

  makeSamplers(decisionModel='RDM', learningModel=learningModel, dctModel=dctModel, dat=simulated_data,
               task=task, samples_dir='./parameter_recovery/samples', n_chains=3)
  samplers_fn = paste0('./parameter_recovery/samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel, '.RData')
  load(samplers_fn)
  progress <- EMC2:::check_progress(samplers, stage='sample', iter=1e3, max_gd=1.1, mean_gd=1.05, min_es=0, min_unique=600, max_trys=20, step_size=100, n_cores=8, verbose=FALSE)
  if(!progress$done) {
    out = run_emc(samplers_fn, cores_per_chain = 8, cores_for_chains=3, iter=1000, max_trys=10, verbose = TRUE, verboseProgress = TRUE)
  }
}



# Visualize recovery ------------------------------------------------------
#task <- 'wagenmakers2004_CS' #'wagenmakers2004_CS'
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
for(task in tasks) {
  simulatedData <- EMC2:::loadRData(paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData'))
  samples <- EMC2:::loadRData(Sys.glob(paste0('./parameter_recovery/samples/', task, '_*.RData'))[1])
  filter <- names(which(chain_n(samples)[1,]>0)); filter <- filter[length(filter)]
  #plot_chains(samples, selection='mu', filter=filter)
  #plot_chains(samples, selection='alpha', filter=filter, subject=names(samples[[1]]$data)[1])


  #subfilter=0
  idx <- which(samples[[1]]$samples$stage==filter)
  #idx <- idx[(length(idx)-subfilter):length(idx)]
  alpha <- abind(lapply(samples, function(x) x$samples$alpha[,,idx]), along=3)
  mu <- abind(lapply(samples, function(x) x$samples$theta_mu[,idx]), along=2)
  fitParsMean <- apply(alpha,1:2,mean)
  fitParsQps <- apply(alpha,1:2,quantile, c(0.025, .5, 0.975))
  fitMuQps <- apply(mu, 1, quantile, c(0.025, 0.5, 0.975))

  truePars <- attr(simulatedData, 'pars')


  pname_mapping = list('v'='$u$',
                       'v_lMd'='$dv$',
                       'B'='b',
                       'B_Ea-n'='$b_{a-n}$',
                       'B_Ea-s'='$b_{a-s}$',
                       'B_lRd'='$b_{l-r}$',
                       't0'='$t_0$',
                       's_lMd'='$ds$',
                       'alpha3'='$\\alpha_{FM}$',
                       'weight3'='$w_{FM}$',
                       'alpha2'='$\\alpha_{AM}$',
                       'weight2'='$w_{AM}$',
                       'alpha1'='$\\alpha_{SM}$',
                       'weight1'='$w_{SM}$',
                       's_errorAccumulatorTRUE'='$ds$',
                       'v_lMTRUE'='$dv$',
                       'v_lMTRUE:Wvlf'='$dv:W_{vlf}$',
                       'v_lMTRUE:Whf'='$dv:W_{hf}$',
                       'v_lMTRUE:Wnonword'='$dv:W_{hf}$')

  for(i in 1:10) {
    pname_mapping[paste0('vcos', i)] = paste0('$\\beta_{',i,'}$')
    pname_mapping[paste0('Bcos', i)] = paste0('$\\beta_{',i,'}$')
  }


  if(task == 'forstmann2008') {
    nrows=3
  } else if(task == 'wagenmakers2004_CS') {
    nrows=3
  } else if(task == 'mileticvanmaanen2019exp2block2') {
    nrows=3
  } else {
    nrows=3
  }

  for(ftype in c('jpeg', 'pdf')) {
    fn = paste0('./figures/parameter_recovery_', task, '.pdf')
    if(ftype == 'pdf') pdf(file=fn, width=8, height=nrows*2)
    if(ftype == 'jpeg') jpeg(gsub('.pdf','.jpeg',fn), width=8, height=nrows*2, units='in', quality=100, res=500)

    par(mfrow=c(nrows,5), mar=c(3,3,2,.5), bty='l', mgp=c(2,.5,0), tck=-0.05)
    for(pname in colnames(truePars)) {
      pname_label = pname_mapping[[pname]]
      x<-truePars[,pname]
      y<-fitParsQps[2,pname,]
      xlim <- ylim <- range(c(x,fitParsQps[,pname,]))

      CI_contains_x <- fitParsQps[3,pname,] > x & fitParsQps[1,pname,] < x

      xmean <- mean(x)

      xlab <- ifelse(which(pname == colnames(truePars))>(nrows-1)*5, 'Data generating', '')
      ylab <- ifelse(which(pname == colnames(truePars)) %% 5 == 1, 'Posterior (95% CI)', '')
      plot(x,y, main=TeX(pname_label), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
      arrows(x0=x, y0=fitParsQps[1,pname,], y1=fitParsQps[3,pname,], code=3, angle=90, length=.025, col=adjustcolor(1, alpha.f=.5))
      legend('topleft', paste0('r=',round(cor(x,y),2)), bty='n', cex=.8)
      legend('bottomright', paste0(round(mean(CI_contains_x)*100,1), '%'), bty='n', cex=.8)
      abline(a=0,b=1)

      # group-level recovery plot
      # if(pname %in% colnames(fitMuQps)) {
      #   usr <- par('usr')
      #   arrows(x0=usr[1],
      #          y0=fitMuQps[1,pname],
      #          y1=fitMuQps[3,pname],
      #          code=3,
      #          angle=90,
      #          length=0.025, col=2, lwd=2,
      #          xpd=TRUE
      #   )
      #   segments(x0=xmean, x1=xmean, y0=usr[3], y1=xmean, col=2, lwd=2)
      #   segments(x0=xmean, x1=usr[1], y0=xmean, y1=xmean, col=2, lwd=2)
      # }
    }
    dev.off()
  }
}




#
# print(load('./parameter_recovery/datasets/wagenmakers2008exp2_simulated_data.RData'))
# print(load('./parameter_recovery/samples/wagenmakers2008exp2_model-RDM-zSMuAHbV-DCT-NULL.RData'))
# pp <- post_predict(samplers, n_cores=20)
#
#
# par(mfrow=c(3,3))
# plot_fit(data=simulated_data, pp=pp, factors=c('S','W'))







# simData_fn <- paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData')
# print(load(paste0('./samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel,'.RData')))
# if(!file.exists(simData_fn)) {
#   simulated_data <- post_predict(samples=samplers, n_post=1, use_par = 'mean')
#   simulated_data <- simulated_data[order(simulated_data$subjects, simulated_data$trials),]
#   save(simulated_data, file=paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData'))
# } else {
#   load(simData_fn)
# }
#
# makeSamplers(decisionModel='RDM', learningModel=learningModel, dctModel=dctModel, dat=simulated_data,
#              task=task, samples_dir='./parameter_recovery/samples', n_chains=3)
# samplers_fn = paste0('./parameter_recovery/samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel, '.RData')
# out = run_emc(samplers_fn, cores_per_chain = 6, cores_for_chains=3, verbose = TRUE, verboseProgress = TRUE)
#
#
# # Dataset 2 ---------------------------------------------------------------
# task <- 'forstmann2008'
# learningModel <- 'zSMuAHbV'
# dctModel <- 'NULL'
#
# print(load(paste0('./samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel,'.RData')))
# if(!file.exists(simData_fn)) {
#   simulated_data <- post_predict(samples=samplers, n_post=1, use_par = 'mean')
#   simulated_data <- simulated_data[order(simulated_data$subjects, simulated_data$trials),]
#   save(simulated_data, file=paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData'))
# } else {
#   load(simData_fn)
# }
#
# makeSamplers(decisionModel='RDM', learningModel=learningModel, dctModel=dctModel, dat=simulated_data,
#              task=task, samples_dir='./parameter_recovery/samples', n_chains=3)
# samplers_fn = paste0('./parameter_recovery/samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel, '.RData')
# out = run_emc(samplers_fn, cores_per_chain = 6, cores_for_chains=3, verbose = TRUE, verboseProgress = TRUE)
#
#
# # Dataset 3 ---------------------------------------------------------------
# task <- 'mileticvanmaanen2019exp2block2'
# learningModel <- 'zSMaAHbV'
# dctModel <- 'NULL'
#
# print(load(paste0('./samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel,'.RData')))
# if(!file.exists(simData_fn)) {
#   simulated_data <- post_predict(samples=samplers, n_post=1, use_par = 'mean')
#   simulated_data <- simulated_data[order(simulated_data$subjects, simulated_data$trials),]
#   save(simulated_data, file=paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData'))
# } else {
#   load(simData_fn)
# }
#
# makeSamplers(decisionModel='RDM', learningModel=learningModel, dctModel=dctModel, dat=simulated_data,
#              task=task, samples_dir='./parameter_recovery/samples', n_chains=3)
# samplers_fn = paste0('./parameter_recovery/samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel, '.RData')
# out = run_emc(samplers_fn, cores_per_chain = 6, cores_for_chains=3, verbose = TRUE, verboseProgress = TRUE)
#
#
# # Dataset 4 ---------------------------------------------------------------
# task <- 'wagenmakers2008exp2'
# learningModel <- 'zSMuAHbV'
# dctModel <- 'NULL'
#
# print(load(paste0('./samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel,'.RData')))
# simulated_data <- post_predict(samples=samplers, n_post=1, use_par = 'mean')
# simulated_data <- simulated_data[order(simulated_data$subjects, simulated_data$trials),]
# save(simulated_data, file=paste0('./parameter_recovery/datasets/', task, '_simulated_data.RData'))
#
# makeSamplers(decisionModel='RDM', learningModel=learningModel, dctModel=dctModel, dat=simulated_data,
#              task=task, samples_dir='./parameter_recovery/samples', n_chains=3)
# samplers_fn = paste0('./parameter_recovery/samples/', task, '_model-RDM-', learningModel, '-DCT-', dctModel, '.RData')
# out = run_emc(samplers_fn, cores_per_chain = 6, cores_for_chains=3, verbose = TRUE, verboseProgress = TRUE)


