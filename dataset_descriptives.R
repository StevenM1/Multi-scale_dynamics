## Descriptives
rm(list=ls())

rsis <- c('550-950 ms (mean 750 ms)', '', '0.8 s', '0.15 s')
forstmann_ISI <- c('2850 ms')
# forstmann:
# cue = 1000 ms
# cue-stim interval = 500 ms
# stim = 1000 ms
# RSI = unclear, but if we assume response triggers feedback, it's 350 millisecond (ie the feedback)
# feedback = 350 ms
# ISI then is about 1000 ms + 500 ms + 1000 ms (NB: mean RT is shorter than stimulus duration, but stimulus duration is fixed to 1000 ms) + 350 ms = 2850 ms


for(task in c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2', 'mileticvanmaanen2019exp2block1')) {
  load(paste0('./datasets/', task, '.RData'))
  print(task)
  nsubs <- length(unique(dat$subjects))
  print(paste0('Number of subjects: ', nsubs))

  rtrange <- range(aggregate(rt~subjects,dat,length)$rt)
  print(paste0('Number of trials (min/max):', rtrange[1], ' -- ', rtrange[2]))

  mrts <- aggregate(rt~subjects,dat,mean)$rt
  print(paste0('Mean RT (SE): ' , round(mean(mrts),3), ' (', round(sd(mrts)/sqrt(nsubs),3), ')'))

  macc <- aggregate(accuracy~subjects,dat,mean)$accuracy
  print(paste0('Mean accuracy (SE): ' , round(mean(macc),3), ' (', round(sd(macc)/sqrt(nsubs),3), ')'))

  cat('\n')
}
