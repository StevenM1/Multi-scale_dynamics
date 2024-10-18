## preps some datasets

# Wagenmakers et al 2004 'CS' experiment -------------------------------------------------------
for(condition in c('CS', 'CL', 'SS', 'SL')) {
  dat <- NULL
  for(subject in 1:6) {
    thisDat <- read.csv(paste0('./datasets/noisedat/S0', subject, condition, '.DAT'), header=FALSE, sep=' ', col.names=c('response_key', 'RT', 'practice_or_real', 'stimulus_is_odd', 'RSI', 'stimulus_number'))
    thisDat$subject <- subject
    dat <- rbind(dat, thisDat)
  }

  dat$subjects <- as.factor(dat$subject)
  dat$rt <- dat$RT/1000                 # seconds
  dat$S <- factor(ifelse(dat$stimulus_is_odd==0, 'even', 'odd'))
  dat$R <- factor(ifelse(dat$response_key=='/', 'even', 'odd'))
  dat$accuracy <- as.numeric(dat$S) == as.numeric(dat$R)     # unfortunately quite high
  dat <- dat[,c('subjects', 'S', 'R', 'rt')]
  dat <- dat[dat$rt>.15,]

  dat$accuracy <- as.numeric(as.numeric(dat$S) == as.numeric(dat$R))
  dat$choice <- as.character(dat$R) == 'odd'      # 1 == odd == left, 0 == even == right.
  save(dat, file=paste0('./datasets/wagenmakers2004_', condition, '_censored.RData'))
}
rm(list=ls())



## Miletic van Maanen exp 2 block 1 and long deadline
print(load('./datasets/timingDM2.RData'))
dat <- choiceDat[choiceDat$block=='Palmer',]
dat <- dat[dat$include,]           # ONLY included
dat <- dat[!is.na(dat$RT),]
dat <- dat[(dat$RT<3) & (dat$RT>.15),]  # not too fast or too slow
dat$S <- NA
dat$S[as.character(dat$response)=='left' & dat$acc==1] <- 'left'
dat$S[as.character(dat$response)=='right' & dat$acc==1] <- 'right'
dat$S[as.character(dat$response)=='left' & dat$acc==0] <- 'right'
dat$S[as.character(dat$response)=='right' & dat$acc==0] <- 'left'
dat$S <- factor(dat$S, levels=c('left', 'right'))
#dat$S <- dat$stimulus
dat$R <- dat$response
dat$rt <- dat$RT
dat$accuracy <- dat$acc
dat$subjects <- dat$pp
dat <- dat[!is.na(dat$rt),]
dat$choice <- as.character(dat$R) == 'left'
dat <- droplevels(dat)
dat <- dat[,c('subjects', 'S', 'R', 'rt', 'stimStrength', 'choice', 'accuracy')]
save(dat, file='./datasets/mileticvanmaanen2019exp2block1.RData')

## Miletic van Maanen exp 2 block 1 and long deadline
print(load('./datasets/timingDM2.RData'))
dat <- choiceDat[choiceDat$block=='lowUrg',]
dat <- dat[dat$include,]           # ONLY included
dat <- dat[(dat$RT<3) & (dat$RT>.15),]  # not too fast or too slow
dat$S <- NA
dat$S[as.character(dat$response)=='left' & dat$acc==1] <- 'left'
dat$S[as.character(dat$response)=='right' & dat$acc==1] <- 'right'
dat$S[as.character(dat$response)=='left' & dat$acc==0] <- 'right'
dat$S[as.character(dat$response)=='right' & dat$acc==0] <- 'left'
dat$S <- factor(dat$S, levels=c('left', 'right'))
dat$R <- dat$response
dat$rt <- dat$RT
dat$accuracy <- dat$acc
dat$subjects <- dat$pp
dat <- dat[!is.na(dat$rt),]
dat$choice <- as.character(dat$R) == 'left'
dat <- droplevels(dat)
dat <- dat[,c('subjects', 'S', 'R', 'rt', 'stimStrength', 'choice', 'accuracy')]
save(dat, file='./datasets/mileticvanmaanen2019exp2block2.RData')


print(load('./datasets/timingDM2.RData'))   # Not used
dat <- choiceDat[choiceDat$block=='highUrg',]
dat <- dat[dat$include,]           # ONLY included
dat <- dat[(dat$RT<3) & (dat$RT>.15),]  # not too fast or too slow
dat$S <- NA
dat$S[as.character(dat$response)=='left' & dat$acc==1] <- 'left'
dat$S[as.character(dat$response)=='right' & dat$acc==1] <- 'right'
dat$S[as.character(dat$response)=='left' & dat$acc==0] <- 'right'
dat$S[as.character(dat$response)=='right' & dat$acc==0] <- 'left'
dat$S <- factor(dat$S, levels=c('left', 'right'))
dat$R <- dat$response
dat$rt <- dat$RT
dat$accuracy <- dat$acc
dat$subjects <- dat$pp
dat <- dat[!is.na(dat$rt),]
dat$choice <- as.character(dat$R) == 'left'
dat <- droplevels(dat)
dat <- dat[,c('subjects', 'S', 'R', 'rt', 'stimStrength', 'choice', 'accuracy')]
save(dat, file='./datasets/mileticvanmaanen2019exp2block3.RData')



##### PNAS data
print(load('./datasets/PNAS.RData'))
data <- dat
head(data)
data$subjects <- data$s
data$rt <- data$RT
data$trials <- data$trial
data$accuracy <- data$C
levels(data$R) <- levels(data$S)
data$choice = data$R == 'right'
data$accuracy <- as.numeric(data$accuracy)
dat <- data[,c('subjects', 'S', 'R', 'rt', 'E', 'trials', 'accuracy', 'choice')]
save(dat, file='./datasets/forstmann2008.RData')



### Wagenmakers2008exp2
expt1 <- read.delim("./datasets/LexDecData/SpeedAccData.txt",header=F)
expt2 <- read.delim("./datasets/LexDecData/PropData.txt",header=F)
names(expt1) <- c("s","blk","practice","E","stim","F","R","RT","censor")
round(prop.table(table(expt1$censor)),3)
# 0     1
# 0.979 0.021
names(expt2) <- c("s","blk","practice","P","stim","F","R","RT","censor")
round(prop.table(table(expt2$censor)),3)
#     0     1
# 0.983 0.017
expt1 <- expt1[expt1$censor!=1,-c(3,9)]
expt2 <- expt2[expt2$censor!=1,-c(3,9)]
expt1$s=factor(expt1$s)
expt2$s=factor(expt2$s)
expt1$E=factor(expt1$E,labels=c("accuracy","speed"))
expt2$P=factor(expt2$P,labels=c("25","75"))
expt1$S <- expt1$F<4
expt1$S <- factor(expt1$S,labels=c("nonword","word"))
expt2$S <- expt2$F<4
expt2$S <- factor(expt2$S,labels=c("nonword","word"))
expt1$R <- factor(expt1$R,labels=c("nonword","word"))
expt2$R <- factor(expt2$R,labels=c("nonword","word"))
expt1$F <- ((expt1$F-1) %% 3)+1
expt2$F <- ((expt2$F-1) %% 3)+1
expt1$F <- factor(expt1$F,labels=c("hf","lf","vlf"))
expt2$F <- factor(expt2$F,labels=c("hf","lf","vlf"))
expt1$C <- toupper(expt1$S)==toupper(expt1$R)
expt2$C <- toupper(expt2$S)==toupper(expt2$R)
# add in combined NW and W factor (4 levels, lump NW together)
expt1$W <- as.character(expt1$F)
expt2$W <- as.character(expt2$F)
expt1$W[expt1$S=="nonword"] <- "nonword"
expt2$W[expt2$S=="nonword"] <- "nonword"
expt1$W <- factor(expt1$W,c("hf","lf","vlf","nonword"))
expt2$W <- factor(expt2$W,c("hf","lf","vlf","nonword"))
dat <- expt2
colnames(dat) <- c('subjects', 'block', 'proportion_word', 'stim_identifier', 'frequency', 'R', 'rt', 'S', 'accuracy', 'W')
head(dat)
dat$W <- factor(dat$W, levels=c('lf', 'vlf', 'hf', 'nonword'))
save(dat, file='datasets/wagenmakers2008exp2.RData')

# remove subs with accuracy >96%
toremove <- which(aggregate(accuracy~subjects,dat,mean)[,2]>.96)
dat <- droplevels(dat[!dat$subjects %in% toremove,])
save(dat, file='datasets/wagenmakers2008exp2_subset.RData')


#save(expt1,expt2,file="ejldt-dat.RData")  # Not used, just saved here
dat <- expt1
colnames(dat) <- c('subjects', 'block', 'SAT', 'stim_identifier', 'frequency', 'R', 'rt', 'S', 'accuracy', 'W')
head(dat)
dat$W <- factor(dat$W, levels=c('lf', 'vlf', 'hf', 'nonword'))
save(dat, file='datasets/wagenmakers2008exp1.RData')



