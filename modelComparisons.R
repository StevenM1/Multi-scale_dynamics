rm(list=ls())
library(EMC2)
library(emcAdapt)
library(openxlsx)
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMsNewEMC')
source(file.path(root_dir, 'extra_EMC2_functions/compare_IC.R'))
source(file.path(root_dir, 'extra_EMC2_functions/adaptive.R'))
source(file.path(root_dir, 'extra_EMC2_functions/model_RDMdynamic.R'))
source(file.path(root_dir, 'extra_EMC2_functions/utils.R'))

model_averaging <- function(IC_for, IC_against) {
  if(is.null(IC_for)) return(NULL)

  if(is.data.frame(IC_for)){
    # Recursive call to make it work with the output of compare
    MD <- model_averaging(IC_for$MD, IC_against$MD)
    BPIC <- model_averaging(IC_for$BPIC, IC_against$BPIC)
    DIC <- model_averaging(IC_for$DIC, IC_against$DIC)
    return(rbind(MD = MD, BPIC = BPIC, DIC = DIC))
  }
  # Combine the IC values from both groups
  all_IC <- c(IC_for, IC_against)

  # Find the smallest IC value (for numerical stability)
  min_IC <- min(all_IC)

  # Compute the unnormalized weights
  unnorm_weights <- exp(-0.5 * (all_IC - min_IC))

  # Normalize weights so they sum to 1
  weights <- unnorm_weights / sum(unnorm_weights)

  # Separate the weights for the two groups
  weight_for <- sum(weights[seq_along(IC_for)])
  weight_against <- sum(weights[(length(IC_for) + 1):length(all_IC)])

  # Compute the Bayes factor: evidence in favor relative to against
  bayes_factor <- weight_for / weight_against

  # Return the results as a data frame
  return(data.frame(
    wFor = weight_for,
    wAgainst = weight_against,
    Factor = bayes_factor
  ))
}

get_model_comparison <- function(fns, subfilter=0, byParticipant=FALSE, hasTrend=FALSE) {
  allSamplers <- list()
  allGds <- c()
  for(fn in fns) {
    samplers <- EMC2:::loadRData(fn)
#    modelName <- learningModel <- gsub('.*?((v|z)SM(v|u|b)AH(b|u)V).*', '\\1', fn, perl=TRUE)
    modelName <- learningModel <- gsub('.*_model-?(.*)_trend-.*', '\\1', fn, perl=TRUE)
    if(hasTrend) modelName <- learningModel <- gsub('.*_model-?(.*)_trend-(.*).RData', '\\2', fn, perl=TRUE)
    nchains <- chain_n(samplers)
    if(nchains[1,4]>=100) {
      allSamplers[[modelName]] <- samplers
      gdmu <- EMC2:::gd_summary(samplers, selection='alpha')
      gdalpha <- EMC2:::gd_summary(samplers, selection='mu')
      allGds <- c(allGds, max(c(gdmu, gdalpha)))
    }
  }
  if(length(allSamplers)>0) {
    if(byParticipant) {
      out <- EMC2:::compare_subject(allSamplers, print_summary = FALSE)
      wDIC <- lapply(out,function(x)x["wDIC"])
      wBPIC <- lapply(out,function(x)x["wBPIC"])
      pDIC <- do.call(rbind,lapply(wDIC,function(x){
        setNames(data.frame(t(x)),paste("wDIC",rownames(x),sep="_"))}))
      pBPIC <- do.call(rbind,lapply(wBPIC,function(x){
        setNames(data.frame(t(x)),paste("wBPIC",rownames(x),sep="_"))}))
      # print(round(cbind(pDIC,pBPIC),3))
      mnams <- unlist(lapply(strsplit(dimnames(pDIC)[[2]],"_"),function(x){x[[2]]}))
      # cat("\nWinners\n")
      mc <- rbind(DIC=table(factor(mnams[apply(pDIC,1,which.max)], levels=mnams)),
                  BPIC=table(factor(mnams[apply(pBPIC,1,which.max)], levels=mnams)))
      print(mc)
    } else {
      mc = EMC2:::compare(allSamplers, print_summary=FALSE, subfilter=subfilter, BayesFactor = FALSE)
      mc <- cbind(mc, max.gd=allGds)
      mc <- get_delta_IC(mc)
      print(round(mc,3))
      mc
    }
  }
}

get_delta_IC <- function(mc) {
  mc$dDIC <- mc$DIC-min(mc$DIC)
  mc$dBPIC <- mc$BPIC-min(mc$BPIC)
  return(mc)
}

getFns <- function(task, samples_dir=c('./samples/', './samples_mc123/', './samples_trends/')) {
  all_fns <- Sys.glob(paste0(samples_dir, task, '*.RData'))
  fns = all_fns[!grepl('_pp', all_fns)]
  return(fns)
}

make_MC_table <- function(fns_mc1, row_names, row_order=c(1,2,3), ...) {
  mc1 <- lapply(fns_mc1, get_model_comparison, byParticipant=FALSE, ...)
  mc1_raw <- mc1
  names(mc1) <- paste0('ds', 1:4)
  mc1 <- do.call(cbind, mc1)
  if(!is.null(row_names)) row.names(mc1) <- row_names
  if(!is.null(row_order)) mc1 <- mc1[row_order,]

  ## By subject..
  mc1P <- lapply(fns_mc1, get_model_comparison, byParticipant=TRUE, ...)
  mc1Psummary <- do.call(cbind, lapply(mc1P, function(x) x[2,]))
  if(!is.null(row_order)) mc1Psummary <- mc1Psummary[row_order,]
  colnames(mc1Psummary) <- paste0('ds', 1:4, '.N')
  mc1combined <- cbind(round(mc1), mc1Psummary)
  mc1combined <- mc1combined[,c('ds1.dBPIC', 'ds1.N', 'ds2.dBPIC', 'ds2.N', 'ds3.dBPIC', 'ds3.N', 'ds4.dBPIC', 'ds4.N')]
  return(list(table=mc1combined, group=mc1, participant=mc1P, group_raw=mc1_raw))
}



datasets <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')
all_fns <- lapply(datasets, getFns)

# MC 1: SM vs null --------------------------------------------------------------------
include_fns_mc1 <- '(vSM|zSM|NULL)_trend-NULL'
fns_mc1 <- lapply(all_fns, function(x) x[grepl(include_fns_mc1,x)])
out1 <- make_MC_table(fns_mc1, row_names=c('NULL', 'dr~SM', 'db~SM'), row_order=1:3)
write.xlsx(out1$table, file='./tables/mc1.xlsx', colNames = TRUE, rowNames = TRUE)

# MC2: SM+AM vs SM -------
include_fns_mc2 <- '(zSM|zSMuAH|zSMbAH|zSMvAH)_trend-NULL'
fns_mc2 <- lapply(all_fns, function(x) x[grepl(include_fns_mc2,x)])
out2 <- make_MC_table(fns_mc2, row_names=c('db~SM, b~AM', 'db~SM', 'db~SM, u~AM', 'db~SM, v~AM'),
                      row_order=c(2,1,3,4))
write.xlsx(out2$table, file='./tables/mc2.xlsx', colNames = TRUE, rowNames = TRUE)

# MC3: SM+AM+FM vs SM+AM -------
include_fns_mc3 <- '(zSMuAH|zSMuAHbV|zSMuAHuV)_trend-NULL'
fns_mc3 <- lapply(all_fns, function(x) x[grepl(include_fns_mc3,x)])
out3 <- make_MC_table(fns_mc3, row_names=c('db~SM, u~AM, b~FM', 'db~SM, u~AM', 'db~SM, u~AM, u~FM'),
                      row_order=c(2,1,3))
write.xlsx(out3$table, file='./tables/mc3.xlsx', colNames = TRUE, rowNames = TRUE)


# MC 4: All permutations ----------
processByParticipant <- function(mcP) {
  ## function to sum all pps for which SM prefers Z or V; AM u or b or v; and FM b or u
  byParticipant <- data.frame(matrix(NA, nrow=2,ncol=7))
  colnames(byParticipant) <- c('zSM', 'vSM', 'uAH', 'bAH', 'vAH', 'uV', 'bV')
  row.names(byParticipant) <- c('DIC', 'BPIC')
  for(coln in colnames(byParticipant)) {
    byParticipant[,coln] <- apply(mcP[,grepl(coln, colnames(mcP))],1,sum)
  }
  return(byParticipant)
}
include_fns_mc4 <- '(v|z)SM(v|u|b)AH(u|b)V_trend-NULL'
fns_mc4 <- lapply(all_fns, function(x) x[grepl(include_fns_mc4,x)])
out4 <- make_MC_table(fns_mc4, row_names=NULL, row_order=NULL)
out4$table
byparticipant <- do.call(rbind, lapply(out4$participant, function(x) processByParticipant(x)[2,]))
write.xlsx(out4$table, file='./tables/mc4.xlsx', colNames = TRUE, rowNames = TRUE)
all_weights <- list() #$vector(mode='list', length=length(c('zSM', 'vSM', 'uAH', 'vAH', 'bAH', 'uV', 'bV')))
for(mechanism in c('zSM', 'vSM', 'uAH', 'vAH', 'bAH', 'uV', 'bV')) {
  weights <- setNames(sapply(1:4, function(x) model_averaging(IC_for = out4$group_raw[[x]][grepl(mechanism, row.names(out4$group_raw[[x]])),],
                                                              IC_against=out4$group_raw[[x]][!grepl(mechanism, row.names(out4$group_raw[[x]])),])[1,1]),
                      paste0('ds', 1:4))
  all_weights[[mechanism]]<-weights

  print(paste0('Weight for ', mechanism))
  print(round(weights,2))
}
summary_table <- data.frame(matrix(NA, ncol=8, nrow=9))
row.names(summary_table) <- c('zSM', 'vSM', 'AH:', 'uAH', 'vAH', 'bAH', 'FM:', 'uV', 'bV')
colnames(summary_table) <- c(paste0('ds.', 1:4, '_wBPIC'), paste0('ds.', 1:4, '_N'))
for(mechanism in c('zSM', 'vSM', 'uAH', 'vAH', 'bAH', 'uV', 'bV')) {
  summary_table[mechanism, 1:4] <- round(all_weights[[mechanism]],2)
  summary_table[mechanism, 5:8] <- byparticipant[[mechanism]]
}
summary_table <- summary_table[,c(1,5,2,6,3,7,4,8)]
row.names(summary_table) <- c('db~SM','dr~SM', 'AM', 'u~AM', 'v~AM', 'b~AM', 'FM','u~FM','b~FM')
write.xlsx(summary_table, file='./tables/mc4b.xlsx', colNames = TRUE, rowNames = TRUE)


# MC 5: Trends vs no trends -----------------------------------------------
include_fns_mc5 <- 'zSMuAHbV_trend-(DCT|NULL)'
fns_mc5 <- lapply(all_fns, function(x) x[grepl(include_fns_mc5,x)])
out5 <- make_MC_table(fns_mc5, row_names=c('MS3-RDM', 'MS3-RDM + B~DCT', 'MS3-RDM + u~DCT', 'MS3-RDM + v~DCT'),
                      row_order=c(1:4), hasTrend=TRUE)
write.xlsx(out5$table, file='./tables/mc5.xlsx', colNames = TRUE, rowNames = TRUE)





# ## MC2 --------------------------------------------------------------------
# fns_mc2 <- lapply(all_fns, function(x) x[grepl(include_fns_mc2,x)])
# mc2 <- lapply(fns_mc2, get_model_comparison, byParticipant=FALSE)
# names(mc2) <- paste0('ds', 1:4)
# mc2 <- do.call(cbind, mc2)
# row.names(mc2) <- c('z~SM, b~AH', 'z~SM', 'z~SM, u~AH', 'z~SM, v~AH')
# mc2 <- mc2[c(2,1,3,4),]
# mc2[,grepl('delta_BPIC', colnames(mc2))]
#
# ## By subject..
# mc2P <- lapply(fns_mc2, get_model_comparison, byParticipant=TRUE)
# mc2Psummary <- do.call(cbind, lapply(mc2P, function(x) x[2,]))
# colnames(mc2Psummary) <- paste0('ds', 1:4)
# mc2Psummary <- mc2Psummary[c(2,1,3,4),]  # align with group
# mc2combined <- cbind(round(mc2), mc2Psummary)
# mc2combined <- mc2combined[,c('ds1.delta_BPIC', 'ds1', 'ds2.delta_BPIC', 'ds2', 'ds3.delta_BPIC', 'ds3', 'ds4.delta_BPIC', 'ds4')]
# write.xlsx(mc2combined, file='./tables/mc2.xlsx', colNames = TRUE, rowNames = TRUE)
#
#
# # MC 3 --------------------------------------------------------------------
# fns_mc3 <- lapply(all_fns, function(x) x[grepl(include_fns_mc3,x)])
# mc3 <- lapply(fns_mc3, get_model_comparison, byParticipant=FALSE)
# names(mc3) <- paste0('ds', 1:4)
# mc3 <- do.call(cbind, mc3)
# row.names(mc3) <- c('z~SM, u~AH, b~V', 'z~SM, u~AH', 'z~SM, u~AH, u~V')
# mc3 <- mc3[c(2,1,3),]
# mc3[,grepl('delta_BPIC', colnames(mc3))]
#
# ## By subject. Terribly slow, but not sure why..
# mc3P <- lapply(fns_mc3, get_model_comparison, byParticipant=TRUE)
# mc3Psummary <- do.call(cbind, lapply(mc3P, function(x) x[2,]))
# colnames(mc3Psummary) <- paste0('ds', 1:4)
# mc3Psummary <- mc3Psummary[c(2,1,3),]
# mc3combined <- cbind(round(mc3), mc3Psummary)
# mc3combined <- mc3combined[,c('ds1.delta_BPIC', 'ds1', 'ds2.delta_BPIC', 'ds2', 'ds3.delta_BPIC', 'ds3', 'ds4.delta_BPIC', 'ds4')]
# write.xlsx(mc3combined, file='./tables/mc3.xlsx', colNames = TRUE, rowNames = TRUE)
#
#
# # MC 4 --------------------------------------------------------------------
# all_fns1 <- lapply(datasets, getFns, samples_dir='./samples_trends/')
# all_fns2 <- lapply(datasets, getFns, samples_dir='./samples/')
# all_fns4 <- lapply(1:4, function(x) c(all_fns1[[x]], all_fns2[[x]]))
# fns_mc4 <- lapply(all_fns4, function(x) x[grepl(include_fns_mc4,x)])
# mc4 <- lapply(fns_mc4, get_model_comparison, byParticipant=FALSE, hasTrend=TRUE)
# names(mc4) <- paste0('ds', 1:4)
# mc4 <- do.call(cbind, mc4)
# row.names(mc4) <- c('B', 'u', 'v', '0')
# mc4[,grepl('delta_BPIC', colnames(mc4))]
#
# ## By subject. Terribly slow, but not sure why..
# mc4P <- lapply(fns_mc4, get_model_comparison, byParticipant=TRUE, hasTrend=TRUE)
# mc4Psummary <- do.call(cbind, lapply(mc4P, function(x) x[2,]))
# colnames(mc4Psummary) <- paste0('ds', 1:4)
# row.names(mc4Psummary) <- c('B', 'u', 'v', '0')
# mc4combined <- cbind(round(mc4), mc4Psummary)
# mc4combined <- mc4combined[,c('ds1.delta_BPIC', 'ds1', 'ds2.delta_BPIC', 'ds2', 'ds3.delta_BPIC', 'ds3', 'ds4.delta_BPIC', 'ds4')]
# colnames(mc4combined) <- gsub('delta_', 'd', colnames(mc4combined))
# mc4combined <- mc4combined[c(4,1,2,3),]
# write.xlsx(mc4combined, file='./tables/mc4.xlsx', colNames = TRUE, rowNames = TRUE)
#
#
# # New MCs: All permutations --------------------------------------------------------
# processByParticipant <- function(mcP) {
#   byParticipant <- data.frame(matrix(NA, nrow=2,ncol=7))
#   colnames(byParticipant) <- c('zSM', 'vSM', 'uAH', 'bAH', 'vAH', 'uV', 'bV')
#   row.names(byParticipant) <- c('DIC', 'BPIC')
#   for(coln in colnames(byParticipant)) {
#     byParticipant[,coln] <- apply(mcP[,grepl(coln, colnames(mcP))],1,sum)
#   }
#   return(byParticipant)
# }
# # processGroup <- function(mc, sum_over_col='BPIC', sum_or_mean='sum') {
# #   out <- data.frame(matrix(NA, nrow=1,ncol=7))
# #   mechanisms <- c('zSM', 'vSM', 'uAH', 'bAH', 'vAH', 'uV', 'bV')
# #   colnames(out) <- mechanisms
# #   for(mechanism in mechanisms) {
# #     if(sum_or_mean == 'sum') out[1,mechanism] <- sum(mc[grepl(mechanism, row.names(mc)),sum_over_col])
# #     if(sum_or_mean == 'mean') out[1,mechanism] <- mean(mc[grepl(mechanism, row.names(mc)),sum_over_col])
# #   }
# #   out[1,c('zSM', 'vSM')] <- out[1,c('zSM', 'vSM')]-out[1,c('zSM')]
# #   out[1,c('uAH', 'bAH', 'vAH')] <- out[1,c('uAH', 'bAH', 'vAH')]-out[1,c('uAH')]
# #   out[1,c('uV', 'bV')] <- out[1,c('uV', 'bV')]-out[1,c('bV')]
# #   out
# # }
#
# # group wise: summed BPIC
# include_fns_mc5 <- '(v|z)SM(v|u|b)AH(u|b)V'
# all_fns <- lapply(datasets, getFns, samples_dir='./samples/')
# fns_mc5 <- lapply(all_fns, function(x) x[grepl(include_fns_mc5,x)])
# mc5 <- lapply(fns_mc5, get_model_comparison, byParticipant=FALSE)
# lapply(mc5, processGroup, sum_or_mean='sum')
#
# ## TMP
# names(mc5) <- paste0('ds', 1:4)
# rn <- mc5
# mc5c <- merge(do.call(cbind, mc5[1:3]), mc5[[4]], by=0, outer=TRUE, all=TRUE)
# row.names(mc5c) <- mc5c$Row.names
# # row.names(mc5) <- c('B', 'u', 'v', '0')
# mc5c[,grepl('delta_BPIC', colnames(mc5c))]
#
# # subject-wise: winners per mechanism. This is quite slow.
# all_fns <- lapply(datasets, getFns, samples_dir='./samples/')
# fns_mc5 <- lapply(all_fns, function(x) x[grepl(include_fns_mc5,x)])
# mc5P <- lapply(fns_mc5, get_model_comparison, byParticipant=TRUE)
# mc5Psummary <- lapply(mc5P, processByParticipant)
#
#
# mc5Psummary <- do.call(cbind, lapply(mc5P[1:3], function(x) x[2,]))
# mc5Psummary <- cbind(mc5Psummary,0)
# mc5Psummary[colnames(t(mc5P[[4]][2,])),4] <- t(mc5P[[4]][2,])
# colnames(mc5Psummary) <- paste0('ds', 1:4)
# # row.names(mc4Psummary) <- c('B', 'u', 'v', '0')
# mc5combined <- cbind(round(mc5c[,!colnames(mc5c)=='Row.names']), mc5Psummary)
# mc5combined <- mc5combined[,c('ds1.delta_BPIC', 'ds1', 'ds2.delta_BPIC', 'ds2', 'ds3.delta_BPIC', 'ds3', 'delta_BPIC', 'ds4')]
# row.names(mc5combined) <- gsub('vSM', 'dr~SM, ', gsub('zSM', 'db~SM, ', gsub('AH', '~AM, ', gsub('V', '~FM', gsub('v~', 'dv~', row.names(mc5combined))))))
# colnames(mc5combined) <- gsub('delta_', 'd', colnames(mc5combined))
# mc5combined
#
# write.xlsx(mc5combined, file='./tables/mc5_full.xlsx', colNames = TRUE, rowNames = TRUE)
#
#
#
# mechanism <- 'bV'
# sapply(1:4, function(x) model_averaging(IC_for = mc5[[x]][grepl(mechanism, row.names(mc5[[x]])),],
#                                         IC_against=mc5[[x]][!grepl(mechanism, row.names(mc5[[x]])),])[1,])
#
#
# uAM <- mc5[[3]][grepl('uAH', row.names(mc5[[1]])),]
# notUAMM <- mc5[[3]][!grepl('uAH', row.names(mc5[[1]])),]
# model_averaging(IC_for = uAM, IC_against=notUAMM)
#
#
# processByParticipant(mc5P[[2]])
#
#
#
#
#
#
#
#
#
#
# # old stuff below ---------------------------------------------------------
#
#
# fns <- getFns('wagenmakers2004_CS')
# fns <- fns[grepl(include_fns_mc5, fns)]
# mc5.1 <- get_model_comparison(fns, byParticipant = FALSE)    # ALL MODELS
# processGroup(mc5.1)
#
#
#
# fns <- getFns('wagenmakers2004_CS')
# fns <- fns[grepl('(v|z)SM(v|u|b)AH(u|b)V', fns)]
# mc <- get_model_comparison(fns, byParticipant = FALSE)    # ALL MODELS
# mcP1 <- get_model_comparison(fns, byParticipant = TRUE)     # ALL MODELS
# mcP1summary <- processByParticipant(mcP1)
#
#
# fns <- getFns('forstmann2008')
# fns <- fns[grepl('(v|z)SM(v|u|b)AH(u|b)V', fns)]
# mc <- get_model_comparison(fns, byParticipant = FALSE)     # ALL MODELS
# mcP2 <- get_model_comparison(fns, byParticipant = TRUE)     # ALL MODELS
# mcP2summary <- processByParticipant(mcP2)
#
# fns <- getFns('mileticvanmaanen2019exp2block2')
# fns <- fns[grepl('(v|z)SM(v|u|b)AH(u|b)V', fns)]
# mc <- get_model_comparison(fns, byParticipant = FALSE)    # ALL MODELS
# mcP3 <- get_model_comparison(fns, byParticipant = TRUE)    #
# mcP3summary <- processByParticipant(mcP3)
#
#
# fns <- getFns('wagenmakers2008exp2')
# fns <- fns[grepl('(v|z)SM(v|u|b)AH(u|b)V', fns)]
# mc <- get_model_comparison(fns, byParticipant = FALSE)    # ALL MODELS
#
#
#
#
#
#
# # out <- EMC2:::compare_subject(sList=list(vSMuAHbV=EMC2:::loadRData('./samples/wagenmakers2004_CS_model-vSMuAHbV_trend-NULL.RData'),
# #                                          zSMuAHbV=EMC2:::loadRData('./samples/wagenmakers2004_CS_model-zSMuAHbV_trend-NULL.RData'),
# #                                          zSMuAHuV=EMC2:::loadRData('./samples/wagenmakers2004_CS_model-zSMuAHuV_trend-NULL.RData')), print_summary = FALSE)
#
#
#
#
#
# # OLD BELOW ---------------------------------------------------------------
#
#
#
# getFns <- function(task, samples_dir='./samples/') {
#   all_fns <- Sys.glob(paste0(samples_dir, task, '*.RData'))
#   fns = all_fns[!grepl('_pp', all_fns)]
#   return(fns)
# }
#
#
# doMC1 <- function(task, byParticipant=FALSE) {
#   all_fns <- getFns(task)
#   # 1. Compare learning models, ignore DCT
#   print(paste0(task, ': MC 1'))
#   select_fns <- all_fns[grepl('RDM-(vSM|zSM|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
#   if(length(select_fns)>0) {
#     mc <- get_model_comparison(select_fns, byParticipant=byParticipant)
#   }
#   return(mc)
# }
#
# doMC2 <- function(task, byParticipant=FALSE) {
#   all_fns <- getFns(task)
#   # 2. Threshold vs urgency as a function of accuracy history
#   print(paste0(task, ': MC 2'))
#   select_fns <- all_fns[grepl('RDM-(zSMuAH|zSMaAH|zSMvAH|zSM|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
#   if(length(select_fns)>0) {
#     mc <- get_model_comparison(select_fns, byParticipant=byParticipant)
#   }
#   return(mc)
# }
#
# doMC3 <- function(task, subfilter=0, byParticipant=FALSE) {
#   all_fns <- getFns(task)
#   # 3. Threshold vs urgency as a function of response speed
#   print(paste0(task, ': MC 3'))
#   #  if(task == 'mileticvanmaanen2019exp2block2') {
#   select_fns <- all_fns[grepl('RDM-(zSMuAH|zSMaAHbV|zSMuAHbV|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
#   #  } else {
#   #    select_fns <- all_fns[grepl('RDM-(zSMuAH|zSMuAHbV|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
#   #  }
#   subfilter <- ifelse(task == 'wagenmakers2008exp2', 500, 0)
#   if(length(select_fns)>0) {
#     mc <- get_model_comparison(select_fns, subfilter=subfilter, byParticipant=byParticipant)
#   }
# }
#
# doMC4a <- function(task, byParticipant=FALSE) {
#   all_fns <- getFns(task, samples_dir='./samples/')
#   # 4a. Winning models with DCT, including null
#   print(paste0(task, ': MC 4'))
#   select_fns <- all_fns[grepl('-zSMuAHbV-trend-DCT-(B|u|v)-3', all_fns) | grepl('-zSMuAHbV-DCT-NULL', all_fns)]
#   if(length(select_fns)>0) {
#     mc <- get_model_comparison(select_fns, byParticipant=byParticipant)
#   }
# }
#
# doMC4b <- function(task) {
#   all_fns <- getFns(task, samples_dir='./samples/')
#   # 4b. Winning models with DCT, excluding null
#   print(paste0(task, ': MC 4'))
#   select_fns <- all_fns[grepl('-zSMuAHbV-trend-DCT-(B|u|v)-3', all_fns)] # | grepl('-zSMuAHbV-DCT-NULL', all_fns)]
#   if(length(select_fns)>0) {
#     mc <- get_model_comparison(select_fns)
#   }
# }
#
# process_by_participant <- function(x, nms) {
#   x <- t(do.call(rbind, lapply(x, function(y) y[2,])))
#   row.names(x) <- nms
#   colnames(x) <- paste('Dataset', 1:ncol(x))
#   x
# }
#
#
#
# ## 1  v~SM; z~SM
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2',  'wagenmakers2008exp2')[c(1,2,3,4)]
# allMCs1 <- setNames(lapply(tasks, doMC1), 1:length(tasks))
# allMCs1 <- lapply(allMCs1, get_delta_IC)
# mc1 <- do.call(cbind, allMCs1) # lapply(allMCs, function(x) round(x[,c('BPIC', 'wBPIC')],3)))
# row.names(mc1) <- c('NULL', 'v~SM', 'z~SM')
# write.csv(mc1, file='./tables/model_comparisons_tables/mc1_full.csv')
# write.xlsx(round(mc1,7), file='./tables/model_comparisons_tables/mc1_full.xlsx', colNames = TRUE, rowNames = TRUE)
#
# # By participant
# allMCs1_participant <- setNames(lapply(tasks, doMC1, byParticipant=TRUE), 1:length(tasks))
# allMCs1_participant_BPICs <- process_by_participant(allMCs1_participant, c('NULL', 'v~SM', 'z~SM'))
#
#
# ## 2  u~AH, b~AH, v~AH
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2',  'wagenmakers2008exp2')[c(1,2,3,4)]
# allMCs2 <- setNames(lapply(tasks, doMC2), 1:length(tasks))
# allMCs2 <- lapply(allMCs2, get_delta_IC)
# mc2 <- do.call(cbind, allMCs2)  #lapply(allMCs, function(x) round(x[,c('BPIC', 'wBPIC')],3)))
# row.names(mc2) <- c('NULL', 'z~SM, b~AH', 'z~SM', 'z~SM, u~AH', 'z~SM, v~AH')
# mc2 <- mc2[c(1,3,2,4,5),]
# write.csv(mc2, file='./tables/model_comparisons_tables/mc2_full.csv')
# write.xlsx(round(mc2,7), file='./tables/model_comparisons_tables/mc2_full.xlsx', colNames = TRUE, rowNames = TRUE)
#
# # By participant
# allMCs2_participant <- setNames(lapply(tasks, doMC2, byParticipant=TRUE), 1:length(tasks))
# allMCs2_participant_BPICs <- process_by_participant(allMCs2_participant, c('NULL', 'z~SM, b~AH', 'z~SM', 'z~SM, u~AH', 'z~SM, v~AH'))
# allMCs2_participant_BPICs <- allMCs2_participant_BPICs[c(1,3,2,4,5),]
#
#
# ## 3 b~RS
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')[c(1,2,3,4)]
# allMCs3 <- setNames(lapply(tasks, doMC3), 1:length(tasks))
# allMCs3 <- lapply(allMCs3, get_delta_IC)
# mc3 <- do.call(cbind, allMCs3) #, function(x) round(x[,c('BPIC', 'wBPIC')],3)))
# row.names(mc3) <- c('NULL', 'z~SM, b~AH, b~RS', 'z~SM, u~AH, b~RS', 'z~SM, u~AH')
# mc3 <- mc3[c(1,4,2,3),]
# write.csv(mc3, file='./tables/model_comparisons_tables/mc3_full.csv')
# write.xlsx(round(mc3,7), file='./tables/model_comparisons_tables/mc3_full.xlsx', colNames = TRUE, rowNames = TRUE)
#
# # By participant
# allMCs3_participant <- setNames(lapply(tasks, doMC3, byParticipant=TRUE), 1:length(tasks))
# allMCs3_participant <- process_by_participant(allMCs3_participant, c('NULL', 'z~SM, b~AH, b~RS', 'z~SM, u~AH, b~RS', 'z~SM, u~AH'))
# allMCs3_participant <- allMCs3_participant[c(1,4,2,3),]
#
#
# ## 4. Cosines based on winning models
# tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2',  'wagenmakers2008exp2')[c(1,2,3,4)]
# allMCs4 <- setNames(lapply(tasks, doMC4a), 1:length(tasks))
# allMCs4 <- lapply(allMCs4, get_delta_IC)
# mc4 <- do.call(cbind, allMCs4)   #lapply(allMCs, function(x) round(x[,c('BPIC', 'wBPIC', 'minD', 'Dmean')],1)))
# row.names(mc4) <- c('z~SM, u~AH, b~RS, b~DCT',
#                     'z~SM, u~AH, b~RS, no DCT',
#                     'z~SM, u~AH, b~RS, u~DCT',
#                     'z~SM, u~AH, b~RS, v~DCT')
# mc4 <- mc4[c(2,1,3,4),]
# write.csv(mc4, file='./tables/model_comparisons_tables/mc4_full.csv')
# write.xlsx(round(mc4,7), file='./tables/model_comparisons_tables/mc4_full.xlsx', colNames = TRUE, rowNames = TRUE)
#
# allMCs4_participant <- setNames(lapply(tasks, doMC4a, byParticipant=TRUE), 1:length(tasks))
# allMCs4_participant <- process_by_participant(allMCs4_participant, c('z~SM, u~AH, b~RS, b~DCT',
#                                                                      'z~SM, u~AH, b~RS, no DCT',
#                                                                      'z~SM, u~AH, b~RS, u~DCT',
#                                                                      'z~SM, u~AH, b~RS, v~DCT'))
# allMCs4_participant <- allMCs4_participant[c(2,1,3,4),]
#
#
#
#
# # Backwards sequential ----------------------------------------------------
# getFns <- function(task, samples_dir='./samples/') {
#   all_fns <- Sys.glob(paste0(samples_dir, task, '*.RData'))
#   fns = all_fns[!grepl('_pp', all_fns)]
#   return(fns)
# }
#
# to_include = 'SM(u|v|b)AH(b|u)V-DCT-'
# fns <- getFns('forstmann2008', './samples_alice/')
# fns <- fns[grepl(to_include, fns)]
# mc <- get_model_comparison(fns, byParticipant = FALSE)
#
#
# to_include = 'SM(u|v|b)AH(b|u)V-DCT-'
# fns <- getFns('wagenmakers2008exp2', './samples_alice/')
# fns <- fns[grepl(to_include, fns)]
# mc <- get_model_comparison(fns)
#
#
# fns <- getFns('wagenmakers2004_CS', './samples/')
# fns <- fns[grepl(to_include, fns)]
# mc <- get_model_comparison(fns)
#
#
#
#
#
# ## manual
# samplers_zSMuAHuV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMuAHuV-DCT-NULL.RData')
# samplers_zSMuAHbV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMuAHbV-DCT-NULL.RData')
#
# samplers_zSMbAHbV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMbAHbV-DCT-NULL.RData')
# samplers_zSMbAHuV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMbAHuV-DCT-NULL.RData')
#
# samplers_zSMvAHbV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMvAHbV-DCT-NULL.RData')
# samplers_zSMvAHuV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMvAHuV-DCT-NULL.RData')
#
# samplers_vSMuAHuV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-vSMuAHuV-DCT-NULL.RData')
# samplers_vSMuAHbV <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-vSMuAHbV-DCT-NULL.RData')
#
# samplers_zSMuAH <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMuAH-DCT-NULL.RData')
# samplers_zSMbAH <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMaAH-DCT-NULL.RData')
#
# slist <-list('zSMuAHuV'=samplers_zSMuAHuV,
#              'zSMuAHbV'=samplers_zSMuAHbV,
#              'zSMbAHbV'=samplers_zSMbAHbV,
#              'zSMbAHuV'=samplers_zSMbAHuV,
#
#              'zSMvAHbV'=samplers_zSMvAHbV,
#              'zSMvAHuV'=samplers_zSMvAHuV,
#
#              #                      'samplers_vSMuAHbV'=samplers_vSMuAHbV,
#              #                      'samplers_vSMuAHbV'=samplers_vSMuAHbV,
#
#              'zSMuAH'=samplers_zSMuAH,
#              'zSMbAH'=samplers_zSMbAH
# )
# compare_IC(sList=slist)
#
# #EMC2:::compare_ICs(sList=slist)
#
#
#
# samplers_zSMuAH <- EMC2:::loadRData('./samples/forstmann2008_model-RDM-zSMuAH-DCT-NULL.RData')
# samplers_zSMuAHuV <- EMC2:::loadRData('./samples/forstmann2008_model-RDM-zSMuAHuV-DCT-NULL.RData')
# samplers_zSMuAHbV <- EMC2:::loadRData('./samples/forstmann2008_model-RDM-zSMuAHbV-DCT-NULL.RData')
# compare_IC(sList=list('zSMuAH'=samplers_zSMuAH,
#                       'zSMuAHuV'=samplers_zSMuAHuV,
#                       'zSMuAHbV'=samplers_zSMuAHbV
# ))
#
# EMC2:::compare_ICs(sList=list('zSMuAHuV'=samplers_zSMuAHuV,
#                               'zSMuAHbV'=samplers_zSMuAHbV
# ))
#
#
#
# getFns <- function(task, samples_dir='./samples/') {
#   all_fns <- Sys.glob(paste0(samples_dir, task, '*.RData'))
#   fns = all_fns[!grepl('_pp', all_fns)]
#   return(fns)
# }
#
# fns <- getFns('forstmann2008', './samples_alice/')
# mc <- get_model_comparison(fns)
#
# fns <- getFns('mileticvanmaanen2019exp2block2', './samples_alice/')
# mc <- get_model_comparison(fns)
#
# fns <- getFns('wagenmakers2004_CS', './samples/')
# fns <- fns[grepl('.*SM(u|b|v)AH.*V-DCT', fns)]
# mc <- get_model_comparison(fns)
#
#
# fns <- getFns('forstmann2008', './samples_alice/')
# fns <- fns[grepl('.*SM(u|b|v)AH.*V-DCT', fns)]
# mc <- get_model_comparison(fns)
#
#
# fns <- getFns('mileticvanmaanen2019exp2block2', './samples_alice/')
# fns <- fns[grepl('.*SM(u|b|v|a)AH.*V-DCT', fns)]
# mc <- get_model_comparison(fns)
#
#
#
# ## kakscheiss
# add_learningmodel_to_env <- function(samplers, learningModel) {
#   list2env(list('learningModel'=learningModel), env=environment(attr(samplers[[1]]$data[[1]], 'adapt')$design$dynamic$output_fun))
# }
# samplers_old <- EMC2:::loadRData('./samples/samples_bigger_bound/wagenmakers2004_CS_model-RDM-zSMuAHuV-DCT-NULL.RData')
# add_learningmodel_to_env(samplers_old, 'zSMuAHuV')
#
# samplers_new <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMuAHuV-DCT-NULL.RData')
# list2env(list('learningModel'='zSMuAHuV'), env=environment(attr(samplers_new[[1]]$data[[1]], 'adapt')$design$dynamic$output_fun))
#
# samplers_bV_old <- EMC2:::loadRData('./samples/samples_bounded/wagenmakers2004_CS_model-RDM-zSMuAHbV-DCT-NULL.RData')
#
# samplers_bV_new <- EMC2:::loadRData('./samples/wagenmakers2004_CS_model-RDM-zSMuAHbV-DCT-NULL.RData')
# list2env(list('learningModel'='zSMuAHbV'),
#          env=environment(attr(samplers_bV_new[[1]]$data[[1]], 'adapt')$design$dynamic$output_fun))
#
# compare_IC(sList=list('old'=samplers_old,
#                       'new'=samplers_new,
#                       'oldbV'=samplers_bV_old,
#                       'newbV'=samplers_bV_new
# ))
#
# undebug(EMC2:::update_pars)
# debug(EMC2:::IC)
# undebug(EMC2:::log_likelihood_race)
# EMC2:::IC(samplers_new)
# #
# #attr(samplers_bV_new[[1]]$data[[1]], 'adapt')$design$dynamic$output_fun
#
# list2env(list('learningModel'='zSMuAHbV'),
#          env=environment(attr(samplers_bV_new[[1]]$data[[1]], 'adapt')$design$dynamic$output_fun))
#
# # plot_chains(samplers_bV_old, selection='LL', layout=c(3,4))
# # plot_chains(samplers_bV_new, selection='LL', layout=c(3,4))
# #
# # apply(EMC2:::parameters_data_frame(samplers_bV_old),2,mean)
# # apply(EMC2:::parameters_data_frame(samplers_bV_new),2,mean)
#
#
#
# ## old stuff
#
# mc1 <- lapply(fns_mc1, get_model_comparison, byParticipant=FALSE)
# names(mc1) <- paste0('ds', 1:4)
# mc1 <- do.call(cbind, mc1)
# row.names(mc1) <- c('NULL', 'dr~SM', 'db~SM')
# round(mc1[,grepl('dBPIC', colnames(mc1))])
#
# ## By subject..
# mc1P <- lapply(fns_mc1, get_model_comparison, byParticipant=TRUE)
# mc1Psummary <- do.call(cbind, lapply(mc1P, function(x) x[2,]))
# colnames(mc1Psummary) <- paste0('ds', 1:4, '.N')
# mc1combined <- cbind(round(mc1), mc1Psummary)
# mc1combined <- mc1combined[,c('ds1.dBPIC', 'ds1.N', 'ds2.dBPIC', 'ds2.N', 'ds3.dBPIC', 'ds3.N', 'ds4.dBPIC', 'ds4.N')]
# write.xlsx(mc1combined, file='./tables/mc1.xlsx', colNames = TRUE, rowNames = TRUE)
#
