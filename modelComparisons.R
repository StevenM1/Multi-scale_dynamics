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
out3 <- make_MC_table(fns_mc3, row_names=c('db~SM, u~AM, b~FM', 'db~SM, u~AM, u~FM', 'db~SM, u~AM'),
                      row_order=c(3,1,2))
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


