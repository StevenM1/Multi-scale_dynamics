## some model comparisons
library(EMC2)
library(emcAdapt)
#library(stringr)
#library(forecast)
library(openxlsx)

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
#       allGds <- c(allGds, max(EMC2:::check_gd(samplers, stage='sample', max_gd=1.1, mean_gd=1.05, trys=1, verbose=FALSE)$gd))
    }
  }
  if(length(allSamplers)>0) {
    mc = EMC2:::compare_IC(allSamplers, print_summary=FALSE, subfilter=subfilter)
    mc <- cbind(mc, max.gd=allGds)
    print(round(mc,3))
    mc
  }
}

getFns <- function(task, samples_dir='./samples/') {
  all_fns <- Sys.glob(paste0(samples_dir, task, '*.RData'))
  fns = all_fns[!grepl('_pp', all_fns)]
  return(fns)
}

get_delta_IC <- function(mc) {
  mc$delta_DIC <- mc$DIC-min(mc$DIC)
  mc$delta_BPIC <- mc$BPIC-min(mc$BPIC)
  return(mc)
}

doMC1 <- function(task) {
  all_fns <- getFns(task)
  # 1. Compare learning models, ignore DCT
  print(paste0(task, ': MC 1'))
  select_fns <- all_fns[grepl('RDM-(vSM|zSM|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
  if(length(select_fns)>0) {
    mc <- get_model_comparison(select_fns)
  }
  return(mc)
}

doMC2 <- function(task) {
  all_fns <- getFns(task)
  # 2. Threshold vs urgency as a function of accuracy history
  print(paste0(task, ': MC 2'))
  select_fns <- all_fns[grepl('RDM-(zSMuAH|zSMaAH|zSMvAH|zSM|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
  if(length(select_fns)>0) {
    mc <- get_model_comparison(select_fns)
  }
  return(mc)
}

doMC3 <- function(task, subfilter=0) {
  all_fns <- getFns(task)
  # 3. Threshold vs urgency as a function of response speed
  print(paste0(task, ': MC 3'))
#  if(task == 'mileticvanmaanen2019exp2block2') {
    select_fns <- all_fns[grepl('RDM-(zSMuAH|zSMaAHbV|zSMuAHbV|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
#  } else {
#    select_fns <- all_fns[grepl('RDM-(zSMuAH|zSMuAHbV|NULL)-', all_fns) & grepl('-DCT-NULL', all_fns)]
#  }
  subfilter <- ifelse(task == 'wagenmakers2008exp2', 500, 0)
  if(length(select_fns)>0) {
    mc <- get_model_comparison(select_fns, subfilter=subfilter)
  }
}

doMC4a <- function(task) {
  all_fns <- getFns(task, samples_dir='./samples/')
  # 4a. Winning models with DCT, including null
  print(paste0(task, ': MC 4'))
  select_fns <- all_fns[grepl('-zSMuAHbV-trend-DCT-(B|u|v)-3', all_fns) | grepl('-zSMuAHbV-DCT-NULL', all_fns)]
  if(length(select_fns)>0) {
    mc <- get_model_comparison(select_fns)
  }
}

doMC4b <- function(task) {
  all_fns <- getFns(task, samples_dir='./samples/')
  # 4b. Winning models with DCT, excluding null
  print(paste0(task, ': MC 4'))
  select_fns <- all_fns[grepl('-zSMuAHbV-trend-DCT-(B|u|v)-3', all_fns)] # | grepl('-zSMuAHbV-DCT-NULL', all_fns)]
  if(length(select_fns)>0) {
    mc <- get_model_comparison(select_fns)
  }
}



## 1  v~SM; z~SM
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2',  'wagenmakers2008exp2')[c(1,2,3,4)]
allMCs1 <- setNames(lapply(tasks, doMC1), 1:length(tasks))
allMCs1 <- lapply(allMCs1, get_delta_IC)
mc1 <- do.call(cbind, allMCs1) # lapply(allMCs, function(x) round(x[,c('BPIC', 'wBPIC')],3)))
row.names(mc1) <- c('NULL', 'v~SM', 'z~SM')
write.csv(mc1, file='./tables/model_comparisons_tables/mc1_full.csv')
write.xlsx(round(mc1,7), file='./tables/model_comparisons_tables/mc1_full.xlsx', colNames = TRUE, rowNames = TRUE)

## 2  u~AH, b~AH, v~AH
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2',  'wagenmakers2008exp2')[c(1,2,3,4)]
allMCs2 <- setNames(lapply(tasks, doMC2), 1:length(tasks))
allMCs2 <- lapply(allMCs2, get_delta_IC)
mc2 <- do.call(cbind, allMCs2)  #lapply(allMCs, function(x) round(x[,c('BPIC', 'wBPIC')],3)))
row.names(mc2) <- c('NULL', 'z~SM, b~AH', 'z~SM', 'z~SM, u~AH', 'z~SM, v~AH')
mc2 <- mc2[c(1,3,2,4,5),]
write.csv(mc2, file='./tables/model_comparisons_tables/mc2_full.csv')
write.xlsx(round(mc2,7), file='./tables/model_comparisons_tables/mc2_full.xlsx', colNames = TRUE, rowNames = TRUE)


## 3 b~RS
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2', 'wagenmakers2008exp2')[c(1,2,3,4)]
allMCs3 <- setNames(lapply(tasks, doMC3), 1:length(tasks))
allMCs3 <- lapply(allMCs3, get_delta_IC)
mc3 <- do.call(cbind, allMCs3) #, function(x) round(x[,c('BPIC', 'wBPIC')],3)))
row.names(mc3) <- c('NULL', 'z~SM, b~AH, b~RS', 'z~SM, u~AH, b~RS', 'z~SM, u~AH')
mc3 <- mc3[c(1,4,2,3),]
write.csv(mc3, file='./tables/model_comparisons_tables/mc3_full.csv')
write.xlsx(round(mc3,7), file='./tables/model_comparisons_tables/mc3_full.xlsx', colNames = TRUE, rowNames = TRUE)


## 4. Cosines based on winning models
tasks <- c('wagenmakers2004_CS', 'forstmann2008', 'mileticvanmaanen2019exp2block2',  'wagenmakers2008exp2')[c(1,2,3,4)]
allMCs4 <- setNames(lapply(tasks, doMC4a), 1:length(tasks))
allMCs4 <- lapply(allMCs4, get_delta_IC)
mc4 <- do.call(cbind, allMCs4)   #lapply(allMCs, function(x) round(x[,c('BPIC', 'wBPIC', 'minD', 'Dmean')],1)))
row.names(mc4) <- c('z~SM, u~AH, b~RS, b~DCT',
                    'z~SM, u~AH, b~RS, no DCT',
                    'z~SM, u~AH, b~RS, u~DCT',
                    'z~SM, u~AH, b~RS, v~DCT')
mc4 <- mc4[c(2,1,3,4),]
write.csv(mc4, file='./tables/model_comparisons_tables/mc4_full.csv')
write.xlsx(round(mc4,7), file='./tables/model_comparisons_tables/mc4_full.xlsx', colNames = TRUE, rowNames = TRUE)

