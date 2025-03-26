## Generate all posterior predictives
rm(list=ls())
library(EMC2)
library(emcAdapt)
root_dir = file.path(Sys.getenv('HOME'), 'Projects', 'dynamicEAMsNewEMC')
source(file.path(root_dir, 'extra_EMC2_functions/adaptive.R'))
source(file.path(root_dir, 'extra_EMC2_functions/model_RDMdynamic.R'))
source(file.path(root_dir, 'extra_EMC2_functions/utils.R'))
source(file.path(root_dir, 'extra_EMC2_functions/make_data.R'))

cores_for_chains <- 3
cores_per_chain <- 10

## Generate posterior predictives of all most-relevant models
samplers_dir = './samples'
samplers_ms3 <- Sys.glob(file.path(samplers_dir, '*NULL.RData'))
samplers_ms3 <- samplers_ms3[!grepl('_pp-', samplers_ms3)]
samplers_ms3 <- samplers_ms3[grepl('(zSMuAHbV|zSMbAHbV|zSMvAHbV)', samplers_ms3)]

# add zSM, uAH, zSMuAH, bV
samplers_null <- Sys.glob(file.path('./samples_mc123', '*NULL.RData'))
samplers_null <- samplers_null[!grepl('_pp-', samplers_null)]
samplers_null <- samplers_null[grepl('model-(zSMuAH|zSM|uAH|bV|NULL)_', samplers_null)]

# add trends
samplers_trends <- Sys.glob(file.path('./samples_trends', '*.RData'))
samplers_trends <- samplers_trends[!grepl('_pp-', samplers_trends)]
samplers_trends <- c(samplers_trends[grepl('(wagenmakers2004_CS|wagenmakers2008exp2)_model-(zSMuAHbV|NULL)_trend-DCTB3', samplers_trends)],
                     samplers_trends[grepl('(mileticvanmaanen2019exp2block2|forstmann2008)_model-(zSMuAHbV|NULL)_trend-DCTv3', samplers_trends)],
                     samplers_trends[grepl('model-(uAH|bV)', samplers_trends)],
                     samplers_trends[grepl('mileticvanmaanen2019exp2block1', samplers_trends)]
)


# initial QAM models for supplementary figure 2
samplers_initialqam <- Sys.glob(file.path('./samples_initialQAM', '*NULL.RData'))
samplers_initialqam <- samplers_initialqam[!grepl('_pp-', samplers_initialqam)]
samplers_initialqam <- samplers_initialqam[grepl('model-(zSMuAHbV)_', samplers_initialqam)]

all_samplers <- c(samplers_ms3, samplers_initialqam, samplers_null, samplers_trends)

check_file_age <- function(file1, file2) {
  # Check if both files exist
  if (!file.exists(file1)) {
    stop(paste("File", file1, "does not exist"))
  }
  if (!file.exists(file2)) {
    stop(paste("File", file2, "does not exist"))
  }

  # Get the modification times of the files
  file1_time <- file.info(file1)$mtime
  file2_time <- file.info(file2)$mtime

  # Compare the modification times
  if (file1_time < file2_time) {
    return(TRUE)  # file1 is older
  } else {
    return(FALSE) # file1 is not older
  }
}

for(samplers_fn in all_samplers) {
  print(samplers_fn)
  pp_fn_root <- gsub('/samples', '/posteriorpredictives', samplers_fn)
  dn <- dirname(pp_fn_root)
  if(!dir.exists(dn)) dir.create(dn, recursive=TRUE)

  ## filenames
  fn_pp1 <- gsub('.RData', '_pp-unconditional.RData', pp_fn_root)
  fn_pp2 <- gsub('.RData', '_pp-conditional.RData', pp_fn_root)

  genPPs <- FALSE
  if(file.exists(fn_pp1)) {
    if(check_file_age(fn_pp1, samplers_fn)) {
      genPPs <- TRUE
    }
  }

  if(!(file.exists(fn_pp1) & file.exists(fn_pp2))) {
    gen_PPs <- TRUE
  }

  if(genPPs) {
    emc <- EMC2:::loadRData(samplers_fn)
    print(fn_pp1)
    if(!file.exists(fn_pp1)) {
      pp_unconditional <- predict(emc, n_post = 100, n_cores=cores_for_chains*cores_per_chain, conditionalOnData=FALSE)
      save(pp_unconditional, file=fn_pp1)
    }
    if(!file.exists(fn_pp2)) {
      pp_conditional <- predict(emc, n_post = 100, n_cores=cores_for_chains*cores_per_chain, conditionalOnData=TRUE)
      save(pp_conditional, file=fn_pp2)
    }
  }
}
