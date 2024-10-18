### SM Playing around with h5
# rm(list=ls())
# install.packages("BiocManager")
# BiocManager::install(c("rhdf5", 'HDF5Array'))
# library(EMC2)
# library(rhdf5)
# library(HDF5Array)

# Write *single* pmwg sampler to file
write_pmwg_samples <- function(samples, fn, groupName) {
  nms <- names(samples)
  for(nm in nms) {

    # matrices have dimensions
    dimensions <- dim(samples[[nm]])
    # vectors don't, use length
    if(is.null(dimensions)) dimensions <- length(samples[[nm]])

    # create dataset with dimensions and infinite max dimensions. This is important for extending later on
    dtype <- ifelse(nm=='stage', 'character', 'double')
    rhdf5::h5createDataset(file=fn, dataset=paste0(groupName, '/', nm), storage.mode=dtype, dims=dimensions, maxdims=rep(1e5, length(dimensions)))  #NB: Inf is better but doesn't work on all clusters :-/
                           #maxdims=rep(Inf, length(dimensions)))

    # then write
    rhdf5::h5write(obj=samples[[nm]], file=fn, name=paste0(groupName, '/', nm))

    # check for dimnames; if these exist, write
    if(!is.null(dimnames(samples[[nm]]))) {
      HDF5Array::h5writeDimnames(dimnames = dimnames(samples[[nm]]), filepath=fn, name=paste0(groupName, '/', nm))
    }
  }
}

write_pmwgs_samples_to_h5 <- function(samplers, fn) {
  # check for existing groups, create these if they don't exist
  n_chains <- length(samplers)
  groupNames <- paste0('chain', 1:n_chains)

  existingGroups <- rhdf5::h5ls(fn)$name
  if((length(existingGroups) == 0) | (!all(groupNames == existingGroups))) {
    sapply(groupNames, function(x) rhdf5::h5createGroup(file=fn, group=x))
    sapply(groupNames, function(x) rhdf5::h5createGroup(file=fn, group=paste0(x, '/samples')))
  }

  # write chains
  lapply(1:length(samplers), function(x) write_pmwg_samples(samplers[[x]]$samples, fn=fn, groupName=paste0('chain',x, '/samples')))

  # check for nuisance, write if exists
  if('sampler_nuis' %in% names(samplers[[1]])) {
    lapply(1:length(samplers), function(x) rhdf5::h5createGroup(file=fn, group=paste0('chain',x, '/sampler_nuis')))
    lapply(1:length(samplers), function(x) rhdf5::h5createGroup(file=fn, group=paste0('chain',x, '/sampler_nuis/samples')))
    lapply(1:length(samplers), function(x) write_pmwg_samples(samplers[[x]]$sampler_nuis$samples, fn=fn, groupName=paste0('chain',x, '/sampler_nuis/samples')))
  }
}


reconstruct_dimnames <- function(x) {
  nms <- names(x)
  nms_dimnames <- nms[grepl('_dimnames', nms)]
  nms_other <- nms[!grepl('_dimnames', nms)]

  for(i in nms_dimnames) {
    element_name <- gsub('(\\.|_dimnames)', '', i)
    dimnames(x[[element_name]]) <- x[[i]]
  }
  x <- x[nms_other]
  return(x)
}

read_pmwg_samples <- function(fn) {
  reconstructed_samples <- rhdf5::h5dump(fn)

  ## fix dimnames
  for(chain in 1:length(reconstructed_samples)) {
    reconstructed_samples[[chain]]$samples <- reconstruct_dimnames(reconstructed_samples[[chain]]$samples) # lapply(reconstructed_samples, reconstruct_dimnames)
  }

  # check for nuisance object, if exists, fix dimnames there as well
  if('sampler_nuis' %in% names(reconstructed_samples[[1]])) {
    for(chain in 1:length(reconstructed_samples)) {
      reconstructed_samples[[chain]]$sampler_nuis$samples <- reconstruct_dimnames(reconstructed_samples[[chain]]$sampler_nuis$samples)
    }
  }

  return(reconstructed_samples)
}


## start clean
shorten_chains <- function(samplers, select, remove = TRUE, last_select = FALSE, filter = NULL) {
  design_list <- attr(samplers, 'design_list')
  data_list <- attr(samplers, 'data_list')
  model_list <- attr(samplers, 'model_list')

  samplers2 <- lapply(samplers, remove_iterations, select, remove = remove, last_select = last_select, filter = filter)
  attr(samplers2, 'design_list') <- design_list
  attr(samplers2, 'data_list') <- data_list
  attr(samplers2, 'model_list') <- model_list

  return(samplers2)
}


get_dimensions_from_h5 <- function(fn, groupName, datasetName) {
  lst <- rhdf5::h5ls(fn)
  groupName <- paste0('/', groupName)
  dims <- lst[lst$group==groupName&lst$name==datasetName,'dim']
  dims <- strsplit(dims, ' x ')[[1]]
  return(as.numeric(dims))
}

# Next, we append
append_samples_to_h5 <- function(samples, nAppend, fn, idx=NULL, rootGroup='chain') {
  for(nm in names(samples)) {
    groupName <- paste0(rootGroup, '/', nm)

    ## for simplciity, assume all objects need to extended on the last dimension
    ## except for last_theta_var_inv and idx
    if(nm %in% c('last_theta_var_inv', 'idx')) {
      ## these can just be overwritten
      rhdf5::h5write(obj=samples[[nm]], file=fn, name=groupName)
    } else {
      currentDimensions <- get_dimensions_from_h5(fn, rootGroup, nm)

      ## set extent to be able to append
      dimensions <- dim(samples[[nm]])
      if(is.null(dimensions)) dimensions <- length(samples[[nm]])
      dimensions[length(dimensions)] <- currentDimensions[length(currentDimensions)] + nAppend #length(iters)

      # we always append to the last dimension
      indexToWrite <- lapply(1:length(dimensions), function(x) return(NULL))
      indexToWrite[[length(indexToWrite)]] <- (currentDimensions[length(currentDimensions)]+1):(currentDimensions[length(currentDimensions)] + nAppend) #length(iters)) #iters

      # set extent, write to index
      rhdf5::h5set_extent(file=fn, dataset=groupName, dims=dimensions)
      rhdf5::h5write(obj=samples[[nm]], file=fn, name=groupName, index=indexToWrite)
    }
  }
}

append_samplers_to_h5 <- function(samplers, select, fn, last_select=TRUE, remove=FALSE) {
  # select ONLY samples to append
  to_append <- shorten_chains(samplers, select=select, last_select = last_select, remove=remove)

  # how many samples are we appending?
  nAppend <- ifelse(length(select)==1, select, length(select))

  # append
  lapply(1:length(to_append), function(x) append_samples_to_h5(to_append[[x]]$samples, nAppend=nAppend, fn=fn, rootGroup = paste0('chain', x, '/samples')))

  # check for nuisance, append if exists
  if('sampler_nuis' %in% names(to_append[[1]])) {
    lapply(1:length(to_append), function(x) append_samples_to_h5(to_append[[x]]$sampler_nuis$samples, nAppend=nAppend, fn=fn, rootGroup = paste0('chain', x, '/sampler_nuis/samples')))
  }
}
# load('./pnassamples.RData')
# samples <- read_pmwg_samples('pnassamples.h5')
# samplers2 <- samplers
# samplers2[[1]]$samples <- samples$chain1$samples
# samplers2[[2]]$samples <- samples$chain2$samples
# samplers2[[3]]$samples <- samples$chain3$samples
# samplers2[[1]]$samples$idx <- rowSums(chain_n(samplers2))[1]+1
# samplers2[[2]]$samples$idx <- rowSums(chain_n(samplers2))[2]+1
# samplers2[[3]]$samples$idx <- rowSums(chain_n(samplers2))[3]+1
# Read HDF5 file and
#' HDF5 file and RData file to re-create single samplers file
#'
#' @return
#' @export
#'
#' @examples
loadH5 <- function(fn) {
  fn_RData <- gsub('\\.h5$', '.RData', fn)
  ### UNTESTED
  samplers <- loadRData(fn_RData)
  samples <- read_pmwg_samples(fn)

  for(i in 1:length(samplers)) {
    samplers[[i]]$samples <- samples[[i]]$samples
    samplers[[i]]$samples$idx <- rowSums(chain_n(samplers))[i]+1
    if(samplers[[i]]$samples$idx > 1) samplers[[i]]$init <- TRUE

    if('sampler_nuis' %in% names(samplers[[i]])) {
      samplers[[i]]$sampler_nuis$samples <- samples[[i]]$sampler_nuis$samples
      samplers[[i]]$sampler_nuis$samples$idx <- samplers[[i]]$samples$idx
    }
  }

  ## always make references to data, not copies of the same data
  for(i in 2:length(samplers)) samplers[[i]]$data <- samplers[[1]]$data

  return(samplers)
}


shorten_chains_lowmem <- function(samplers, max_burn_iters=500, max_sample_iters=1000, keep_single_sample_per_stage=TRUE, verbose=TRUE) {
  # max_burn_iters <- max_burn_iters#-already_removed
  # max_sample_iters <- max_sample_iters#-already_removed

  stage <- samplers[[1]]$samples$stage[length(samplers[[1]]$samples$stage)]
  nsamples <- rowSums(chain_n(samplers))[1]
  to_remove = c()

  max_iters <- ifelse(stage %in% c('preburn', 'burn', 'adapt'), max_burn_iters, max_sample_iters)
  if(nsamples < max_iters) {
    to_remove <- c()
  } else {
    idx <- 1:nsamples
    potentially_remove <- idx[1:(nsamples-max_iters)]
    # never remove adapt
    is_adapt <- which(samplers[[1]]$samples$stage=='adapt')
    to_remove <- setdiff(potentially_remove, is_adapt)
  }

  ## shorten chains
  # if(stage %in% c('preburn', 'burn')) {
  #   ## Keep maximally 500 samples of c('preburn', 'burn') in memory
  #   is_stage <- samplers[[1]]$samples$stage %in% c('preburn', 'burn')
  #   nsample <- sum(is_stage)
  #   if(nsample > max_burn_iters) {
  #     to_remove = which(is_stage)[1:(nsample-max_burn_iters)]
  #   }
  # } else if(stage == 'adapt') {
  #   ## always keep adapt, don't shorten this part
  #   ## but we can get rid of preburn now
  #   is_preburn <- samplers[[1]]$samples$stage == 'preburn'
  #
  #   if(sum(is_preburn) > 0) {
  #     to_remove <- c(to_remove, which(is_preburn))
  #   }
  #
  #   ## and we can shorten burn if (n_burn+n_preburn+n_adapt) > max_n_burn
  #   is_stage <- samplers[[1]]$samples$stage %in% c('preburn', 'burn', 'adapt')
  #   nsample <- sum(is_stage)
  #   if(nsample > max_burn_iters) {
  #     to_remove = which(is_stage)[1:(nsample-max_burn_iters)]
  #   }
  #
  # } else if(stage == 'sample') {
  #   ## remove preburn if it's still in (eg via a restart)
  #   is_preburn <- samplers[[1]]$samples$stage == 'preburn'
  #   if(sum(is_preburn) > 0) {
  #     to_remove <- c(to_remove, which(is_preburn))
  #   }
  #   ## Remove burn if that's still in (eg via a restart)
  #   is_burn <- samplers[[1]]$samples$stage == 'burn'
  #   if(sum(is_burn) > 0) {
  #     to_remove <- c(to_remove, which(is_burn))
  #   }
  #
  #   ## keep adapt and last 1e3 iterations of sample
  #   is_sample <- samplers[[1]]$samples$stage=='sample'
  #   nsample <- sum(is_sample)
  #   if(nsample > max_sample_iters) {
  #     to_remove = c(to_remove, which(is_sample)[1:(nsample-max_sample_iters)])
  #   }
  # }

  if(keep_single_sample_per_stage) {
    keep <- c()
    if(stage %in% c('burn', 'adapt', 'sample')) {
      in_stage <- which(samplers[[1]]$samples$stage == 'preburn')
      keep <- c(keep, in_stage[length(in_stage)])
    }
    if(stage %in% c('adapt', 'sample')) {
      in_stage <- which(samplers[[1]]$samples$stage == 'burn')
      keep <- c(keep, in_stage[length(in_stage)])
    }
    to_remove <- setdiff(to_remove, keep)
  }

  if(length(to_remove) == 0) {
    if(verbose) {
      print('Low-memory mode. Nothing to remove')
    }
  } else {
    if(verbose) {
      print(paste0('Low-memory mode. Removing ', length(to_remove), ' samples'))
    }
    samplers <- shorten_chains(samplers, select=to_remove, remove=TRUE, last_select=FALSE)
  }
  if(verbose) {
    print('After removal:')
    print(chain_n(samplers))
  }

  # if(verbose) {
  #   print(paste0('Removing ', length(to_remove), ' samples'))
  #   print('AFTER REMOVAL: ')
  #   print(chain_n(samplers))
  # }

  return(samplers)
}

chain_n_h5 <- function(fn) {
  if(!file.exists(fn)) return(NULL)
  if(!rhdf5::H5Fis_hdf5(fn)) return(NULL)
  lst <- rhdf5::h5ls(fn)

  groupName <- '/chain1/samples'
  datasetName <- 'stage'
  chains <- lst$name[grepl('chain(\\d+)', lst$name)]

  do.call(rbind, lapply(chains, function(x) table(factor(rhdf5::h5read(fn, paste0(x, '/samples/stage')), levels=c("preburn", "burn","adapt","sample")))))
}

append_and_reload <- function(samplers, fileNameH5, fileNameRData, step_size) {
  if(!file.exists(fileNameH5)) {
    ## make new, this is easy
    rhdf5::h5createFile(fileNameH5)
    write_pmwgs_samples_to_h5(samplers, fileNameH5)
  } else {
    if(!rhdf5::H5Fis_hdf5(fileNameH5)) {
      stop('Something is wrong, h5 file not recognized as such, corrupt? Blame Steven for this one')
    }
    ## file already exists, we can append the last step_size iterations
    # there's a very weird edge case where the first 1/3rd of the chains is removed (say, 33 samples of 100) on the very first
    # iteration. In that case, we try to add the last step_size (100) iterations, the first 33 of those are actually adapt, which
    # are already added... hence this weird thing
    append_samplers_to_h5(samplers, select=step_size, fn=fileNameH5, last_select=TRUE)

    # # here, we can load H5 & save as RData to update the RData file, which will then always contain the full samples
    samplers <- loadH5(fileNameH5)
    save(samplers, file=fileNameRData)
  }
  return(samplers)
}


# out <- shorten_chains_lowmem(samplers)
#chain_n(out)

#

# print(load('./samples/murphy2014_model-RDM-accOnV-DCT-B.RData'))
# fn = './h5samplers_testing.h5'
# file.remove(fn)            # start afresh
# h5createFile(fn)
#
#
# samplers_short <- shorten_chains(samplers, select=10:2584, remove=TRUE)
# write_pmwgs_samples_to_h5(samplers_short, fn)
# out <- read_pmwg_samples(fn)
#
#
# #debug(append_samplers_to_h5)
# append_samplers_to_h5(samplers, select=500:600, fn=fn)
# append_samplers_to_h5(samplers, select=601:700, fn=fn)
#
# out2 <- read_pmwg_samples(fn)
# str(out2$chain1$sampler_nuis$samples, max.level=1)
#
##
# to_append <- shorten_chains(samplers, select=11:20, remove=FALSE, overwriteIdx=FALSE)
# to_append[[1]]$samples$idx <- 20
# to_append[[2]]$samples$idx <- 20
# to_append[[3]]$samples$idx <- 20
# chain_n(to_append)


# samplers_short <- shorten_chains(samplers, select=100:2584)
# undebug(append_to_h5)
# append_to_h5(to_append, iters=11:20, fn=fn)



# ## write example array
# h5write(obj=samplers[[2]]$samples$theta_mu, file=fn, name='chain2/theta_mu')
# h5writeDimnames(dimnames = dimnames(samplers[[2]]$samples$theta_mu), filepath=fn, name='chain2/theta_mu')
#
# # h5writeDimnames
#
# ## load back
# theta_mu <- h5read(fn, 'chain2/theta_mu')
# dimnames(theta_mu) <- h5readDimnames(fn, 'chain2/theta_mu')



#str(samplers[[1]]$samples$theta_mu)
#str(h5read(fn, 'chain2/theta_mu'))
#str(h5readDimnames(fn, 'chain2/theta_mu'))
#
# h5createFile('myhdf5file.h5')
# h5createGroup('myhdf5file.h5', 'foo')
# h5createDataset("myhdf5file.h5", "foo/S", c(5,8),
#                 storage.mode = "integer", chunk=c(5,1), level=7)
# h5write(matrix(1:5,nr=5,nc=1), file="myhdf5file.h5",
#         name="foo/S", index=list(NULL,1))
# h5read("myhdf5file.h5", "foo/S")
# h5write(6:10, file="myhdf5file.h5",
#         name="foo/S", index=list(1,2:6))
# h5read("myhdf5file.h5", "foo/S")
# h5write(6:10, file="myhdf5file.h5",
#         name="foo/S", index=list(1,2:6+5))
# h5read("myhdf5file.h5", "foo/S")
#
#
#
#
#
#
#
# ##### appending
# fid <- H5Fcreate('test.h5')
# h5createGroup(fid,'1')
# h5createDataset(fid,'1/1', dims = c(2,2,2), maxdims = c(Inf,Inf,Inf), fillValue = NA)
#
# arr <- array(c(1:16),c(4,2,2))
# h5write(arr,fid,'1/1')
#
# h5set_extent(fid, '1/1', c(2,2,4))
#
# h5dump(fid)
