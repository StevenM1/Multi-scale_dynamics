## MINOR BUGFIX
compare <- function(sList,stage="sample",filter=NULL,use_best_fit=TRUE,
         BayesFactor = TRUE, cores_for_props =4, cores_per_prop = 1,
         print_summary=TRUE,digits=0,digits_p=3, ...) {
  if(is(sList, "emc")) sList <- list(sList)
  getp <- function(IC) {
    IC <- -(IC - min(IC))/2
    exp(IC)/sum(exp(IC))
  }
  if(!is.null(list(...)$subject)) subject <- list(...)$subject else subject <- NULL

  if (is.numeric(filter)) defaultsf <- filter[1] else defaultsf <- 0
  sflist <- as.list(setNames(rep(defaultsf,length(sList)),names(sList)))
  if (is.list(filter)) for (i in names(filter))
    if (i %in% names(sflist)) sflist[[i]] <- filter[[i]]
  dots <- EMC2:::add_defaults(list(...), group_only = FALSE)
  ICs <- setNames(vector(mode="list",length=length(sList)),names(sList))
  for (i in 1:length(ICs)) ICs[[i]] <- EMC2:::IC(sList[[i]],stage=stage,
                                                filter=sflist[[i]],use_best_fit=use_best_fit,subject=subject,print_summary=FALSE,
                                                group_only = dots$group_only)
  ICs <- data.frame(do.call(rbind,ICs))
  DICp <- getp(ICs$DIC)
  BPICp <- getp(ICs$BPIC)
  out <- cbind.data.frame(DIC=ICs$DIC,wDIC=DICp,BPIC=ICs$BPIC,wBPIC=BPICp,ICs[,-c(1:2)])

  if(BayesFactor){
    MLLs <- numeric(length(sList))
    for(i in 1:length(MLLs)){
      MLLs[i] <- EMC2:::run_bridge_sampling(sList[[i]], stage = stage, filter = sflist[[i]], both_splits = FALSE,
                                     cores_for_props = cores_for_props, cores_per_prop = cores_per_prop)
    }
    MD <- -2*MLLs
    modelProbability <- getp(MD)
    out <- cbind.data.frame(MD = MD, wMD = modelProbability, out)
  }
  if (print_summary) {
    tmp <- out
    tmp$wDIC <- round(tmp$wDIC,digits_p)
    tmp$wBPIC <- round(tmp$wBPIC,digits_p)
    if(BayesFactor){
      tmp$wMD <- round(tmp$wMD, digits_p)
      tmp[,-c(2,4,6)] <- round(tmp[,-c(2,4,6)],digits=digits)
    } else{
      tmp[,-c(2,4)] <- round(tmp[,-c(2,4)],digits=digits)
    }
    print(tmp)
  }
  invisible(out)
}

assignInNamespace("compare", compare, ns="EMC2")

IC <- function (emc, stage = "sample", filter = 0, use_best_fit = TRUE,
          print_summary = TRUE, digits = 0, subject = NULL, group_only = FALSE)
{
  ll <- get_pars(emc, stage = stage, filter = filter, selection = "LL",
                 merge_chains = TRUE)
  minDs <- -2 * apply(ll[[1]][[1]], 2, min)
  mean_lls <- apply(ll[[1]][[1]], 2, mean)
  alpha <- get_pars(emc, selection = "alpha", stage = stage,
                    filter = filter, by_subject = TRUE, merge_chains = TRUE)
  mean_pars <- lapply(alpha, function(x) {
    apply(do.call(rbind, x), 2, mean)
  })
  data <- emc[[1]]$data
  mean_pars_lls <- setNames(numeric(length(mean_pars)), names(mean_pars))
  if(is.null(subject)) {
    for (sub in names(mean_pars)) {
      mean_pars_lls[sub] <- EMC2:::calc_ll_manager(t(mean_pars[[sub]]),
                                            dadm = data[[sub]], emc[[1]]$model)
    }
  } else {
    mean_pars_lls[subject[1]] <- EMC2:::calc_ll_manager(t(mean_pars[[subject[1]]]),
                                          dadm = data[[subject[1]]], emc[[1]]$model)
  }
  Dmeans <- -2 * mean_pars_lls
  if (!is.null(subject)) {
    Dmeans <- Dmeans[subject[1]]
    mean_lls <- mean_lls[subject[1]]
    minDs <- minDs[subject[1]]
  }
  else {
    group_stats <- EMC2:::group_IC(emc, stage = stage, filter = filter,
                            type = emc[[1]]$type)
    if (group_only) {
      mean_lls <- group_stats$mean_ll
      minDs <- group_stats$minD
      Dmeans <- group_stats$Dmean
    }
    else {
      mean_lls <- c(mean_lls, group_stats$mean_ll)
      minDs <- c(minDs, group_stats$minD)
      Dmeans <- c(Dmeans, group_stats$Dmean)
    }
  }
  if (use_best_fit)
    minDs <- pmin(minDs, Dmeans)
  mD <- sum(-2 * mean_lls)
  Dmean <- sum(Dmeans)
  minD <- sum(minDs)
  if (!use_best_fit)
    Dm <- Dmean
  else Dm <- minD
  pD <- mD - Dm
  DIC <- mD + pD
  BPIC <- mD + 2 * pD
  out <- c(DIC = DIC, BPIC = BPIC, EffectiveN = pD, meanD = mD,
           Dmean = Dmean, minD = minD)
  names(out) <- c("DIC", "BPIC", "EffectiveN", "meanD", "Dmean",
                  "minD")
  if (print_summary)
    print(round(out, digits))
  invisible(out)
}

assignInNamespace("IC", IC, ns="EMC2")
