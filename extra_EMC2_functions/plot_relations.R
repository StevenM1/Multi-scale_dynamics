plot_relations <- function(emc = NULL, stage = "sample", plot_cred = TRUE, plot_means = TRUE,
                           only_cred = FALSE, nice_names = NULL, ...)
{
  loadings <- NULL
  standardize <- TRUE
  do_corr <- TRUE
  corrs <- NULL
  optionals <- list(...)
  for (name in names(optionals)) {
    assign(name, optionals[[name]])
  }
  if (!is.null(loadings))
    do_corr <- F
  addCoef.col <- "black"
  if (!plot_means)
    addCoef.col <- NULL
  if (!is.null(emc))
    sampled <- merge_chains(emc)
  if (do_corr || !is.null(corrs)) {
    if (is.null(corrs)) {
      values <- sampled$samples$theta_var[, , sampled$samples$stage ==
                                            stage, drop = F]
    }
    else {
      values <- corrs
    }
    means <- cov2cor(apply(values, 1:2, mean))
  }
  else {
    if (is.null(loadings)) {
      if (standardize) {
        loadings <- standardize_loadings(emc, stage = stage)
      }
      else {
        samples <- merge_chains(emc)
        loadings <- samples$samples$theta_lambda[, ,
                                                 samples$samples$stage == stage]
      }
    }
    values <- loadings
    means <- apply(values, 1:2, mean)
  }
  if (plot_cred || only_cred) {
    if (do_corr) {
      for (i in 1:dim(values)[3]) {
        values[, , i] <- cov2cor(values[, , i])
      }
    }
    cred <- aperm(apply(values, 1:2, quantile, probs = c(0.025,
                                                         0.975)))
    conf <- paste0("[", format(cred[, , 1, drop = F], digits = 1),
                   ":", format(cred[, , 2, drop = F], digits = 1), "]")
    if (only_cred) {
      is_cred <- unique(c(which(cred[, , 1] > 0), which(cred[,
                                                             , 2] < 0)))
      conf[!1:length(conf) %in% is_cred] <- ""
      is_cred <- unique(c(which(t(cred[, , 1] > 0)), which(t(cred[,
                                                                  , 2] < 0))))
      means[!1:length(means) %in% is_cred] <- NA
    }
    cred <- round(cred, 2)
  }
  if (do_corr || !is.null(corrs)) {
    if (!is.null(nice_names))
      colnames(means) <- rownames(means) <- nice_names
    if (only_cred) {
      p.mat <- means
      p.mat[is.na(means)] <- 1
      p.mat[!is.na(means)] <- 0
      means[is.na(means)] <- 0
      corrplot(means, addCoef.col = addCoef.col, number.cex = 0.75,
               tl.col = "black", p.mat = p.mat, insig = "blank",
               sig.level = 0.05)
    }
    else {
      corrplot(means, addCoef.col = addCoef.col, number.cex = 0.75,
               tl.col = "black")
    }
  }
  else {
    cols <- diverging_hcl(200, palette = "Red-Green")
    if (!is.null(nice_names))
      rownames(means) <- nice_names
    if (only_cred) {
      p.mat <- means
      p.mat[is.na(means)] <- 1
      p.mat[!is.na(means)] <- 0
      means[is.na(means)] <- 0
      colnames(p.mat) <- colnames(means) <- 1:ncol(means)
      p <- corrplot(means, is.corr = F, col = cols, col.lim = c(-1,
                                                                1), cl.pos = "n", addCoef.col = addCoef.col,
                    number.cex = 0.75, tl.col = "black", p.mat = p.mat,
                    insig = "blank", sig.level = 0.05)
    }
    else {
      p <- corrplot(means, is.corr = F, col = cols, col.lim = c(-1,
                                                                1), cl.pos = "n", addCoef.col = addCoef.col,
                    number.cex = 0.75, tl.col = "black")
    }
    max_x <- max(p$corrPos$x)
    max_y <- max(p$corrPos$y)
    colorlegend(xlim = c(max_x + 1, max_x + 3), ylim = c(max_y/2 -
                                                           max_y/5, max_y/2 + max_y/5), cols, c(seq(-1, 1, 0.25)),
                align = "l", vertical = TRUE, addlabels = TRUE)
  }
  if (plot_cred) {
    xs <- row(matrix(cred[, , 1], nrow = ncol(means)))
    ys <- (ncol(matrix(cred[, , 1], nrow = ncol(means))) +
             1) - col(matrix(cred[, , 2], nrow = ncol(means))) -
      0.05
    if (plot_means) {
      ys <- ys - 0.1
      text(xs, ys, conf, pos = 1, cex = 0.55, font = 2)
    }
    else {
      text(xs, ys, conf, pos = 1, cex = 0.7, font = 2)
    }
  }
}
