#' @include penalties.R

unstandardize <- function(intercepts,
                          betas,
                          x,
                          y,
                          fit_intercept,
                          x_center,
                          x_scale,
                          y_center,
                          y_scale) {

  p <- NROW(betas)
  m <- NCOL(betas)
  n_penalties <- NSLICE(betas)

  for (k in seq_len(m)) {
    x_bar_beta_sum <- double(n_penalties)

    for (j in seq_len(p)) {
      betas[j, k, ] <- betas[j, k, ] * y_scale[k]/x_scale[j]
      x_bar_beta_sum <- x_bar_beta_sum + x_center[j] * betas[j, k, ]
    }

    if (fit_intercept)
      intercepts[, k, ] <-
        intercepts[, k, ]*y_scale[k] + y_center[k] - x_bar_beta_sum
  }

  list(intercepts = intercepts, betas = betas)
}

setGeneric(
  "postProcess",
  function(penalty,
           intercepts,
           betas,
           x,
           y,
           fit_intercept,
           x_center,
           x_scale,
           y_center,
           y_scale)
    standardGeneric("postProcess")
)

setMethod(
  "postProcess",
  "Penalty",
  function(penalty,
           intercepts,
           betas,
           x,
           y,
           fit_intercept,
           x_center,
           x_scale,
           y_center,
           y_scale) {

    c(intercepts, betas) %<-% unstandardize(intercepts,
                                            betas,
                                            x,
                                            y,
                                            fit_intercept,
                                            x_center,
                                            x_scale,
                                            y_center,
                                            y_scale)
    selected <- which(betas > 0)
    list(intercepts = intercepts, betas = betas, selected = selected)
  }
)

setMethod(
  "postProcess",
  "GroupSlope",
  function(penalty,
           intercepts,
           betas,
           x,
           y,
           fit_intercept,
           x_center,
           x_scale,
           y_center,
           y_scale) {

    group_id <- penalty@group_id
    n_groups <- length(group_id)
    ortho_group_id <- penalty@ortho_group_id
    orthogonalize <- penalty@orthogonalize

    if (orthogonalize) {

      group_lengths <- lengths(group_id)
      ortho_group_lengths <- lengths(ortho_group_id)

      ortho <- penalty@ortho

      if (all(group_lengths == ortho_group_lengths)) {

        beta_tilde <- betas

        for (i in seq_len(n_groups)) {
          # c corresponds to the (reordered) group
          # structure in orthogonalized version of X
          ci <- betas[ortho_group_id[[i]]]
          li <- ortho_group_lengths[i]

          if (inherits(ortho[[i]]$R, "sparseMatrix")) {
            bi <- tryCatch({
              Matrix::solve(ortho[[i]]$R, ci)
            }, error = function(err) {
              warning(paste("golem caught an error:", err))
              rep(NA, li)
            })
          } else {
            bi <- tryCatch({
              backsolve(ortho[[i]]$R, ci)
            }, error = function(err) {
              warning(paste("golem caught an error:", err))
              rep(NA, li)
            })
          }

          or <- double(li)

          for (j in seq_len(li))
            or[j] <- which(ortho[[i]]$P == j)

          # beta corresponds to the group structure in the original matrix
          beta_tilde[group_id[[i]]] <- bi[or]
        }
      } else {
        beta_tilde <- NULL
      }
    } else {
      beta_tilde <- betas
    }

    # compute group norms ||X_I beta_I||
    group_norms <- double(n_groups)

    for (i in seq_len(n_groups)) {
      if (orthogonalize) {
        group_norms[i] <- norm(as.matrix(betas[ortho_group_id[[i]]]), "f")
      } else {
        xbetai <- x[, group_id[[i]]] %*% as.matrix(beta_tilde[group_id[[i]]])
        group_norms[i] <- norm(as.matrix(xbetai), "f")
      }
    }

    group_names <- names(group_id)
    names(group_norms) <- group_names

    res <- unstandardize(intercepts,
                         beta_tilde,
                         x,
                         y,
                         fit_intercept,
                         x_center,
                         x_scale,
                         y_center,
                         y_scale)
    res$selected <- which(group_norms > 0)
    res
  }
)

