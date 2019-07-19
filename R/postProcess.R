postProcess <- function(object, ...) {
  UseMethod("postProcess", object)
}

postProcess.Penalty <- function(object,
                                intercepts,
                                betas,
                                x,
                                y,
                                fit_intercept,
                                x_center,
                                x_scale,
                                y_center,
                                y_scale,
                                groups) {

  res <- unstandardize(intercepts,
                       betas,
                       x,
                       y,
                       fit_intercept,
                       x_center,
                       x_scale,
                       y_center,
                       y_scale)

  intercepts <- res$intercepts
  betas <- res$betas
  nonzeros <- apply(betas, c(2, 3), function(x) abs(x) > 0)

  list(intercepts = intercepts,
       betas = betas,
       nonzeros = nonzeros)
}

postProcess.GroupSlope <- function(object,
                                   intercepts,
                                   betas,
                                   x,
                                   y,
                                   fit_intercept,
                                   x_center,
                                   x_scale,
                                   y_center,
                                   y_scale,
                                   groups) {

  # compute group norms ||X_I beta_I||
  group_id <- groups$group_id
  ortho_group_id <- groups$ortho_group_id
  n_groups <- length(group_id)
  orthogonalize <- groups$orthogonalize

  group_norms <- apply(betas, c(2, 3), function(beta_i) {
    group_norms <- double(n_groups)

    for (i in seq_len(n_groups)) {
      if (orthogonalize) {
        group_norms[i] <- norm(as.matrix(beta_i[ortho_group_id[[i]]]), "f")
      } else {
        xbetai <- x[, group_id[[i]]] %*% as.matrix(beta_i[group_id[[i]]])
        group_norms[i] <- norm(as.matrix(xbetai), "f")
      }
    }
    group_norms
  })

  dimnames(group_norms) <- list(names(group_id), NULL, NULL)

  res <- unstandardize(intercepts,
                       betas,
                       x,
                       y,
                       fit_intercept,
                       x_center,
                       x_scale,
                       y_center,
                       y_scale)

  res$nonzeros <- apply(group_norms, c(2, 3), function(x) abs(x) > 0)
  res
}
