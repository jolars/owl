standardize <- function(x, standardize_features) {
    p <- NCOL(x)

    if (standardize_features) {
      x_center <- Matrix::colMeans(x)

      if (!inherits(x, "sparseMatrix")) {
        x <- sweep(x, 2, x_center)
        x_scale <- apply(x, 2, norm, type = "2")
        x_scale[x_scale == 0] <- 1
        x <- sweep(x, 2, x_scale, "/")
      } else {
        x_scale <- standardizedSparseColNorms(x, x_center)
        x_scale[x_scale == 0] <- 1
        for (j in seq_len(p))
          x[, j] <- x[, j]/x_scale[j]
      }
    } else {
      x_center <- rep.int(0, p)
      x_scale  <- rep.int(1, p)
    }

    list(x = x,
         x_center = x_center,
         x_scale = x_scale)
}

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

  list(intercepts = intercepts,
       betas = betas)
}
