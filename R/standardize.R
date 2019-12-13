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
  n_penalties <- dim(betas)[3]

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
