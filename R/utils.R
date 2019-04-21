# Compute the usual unbiased estimate of the variance in a linear model.
estimate_noise <- function(X, y, intercept = TRUE) {
  n <- nrow(X)
  if (intercept)
    X <- cbind(rep(1, n), X)
  p <- ncol(X)
  stopifnot(n > p)

  fit <- stats::lm.fit(X, y)
  sqrt(sum(fit$residuals^2) / (n-p))
}

# Generate a random synthetic model and data. Used for testing.
random_problem <- function(n, p, k = NULL, amplitude = 3, sigma = 1) {
  if (is.null(k))
    k <- max(1, as.integer(p/5))

  X <- matrix(stats::rnorm(n*p), n, p)
  nonzero <- sample(p, k)
  beta <- amplitude * (1:p %in% nonzero)
  y <- X %*% beta + stats::rnorm(n, sd = sigma)
  list(X = X, y = y, beta = beta, nonzero = nonzero)
}

firstUpper <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

logSeq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

unstandardize <- function(intercepts,
                          betas,
                          x_center,
                          x_scale,
                          y_center,
                          y_scale,
                          fit_intercept) {
  p <- NROW(betas)
  m <- NCOL(betas)
  K <- dim(betas)[3]

  for (k in seq_len(m)) {
    x_bar_beta_sum <- double(K)

    for (j in seq_len(p)) {
      betas[j, k, ] <- betas[j, k, ] * y_scale[k]/x_scale[j]
      x_bar_beta_sum <- x_bar_beta_sum + x_center[j] * betas[j, k, ]
    }

    if (fit_intercept) {
      intercepts[m, k, ] <-
        intercepts[m, k, ]*y_scale[k] + y_center[k] - x_bar_beta_sum
    }
  }

  list(intercepts = intercepts, betas = betas)
}
