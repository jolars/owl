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
  list(X = X, y = y, beta = beta)
}

firstUpper <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
