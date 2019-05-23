# updates attributes given a list of attributes
updateAttributes <- function(x, val) {
  attributes(x) <- utils::modifyList(attributes(x), val, keep.null = TRUE)
  x
}

# Compute the usual unbiased estimate of the variance in a linear model.
estimate_noise <- function(x, y, intercept = TRUE) {
  n <- nrow(x)
  if (intercept)
    x <- cbind(rep(1, n), x)
  p <- ncol(x)
  stopifnot(n > p)

  fit <- stats::lm.fit(x, y)
  sqrt(sum(fit$residuals^2) / (n-p))
}

# Generate a random synthetic model and data. Used for testing.
random_problem <- function(n, p, k = NULL, amplitude = 3, sigma = 1) {
  if (is.null(k))
    k <- max(1, as.integer(p/5))

  x <- matrix(stats::rnorm(n*p), n, p)
  nonzero <- sample(p, k)
  beta <- amplitude * (1:p %in% nonzero)
  y <- x %*% beta + stats::rnorm(n, sd = sigma)
  list(x = x, y = y, beta = beta, nonzero = nonzero)
}

firstUpper <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

logSeq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

groupID <- function(group) {
  group <- as.integer(group)
  id <- unique(group)

  members <- lapply(id, function(id_i) which(group %in% id_i))
  names(members) <- id

  members
}

NSLICE <- function(x) {
  d <- dim(x)
  n_slices <- if (length(d) == 3) d[3] else 1
  n_slices
}
