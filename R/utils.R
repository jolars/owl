# Compute the usual unbiased estimate of the variance in a linear model.
estimate_noise <- function(x, y, intercept = TRUE) {
  n_samples <- NROW(x)

  if (intercept)
    x <- cbind(rep(1, n_samples), x)

  n_features <- NCOL(x)
  stopifnot(n_samples > n_features)

  fit <- stats::lm.fit(x, y)
  sqrt(sum(fit$residuals^2) / (n_samples - n_features))
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

colNorms <- function(x, norm_type = 2) {
  if (inherits(x, "sparseMatrix"))
    colNormsSparse(x, norm_type)
  else
    colNormsDense(x, norm_type)
}

rowNorms <- function(x, norm_type = 2) {
  if (inherits(x, "sparseMatrix"))
    rowNormsSparse(x, norm_type)
  else
    rowNormsDense(x, norm_type)
}

randomProblem <- function(n = 1000,
                          p = 100,
                          q = 0.2,
                          n_groups = NULL,
                          density = 1,
                          amplitude = 3,
                          sigma = 1,
                          response = c("gaussian", "binomial")) {
  if (density == 1) {
    x <- matrix(stats::rnorm(n*p), n)
  } else {
    x <- Matrix::rsparsematrix(n, p, density)
  }

  if (!is.null(n_groups)) {
    groups <- rep(seq_len(n_groups), each = ceiling(p/n_groups),
                  length.out = p)
    nonzero <- which(groups %in% seq_len(max(floor(n_groups*q), 1)))
  } else {
    groups <- NA
    nonzero <- sample(p, max(floor(q*p), 1))
  }

  beta <- amplitude * (1:p %in% nonzero)

  y <- switch(match.arg(response),
              gaussian = x %*% beta + stats::rnorm(n, sd = sigma),
              binomial = {
                y <- x %*% beta + stats::rnorm(n, sd = sigma)
                (sign(y) + 1)/2
              })

  list(x = x,
       y = as.double(y),
       beta = beta,
       groups = groups,
       nonzero = nonzero,
       q = q)
}
