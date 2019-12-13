# Compute the usual unbiased estimate of the variance in a linear model.
estimate_noise <- function(x, y, intercept = TRUE) {
  n <- NROW(x)

  if (intercept)
    x <- cbind(rep(1, n), x)

  p <- NCOL(x)
  stopifnot(n > p)

  fit <- stats::lm.fit(x, y)
  sqrt(sum(fit$residuals^2) / (n - p))
}

firstUpper <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

camelCase <- function(x) {
  s <- strsplit(x, "[^[:alnum:]]")

  sapply(s, function(y) {
    first <- toupper(substring(y, 1, 1))
    paste(first, substring(y, 2), sep = "", collapse = "")
  })
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

# rowNorms <- function(x, norm_type = 2) {
#   if (inherits(x, "sparseMatrix"))
#     rowNormsSparse(x, norm_type)
#   else
#     rowNormsDense(x, norm_type)
# }

randomProblem <-
  function(n = 1000,
           p = 100,
           q = 0.2,
           n_groups = NULL,
           n_targets = if (match.arg(response) == "multinomial") 3 else 1,
           density = 1,
           amplitude = if (match.arg(response) == "poisson") 1 else 3,
           sigma = 1,
           response = c("gaussian", "binomial", "poisson", "multinomial")) {
  m <- n_targets

  if (density == 1) {
    x <- matrix(stats::rnorm(n*p), n)
  } else {
    x <- Matrix::rsparsematrix(n, p, density)
  }

  if (!is.null(n_groups)) {
    groups <- rep(seq_len(n_groups), each = ceiling(m*p/n_groups),
                  length.out = p*m)
    nonzero <- which(groups %in% seq_len(max(floor(n_groups*q), 1)))
  } else {
    groups <- NA
    nonzero <- sample(p*m, max(floor(q*p*m), 1))
  }

  signs <- sample(c(-1, 1), p*m, replace = TRUE)

  beta <- signs * amplitude * (1:(p*m) %in% nonzero)

  y <- switch(match.arg(response),
              gaussian = x %*% beta + stats::rnorm(n, sd = sigma),
              binomial = {
                y <- x %*% beta + stats::rnorm(n, sd = sigma)
                (sign(y) + 1)/2
              },
              multinomial = {
                beta <- matrix(beta, p, m)
                lin_pred_exp <- exp(x %*% beta)
                prob <- lin_pred_exp/rowSums(lin_pred_exp)
                y <- apply(prob, 1, function(x) sample(1:m, 1, prob = x))
              },
              poisson = {
                lambda <- as.double(exp(x %*% beta))
                y <- stats::rpois(n, lambda)
              })

  dimnames(x) <- list(seq_len(nrow(x)),
                      paste0("V", seq_len(ncol(x))))

  list(x = x,
       y = as.double(y),
       beta = beta,
       groups = groups,
       nonzero = nonzero,
       q = q)
}
