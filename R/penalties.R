Slope <- function(x,
                  y,
                  y_scale,
                  lambda = c("gaussian", "bhq"),
                  sigma = c("sequence", "estimate"),
                  sigma_min_ratio = NULL,
                  n_sigma = 100,
                  fdr = 0.2,
                  family) {

  n <- NROW(x)
  p <- NCOL(x)

  if (is.null(sigma_min_ratio))
    sigma_min_ratio <- if (n < p) 0.01 else 0.0001

  if (is.null(lambda))
    lambda <- "gaussian"

  # noise estimate
  if (is.character(sigma)) {
    sigma_type <- match.arg(sigma)
    sigma <- NA_real_
  } else {
    stopifnot(length(sigma) > 0, sigma >= 0, is.finite(sigma))
    sigma_type <- "user"
  }

  # regularization strength
  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)

    if (lambda_type %in% c("bhq", "gaussian")) {
      q <- 1:p * fdr/(2*p)
      lambda <- stats::qnorm(1 - q)

      if (lambda_type == "gaussian" && p > 1) {
        sum_sq <- 0
        for (i in 2:p) {
          sum_sq <- sum_sq + lambda[i - 1]^2
          w <- max(1, n - i)
          lambda[i] <- lambda[i]*sqrt(1 + sum_sq/w)
        }
      }

      # ensure non-increasing lambdas
      lambda[which.min(lambda):p] <- min(lambda)

    }
  } else {
    lambda <- as.double(lambda)

    if (length(lambda) != p)
      stop("lambda sequence must be as long as there are variables")

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")

    if (any(lambda < 0))
      stop("lambda sequence cannot contain negative values")
  }

  if (sigma_type == "sequence") {
    lambda_max <- lambdaMax(family, x, y, y_scale)*NROW(x)

    sigma <- lambda_max/min(lambda)
    sigma <- logSeq(sigma, sigma*sigma_min_ratio, n_sigma)
  }

  lambda <- matrix(lambda, p, 1)
  sigma  <- sigma

  structure(list(name = "slope",
                 tuning_parameters = c("sigma", "fdr"),
                 sigma = sigma,
                 fdr = fdr,
                 lambda = lambda),
            class = c("Slope", "Penalty"))
}

GroupSlope <- function(x,
                       y,
                       y_scale,
                       groups,
                       lambda = c("corrected", "mean", "max"),
                       sigma = c("sequence", "estimate"),
                       sigma_min_ratio = NULL,
                       n_sigma = 100,
                       fdr = 0.2,
                       family) {

  group_id <- groups$group_id
  ortho_group_id <- groups$ortho_group_id
  orthogonalize <- groups$orthogonalize
  wt <- groups$wt

  n <- NROW(x)
  p <- NCOL(x)

  n_groups <- length(group_id)

  if (is.null(lambda))
    lambda <- "corrected"

  if (is.null(sigma_min_ratio))
    sigma_min_ratio <- if (n < p) 0.01 else 0.0001

  group_sizes <- if (orthogonalize)
    lengths(ortho_group_id)
  else
    lengths(group_id)

  # noise estimate
  if (is.character(sigma)) {
    sigma_type <- match.arg(sigma)
    sigma <- NA_real_
  } else {
    stopifnot(length(sigma) > 0, sigma >= 0, is.finite(sigma))
    sigma_type <- "user"
  }

  # regularization strength
  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)

    if (lambda_type %in% c("max", "mean")) {

      lambda <- lambdaChiOrtho(fdr = fdr,
                               n.group = n_groups,
                               wt = wt,
                               group.sizes = group_sizes,
                               method = lambda_type)

    } else if (lambda_type == "corrected") {

      # Check for equal group sizes and equal weights
      if ((length(unique(group_sizes)) == 1) & (length(unique(wt)) == 1)) {
        # lambdas of Procedure 6 in Brzyski et. al. (2016)
        m <- unique(group_sizes)
        w <- unique(wt)

        lambda <- lambdaChiEqual(fdr = fdr,
                                 n.obs = n,
                                 n.group = n_groups,
                                 m = m,
                                 w = w)
      } else {
        # lambdas of Procedure 1 in Brzyski et. al. (2016)
        lambda <- lambdaChiMean(fdr = fdr,
                                n.obs = n,
                                n.group = n_groups,
                                group.sizes = group_sizes,
                                wt = wt)
      }
    }
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")

    if (any(lambda < 0))
      stop("lambda sequence cannot contain negative values")
  }

  if (sigma_type == "sequence") {
    lambda_max <- lambdaMax(family, x, y, y_scale)*NROW(x)

    sigma <- lambda_max/min(lambda)
    sigma <- logSeq(sigma, sigma*sigma_min_ratio, n_sigma)
  }

  lambda <- matrix(lambda, n_groups, 1)
  sigma <- sigma

  structure(list(name = "group_slope",
                 tuning_parameters = c("sigma", "fdr"),
                 sigma = sigma,
                 fdr = fdr,
                 lambda = lambda),
            class = c("GroupSlope", "Penalty"))
}

Lasso <- function(x,
                  y,
                  y_scale,
                  family,
                  lambda = NULL,
                  lambda_min_ratio = NULL,
                  n_lambda = 100) {

  if (is.null(lambda_min_ratio))
    lambda_min_ratio <- ifelse(NROW(x) < NCOL(x), 0.01, 0.0001)

  # lambda (regularization strength)
  if (is.null(lambda)) {
    lambda_max <- lambdaMax(family, x, y, y_scale)

    lambda <- logSeq(lambda_max,
                     lambda_max*lambda_min_ratio,
                     n_lambda)
    lambda_scale <- max(y_scale)

  } else {
    lambda <- as.double(lambda)
    stopifnot(length(lambda) > 0, all(lambda >= 0))
    lambda_scale <- 1
  }

  lambda_scale <- lambda_scale/NROW(x)
  lambda <- matrix(lambda, 1, length(lambda))

  structure(list(name = "lasso",
                 tuning_parameters = c("lambda"),
                 lambda = lambda,
                 lambda_scale = lambda_scale),
            class = c("Lasso", "Penalty"))
}
