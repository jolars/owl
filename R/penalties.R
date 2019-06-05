#' @include families.R utils.R

setClass("Penalty",
         slots = c(name = "character",
                   lambda = "numeric",
                   lambda_type = "character"))

setClass("Slope",
         contains = "Penalty",
         slots = c(sigma = "numeric",
                   fdr = "numeric"))

Slope <- function(lambda = c("gaussian", "bhq"),
                  sigma = NULL,
                  fdr = 0.2) {

  if (is.null(lambda))
    lambda <- "gaussian"

  stopifnot(length(fdr) == 1,
            fdr >= 0 && fdr <= 1)

  # noise estimate
  if (is.null(sigma)) {
    sigma_type <- "auto"
    sigma <- NA_real_
  } else {
    stopifnot(length(sigma) == 1, sigma >= 0, is.finite(sigma))
    sigma_type <- "user"
  }

  # regularization strength
  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)
    lambda <- NA_real_
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")

    if (any(lambda < 0))
      stop("lambda sequence cannot contain negative values")
  }

  new("Slope",
      name = "slope",
      sigma = sigma,
      lambda = lambda,
      lambda_type = lambda_type,
      fdr = fdr)
}

setClass("GroupSlope",
         contains = "Penalty",
         slots = c(sigma = "numeric",
                   fdr = "numeric",
                   orthogonalize = "logical",
                   n_obs = "numeric",
                   groups = "numeric",
                   group_id = "list",
                   ortho = "list",
                   ortho_groups = "numeric",
                   ortho_group_id = "list",
                   ortho_group_length = "numeric",
                   wt = "numeric",
                   wt_per_coef = "numeric"))

GroupSlope <- function(groups,
                       lambda = c("corrected", "mean", "max"),
                       sigma = NULL,
                       fdr = 0.2,
                       orthogonalize = TRUE) {

  if (is.null(lambda))
    lambda <- "corrected"

  stopifnot(length(fdr) == 1, fdr >= 0, fdr <= 1)

  if (anyNA(groups))
    stop("NA values not allowed in 'groups'")

  # noise estimate
  if (is.null(sigma)) {
    sigma_type <- "auto"
    sigma <- NA_real_
  } else {
    stopifnot(length(sigma) == 1, sigma >= 0, is.finite(sigma))
    sigma_type <- "user"
  }

  # regularization strength
  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)
    lambda <- NA_real_
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")

    if (any(lambda < 0))
      stop("lambda sequence cannot contain negative values")
  }

  group_id <- groupID(groups)

  new("GroupSlope",
      name = "group_slope",
      sigma = sigma,
      lambda = lambda,
      lambda_type = lambda_type,
      fdr = fdr,
      n_obs = NA_real_,
      groups = groups,
      group_id = group_id,
      ortho_groups = groups,
      ortho_group_id = list(),
      orthogonalize = orthogonalize)
}

setClass("Lasso",
         contains = "Penalty",
         slots = c(lambda_min_ratio = "numeric",
                   n_lambda = "numeric",
                   lambda_scale = "numeric"))

Lasso <- function(lambda = NULL,
                  lambda_min_ratio = 0.0001,
                  n_lambda = 100) {
  # lambda (regularization strength)
  if (is.null(lambda)) {
    lambda_type <- "auto"
    lambda <- NA_real_

  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)
    n_lambda <- length(lambda)
    stopifnot(n_lambda > 0,
              all(lambda >= 0))
  }

  new("Lasso",
      name = "lasso",
      lambda = lambda,
      lambda_min_ratio = lambda_min_ratio,
      lambda_type = lambda_type,
      n_lambda = n_lambda)

}

setGeneric("setup",
           function(object, family, x, y, y_scale) standardGeneric("setup"))

setMethod(
  "setup",
  "Slope",
  function(object, family, x, y, y_scale) {

    lambda      <- object@lambda
    lambda_type <- object@lambda_type
    sigma       <- object@sigma
    fdr         <- object@fdr

    n  <- NROW(x)
    p <- NCOL(x)

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
    } else {
      if (length(lambda) != p)
        stop("lambda sequence must be as long as there are variables")
    }

    object@lambda <- lambda
    object@sigma  <- sigma
    object
  }
)

setMethod(
  "setup",
  "GroupSlope",
  function(object, family, x, y, y_scale) {

    lambda        <- object@lambda
    lambda_type   <- object@lambda_type
    sigma         <- object@sigma
    fdr           <- object@fdr
    orthogonalize <- object@orthogonalize
    group_id      <- object@group_id

    n_groups      <- length(group_id)

    n  <- NROW(x)
    wt <- object@wt

    group_sizes <- if (orthogonalize)
      lengths(object@ortho_group_id)
    else
      lengths(object@group_id)

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

    } else {
      if (length(lambda) != n_groups)
        stop("lambda sequence must be as long as there are variables")
    }

    object@lambda         <- lambda
    object@sigma          <- sigma
    # object@ortho_groups   <- ortho_groups
    # object@ortho_group_id <- ortho_group_id

    object
  }
)

setMethod(
  "setup",
  "Lasso",
  function(object, family, x, y, y_scale) {

    lambda           <- object@lambda
    lambda_min_ratio <- object@lambda_min_ratio
    lambda_type      <- object@lambda_type
    n_lambda         <- object@n_lambda

    # much of the scaling here is done to ensure that output
    # is equivalent to glmnet
    if (lambda_type == "auto") {
      lambda_max <- lambdaMax(family, x, y, y_scale)

      lambda <- logSeq(lambda_max, lambda_max*lambda_min_ratio, n_lambda)
      lambda_scale <- max(y_scale)
    } else {
      lambda_scale <- 1
    }

    object@lambda       <- lambda
    object@lambda_scale <- lambda_scale/NROW(x)
    object
  }
)

setGeneric(
  "getWeights",
  function(penalty, x) standardGeneric("getWeights")
)

setMethod(
  "getWeights",
  "Penalty",
  function(penalty, x) {
    rep.int(1, NCOL(x))
  }
)

setMethod(
  "getWeights",
  "GroupSlope",
  function(penalty, x) {
    penalty@wt_per_coef
  }
)
