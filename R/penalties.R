setClass("Penalty",
         slots = c(name = "character",
                   lambda = "numeric",
                   lambda_type = "character"),
         prototype = list(name = NA_character_,
                          lambda = NA_real_,
                          lambda_type = NA_character_))

setClass("Slope",
         contains = "Penalty",
         slots = c(sigma = "numeric",
                   fdr = "numeric"),
         prototype = list(sigma = NA_real_,
                          fdr = NA_real_))

setClass("Lasso",
         contains = "Penalty",
         slots = c(lambda_min_ratio = "numeric",
                   n_lambda = "numeric",
                   lambda_scale = "numeric"),
         prototype = list(lambda_min_ratio = NA_real_,
                          n_lambda = NA_integer_,
                          lambda_scale = NA_real_))

#' Sorted L-One Penalized Estimation (SLOPE)
#'
#' @param lambda penalty strength
#' @param sigma estimate of signal noise
#' @param fdr target for false discovery rate (FDR)
#'
#' @return A parameter pack for the SLOPE penalty.
#'
#' @export
Slope <- function(lambda = c("gaussian", "bhq"),
                  sigma = NULL,
                  fdr = 0.2) {

  stopifnot(length(fdr) == 1,
            fdr >= 0 && fdr <= 1)

  # noise estimate
  if (is.null(sigma)) {
    sigma <- 1 # temporarily set sigma to 1
    sigma_type <- "auto"
  } else {
    stopifnot(length(sigma) == 1)
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
  }

  new("Slope",
      name = "slope",
      sigma = sigma,
      lambda = lambda,
      lambda_type = lambda_type,
      fdr = fdr)
}

#' Lasso
#'
#' The lasso penalty penalized coefficients via the L1 norm, which induces
#' sparse solutions if the regularization strength (lambda) is sufficiently
#' strong.
#'
#' @param lambda the regularization strength, which can either be a
#'   user-supplied vector of (theoreticall) any length, or `NULL`, in which
#'   case golem automatically computes a sequence so that the first value
#'   leads to the null (intercept-only) model.
#' @param lambda_min_ratio the lowest permissible value of the
#'   lambda penalty
#' @param n_lambda the length of the \eqn{\lambda} sequence -- ignored if
#'   a value is supplied to `lambda`.
#'
#' @return A paramter pack for the lasso penalty.
#' @export
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

setMethod("setup",
          "Slope",
          function(object, family, x, y, y_scale) {

  lambda      <- object@lambda
  lambda_type <- object@lambda_type
  sigma       <- object@sigma
  fdr         <- object@fdr

  n <- NROW(x)
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
})

setMethod("setup",
          "Lasso",
          function (object, family, x, y, y_scale) {

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
})
