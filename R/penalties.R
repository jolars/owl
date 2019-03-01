#' Sorted L-One Penalized Estimation (SLOPE)
#'
#' @param lambda penalty strength
#' @param sigma estimate of signal noise
#' @param fdr target for false discovery rate (FDR)
#'
#' @return A parameter pack for the SLOPE penalty.
#'
#' @export
slope <- function(lambda = c("gaussian", "bhq"),
                  sigma = NULL,
                  fdr = 0.2) {

  stopifnot(length(fdr) == 1,
            fdr >= 0 && fdr <= 1)

  # lambda (regularization strength)
  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)
    lambda <- 0
    #lambda <- rep.int(0, p)
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)
    #stopifnot(length(lambda) == p)

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")
  }

  if (is.null(sigma)) {
    sigma <- 1 # temporarily set sigma to 1
    sigma_type <- "auto"
  } else {
    stopifnot(length(sigma) == 1)
    sigma_type <- "user"
  }

  structure(list(name = "slope",
                 lambda = lambda,
                 lambda_type = lambda_type,
                 sigma = sigma,
                 sigma_type = sigma_type,
                 fdr = fdr),
            class = "penalty")
}

#' Elastic Net (Lasso and Ridge)
#'
#' The elastic net penalty is a mix of the lasso (L1) and ridge (L2)
#' regularization penalties. This mix is controlled by `alpha`, where
#' a value of 1 imposes the lasso penalty and 0 the ridge penalty, whereas
#' anything in between serves as a combination of the two.
#'
#' @param lambda the regularization strength, which can either be a
#'   user-supplied vector of (theoreticall) any length, or `NULL`, in which
#'   case golem automatically computes a sequence so that the first value
#'   leads to the null (intercept-only) model.
#' @param n_lambda the length of the \eqn{\lambda} sequence -- ignored if
#'   a value is supplied to `lambda`.
#' @param alpha the elastic net mix, 1 for the lasso penalty, 0 for the ridge,
#'   and anything in between for a mix of the two (the elastic net)
#'
#' @return A paramter pack for the elastic net penalty.
#' @export
elasticNet <- function(lambda = NULL,
                       lambda_min_ratio = 0.0001,
                       n_lambda = 100,
                       alpha = 1) {

  stopifnot(length(alpha) == 1,
            alpha >= 0 && alpha <= 1)

  # lambda (regularization strength)
  if (is.null(lambda)) {
    lambda_type <- "auto"
    lambda <- rep.int(0, n_lambda)
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)
    stopifnot(length(lambda) > 0)
  }

  structure(list(name = "elasticNet",
                 lambda = lambda,
                 lambda_min_ratio = 0.0001,
                 lambda_type = lambda_type,
                 n_lambda = length(lambda),
                 alpha = alpha),
            class = "penalty")
}
