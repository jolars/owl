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
