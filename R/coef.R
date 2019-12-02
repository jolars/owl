#' Obtain Coefficients from Model fit by Golem
#'
#' This function is equivalent to simply calling `drop(object$coefficients)`
#' on a model fit by [golem()].
#'
#' @param object an object of class `'Golem'`.
#' @param ... arguments that are passed on to [stats::update()] (and therefore
#'   also to [golem()]) if `exact = TRUE` and the given penalty
#'   is not in `object`
#' @inheritParams predict.Golem
#'
#' @return Coefficients from the model after having dropped extraneous
#'   dimensions by calling drop.
#'
#' @export
#' @examples
#' fit <- golem(mtcars$mpg, mtcars$vs, n_sigma = 1)
#' coef(fit)
coef.Golem <- function(object,
                       lambda = NULL,
                       sigma = NULL,
                       exact = FALSE,
                       simplify = TRUE,
                       ...) {
  beta <- object$coefficients

  # p <- NROW(beta)
  # m <- NCOL(beta)
  n_penalties <- dim(beta)[3]

  if (object$penalty$name %in% c("slope", "group_slope")) {
    penalty <- object$sigma
    value <- sigma
  } else {
    penalty <- object$lambda
    value <- lambda
  }

  if (is.null(value)) {
    n_penalties <- length(penalty)
  } else if (all(value %in% penalty)) {
    n_penalties <- length(value)
    beta <- beta[, , penalty %in% value, drop = FALSE]
  } else if (exact) {
    object <- stats::update(object, lambda = lambda, sigma = sigma, ...)
    beta <- object$coefficients
  } else {
    stopifnot(value >= 0)
    interpolation_list <- interpolatePenalty(penalty, value)
    beta <- interpolateCoefficients(beta, interpolation_list)
    n_penalties <- length(value)
  }

  if (simplify)
    beta <- drop(beta)

  beta
}
