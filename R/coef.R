#' Obtain Coefficients from Model fit by Owl
#'
#' This function is equivalent to simply calling `drop(object$coefficients)`
#' on a model fit by [owl()].
#'
#' @param object an object of class `'Owl'`.
#' @param ... arguments that are passed on to [stats::update()] (and therefore
#'   also to [owl()]) if `exact = TRUE` and the given penalty
#'   is not in `object`
#' @inheritParams predict.Owl
#'
#' @return Coefficients from the model after having dropped extraneous
#'   dimensions by calling drop.
#'
#' @export
#' @examples
#' fit <- owl(mtcars$mpg, mtcars$vs, n_sigma = 1)
#' coef(fit)
coef.Owl <- function(object,
                     lambda = NULL,
                     sigma = NULL,
                     exact = FALSE,
                     simplify = TRUE,
                     ...) {
  beta <- object$coefficients

  n_penalties <- dim(beta)[3]

  penalty <- object$sigma
  value <- sigma

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
