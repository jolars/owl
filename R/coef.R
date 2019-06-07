#' Obtain Coefficients from Model fit by Golem
#'
#' This function is equivalent to calling `$coef()` on a model fit by
#' [golem()].
#'
#' @param object an object of class `'Golem'`.
#' @param ... ignored
#'
#' @return Coefficients from the model after having dropped extraneous
#'   dimensions by calling drop.
#' @export
coef.Golem <- function(object, ...) {
  object$coef()
}
