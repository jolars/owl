#' Obtain Coefficients from Model fit by Golem
#'
#' This function is equivalent to simply calling `drop(object$coefficients)`
#' on a model fit by [golem()].
#'
#' @param object an object of class `'Golem'`.
#' @param ... ignored
#'
#' @return Coefficients from the model after having dropped extraneous
#'   dimensions by calling drop.
#'
#' @export
#' @examples
#' fit <- golem(mtcars$mpg, mtcars$vs, n_sigma = 1)
#' coef(fit)
coef.Golem <- function(object, ...) {
  drop(object$coefficients)
}
