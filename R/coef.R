#' Return Coefficients from a Golem Model
#'
#' Returns the coefficients for a model fit via [golem::golem()].
#' Unlike retrieving the coefficients using simple indexing, such as
#' `fit$coefficients`, this function attempts to reduce the dimensions of
#' the coefficients---that are otherwise always a three-dimensional array---by
#' calling [base::drop()].
#'
#' @param object an object of class `"Golem"`
#' @param ... parameters passed onto [stats::coef()]
#'
#' @return The coefficients, including the intercept if it was fit.
#' @export
#' @include golem.R
setMethod(
  "coef",
  "Golem",
  function(object, ...)
    drop(object@coefficients)
)
