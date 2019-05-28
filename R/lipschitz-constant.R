#' @include families.R penalties.R

setGeneric(
  "lipschitzConstant",
  function(family, penalty, x, fit_intercept)
    standardGeneric("lipschitzConstant")
)

setMethod(
  "lipschitzConstant",
  c("Gaussian", "Penalty"),
  function(family, penalty, x, fit_intercept) {
    max(rowSums(x^2)) + fit_intercept
  }
)

setMethod(
  "lipschitzConstant",
  c("Binomial", "Penalty"),
  function(family, penalty, x, fit_intercept) {
    0.25 * (max(rowSums(x^2)) + fit_intercept)
  }
)
