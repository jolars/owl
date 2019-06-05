#' @include families.R penalties.R

setGeneric(
  "lipschitzConstant",
  function(family,
           penalty,
           x,
           fit_intercept,
           x_center,
           x_scale,
           standardize_features)
    standardGeneric("lipschitzConstant")
)

setMethod(
  "lipschitzConstant",
  c("Gaussian", "Penalty"),
  function(family,
           penalty,
           x,
           fit_intercept,
           x_center,
           x_scale,
           standardize_features) {

    maxSquaredRowNorm(x, x_center/x_scale, standardize_features) + fit_intercept
  }
)

setMethod(
  "lipschitzConstant",
  c("Binomial", "Penalty"),
  function(family,
           penalty,
           x,
           fit_intercept,
           x_center,
           x_scale,
           standardize_features) {

    0.25 * (maxSquaredRowNorm(x,
                              x_center/x_scale,
                              standardize_features) + fit_intercept)
  }
)

