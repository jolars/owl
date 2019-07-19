lipschitzConstant <- function(object, ...) {
  UseMethod("lipschitzConstant", object)
}

lipschitzConstant.Gaussian <- function(object,
                                       x,
                                       fit_intercept,
                                       x_center,
                                       x_scale,
                                       standardize_features) {

  maxSquaredRowNorm(x, x_center/x_scale, standardize_features) +
    fit_intercept
}

lipschitzConstant.Binomial <- function(object,
                                       x,
                                       fit_intercept,
                                       x_center,
                                       x_scale,
                                       standardize_features) {
    0.25 * (maxSquaredRowNorm(x, x_center/x_scale, standardize_features) +
              fit_intercept)
}
