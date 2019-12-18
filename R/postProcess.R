postProcess <- function(object, ...) {
  UseMethod("postProcess", object)
}

postProcess.Penalty <- function(object,
                                intercepts,
                                betas,
                                x,
                                y,
                                fit_intercept,
                                x_center,
                                x_scale,
                                y_center,
                                y_scale) {

  res <- unstandardize(intercepts,
                       betas,
                       x,
                       y,
                       fit_intercept,
                       x_center,
                       x_scale,
                       y_center,
                       y_scale)

  intercepts <- res$intercepts
  betas <- res$betas
  nonzeros <- apply(betas, c(2, 3), function(x) abs(x) > 0)

  list(intercepts = intercepts,
       betas = betas,
       nonzeros = nonzeros)
}

