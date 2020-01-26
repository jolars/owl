#' Model deviance
#'
#' @param object an object of class `'Owl'`.
#' @param ... ignored
#'
#' @return For Gaussian models this is twice the residual sums of squares. For
#'   all other models, two times the negative loglikelihood is returned.
#' @export
#'
#' @examples
#' fit <- owl(heart$x, heart$y, family = "binomial")
#' deviance(fit)
deviance.Owl <- function(object, ...) {
  deviance_ratio <- object$deviance_ratio
  null_deviance <- object$null_deviance

  (1 - deviance_ratio) * null_deviance
}
