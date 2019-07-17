link <- function(object, y) {
  UseMethod("link", object)
}

link.Gaussian <- function(object, x) {
  x %*% object$beta + object$intercept
}

link.Binomial <- function(object, x) {
  lin_pred <- x %*% object$beta + object$intercept

  1 / (1 + exp(-lin_pred))
}
