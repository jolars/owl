#' Return linear predictors
#'
#' @param object an ojbect of class `"Golem"`
#' @param x a feature matrix
#'
#' @return Linear predictors
#'
#' @keywords internal
linearPredictors <- function(object, x) {

  beta <- object$coefficients

  if (is.numeric(which)) {
    beta <- beta[, , which, drop = FALSE]
  }

  p <- NROW(beta)
  n <- NROW(x)
  m <- NCOL(beta)
  n_penalties <- dim(beta)[3]

  if (inherits(x, "sparseMatrix"))
    x <- methods::as(x, "dgCMatrix")

  if (inherits(x, "data.frame"))
    x <- as.matrix(x)

  if (names(beta[, , 1])[1] == "(Intercept)") {
    x <- methods::cbind2(1, x)
  }

  stopifnot(p == NCOL(x))

  lin_pred <- array(dim = c(n, m, n_penalties),
                    dimnames = list(NULL,
                                    dimnames(beta)[[2]],
                                    dimnames(beta)[[3]]))

  for (i in seq_len(n_penalties)) {
    lin_pred[, , i] <- as.matrix(x %*% beta[, , i])
  }

  lin_pred
}

#' Generate predictions from golem models
#'
#' Return predictions from models fit by [golem()] based on
#' new data.
#'
#' @param object an object of class `"Golem"`, typically the result of
#'   a call to [golem()]
#' @param x new data
#' @param type type of prediction; `"link"` returns the linear predictors,
#'   `"response"` returns the result of applying the link function,
#'    and `"class"` returns class predictions.
#' @param ... ignored and only here for method consistency
#'
#' @seealso [stats::predict()], [stats::predict.glm()]
#'
#' @return Predictions from the model with scale determined by `type`.
#'
#'
#' @examples
#' fit <- with(mtcars, golem(cbind(mpg, hp), vs, family = "binomial"))
#' predict(fit, with(mtcars, cbind(mpg, hp)), type = "class")
#'
#' @export
predict.Golem <- function(object, x, type, ...) {
  NULL
}

#' @rdname predict.Golem
#' @export
predict.GolemGaussian <- function(object,
                                  x,
                                  type = c("link", "response"),
                                  ...) {
  type <- match.arg(type)

  linearPredictors(object, x)
}

#' @rdname predict.Golem
#' @export
predict.GolemBinomial <- function(object,
                                  x,
                                  type = c("link", "response", "class"),
                                  ...) {

  type <- match.arg(type)

  lin_pred <- linearPredictors(object, x)

  switch(
    type,
    response = 1 / (1 + exp(-lin_pred)),
    class = {
      cnum <- ifelse(lin_pred > 0, 2, 1)
      clet <- object$class_names[cnum]

      if (is.matrix(cnum))
        clet <- array(clet, dim(cnum), dimnames(cnum))

      clet
    }
  )
}
