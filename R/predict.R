#' Generate predictions from owl models
#'
#' Return predictions from models fit by [owl()].
#'
#' @param object an object of class `"owl"`, typically the result of
#'   a call to [owl()]
#' @param x new data
#' @param type type of prediction; `"link"` returns the linear predictors,
#'   `"response"` returns the result of applying the link function,
#'    and `"class"` returns class predictions.
#' @param ... ignored and only here for method consistency
#' @param lambda penalty parameter for Lasso-type models; if `NULL`, the
#'   values used in the original fit will be used
#' @param sigma penalty parameter for SLOPE models; if `NULL`, the
#'   values used in the original fit will be used
#' @param simplify if `TRUE`, [base::drop()] will be called before returning
#'   the coefficients to drop extraneous dimensions
#' @param exact if `TRUE` and the given parameter values differ from those in
#'   the original fit, the model will be refit by calling [stats::update()] on
#'   the object with the new parameters. If `FALSE`, the predicted values
#'   will be based on interpolated coefficients from the original
#'   penalty path.
#'
#' @seealso [stats::predict()], [stats::predict.glm()]
#'
#' @return Predictions from the model with scale determined by `type`.
#'
#'
#' @examples
#' fit <- with(mtcars, owl(cbind(mpg, hp), vs, family = "binomial"))
#' predict(fit, with(mtcars, cbind(mpg, hp)), type = "class")
#'
#' @export
predict.Owl <- function(object,
                        x,
                        lambda = NULL,
                        sigma = NULL,
                        type = "link",
                        exact = FALSE,
                        simplify = TRUE,
                        ...) {
  # This method (the base method) only generates linear predictors

  if (inherits(x, "sparseMatrix"))
    x <- methods::as(x, "dgCMatrix")

  if (inherits(x, "data.frame"))
    x <- as.matrix(x)

  beta <- stats::coef(object, lambda = lambda, sigma = sigma, simplify = FALSE)

  if (names(beta[, , 1])[1] == "(Intercept)")
    x <- methods::cbind2(1, x)

  n <- NROW(x)
  p <- NROW(beta)
  m <- NCOL(beta)
  n_penalties <- dim(beta)[3]

  stopifnot(p == NCOL(x))

  lin_pred <- array(dim = c(n, m, n_penalties),
                    dimnames = list(rownames(x),
                                    dimnames(beta)[[2]],
                                    dimnames(beta)[[3]]))

  for (i in seq_len(n_penalties))
    lin_pred[, , i] <- as.matrix(x %*% beta[, , i])

  lin_pred
}

#' @rdname predict.Owl
#' @export
predict.OwlGaussian <- function(object,
                                x,
                                lambda = NULL,
                                sigma = NULL,
                                type = c("link", "response"),
                                exact = FALSE,
                                simplify = TRUE,
                                ...) {
  type <- match.arg(type)

  out <- NextMethod(object, type = type) # always linear predictors

  if (simplify)
    out <- drop(out)

  out
}

#' @rdname predict.Owl
#' @export
predict.OwlBinomial <- function(object,
                                x,
                                lambda = NULL,
                                sigma = NULL,
                                type = c("link", "response", "class"),
                                exact = FALSE,
                                simplify = TRUE,
                                ...) {

  type <- match.arg(type)

  lin_pred <- NextMethod(object, type = type)

  out <- switch(
    type,
    link = lin_pred,
    response = 1 / (1 + exp(-lin_pred)),
    class = {
      cnum <- ifelse(lin_pred > 0, 2, 1)
      clet <- object$class_names[cnum]

      if (is.matrix(cnum))
        clet <- array(clet, dim(cnum), dimnames(cnum))

      clet
    }
  )

  if (simplify)
    out <- drop(out)

  out
}

#' @rdname predict.Owl
#' @export
predict.OwlPoisson <- function(object,
                               x,
                               sigma = NULL,
                               type = c("link", "response"),
                               exact = FALSE,
                               simplify = TRUE,
                               ...) {

  type <- match.arg(type)

  lin_pred <- NextMethod(object, type = type)

  out <- switch(
    type,
    link = lin_pred,
    response = exp(lin_pred)
  )

  if (simplify)
    out <- drop(out)

  out
}



