#' Make Predictions from Golem Models
#'
#' @param object an object of class "Golem"
#' @param newdata new data to make predictions for
#' @param type type of predictions. `"link"` returns the linear predictors,
#'   `"response"` returns the response after applying the link functions,
#'   and `"class"` returns the predicted class (for the binomial family).
#' @param ... ignored.
#'
#' @return A vector of predictions.
#' @export
#'
#' @examples
#' x <- as.matrix(with(Puromycin, cbind(conc, rate)))
#' y <- Puromycin$state
#'
#' fit <- golem(x, y, family = "binomial")
#' predict(fit, x, "class")
predict.Golem <- function(object,
                          newdata,
                          type = c("link", "response"),
                          ...) {
  if (is.null(newdata)) {
    if (type %in% c("link", "response", "class"))
      stop("you need to supply a value for 'newdata' for type = '", type, "'")
  }

  beta <- stats::coef(object)

  if (inherits(newdata, "sparseMatrix"))
    newdata <- methods::as(newdata, "dgCMatrix")
  if (inherits(newdata, "data.frame"))
    newdata <- as.matrix(newdata)

  if (names(beta)[1] == "(Intercept)") {
    newdata <- methods::cbind2(1, newdata)
  }

  lin_pred <- newdata %*% beta

  switch(type,
         link = lin_pred,
         response = lin_pred,
         lin_pred)
}

#' @inherit predict.Golem
#'
#' @export
#' @rdname predict.Golem
predict.GolemGaussian <- function(object,
                                  newdata = NULL,
                                  type = c("link",
                                           "response"),
                                  ...) {
  type <- match.arg(type)
  NextMethod("predict", type = type)
}

#' @inherit predict.Golem
#'
#' @export
#' @rdname predict.Golem
predict.GolemBinomial <- function(object,
                                  newdata = NULL,
                                  type = c("link",
                                           "response",
                                           "class"),
                                  ...) {
  type <- match.arg(type)
  fit <- NextMethod("predict", type = type)
  switch(
    type,
    link = fit,
    response = 1 / (1 + exp(-fit)),
    class = {
      cnum <- ifelse(fit > 0, 2, 1)
      clet <- object$class_names[cnum]
      if (is.matrix(cnum))
        clet <- array(clet, dim(cnum), dimnames(cnum))
      clet
    },
    fit
  )
}
