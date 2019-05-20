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
setMethod(
  "predict",
  "Golem",
  function(object, newdata, type = c("link", "response", "class"), ...) {
    if (is.null(newdata)) {
      if (type %in% c("link", "response", "class"))
        stop("you need to supply a value for 'newdata' for type = '", type, "'")
    }

    family <- object@family

    beta <- coef(object)

    if (inherits(newdata, "sparseMatrix"))
      newdata <- methods::as(newdata, "dgCMatrix")

    if (inherits(newdata, "data.frame"))
      newdata <- as.matrix(newdata)

    if (names(beta)[1] == "(Intercept)") {
      newdata <- methods::cbind2(1, newdata)
    }

    lin_pred <- newdata %*% beta

    switch(
      match.arg(type),
      link = lin_pred,
      response = link(family, lin_pred),
      class = predictClass(family, lin_pred, object@class_names),
      lin_pred
    )
  }
)
