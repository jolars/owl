#' Title
#'
#' @param object
#' @param newdata
#' @param type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.Golem <- function(object,
                          newdata = NULL,
                          type = c("link", "response", "terms"),
                          ...) {

  beta <- object$beta
  intercept <- object$intercept
  if (is.null(intercept))
    intercept <- 0

  if (is.null(newdata)) {
    if (type %in% c("link", "response", "class"))
      stop("you need to supply a value for 'newdata' for type = '", type, "'")
  } else {
    if (inherits(newdata, "sparseMatrix"))
      newdata <- methods::as(newdata, "dgCMatrix")
    if (inherits(newdata, "data.frame"))
      newdata <- as.matrix(newdata)

    lin_pred <- newdata %*% beta + intercept
  }

  switch(type,
         link = lin_pred,
         response = lin_pred,
         terms = rbind(intercept, beta),
         fit)
}

#' @inherit predict.Golem
#'
#' @export
#' @rdname predict.Golem
predict.GolemGaussian <- function(object,
                                  newdata = NULL,
                                  type = c("link",
                                           "response",
                                           "terms"),
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
                                           "terms",
                                           "class"),
                                  ...) {
  type <- match.arg(type)
  fit <- NextMethod("predict", type = type)
  switch(
    type,
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
