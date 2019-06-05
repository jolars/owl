setClass("Family",
         slots = c(name = "character"),
         prototype = list(name = NA_character_))

setClass("Gaussian", contains = "Family")
setClass("Binomial", contains = "Family")

setGeneric(
  "lambdaMax",
  function(object, x, y, y_scale)
    standardGeneric("lambdaMax")
)

setMethod(
  "lambdaMax",
  "Gaussian",
  function(object, x, y, y_scale)
    y_scale*max(abs(crossprod(x, y)))/NROW(x)
)

setMethod(
  "lambdaMax",
  "Binomial",
  function(object, x, y, y_scale) {
    # convert y from {-1, 1} to {0, 1}
    y <- (y+1) / 2

    # standardize
    y_bar <- mean(y)
    y_sd <- stats::sd(y)

    y <- (y-y_bar) / y_sd

    y_sd*max(abs(crossprod(x, y))) / NROW(x)
  }
)

setGeneric("preprocessResponse",
           function(object, y) standardGeneric("preprocessResponse"))

setMethod(
  "preprocessResponse",
  "Gaussian",
  function(object, y) {

    y <- as.numeric(y)

    if (NCOL(y) > 1)
      stop("response for Gaussian regression must be one-dimensional.")

    y_center <- mean(y)
    y_scale  <- 1

    y <- as.matrix(y - y_center)

    list(y, y_center, y_scale, n_classes = 1L, class_names = NA_character_)
  }
)

setMethod(
  "preprocessResponse",
  "Binomial",
  function(object, y) {
    if (NCOL(y) > 1)
      stop("response for binomial regression must be one-dimensional.")

    if (length(unique(y)) > 2)
      stop("more than two classes in response")

    if (length(unique(y)) == 1)
      stop("only one class in response.")

    y_table <- table(y)
    min_class <- min(y_table)

    if (min_class <= 1)
      stop("one class only has ", min_class, " observations.")

    class_names <- names(y_table)

    # Transform response to {-1, 1}, which is used internally
    y <- as.matrix(ifelse(as.numeric(as.factor(y)) == 1, -1, 1))

    list(y,
         y_center = 0,
         y_scale = 1,
         n_classes = 1L,
         class_names = class_names)
  }
)

setGeneric(
  "link",
  function(object, lin_pred)
    standardGeneric("link")
)

setMethod("link",
          "Family",
          function(object, lin_pred) {
            lin_pred
          })

setMethod("link",
          "Binomial",
          function(object, lin_pred) {
            1 / (1 + exp(-lin_pred))
          })

setGeneric("predictClass",
           function(object, lin_pred, class_names)
             standardGeneric("predictClass"))

setMethod("predictClass",
          "Gaussian",
          function(object, lin_pred, class_names) {
            stop("class predictions do not make sense for Gaussian responses.")
          })

setMethod(
  "predictClass",
  "Binomial",
  function(object, lin_pred, class_names) {

    cnum <- ifelse(lin_pred > 0, 2, 1)
    clet <- class_names[cnum]

    if (is.matrix(cnum))
      clet <- array(clet, dim(cnum), dimnames(cnum))

    clet
  }
)

Gaussian <- function() {
  new("Gaussian",
      name = "gaussian")
}

Binomial <- function() {
  new("Binomial",
      name = "binomial")
}
