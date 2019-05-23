setClass("Family",
         slots = c(name = "character"),
         prototype = list(name = NA_character_))

setClass("Gaussian", contains = "Family")
setClass("Binomial", contains = "Family")

setGeneric("lambdaMax",
           function(object, x, y, y_scale) standardGeneric("lambdaMax"))

setMethod("lambdaMax",
          "Gaussian",
          function(object, x, y, y_scale)
            y_scale*max(abs(crossprod(x, y)))/NROW(x))

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

setGeneric("lipschitzConstant",
           function(object, x, fit_intercept)
             standardGeneric("lipschitzConstant"))

setMethod("lipschitzConstant",
          "Gaussian",
          function(object, x, fit_intercept) {
            max(rowSums(x^2)) + fit_intercept
          })

setMethod("lipschitzConstant",
          "Binomial",
          function(object, x, fit_intercept) {
            0.25 * (max(rowSums(x^2)) + fit_intercept)
          })

setGeneric("preprocessResponse",
           function(object, y) standardGeneric("preprocessResponse"))

setMethod("preprocessResponse",
          "Gaussian",
          function(object, y) {
            m <- NCOL(y)

            if (!is.numeric(y))
              stop("non-numeric response.")

            if (m > 1)
              stop("response for Gaussian regression must be one-dimensional.")

            y <- as.numeric(y)

            y_center <- mean(y)
            #y_scale  <- sd(y)
            y_scale  <- 1

            y <- as.matrix(y - y_center)
            attr(y, "center") <- y_center
            attr(y, "scale") <- y_scale
            attr(y, "n_classes") <- 1
            attr(y, "class_names") <- NA_character_
            y
          })

setMethod(
  "preprocessResponse",
  "Binomial",
  function(object, y) {
    # m <- NCOL(y)

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

    attr(y, "scale") <- 1
    attr(y, "center") <- 0
    attr(y, "n_classes") <- 1
    attr(y, "class_names") <- class_names
    y
  }
)

setGeneric("link",
           function(object, lin_pred) standardGeneric("link"))

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
