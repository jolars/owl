Family <- R6::R6Class(
  "Family",
  public = list(
    name = NULL,

    predictClass = function(lin_pred, class_names) {
      stop("class predictions are not available for this family")
    },

    link = function(lin_pred) {
      lin_pred
    }
  )
)

Gaussian <- R6::R6Class(
  "Gaussian",
  inherit = Family,
  public = list(
    name = "gaussian",

    lambdaMax = function(x, y, y_scale) {
      y_scale * max(abs(crossprod(x, y))) / NROW(x)
    },

    preprocessResponse = function(y) {
      y <- as.numeric(y)

      if (NCOL(y) > 1)
        stop("response for Gaussian regression must be one-dimensional.")

      y_center <- mean(y)
      y_scale  <- 1

      y <- as.matrix(y - y_center)

      list(y, y_center, y_scale, n_classes = 1L, class_names = NA_character_)
    },

    lipschitzConstant =
      function(x,
               fit_intercept,
               x_center,
               x_scale,
               standardize_features) {
      maxSquaredRowNorm(x, x_center/x_scale, standardize_features) +
          fit_intercept
    }
  )
)

Binomial <- R6::R6Class(
  "Binomial",
  inherit = Family,
  public = list(
    name = "binomial",

    lambdaMax = function(x, y, y_scale) {
      # convert y from {-1, 1} to {0, 1}
      y <- (y + 1) / 2

      # standardize
      y_bar <- mean(y)
      y_sd <- stats::sd(y)

      y <- (y - y_bar) / y_sd

      y_sd * max(abs(crossprod(x, y))) / NROW(x)
    },

    preprocessResponse = function(y) {
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
    },

    link = function(lin_pred) {
      1 / (1 + exp(-lin_pred))
    },

    predictClass = function(lin_pred, class_names) {

      cnum <- ifelse(lin_pred > 0, 2, 1)
      clet <- class_names[cnum]

      if (is.matrix(cnum))
        clet <- array(clet, dim(cnum), dimnames(cnum))

      clet
    },

    lipschitzConstant =
      function(x,
               fit_intercept,
               x_center,
               x_scale,
               standardize_features) {
        0.25 * (maxSquaredRowNorm(x, x_center/x_scale, standardize_features) +
          fit_intercept)
      }
  )
)
