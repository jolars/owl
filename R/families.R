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

    lipschitzConstant = function(x,
                                 fit_intercept,
                                 x_center,
                                 x_scale,
                                 standardize_features) {
      maxSquaredRowNorm(x, x_center/x_scale, standardize_features) +
          fit_intercept
    },

    score = function(fit, x, y, measure = c("deviance", "mse", "mae")) {
      measure <- match.arg(measure)

      y <- as.vector(y)
      y_hat <- fit$predict(x)

      switch(measure,
             deviance = apply((y_hat - y)^2, 3, mean),
             mse = apply((y_hat - y)^2, 3, mean),
             mae = apply(abs(y_hat - y), 3, mean))
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
      },

    score = function(fit,
                     x,
                     y,
                     measure = c("deviance",
                                 "mse",
                                 "mae",
                                 "misclass_error",
                                 "auc")) {
      measure <- match.arg(measure)

      prob_min <- 1e-05
      prob_max <- 1 - prob_min

      y <- as.factor(y)
      y <- diag(2)[as.numeric(y), ]

      y_hat <- fit$predict(x, type = "response")

      switch(
        measure,
        auc = {
          apply(y_hat, 3, function(y_hat_i) auc(y, y_hat_i))
        },
        mse = apply((y_hat + y[, 1] - 1)^2 + (y_hat - y[, 2])^2, 3, mean),
        mae = apply(abs(y_hat + y[, 1] - 1) + abs(y_hat - y[, 2]), 3, mean),
        deviance = {
          y_hat <- pmin(pmax(y_hat, prob_min), prob_max)
          lp <- y[, 1] * log(1 - y_hat) + y[, 2] * log(y_hat)
          ly <- log(y)
          ly[y == 0] <- 0
          ly <- drop((y * ly) %*% c(1, 1))
          apply(2 * (ly - lp), 3, mean)
        },
        misclass_error =
          apply(y[, 1] * (y_hat > 0.5) + y[, 2] * (y_hat <= 0.5), 3, mean)
      )
    }
  )
)
