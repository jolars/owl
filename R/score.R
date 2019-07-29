#' Return the loss of a golem model
#'
#' This function is a unified interface to return various types of loss for a
#' model fit with [golem()].
#'
#' @param object an object of class `"Golem"`
#' @param x feature matrix
#' @param y response
#' @param measure type of target measure. The default, `"deviance"`,
#'   returns mean squared error for Gaussian models and deviance for
#'   Binomial models. `"mse"` returns mean squared error. `"mae"` returns
#'   mean absolute error, `"accuracy"` returns classification rate accuracy,
#'   and `"auc"` returns area under the ROC curve.
#'
#' @return The measure along the regularization path depending on the
#'   value in `measure`.
#'
#' @export
#'
#' @examples
#' x <- subset(infert, select = c("induced", "age", "pooled.stratum"))
#' y <- infert$case
#'
#' fit <- golem(x, y, family = "binomial")
#' score(fit, x, y, measure = "auc")
score <- function(object, x, y, measure)
  UseMethod("score")

#' @rdname score
#' @export
score.GolemGaussian <- function(object,
                                x,
                                y,
                                measure = c("deviance", "mse", "mae")) {
  measure <- match.arg(measure)

  y <- as.vector(y)
  y_hat <- stats::predict(object, x, simplify = FALSE)

  switch(measure,
         deviance = apply((y_hat - y)^2, 3, mean),
         mse = apply((y_hat - y)^2, 3, mean),
         mae = apply(abs(y_hat - y), 3, mean))
}

#' @rdname score
#' @export
score.GolemBinomial <- function(object,
                                x,
                                y,
                                measure = c("deviance",
                                            "mse",
                                            "mae",
                                            "accuracy",
                                            "auc")) {
  measure <- match.arg(measure)

  prob_min <- 1e-05
  prob_max <- 1 - prob_min

  y <- as.factor(y)
  y <- diag(2)[as.numeric(y), ]

  y_hat <- stats::predict(object, x, type = "response", simplify = FALSE)

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
    accuracy = 1 -
      apply(y[, 1] * (y_hat > 0.5) + y[, 2] * (y_hat <= 0.5), 3, mean)
  )
}

auc <- function(y, prob, weights = rep.int(1, nrow(y))) {
  if (is.matrix(y) || is.data.frame(y)) {

    ny <- nrow(y)
    auc(rep(c(0, 1), c(ny, ny)), c(prob,prob), as.vector(weights*y))

  } else {

    if (is.null(weights)) {
      rprob <- rank(prob)
      n1 <- sum(y)
      n0 <- length(y) - n1
      u <- sum(rprob[y == 1]) - n1*(n1 + 1)/2
      exp(log(u) - log(n1) - log(n0))
    } else {
      # randomize ties
      rprob <- stats::runif(length(prob))
      op <- order(prob, rprob)
      y <- y[op]
      weights <- weights[op]
      cw <- cumsum(weights)
      w1 <- weights[y == 1]
      cw1 <- cumsum(w1)
      wauc <- log(sum(w1 * (cw[y == 1] - cw1)))
      sumw1 <- cw1[length(cw1)]
      sumw2 <- cw[length(cw)] - sumw1
      exp(wauc - log(sumw1) - log(sumw2))
    }
  }
}
