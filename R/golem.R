#' @include families.R penalties.R solvers.R
setClass("Golem",
         slots = c(coefficients = "array",
                   nonzeros = "numeric",
                   family = "Family",
                   penalty = "Penalty",
                   class_names = "character",
                   n = "numeric",
                   p = "numeric",
                   m = "numeric",
                   diagnostics = "data.frame",
                   passes = "numeric",
                   call = "call"))

#' Regularized Generalized Linear Models
#'
#' This functions fits a generalized linear model (GLM) using efficient
#' optimization routines suitable to big data problems.
#'
#' The objective for each model is simply the loss function for
#' each family plus a penalty term.
#'
#' @section Families:
#'
#' **Gaussian**
#'
#' The Gaussian model (Ordinary Least Squares) minimizes the following
#' objective.
#'
#' \deqn{
#'   ||\boldsymbol{y} - \boldsymbol{X\beta}||_2^2
#' }{
#'   ||y - X\beta||_2^2
#' }
#'
#' **Binomial**
#'
#' The binomial model (logistic regression) has the following objective.
#'
#' \deqn{
#'   \sum_{i=1}^n \log\left[1+ \exp\left(- y_i\boldsymbol{x}_i'\boldsymbol{\beta} \right) \right]
#' }{
#'   \sum log[1+ exp(- y_i x_i' \beta)]
#' }
#'
#' @section Penalties:
#' Models fit by [golem::golem()] can be regularized via several
#' penalties.
#'
#' **Lasso**
#'
#' The Lasso penalizes coefficients using the L1 norm. The penalty
#' term in the lagrangian form of the loss function is
#'
#' \deqn{
#'   \lambda \sum_{j=1}^p |\beta_j|
#' }{
#'   \lambda \sum |\beta|
#' }
#'
#' **SLOPE**
#'
#' SLOPE (Sorted L-One Penalized Estimation) is an extension of the Lasso.
#' Unlike the latter, however, SLOPE uses a non-increasing
#' sequence of \eqn{\lambda}---one
#' for each coefficient. The penalty term looks like
#'
#' \deqn{
#'   \sigma \sum_{i=j}^p \lambda_j |\beta|_{(j)}
#' }{
#'   \sigma \sum \lambda |\beta|(j)
#' }
#'
#' **Group SLOPE**
#'
#' Group SLOPE is an extension of Group LASSO. It applies the following
#' penalty
#'
#' \deqn{
#'   J_\lambda(\beta) = \sum_{j=1}^p \lambda_j |W||\beta||_{I,X}|_{(j)}
#' }{
#'   J_\lambda(\beta) = \sum \lambda_j |W||\beta||_{I,X}|_(j)
#' }
#'
#' @section Solvers:
#' There is currently a single solver available for [golem::golem].
#'
#' **FISTA**
#'
#' FISTA (Fast Iterative Shrinking-Tresholding Algorithm) is an extension
#' of the classical gradient algorithm.
#'
#' @param x feature matrix
#' @param y response
#' @param family response type. See **Families** for details.
#' @param penalty the regularization penalty to use. See **Penalties** for
#'   details.
#' @param solver the numerical solver to use. See **Solvers** for details.
#' @param intercept whether to fit an intercept
#' @param standardize_features whether to standardize features (predictors)
#' @param ... currently ignored
#' @param groups vector of integers to indicate group membership of each
#'   feature (only applies to Group SLOPE)
#' @param sigma noise estimate (only applies to SLOPE and Group SLOPE)
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or
#' @param fdr target false discovery rate (only applies to SLOPE and
#'   Group SLOPE)
#' @param n_lambda length of regularization path (only relevant for lasso)
#' @param lambda_min_ratio smallest value for `lambda` as a fraction of
#'   \eqn{\lambda_\mathrm{max}}{\lambda_max} (only applies to lasso)
#' @param tol tolerance for optimizer
#' @param max_passes maximum number of passes for optimizer
#' @param diagnostics should diagnostics be saved for the model fit (timings,
#'   primal and dual objectives, and infeasibility)
#' @param orthogonalize whether `x` should be orthogonalized. Note that
#'   setting this to TRUE when `x` is sparse will through an error.
#'   (only applies to Group SLOPE)
#' @return An object of class `"Golem"`.
#' @export
#'
#' @examples
#'
#' # Gaussian response, slope penalty (default)
#' gaussian_fit <- golem(abalone$x, abalone$y, family = "gaussian")
#'
#' # Binomial response, lasso penalty
#' binomial_fit <- golem(heart$x, heart$y, family = "binomial",
#'                       penalty = "lasso")
#'
golem <- function(x,
                  y,
                  groups = NULL,
                  family = c("gaussian", "binomial"),
                  penalty = c("slope", "group_slope", "lasso"),
                  solver = "fista",
                  intercept = TRUE,
                  standardize_features = TRUE,
                  orthogonalize = TRUE,
                  sigma = NULL,
                  lambda = NULL,
                  fdr = 0.2,
                  n_lambda = 100,
                  lambda_min_ratio = ifelse(NROW(x) < NCOL(x), 0.01, 0.0001),
                  tol = 1e-06,
                  max_passes = 1e4,
                  diagnostics = FALSE,
                  ...) {

  stopifnot(is.logical(intercept),
            is.logical(standardize_features),
            is.character(family),
            is.character(solver),
            is.character(penalty))

  # collect the call so we can use it in update() later on
  ocall <- match.call()

  family <- match.arg(family)
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)

  fit_intercept <- intercept

  if (NROW(y) != NROW(x))
    stop("the number of samples in 'x' and 'y' must match")

  if (NROW(y) == 0)
    stop("the response (y) is empty")

  if (NROW(x) == 0)
    stop("the feature matrix (x) is empty")

  if (anyNA(y) || anyNA(x))
    stop("missing values are not allowed")

  n <- NROW(x)
  p <- NCOL(x)
  m <- NCOL(y)

  # convert sparse x to dgCMatrix class from package Matrix.
  is_sparse <- inherits(x, "sparseMatrix")

  if (is_sparse) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  if (penalty == "group_slope" && standardize_features && is_sparse &&
      orthogonalize)
    stop("orthogonalization is currently not implemented for sparse data ",
         "when standardization is required")

  # setup penalty settings
  penalty <- switch(penalty,
                    slope = Slope(lambda = lambda,
                                  sigma = sigma,
                                  fdr = fdr),
                    group_slope = GroupSlope(groups = groups,
                                             lambda = lambda,
                                             sigma = sigma,
                                             fdr = fdr,
                                             orthogonalize = orthogonalize),
                    lasso = Lasso(lambda = lambda,
                                  lambda_min_ratio = lambda_min_ratio,
                                  n_lambda = n_lambda))

  # setup family
  family <- switch(family,
                   gaussian = Gaussian(),
                   binomial = Binomial())

  y_center <- y_scale <- double(m)
  n_classes <- integer(0)
  class_names <- character(0)

  c(y, y_center, y_scale, n_classes, class_names) %<-%
    preprocessResponse(family, y)

  x_center <- x_scale <- double(p)
  c(x, x_center, x_scale) %<-% standardize(x, standardize_features)

  c(penalty, x) %<-% preprocessFeatures(penalty,
                                        x,
                                        x_center,
                                        x_scale,
                                        standardize_features)

  penalty <- setup(penalty, family, x, y, y_scale)

  is_slope <- inherits(penalty, "Slope") || inherits(penalty, "GroupSlope")

  # setup solver settings
  solver <- Fista(tol = tol,
                  max_passes = max_passes,
                  diagnostics = diagnostics)

  # collect response and variable names (if they are given) and otherwise
  # make new
  response_names <- colnames(y)
  variable_names <- colnames(x)

  if (is.null(variable_names))
    variable_names <- paste0("V", seq_len(p))
  if (is.null(response_names))
    response_names <- paste0("y", seq_len(m))

  weights <- getWeights(penalty, x)

  for (j in seq_len(p))
    x[, j] <- x[, j]/weights[j]

  x_center <- x_center/weights

  lipschitz_constant <- lipschitzConstant(family,
                                          penalty,
                                          x,
                                          fit_intercept,
                                          x_center,
                                          x_scale,
                                          standardize_features)

  control <- list(family = family,
                  penalty = penalty,
                  solver = solver,
                  fit_intercept = fit_intercept,
                  is_sparse = is_sparse,
                  weights = weights,
                  standardize_features = standardize_features,
                  x_scaled_center = x_center/x_scale,
                  lipschitz_constant = lipschitz_constant)

  golemFit <- if (is_sparse) golemSparse else golemDense

  if (is_sparse)
    x <- methods::as(x, "dgCMatrix")

  if (is_slope && is.null(sigma)) {
    if (inherits(family, "Gaussian")) {
      control$penalty@sigma <- penalty@sigma <- stats::sd(y)
      res <- golemFit(x, y, control)

      S_new <- which(res$beta != 0)
      S <- c()

      while(!isTRUE(all.equal(S, S_new)) && (length(S_new) > 0)) {
        S <- S_new
        if (length(S) > n) {
          stop("sigma estimation fails because more predictors got ",
               "selected than there are observations.")
        }

        new_x <- x[, S, drop = FALSE]

        if (fit_intercept)
          new_x <- Matrix::cbind2(1, new_x)

        OLS <- stats::lm.fit(as.matrix(new_x), y)
        if (standardize_features) {
          sigma <- sqrt(sum(OLS$residuals^2) / (n - length(S) - 1))
        } else {
          sigma <- sqrt(sum(OLS$residuals^2) / (n - length(S)))
        }

        control$penalty@sigma <- penalty@sigma <- sigma

        res <- golemFit(x, y, control)
        S_new <- which(res$beta != 0)
      }
    } else if (inherits(family, "Binomial")) {
      control$penalty@sigma <- penalty@sigma <- 0.5
      res <- golemFit(x, y, control)
    }
  } else {
    res <- golemFit(x, y, control)
  }

  beta <- sweep(res$beta, 1, weights, "/")
  x_center <- x_center*weights

  nonzeros <- integer(0)

  c(intercept, beta, nonzeros) %<-% postProcess(penalty,
                                                res$intercept,
                                                beta,
                                                x,
                                                y,
                                                fit_intercept,
                                                x_center,
                                                x_scale,
                                                y_center,
                                                y_scale)

  n_penalties <- dim(beta)[3]

  if (fit_intercept) {
    coefficients <- array(NA, dim = c(p + 1, m, n_penalties))

    for (i in seq_len(n_penalties)) {
      coefficients[1, , i] <- intercept[, , i]
      coefficients[-1, , i] <- beta[, , i]
    }
    dimnames(coefficients) <- list(c("(Intercept)", variable_names),
                                   response_names,
                                   paste0("p", seq_len(n_penalties)))
  } else {
    coefficients <- beta
    dimnames(coefficients) <- list(variable_names,
                                   response_names,
                                   paste0("p", seq_len(n_penalties)))
  }

  diagnostics <- solver@diagnostics

  if (diagnostics) {
    nl <- length(res$time)
    nn <- lengths(res$time)
    time <- unlist(res$time)
    primal <- unlist(res$primals)
    dual <- unlist(res$duals)
    infeasibility <- unlist(res$infeasibilities)

    diag <- data.frame(time = time,
                       primal = primal,
                       dual = dual,
                       infeasibility = infeasibility,
                       penalty = rep(seq_len(nl), nn))
  } else {
    diag <- data.frame()
  }

  new("Golem",
      coefficients = coefficients,
      nonzeros = nonzeros,
      family = family,
      penalty = penalty,
      class_names = class_names,
      n = n,
      p = p,
      m = m,
      diagnostics = diag,
      passes = as.integer(res$passes),
      call = ocall)
}

