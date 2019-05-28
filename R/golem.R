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
#' @section Regularization Penalties:
#' There is a multitude of ways to penalize the models created by
#' [golem::golem()], currently they are:
#'
#' * SLOPE
#' * Group SLOPE
#' * LASSO
#'
#' These functions are in fact only parameter packs for the actual
#' implementations of the penalties and will be passed on to
#' the respective C++ creator functions, where the magic happens.
#'
#' Do *not* attempt to create your own penalty functions using this interface.
#' Such attempts will most likely be caught in assertions before anything bad
#' happens, but all bets are off if you are able to sneak them through
#' the various checks.
#'
#' @section Solvers:
#' There is currently a single solver available for [golem::golem], namely
#'
#' * FISTA
#'
#' @param x input matrix
#' @param y response variable
#' @param family response type, one of `'gaussian'`, `'binomial'`,
#'   `'multinomial'`, or `'mgaussian'`. See **Supported families** for details.
#' @param penalty the regularization penalty to use, either in the
#'   form of the output from one of this package's penalty functions,
#'   the function itself, or a character vector specifying one such function.
#'   Each function has its respective set of parameters, such as the
#'   regularization strength. Please see
#'   *Regularization Penalties* for more information.
#' @param solver the solver to use to optimize the loss function (objective).
#'   Just like the `penalty` parameter, this argument may be
#'   either a function, the function's output, or a character vector.
#'   Control arguments (such as convergence threshold) are set in the
#'   solver function itself. Please see **Solvers** for more information.
#' @param intercept whether to fit an intercept or not
#' @param standardize the type of standardization to carry out. Note that
#'   currently, standardization of response has no real effect. The
#'   response is always standardized for Gaussian responses and never
#'   for binomial
#' @param ... currently ignored
#' @param groups vector of integers to indicate group membership of each
#'   feature (only relevant for Group SLOPE)
#' @param sigma noise estimate (only relevant for SLOPE and Group SLOPE)
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or
#' @param fdr target false discovery rate (only relevant for SLOPE and
#'   Group SLOPE)
#' @param n_lambda length of regularization path (only relevant for lasso)
#' @param lambda_min_ratio smallest value for `lambda` as a fraction of
#'   \eqn{\lambda_\text{max}}{\lambda_max}#'
#' @param tol tolerance for optimizer
#' @param max_passes maximum number of passes for optimizer
#' @param diagnostics should diagnostics be saved for the model fit (timings,
#'   primal and dual objectives, and infeasibility)
#' @param orthogonalize whether `x` should be orthogonalized
#' @return An object of class `"Golem"`.
#' @export
#'
#' @examples
#' X <- with(mtcars, cbind(cyl, wt, disp, hp, drat))
#' y <- mtcars$mpg
#'
#' golem_fit <- golem::golem(X, y, family = "gaussian")
golem <- function(x,
                  y,
                  groups = NULL,
                  family = c("gaussian", "binomial"),
                  penalty = c("slope", "group_slope", "lasso"),
                  solver = "fista",
                  intercept = TRUE,
                  standardize = c("features", "response", "both", "none"),
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
            is.character(family),
            is.character(solver),
            is.character(penalty))

  # collect the call so we can use it in update() later on
  ocall <- match.call()

  family <- match.arg(family)
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)

  fit_intercept <- intercept
  n <- NROW(x)
  p <- NCOL(x)
  m <- NCOL(y)

  if (NROW(y) != NROW(x))
    stop("the number of samples in 'x' and 'y' must match")

  if (NROW(y) == 0)
    stop("the response (y) is empty")

  if (NROW(x) == 0)
    stop("the feature matrix (x) is empty")

  if (anyNA(y) || anyNA(x))
    stop("missing values are not allowed")

  # convert sparse x to dgCMatrix class from package Matrix.
  if (is_sparse <- inherits(x, "sparseMatrix")) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  # collect settings
  if (isFALSE(standardize)) {
    standardize <- "none"
  } else if (isTRUE(standardize)) {
    standardize <- "features"
  } else {
    standardize <- match.arg(standardize)
  }

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

  y <- preprocessResponse(family, y)
  res <- preprocessFeatures(penalty, x, standardize)
  x <- res$x
  penalty <- res$penalty

  y_center <- attr(y, "center")
  y_scale  <- attr(y, "scale")
  x_center <- attr(x, "center")
  x_scale  <- attr(x, "scale")

  class_names <- attr(y, "class_names")

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

  x <- sweep(x, 2, weights, "/")

  lipschitz_constant <- lipschitzConstant(family, penalty, x, fit_intercept)

  control <- list(family = family,
                  penalty = penalty,
                  solver = solver,
                  fit_intercept = fit_intercept,
                  is_sparse = is_sparse,
                  weights = weights,
                  standardize = standardize,
                  x_center = x_center,
                  x_scale = x_scale,
                  y_center = y_center,
                  y_scale = y_scale,
                  lipschitz_constant = lipschitz_constant)

  golemFit <- if (is_sparse) {
    stop("sparse feature matrices are not yet supported.")
  } else {
    golemDense
  }

  if (is_slope && is.null(sigma)) {
    if (inherits(family, "Gaussian")) {
      control$penalty@sigma <- penalty@sigma <- sd(y)
      res <- golemDense(x, y, control)

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
          new_x <- cbind(1, new_x)

        OLS <- lm.fit(new_x, y)
        if (standardize %in% c("features", "both")) {
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

  unstandardized_coefs <-
    postProcess(penalty, res$intercept, beta, x, y, fit_intercept)

  beta <- unstandardized_coefs$betas
  intercept <- unstandardized_coefs$intercepts

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

  nonzeros <- unstandardized_coefs$selected

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

