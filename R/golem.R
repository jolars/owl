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
#' * [golem::Slope()]
#' * [golem::Lasso()]
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
#' * [golem::Fista()]
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
#'
#' @return The result of fitting
#' @export
#'
#' @examples
#' X <- with(mtcars, cbind(cyl, wt, disp, hp, drat))
#' y <- mtcars$mpg
#'
#' golem_fit <- golem::golem(X, y, family = "gaussian")
golem <- function(x,
                  y,
                  family = c("gaussian", "binomial"),
                  penalty = golem::Slope(),
                  solver = golem::Fista(),
                  intercept = TRUE,
                  standardize = c("features", "response", "both", "none"),
                  ...) {

  stopifnot(is.logical(intercept))

  # collect the call so we can use it in update() later on
  ocall <- match.call()

  fit_intercept <- intercept
  n_samples  <- n <- NROW(x)
  n_features <- p <- NCOL(x)
  n_targets  <- m <- NCOL(y)

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

  # setup family
  family <- switch(match.arg(family),
                   gaussian = Gaussian(),
                   binomial = Binomial())

  y <- preprocessResponse(family, y)
  x <- preprocessFeatures(x, standardize)

  y_center <- attr(y, "center")
  y_scale  <- attr(y, "scale")
  x_center <- attr(x, "center")
  x_scale  <- attr(x, "scale")

  class_names <- attr(y, "class_names")

  # setup penalty settings
  if (is.character(penalty))
    penalty <- match.fun(penalty)()
  else if (is.function(penalty))
    penalty <- do.call(penalty, list())

  penalty <- setup(penalty, family, x, y, y_scale)

  # setup solver settings
  if (is.character(solver))
    solver <- match.fun(solver)()
  else if (is.function(solver))
    solver <- do.call(solver, list())

  stopifnot(isClass(solver, "Solver"),
            isClass(penalty, "Penalty"),
            isClass(family, "Family"))

  # collect response and variable names (if they are given) and otherwise
  # make new
  response_names <- colnames(y)
  variable_names <- colnames(x)

  if (is.null(variable_names))
    variable_names <- paste0("V", seq_len(p))
  if (is.null(response_names))
    response_names <- paste0("y", seq_len(m))

  lipschitz_constant <- lipschitzConstant(family, x, fit_intercept)

  control <- list(debug = debug,
                  family = family,
                  penalty = penalty,
                  solver = solver,
                  fit_intercept = fit_intercept,
                  is_sparse = is_sparse,
                  standardize = standardize,
                  x_center = x_center,
                  x_scale = x_scale,
                  y_center = y_center,
                  y_scale = y_scale,
                  lipschitz_constant = lipschitz_constant)

  # run the solver, perhaps iteratively
  # if (is.null(sigma)) {
  #   # run Algorithm 5 of Section 3.2.3.
  #   selected <- NULL
  #   repeat {
  #     selected.prev <- selected
  #     sigma <- c(sigma, estimate_noise(X[, selected, drop = FALSE], y))
  #     result <- SLOPE_solver_call(X, y, tail(sigma,1) * lambda)
  #     selected <- result$selected
  #     if (identical(selected, selected.prev))
  #       break
  #     if (length(selected) + 1 >= n)
  #       stop("selected >= n-1 variables. Cannot estimate variance.")
  #   }
  # } else {
  #   result = SLOPE_solver_call(X, y, sigma * lambda)
  # }

  if (is_sparse) {
    stop("sparse feature matrices are not yet supported.")
  } else {
    res <- golemDense(x, y, control)
  }

  unstandardized_coefs <- unstandardize(res$intercept,
                                        res$beta,
                                        x_center,
                                        x_scale,
                                        y_center,
                                        y_scale,
                                        fit_intercept)

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


  nonzeros <- apply(beta, 3, function(x) colSums(x > 0))

  diagnostics <- solver@diagnostics

  if (diagnostics) {
    nl <- length(res$time)
    nn <- lengths(res$time)
    time <- unlist(res$time)
    primal <- unlist(res$primals)
    dual <- unlist(res$duals)

    diag <- data.frame(time = time,
                       primal = primal,
                       dual = dual,
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
      passes = as.double(res$passes),
      call = ocall)
}

