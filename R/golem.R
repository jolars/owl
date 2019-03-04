#' Regularized Generalized Linear Models
#'
#' This functions fits a generalized linear model (GLM) using efficient
#' optimization routines suitable to big data problems.
#'
#' @section Regularization Penalties:
#' There is a multitude of ways to penalize the models created by
#' [golem::golem()], currently they are:
#'
#' * [golem::slope()]
#' * [golem::elasticNet()]
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
#' * [golem::fista()]
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
#' @param diagnostics whether to collect diagnostics from the solver,
#'   which currently includes loss and timing estimates at each iteration
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
                  penalty = golem::slope(),
                  solver = golem::fista(),
                  intercept = TRUE,
                  standardize = c("features", "response", "both", "none"),
                  diagnostics = FALSE,
                  ...) {

  stopifnot(is.logical(diagnostics),
            is.logical(intercept))

  # collect the call so we can use it in update() later on
  ocall <- match.call()

  n <- NROW(x)
  p <- NCOL(x)
  m <- NCOL(y)

  # collect settings
  if (isFALSE(standardize)) {
    standardize <- "none"
  } else if (isTRUE(standardize)) {
    standardize <- "features"
  } else {
    standardize <- match.arg(standardize)
  }

  # setup penalty settings
  if (is.character(penalty))
    penalty_args <- match.fun(penalty)()
  else if (is.function(penalty))
    penalty_args <- do.call(penalty, list())
  else
    penalty_args <- penalty

  penalty_args <- utils::modifyList(penalty_args,
                                    list(n = n,
                                         p = p))

  # setup solver settings
  if (is.character(solver))
    solver_args <- match.fun(solver)()
  else if (is.function(solver))
    solver_args <- do.call(solver, list())
  else
    solver_args <- solver

  stopifnot(inherits(solver_args, "solver"),
            inherits(penalty_args, "penalty"))

  solver_args$diagnostics <- diagnostics

  family <- match.arg(family)
  fit_intercept <- intercept

  # collect options
  debug <- isTRUE(getOption("golem.debug"))

  n_samples <- NROW(x)
  n_features <- NCOL(x)
  n_targets <- NCOL(y)

  stopifnot(is.logical(intercept),
            is.logical(debug))

  if (NROW(y) != NROW(x))
    stop("the number of samples in 'x' and 'y' must match")

  if (NROW(y) == 0)
    stop("the response (y) is empty.")

  if (NROW(x) == 0)
    stop("the feature matrix (x) is empty.")

  # convert sparse x to dgCMatrix class from package Matrix.
  if (is_sparse <- inherits(x, "sparseMatrix")) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  # collect response and variable names (if they are given) and otherwise
  # make new
  response_names <- colnames(y)
  variable_names <- colnames(x)
  class_names    <- NULL

  if (is.null(variable_names))
    variable_names <- paste0("V", seq_len(p))
  if (is.null(response_names))
    response_names <- paste0("y", seq_len(m))

  if (any(is.na(y)) || any(is.na(x)))
    stop("NA values are not allowed.")

  switch(
    family,
    gaussian = {
      if (m > 1)
        stop("response for Gaussian regression must be one-dimensional.")

      if (!is.numeric(y))
        stop("non-numeric response.")

      n_classes <- 1L
      y <- as.numeric(y)
    },
    binomial = {
      if (length(unique(y)) > 2)
        stop("more than two classes in response")

      if (length(unique(y)) == 1)
        stop("only one class in response.")

      y_table <- table(y)
      min_class <- min(y_table)
      n_classes <- 1L

      if (min_class <= 1)
        stop("one class only has ", min_class, " observations.")

      class_names <- names(y_table)

      # Transform response to {-1, 1}, which is used internally
      y <- ifelse(as.numeric(as.factor(y)) == 1, -1, 1)
    }
  )

  y <- as.matrix(y)

  control <- list(debug = debug,
                  family_choice = family,
                  penalty_args = penalty_args,
                  solver_args = solver_args,
                  fit_intercept = intercept,
                  is_sparse = is_sparse,
                  standardize = standardize)

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

  beta <- res$beta
  intercept <- res$intercept
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

  out <- structure(list(coefficients = coefficients,
                        penalty = res$penalty,
                        nonzeros = nonzeros,
                        passes = res$passes,
                        class_names = class_names,
                        n = n,
                        p = p,
                        m = m,
                        call = ocall),
                   class = c(paste0("Golem", firstUpper(family)),
                             "Golem"))

  if (diagnostics) {
    nl <- length(res$time)
    nn <- lengths(res$time)
    time <- unlist(res$time)
    loss <- unlist(res$loss)

    diag <- data.frame(time = time,
                       loss = loss,
                       penalty = rep(seq_len(nl), nn))

    attr(out, "diagnostics") <- diag
  }

  out
}
