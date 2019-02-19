#' Golem
#'
#' @param x asdf
#' @param y
#' @param family
#' @param penalty
#' @param solver
#' @param intercept
#' @param lambda
#' @param sigma
#' @param fdr
#' @param standardize
#' @param max_passes
#' @param tol
#' @param ...
#'
#' @return wef
#' @export
#'
#' @examples
#' fit(3)
golem <- function(x,
                  y,
                  family = c("gaussian",
                             "binomial"),
                  penalty = "SLOPE",
                  solver = "auto",
                  intercept = TRUE,
                  lambda = "gaussian",
                  sigma = NULL,
                  fdr = 0.2,
                  standardize = c("features", "response", "both", "none"),
                  max_passes = 1000,
                  tol = 1e-6,
                  ...) {
  # collect the call so we can use it in update() later on
  ocall <- match.call()

  # collect settings
  if (isFALSE(standardize)) {
    standardize <- "none"
  } else if (isTRUE(standardize)) {
    standardize <- "features"
  } else {
    standardize <- match.arg(standardize)
  }

  family <- match.arg(family)
  fit_intercept <- intercept
  lambda_type <- lambda

  # collect options
  debug <- isTRUE(getOption("golem.debug"))

  n_samples <- NROW(x)
  n_features <- NCOL(x)
  n_targets <- NCOL(y)

  n <- NROW(x)
  p <- NCOL(x)
  m <- NCOL(y)

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

  # lambda (regularization strength)
  if (is.character(lambda)) {
    lambda_type <- lambda
    lambda <- rep.int(0, p)
    stopifnot(lambda_type %in% c("gaussian", "bhq"))
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)
    stopifnot(length(lambda) == p)

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")
  }

  if (is.null(sigma)) {
    sigma <- 1 # temporarily set sigma to 1
    if (n >= p + 30) {
      sigma_type <- "residuals"
    } else {
      sigma_type <- "iterative"
    }
  } else {
    sigma_type <- "user"
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

      # Transform response to {0, 1}, which is used internally
      y <- as.numeric(as.factor(y)) - 1
    }
  )

  y <- as.matrix(y)

  control <- list(debug = debug,
                  sigma = sigma,
                  sigma_type = sigma_type,
                  family_choice = family,
                  penalty_choice = penalty,
                  solver_choice = solver,
                  fit_intercept = intercept,
                  is_sparse = is_sparse,
                  lambda = lambda,
                  lambda_type = lambda_type,
                  fdr = fdr,
                  max_passes = max_passes,
                  standardize = standardize,
                  tol = tol)

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

  if (family %in% c("gaussian", "binomial")) {
    beta <- as.double(res$beta)

    if (fit_intercept) {
      coefficients <- c(res$intercept, beta)
      names(coefficients) <- c("(Intercept)", variable_names)
    } else {
      coefficients <- beta
      names(coefficients) <- variable_names
    }

    df <- sum(beta > 0)
  }

  out <- structure(list(coefficients = coefficients,
                        lambda = drop(res$lambda),
                        df = df,
                        passes = res$passes,
                        class_names = class_names,
                        sigma = res$sigma,
                        n = n,
                        p = p,
                        call = ocall),
                   class = c(paste0("Golem", firstUpper(family)),
                             "Golem"))

  # if (debug)
  #   attr(out, "diagnostics") <- list(loss = res$losses)
  out
}
