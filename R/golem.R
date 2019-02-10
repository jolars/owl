#' Title
#'
#' @param x asdf
#'
#' @return wef
#' @export
#'
#' @examples
#' fit(3)
fit <- function(X,
                y,
                fdr = 0.20,
                lambda = "gaussian",
                sigma = NULL,
                normalize = TRUE,
                ...) {
  # Validate input types.
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X))
    X.names = colnames(X)
  else
    stop("Input X must be a matrix or data frame")
  n <- nrow(X)
  p <- ncol(X)
  y <- as.numeric(y)

  # Normalize input, if necessary.
  if (normalize) {
    X = normalize(X)
    y = y - mean(y)
  }

  # Create lambda sequence.
  if (is.character(lambda)) {
    lambda_method = lambda
    lambda = create_lambda(n, p, fdr, lambda)
  } else {
    lambda_method = 'user'
    lambda = as.numeric(lambda)
  }

  # Validate input constraints.
  stopifnot(length(y) == n, length(lambda) == p)
  if (is.unsorted(rev(lambda)))
    stop('Lambda sequence must be non-increasing');

  # Estimate the noise level, if possible.
  if (is.null(sigma) && n >= p + 30)
    sigma = estimate_noise(X, y)

  # Run the solver, iteratively if necessary.
  if (is.null(sigma)) {
    # Run Algorithm 5 of Section 3.2.3.
    selected = NULL
    repeat {
      selected.prev = selected
      sigma = c(sigma, estimate_noise(X[,selected,drop=F], y))
      result = SLOPE_solver_call(X, y, tail(sigma,1) * lambda)
      selected = result$selected
      if (identical(selected, selected.prev))
        break
      if (length(selected)+1 >= n)
        stop('Selected >= n-1 variables. Cannot estimate variance.')
    }
  } else {
    result = SLOPE_solver_call(X, y, sigma * lambda)
  }

  # Package up the results.
  beta = result$beta
  selected = result$selected
  if (!is.null(X.names))
    names(selected) = X.names[selected]

  structure(list(call = match.call(),
                 lambda = lambda,
                 lambda_method = lambda_method,
                 sigma = sigma,
                 beta = beta,
                 selected = selected),
            class = 'SLOPE.result')
}


# Helper function for invoking the SLOPE solver.
SLOPE_solver_call <- function(X, y, lambda) {

  result = slope_solver(X, y, lambda)
  beta = result$x
  tol = 0 # Our solver sets un-selected beta's to *exactly* zero.

  selected = which(abs(beta) > tol)
  list(beta = beta, selected = selected)
}

# Print method for SLOPE results.
#' @export
#' @keywords internal
print.SLOPE.result <- function(x, ...) {
  result = x
  cat("\nCall:\n", paste(deparse(result$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  selected = result$selected
  if (length(selected)) {
    beta.selected <- result$beta[selected]
    if (is.null(names(selected)))
      names(beta.selected) <- as.character(selected)
    else
      names(beta.selected) <- names(selected)
    cat("Selected", length(beta.selected), "variables with coefficients:\n")
    print.default(beta.selected, print.gap = 2L, quote = FALSE)
  } else {
    cat("No selected variables\n")
  }
  cat("\n")
  invisible(result)
}
