#' Print results from owl fit
#'
#' @param x an object of class `'Owl'` or `'TrainedOwl'`
#' @param digits number of significant digits to print in results
#' @param ... other arguments passed to [print()]
#'
#' @return Prints output on the screen
#' @export
print.Owl <- function(x, digits = 2, ...) {
  sigma <- x$sigma
  n_nonzero <- apply(x$nonzeros, 3, sum)
  deviance_ratio <- x$deviance_ratio

  out <- data.frame(sigma = sigma,
                    deviance_ratio = deviance_ratio,
                    n_nonzero = n_nonzero)

  # print call
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  # print path summary
  cat("\nPath summary:\n")
  print(out, digits = digits, ...)

  # print lambda path
  cat("\nLambda sequence:\n")
  print(as.vector(x$lambda), digits = digits, ...)
}

#' @rdname print.Owl
#' @export
print.TrainedOwl <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Optimum values:\n")

  print(x$optima, ...)
}
