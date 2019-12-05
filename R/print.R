print.Owl <- function(x, ...) {

}

print.TrainedOwl <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Optimum values:\n")

  print(x$optima)
}
