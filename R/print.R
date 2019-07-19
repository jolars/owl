print.Golem <- function(x, ...) {

}

print.TrainedGolem <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Optimum values:\n")

  print(x$optima)
}
