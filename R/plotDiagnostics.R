#' Plot results from diagnostics collected during model fitting
#'
#' This function plots various diagnostics collected during
#' the model fitting resulting from a call to [owl()] *provided that
#' `diagnostics = TRUE`*.
#'
#' @param object an object of class `"Owl"`.
#' @param ind either "last"
#' @param xvar what to place on the x axis. `iteration` plots each iteration, `time`
#'   plots the wall-clock time.
#' @param yvar what to place on the y axis. `objectives` returns the
#'   primal and dual objectives whereas `infeasibility` returns
#'   the infeasibility metric.
#' @param ... other arguments that will be used to modify the call to
#'   [lattice::xyplot()]
#'
#' @return An object of class `"trellis"`, which, unless stored in a variable,
#'   will be plotted when its default `print()` method is called.
#' @export
#'
#' @examples
#' x <- owl(abalone$x, abalone$y, sigma = 2, diagnostics = TRUE)
#' plotDiagnostics(x)
plotDiagnostics <- function(object,
                            ind = max(object$diagnostics$penalty),
                            xvar = c("time", "iteration"),
                            yvar = c("objectives", "infeasibility"),
                            ...) {

  stopifnot(inherits(object, "Owl"),
            is.numeric(ind),
            length(ind) == 1)

  xvar <- match.arg(xvar)
  yvar <- match.arg(yvar)

  if (is.null(object$diagnostics))
    stop("no diagnostics found in fit;",
         "did you call owl() with diagnostics = TRUE?")

  d <- object$diagnostics

  d <- subset(d, subset = d$penalty == ind)

  # setup a list of arguments to be provided to lattice::xyplot()
  args <- list(data = d, type = "l")

  if (nrow(d) > 1)
    args$grid <- TRUE

  if (yvar == "objectives") {
    args$x <- "primal + dual"
    args$ylab <- "Objective"
    args$auto.key <- list(space = "inside",
                          corner = c(0.95, 0.95),
                          lines = TRUE,
                          points = FALSE)

  } else if (yvar == "infeasibility") {
    args$x <- "infeasibility"
    args$ylab <- "Infeasibility"
  }

  if (xvar == "time") {
    args$x <- paste(args$x, "~ time")
    args$xlab <- "Time (seconds)"
  } else if (xvar == "iteration") {
    args$x <- paste(args$x, "~ iteration")
    args$xlab <- "Iteration"
  }

  args$x <- stats::as.formula(args$x)

  args <- utils::modifyList(args,
                            list(...))

  do.call(lattice::xyplot, args)
}
