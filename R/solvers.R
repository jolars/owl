setClass("Solver",
         slots = c(name = "character",
                   tol = "numeric",
                   max_passes = "numeric",
                   diagnostics = "logical"),
         prototype = c(name = NA_character_,
                       tol = NA_real_,
                       max_passes = NA_real_,
                       diagnostics = FALSE))

setClass("Fista",
         contains = "Solver")

#' FISTA solver
#'
#' @param tol tolerance threshold for convergence
#' @param max_passes maximum number of passes for the solver
#' @param diagnostics whether to save diagnostics (wall-clock time and loss)
#'
#' @return A set of parameters for the FISTA solver.
#' @export
Fista <- function(tol = 1e-5, max_passes = 1e4, diagnostics = FALSE) {

  stopifnot(tol > 0,
            max_passes > 0)

  new("Fista",
      name = "fista",
      tol = tol,
      max_passes = max_passes,
      diagnostics = diagnostics)
}

