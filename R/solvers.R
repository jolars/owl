#' FISTA solver
#'
#' @param tol tolerance threshold for convergence
#' @param max_passes maximum number of passes for the solver
#'
#' @return A set of parameters for the FISTA solver.
#' @export
fista <- function(tol = 1e-3, max_passes = 1e4) {

  stopifnot(tol > 0,
            max_passes > 0)

  structure(list(name = "fista",
                 tol = tol,
                 max_passes = max_passes),
            class = "solver")
}
