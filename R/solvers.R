setClass("Solver",
         slots = c(name = "character",
                   tol_rel_gap = "numeric",
                   tol_infeas = "numeric",
                   max_passes = "numeric",
                   diagnostics = "logical"))

setClass("Fista",
         contains = "Solver")

Fista <- function(tol_rel_gap = 1e-6,
                  tol_infeas = 1e-6,
                  max_passes = 1e4,
                  diagnostics = FALSE) {

  stopifnot(tol_rel_gap >= 0,
            tol_infeas >= 0,
            max_passes > 0)

  new("Fista",
      name = "fista",
      tol_rel_gap = tol_rel_gap,
      tol_infeas = tol_infeas,
      max_passes = max_passes,
      diagnostics = diagnostics)
}

