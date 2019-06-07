Solver <- R6::R6Class(
  "Solver",
  list(
    name = NULL,
    max_passes = 1e4,
    diagnostics = FALSE
  )
)

Fista <- R6::R6Class(
  "Fista",
  inherit = Solver,
  list(
    tol_rel_gap = NULL,
    tol_infeas = NULL,

    initialize = function(tol_rel_gap = 1e-6,
                          tol_infeas = 1e-6,
                          max_passes = 1e4,
                          diagnostics = FALSE) {

      self$name <- "fista"

      stopifnot(tol_rel_gap >= 0,
                tol_infeas >= 0,
                max_passes > 0,
                is.logical(diagnostics),
                is.finite(tol_rel_gap),
                is.finite(tol_infeas),
                is.finite(max_passes))

      self$tol_rel_gap <- tol_rel_gap
      self$tol_infeas <- tol_infeas
      self$max_passes <- max_passes
      self$diagnostics <- diagnostics
    }
  )
)

