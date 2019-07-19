Fista <- function(tol_rel_gap = 1e-6,
                  tol_infeas = 1e-6,
                  max_passes = 1e4,
                  diagnostics = FALSE) {

  stopifnot(tol_rel_gap >= 0,
            tol_infeas >= 0,
            max_passes > 0,
            is.logical(diagnostics),
            is.finite(tol_rel_gap),
            is.finite(tol_infeas),
            is.finite(max_passes))

  structure(list(name = "fista",
                 tol_rel_gap = tol_rel_gap,
                 tol_infeas = tol_infeas,
                 max_passes = max_passes,
                 diagnostics = diagnostics),
            class = c("Fista", "Solver"))
}
