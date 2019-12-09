#' FISTA (Fast Iterative Shrinking--Thresholding Algorithm)
#'
#' FISTA is a first-order gradient method.
#'
#' @param tol_rel_gap relative tolerance criterion for the duality gap
#' @param tol_infeas tolerance criterion for infeasibility of the
#'   dual objective
#'
#' @return An object of class `c("Fista", "Solver")`
#' @export
FISTA <- function(tol_rel_gap = 1e-6,
                  tol_infeas = 1e-6) {

  stopifnot(tol_rel_gap >= 0,
            tol_infeas >= 0,
            is.finite(tol_rel_gap),
            is.finite(tol_infeas))

  structure(list(name = "fista",
                 tol_rel_gap = tol_rel_gap,
                 tol_infeas = tol_infeas),
            class = c("FISTA", "Solver"))
}

#' ADMM (Alternating Direction Method of Multipliers)
#'
#' @param tol_rel relative tolerance criterion for convergence
#' @param tol_abs absolute tolerance criterion for convergence
#' @param alpha over-regularization paramter, must be in `(1.0, 1.8`)
#'
#' @return An object of clas `c("ADMM", "Solver"`
#' @export
ADMM <- function(tol_rel = 1e-5,
                 tol_abs = 1e-4,
                 alpha = 1.6) {

  stopifnot(tol_rel >= 0,
            tol_abs >= 0,
            alpha >= 1 && alpha <= 1.8,
            is.finite(tol_rel),
            is.finite(tol_abs))

  structure(list(name = "admm",
                 tol_rel = tol_rel,
                 tol_abs = tol_abs,
                 alpha = alpha),
            class = c("ADMM", "Solver"))
}
