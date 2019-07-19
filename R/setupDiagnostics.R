#' Setup a data.frame of diagnostics
#'
#' @param res the result from calling the C++ routine used to fit a model
#'   in golem
#'
#' @return A data.frame
#'
#' @keywords internal
setupDiagnostics <- function(res) {
  time <- res$time
  primals <- res$primals
  duals <- res$duals
  infeasibilities <- res$infeasibilities

  nl <- length(time)
  nn <- lengths(time)
  time <- unlist(time)
  primal <- unlist(primals)
  dual <- unlist(duals)
  infeasibility <- unlist(infeasibilities)

  data.frame(iteration = seq_along(time),
             time = time,
             primal = primal,
             dual = dual,
             infeasibility = infeasibility,
             penalty = rep(seq_len(nl), nn))
}

