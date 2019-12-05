#' Setup a data.frame of diagnostics
#'
#' @param res the result from calling the C++ routine used to fit a model
#'   in owl
#'
#' @return A data.frame
#'
#' @keywords internal
setupDiagnostics <- function(res) {
  time <- res$time
  primals <- res$primals
  duals <- res$duals
  infeasibilities <- res$infeasibilities
  line_searches <- res$line_searches

  nl <- length(time)
  nn <- lengths(time)
  time <- unlist(time)
  primal <- unlist(primals)
  dual <- unlist(duals)
  infeasibility <- unlist(infeasibilities)
  line_searches <- unlist(line_searches)

  data.frame(iteration = unlist(lapply(nn, seq_len)),
             time = time,
             primal = primal,
             dual = dual,
             infeasibility = infeasibility,
             line_searches = line_searches,
             penalty = rep(seq_len(nl), nn))
}

