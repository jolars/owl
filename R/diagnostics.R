#' Diagnostics
#'
#' This class encapsulates diagnostics from model fit with
#' [golem()].
#'
#' @docType class
#'
#' @format NULL
#' @keywords NULL
#'
#' @field data contains a data frame with variables
#'   \describe{
#'     \item{`time`}{time}
#'     \item{`primal`}{primal objective}
#'     \item{`dual`}{dual objective}
#'     \item{`infeasibility`}{infeasibility}
#'     \item{`penalty`}{the index along the regularization path}
#'   }
#'
#' @section Methods:
#'
#' \describe{
#'   \item{`plot()`}{
#'     Plot objectives and infeasibility side-by-side using
#'     trellis graphics.
#'   }
#' }
#'
#' @export
Diagnostics <- R6::R6Class(
  "Diagnostics",
  list (
    data = NULL,

    initialize = function(time, primals, duals, infeasibilities) {
      nl <- length(time)
      nn <- lengths(time)
      time <- unlist(time)
      primal <- unlist(primals)
      dual <- unlist(duals)
      infeasibility <- unlist(infeasibilities)

      self$data <- data.frame(time = time,
                              primal = primal,
                              dual = dual,
                              infeasibility = infeasibility,
                              penalty = rep(seq_len(nl), nn))
    },

    plot = function() {
      d <- self$data

      p1 <- lattice::xyplot(primal + dual ~ time,
                            data = d,
                            type = "l",
                            ylab = "Objetive",
                            xlab = "Time (Seconds)",
                            grid = TRUE,
                            auto.key = list(space = "inside",
                                            lines = TRUE,
                                            points = FALSE))

      p2 <- lattice::xyplot(infeasibility ~ time,
                            data = d,
                            type = "l",
                            grid = TRUE,
                            xlab = "Time (Seconds)",
                            ylab = "Infeasibility")

      gridExtra::grid.arrange(p1, p2, ncol = 2)
    }
  )
)
