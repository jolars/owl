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
#'   \item{`plot(ind = "last", what = c("objectives", "infeasibility"))`}{
#'     Plot objectives and infeasibility side-by-side using
#'     trellis graphics.
#'     \describe{
#'       \item{`ind`}{
#'         index of the fit for which diagnostics should be plotted
#'       }
#'       \item{`what`}{
#'         type of diagnostics to plot. `objectives` returns the
#'         primal and dual objectives whereas `infeasibility` returns
#'         the infeasibility metric. If both are provided, all plots
#'         will be arranged using [gridExtra::grid.arrange()]
#'       }
#'     }
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

    plot = function(ind = "last", what = c("objectives", "infeasibility")) {
      d <- self$data

      n_penalties <- length(unique(d$penalty))

      if (ind == "last") {
        ind <- unique(d$penalty)[n_penalties]
      } else {
        stopifnot(ind <= n_penalties,
                  ind >= 1)
      }

      d <- subset(d, subset = d$penalty == ind)

      p <- vector("list", length(what))
      i <- 1

      if ("objectives" %in% what) {
        p[[i]] <- lattice::xyplot(primal + dual ~ time,
                                  data = d,
                                  type = "l",
                                  ylab = "Objective",
                                  xlab = "Time (Seconds)",
                                  grid = TRUE,
                                  auto.key = list(space = "inside",
                                                  lines = TRUE,
                                                  points = FALSE))
        i <- i + 1
      }

      if ("infeasibility" %in% what) {
        p[[i]] <- lattice::xyplot(infeasibility ~ time,
                                  data = d,
                                  type = "l",
                                  grid = TRUE,
                                  xlab = "Time (Seconds)",
                                  ylab = "Infeasibility")
      }

      if (length(what) > 1)
        gridExtra::grid.arrange(grobs = p, ncol = length(p))
      else
        p[[1]]
    }
  )
)
