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
#'     \item{`iteration`}{iteration}
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
#'       \item{`x_var`}{
#'         what to place on the x axis. `iteration` plots each iteration, `time`
#'         plots the wall-clock time.
#'       }
#'       \item{`y_var`}{
#'         what to place on the y axis. `objectives` returns the
#'         primal and dual objectives whereas `infeasibility` returns
#'         the infeasibility metric.
#'       }
#'       \item{`\dots`}{
#'         other arguments that will be used to modify the call to
#'         [lattice::xyplot()]
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

      self$data <- data.frame(iteration = seq_along(time),
                              time = time,
                              primal = primal,
                              dual = dual,
                              infeasibility = infeasibility,
                              penalty = rep(seq_len(nl), nn))
    },

    plot = function(ind = "last",
                    x_var = c("time", "iteration"),
                    y_var = c("objectives", "infeasibility"),
                    ...) {
      d <- self$data

      n_penalties <- length(unique(d$penalty))

      x_var <- match.arg(x_var)
      y_var <- match.arg(y_var)

      if (ind == "last") {
        ind <- unique(d$penalty)[n_penalties]
      } else {
        stopifnot(ind <= n_penalties,
                  ind >= 1)
      }

      d <- subset(d, subset = d$penalty == ind)

      args <- list(data = d,
                   type = "l",
                   grid = TRUE)

      if (y_var == "objectives") {
        args$x <- "primal + dual"
        args$ylab <- "Objective"
        args$auto.key <- list(space = "right",
                              lines = TRUE,
                              points = FALSE)

      } else if (y_var == "infeasibility") {
        args$x <- "infeasibility"
        args$ylab <- "Infeasibility"
      }

      if (x_var == "time") {
        args$x <- paste(args$x, "~ time")
        args$xlab <- "Time (seconds)"
      } else if (x_var == "iteration") {
        args$x <- paste(args$x, "~ iteration")
        args$xlab <- "Iteration"
      }

      args$x <- as.formula(args$x)

      args <- utils::modifyList(args,
                                list(...))

      do.call(lattice::xyplot, args)
    }
  )
)
