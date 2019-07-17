#' Plot coefficients for a model fit by golem
#'
#' Plot the model's coefficient along the regularization path.
#'
#' @param x an object of class `"Golem"`
#' @param ... parameters that will be used to modify the call to
#'   [lattice::xyplot()]
#'
#' @seealso [lattice::xyplot()], [golem()], [plotDiagnostics()]
#'
#' @return An object of class `"trellis"`, which will be plotted on the
#'   current device unless stored in a variable.
#' @export
#'
#' @examples
#' fit <- golem(heart$x, heart$y, penalty = "lasso")
#' plot(fit)
plot.Golem = function(x, ...) {
  object <- x

  coefs <- object$coefficients

  if (is.null(coefs))
    stop("nothing to plot since model is yet to be fit")

  p <- NROW(coefs) # number of features
  m <- NCOL(coefs) # number of responses

  is_slope <- inherits(object$penalty, c("Slope", "GroupSlope"))

  args <- list()

  if (is_slope) {
    x <- object$penalty$sigma
    xlab <- expression(sigma)
  } else {
    x <- object$penalty$lambda
    xlab <- expression(lambda)
  }

  n_x <- length(x)
  d <- as.data.frame(as.table(coefs))
  d$x <- rep(x, each = p*m)

  args <- list(
    x = if (m > 1)
      quote(Freq ~ x | Var2)
    else
      quote(Freq ~ x),
    type = if (n_x == 1) "p" else "l",
    groups = quote(Var1),
    data = quote(d),
    ylab = expression(hat(beta)),
    xlab = xlab,
    # scales = list(x = list(log = "e")),
    # xscale.components = function(lim, ...) {
    #   x <- lattice::xscale.components.default(lim, ...)
    #   x$bottom$labels$labels <- parse(text = x$bottom$labels$labels)
    #   x
    # },
    auto.key = if (p <= 10)
      list(space = "right", lines = TRUE, points = FALSE)
    else FALSE,
    abline = within(lattice::trellis.par.get("reference.line"), {h = 0})
  )

  # switch(match.arg(xvar),
  #        norm = {
  #          plot_args$xlab <-
  #            expression(group("|", group("|", hat(beta), "|"), "|")[1])
  #          plot_data$xval <- if (is.list(beta))
  #            rowSums(vapply(beta,
  #                           function(x) colSums(abs(as.matrix(x))),
  #                           double(ncol(beta[[1]]))))
  #          else
  #            colSums(abs(as.matrix(beta)))
  #        },
  #        lambda = {
  #          plot_args$xlab <- expression(lambda)
  #          plot_args$scales <- list(x = list(log = "e"))
  #          plot_data$xval <- x$lambda
  #
  #          # Prettier x scale
  #          plot_args$xscale.components <- function(lim, ...) {
  #            x <- lattice::xscale.components.default(lim, ...)
  #            x$bottom$labels$labels <- parse(text = x$bottom$labels$labels)
  #            x
  #          }
  #        },
  #        dev = {
  #          plot_args$xlab <- "Fraction of deviance explained"
  #          plot_data$xval <- x$dev.ratio
  #        })

  # Let the user modify the plot parameters
  do.call(lattice::xyplot, utils::modifyList(args, list(...)))
}
