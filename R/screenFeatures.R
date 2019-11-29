#' Screen features
#'
#' @param family family
#' @param penalty penalty
#' @param x features
#' @param y response
#' @param method type of screening
#' @param lambda regularization sequence of current fit
#' @param lambda_prev regularization sequence of previous fit
#' @param beta_prev coefficients from previous fit
#' @param intercept_prev intercept from previous fit
#' @param gradient_prev gradient from previous fit
#'
#' @return A logical vector indicating whether to drop a feature or not
#' @export
screenFeatures <- function(family,
                           penalty,
                           x,
                           y,
                           lambda,
                           lambda_prev,
                           beta_prev,
                           intercept_prev,
                           gradient_prev,
                           method = c("none", "strong", "safe")) {
  method <- match.arg(method)
  n <- NROW(x)
  p <- NCOL(x)

  switch(
    method,
    none = {
      rep(FALSE, p)
    },

    strong = {
      ord <- order(abs(gradient_prev), decreasing = TRUE)

      out <- logical(p)
      abs_gradient_prev <- abs(gradient_prev)
      abs_gradient_prev_sorted <- sort(abs_gradient_prev)

      for (i in seq_len(p)) {
        tmp <- abs_gradient_prev[i] + abs(lambda_prev[ord][i] - lambda)
        neword <- (p:1)[findInterval(tmp,
                                     abs_gradient_prev_sorted,
                                     all.inside = TRUE)]
        # out[i] <- all(tmp < lambda[neword])
        out[i] <- all(cumsum(tmp) < cumsum(lambda[neword]))
      }

      !out
    },

    safe = {
      linear_predictor <- x %*% beta_prev + intercept_prev
      pseudo_gradient_prev <- -(y - linear_predictor)

      ord <- order(abs(gradient_prev), decreasing = TRUE)

      lh <- abs(gradient_prev)

      rh <- lambda[ord] -
        apply(x, 2, norm, "2") *
        norm(pseudo_gradient_prev, "2") *
        (lambda_prev[ord] - lambda[ord])/lambda_prev[ord]

      lh >= rh
    }
  )
}
