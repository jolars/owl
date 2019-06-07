#' @include families.R utils.R

Penalty <- R6::R6Class(
  "Penalty",
  list(
    name = NULL,
    lambda = NULL,
    lambda_type = NULL,

    # by default return coefficients on original scale
    postProcess = function(intercepts,
                           betas,
                           x,
                           y,
                           fit_intercept,
                           x_center,
                           x_scale,
                           y_center,
                           y_scale,
                           groups) {

      c(intercepts, betas) %<-% unstandardize(intercepts,
                                              betas,
                                              x,
                                              y,
                                              fit_intercept,
                                              x_center,
                                              x_scale,
                                              y_center,
                                              y_scale)
      selected <- which(betas > 0)
      list(intercepts = intercepts, betas = betas, selected = selected)
    }
  )
)

Slope <- R6::R6Class(
  "Slope",
  inherit = Penalty,
  public = list(
    sigma = NULL,

    initialize = function(x,
                          y,
                          lambda = c("gaussian", "bhq"),
                          sigma = NULL,
                          fdr = 0.2) {

      stopifnot(length(fdr) == 1L, fdr >= 0, fdr <= 1)

      self$name <- "slope"

      n <- NROW(x)
      p <- NCOL(x)

      if (is.null(lambda))
        lambda <- "gaussian"

      # noise estimate
      if (is.null(sigma)) {
        sigma_type <- "auto"
        sigma <- NA_real_
      } else {
        stopifnot(length(sigma) == 1, sigma >= 0, is.finite(sigma))
        sigma_type <- "user"
      }

      # regularization strength
      if (is.character(lambda)) {
        lambda_type <- match.arg(lambda)

        if (lambda_type %in% c("bhq", "gaussian")) {
          q <- 1:p * fdr/(2*p)
          lambda <- stats::qnorm(1 - q)

          if (lambda_type == "gaussian" && p > 1) {
            sum_sq <- 0
            for (i in 2:p) {
              sum_sq <- sum_sq + lambda[i - 1]^2
              w <- max(1, n - i)
              lambda[i] <- lambda[i]*sqrt(1 + sum_sq/w)
            }
          }

          # ensure non-increasing lambdas
          lambda[which.min(lambda):p] <- min(lambda)

        }
      } else {
        lambda <- as.double(lambda)

        if (length(lambda) != p)
          stop("lambda sequence must be as long as there are variables")

        if (is.unsorted(rev(lambda)))
          stop("lambda sequence must be non-increasing")

        if (any(lambda < 0))
          stop("lambda sequence cannot contain negative values")
      }

      self$lambda <- lambda
      self$sigma  <- sigma
    }
  )
)

GroupSlope <- R6::R6Class(
  "GroupSlope",
  inherit = Slope,
  public = list(

    initialize = function(x,
                          y,
                          groups,
                          lambda = c("corrected", "mean", "max"),
                          sigma = NULL,
                          fdr = 0.2) {

      stopifnot(length(fdr) == 1, fdr >= 0, fdr <= 1)

      self$name <- "group_slope"
      group_id <- groups$group_id
      ortho_group_id <- groups$ortho_group_id
      orthogonalize <- groups$orthogonalize
      wt <- groups$wt

      n <- NROW(x)

      n_groups <- length(group_id)

      if (is.null(lambda))
        lambda <- "corrected"

      group_sizes <- if (orthogonalize)
        lengths(ortho_group_id)
      else
        lengths(group_id)

      # noise estimate
      if (is.null(sigma)) {
        sigma_type <- "auto"
        sigma <- NA_real_
      } else {
        stopifnot(length(sigma) == 1, sigma >= 0, is.finite(sigma))
        sigma_type <- "user"
      }

      # regularization strength
      if (is.character(lambda)) {
        lambda_type <- match.arg(lambda)

        if (lambda_type %in% c("max", "mean")) {

          lambda <- lambdaChiOrtho(fdr = fdr,
                                   n.group = n_groups,
                                   wt = wt,
                                   group.sizes = group_sizes,
                                   method = lambda_type)

        } else if (lambda_type == "corrected") {

          # Check for equal group sizes and equal weights
          if ((length(unique(group_sizes)) == 1) & (length(unique(wt)) == 1)) {
            # lambdas of Procedure 6 in Brzyski et. al. (2016)
            m <- unique(group_sizes)
            w <- unique(wt)

            lambda <- lambdaChiEqual(fdr = fdr,
                                     n.obs = n,
                                     n.group = n_groups,
                                     m = m,
                                     w = w)
          } else {
            # lambdas of Procedure 1 in Brzyski et. al. (2016)
            lambda <- lambdaChiMean(fdr = fdr,
                                    n.obs = n,
                                    n.group = n_groups,
                                    group.sizes = group_sizes,
                                    wt = wt)
          }
        }
      } else {
        lambda_type <- "user"
        lambda <- as.double(lambda)

        if (is.unsorted(rev(lambda)))
          stop("lambda sequence must be non-increasing")

        if (any(lambda < 0))
          stop("lambda sequence cannot contain negative values")
      }

      self$lambda <- lambda
      self$sigma <- sigma
    },

    postProcess = function(intercepts,
                           betas,
                           x,
                           y,
                           fit_intercept,
                           x_center,
                           x_scale,
                           y_center,
                           y_scale,
                           groups) {

      # compute group norms ||X_I beta_I||
      group_id <- groups$group_id
      ortho_group_id <- groups$ortho_group_id
      n_groups <- length(group_id)
      orthogonalize <- groups$orthogonalize

      group_norms <- double(n_groups)

      for (i in seq_len(n_groups)) {
        if (orthogonalize) {
          group_norms[i] <- norm(as.matrix(betas[ortho_group_id[[i]]]), "f")
        } else {
          xbetai <- x[, group_id[[i]]] %*% as.matrix(betas[group_id[[i]]])
          group_norms[i] <- norm(as.matrix(xbetai), "f")
        }
      }

      names(group_norms) <- names(group_id)

      res <- unstandardize(intercepts,
                           betas,
                           x,
                           y,
                           fit_intercept,
                           x_center,
                           x_scale,
                           y_center,
                           y_scale)

      res$selected <- which(group_norms > 0)
      res
    }
  )
)

Lasso <- R6::R6Class(
  "Lasso",
  inherit = Penalty,
  public = list(
    name = "lasso",
    lambda_scale = NULL,

    initialize = function(x,
                          y,
                          y_scale,
                          family,
                          lambda = NULL,
                          lambda_min_ratio = "auto",
                          n_lambda = 100) {

      if (lambda_min_ratio == "auto")
        lambda_min_ratio <- ifelse(NROW(x) < NCOL(x), 0.01, 0.0001)

      # lambda (regularization strength)
      if (is.null(lambda)) {
        lambda_max <- family$lambdaMax(x, y, y_scale)

        lambda <- logSeq(lambda_max,
                         lambda_max*lambda_min_ratio,
                         n_lambda)
        lambda_scale <- max(y_scale)

      } else {
        lambda <- as.double(lambda)
        stopifnot(length(lambda) > 0, all(lambda >= 0))
        lambda_scale <- 1
      }

      self$lambda_scale <- lambda_scale/NROW(x)
      self$lambda <- lambda
    }
  )
)
