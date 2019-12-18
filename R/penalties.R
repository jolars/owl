Slope <- function(x,
                  y,
                  x_center,
                  x_scale,
                  y_scale,
                  standardize_features,
                  lambda = c("gaussian", "bhq"),
                  sigma = NULL,
                  lambda_min_ratio = NULL,
                  n_sigma = 100,
                  fdr = 0.2,
                  family,
                  n_targets) {

  n <- NROW(x)
  p <- NCOL(x)

  if (is.null(lambda_min_ratio))
    lambda_min_ratio <- if (n < p) 0.01 else 0.0001

  if (is.null(lambda))
    lambda <- "gaussian"

  n_lambda <- p*n_targets

  # regularization strength
  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)

    if (lambda_type %in% c("bhq", "gaussian")) {
      q <- 1:n_lambda * fdr/(2*n_lambda)
      lambda <- stats::qnorm(1 - q)

      if (lambda_type == "gaussian" && n_lambda > 1) {
        sum_sq <- 0
        for (i in 2:n_lambda) {
          sum_sq <- sum_sq + lambda[i - 1]^2
          w <- max(1, n - i)
          lambda[i] <- lambda[i]*sqrt(1 + sum_sq/w)
        }
      }

      # ensure non-increasing lambdas
      lambda[which.min(lambda):n_lambda] <- min(lambda)

    }
  } else {
    lambda <- as.double(lambda)

    if (length(lambda) != n_lambda)
      stop("lambda sequence must be as long as there are variables")

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")

    if (any(lambda < 0))
      stop("lambda sequence cannot contain negative values")
  }

  if (is.null(sigma)) {
    lambda_max <- lambdaMax(x,
                            y,
                            x_center,
                            x_scale,
                            y_scale,
                            n_targets,
                            family$name,
                            standardize_features)

    start <- max(sort(lambda_max, decreasing = TRUE)/lambda)

    sigma <- exp(seq(log(start),
                     log(start*lambda_min_ratio),
                     length.out = n_sigma))
  }

  structure(list(name = "slope",
                 tuning_parameters = c("sigma", "fdr"),
                 sigma = sigma,
                 fdr = fdr,
                 lambda = lambda),
            class = c("Slope", "Penalty"))
}
