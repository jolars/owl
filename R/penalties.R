#' @include families.R utils.R

Penalty <- R6::R6Class(
  "Penalty",
  list(
    name = NULL,
    lambda = NULL,
    lambda_type = NULL,


    # no additional preprocessing by default
    preprocessFeatures = function(x,
                                  x_center,
                                  x_scale,
                                  standardize_features) {
      x
    },

    # by default return coefficients on original scale
    postProcess = function(intercepts,
                           betas,
                           x,
                           y,
                           fit_intercept,
                           x_center,
                           x_scale,
                           y_center,
                           y_scale) {

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
    },

    getWeights = function(x) {
      rep.int(1, NCOL(x))
    }
  )
)

Slope <- R6::R6Class(
  "Slope",
  inherit = Penalty,
  public = list(
    sigma = NULL,
    sigma_type = NULL,
    fdr = NULL,

    initialize = function(lambda = c("gaussian", "bhq"),
                          sigma = NULL,
                          fdr = 0.2) {

      self$name <- "slope"

      if (is.null(lambda))
        lambda <- "gaussian"

      stopifnot(length(fdr) == 1, fdr >= 0, fdr <= 1)

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
        lambda <- NA_real_
      } else {
        lambda_type <- "user"
        lambda <- as.double(lambda)

        if (is.unsorted(rev(lambda)))
          stop("lambda sequence must be non-increasing")

        if (any(lambda < 0))
          stop("lambda sequence cannot contain negative values")
      }

      self$sigma <- sigma
      self$sigma_type <- sigma_type
      self$lambda <- lambda
      self$lambda_type <- lambda_type
      self$fdr <- fdr
    },

    setup = function(family, x, y, y_scale) {

      lambda      <- self$lambda
      lambda_type <- self$lambda_type
      sigma       <- self$sigma
      fdr         <- self$fdr

      n <- NROW(x)
      p <- NCOL(x)

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
      } else {
        if (length(lambda) != p)
          stop("lambda sequence must be as long as there are variables")
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
    do_orthogonalize = NULL,
    groups = NULL,
    group_id = NULL,
    ortho = NULL,
    ortho_groups = NULL,
    ortho_group_id = NULL,
    ortho_group_length = NULL,
    wt = NULL,
    wt_per_coef = NULL,

    initialize = function(groups,
                          lambda = c("corrected", "mean", "max"),
                          sigma = NULL,
                          fdr = 0.2,
                          orthogonalize = TRUE) {

      self$name <- "group_slope"

      if (is.null(lambda))
        lambda <- "corrected"

      stopifnot(length(fdr) == 1, fdr >= 0, fdr <= 1)

      if (anyNA(groups))
        stop("NA values are not allowed in 'groups'")

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
        lambda <- NA_real_
      } else {
        lambda_type <- "user"
        lambda <- as.double(lambda)

        if (is.unsorted(rev(lambda)))
          stop("lambda sequence must be non-increasing")

        if (any(lambda < 0))
          stop("lambda sequence cannot contain negative values")
      }

      self$sigma <- sigma
      self$sigma_type <- sigma_type
      self$lambda <- lambda
      self$lambda_type <- lambda_type
      self$fdr <- fdr
      self$groups <- groups
      self$group_id <- groupID(groups)
      self$ortho_groups <- groups
      self$ortho_group_id <- list()
      self$do_orthogonalize <- orthogonalize
    },

    orthogonalize = function(x,
                             x_center,
                             x_scale,
                             standardize_features) {

      orthogonalize <- self$do_orthogonalize
      group_id      <- self$group_id

      n_groups <- length(group_id)
      group_names <- names(group_id)

      n <- NROW(x)
      # p <- NCOL(x)

      if (orthogonalize) {
        ortho <- orthogonalizeGroups(x,
                                     group_id,
                                     x_center,
                                     x_scale,
                                     standardize_features)

        # determine sizes of orthogonalized groups:
        ortho_group_length <- rep(NA, n_groups)
        names(ortho_group_length) <- group_names

        for (i in 1:n_groups) {
          ortho_group_length[i] <- ncol(ortho[[i]]$Q)
        }

        # overwrite x with a matrix that contains orthogonalizations
        # of the original blocks of x
        # and set up grouping info for the orthogonalized version of x
        # to be used in the optimization
        x <- matrix(nrow = n, ncol = sum(ortho_group_length))
        grp <- rep(NA, sum(ortho_group_length))
        block_end <- cumsum(ortho_group_length)
        block_start <- utils::head(c(1, block_end + 1), n_groups)

        for (i in seq_len(n_groups)) {
          ind <- block_start[i]:block_end[i]
          grp[ind] <- group_names[i]
          x[, ind] <- as.matrix(ortho[[i]]$Q)
        }

        ortho_group_id <- groupID(grp)
        # set prior weights per group:
        wt <- sqrt(ortho_group_length)
        wt_per_coef <- rep(NA, ncol(x))

        for (i in 1:n_groups)
          wt_per_coef[ortho_group_id[[i]]] <- wt[i]

        self$ortho <- ortho
        self$ortho_group_length <- ortho_group_length
        self$ortho_group_id <- ortho_group_id

      } else {
        # set prior weights per group:
        wt <- sqrt(lengths(group_id))
        wt_per_coef <- rep(NA, ncol(x))

        for (i in seq_len(n_groups))
          wt_per_coef[group_id[[i]]] <- wt[i]
      }

      self$wt <- wt
      self$wt_per_coef <- wt_per_coef
      x
    },

    preprocessFeatures = function(x,
                                  x_center,
                                  x_scale,
                                  standardize_features) {
      # orthogonalize (if required)
      x <- self$orthogonalize(x, x_center, x_scale, standardize_features)

      x
    },

    setup = function(family, x, y, y_scale) {

      lambda        <- self$lambda
      lambda_type   <- self$lambda_type
      sigma         <- self$sigma
      fdr           <- self$fdr
      orthogonalize <- self$do_orthogonalize
      group_id      <- self$group_id

      n_groups      <- length(group_id)

      n  <- NROW(x)
      wt <- self$wt

      group_sizes <- if (orthogonalize)
        lengths(self$ortho_group_id)
      else
        lengths(self$group_id)

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

      } else {
        if (length(lambda) != n_groups)
          stop("lambda sequence must be as long as there are variables")
      }

      self$lambda         <- lambda
      self$sigma          <- sigma
      # self$ortho_groups   <- ortho_groups
      # self$ortho_group_id <- ortho_group_id
    },

    postProcess = function(intercepts,
                           betas,
                           x,
                           y,
                           fit_intercept,
                           x_center,
                           x_scale,
                           y_center,
                           y_scale) {

      group_id <- self$group_id
      n_groups <- length(group_id)
      ortho_group_id <- self$ortho_group_id
      orthogonalize <- self$do_orthogonalize

      if (orthogonalize) {

        group_lengths <- lengths(group_id)
        ortho_group_lengths <- lengths(ortho_group_id)

        ortho <- self$ortho

        if (all(group_lengths == ortho_group_lengths)) {

          beta_tilde <- betas

          for (i in seq_len(n_groups)) {
            # c corresponds to the (reordered) group
            # structure in orthogonalized version of X
            ci <- betas[ortho_group_id[[i]]]
            li <- ortho_group_lengths[i]

            if (inherits(ortho[[i]]$R, "sparseMatrix")) {
              bi <- tryCatch({
                Matrix::solve(ortho[[i]]$R, ci)
              }, error = function(err) {
                warning(paste("golem caught an error:", err))
                rep(NA, li)
              })
            } else {
              bi <- tryCatch({
                backsolve(ortho[[i]]$R, ci)
              }, error = function(err) {
                warning(paste("golem caught an error:", err))
                rep(NA, li)
              })
            }

            or <- double(li)

            for (j in seq_len(li))
              or[j] <- which(ortho[[i]]$P == j)

            # beta corresponds to the group structure in the original matrix
            beta_tilde[group_id[[i]]] <- bi[or]
          }
        } else {
          beta_tilde <- NULL
        }
      } else {
        beta_tilde <- betas
      }

      # compute group norms ||X_I beta_I||
      group_norms <- double(n_groups)

      for (i in seq_len(n_groups)) {
        if (orthogonalize) {
          group_norms[i] <- norm(as.matrix(betas[ortho_group_id[[i]]]), "f")
        } else {
          xbetai <- x[, group_id[[i]]] %*% as.matrix(beta_tilde[group_id[[i]]])
          group_norms[i] <- norm(as.matrix(xbetai), "f")
        }
      }

      group_names <- names(group_id)
      names(group_norms) <- group_names

      res <- unstandardize(intercepts,
                           beta_tilde,
                           x,
                           y,
                           fit_intercept,
                           x_center,
                           x_scale,
                           y_center,
                           y_scale)
      res$selected <- which(group_norms > 0)
      res
    },

    getWeights = function(x) {
      self$wt_per_coef
    }
  )
)

Lasso <- R6::R6Class(
  "Lasso",
  inherit = Penalty,
  public = list(
    name = "lasso",
    lambda_scale = NULL,
    lambda_min_ratio = NULL,
    lambda_type = NULL,
    n_lambda = NULL,

    initialize = function(lambda = NULL,
                          lambda_min_ratio = 1e-4,
                          n_lambda = 100) {

      # lambda (regularization strength)
      if (is.null(lambda)) {
        lambda_type <- "auto"
        lambda <- NA_real_

      } else {
        lambda_type <- "user"
        lambda <- as.double(lambda)
        n_lambda <- length(lambda)
        stopifnot(n_lambda > 0, all(lambda >= 0))
      }

      self$lambda <- lambda
      self$lambda_min_ratio <- lambda_min_ratio
      self$lambda_type <- lambda_type
      self$n_lambda <- n_lambda
    },

    setup = function(family, x, y, y_scale) {

      if (self$lambda_min_ratio == "auto") {
        self$lambda_min_ratio <- ifelse(NROW(x) < NCOL(x), 0.01, 0.0001)
      }

      # much of the scaling here is done to ensure that output
      # is equivalent to glmnet
      if (self$lambda_type == "auto") {
        lambda_max <- family$lambdaMax(x, y, y_scale)

        self$lambda <- logSeq(lambda_max,
                              lambda_max*self$lambda_min_ratio,
                              self$n_lambda)
        lambda_scale <- max(y_scale)
      } else {
        lambda_scale <- 1
      }

      self$lambda_scale <- lambda_scale/NROW(x)
    }
  )
)
