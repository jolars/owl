#' @include families.R penalties.R solvers.R
Golem <- R6::R6Class(
  "Golem",
  private = list(
    n = NULL,
    p = NULL,
    m = NULL,

    intercept = NULL,
    beta = NULL
  ),

  public = list(
    args = list(),

    family = NULL,
    penalty = NULL,
    solver = NULL,

    coefficients = NULL,
    nonzeros = NULL,
    class_names = NULL,
    diagnostics = NULL,
    passes = 0L,

    fit = function(x, y, groups = NULL, warm_start = TRUE, ...) {

      stopifnot(is.logical(warm_start))

      self$args <- utils::modifyList(self$args, list(...))

      self$args <- checkArgs(self$args)
      args <- self$args

      fit_intercept <- args$intercept
      orthogonalize <- args$orthogonalize && args$penalty == "group_slope"
      group_penalty <- args$penalty == "group_slope"

      if (NROW(y) != NROW(x))
        stop("the number of samples in 'x' and 'y' must match")

      if (NROW(y) == 0)
        stop("the response (y) is empty")

      if (NROW(x) == 0)
        stop("the feature matrix (x) is empty")

      if (anyNA(y) || anyNA(x))
        stop("missing values are not allowed")

      if (anyNA(groups))
        stop("NA values are not allowed in 'groups'")

      private$n <- n <- NROW(x)
      private$p <- p <- NCOL(x)
      private$m <- m <- NCOL(y)

      # use warm starts if applicable
      if (!warm_start || is.null(private$beta)) {

        intercept_init <- double(m)
        beta_init <- matrix(0, p, m)

      } else if (warm_start &&
                 p == NROW(private$beta) &&
                 m == NCOL(private$beta)) {

        # use the coefficients at the end of the regularization path
        # TODO(jolars): perhaps this could be set more intellegently
        #               based on the penalty strength and penalty?
        last <- dim(private$intercept)[3]
        intercept_init <- private$intercept[, , last]
        beta_init <- matrix(private$beta[, , last], p, m)
      }

      # convert sparse x to dgCMatrix class from package Matrix.
      is_sparse <- inherits(x, "sparseMatrix")

      if (is_sparse) {
        x <- methods::as(x, "dgCMatrix")
      } else {
        x <- as.matrix(x)
      }

      if (args$penalty == "group_slope" &&
          args$standardize_features &&
          is_sparse &&
          args$orthogonalize)
        stop("orthogonalization is currently not implemented for sparse data ",
             "when standardization is required")

      # setup family
      family <- switch(args$family,
                       gaussian = Gaussian$new(),
                       binomial = Binomial$new())

      y_center <- y_scale <- double(m)
      n_classes <- integer(0)
      class_names <- character(0)

      c(y, y_center, y_scale, n_classes, class_names) %<-%
        family$preprocessResponse(y)

      x_center <- x_scale <- double(p)
      c(x, x_center, x_scale) %<-% standardize(x, args$standardize_features)

      if (group_penalty) {
        c(x, groups) %<-% groupify(x,
                                   x_center,
                                   x_scale,
                                   groups,
                                   args$standardize_features,
                                   orthogonalize)
      }

      if (group_penalty) {
        weights <- groups$wt_per_coef

        for (j in seq_len(p))
          x[, j] <- x[, j]/weights[j]

        x_center <- x_center/weights
      }

      # setup penalty settings
      penalty <- switch(
        args$penalty,

        slope = Slope$new(x = x,
                          y = y,
                          y_scale = y_scale,
                          lambda = args$lambda,
                          sigma = args$sigma,
                          sigma_min_ratio = args$sigma_min_ratio,
                          n_sigma = args$n_sigma,
                          fdr = args$fdr,
                          family = family),

        group_slope = GroupSlope$new(x = x,
                                     y = y,
                                     y_scale = y_scale,
                                     groups = groups,
                                     lambda = args$lambda,
                                     sigma = args$sigma,
                                     sigma_min_ratio = args$sigma_min_ratio,
                                     n_sigma = args$n_sigma,
                                     fdr = args$fdr,
                                     family = family),

        lasso = Lasso$new(x = x,
                          y = y,
                          y_scale = y_scale,
                          family = family,
                          lambda = args$lambda,
                          lambda_min_ratio = args$lambda_min_ratio,
                          n_lambda = args$n_lambda)
      )

      is_slope <- inherits(penalty, "Slope") || inherits(penalty, "GroupSlope")

      n_penalties <- if (is_slope)
        length(penalty$sigma)
      else
        NCOL(penalty$lambda)

      # setup solver settings
      solver <- Fista$new(tol_rel_gap = args$tol_rel_gap,
                          tol_infeas = args$tol_infeas,
                          max_passes = args$max_passes,
                          diagnostics = args$diagnostics)

      # collect response and variable names (if they are given) and otherwise
      # make new
      response_names <- colnames(y)
      variable_names <- colnames(x)

      if (is.null(variable_names))
        variable_names <- paste0("V", seq_len(p))
      if (is.null(response_names))
        response_names <- paste0("y", seq_len(m))

      lipschitz_constant <-
        family$lipschitzConstant(x,
                                 args$intercept,
                                 x_center,
                                 x_scale,
                                 args$standardize_features)

      control <- list(intercept_init = intercept_init,
                      beta_init = beta_init,
                      family = family,
                      penalty = penalty,
                      solver = solver,
                      groups = groups,
                      fit_intercept = fit_intercept,
                      is_sparse = is_sparse,
                      weights = weights,
                      standardize_features = args$standardize_features,
                      x_scaled_center = x_center/x_scale,
                      lipschitz_constant = lipschitz_constant,
                      n_penalties = n_penalties)

      golemFit <- if (is_sparse) golemSparse else golemDense

      if (is_sparse)
        x <- methods::as(x, "dgCMatrix")

      if (is_slope && args$sigma == "estimate") {
        if (inherits(family, "Gaussian")) {
          control$penalty$sigma <- penalty$sigma <- stats::sd(y)
          res <- golemFit(x, y, control)

          S_new <- which(res$beta != 0)
          S <- c()

          while(!isTRUE(all.equal(S, S_new)) && (length(S_new) > 0)) {
            S <- S_new
            if (length(S) > n) {
              stop("sigma estimation fails because more predictors got ",
                   "selected than there are observations.")
            }

            new_x <- x[, S, drop = FALSE]

            if (args$intercept)
              new_x <- Matrix::cbind2(1, new_x)

            OLS <- stats::lm.fit(as.matrix(new_x), y)
            if (args$standardize_features) {
              sigma <- sqrt(sum(OLS$residuals^2) / (n - length(S) - 1))
            } else {
              sigma <- sqrt(sum(OLS$residuals^2) / (n - length(S)))
            }

            control$penalty$sigma <- penalty$sigma <- sigma
            control$intercept_init <- res$intercept
            control$beta_init <- as.matrix(res$beta)

            res <- golemFit(x, y, control)
            S_new <- which(res$beta != 0)
          }
        } else if (inherits(family, "Binomial")) {
          control$penalty$sigma <- penalty$sigma <- 0.5
          res <- golemFit(x, y, control)
        }
      } else {
        res <- golemFit(x, y, control)
      }

      intercept <- res$intercept
      beta <- res$beta

      # store coefficients for warm starts
      private$intercept <- intercept
      private$beta <- beta

      # reverse scaling when using group penalty type
      if (group_penalty) {
        beta <- sweep(beta, 1, weights, "/")
        x_center <- x_center*weights
      }

      nonzeros <- integer(0)

      if (group_penalty && orthogonalize) {
        beta <- unorthogonalize(beta, groups)
      }

      c(intercept, beta, nonzeros) %<-% penalty$postProcess(intercept,
                                                            beta,
                                                            x,
                                                            y,
                                                            fit_intercept,
                                                            x_center,
                                                            x_scale,
                                                            y_center,
                                                            y_scale,
                                                            groups)

      n_penalties <- dim(beta)[3]

      if (args$intercept) {
        coefficients <- array(NA, dim = c(p + 1, m, n_penalties))

        for (i in seq_len(n_penalties)) {
          coefficients[1, , i] <- intercept[, , i]
          coefficients[-1, , i] <- beta[, , i]
        }
        dimnames(coefficients) <- list(c("(Intercept)", variable_names),
                                       response_names,
                                       paste0("p", seq_len(n_penalties)))
      } else {
        coefficients <- beta
        dimnames(coefficients) <- list(variable_names,
                                       response_names,
                                       paste0("p", seq_len(n_penalties)))
      }

      if (args$diagnostics) {
        self$diagnostics <- Diagnostics$new(res$time,
                                            res$primals,
                                            res$duals,
                                            res$infeasibilities)
      }

      self$coefficients <- coefficients
      self$nonzeros <- nonzeros

      self$family <- family
      self$penalty <- penalty
      self$solver <- solver

      self$class_names <- class_names
      self$passes <- as.integer(res$passes)

      invisible(self)
    },

    coef = function() {
      drop(self$coefficients)
    },

    predict = function(x, type = c("link", "response", "class")) {
      beta <- self$coefficients

      # p <- NROW(beta)
      n <- NROW(x)
      m <- NCOL(beta)
      n_penalties <- dim(beta)[3]

      if (inherits(x, "sparseMatrix"))
        x <- methods::as(x, "dgCMatrix")

      if (inherits(x, "data.frame"))
        x <- as.matrix(x)

      if (names(beta[, , 1])[1] == "(Intercept)") {
        x <- methods::cbind2(1, x)
      }

      lin_pred <- array(dim = c(n, m, n_penalties),
                        dimnames = list(NULL,
                                        dimnames(beta)[[2]],
                                        dimnames(beta)[[3]]))

      for (i in seq_len(n_penalties)) {
        lin_pred[, , i] <- x %*% beta[, , i]
      }

      switch(
        match.arg(type),
        link = lin_pred,
        response = self$family$link(lin_pred),
        class = self$family$predictClass(lin_pred, self$class_names)
      )
    },

    plot = function(...) {
      coefs <- self$coefficients

      if (is.null(coefs))
        stop("nothing to plot since model is yet to be fit")

      p <- NROW(coefs) # number of features
      m <- NCOL(coefs) # number of responses

      is_slope <- inherits(self$penalty, c("Slope", "GroupSlope"))

      args <- list()

      if (is_slope) {
        x <- self$penalty$sigma
        xlab <- expression(sigma)
      } else {
        x <- self$penalty$lambda
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
  )
)

#' Regularized Generalized Linear Models
#'
#' This functions fits a generalized linear model (GLM) using efficient
#' optimization routines suitable to big data problems.
#'
#' The objective for each model is simply the loss function for
#' each family plus a penalty term.
#'
#' @section Methods:
#'
#' \describe{
#'   \item{`fit(x, y, groups = NULL, warm_start = TRUE, ...)`}{
#'     This method fits models specified by [golem()].
#'     \describe{
#'       \item{`x`}{
#'         the feature matrix, which can be either a dense
#'         matrix of the standard *matrix* class, or a sparse matrix
#'         inheriting from [Matrix::sparseMatrix] Data frames will
#'         be converted to matrices internally.
#'       }
#'       \item{`y`}{
#'         the response. For Gaussian models, this must be numeric; for
#'         binomial models, it can be a factor.
#'       }
#'       \item{`groups`}{
#'         a vector of integers giving the group membership of each
#'         feature (only applies to Group SLOPE)
#'       }
#'       \item{`warm_start`}{
#'         whether to use the coefficient estimates from the previous
#'         fit when refitting the model using new data. Provided that
#'         the same penalty (with approximately the same parameters)
#'         is used, setting this to true might lead to
#'         substantial performance boosts.
#'       }
#'       \item{`\dots`}{
#'         arguments that will be used to modify the original model
#'         specification from the call to [golem()]. Note that arguments
#'         pushed through this interface will modify the model, which
#'         might have unintended consequences.
#'       }
#'     }
#'   }
#'   \item{`coef()`}{
#'     Return the coefficients from the model fit (after dropping extraneous
#'     dimensions). If you prefer to always return a three-dimensional array,
#'     call `model$coefficients` instead.
#'   }
#'   \item{`predict(x, type = c("link", "response", "class"))`}{
#'     Return predictions from models fit by [golem()] based on
#'     new data.
#'     \describe{
#'       \item{`x`}{new data to make predictions for.}
#'       \item{`type`}{
#'         type of predictions to make. `"link"` gives the linear
#'         predictors, `"response"` gives the result of applying the link
#'         function, and `"class"` gives class predictions.
#'       }
#'     }
#'   }
#'   \item{`plot(...)`}{
#'     Plot the model's coefficient along the regularization path.
#'     \describe{
#'       \item{`\dots`}{
#'         graphical parameters for the plot passed on
#'         to [lattice::xyplot()].
#'       }
#'     }
#'   }
#' }
#'
#' @section Families:
#'
#' **Gaussian**
#'
#' The Gaussian model (Ordinary Least Squares) minimizes the following
#' objective.
#'
#' \deqn{
#'   ||\boldsymbol{y} - \boldsymbol{X\beta}||_2^2
#' }{
#'   ||y - X\beta||_2^2
#' }
#'
#' **Binomial**
#'
#' The binomial model (logistic regression) has the following objective.
#'
#' \deqn{
#'   \sum_{i=1}^n \log\left[1+ \exp\left(- y_i\boldsymbol{x}_i'\boldsymbol{\beta} \right) \right]
#' }{
#'   \sum log[1+ exp(- y_i x_i' \beta)]
#' }
#'
#' @section Penalties:
#' Models fit by [golem::golem()] can be regularized via several
#' penalties.
#'
#' **Lasso**
#'
#' The Lasso penalizes coefficients using the L1 norm. The penalty
#' term in the lagrangian form of the loss function is
#'
#' \deqn{
#'   \lambda \sum_{j=1}^p |\beta_j|
#' }{
#'   \lambda \sum |\beta|
#' }
#'
#' **SLOPE**
#'
#' SLOPE (Sorted L-One Penalized Estimation) is an extension of the Lasso.
#' Unlike the latter, however, SLOPE uses a non-increasing
#' sequence of \eqn{\lambda}---one
#' for each coefficient. The penalty term looks like
#'
#' \deqn{
#'   \sigma \sum_{i=j}^p \lambda_j |\beta|_{(j)}
#' }{
#'   \sigma \sum \lambda |\beta|(j)
#' }
#'
#' **Group SLOPE**
#'
#' Group SLOPE is an extension of Group LASSO. It applies the following
#' penalty
#'
#' \deqn{
#'   J_\lambda(\beta) = \sum_{j=1}^p \lambda_j |W||\beta||_{I,X}|_{(j)}
#' }{
#'   J_\lambda(\beta) = \sum \lambda_j |W||\beta||_{I,X}|_(j)
#' }
#'
#' @section Solvers:
#' There is currently a single solver available for [golem::golem].
#'
#' **FISTA**
#'
#' FISTA (Fast Iterative Shrinking-Tresholding Algorithm) is an extension
#' of the classical gradient algorithm.
#'
#' @param family response type. See **Families** for details.
#' @param penalty the regularization penalty to use. See **Penalties** for
#'   details.
#' @param solver the numerical solver to use. See **Solvers** for details.
#' @param intercept whether to fit an intercept
#' @param standardize_features whether to standardize features (predictors)
#' @param sigma noise estimate (only applies to SLOPE and Group SLOPE)
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or
#' @param fdr target false discovery rate (only applies to SLOPE and
#'   Group SLOPE)
#' @param n_lambda length of regularization path (only relevant for lasso)
#' @param lambda_min_ratio smallest value for `lambda` as a fraction of
#'   \eqn{\lambda_\mathrm{max}}{\lambda_max} (only applies to lasso)
#' @param max_passes maximum number of passes for optimizer
#' @param diagnostics should diagnostics be saved for the model fit (timings,
#'   primal and dual objectives, and infeasibility)
#' @param orthogonalize whether `x` should be orthogonalized. Note that
#'   setting this to TRUE when `x` is sparse will through an error. (only
#'   applies to Group SLOPE)
#' @param tol_rel_gap relative tolerance threshold for duality gap check
#' @param n_sigma number of sigmas to generate (only relevant for
#'   SLOPE and Group SLOPE)
#' @param sigma_min_ratio smallest value for `sigma` as a fraction of
#    \eqn{\sigma_\mathrm{max}}{\sigma_max}
#' @param tol_infeas tolerance threshold for infeasibility
#'
#' @return An object of class `"Golem"`.
#' @export
#'
#' @examples
#'
#' # Gaussian response, slope penalty (default) --------------------------------
#'
#' # Specify the model
#' gaussian_model <- golem(family = "gaussian")
#'
#' # Fit the model
#' gaussian_model$fit(abalone$x, abalone$y)
#'
#' # Get coefficients from the model fit
#' gaussian_model$coef()
#'
#' # Binomial response, lasso penalty ------------------------------------------
#'
#' binomial_model <- golem(family = "binomial", penalty = "lasso")
#'
#' binomial_model$fit(heart$x, heart$y)
#'
golem <- function(family = c("gaussian", "binomial"),
                  penalty = c("slope", "group_slope", "lasso"),
                  solver = "fista",
                  intercept = TRUE,
                  standardize_features = TRUE,
                  orthogonalize = TRUE,
                  sigma = c("sequence", "estimate"),
                  n_sigma = 100,
                  sigma_min_ratio = NULL,
                  lambda = NULL,
                  n_lambda = 100,
                  lambda_min_ratio = NULL,
                  fdr = 0.2,
                  tol_rel_gap = 1e-6,
                  tol_infeas = 1e-6,
                  max_passes = 1e4,
                  diagnostics = FALSE) {

  out <- Golem$new()

  args <- list(family = family,
               penalty = penalty,
               solver = solver,
               intercept = intercept,
               standardize_features = standardize_features,
               orthogonalize = orthogonalize,
               sigma = sigma,
               n_sigma = n_sigma,
               sigma_min_ratio = sigma_min_ratio,
               lambda = lambda,
               n_lambda = n_lambda,
               lambda_min_ratio = lambda_min_ratio,
               fdr = fdr,
               tol_rel_gap = tol_rel_gap,
               tol_infeas = tol_infeas,
               max_passes = max_passes,
               diagnostics = diagnostics)

  args <- checkArgs(args)

  out$args <- args

  out
}

checkArgs <- function(args) {
  args$family <- match.arg(args$family, c("gaussian", "binomial"))
  args$penalty <- match.arg(args$penalty, c("slope", "group_slope", "lasso"))

  if (is.character(args$sigma))
    args$sigma <- match.arg(args$sigma, c("sequence", "estimate"))

  with(
    args,
    stopifnot(
      is.null(lambda_min_ratio) ||
        (lambda_min_ratio > 0 && lambda_min_ratio < 1),
      tol_rel_gap > 0,
      tol_infeas > 0,
      max_passes > 0,
      n_lambda >= 1,
      fdr > 0,
      fdr < 1,
      is.character(sigma) ||
        (all(sigma >= 0 & is.finite(sigma))),
      is.null(sigma_min_ratio) ||
        (sigma_min_ratio > 0 && sigma_min_ratio < 1),
      length(n_sigma) == 1,
      n_sigma >= 1,
      is.null(lambda) ||
        is.character(lambda) || is.numeric(lambda),
      is.finite(max_passes),
      is.finite(tol_rel_gap),
      is.finite(tol_infeas),
      is.finite(n_lambda),
      n_lambda >= 1,
      is.logical(diagnostics),
      is.logical(intercept),
      is.logical(standardize_features),
      is.logical(orthogonalize)
    )
  )

  args
}
