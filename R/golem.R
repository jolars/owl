#' @include families.R penalties.R solvers.R
Golem <- R6::R6Class(
  "Golem",
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

    initialize = function(family,
                          penalty,
                          solver,
                          intercept,
                          standardize_features,
                          orthogonalize,
                          sigma,
                          lambda,
                          fdr,
                          n_lambda,
                          lambda_min_ratio,
                          tol_rel_gap,
                          tol_infeas,
                          max_passes,
                          diagnostics) {

      self$args <- list(family = family,
                        penalty = penalty,
                        solver = solver,
                        intercept = intercept,
                        standardize_features = standardize_features,
                        orthogonalize = orthogonalize,
                        sigma = sigma,
                        lambda = lambda,
                        fdr = fdr,
                        n_lambda = n_lambda,
                        lambda_min_ratio = lambda_min_ratio,
                        tol_rel_gap = tol_rel_gap,
                        tol_infeas = tol_infeas,
                        max_passes = max_passes,
                        diagnostics = diagnostics)
    },

    fit = function(x, y, groups = NULL, ...) {

      self$args <- utils::modifyList(self$args, list(...))
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

      # setup penalty settings
      penalty <- switch(
        args$penalty,

        slope = Slope$new(x = x,
                          y = y,
                          lambda = args$lambda,
                          sigma = args$sigma,
                          fdr = args$fdr),

        group_slope = GroupSlope$new(x = x,
                                     y = y,
                                     groups = groups,
                                     lambda = args$lambda,
                                     sigma = args$sigma,
                                     fdr = args$fdr),

        lasso = Lasso$new(x = x,
                          y = y,
                          y_scale = y_scale,
                          family = family,
                          lambda = args$lambda,
                          lambda_min_ratio = args$lambda_min_ratio,
                          n_lambda = args$n_lambda)
      )

      n_penalties <- NCOL(penalty$lambda)

      is_slope <- inherits(penalty, "Slope") || inherits(penalty, "GroupSlope")

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

      if (group_penalty) {
        weights <- groups$wt_per_coef

        for (j in seq_len(p))
          x[, j] <- x[, j]/weights[j]

        x_center <- x_center/weights
      }

      lipschitz_constant <-
        family$lipschitzConstant(x,
                                 args$intercept,
                                 x_center,
                                 x_scale,
                                 args$standardize_features)

      control <- list(family = family,
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

      if (is_slope && is.null(args$sigma)) {
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

      beta <- res$beta

      # reverse scaling when using group penalty type
      if (group_penalty) {
        beta <- sweep(beta, 1, weights, "/")
        x_center <- x_center*weights
      }

      nonzeros <- integer(0)

      if (group_penalty && orthogonalize) {
        beta <- unorthogonalize(beta, groups)
      }

      c(intercept, beta, nonzeros) %<-% penalty$postProcess(res$intercept,
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

      diagnostics <- solver$diagnostics

      if (diagnostics) {
        nl <- length(res$time)
        nn <- lengths(res$time)
        time <- unlist(res$time)
        primal <- unlist(res$primals)
        dual <- unlist(res$duals)
        infeasibility <- unlist(res$infeasibilities)

        diag <- data.frame(time = time,
                           primal = primal,
                           dual = dual,
                           infeasibility = infeasibility,
                           penalty = rep(seq_len(nl), nn))
      } else {
        diag <- data.frame()
      }

      self$diagnostics <- diag

      self$coefficients <- coefficients
      self$nonzeros <- nonzeros

      self$family <- family
      self$penalty <- penalty
      self$solver <- solver

      self$class_names <- class_names
      self$passes <- as.integer(res$passes)

      private$n <- n
      private$p <- p
      private$m <- m

      invisible(self)
    },

    coef = function() {
      drop(self$coefficients)
    },

    predict = function(x, type = c("link", "response", "class")) {
      beta <- self$coef()

      if (inherits(x, "sparseMatrix"))
        x <- methods::as(x, "dgCMatrix")

      if (inherits(x, "data.frame"))
        x <- as.matrix(x)

      if (names(beta)[1] == "(Intercept)") {
        x <- methods::cbind2(1, x)
      }

      lin_pred <- x %*% beta

      switch(
        match.arg(type),
        link = lin_pred,
        response = self$family$link(lin_pred),
        class = self$family$predictClass(lin_pred, self$class_names)
      )
    }
  ),

  private = list(
    n = NULL,
    p = NULL,
    m = NULL
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
#'   \item{`fit(x, y, groups = NULL, ...)`}{
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
#' @param tol_infeas tolerance threshold for infeasibility
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
                  sigma = NULL,
                  lambda = NULL,
                  fdr = 0.2,
                  n_lambda = 100,
                  lambda_min_ratio = NULL,
                  tol_rel_gap = 1e-6,
                  tol_infeas = 1e-6,
                  max_passes = 1e4,
                  diagnostics = FALSE) {

  family <- match.arg(family)
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)

  stopifnot(is.null(lambda_min_ratio) ||
              (lambda_min_ratio > 0 && lambda_min_ratio < 1),
            tol_rel_gap > 0,
            tol_infeas > 0,
            max_passes > 0,
            n_lambda >= 1,
            fdr > 0,
            fdr < 1,
            is.null(sigma) || (sigma >= 0 && is.finite(sigma)),
            is.null(lambda) || is.character(lambda) || is.numeric(lambda),
            is.finite(max_passes),
            is.finite(tol_rel_gap),
            is.finite(tol_infeas),
            is.finite(n_lambda),
            is.logical(diagnostics),
            is.logical(intercept),
            is.logical(standardize_features),
            is.logical(orthogonalize))

  Golem$new(family = family,
            penalty = penalty,
            solver = solver,
            intercept = intercept,
            standardize_features = standardize_features,
            orthogonalize = orthogonalize,
            sigma = sigma,
            lambda = lambda,
            fdr = fdr,
            n_lambda = n_lambda,
            lambda_min_ratio = lambda_min_ratio,
            tol_rel_gap = tol_rel_gap,
            tol_infeas = tol_infeas,
            max_passes = max_passes,
            diagnostics = diagnostics)
}

