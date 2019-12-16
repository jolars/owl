#' Regularized Generalized Linear Models
#'
#' This functions fits a generalized linear model (GLM) using efficient
#' optimization routines suitable to big data problems.
#'
#' The objective for each model is simply the loss function for
#' each family plus a penalty term.
#'
#' @section Families:
#'
#' **Gaussian**
#'
#' The Gaussian model (Ordinary Least Squares) minimizes the following
#' objective.
#'
#' \deqn{
#'   ||y - X\beta||_2^2
#' }{
#'   ||y - X\beta||_2^2
#' }
#'
#' **Binomial**
#'
#' The binomial model (logistic regression) has the following objective.
#'
#' \deqn{
#'   \sum_{i=1}^n \log\left(1+ \exp\left(- y_i \left(x_i^T\beta + \alpha \right) \right) \right)
#' }{
#'   \sum log(1+ exp(- y_i x_i^T \beta))
#' }
#'
#' **Poisson**
#'
#' In poisson regression, we use the following objective.
#'
#' \deqn{
#'   -\sum_{i=1}^n \left(y_i\left(x_i^T\beta + \alpha\right) - \exp\left(x_i^T\beta + \alpha\right)\right)
#' }{
#'   -\sum (y_i(x_i^T\beta + \alpha) - exp(x_i^T\beta + \alpha))
#' }
#'
#'
#' **Multinomial**
#'
#' In multinomial regression, we use the following objective.
#'
#' \deqn{
#'   -\sum_{i=1}^n\left( \sum_{k=1}^m y_{ik}(x_i^T\beta_k + \alpha_k)
#'                      - \log\sum_{k=1}^m \exp(x_i^T\beta_k + \alpha_k) \right)
#' }{
#'   -\sum(y_ik(x_i^T\beta_k + \alpha_k) - log(\sum exp(x_i^T\beta_k + \alpha_k)))
#' }
#'
#' @section Penalties:
#' Models fit by [owl()] can be regularized via several
#' penalties.
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
#' @param x the feature matrix, which can be either a dense
#'   matrix of the standard *matrix* class, or a sparse matrix
#'   inheriting from [Matrix::sparseMatrix] Data frames will
#'   be converted to matrices internally.
#' @param y the response. For Gaussian models this must be numeric; for
#'   binomial models, it can be a factor.
#' @param groups a vector of integers giving the group membership of each
#'   feature (only applies to Group SLOPE)
#' @param family response type. See **Families** for details.
#' @param penalty the regularization penalty to use. See **Penalties** for
#'   details.
#' @param solver the numerical solver to use. See **Solvers** for details.
#' @param intercept whether to fit an intercept
#' @param standardize_features whether to standardize features (predictors)
#' @param orthogonalize whether `x` should be orthogonalized. Note that
#'   setting this to TRUE when `x` is sparse will through an error. (only
#'   applies to Group SLOPE)
#' @param sigma noise estimate (only applies to SLOPE and Group SLOPE)
#' @param n_sigma length of regularization path
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or the a vector or matrix
#' @param lambda_min_ratio smallest value for `lambda` as a fraction of
#'   `lambda_max`
#' @param fdr target false discovery rate (only applies to SLOPE and
#'   Group SLOPE)
#' @param max_passes maximum number of passes for optimizer
#' @param diagnostics should diagnostics be saved for the model fit (timings,
#'   primal and dual objectives, and infeasibility)
#' @param screening_rule type of screening rule to use
#' @param ... arguments passed on to the solver (see [FISTA()], and [ADMM()])
#' @param verbosity level of verbosity for displaying output from the
#'   program. Setting this to 1 displays information on the path level,
#'   while setting it to 2 displays information also from inside the solver.
#'
#' @return An object of class `"Owl"` with the following slots:
#' \item{coefficients}{
#'   a three-dimensional array of the coefficients from the
#'   model fit, including the intercept if it was fit.
#'   There is one row for each coefficient, one column
#'   for each target (dependent variable), and
#'   one slice for each penalty.
#' }
#' \item{nonzeros}{
#'   a three-dimensional boolean array indicating whether a
#'   coefficient was zero or not
#' }
#' \item{lambda}{
#'   the lambda vector that when multiplied by a value in `sigma`
#'   gives the penalty vector at that point along the regularization
#'   path
#' }
#' \item{sigma}{the vector of sigma, indicating the scale of the lambda vector}
#' \item{class_names}{
#'   a character vector giving the names of the classes for binomial and
#'   multinomial families
#' }
#' \item{passes}{the number of passes the solver took at each path}
#' \item{violations}{the number of violations of the screening rule}
#' \item{active_sets}{
#'   a list where each element indicates the indices of the
#'   coefficients that were active at that point in the regularization path
#' }
#' \item{lipschitz_constants}{
#'   the computed Lipschitz constants along the path
#' }
#' \item{diagnostics}{
#'   a `data.frame` of objective values for the primal and dual problems, as
#'   well as a measure of the infeasibility, time, and iteration. Only
#'   available if `diagnostics = TRUE` in the call to [owl()].
#' }
#' \item{call}{the call used for fitting the model}
#' @export
#'
#' @seealso [plot.Owl()], [plotDiagnostics()], [score()], [predict.Owl()],
#'   [trainOwl()]
#'
#' @examples
#'
#' # Gaussian response
#'
#' fit <- owl(abalone$x, abalone$y)
#'
#' # Binomial response
#'
#' fit <- owl(heart$x, heart$y, family = "binomial")
#'
#' # Poisson response
#'
#' fit <- owl(abalone$x, abalone$y, family = "poisson")
#'
owl <- function(x,
                y,
                groups = NULL,
                family = c("gaussian", "binomial", "multinomial", "poisson"),
                penalty = c("slope", "group_slope"),
                solver = c("fista", "admm"),
                intercept = TRUE,
                standardize_features = TRUE,
                orthogonalize = TRUE,
                sigma = c("sequence", "estimate"),
                lambda = NULL,
                lambda_min_ratio = if (NROW(x) < NCOL(x)) 1e-2 else 1e-4,
                n_sigma = 100,
                fdr = 0.2,
                screening_rule = c("none", "strong"),
                max_passes = 1e6,
                diagnostics = FALSE,
                verbosity = 0,
                ...) {

  ocall <- match.call()
  solver_options <- list(...)

  family <- match.arg(family)
  penalty <- match.arg(penalty)
  solver <- match.arg(solver)
  screening_rule <- match.arg(screening_rule)

  if (is.character(sigma)) {
    sigma <- match.arg(sigma)
    sigma_type <- sigma
  } else {
    sigma <- as.double(sigma)
    sigma_type <- "user"
  }

  stopifnot(
    is.null(lambda_min_ratio) ||
      (lambda_min_ratio > 0 && lambda_min_ratio < 1),
    max_passes > 0,
    fdr > 0,
    fdr < 1,
    is.character(sigma) ||
      (all(sigma >= 0 & is.finite(sigma))),
    length(n_sigma) == 1,
    n_sigma >= 1,
    is.null(lambda) ||
      is.character(lambda) || is.numeric(lambda),
    is.finite(max_passes),
    is.logical(diagnostics),
    is.logical(intercept),
    is.logical(standardize_features),
    is.logical(orthogonalize)
  )

  fit_intercept <- intercept
  orthogonalize <- orthogonalize && penalty == "group_slope"
  group_penalty <- penalty == "group_slope"

  n <- NROW(x)
  p <- NCOL(x)

  # convert sparse x to dgCMatrix class from package Matrix.
  is_sparse <- inherits(x, "sparseMatrix")

  if (group_penalty && !(screening_rule == "none"))
    stop("screening rules are only implemented for standard slope penalty")

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

  if (!standardize_features && intercept && solver == "admm")
    stop("'standardize_features' must be set to true when using ADMM")

  if (family != "gaussian" && solver == "admm")
    stop("the ADMM solver currently only works with Gaussian responses")

  if (is_sparse && standardize_features && solver == "admm")
    stop("the ADMM solver does not (yet) work with sparse features",
         "and standardization")

  if (is_sparse && solver == "admm")
    stop("the ADMM solver does not (yet) work with sparse features")

  if (penalty == "group_slope" &&
      standardize_features &&
      is_sparse &&
      orthogonalize)
    stop("orthogonalization is currently not implemented for sparse data ",
         "when standardization is required")

  if (identical(solver_options$line_search, FALSE) && family == "poisson")
    stop("line search must be enabled for poisson models")

  if (identical(solver_options$line_search, FALSE) && penalty == "group_slope")
    stop("line search must be enabled for group slope")

  if (penalty == "group_slope" && family == "multinomial")
    stop("cannot use group slope penalty with multinomial response")

  if (is_sparse) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  # setup response
  family <- switch(family,
                   gaussian = Gaussian(),
                   binomial = Binomial(),
                   multinomial = Multinomial(),
                   poisson = Poisson())

  res <- preProcessResponse(family, y)
  y <- as.matrix(res$y)
  y_center <- res$y_center
  y_scale <- res$y_scale
  class_names <- res$class_names
  m <- n_targets <- res$n_targets
  response_names <- res$response_names
  variable_names <- colnames(x)

  if (is.null(variable_names))
    variable_names <- paste0("V", seq_len(p))
  if (is.null(response_names))
    response_names <- paste0("y", seq_len(NCOL(y)))

  # prepare responses
  intercept_init <- double(m)
  beta_init <- matrix(0, p, m)

  # setup feature matrix
  if (standardize_features) {
    res <- standardize(x)
    x <- res$x
    x_center <- res$x_center
    x_scale <- res$x_scale
  } else {
    x_center <- rep(0, p)
    x_scale <- rep(1, p)
  }

  names(x_center) <- variable_names
  names(x_scale) <- variable_names

  if (group_penalty) {
    res <- groupify(x,
                    x_center,
                    x_scale,
                    groups,
                    standardize_features,
                    orthogonalize)
    x <- res$x
    groups <- res$groups
  }

  if (group_penalty) {
    weights <- groups$wt_per_coef[1:p]

    for (j in seq_len(p))
      x[, j] <- x[, j]/weights[j]

    x_center <- x_center/weights
  } else {
    weights <- rep(1, p)
  }

  if (is_sparse)
    x <- methods::as(x, "dgCMatrix")

  # setup penalty settings
  penalty <- switch(
    penalty,

    slope = Slope(x = x,
                  y = y,
                  x_center = x_center,
                  x_scale = x_scale,
                  y_scale = y_scale,
                  standardize_features = standardize_features,
                  lambda = lambda,
                  sigma = sigma,
                  lambda_min_ratio = lambda_min_ratio,
                  n_sigma = n_sigma,
                  fdr = fdr,
                  family = family,
                  n_targets = n_targets),

    group_slope = GroupSlope(x = x,
                             y = y,
                             x_center = x_center,
                             x_scale = x_scale,
                             y_scale = y_scale,
                             standardize_features = standardize_features,
                             groups = groups,
                             lambda = lambda,
                             sigma = sigma,
                             lambda_min_ratio = lambda_min_ratio,
                             n_sigma = n_sigma,
                             fdr = fdr,
                             family = family,
                             n_targets)
  )

  sigma <- penalty$sigma
  lambda <- penalty$lambda
  n_sigma <- length(penalty$sigma)

  if (family == "gaussian" && !(inherits(penalty, "GroupSlope")))
    solver_options$line_search <- FALSE
  else
    solver_options$line_search <- TRUE

  # setup solver settings
  solver <- switch(solver,
                   fista = do.call(FISTA, solver_options),
                   admm = do.call(ADMM, solver_options))

  control <- list(intercept_init = intercept_init,
                  beta_init = beta_init,
                  family = family,
                  penalty = penalty,
                  solver = solver,
                  groups = groups,
                  fit_intercept = fit_intercept,
                  is_sparse = is_sparse,
                  weights = weights,
                  standardize_features = standardize_features,
                  x_center = x_center,
                  x_scale = x_scale,
                  n_sigma = n_sigma,
                  n_targets = n_targets,
                  sigma_type = sigma_type,
                  screening_rule = screening_rule,
                  sigma = sigma,
                  lambda = lambda,
                  max_passes = max_passes,
                  diagnostics = diagnostics,
                  verbosity = verbosity)

  owlFit <- if (is_sparse) owlSparse else owlDense

  if (sigma_type == "estimate") {
    if (inherits(family, "Gaussian")) {
      control$sigma <- sigma <- stats::sd(y)
      fit <- owlFit(x, y, control)

      S_new <- which(fit$beta != 0)
      S <- c()

      while(!isTRUE(all.equal(S, S_new)) && (length(S_new) > 0)) {
        S <- S_new
        if (length(S) > n) {
          stop("sigma estimation fails because more predictors got ",
               "selected than there are observations.")
        }

        new_x <- x[, S, drop = FALSE]

        if (intercept)
          new_x <- Matrix::cbind2(1, new_x)

        OLS <- stats::lm.fit(as.matrix(new_x), y)
        if (standardize_features) {
          sigma <- sqrt(sum(OLS$residuals^2) / (n - length(S) - 1))
        } else {
          sigma <- sqrt(sum(OLS$residuals^2) / (n - length(S)))
        }

        control$sigma <- penalty$sigma <- sigma
        control$intercept_init <- fit$intercept
        control$beta_init <- as.matrix(fit$beta)

        fit <- owlFit(x, y, control)
        S_new <- which(fit$beta != 0)
      }
    } else if (inherits(family, "Binomial")) {
      control$sigma <- penalty$sigma <- sigma <- 0.5
      fit <- owlFit(x, y, control)
    } else if (inherits(family, "Poisson")) {
      control$sigma <- penalty$sigma <- sigma <- 1
      fit <- owlFit(x, y, control)
    } else if (inherits(family, "Multinomial")) {
      control$sigma <- penalty$sigma <- sigma <- 0.5
      fit <- owlFit(x, y, control)
    }
  } else {
    fit <- owlFit(x, y, control)
  }

  # c++ to R indexing
  active_sets <- lapply(drop(fit$active_sets),
                        function(x) drop(x) + 1)

  intercept <- fit$intercept
  beta <- fit$beta

  # reverse scaling when using group penalty type
  if (group_penalty) {
    beta <- sweep(beta, 1, weights, "/")
    x_center <- x_center*weights
  }

  if (group_penalty && orthogonalize) {
    beta <- unorthogonalize(beta, groups)
  }

  # post-processing to get back non-standardized coefficients
  res <- postProcess(penalty,
                     intercept,
                     beta,
                     x,
                     y,
                     fit_intercept,
                     x_center,
                     x_scale,
                     y_center,
                     y_scale,
                     groups)

  intercept <- res$intercepts
  beta <- res$betas
  nonzeros <- res$nonzeros

  # make sure intercepts sum to zero in multinomial model
  # if (family$name == "multinomial")
  #   intercept <- sweep(intercept, 3, apply(intercept, 3, mean), "-")

  n_sigma <- dim(beta)[3]

  if (fit_intercept) {
    coefficients <- array(NA, dim = c(p + 1, m, n_sigma))

    for (i in seq_len(n_sigma)) {
      coefficients[1, , i] <- intercept[, , i]
      coefficients[-1, , i] <- beta[, , i]
    }
    dimnames(coefficients) <- list(c("(Intercept)", variable_names),
                                   response_names[1:n_targets],
                                   paste0("p", seq_len(n_sigma)))
  } else {
    coefficients <- beta
    dimnames(coefficients) <- list(variable_names,
                                   response_names[1:n_targets],
                                   paste0("p", seq_len(n_sigma)))
  }

  diagnostics <- if (diagnostics) setupDiagnostics(fit) else NULL

  structure(list(coefficients = coefficients,
                 nonzeros = nonzeros,
                 lambda = lambda,
                 sigma = sigma,
                 class_names = class_names,
                 passes = fit$passes,
                 violations = fit$violations,
                 active_sets = active_sets,
                 lipschitz_constants = fit$lipschitz_constants,
                 diagnostics = diagnostics,
                 call = ocall),
            class = c(paste0("Owl", camelCase(family$name)),
                      paste0("Owl", camelCase(penalty$name)),
                      "Owl"))
}
