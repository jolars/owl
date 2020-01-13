#' Generalized Linear Models Penalized with the SLOPE (OWL) Norm
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
#' @param x the feature matrix, which can be either a dense
#'   matrix of the standard *matrix* class, or a sparse matrix
#'   inheriting from [Matrix::sparseMatrix] Data frames will
#'   be converted to matrices internally.
#' @param y the response. For Gaussian models this must be numeric; for
#'   binomial models, it can be a factor.
#' @param family response type. See **Families** for details.
#' @param intercept whether to fit an intercept
#' @param standardize_features whether to standardize features (predictors)
#' @param sigma scale of lambda sequence
#' @param n_sigma length of regularization path
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or the a vector or matrix
#' @param lambda_min_ratio smallest value for `lambda` as a fraction of
#'   `lambda_max`
#' @param q shape of lambda sequence
#' @param max_passes maximum number of passes for optimizer
#' @param diagnostics should diagnostics be saved for the model fit (timings,
#'   primal and dual objectives, and infeasibility)
#' @param screening whether the strong rule for SLOPE be used to screen
#'   variables for inclusion
#' @param verbosity level of verbosity for displaying output from the
#'   program. Setting this to 1 displays information on the path level,
#'   while setting it to 2 displays information also from inside the solver.
#' @param tol_dev_change the regularization path is stopped if the
#'   fractional change in deviance falls below this value. Note that this is
#'   automatically set to 0 if a sigma is manually entered
#' @param tol_dev_ratio the regularization path is stopped if the
#'   deviance ratio
#'   \eqn{1 - \mathrm{deviance}/\mathrm{(null-deviance)}}{1 - deviance/(null deviance)}
#'   is above this threshold
#' @param max_variables criterion for stoping the path in terms of the
#'   maximum number of unique, nonzero
#'   coefficients in absolute value in model
#' @param tol_rel_gap stopping criterion for the duality gap
#' @param tol_infeas stopping criterion for the level of infeasibility
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
#' fit <- owl(bodyfat$x, bodyfat$y)
#'
#' # Binomial response
#'
#' fit <- owl(heart$x, heart$y, family = "binomial")
#'
#' # Poisson response
#'
#' fit <- owl(abalone$x, abalone$y, family = "poisson")
#'
#' # Multinomial response
#'
#' fit <- owl(wine$x, wine$y, family = "multinomial")
#'
owl <- function(x,
                y,
                family = c("gaussian", "binomial", "multinomial", "poisson"),
                intercept = TRUE,
                standardize_features = TRUE,
                sigma = NULL,
                lambda = c("gaussian", "bhq"),
                lambda_min_ratio = if (n < p) 1e-2 else 1e-4,
                n_sigma = 100,
                q = 0.1*min(1, n/p),
                screening = FALSE,
                tol_dev_change = 1e-5,
                tol_dev_ratio = 0.995,
                max_variables = n*m,
                max_passes = 1e6,
                tol_rel_gap = 1e-5,
                tol_infeas = 1e-4,
                diagnostics = FALSE,
                verbosity = 0) {

  ocall <- match.call()

  family <- match.arg(family)

  n <- NROW(x)
  p <- NCOL(x)

  stopifnot(
    is.null(lambda_min_ratio) ||
      (lambda_min_ratio > 0 && lambda_min_ratio < 1),
    max_passes > 0,
    q > 0,
    q < 1,
    length(n_sigma) == 1,
    n_sigma >= 1,
    is.null(lambda) || is.character(lambda) || is.numeric(lambda),
    is.finite(max_passes),
    is.logical(diagnostics),
    is.logical(intercept),
    is.logical(standardize_features),
    tol_rel_gap >= 0,
    tol_infeas >= 0
  )

  fit_intercept <- intercept

  # convert sparse x to dgCMatrix class from package Matrix.
  is_sparse <- inherits(x, "sparseMatrix")

  if (NROW(y) != NROW(x))
    stop("the number of samples in 'x' and 'y' must match")

  if (NROW(y) == 0)
    stop("y is empty")

  if (NROW(x) == 0)
    stop("x is empty")

  if (anyNA(y) || anyNA(x))
    stop("missing values are not allowed")

  if (is_sparse) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  # setup response
  family_choice <- family
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
    response_names <- paste0("y", seq_len(m))

  if (is.null(sigma)) {
    sigma_type <- "auto"
    sigma <- double(n_sigma)
  } else {
    sigma_type <- "user"

    sigma <- as.double(sigma)
    n_sigma <- length(sigma)

    stopifnot(n_sigma > 0)

    # do not stop path early if user requests specific sigma
    tol_dev_change <- 0
    tol_dev_ratio <- 1
    max_variables <- (NCOL(x) + intercept)*m
  }

  n_lambda <- m*p

  if (is.null(lambda)) {
    lambda_type <- "bhq"
    lambda <- double(n_lambda)
  } else if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)
    lambda <- double(n_lambda)
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)

    if (length(lambda) != n_lambda)
      stop("lambda sequence must be as long as there are variables")

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")

    if (any(lambda < 0))
      stop("lambda sequence cannot contain negative values")
  }

  control <- list(family = family$name,
                  fit_intercept = fit_intercept,
                  is_sparse = is_sparse,
                  standardize_features = standardize_features,
                  n_sigma = n_sigma,
                  n_targets = n_targets,
                  screening = screening,
                  sigma = sigma,
                  sigma_type = sigma_type,
                  lambda = lambda,
                  lambda_type = lambda_type,
                  lambda_min_ratio = lambda_min_ratio,
                  q = q,
                  y_center = y_center,
                  y_scale = y_scale,
                  max_passes = max_passes,
                  diagnostics = diagnostics,
                  verbosity = verbosity,
                  max_variables = max_variables,
                  tol_dev_change = tol_dev_change,
                  tol_dev_ratio = tol_dev_ratio,
                  tol_rel_gap = tol_rel_gap,
                  tol_infeas = tol_infeas)

  owlFit <- if (is_sparse) owlSparse else owlDense

  fit <- owlFit(x, y, control)

  lambda <- fit$lambda
  sigma <- fit$sigma
  n_sigma <- length(sigma)
  active_sets <- lapply(drop(fit$active_sets), function(x) drop(x) + 1)
  intercept <- fit$intercepts
  beta <- fit$betas
  nonzeros <- apply(beta, c(2, 3), function(x) abs(x) > 0)

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
                 unique = fit$n_unique,
                 deviance_ratio = as.vector(fit$deviance_ratio),
                 null_deviance = fit$null_deviance,
                 family = family_choice,
                 diagnostics = diagnostics,
                 call = ocall),
            class = c(paste0("Owl", camelCase(family$name)),
                      "Owl"))
}
