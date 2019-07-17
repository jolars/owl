#' Train (tune) model parameter for a golem model
#'
#' @inheritParams golem
#' @param method type of model tuning method
#' @param number number of folds (cross-validation)
#' @param repeats number of repeats for each fold (for repeated *k*-fold
#'   cross validation)
#' @param measure measure to try to optimize
#' @param n_cores number of cores to use for parallel processing
#' @param ... other arguments to pass on to [golem()]
#'
#' @return An object of class `"TrainedGolem"`
#'
#' @export
#'
#' @seealso [parallel::parallel]
trainGolem <- function(x,
                       y,
                       groups = NULL,
                       method = "cv",
                       number = 10,
                       repeats = 1,
                       measure = c("deviance",
                                   "mse",
                                   "mae",
                                   "misclass_error",
                                   "auc"),
                       n_cores = 1,
                       ...) {

  measure <- match.arg(measure)
  method <- match.arg(method)

  n <- NROW(x)
  p <- NCOL(x)

  y <- as.matrix(y)

  stopifnot(NROW(x) > number)

  # get initial lambda sequence
  fit <- golem(x, y, groups = groups, ...)

  name <- fit$penalty$name
  lambda <- fit$penalty$lambda

  is_slope <- name == "slope" || name == "group_slope"

  fold_id <- as.numeric(cut(sample(n), number))

  if (is_slope) {
    sigma <- fit$penalty$sigma
    # fdr <- fit$penalty$fdr

    params <- expand.grid(sigma = sigma)

    result <- expand.grid(fold = seq_len(number),
                          sigma = sigma)

    dim <- c(number, length(sigma), 1)
  } else {
    params <- expand.grid(lambda = lambda)

    dim <- c(number, length(lambda), 1)
  }

  result <- array(dim = dim)

  for (i in seq_len(dim[3])) {
    for (j in seq_len(dim[1])) {
      train_ind <- j == fold_id
      test_ind <- !train_ind

      x_train <- x[train_ind, , drop = FALSE]
      y_train <- y[train_ind, , drop = FALSE]

      if (!is.null(groups)) {
        g_train <- groups[train_ind]
        g_test <- groups[test_ind]
      } else {
        g_train <- NULL
        g_test <- NULL
      }

      x_test <- x[test_ind, , drop = FALSE]
      y_test <- y[test_ind, , drop = FALSE]

      if (is_slope)
        fit <- golem(x_train,
                     y_train,
                     groups = g_train,
                     sigma = sigma,
                     lambda = lambda)
      else
        fit <- golem(x_train,
                     y_train,
                     groups = g_train,
                     lambda = lambda)

      result[j, , i] <-
        score(fit, x_test, y_test, measure)
    }
  }
  fit$tuning_result <- result
}
