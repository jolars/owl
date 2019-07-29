#' Train a golem model
#'
#' This function trains a model fit by [golem()] by tuning its parameters
#' through the resampling method chosen by the `method` argument.
#'
#' Note that by default this method matches all of the available metrics
#' for the given model family against those provided in the argumet
#' `measure`. Collecting these measures is not particularly demanding
#' computationally so it is almost always best to leave this argument
#' as it is and then choose which argument to focus on in the call
#' to [plot.TrainedGolem()].
#'
#' @inheritParams golem
#' @param method type of model tuning method
#' @param number number of folds (cross-validation)
#' @param repeats number of repeats for each fold (for repeated *k*-fold
#'   cross validation)
#' @param measure measure to try to optimize; note that you may
#'   supply *multiple* values here and that, by default,
#'   all the possible measures for the given model will be used.
#' @param n_cores number of cores to use for parallel processing
#' @param ... other arguments to pass on to [golem()]
#'
#' @return An object of class `"TrainedGolem"`, with the following slots:
#' \item{summary}{a summary of the resutls with means, standard errors,
#'                and 0.95 confidence levels}
#' \item{data}{the raw data from the model training}
#' \item{optima}{a `data.frame` of the best (mean) values for the different metrics and their corresponding parameter values}
#' \item{measure}{a `data.frame` listing the used metrics and their labels}
#' \item{model}{the model fit to the entire data set}
#' \item{call}{the call}
#'
#' @export
#'
#' @seealso [parallel::parallel], [plot.TrainedGolem()]
#'
#' @examples
#' # 8-fold cross-validation repeated 5 times
#' tune <- trainGolem(subset(mtcars, select = c("mpg", "drat", "wt")),
#'                    mtcars$hp,
#'                    fdr = c(0.1, 0.2),
#'                    penalty = "slope",
#'                    number = 8,
#'                    repeats = 5)
trainGolem <- function(x,
                       y,
                       groups = NULL,
                       fdr = 0.2,
                       method = "cv",
                       number = 10,
                       repeats = 1,
                       measure = c("deviance",
                                   "mse",
                                   "mae",
                                   "accuracy",
                                   "auc"),
                       n_cores = 1,
                       ...) {
  ocall <- match.call()

  # dots <- list(...)

  measure <- match.arg(measure, several.ok = TRUE)
  method <- match.arg(method)

  n <- NROW(x)

  y <- as.matrix(y)

  stopifnot(NROW(x) > number,
            number > 1,
            repeats >= 1)

  # get initial penalty sequence
  fit <- golem(x, y, groups = groups, ...)

  # match measure against accepted measure for the given family
  ok <- switch(fit$family$name,
               gaussian = c("deviance", "mse", "mae"),
               binomial = c("deviance", "mse", "mae", "accuracy", "auc"))
  measure <- measure[measure %in% ok]

  is_slope <- fit$penalty$name %in% c("slope", "group_slope")

  foldnames <- paste("fold", seq_len(number), sep = "_")

  if (repeats > 1) {
    foldnames <- paste(foldnames,
                       "rep",
                       rep(seq_len(repeats), each = number),
                       sep = "_")
  }

  if (is_slope) {
    sigma <- fit$penalty$sigma

    result <- array(
      NA,
      dim = c(number*repeats, length(sigma), length(fdr)),
      dimnames = list(
        foldnames,
        paste("sigma", signif(sigma, 2), sep = "_"),
        paste("fdr", signif(fdr, 2), sep = "_")
      )
    )
  } else {
    lambda <- fit$penalty$lambda
    result <- array(
      NA,
      dim = c(number*repeats, length(lambda), 1),
      dimnames = list(foldnames,
                      paste("lambda", signif(lambda, 2), sep = "_"),
                      NULL)
    )
  }

  d <- dim(result)

  result_list <- lapply(seq_along(measure), function(x) result)
  names(result_list) <- measure

  # repeat each k-fold cv
  for (i in seq_len(repeats)) {
    # sample fold indicies for each data point
    fold_id <- as.numeric(cut(sample(n), number))

    # loop over each fold
    for (j in seq_len(number)) {

      train_ind <- j == fold_id
      test_ind <- !train_ind

      x_train <- x[train_ind, , drop = FALSE]
      y_train <- y[train_ind, , drop = FALSE]
      x_test <- x[test_ind, , drop = FALSE]
      y_test <- y[test_ind, , drop = FALSE]

      args <- utils::modifyList(list(x = x_train,
                                     y = y_train,
                                     groups = groups), list(...))

      for (k in seq_len(d[3])) {

        if (is_slope) {
          args$fdr <- fdr[k]
        } else {
          args$lambda <- lambda
        }

        # collect each measure
        for (l in seq_along(measure)) {
          result_list[[l]][j + (i-1)*number, , k] <-
            score(do.call(golem, args), x_test, y_test, measure[l])
        }
      } # loop over each outer parameter (e.g. fdr)
    } # loop over each fold
  } # repeated cv

  if (is_slope) {
    grid <- expand.grid(sigma = sigma, fdr = fdr)
  } else {
    grid <- expand.grid(lambda = lambda)
  }

  n <- number*repeats

  arrange_results <- function(result, grid) {
    out <- grid
    out$mean <- as.vector(apply(result, 3, colMeans))
    out$se <- as.vector(apply(result, c(2, 3), stats::sd))/sqrt(n)
    out$lo <- out$mean - stats::qt(0.975, n - 1)*out$se
    out$hi <- out$mean + stats::qt(0.975, n - 1)*out$se
    out
  }

  summary <- lapply(result_list, arrange_results, grid = grid)

  optima <- vector("list", length(measure))

  for (i in seq_along(measure)) {
    s <- summary[[i]]

    if (measure[i] %in% c("auc", "accuracy"))
      best_ind <- which.max(s$mean)
    else
      best_ind <- which.min(s$mean)

    if (is_slope) {
      optima[[i]] <- data.frame(sigma = s[best_ind, "sigma"],
                                fdr = s[best_ind, "fdr"],
                                measure = measure[i],
                                mean = s$mean[best_ind],
                                lo = s$lo[best_ind],
                                hi = s$hi[best_ind])
    } else {
      optima[[i]] <- data.frame(lambda = summary[[i]][best_ind, "lambda"])
    }
  }

  optima <- do.call(rbind, optima)

  labels <- vapply(measure, function(m) {
    switch(
      m,
      deviance = {
        if (inherits(fit, "GolemGaussian"))
          "Mean-Squared Error"
        else if (inherits(fit, "GolemBinomial"))
          "Binomial Deviance"
      },
      mse = "Mean Squared Error",
      mae = "Mean Absolute Error",
      accuracy = "Accuracy",
      auc = "AUC"
    )
  }, FUN.VALUE = character(1))

  structure(list(summary = summary,
                 data = result_list,
                 optima = optima,
                 measure = data.frame(measure = measure,
                                      label = labels,
                                      row.names = NULL,
                                      stringsAsFactors = FALSE),
                 model = fit,
                 call = ocall),
            class = "TrainedGolem")
}

