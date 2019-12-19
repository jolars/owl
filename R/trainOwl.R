#' Train a owl model
#'
#' This function trains a model fit by [owl()] by tuning its parameters
#' through the resampling method chosen by the `method` argument.
#'
#' Note that by default this method matches all of the available metrics
#' for the given model family against those provided in the argument
#' `measure`. Collecting these measures is not particularly demanding
#' computationally so it is almost always best to leave this argument
#' as it is and then choose which argument to focus on in the call
#' to [plot.TrainedOwl()].
#'
#' @inheritParams owl
#' @param method type of model tuning method
#' @param number number of folds (cross-validation)
#' @param repeats number of repeats for each fold (for repeated *k*-fold
#'   cross validation)
#' @param measure measure to try to optimize; note that you may
#'   supply *multiple* values here and that, by default,
#'   all the possible measures for the given model will be used.
#' @param n_cores number of cores to use for parallel processing
#' @param ... other arguments to pass on to [owl()]
#'
#' @return An object of class `"TrainedOwl"`, with the following slots:
#' \item{summary}{a summary of the results with means, standard errors,
#'                and 0.95 confidence levels}
#' \item{data}{the raw data from the model training}
#' \item{optima}{a `data.frame` of the best (mean) values for the different metrics and their corresponding parameter values}
#' \item{measure}{a `data.frame` listing the used metrics and their labels}
#' \item{model}{the model fit to the entire data set}
#' \item{call}{the call}
#'
#' @export
#'
#' @seealso [parallel::parallel], [plot.TrainedOwl()]
#'
#' @examples
#' # 8-fold cross-validation repeated 5 times
#' tune <- trainOwl(subset(mtcars, select = c("mpg", "drat", "wt")),
#'                    mtcars$hp,
#'                    q = c(0.1, 0.2),
#'                    number = 8,
#'                    repeats = 5)
trainOwl <- function(x,
                     y,
                     q = 0.2,
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

  measure <- match.arg(measure, several.ok = TRUE)
  method <- match.arg(method)

  n <- NROW(x)

  y <- as.matrix(y)

  stopifnot(NROW(x) > number,
            number > 1,
            repeats >= 1)

  # get initial penalty sequence
  fit <- owl(x, y, ...)

  # match measure against accepted measure for the given family
  family <- if (inherits(fit, "OwlGaussian"))
    "gaussian"
  else if (inherits(fit, "OwlBinomial"))
    "binomial"
  else if (inherits(fit, "OwnPoisson"))
    "poisson"

  ok <- switch(family,
               gaussian = c("deviance", "mse", "mae"),
               binomial = c("deviance", "mse", "mae", "accuracy", "auc"))
  measure <- measure[measure %in% ok]

  foldnames <- paste("fold", seq_len(number), sep = "_")

  if (repeats > 1) {
    foldnames <- paste(foldnames,
                       "rep",
                       rep(seq_len(repeats), each = number),
                       sep = "_")
  }

  sigma <- fit$sigma

  result <- array(
    NA,
    dim = c(number*repeats, length(sigma), length(q)),
    dimnames = list(
      foldnames,
      paste("sigma", signif(sigma, 2), sep = "_"),
      paste("q", signif(q, 2), sep = "_")
    )
  )

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
                                     sigma = sigma), list(...))

      for (k in seq_len(d[3])) {

        args$q <- q[k]

        # collect each measure
        for (l in seq_along(measure)) {
          result_list[[l]][j + (i-1)*number, , k] <-
            score(do.call(owl, args), x_test, y_test, measure[l])
        }
      } # loop over each outer parameter (e.g. q)
    } # loop over each fold
  } # repeated cv

  grid <- expand.grid(sigma = sigma, q = q)

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

    optima[[i]] <- data.frame(sigma = s[best_ind, "sigma"],
                              q = s[best_ind, "q"],
                              measure = measure[i],
                              mean = s$mean[best_ind],
                              lo = s$lo[best_ind],
                              hi = s$hi[best_ind])
  }

  optima <- do.call(rbind, optima)

  labels <- vapply(measure, function(m) {
    switch(
      m,
      deviance = {
        if (inherits(fit, "OwlGaussian"))
          "Mean-Squared Error"
        else if (inherits(fit, "OwlBinomial"))
          "Binomial Deviance"
        else if (inherits(fit, "OwlPoisson"))
          "Mean-Squared Error"
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
            class = "TrainedOwl")
}

