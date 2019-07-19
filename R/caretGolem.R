#' Model objects for model tuning with caret
#'
#' This function can be used in a call to [caret::train()] to enable
#' model tuning using caret. Note that this function does not properly work
#' with sparse feature matrices and standardization due to the way
#' resampling is implemented in caret. So for these cases, please
#' check out [trainGolem()] instead.
#'
#' @return A model description list to be used in the `method` argument
#' in [caret::train()].
#'
#' @seealso [caret::train()], [trainGolem()], [golem()]
#'
#' @export
caretSlopeGolem <- function() {
  list(
    label = "golem",
    library = c("golem", "Matrix"),
    type = c("Regression", "Classification"),

    parameters = data.frame(parameter = c("sigma", "fdr"),
                            class = rep("numeric", 2),
                            label = c("sigma", "False discovery rate")),

    grid = function(x, y, len = NULL, search = "grid") {

      if (length(len) == 1)
        len <- c(len, 1)

      s <- 0.2*(1 - 1/len[2])
      fdr <- seq(-s, s, length.out = len[2]) + 0.2

      numLev <- if (is.character(y) | is.factor(y)) length(levels(y)) else NA

      if (!is.na(numLev))
        fam <- ifelse(numLev > 2, "multinomial", "binomial")
      else
        fam <- "gaussian"

      fit <- golem::golem(x,
                          y,
                          family = fam,
                          penalty = "slope",
                          n_sigma = len[1],
                          standardize_features = FALSE)

      if (search == "grid") {
        out <- expand.grid(sigma = fit$penalty$sigma,
                           fdr = fit$penalty$fdr)
      } else {
        fdr <- stats::runif(0.01, 0.4, len[2])

        sigma <- exp(stats::runif(log(min(fit$sigma)),
                                  log(max(fit$sigma)),
                                  len[1]))

        out <- data.frame(sigma = sigma,
                          fdr = fdr)
      }

      out
    },

    loop = function(grid) {

      fdr <- unique(grid$fdr)
      loop <- data.frame(fdr = fdr)
      loop$sigma <- NA

      submodels <- vector("list", length = length(fdr))

      for (i in seq_along(fdr)) {
        np <- grid[grid$fdr == fdr[i], "sigma"]
        loop$sigma[loop$fdr == fdr[i]] <- np[which.max(np)]
        submodels[[i]] <- data.frame(sigma = np[-which.max(np)])
      }

      list(loop = loop, submodels = submodels)
    },

    fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {

      dots <- list(...)

      numLev <- if(is.character(y) | is.factor(y)) length(levels(y)) else NA

      if (all(names(dots) != "family")) {
        if (!is.na(numLev)) {
          fam <- ifelse(numLev > 2, "multinomial", "binomial")
        } else fam <- "gaussian"
        dots$family <- fam
      }

      dots$x <- x
      dots$y <- y
      dots$sigma <- param$sigma
      dots$fdr <- param$fdr
      dots$standardize_features = FALSE

      do.call(golem::golem, dots)
    },

    predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      library(golem)
      stats::predict(modelFit, x = newdata, type = "class")

      # if (!is.matrix(newdata))
      #   newdata <- Matrix::as.matrix(newdata)
      #
      # if (length(modelFit$obsLevels) < 2) {
      #   out <- stats::predict(modelFit,
      #                         newdata,
      #                         s = modelFit$lambdaOpt)
      # } else {
      #   out <- stats::predict(modelFit,
      #                         newdata,
      #                         s = modelFit$lambdaOpt,
      #                         type = "class")
      # }
      #
      # if (is.matrix(out)) out <- out[, 1]
      #
      # if (!is.null(submodels)) {
      #   if (length(modelFit$obsLevels) < 2) {
      #     tmp <- as.list(as.data.frame(predict(modelFit,
      #                                          newdata,
      #                                          penalty = submodels$lambda)))
      #   } else {
      #     tmp <- predict(modelFit,
      #                    newdata,
      #                    s = submodels$lambda,
      #                    type = "class")
      #     tmp <- if (is.matrix(tmp))
      #       as.data.frame(tmp, stringsAsFactors = FALSE)
      #     else
      #       as.character(tmp)
      #     tmp <- as.list(tmp)
      #   }
      #   out <- c(list(out), tmp)
      # }
      # out

    },

    prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      library(golem)
      stats::predict(modelFit, x = newdata, type = "response")
    },

    sort = function(x) x[order(-x$sigma, x$fdr),],

    tags = c("Generalized Linear Model", "Implicit Feature Selection",
             "L1 Regularization", "Linear Classifier",
             "Linear Regression"),
  )
}


