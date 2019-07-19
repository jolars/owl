#' Model objects for model tuning with caret
#'
#' This function can be used in a call to [caret::train()] to enable
#' model tuning using caret. This function is currently not exported
#' pending investigation.
#'
#' @return Please see [caret::train()]
#' @keywords internal
caretGolem <- function() {
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

      fit <- golem(x,
                   y,
                   family = fam)

      if (search == "grid") {
        res <- Slope(x,
                     as.numeric(as.factor(y)) - 1,
                     1,
                     "gaussian",
                     "sequence",
                     0.001,
                     len[1],
                     fdr,
                     Binomial())
        out <- expand.grid(sigma = res$sigma,
                           fdr = fdr)
      } else {
        fdr <- stats::runif(0.01, 0.4, len[2])
        res <- Slope(x,
                     as.numeric(as.factor(y)) - 1,
                     1,
                     "gaussian",
                     "sequence",
                     0.001,
                     len[1],
                     fdr,
                     Binomial())

        sigma <- exp(stats::runif(log(min(res$sigma)),
                                  log(max(res$sigma)),
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
      golem::golem(as.matrix(x),
                   y,
                   sigma = param$sigma,
                   fdr = param$fdr,
                   standardize_features = FALSE,
                   ...)
    },

    predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      library(golem)
      stats::predict(modelFit, x = newdata, type = "class")
    },

    prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      library(golem)
      stats::predict(modelFit, x = newdata, type = "response")
    }
  )
}


