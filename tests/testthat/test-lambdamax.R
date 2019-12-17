test_that("lambda_max for sparse input, standardizes vs not, and different families", {
  library(Matrix)

  set.seed(-45)

  for (family in c("binomial", "gaussian", "poisson", "multinomial")) {
    xy <- owl:::randomProblem(100, 10, density = 0.3, response = family, amplitude = 1)
    x <- xy$x
    y <- xy$y

    n <- nrow(x)

    x_center <- colMeans(x)
    x_scale <- apply(x, 2, sd)
    xc <- matrix(x_center/x_scale, nrow(x), ncol(x), byrow = TRUE)

    y_scale <- 1

    # make sure y is as lambdaMax expects it
    if (family == "gaussian")
      y <- y - mean(y)
    if (family == "binomial")
      y <- y*2 - 1

    n_targets <- ifelse(family == "multinomial", length(unique(y)) - 1, 1)

    y <- switch(
      family,
      gaussian = y - mean(y),
      binomial = y*2 - 1,
      multinomial = {
        y_map <- matrix(NA, nrow = n, ncol = n_targets)
        for (k in 1:n_targets) {
          y_map[, k] <- as.integer(y == k)
        }
        y_map
      },
      poisson = y
    )

    for (standardize in c(TRUE, FALSE)) {
      if (standardize) {
        x_sparse <- as(sweep(x, 2, x_scale, "/"), "dgCMatrix")
        x_dense <- scale(as.matrix(x))
      } else {
        x_sparse <- x
        x_dense <- as.matrix(x)
      }

      y <- as.matrix(y)

      lsparse <- owl:::lambdaMax(x_sparse, y, x_center, x_scale, y_scale,
                                 n_targets, family, standardize)
      ldense <- owl:::lambdaMax(x_dense, y, x_center, x_scale, y_scale,
                                n_targets, family, standardize)

      expect_equal(lsparse, ldense)
    }
  }
})
