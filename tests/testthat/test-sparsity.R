test_that("sparse and dense implementations give equivalent results", {
  set.seed(2)
  n <- 1000
  p <- 2

  for (family in c("gaussian", "binomial", "poisson")) {
    for (standardize in c(TRUE, FALSE)) {
      d <- owl:::randomProblem(n, p, 0.5, density = 0.5, response = family)
      x <- d$x
      y <- d$y
      beta <- d$beta

      sparse_fit <-
        owl(x, y, family = family, standardize_features = standardize,
            intercept = FALSE,
            diagnostics = TRUE)

      sparse_coefs <- coef(sparse_fit)

      dense_fit <- owl(as.matrix(x), y, family = family,
                       diagnostics = TRUE,
                       intercept = FALSE,
                       standardize_features = standardize)

      dense_coefs <- coef(dense_fit)

      expect_equal(sparse_coefs, dense_coefs, tol = 1e-6)
    }
  }

  p <- 10
  d <- owl:::randomProblem(n, p, 0.5, density = 0.5, response = "gaussian",
                           n_groups = 5)
  x <- d$x
  y <- d$y
  beta <- d$beta

  for (orthogonalize in c(TRUE, FALSE)) {
    sparse_fit <- owl(x, y,
                      sigma = 1,
                      penalty = "group_slope",
                      groups = d$groups,
                      standardize_features = FALSE,
                      diagnostics = TRUE,
                      orthogonalize = orthogonalize)
    sparse_coefs <- coef(sparse_fit)
    dense_fit <- owl(as.matrix(x), y,
                     sigma = 1,
                     penalty = "group_slope",
                     groups = d$groups,
                     standardize_features = FALSE,
                     diagnostics = TRUE,
                     orthogonalize = orthogonalize)
    dense_coefs <- coef(dense_fit)
    expect_equal(sparse_coefs, dense_coefs, tol = 1e-3)
  }
})
