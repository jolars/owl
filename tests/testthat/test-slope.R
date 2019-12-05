test_that("SLOPE and owl agree for gaussian designs", {
  set.seed(1)

  problem <- owl:::randomProblem(1000, 50, sigma = 1)

  x <- scale(problem$x)
  y <- problem$y

  set.seed(0)
  res_slope <- SLOPE::SLOPE(x, y, sigma = 1, normalize = FALSE)
  set.seed(0)
  g <- owl(x, y,
           intercept = FALSE,
           penalty = "slope",
           sigma = 1,
           diagnostics = TRUE,
           standardize_features = FALSE)

  expect_equivalent(res_slope$beta, coef(g), tol = 1e-5)
})

test_that("SLOPE and owl agree when computing lambda sequences", {
  set.seed(1)

  problem <- owl:::randomProblem(10, 5, sigma = 1)

  x <- problem$x
  y <- problem$y

  for (lambda in c("bhq", "gaussian")) {
    slope_lambda <- SLOPE::SLOPE(x, y, sigma = 1, lambda = lambda)$lambda
    owl_lambda <- owl(x, y, sigma = 1, lambda = lambda,
                      penalty = "slope")$lambda
    expect_equivalent(owl_lambda, slope_lambda)
  }
})


test_that("sigma estimation for SLOPE", {
  set.seed(1)

  problem <- owl:::randomProblem(10, 5, sigma = 1)

  x <- problem$x
  y <- problem$y

  slope_fit <- SLOPE::SLOPE(x, y)
  owl_fit <- owl(x, y, penalty = "slope", intercept = FALSE,
                 sigma = "estimate")

  expect_equal(slope_fit$sigma[length(slope_fit$sigma)],
               owl_fit$sigma)
})

test_that("fitting SLOPE path works", {
  set.seed(239)

  problem <- owl:::randomProblem(100, 5, sigma = 1)

  x <- problem$x
  y <- problem$y

  expect_silent(owl(x, y, penalty = "slope", sigma = "sequence"))
})
