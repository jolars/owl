context("test-slope")

test_that("multiplication works", {
  if (!requireNamespace("SLOPE", quietly = TRUE))
    skip("SLOPE package not available")

  set.seed(1)

  problem <- golem:::random_problem(10000, 50, sigma=1)

  X <- problem$X
  y <- problem$y

  set.seed(0)
  res_slope <- SLOPE::SLOPE(X, y, sigma = 1)
  set.seed(0)
  res_golem <- golem:::fit(X, y, sigma = 1)

  expect_equivalent(res_slope$beta, res_golem$beta, tol = 0.05)
})

