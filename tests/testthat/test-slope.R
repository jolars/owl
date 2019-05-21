test_that("SLOPE and golem agree for gaussian designs", {
  set.seed(1)

  problem <- golem:::random_problem(1000, 50, sigma = 1)

  x <- scale(problem$x)
  y <- problem$y

  set.seed(0)
  res_slope <- SLOPE::SLOPE(x, y, sigma = 1, normalize = FALSE)
  set.seed(0)
  res_golem <- golem::golem(x, y, intercept = FALSE,
                            penalty = Slope(sigma = 1),
                            solver = Fista(diagnostics = TRUE),
                            standardize = FALSE)

  expect_equivalent(res_slope$beta, coef(res_golem), tol = 0.01)
})

test_that("SLOPE and golem agree when computing lambda sequences", {
  set.seed(1)

  problem <- golem:::random_problem(10, 5, sigma = 1)

  x <- problem$x
  y <- problem$y

  for (lambda in c("bhq", "gaussian")) {
    slope_lambda <- SLOPE::SLOPE(x, y, sigma = 1, lambda = lambda)$lambda
    golem_lambda <- golem::golem(x, y, sigma = 1, penalty = Slope(lambda = lambda))@penalty@lambda
    expect_equivalent(golem_lambda, slope_lambda)
  }
})
