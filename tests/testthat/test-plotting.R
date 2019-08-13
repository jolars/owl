test_that("plotting works", {
  set.seed(1)
  xy <- prague:::randomProblem(100, 2)

  # one parameter
  fit <- golem(xy$x, xy$y, sigma = 0.2)
  expect_silent(dont_plot(fit))

  # more parameters
  fit <- golem(xy$x, xy$y, n_sigma = 10)
  expect_silent(dont_plot(fit))

  # lasso
  fit <- golem(xy$x, xy$y, penalty = "lasso")
  expect_silent(dont_plot(fit))
})
