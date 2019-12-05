test_that("plotting works", {
  set.seed(1)
  xy <- owl:::randomProblem(100, 2)

  # one parameter
  fit <- owl(xy$x, xy$y, sigma = 0.2)
  expect_silent(dont_plot(fit))

  # more parameters
  fit <- owl(xy$x, xy$y, n_sigma = 10)
  expect_silent(dont_plot(fit))
})
