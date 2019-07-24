test_that("plotting works", {
  xy <- golem:::randomProblem(100, 2)

  # one parameter
  fit <- golem(xy$x, xy$y, sigma = 0.2)
  p <- plot(fit)
  expect_s3_class(p, "trellis")

  # more parameters
  fit <- golem(xy$x, xy$y, n_sigma = 10)
  p <- plot(fit)
  expect_s3_class(p, "trellis")

  # lasso
  fit <- golem(xy$x, xy$y, penalty = "lasso")
  expect_silent(dont_plot(plot(fit)))
})
