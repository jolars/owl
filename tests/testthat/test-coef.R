test_that("intepolating coefficients works properly", {
  xy <- golem:::randomProblem(100, 10)

  # check for slope
  fit <- golem(xy$x, xy$y)

  expect_type(coef(fit), "double")
  expect_silent(coef(fit, sigma = c(0.001, 0.04)))

  # check for lasso
  fit <- golem(xy$x, xy$y, penalty = "lasso")

  expect_silent(coef(fit))
  expect_silent(coef(fit, lambda = c(0.2, 20)))

  # check simplify

})
