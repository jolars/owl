test_that("intepolating coefficients works properly", {
  xy <- prague:::randomProblem(100, 10)

  # check for slope
  fit <- golem(xy$x, xy$y)

  expect_type(coef(fit), "double")
  expect_silent(coef(fit, sigma = c(0.001, 0.04)))

  # check for lasso
  fit <- golem(xy$x, xy$y, penalty = "lasso")

  expect_silent(coef(fit))
  expect_silent(coef(fit, lambda = c(0.2, 20)))

  # penalties are in the path already
  expect_silent(coef(fit, lambda = fit$penalty$lambda[c(2, 3)]))
})

test_that("simplify argument in coef() works as expected", {
  xy <- prague:::randomProblem(100, 10)

  # check for slope
  fit <- golem(xy$x, xy$y, penalty = "lasso")

  # check simplify
  coefs <- coef(fit, simplify = TRUE)
  expect_true(length(dim(coefs)) == 2)

  coefs <- coef(fit, simplify = FALSE)
  expect_true(length(dim(coefs)) == 3)
})

test_that("refitting works if exact = TRUE", {
  xy <- prague:::randomProblem(100, 10)

  # check for slope
  fit <- golem(xy$x, xy$y)

  # check simplify
  coefs <- coef(fit, sigma = 0.4, exact = TRUE, x = xy$x, y = xy$y)

  expect_type(coefs, "double")
})
