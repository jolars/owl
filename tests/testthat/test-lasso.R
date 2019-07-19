test_that("lasso induces sparse models", {
  set.seed(1)
  x <- with(mtcars, cbind(mpg, disp))
  y <- mtcars$hp

  library(glmnet)
  glmnet.control(fdev = 0)
  glm_fit <- glmnet::glmnet(x, y, standardize = FALSE, nlambda = 10)
  fit <- golem(x,
               y,
               penalty = "lasso",
               standardize_features = FALSE,
               n_lambda = 10)

  expect_equivalent(glm_fit$lambda, fit$penalty$lambda, tol = 1e-4)
})

test_that("lasso and slope fits are equivalent if all lambda are equal", {
  set.seed(1)
  xy <- golem:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  model <- golem(x, y, penalty = "lasso", lambda = 0.2)
  lasso_coef <- coef(model)

  model <- golem(x, y, penalty = "slope",
                 lambda = rep(0.2, NCOL(x))*NROW(x), sigma = 1)
  slope_coef <- coef(model)

  expect_equal(lasso_coef, slope_coef)
})
