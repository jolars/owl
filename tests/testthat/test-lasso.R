test_that("lasso induces sparse models", {
  x <- with(mtcars, cbind(mpg, disp))
  y <- mtcars$hp

  library(glmnet)
  glmnet.control(fdev = 0)
  glm_fit <- glmnet::glmnet(x, y, standardize = FALSE, nlambda = 10)
  fit <- golem(penalty = "lasso",
               standardize_features = FALSE,
               n_lambda = 10)$fit(x, y)

  expect_equivalent(glm_fit$lambda, fit$penalty$lambda, tol = 1e-4)
})

test_that("lasso and slope fits are equivalent if all lambda are equal", {
  set.seed(1)
  xy <- golem:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  model <- golem(penalty = "lasso", lambda = 0.2)
  model$fit(x, y)
  lasso_coef <- model$coef()

  model$fit(x, y, penalty = "slope",
            lambda = rep(0.2, NCOL(x))*NROW(x), sigma = 1)
  slope_coef <- model$coef()

  expect_equal(lasso_coef, slope_coef)
})
