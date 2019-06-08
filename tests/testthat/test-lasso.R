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
