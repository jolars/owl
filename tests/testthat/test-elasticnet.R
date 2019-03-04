context("elastic net penalty")

test_that("elastic net induces sparse models", {
  x <- scale(as.matrix(abalone$x))
  y <- abalone$y
  library(glmnet)
  glmnet.control(fdev = 0)
  glm_fit <- glmnet::glmnet(x, y, standardize = FALSE)

  lambda <- glm_fit$lambda

  fit <- golem(x, y, penalty = elasticNet, standardize = FALSE)

  expect_equivalent(lambda, fit$penalty$lambda)
})
