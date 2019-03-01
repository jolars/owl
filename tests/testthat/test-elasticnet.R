context("elastic net penalty")

test_that("elastic net induces sparse models", {
  x <- as.matrix(abalone$x)
  library(glmnet)
  glmnet.control(fdev = 0)
  glm_fit <- glmnet::glmnet(x, abalone$y)

  lambda <- glm_fit$lambda

  fit <- golem(abalone$x, abalone$y,
               penalty = elasticNet)

  expect_equivalent(lambda, fit$penalty$lambda)
})
