context("gaussian family")

test_that("multiplication works", {
  X <- with(mtcars, cbind(cyl, wt, disp, hp, drat))
  y <- mtcars$mpg
  X <- scale(X)

  lm_fit <- lm(y ~ X)
  glm_fit <- glm(y ~ X)

  golem_fit <- golem::golem(X, y, family = "gaussian",
                            lambda = rep(0, ncol(X)),
                            standardize = T,
                            max_passes = 1e6,
                            sigma = 1, intercept = T)

  expect_equivalent(coef(lm_fit),
                    coef(golem_fit),
                    tol = 0.05)
})
