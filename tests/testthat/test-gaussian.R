context("gaussian family")

test_that("unregularized gaussian models work as expected", {
  X <- with(mtcars, cbind(cyl, wt, disp, hp, drat))
  y <- mtcars$mpg

  for (standardize in c(TRUE, FALSE)) {
    lm_fit <- lm(y ~ X)
    golem_fit <- golem::golem(X, y, family = "gaussian",
                              standardize = standardize,
                              penalty = slope(sigma = 0),
                              solver = fista(max_passes = 1e6, tol = 1e-8),
                              sigma = 1,
                              intercept = T)

    expect_equivalent(coef(lm_fit),
                      coef(golem_fit),
                      tol = 0.05)
  }
})
