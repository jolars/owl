test_that("output from unregularized poisson model matches glm", {

  set.seed(654)

  for (intercept in c(TRUE, FALSE)) {
    xy <- owl:::randomProblem(100, 10, response = "poisson", amplitude = 1)
    x <- xy$x
    y <- xy$y

    x <- scale(x)

    owl_fit <- owl(x, y,
                   family = "poisson",
                   sigma = 1e-8,
                   intercept = intercept,
                   center = FALSE,
                   scale = "none")

    glm_fit <- if (intercept)
      glm(y ~ 1 + ., data = data.frame(y = y, x), family = "poisson")
    else
      glm(y ~ 0 + ., data = data.frame(y = y, x), family = "poisson")

    expect_equivalent(coef(glm_fit), coef(owl_fit), tol = 1e-5)
  }

})

test_that("owl reproduces lasoso fit when all lambda are equal", {

  set.seed(0978213)

  xy <- owl:::randomProblem(100, 10, response = "poisson", amplitude = 1)
  x <- xy$x
  y <- xy$y

  x <- scale(x)

  n <- nrow(x)
  p <- ncol(x)

  library(glmnet)

  sigma <- 1

  gnt_fit <- glmnet(x, y, family = "poisson", lambda = sigma, standardize = FALSE)
  owl_fit <- owl(x, y,
                 family = "poisson",
                 sigma = sigma,
                 lambda = rep(1, p),
                 scale = "none",
                 center = FALSE)

  expect_equivalent(as.double(coef(gnt_fit)), coef(owl_fit), tol = 1e-3)
})
