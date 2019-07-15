test_that("predictions work for all models", {
  set.seed(1)

  for (family in c("gaussian", "binomial")) {
    xy <- golem:::randomProblem(100, 10, response = family)
    x <- xy$x
    y <- xy$y

    fit <- golem(family = family)$fit(x, y)

    for (type in c("link", "response", "class")) {
      if (type == "class" && family == "gaussian")
        next

      expect_silent(fit$predict(x, type = type))
    }
  }
})
