test_that("predictions work for all models", {

  for (family in c("gaussian", "binomial")) {
    xy <- golem:::randomProblem(100, 10, response = family)
    x <- xy$x
    y <- xy$y

    fit <- golem(x, y, family = family)

    for (type in c("link", "response", "class")) {
      if (type == "class" && family == "gaussian")
        next

      expect_silent(predict(fit, x, type = type))
    }
  }
})
