test_that("predictions work for all models", {
  set.seed(1)

  for (family in c("gaussian", "binomial")) {
    for (penalty in c("slope", "group_slope")) {
      xy <- owl:::randomProblem(100, 10, response = family, n_groups = 2)
      x <- xy$x
      y <- xy$y

      fit <- owl(x,
                 y,
                 groups = xy$groups,
                 family = family,
                 penalty = penalty,
                 n_sigma = 10)

      for (type in c("link", "response", "class")) {
        if (type == "class" && family == "gaussian")
          next

        expect_silent(predict(fit, x, type = type))
      }
    }
  }
})

