test_that("predictions work for all models", {
  set.seed(1)

  for (family in c("gaussian", "binomial", "poisson")) {
    for (penalty in c("slope", "group_slope")) {
      xy <- owl:::randomProblem(100, 10, response = family, n_groups = 3)
      x <- xy$x
      y <- xy$y

      fit <- owl(x,
                 y,
                 groups = xy$groups,
                 family = family,
                 penalty = penalty,
                 orthogonalize = FALSE,
                 diagnostics = TRUE,
                 n_sigma = 10)

      for (type in c("link", "response", "class")) {
        if (type == "class" && family %in% c("gaussian", "poisson"))
          next

        expect_silent(predict(fit, x, type = type))
      }
    }
  }
})

