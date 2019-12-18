test_that("regularization path correctly stops if max_variables reached", {

  fit <- owl(heart$x, heart$y, family = "binomial", max_variables = 10)
  n_var <- max(apply(coef(fit) != 0, 2, sum))

  expect_lte(n_var, 10)
})

