test_that("regularization path correctly stops if max_variables reached", {

  x <- scale(heart$x)
  y <- heart$y

  fit <- owl(x, y,
             family = "binomial",
             max_variables = 10,
             intercept = FALSE,
             lambda = "bh",
             standardize_features = FALSE)

  n_var <- max(apply(coef(fit), 2, function(x) {
    length(unique(abs(x[x != 0])))
  }))

  expect_lte(n_var, 10)
})

