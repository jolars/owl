test_that("glmnet and owl return same coefficients for multinomial model", {
  set.seed(41097)

  x <- scale(wine$x)
  y <- wine$y

  gfit <- glmnet::glmnet(x, y,
                         family = "multinomial",
                         thresh = 1e-10,
                         standardize = FALSE, lambda = 0.2,
                         intercept = TRUE)

  ofit <- owl(x, y,
              family = "multinomial",
              sigma = 0.2*nrow(x),
              tol_rel = 1e-7,
              lambda = rep(1, 3*ncol(x)),
              standardize_features = FALSE)

  expect_equivalent(as.matrix(do.call(cbind, coef(gfit))), coef(ofit),
                    tol = 1e-3)

})
