test_that("glmnet and owl return same coefficients for lasso-penalized multinomial model", {
  set.seed(41097)

  x <- scale(wine$x)
  y <- wine$y

  gfit <- glmnet::glmnet(x, y,
                         family = "multinomial",
                         thresh = 1e-10,
                         standardize = FALSE, lambda = 0.2*2/3,
                         intercept = TRUE)

  g_coef <- as.matrix(do.call(cbind, coef(gfit)))

  g_coef <- g_coef - g_coef[, 3]
  g_coef <- g_coef[, 1:2]

  ofit <- owl(x, y,
              family = "multinomial",
              sigma = 0.2*nrow(x),
              tol_rel = 1e-7,
              lambda = rep(1, 2*ncol(x)),
              standardize_features = FALSE)

  g_coef

  expect_equivalent(as.matrix(do.call(cbind, coef(gfit))), coef(ofit),
                    tol = 1e-3)

})

test_that("glmnet and owl return same unpenalized model", {
  library(glmnet)

  set.seed(-329)

  n <- 1000
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  prob <- matrix(c(rep(1, n),
                   exp(3 + 2*x1 + x2),
                   exp(-1 + x1 - 3*x2)),
                 ncol = 3)
  prob <- sweep(prob, 1, apply(prob, 1, sum), "/")

  y = double(n)

  for (i in 1:n)
    y[i] <- sample(3, 1, replace = TRUE, prob = prob[i, ])

  x <- cbind(x1, x2)
  y <- factor(y)

  fit <- glmnet(x, y, family = "multinomial", lambda = 0, thresh = 1e-10)
  glmnet_coefs <- as.matrix(do.call(cbind, coef(fit)))
  glmnet_coefs[, 1:3] <- glmnet_coefs[, 1:3] - glmnet_coefs[, 3]
  glmnet_coefs <- glmnet_coefs[, 1:2]

  ofit <- owl(x, y, family = "multinomial", sigma = 0)

  expect_equivalent(glmnet_coefs, coef(ofit), tol = 1e-4)
})


