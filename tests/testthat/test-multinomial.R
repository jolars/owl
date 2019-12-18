test_that("glmnet and owl return same coefficients for lasso-penalized multinomial model", {
  set.seed(41097)

  x <- scale(wine$x)
  y <- wine$y

  gfit <- glmnet::glmnet(x, y,
                         family = "multinomial",
                         thresh = 1e-10,
                         standardize = FALSE,
                         lambda = 0.2)

  g_coef <- as.matrix(do.call(cbind, coef(gfit)))

  ofit <- owl(x, y,
              family = "multinomial",
              sigma = 0.2*nrow(x),
              tol_rel_gap = 1e-7,
              lambda = rep(1, 3*ncol(x)),
              standardize_features = FALSE)

  expect_equivalent(g_coef, coef(ofit), tol = 1e-4)
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

  x <- scale(cbind(x1, x2), scale = c(norm(x1, "2"), norm(x2, "2")))
  y <- factor(y)

  fit <- glmnet(x, y, family = "multinomial", lambda = 0, thresh = 1e-10,
                standardize = FALSE)
  glmnet_coefs <- as.matrix(do.call(cbind, coef(fit)))

  ofit <- owl(x, y, family = "multinomial", sigma = 1e-9,
              standardize_features = FALSE)

  expect_equivalent(glmnet_coefs, coef(ofit), tol = 1e-4)
})


