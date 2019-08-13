test_that("model training with caret works", {
  set.seed(432)

  library(caret)
  library(prague)

  xy <- prague:::randomProblem(1000, 10, q = 0.5)

  ctrl <- trainControl(method = "cv", number = 3)

  set.seed(849)
  train <- train(xy$x,
                 xy$y,
                 method = caretSlopeGolem(),
                 preProc = c("center", "scale"),
                 tuneLength = 2,
                 trControl = ctrl)
  expect_s3_class(train, "train")
})
