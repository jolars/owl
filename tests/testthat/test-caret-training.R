test_that("model training with caret works", {
  set.seed(432)

  library(caret)
  library(owl)

  xy <- owl:::randomProblem(1000, 2, q = 1)

  ctrl <- trainControl(method = "cv", number = 3)

  train <- train(xy$x,
                 xy$y,
                 method = caretSlopeOwl(),
                 preProc = c("center", "scale"),
                 tuneLength = 2,
                 trControl = ctrl)
  expect_s3_class(train, "train")
})
