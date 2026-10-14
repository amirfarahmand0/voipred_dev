test_that("glm: predict_risk() equals predict(type = 'response')", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  expect_equal(predict_risk(fit, s$val),
               unname(predict(fit, s$val, type = "response")))
})

test_that("a glm whose predictions are not probabilities is rejected", {
  d <- data.frame(y = c(1, 2, 3, 4), x = 1:4)
  expect_error(predict_risk(glm(y ~ x, data = d), d), "outside \\[0, 1\\].*not probabilities")
})

test_that("ranger: the risk is the probability of class 1", {
  skip_if_not_installed("ranger")
  s <- gusto_split()
  rf <- ranger::ranger(factor(day30) ~ age + Killip + sysbp + pulse, data = s$dev,
                       probability = TRUE, num.trees = 100, seed = 1)
  p <- predict_risk(rf, s$val)
  expect_equal(p, unname(predict(rf, data = s$val)$predictions[, "1"]))
})

test_that("ranger: numeric 0/1 outcome (unnamed columns) gives the second column", {
  skip_if_not_installed("ranger")
  s <- gusto_split()
  rf <- ranger::ranger(day30 ~ age + Killip + sysbp + pulse, data = s$dev,
                       probability = TRUE, num.trees = 100, seed = 1)
  expect_equal(predict_risk(rf, s$val),
               unname(predict(rf, data = s$val)$predictions[, 2]))
})

test_that("ranger: a class-only forest is rejected with advice", {
  skip_if_not_installed("ranger")
  s <- gusto_split()
  rf <- ranger::ranger(factor(day30) ~ age + Killip, data = s$dev,
                       num.trees = 20, seed = 1)
  expect_error(predict_risk(rf, s$val), "probability ?= ?TRUE")
})

test_that("randomForest: the risk is the probability of class 1", {
  skip_if_not_installed("randomForest")
  s <- gusto_split()
  set.seed(1)
  rf <- randomForest::randomForest(factor(day30) ~ age + Killip + sysbp + pulse,
                                   data = s$dev, ntree = 100)
  p <- predict_risk(rf, s$val)
  expect_equal(p, unname(predict(rf, s$val, type = "prob")[, "1"]))
})

test_that("a forest with more than two classes is rejected", {
  skip_if_not_installed("randomForest")
  s <- gusto_split()
  set.seed(1)
  rf <- randomForest::randomForest(Killip ~ age + sysbp, data = s$dev, ntree = 20)
  expect_error(predict_risk(rf, s$val), "must be binary")
})

test_that("bayespred: the risk is the posterior mean of risk", {
  skip_if_not_installed("bayespred")
  s <- gusto_split()
  fit <- bayespred::bpmfit(gusto_formula, data = s$dev)
  expect_equal(predict_risk(fit, s$val),
               unname(as.vector(predict(fit, s$val, type = "response", method = "pm"))))
  post <- bayespred::posterior(fit)
  expect_equal(predict_risk(post, s$val), predict_risk(fit, s$val))
  expect_error(nb(fit, s$val, "day30", bayespred_method = "pe"), "unused argument")
})

test_that("a numeric or logical vector is rejected", {
  d <- data.frame(x = 1:3)
  expect_error(predict_risk(c(0.1, 0.5, 0.9), d), "give them with `predicted_risk =`")
  expect_error(predict_risk(c(TRUE, FALSE, TRUE), d), "give them with `predicted_risk =`")
})

test_that("an object without a usable predict() gets a clear error", {
  obj <- structure(list(), class = "not_a_model")
  expect_error(predict_risk(obj, data.frame(x = 1)), "Could not get predicted risks")
})

test_that("default method works for a model with predict(type = 'response')", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  class(fit) <- "my_model"
  predict.my_model <- function(object, newdata, type, ...) {
    class(object) <- c("glm", "lm")
    predict(object, newdata, type = type)
  }
  registerS3method("predict", "my_model", predict.my_model)
  expect_equal(predict_risk(fit, s$val),
               unname(predict.my_model(fit, s$val, type = "response")))
})

test_that("missing covariates are caught", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  v <- s$val
  v$age[1:3] <- NA
  expect_error(predict_risk(fit, v), "NA")
})
