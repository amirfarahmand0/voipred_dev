# The model must return probabilities. Classes coded as 0/1 numbers pass the
# [0, 1] range check, so they are caught separately.

test_that("only 0s and 1s are rejected as classes, not probabilities", {
  d <- data.frame(y = c(1, 0, 1, 0))
  expect_error(nb(fixed_model(c(1, 0, 0, 0)), d, "y"), "look like predicted classes")
  expect_error(predict_risk(fixed_model(c(0, 0, 0, 0)), d), "look like predicted classes")
  expect_error(predict_risk(fixed_model(c(1, 1, 1, 1)), d), "look like predicted classes")
})

test_that("integer 0/1 predictions are rejected too", {
  d <- data.frame(y = c(1, 0, 1, 0))
  expect_error(predict_risk(fixed_model(c(1L, 0L, 1L, 0L)), d), "look like predicted classes")
})

test_that("probabilities that include some exact 0s and 1s are accepted", {
  d <- data.frame(y = c(1, 0, 1, 0))
  p <- c(1, 0, 0.3, 0.05)
  expect_equal(predict_risk(fixed_model(p), d), p)
  expect_s3_class(nb(fixed_model(p), d, "y", thresholds = 0.1), "voipred_nb")
})

test_that("a random forest regression on a 0/1 outcome still works", {
  skip_if_not_installed("randomForest")
  s <- gusto_split()
  set.seed(1)
  rf <- suppressWarnings(  # randomForest warns that a 0/1 outcome has few values
    randomForest::randomForest(day30 ~ age + Killip + sysbp, data = s$dev, ntree = 50))
  p <- predict_risk(rf, s$val)
  expect_true(all(p >= 0 & p <= 1))
  expect_false(all(p %in% c(0, 1)))
})

test_that("a glm asked for log-odds would be outside [0, 1]: predict_risk uses probabilities", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  expect_true(any(predict(fit, s$val) < 0))           # default type = "link"
  p <- predict_risk(fit, s$val)                        # type = "response"
  expect_true(all(p > 0 & p < 1))
})
