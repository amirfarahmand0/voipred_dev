# The model must return probabilities. Classes coded as 0/1 numbers pass the
# [0, 1] range check, so they are caught separately.

test_that("only 0s and 1s are rejected as classes, not probabilities", {
  d <- data.frame(y = c(1, 0, 1, 0))
  expect_error(nb_val(d, fixed_model(c(1, 0, 0, 0)), "y"), "look like predicted classes")
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
  expect_s3_class(nb_val(d, fixed_model(p), "y", thresholds = 0.1), "voipred_nb")
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

test_that("values above 1 or below 0 are rejected as not probabilities", {
  d <- data.frame(y = c(1, 0, 1, 0))
  expect_error(predict_risk(fixed_model(c(0.2, 1.3, 0.5, 0.1)), d), "not probabilities")
  expect_error(predict_risk(fixed_model(c(0.2, -0.1, 0.5, 0.1)), d), "not probabilities")
  expect_error(nb_val(d, fixed_model(c(-2.1, 0.4, 1.7, -0.3)), "y"),
               "outside \\[0, 1\\] \\(from -2.1 to 1.7\\).*scores")
})

test_that("log-odds from a glm would be rejected as scores", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  log_odds <- predict(fit, s$val)                     # type = "link"
  expect_error(predict_risk(fixed_model(log_odds), s$val), "scores")
})

test_that("a linear probability model (lm) with predictions below 0 is rejected", {
  s <- gusto_split()
  fit <- lm(gusto_formula, data = s$dev)
  expect_true(any(predict(fit, s$val) < 0))
  expect_error(nb_val(s$val, fit, "day30"), "not probabilities")
})
