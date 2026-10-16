# If the validation data lacks covariates the model was developed with,
# predict_risk() and nb() must stop and name every missing covariate.

drop_cols <- function(d, cols) d[, setdiff(names(d), cols), drop = FALSE]

test_that("glm: missing covariates are named", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  v <- drop_cols(s$val, c("sysbp", "pulse"))
  expect_error(predict_risk(fit, v), "missing 2 covariates.*sysbp, pulse")
  expect_error(nb(v, fit, "day30"), "missing 2 covariates.*sysbp, pulse")
})

test_that("one missing covariate uses the singular", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  expect_error(predict_risk(fit, drop_cols(s$val, "age")),
               "missing 1 covariate that the model was developed with: age\\.")
})

test_that("transformed covariates are reported by their column name", {
  s <- gusto_split()
  fit <- glm(day30 ~ log(age) + I(sysbp / 10) + Killip, family = binomial,
             data = s$dev)
  expect_error(predict_risk(fit, drop_cols(s$val, c("age", "sysbp"))),
               "age, sysbp")
})

test_that("a model fitted with y ~ . knows all its covariates", {
  s <- gusto_split()
  cols <- c("day30", "age", "sex", "Killip")
  fit <- glm(day30 ~ ., family = binomial, data = s$dev[, cols])
  expect_error(predict_risk(fit, drop_cols(s$val, "sex")), "sex")
})

test_that("ranger: covariates cannot be read generically, so ranger's own error is shown", {
  skip_if_not_installed("ranger")
  s <- gusto_split()
  rf <- ranger::ranger(factor(day30) ~ age + Killip + sysbp + pulse, data = s$dev,
                       probability = TRUE, num.trees = 50, seed = 1)
  v <- drop_cols(s$val, c("Killip", "pulse"))
  expect_error(predict_risk(rf, v),
               "Could not get predicted risks from this ranger model \\(.+\\)\\. If a covariate is missing")
  expect_error(nb(v, rf, "day30"), "Could not get predicted risks from this ranger model")
})

test_that("randomForest (formula and x/y interface): missing covariates are named", {
  skip_if_not_installed("randomForest")
  s <- gusto_split()
  set.seed(1)
  rf <- randomForest::randomForest(factor(day30) ~ age + Killip + sysbp,
                                   data = s$dev, ntree = 50)
  expect_error(predict_risk(rf, drop_cols(s$val, "sysbp")), "sysbp")
  # the x/y interface keeps no formula: randomForest's own error is shown
  rf_xy <- randomForest::randomForest(x = s$dev[, c("age", "sysbp")],
                                      y = factor(s$dev$day30), ntree = 50)
  expect_error(predict_risk(rf_xy, drop_cols(s$val, "age")),
               "Could not get predicted risks from this randomForest model")
})

test_that("bayespred (bpmfit and bpm): missing covariates are named", {
  skip_if_not_installed("bayespred")
  s <- gusto_split()
  fit <- bayespred::bpmfit(gusto_formula, data = s$dev)
  v <- drop_cols(s$val, c("miloc", "sex"))
  expect_error(predict_risk(fit, v), "miloc, sex")
  expect_error(predict_risk(bayespred::posterior(fit), v), "miloc, sex")
  expect_error(nb(v, fit, "day30"), "miloc, sex")
})

test_that("extra columns and a missing outcome column are fine for prediction", {
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  v <- s$val
  v$something_new <- 1
  expect_silent(predict_risk(fit, v))
  expect_silent(predict_risk(fit, drop_cols(s$val, "day30")))
})
