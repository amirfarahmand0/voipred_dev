# nb() must reproduce dcurves::dca() on GUSTO-I for every kind of model:
# model, treat all and treat none, at thresholds 0, 0.01, ..., 0.99.

expect_same_as_dcurves <- function(res, val) {
  ref <- dcurves_nb(attr(res, "risk"), val, "day30", all_z)
  expect_equal(res$threshold, ref$threshold)
  expect_equal(res$NB_model, ref$NB_model, tolerance = 1e-10)
  expect_equal(res$NB_all, ref$NB_all, tolerance = 1e-10)
  expect_equal(res$NB_none, ref$NB_none, tolerance = 1e-10)
}

test_that("logistic regression (glm) agrees with dcurves", {
  skip_if_not_installed("dcurves")
  s <- gusto_split()
  fit <- glm(gusto_formula, family = binomial, data = s$dev)
  expect_same_as_dcurves(nb(fit, s$val, "day30", all_z), s$val)
})

test_that("random forest (ranger) agrees with dcurves, including tied risks", {
  skip_if_not_installed("dcurves")
  skip_if_not_installed("ranger")
  s <- gusto_split()
  rf <- ranger::ranger(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                       data = s$dev, probability = TRUE, num.trees = 200, seed = 1)
  res <- nb(rf, s$val, "day30", all_z)
  expect_same_as_dcurves(res, s$val)
})

test_that("random forest (randomForest) agrees with dcurves", {
  skip_if_not_installed("dcurves")
  skip_if_not_installed("randomForest")
  s <- gusto_split()
  set.seed(1)
  rf <- randomForest::randomForest(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                                   data = s$dev, ntree = 200)
  res <- nb(rf, s$val, "day30", all_z)
  # randomForest votes are multiples of 1/200, so many risks sit exactly on a
  # threshold: this checks the tie rule against dcurves.
  expect_true(any(attr(res, "risk") %in% all_z))
  expect_same_as_dcurves(res, s$val)
})

test_that("bayespred (posterior mean and point estimate) agrees with dcurves", {
  skip_if_not_installed("dcurves")
  skip_if_not_installed("bayespred")
  s <- gusto_split()
  fit <- bayespred::bpmfit(gusto_formula, data = s$dev,
                           prior = bayespred::log_f(m = 2))
  expect_same_as_dcurves(nb(fit, s$val, "day30", all_z), s$val)
  expect_same_as_dcurves(nb(fit, s$val, "day30", all_z, bayespred_method = "pe"), s$val)
})
