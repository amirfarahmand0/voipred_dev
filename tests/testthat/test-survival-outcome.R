# Survival and competing-risks outcomes: parsing, and an error for every
# missing or inconsistent argument.

fit_cox <- function(s) survival::coxph(surv_formula, data = s$dev)

test_that("the survival outcome can be given as a formula, a Surv object or a column", {
  s <- rott_split(); cx <- fit_cox(s)
  r1 <- nb_val(s$val, cx, surv_outcome, horizon = 5, thresholds = c(0.1, 0.3))
  r2 <- nb_val(s$val, cx, survival::Surv(s$val$dtimey, s$val$death), horizon = 5,
               thresholds = c(0.1, 0.3))
  v <- s$val; v$S <- survival::Surv(v$dtimey, v$death)
  r3 <- nb_val(v, cx, "S", horizon = 5, thresholds = c(0.1, 0.3))
  expect_equal(r1$estimates$NB_model, r2$estimates$NB_model)
  expect_equal(r1$estimates$NB_model, r3$estimates$NB_model)
  expect_identical(r1$outcome$type, "survival")
  expect_identical(r1$outcome$horizon, 5)
  expect_identical(r1$outcome$censoring, "km")
})

test_that("a binary outcome can also be given as a formula", {
  toy <- data.frame(y = c(1, 0, 1, 0, 0), p = c(0.9, 0.2, 0.4, 0.05, 0.6))
  expect_equal(nb_val(toy, fixed_model(toy$p), ~ y, thresholds = 0.3)$estimates$NB_model,
               nb_val(toy, fixed_model(toy$p), "y", thresholds = 0.3)$estimates$NB_model)
})

test_that("horizon is required and must be a positive number within follow-up", {
  s <- rott_split(); cx <- fit_cox(s)
  expect_error(nb_val(s$val, cx, surv_outcome), "needs `horizon`")
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = -1), "single positive number")
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = c(2, 5)), "single positive number")
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = NA), "single positive number")
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = "5"), "single positive number")
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = 100),
               "beyond the follow-up of every patient")
})

test_that("survival-only arguments are rejected for a binary outcome", {
  toy <- data.frame(y = c(1, 0, 1, 0), p = c(0.9, 0.2, 0.4, 0.05))
  m <- fixed_model(toy$p)
  expect_error(nb_val(toy, m, "y", horizon = 5), "only used with a survival outcome")
  expect_error(nb_val(toy, m, "y", cause = "x"), "only used with a competing-risks")
  expect_error(nb_val(toy, m, "y", cipcw = ~ p), "only used with a survival outcome")
  expect_error(nb_val(toy, m, "y", cipcw = cipcw_control(~ p, min_prob = 0.05)),
               "only used with a survival outcome")
})

test_that("cause is only for competing risks and must be an event level", {
  s <- rott_split(); cx <- fit_cox(s)
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = 5, cause = "death"),
               "only used with a competing-risks outcome")
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  expect_error(nb_val(s$val, csc, cr_outcome, horizon = 5, cause = "relapse"),
               "must be one of the event types of the outcome: recurrence, death")
})

test_that("only right-censored outcomes are supported", {
  s <- rott_split(); cx <- fit_cox(s)
  v <- s$val; v$start <- 0
  expect_error(nb_val(v, cx, ~ survival::Surv(start, dtimey, death), horizon = 5),
               "Only right-censored")
  expect_error(nb_val(v, cx, ~ survival::Surv(dtimey, dtimey + 1, type = "interval2"),
                      horizon = 5), "Only right-censored")
})

test_that("missing or negative times and missing statuses are rejected", {
  s <- rott_split(); cx <- fit_cox(s)
  v <- s$val; v$dtimey[3] <- NA
  expect_error(nb_val(v, cx, surv_outcome, horizon = 5), "missing times or statuses")
  v <- s$val; v$dtimey[3] <- -1
  expect_error(nb_val(v, cx, surv_outcome, horizon = 5), "negative times")
  expect_error(nb_val(s$val, cx, ~ survival::Surv(dtimey[-1], death[-1]), horizon = 5),
               "rows")
  expect_error(nb_val(s$val, cx, ~ not_a_function(dtimey, death), horizon = 5),
               "Could not evaluate the `outcome` formula")
  expect_error(nb_val(s$val, cx, y ~ survival::Surv(dtimey, death), horizon = 5),
               "must be one-sided")
})

test_that("`cipcw` decides the censoring method, and is checked", {
  s <- rott_split(); cx <- fit_cox(s)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age, data = s$val)
  km <- nb_val(s$val, cx, surv_outcome, horizon = 5, thresholds = 0.2)
  expect_equal(km$outcome$censoring, "km")
  expect_null(km$cipcw)
  ci <- nb_val(s$val, cx, surv_outcome, horizon = 5, thresholds = 0.2, cipcw = ~ age)
  expect_equal(ci$outcome$censoring, "cipcw")
  # there is no `censoring` argument any more
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = 5, censoring = "cipcw"),
               "unused argument")
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = 5, cipcw = list(model = cens)),
               "must be a formula of the covariates of the censoring model")
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = 5, cipcw = "age"),
               "must be a formula")
})

test_that("cipcw_control() checks its inputs", {
  cens <- survival::coxph(survival::Surv(time, status) ~ x, data = data.frame(
    time = 1:10, status = rep(0:1, 5), x = rnorm(10)))
  expect_error(cipcw_control(model = cens, weight_event = 1, weight_horizon = 1),
               "not both")
  expect_error(cipcw_control(), "needs the covariates of the censoring model")
  expect_error(cipcw_control(min_prob = 0.05), "needs the covariates")
  expect_error(cipcw_control(~ x, model = cens), "not several")
  expect_error(cipcw_control(~ x, weight_event = 1, weight_horizon = 1), "not several")
  expect_error(cipcw_control(~ x, refit_function = function(data, weights) cens),
               "are for your own censoring `model`")
  expect_error(cipcw_control(formula = "x"), "must be a formula of the covariates")
  # a model given first is the model
  expect_identical(cipcw_control(cens)$model, cens)
  expect_null(cipcw_control(~ x, min_prob = NULL)$min_prob)
  expect_equal(cipcw_control(~ x)$min_prob, 0.02)
  expect_error(cipcw_control(weight_event = 1:3), "must be given together")
  expect_error(cipcw_control(weight_event = c(1, -1), weight_horizon = c(1, 1)),
               "positive, finite")
  expect_error(cipcw_control(weight_event = c(1, NA), weight_horizon = c(1, 1)),
               "positive, finite")
  expect_error(cipcw_control(weight_event = c(1, Inf), weight_horizon = c(1, 1)),
               "positive, finite")
  expect_error(cipcw_control(weight_event = 1:3, weight_horizon = 1:2), "same length")
  expect_error(cipcw_control(~ x, min_prob = 0), "strictly between 0 and 1")
  expect_error(cipcw_control(~ x, min_prob = 1.2), "strictly between 0 and 1")
  expect_error(cipcw_control(model = c(0.2, 0.3)), "must be a fitted censoring model")
  expect_error(cipcw_control(ties = "middle"), "should be one of")
  expect_output(print(cipcw_control(model = cens, min_prob = 0.02)),
                "censoring model \\(coxph\\), truncation: min_prob = 0.02")
  expect_output(print(cipcw_control(~ x + z)),
                "Cox censoring model ~x \\+ z, fitted by voipred, truncation: min_prob = 0.02")
  expect_output(print(cipcw_control(model = cens, min_prob = NULL)), "truncation: off")
})

test_that("cIPCW weights must have one value per row", {
  s <- rott_split(); cx <- fit_cox(s)
  expect_error(nb_val(s$val, cx, surv_outcome, horizon = 5, cipcw = cipcw_control(weight_event = rep(1, 10), weight_horizon = rep(1, 10))),
               "length 10 but `data` has 1482 rows")
})

test_that("print shows the outcome type, horizon and estimator", {
  s <- rott_split(); cx <- fit_cox(s)
  r <- nb_val(s$val, cx, surv_outcome, horizon = 5, thresholds = c(0.1, 0.2))
  expect_output(print(r), "survival, horizon 5, .*censoring KM")
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  r2 <- nb_val(s$val, csc, cr_outcome, horizon = 5, thresholds = c(0.1, 0.2), cause = "death")
  expect_output(print(r2), "competing risks \\(cause death\\), horizon 5")
})
