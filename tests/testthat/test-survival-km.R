# Kaplan-Meier / Aalen-Johansen net benefit (the default: `cipcw` not given).

test_that("KM net benefit equals dcurves at every threshold, including NA thresholds", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  for (h in c(2, 5, 10)) {
    res <- quiet_nb(s$val, cx, surv_outcome, horizon = h, thresholds = surv_z)
    expect_same_as_dcurves_surv(res, s$val, "dtimey", "death", h)
  }
})

test_that("the risk by the horizon (treat all) is 1 - Kaplan-Meier", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  res <- nb_val(s$val, cx, surv_outcome, horizon = 5, thresholds = 0.1)
  km <- survival::survfit(survival::Surv(dtimey, death) ~ 1, data = s$val)
  expect_equal(res$outcome$prevalence,
               1 - summary(km, times = 5)$surv, tolerance = 1e-12)
})

test_that("the internal KM and Aalen-Johansen match survival::survfit, with weights", {
  set.seed(3)
  n <- 300
  tt <- round(rexp(n, 0.2), 1)                      # ties on purpose
  st <- sample(0:2, n, replace = TRUE, prob = c(0.4, 0.4, 0.2))
  w <- runif(n, 0.2, 3)
  for (tau in c(1, 3, 6)) {
    km <- survival::survfit(survival::Surv(tt, st > 0) ~ 1, weights = w)
    expect_equal(km_cuminc(tt, as.integer(st > 0), 1, w, tau),
                 1 - summary(km, times = tau, extend = TRUE)$surv, tolerance = 1e-12)
    aj <- survival::survfit(survival::Surv(tt, factor(st, 0:2)) ~ 1, weights = w)
    ps <- summary(aj, times = tau, extend = TRUE)$pstate
    expect_equal(km_cuminc(tt, st, 1, w, tau), ps[, 2], tolerance = 1e-12)
    expect_equal(km_cuminc(tt, st, 2, w, tau), ps[, 3], tolerance = 1e-12)
  }
})

test_that("no follow-up to the horizon among the treated gives NA and a warning", {
  d <- data.frame(time = c(1, 2, 3, 8, 9, 10), status = c(1, 0, 1, 1, 0, 0))
  m <- fixed_model(c(0.9, 0.8, 0.7, 0.2, 0.1, 0.1))  # the high-risk are followed < 5 years
  expect_warning(
    res <- withCallingHandlers(
      nb_val(d, m, ~ survival::Surv(time, status), horizon = 5, thresholds = c(0.05, 0.5)),
      warning = function(w) if (grepl("not a survival model", conditionMessage(w)))
        invokeRestart("muffleWarning")),
    "threshold 0.5 is followed up to the horizon")
  expect_true(is.na(res$estimates$NB_model[2]))
  expect_false(is.na(res$estimates$NB_model[1]))
})

test_that("nobody treated gives NB_model = 0", {
  d <- data.frame(time = c(1, 2, 6, 8), status = c(1, 0, 1, 0))
  res <- suppressWarnings(nb_val(d,
                                 fixed_model(c(0.1, 0.2, 0.15, 0.05)), ~ survival::Surv(time, status), horizon = 5,
                                 thresholds = 0.5))
  expect_equal(res$estimates$NB_model, 0)
})

test_that("a risk at the threshold is treated and a risk of 0 never is", {
  d <- data.frame(time = c(1, 2, 6, 7), status = c(1, 1, 0, 0))
  r1 <- suppressWarnings(nb_val(d, fixed_model(c(0.3, 0.1, 0.3, 0.1)),
                                ~ survival::Surv(time, status), horizon = 5, thresholds = 0.3))
  # treated: rows 1 and 3, one event by 5 of two -> q = 0.5, F = 0.5
  expect_equal(r1$estimates$NB_model, 0.5 * 0.5 - 0.5 * 0.5 * 0.3 / 0.7)
  r0 <- suppressWarnings(nb_val(d, fixed_model(c(0.5, 0, 0.4, 0)),
                                ~ survival::Surv(time, status), horizon = 5, thresholds = 0))
  expect_equal(r0$estimates$NB_model, 0.5 * 0.5)           # rows with risk 0 not treated
})

test_that("weights: all ones equal no weights; integer weights equal duplicated rows", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val[1:400, ]
  z <- seq(0, 0.6, 0.05)
  a <- quiet_nb(v, cx, surv_outcome, horizon = 5, thresholds = z)
  b <- quiet_nb(v, cx, surv_outcome, horizon = 5, thresholds = z, weights = rep(1, 400))
  expect_equal(a$estimates$NB_model, b$estimates$NB_model)
  set.seed(2); w <- rpois(400, 1)
  wtd <- quiet_nb(v, cx, surv_outcome, horizon = 5, thresholds = z, weights = w)
  dup <- quiet_nb(v[rep(seq_len(400), w), ], cx, surv_outcome, horizon = 5, thresholds = z)
  expect_equal(wtd$estimates$NB_model, dup$estimates$NB_model, tolerance = 1e-12)
  expect_equal(wtd$estimates$NB_all, dup$estimates$NB_all, tolerance = 1e-12)
})

test_that("Aalen-Johansen net benefit equals dcurves for both causes", {
  s <- rott_split()
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + size + nodes +
                               grade + hormon, data = s$dev)
  r1 <- quiet_nb(s$val, csc, cr_outcome, horizon = 5, thresholds = surv_z)
  expect_identical(r1$outcome$cause, "recurrence")
  expect_same_as_dcurves_surv(r1, s$val, "cr_time", "cr_status", 5)
  r2 <- quiet_nb(s$val, csc, cr_outcome, horizon = 5, thresholds = surv_z, cause = "death")
  expect_same_as_dcurves_surv(r2, s$val, "cr_time", "cr_status", 5,
                              levels = c("censored", "death", "recurrence"))
})

test_that("with competing risks, weights equal duplicated rows", {
  s <- rott_split()
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  v <- s$val[1:400, ]
  set.seed(5); w <- rpois(400, 1)
  z <- seq(0, 0.6, 0.05)
  wtd <- quiet_nb(v, csc, cr_outcome, horizon = 5, thresholds = z, weights = w)
  dup <- quiet_nb(v[rep(seq_len(400), w), ], csc, cr_outcome, horizon = 5, thresholds = z)
  expect_equal(wtd$estimates$NB_model, dup$estimates$NB_model, tolerance = 1e-12)
})

test_that("KM is unbiased with independent censoring and biased with informative censoring", {
  tau <- 5; z <- c(0.1, 0.2, 0.3, 0.4)
  for (inf in c(FALSE, TRUE)) {
    d <- sim_surv(20000, informative = inf, seed = 11)
    p <- 1 - exp(-0.15 * exp(d$lp1) * tau^1.3)          # true risk by tau
    truth <- true_nb(d, p, z, tau)
    res <- nb_val(d, fixed_model(p), ~ survival::Surv(time, status), horizon = tau,
                  thresholds = z) |> suppressWarnings()
    err <- max(abs(res$estimates$NB_model - truth))
    if (!inf) expect_lt(err, 0.01) else expect_gt(err, 0.02)
  }
})
