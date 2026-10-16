# Competing risks: cumulative incidence models, cause labels, and the checks
# on mismatched models and outcomes. Net benefit vs dcurves (Aalen-Johansen).

local_base_formula()

Hist <- prodlim::Hist      # FGR's predict() looks up Hist() in the formula's environment
cr_form <- Hist(cr_time, cr_code) ~ age + size + nodes + grade + hormon
# riskRegression's predict.FGR prints the covariate names; keep the log clean.
fgr_risk <- function(m, nd, cause) { utils::capture.output(
  r <- riskRegression::predictRisk(m, nd, times = 5, cause = cause)); r }
num_vars <- c("age", "nodes", "grade", "hormon", "pgr")

check_cr <- function(model, s, own_risk, cause = NULL, horizon = 5) {
  expect_no_warning(res <- suppress_matching(
    nb(s$val, model, cr_outcome, horizon = horizon, thresholds = surv_z, cause = cause),
    "^No patient treated"))
  r <- res$risk
  expect_true(all(r >= 0 & r <= 1))
  expect_equal(r, as.vector(own_risk), tolerance = 1e-6)
  lev <- if (identical(cause, "death")) c("censored", "death", "recurrence") else NULL
  expect_same_as_dcurves_surv(res, s$val, "cr_time", "cr_status", horizon, levels = lev)
  invisible(res)
}

test_that("riskRegression::CSC for both causes, with numeric or labelled causes (standard route)", {
  s <- rott_split()
  csc <- riskRegression::CSC(cr_form, data = s$dev)
  r1 <- check_cr(csc, s, riskRegression::predictRisk(csc, s$val, times = 5, cause = 1))
  expect_identical(r1$outcome$cause, "recurrence")
  r2 <- check_cr(csc, s, riskRegression::predictRisk(csc, s$val, times = 5, cause = 2),
                 cause = "death")
  expect_identical(r2$outcome$cause, "death")
  # the cumulative incidences of the two causes cannot add up to more than 1
  expect_true(all(r1$risk + r2$risk <= 1 + 1e-12))
  # causes coded by label in the model
  dev <- s$dev
  csl <- riskRegression::CSC(prodlim::Hist(cr_time, cr_status, cens.code = "censored") ~
                               age + size + nodes + grade + hormon, data = dev)
  check_cr(csl, s, riskRegression::predictRisk(csc, s$val, times = 5, cause = 2),
           cause = "death")
})

test_that("riskRegression::FGR (standard route) and cmprsk::crr (voipred method) agree", {
  s <- rott_split()
  fg <- riskRegression::FGR(Hist(cr_time, cr_code) ~ age + nodes + grade + hormon + pgr,
                            data = s$dev, cause = 1)
  X <- as.matrix(s$dev[, num_vars])
  cr <- cmprsk::crr(s$dev$cr_time, s$dev$cr_code, cov1 = X, failcode = 1, cencode = 0)
  own <- fgr_risk(fg, s$val, 1)
  check_cr(fg, s, own)
  check_cr(cr, s, own)
  fg2 <- riskRegression::FGR(Hist(cr_time, cr_code) ~ age + nodes + grade + hormon + pgr,
                             data = s$dev, cause = 2)
  check_cr(fg2, s, fgr_risk(fg2, s$val, 2),
           cause = "death")
})

test_that("randomForestSRC competing-risks forest (standard route, cause by position)", {
  s <- rott_split()
  rf <- randomForestSRC::rfsrc(Surv(cr_time, cr_code) ~ age + size + nodes + grade + hormon,
                               data = s$dev, ntree = 200, seed = -1)
  check_cr(rf, s, riskRegression::predictRisk(rf, s$val, times = 5, cause = 1))
  check_cr(rf, s, riskRegression::predictRisk(rf, s$val, times = 5, cause = 2),
           cause = "death")
})

test_that("survival's multi-state coxph: Aalen-Johansen from cause-specific Cox models", {
  s <- rott_split()
  dev <- s$dev
  dev$st <- factor(dev$cr_code, 0:2, c("censored", "recurrence", "death"))
  ms <- survival::coxph(survival::Surv(cr_time, st) ~ age + size + nodes + grade + hormon,
                        data = dev, id = pid, ties = "breslow", model = TRUE)
  # riskRegression cannot predict from it
  expect_error(riskRegression::predictRisk(ms, s$val, times = 5, cause = 1))
  # reference 1: survfit's own Aalen-Johansen state probabilities
  ps <- summary(survival::survfit(ms, newdata = s$val, ctype = 1, stype = 1),
                times = 5)$pstate
  r1 <- check_cr(ms, s, ps[1, , 2])
  r2 <- check_cr(ms, s, ps[1, , 3], cause = "death")
  expect_true(all(r1$risk + r2$risk <= 1 + 1e-12))
  # reference 2: the same cause-specific Cox models fitted by riskRegression::CSC
  csc <- riskRegression::CSC(Hist(cr_time, cr_code) ~ age + size + nodes + grade + hormon,
                             data = s$dev, ties = "breslow")
  expect_equal(r1$risk, as.vector(riskRegression::predictRisk(csc, s$val, times = 5,
                                                                      cause = 1)), tolerance = 1e-7)
  expect_equal(r2$risk, as.vector(riskRegression::predictRisk(csc, s$val, times = 5,
                                                                      cause = 2)), tolerance = 1e-7)
  # Efron ties: the baseline hazards differ slightly between the two packages
  mse <- survival::coxph(survival::Surv(cr_time, st) ~ age + nodes, data = dev, id = pid,
                         model = TRUE)
  cse <- riskRegression::CSC(Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  re <- nb(s$val, mse, cr_outcome, horizon = 5, thresholds = 0.2)$risk
  expect_lt(max(abs(re - riskRegression::predictRisk(cse, s$val, times = 5, cause = 1))), 1e-3)
})

test_that("multi-state coxph: causes matched by label or position, and the checks", {
  s <- rott_split()
  dev <- s$dev
  dev$st <- factor(dev$cr_code, 0:2, c("censored", "recurrence", "death"))
  ms <- survival::coxph(survival::Surv(cr_time, st) ~ age + nodes, data = dev, id = pid,
                        model = TRUE)
  ref <- nb(s$val, ms, cr_outcome, horizon = 5, thresholds = 0.2, cause = "death")$risk
  # the model's states named differently from the outcome's: matched by position
  dev$st2 <- factor(dev$cr_code, 0:2, c("cens", "rec", "dth"))
  ms2 <- survival::coxph(survival::Surv(cr_time, st2) ~ age + nodes, data = dev, id = pid,
                         model = TRUE)
  expect_equal(nb(s$val, ms2, cr_outcome, horizon = 5, thresholds = 0.2, cause = "death")$risk, ref)
  # with a single-event or binary outcome
  expect_error(suppress_matching(nb(s$val, ms, ~ survival::Surv(cr_time, cr_code > 0),
                                     horizon = 5), "^No patient treated"),
               "competing-risks model")
  expect_error(nb(s$val, ms, "death"), "use it with a (survival|competing-risks) outcome")
  # missing covariates
  v <- s$val; v$nodes <- NULL
  expect_error(nb(v, ms, cr_outcome, horizon = 5), "missing 1 covariate .*: nodes")
  # also when its data is no longer reachable and no model frame was kept
  fit <- function(d) survival::coxph(survival::Surv(cr_time, st) ~ age + nodes, data = d, id = pid)
  far <- fit(dev)
  expect_equal(nb(s$val, far, cr_outcome, horizon = 5, thresholds = 0.2, cause = "death")$risk, ref)
})

test_that("a Fine-Gray model for another cause is rejected (riskRegression would not notice)", {
  s <- rott_split()
  fg <- riskRegression::FGR(Hist(cr_time, cr_code) ~ age + nodes, data = s$dev,
                            cause = 1)
  expect_error(nb(s$val, fg, cr_outcome, horizon = 5, cause = "death"),
               "fitted for cause '1'.*cause '2'")
  # riskRegression itself returns the cause-1 risk when asked for cause 2
  utils::capture.output(r1 <- riskRegression::predictRisk(fg, s$val[1:5, ], times = 5, cause = 1))
  utils::capture.output(r2 <- riskRegression::predictRisk(fg, s$val[1:5, ], times = 5, cause = 2))
  expect_equal(r1, r2)
  X <- as.matrix(s$dev[, c("age", "nodes")])
  cr <- cmprsk::crr(s$dev$cr_time, s$dev$cr_code, cov1 = X, failcode = 2, cencode = 0)
  expect_error(nb(s$val, cr, cr_outcome, horizon = 5), "fitted for failcode 2")
})

test_that("a standard survival model with a competing outcome works, with a warning", {
  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(survival::Surv(cr_time, cr_code == 1) ~ age + size + nodes + grade,
                        data = s$dev)
  expect_warning(res <- withCallingHandlers(
    nb(v, cx, cr_outcome, horizon = 5, thresholds = surv_z),
    warning = function(w) if (grepl("^No patient treated", conditionMessage(w)))
      invokeRestart("muffleWarning")),
    "fitted to a standard survival outcome.*not the cumulative incidence")
  expect_same_as_dcurves_surv(res, v, "cr_time", "cr_status", 5)
})

test_that("a competing-risks model needs a competing-risks outcome", {
  s <- rott_split()
  csc <- riskRegression::CSC(cr_form, data = s$dev)
  fg <- riskRegression::FGR(Hist(cr_time, cr_code) ~ age + nodes, data = s$dev,
                            cause = 1)
  rf <- randomForestSRC::rfsrc(Surv(cr_time, cr_code) ~ age + nodes, data = s$dev, ntree = 20)
  X <- as.matrix(s$dev[, c("age", "nodes")])
  cr <- cmprsk::crr(s$dev$cr_time, s$dev$cr_code, cov1 = X, failcode = 1, cencode = 0)
  csi <- riskRegression::CSC(Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  for (m in list(csc, csi, fg, rf, cr)) {
    expect_error(nb(s$val, m, ~ survival::Surv(cr_time, cr_code > 0), horizon = 5),
                 "competing-risks model")
    expect_error(nb(s$val, m, "death"), "use it with a (survival|competing-risks) outcome")
  }
})

test_that("missing covariates are reported for competing-risks models", {
  s <- rott_split()
  v <- s$val; v$nodes <- NULL
  csc <- riskRegression::CSC(cr_form, data = s$dev)
  fg <- riskRegression::FGR(Hist(cr_time, cr_code) ~ age + nodes, data = s$dev,
                            cause = 1)
  rf <- randomForestSRC::rfsrc(Surv(cr_time, cr_code) ~ age + nodes, data = s$dev, ntree = 20)
  X <- as.matrix(s$dev[, c("age", "nodes")])
  cr <- cmprsk::crr(s$dev$cr_time, s$dev$cr_code, cov1 = X, failcode = 1, cencode = 0)
  for (m in list(csc, fg, rf, cr))
    expect_error(nb(v, m, cr_outcome, horizon = 5), "missing 1 covariate .*: nodes")
})

test_that("numeric competing-risks codes get a hint to use a factor status", {
  s <- rott_split()
  csc <- riskRegression::CSC(cr_form, data = s$dev)
  # Surv() turns status 2 into NA (with its own warning)
  expect_error(suppressWarnings(nb(s$val, csc, ~ survival::Surv(cr_time, cr_code),
                                   horizon = 5)),
               "missing times or statuses.*factor whose first level is censoring")
})
