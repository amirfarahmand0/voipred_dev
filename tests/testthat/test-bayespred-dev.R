# bayespred.dev Bayesian survival models (Cox and parametric AFT): they work in
# nb() without predict_function. The risk by the horizon is the posterior mean
# of risk from posterior(fit, time = horizon), as bayespred does for binary
# models.

skip_if_not_installed("bayespred.dev")

bp_formula <- survival::Surv(dtimey, death) ~ age + size + nodes + grade + hormon
bp_fit <- function(dev, family, ...)
  bayespred.dev::bpmfit(bp_formula, data = dev, family = family,
                        prior = bayespred.dev::flat(), ...)
bp_risk <- function(fit, nd, t, method = "pm")
  as.vector(predict(bayespred.dev::posterior(fit, time = t), nd, method = method))

check_bp <- function(model, s, own_risk, horizon = 5) {
  expect_no_warning(res <- suppress_matching(
    nb(model, s$val, surv_outcome, horizon = horizon, thresholds = surv_z),
    "^No patient treated"))
  expect_equal(attr(res, "risk"), as.vector(own_risk), tolerance = 1e-12)
  expect_same_as_dcurves_surv(res, s$val, "dtimey", "death", horizon)
  invisible(res)
}

test_that("Bayesian Cox (bpmfit, family = 'cox'): posterior mean of risk by the horizon", {
  s <- rott_split()
  fit <- bp_fit(s$dev, "cox")
  for (h in c(2, 5)) check_bp(fit, s, bp_risk(fit, s$val, h), horizon = h)
  # the horizon is used as intended: the plug-in risk equals a Breslow coxph
  cx <- survival::coxph(bp_formula, data = s$dev, ties = "breslow")
  expect_equal(bp_risk(fit, s$val, 5, method = "pe"), predict_risk(cx, s$val, horizon = 5),
               tolerance = 1e-8)
  # and the posterior mean is close to it, not equal
  pm <- attr(nb(fit, s$val, surv_outcome, horizon = 5, thresholds = 0.2), "risk")
  expect_lt(max(abs(pm - bp_risk(fit, s$val, 5, method = "pe"))), 0.05)
  expect_gt(max(abs(pm - bp_risk(fit, s$val, 5, method = "pe"))), 0)
})

test_that("Bayesian parametric survival (Weibull, log-normal, log-logistic, exponential)", {
  s <- rott_split()
  for (fam in c("weibull", "lognormal", "loglogistic", "exponential")) {
    fit <- bp_fit(s$dev, fam)
    check_bp(fit, s, bp_risk(fit, s$val, 5))
    sr <- survival::survreg(bp_formula, data = s$dev, dist = fam)
    expect_equal(bp_risk(fit, s$val, 5, method = "pe"), predict_risk(sr, s$val, horizon = 5),
                 tolerance = 1e-6, label = fam)
  }
})

test_that("a posterior packaged at the horizon (bpm) works at that horizon only", {
  s <- rott_split()
  fit <- bp_fit(s$dev, "cox")
  p5 <- bayespred.dev::posterior(fit, time = 5)
  check_bp(p5, s, bp_risk(fit, s$val, 5))
  expect_error(nb(p5, s$val, surv_outcome, horizon = 3),
               "gives risks at time = 5, but `horizon` is 3.*posterior\\(fit, time = 3\\)")
  expect_error(nb(p5, s$val, "death"), "packaged at time = 5.*horizon = 5")
  # without a horizon there is no baseline hazard, hence no risk
  expect_error(nb(bayespred.dev::posterior(fit), s$val, surv_outcome, horizon = 5),
               "no baseline hazard.*posterior\\(fit, time = 5\\)")
  # a fit kept without its data can only be used at the horizons stored at fit time
  fit_nd <- bp_fit(s$dev, "cox", model = FALSE, time = 5)
  check_bp(fit_nd, s, bp_risk(fit, s$val, 5))
  expect_error(nb(fit_nd, s$val, surv_outcome, horizon = 3), "model = FALSE|individual-level")
})

test_that("checks: binary outcome, missing covariates, predict_function gives the same", {
  s <- rott_split()
  fit <- bp_fit(s$dev, "weibull")
  expect_error(nb(fit, s$val, "death"), "survival outcome made with Surv\\(\\) and give `horizon`")
  v <- s$val; v$nodes <- NULL
  expect_error(nb(fit, v, surv_outcome, horizon = 5), "missing 1 covariate .*: nodes")
  expect_error(nb(bayespred.dev::posterior(fit, time = 5), v, surv_outcome, horizon = 5),
               "missing 1 covariate .*: nodes")
  a <- nb(fit, s$val, surv_outcome, horizon = 5, thresholds = c(0.1, 0.3))
  b <- nb(fit, s$val, surv_outcome, horizon = 5, thresholds = c(0.1, 0.3),
          predict_function = function(model, newdata, horizon, ...)
            predict(bayespred.dev::posterior(model, time = horizon), newdata))
  expect_equal(a$NB_model, b$NB_model)
  # the horizon beyond the model's own follow-up is flagged
  dev <- s$dev; cut <- dev$dtimey > 4
  dev$death[cut] <- 0; dev$dtimey[cut] <- 4
  w <- character()
  withCallingHandlers(nb(bp_fit(dev, "cox"), s$val, surv_outcome, horizon = 5, thresholds = 0.1),
                      warning = function(cond) { w <<- c(w, conditionMessage(cond))
                        invokeRestart("muffleWarning") })
  expect_true(any(grepl("beyond the last time in the bpmfit model's own data \\(4\\)", w)))
})

test_that("the binary bayespred.dev model still goes through the standard route", {
  g <- gusto_split()
  fit <- bayespred.dev::bpmfit(gusto_formula, data = g$dev, family = "binomial")
  expect_equal(predict_risk(fit, g$val), as.vector(predict(fit, g$val)))
  expect_equal(predict_risk(bayespred.dev::posterior(fit), g$val), as.vector(predict(fit, g$val)))
})

test_that("a Bayesian Cox model can be the censoring model of cIPCW", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val
  cens_f <- survival::Surv(dtimey, death == 0) ~ age + nodes + grade
  bcens <- bayespred.dev::bpmfit(cens_f, data = v, family = "cox", prior = bayespred.dev::flat())
  ccens <- survival::coxph(cens_f, data = v, ties = "breslow")
  a <- nb(cx, v, surv_outcome, horizon = 5, censoring = "cipcw",
          cipcw = cipcw_control(model = bcens))
  b <- nb(cx, v, surv_outcome, horizon = 5, censoring = "cipcw",
          cipcw = cipcw_control(model = ccens))
  # posterior mean of the censoring survival vs the plug-in: close
  expect_lt(max(abs(a$NB_model - b$NB_model), na.rm = TRUE), 2e-3)
  # the censoring probabilities at the horizon are the posterior means
  expect_equal(surv_at(bcens, v[1:20, ], rep(5, 20)),
               1 - bp_risk(bcens, v[1:20, ], 5), tolerance = 1e-12)
  # fitted to the event instead of censoring: stopped
  bev <- bayespred.dev::bpmfit(survival::Surv(dtimey, death) ~ age + nodes, data = v,
                               family = "cox", prior = bayespred.dev::flat())
  expect_error(nb(cx, v, surv_outcome, horizon = 5, censoring = "cipcw",
                  cipcw = cipcw_control(model = bev)), "fitted to the event")
})
