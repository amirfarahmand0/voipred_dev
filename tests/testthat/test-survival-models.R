# Survival models: nb() gets each model's predicted risk by the horizon
# (checked against an independent calculation) and the net benefit equals
# dcurves (Kaplan-Meier) on those risks. Models are grouped by the route that
# gives the risk: the standard route (riskRegression::predictRisk), or a
# predict_risk() method shipped by voipred for models it cannot handle.

local_base_formula()

check_model <- function(model, s, own_risk, horizon = 5, outcome = surv_outcome,
                        time = "dtimey", status = "death", ...) {
  # The only warning allowed: no follow-up to the horizon at high thresholds.
  expect_no_warning(res <- suppress_matching(
    nb(s$val, model, outcome, horizon = horizon, thresholds = surv_z, ...),
    "^No patient treated"))
  r <- res$risk
  expect_true(all(r >= 0 & r <= 1))
  expect_equal(r, as.vector(own_risk), tolerance = 1e-6)
  expect_same_as_dcurves_surv(res, s$val, time, status, horizon)
  invisible(res)
}

km_at <- function(sf, t) as.vector(summary(sf, times = t, extend = TRUE)$surv)
# reference for Cox models: survival::survfit on a copy that keeps its model frame
cox_ref <- function(m, nd, data, t = 5)
  1 - km_at(survival::survfit(stats::update(m, data = data, model = TRUE), newdata = nd), t)


# ---- standard route: riskRegression::predictRisk, no extra code ------------

test_that("coxph fitted with x = TRUE (Efron, Breslow, case weights, ties, strata, splines)", {
  s <- rott_split()
  dev <- s$dev
  cx <- survival::coxph(surv_formula, data = dev, x = TRUE)
  check_model(cx, s, cox_ref(cx, s$val, dev))
  br <- survival::coxph(surv_formula, data = dev, ties = "breslow", x = TRUE)
  check_model(br, s, cox_ref(br, s$val, dev))
  dev2 <- dev; dev2$dtimey <- round(dev2$dtimey, 1)
  set.seed(1); dev2$w <- runif(nrow(dev2), 0.5, 2)
  cw <- survival::coxph(surv_formula, data = dev2, weights = w, x = TRUE)
  check_model(cw, s, cox_ref(cw, s$val, dev2))
  st <- survival::coxph(survival::Surv(dtimey, death) ~ age + nodes + size + strata(grade),
                        data = dev, x = TRUE)
  check_model(st, s, cox_ref(st, s$val, dev))
  sp <- survival::coxph(survival::Surv(dtimey, death) ~ splines::ns(age, 3) * hormon +
                          size + nodes, data = dev, x = TRUE)
  check_model(sp, s, cox_ref(sp, s$val, dev))
})

test_that("rms::cph and rms::psm", {
  s <- rott_split()
  cp <- rms::cph(surv_formula, data = s$dev, x = TRUE, y = TRUE, surv = TRUE)
  check_model(cp, s, 1 - rms::survest(cp, newdata = s$val, times = 5)$surv)
  ps <- rms::psm(surv_formula, data = s$dev, dist = "weibull")
  check_model(ps, s, 1 - rms::survest(ps, newdata = s$val, times = 5)$surv)
})

test_that("survreg (Weibull, log-normal, log-logistic, exponential, strata)", {
  s <- rott_split()
  for (dist in c("weibull", "lognormal", "loglogistic", "exponential")) {
    m <- survival::survreg(surv_formula, data = s$dev, dist = dist)
    lp <- predict(m, newdata = s$val, type = "lp"); u <- (log(5) - lp) / m$scale
    own <- switch(dist, weibull = , exponential = 1 - exp(-exp(u)),
                  lognormal = pnorm(u), loglogistic = plogis(u))
    check_model(m, s, own)
  }
  # one scale per stratum
  m <- survival::survreg(survival::Surv(dtimey, death) ~ age + nodes + strata(grade),
                         data = s$dev)
  lp <- predict(m, newdata = s$val, type = "lp")
  sc <- m$scale[paste0("grade=", s$val$grade)]
  # the per-row scale, independently: survreg's own median, S(median) = 1/2
  q <- predict(m, newdata = s$val, type = "quantile", p = 0.5)
  expect_equal(unname((log(q) - lp) / log(log(2))), unname(sc), tolerance = 1e-8)
  check_model(m, s, 1 - exp(-exp((log(5) - lp) / sc)))
  # riskRegression gets the stratum scales wrong for some rows: voipred does not use it here
  expect_gt(max(abs(riskRegression::predictRisk(m, s$val, times = 5) -
                      (1 - exp(-exp((log(5) - lp) / sc))))), 0.01)
})

test_that("flexsurvreg and flexsurvspline (formula written in the call or stored in a variable)", {
  s <- rott_split()
  own <- function(m) 1 - predict(m, newdata = s$val, type = "survival", times = 5)$.pred_survival
  gg <- flexsurv::flexsurvreg(survival::Surv(dtimey, death) ~ age + size + nodes + grade,
                              data = s$dev, dist = "gengamma")
  check_model(gg, s, own(gg))
  sp <- flexsurv::flexsurvspline(survival::Surv(dtimey, death) ~ age + size + nodes + grade,
                                 data = s$dev, k = 2)
  check_model(sp, s, own(sp))
  # riskRegression fails with a formula variable: voipred's fallback takes over
  wb <- flexsurv::flexsurvreg(surv_formula, data = s$dev, dist = "weibull")
  expect_error(riskRegression::predictRisk(wb, s$val, times = 5))
  check_model(wb, s, own(wb))
})

test_that("random survival forests: ranger and randomForestSRC", {
  s <- rott_split()
  set.seed(1)
  rg <- ranger::ranger(surv_formula, data = s$dev, num.trees = 200)
  pr <- predict(rg, data = s$val)
  i <- findInterval(5, pr$unique.death.times)
  check_model(rg, s, 1 - pr$survival[, i])
  rf <- randomForestSRC::rfsrc(rsf_formula, data = s$dev, ntree = 200, seed = -1)
  pr <- predict(rf, newdata = s$val)
  i <- findInterval(5, rf$time.interest)
  check_model(rf, s, 1 - pr$survival[, i])
})

test_that("prodlim", {
  s <- rott_split()
  pl <- prodlim::prodlim(prodlim::Hist(dtimey, death) ~ grade, data = s$dev)
  km <- survival::survfit(survival::Surv(dtimey, death) ~ grade, data = s$dev)
  by_grade <- 1 - summary(km, times = 5)$surv
  names(by_grade) <- sub("grade=", "", summary(km, times = 5)$strata)
  check_model(pl, s, by_grade[as.character(s$val$grade)])
})


# ---- voipred's own methods, for models the standard route cannot handle -----

test_that("coxph without x = TRUE (riskRegression refuses it): own baseline hazard", {
  s <- rott_split()
  dev <- s$dev
  cx <- survival::coxph(surv_formula, data = dev)
  expect_error(riskRegression::predictRisk(cx, s$val, times = 5), "x=TRUE")
  check_model(cx, s, cox_ref(cx, s$val, dev))
  # identical to the standard route on the same model fitted with x = TRUE
  cxx <- survival::coxph(surv_formula, data = dev, x = TRUE)
  expect_equal(predict_risk(cx, s$val, horizon = 5), predict_risk(cxx, s$val, horizon = 5),
               tolerance = 1e-12)
  br <- survival::coxph(surv_formula, data = dev, ties = "breslow")
  check_model(br, s, cox_ref(br, s$val, dev))
  dev2 <- dev; dev2$dtimey <- round(dev2$dtimey, 1)
  set.seed(1); dev2$w <- runif(nrow(dev2), 0.5, 2)
  cw <- survival::coxph(surv_formula, data = dev2, weights = w)
  check_model(cw, s, cox_ref(cw, s$val, dev2))
  # without the stored outcome, or stratified: survfit (needs the model frame)
  ny <- survival::coxph(surv_formula, data = dev, y = FALSE, model = TRUE)
  check_model(ny, s, cox_ref(cx, s$val, dev))
  st <- survival::coxph(survival::Surv(dtimey, death) ~ age + nodes + size + strata(grade),
                        data = dev, model = TRUE)
  check_model(st, s, cox_ref(st, s$val, dev))
})

test_that("a coxph model that needs its (unavailable) data gives a clear error", {
  s <- rott_split()
  fit <- function(d) survival::coxph(surv_formula, data = d, y = FALSE)
  ny <- fit(s$dev)                                  # data not reachable afterwards
  expect_error(nb(s$val, ny, surv_outcome, horizon = 5),
               "Could not compute survival from this coxph model.*x = TRUE")
})

test_that("the own Cox baseline hazard equals survival::basehaz", {
  s <- rott_split()
  dev <- s$dev; dev$dtimey <- round(dev$dtimey, 1)
  set.seed(4); dev$w <- runif(nrow(dev), 0.5, 2)
  for (ties in c("efron", "breslow")) {
    cx <- survival::coxph(surv_formula, data = dev, weights = w, ties = ties, model = TRUE)
    ref <- survival::basehaz(cx, centered = TRUE)
    ref <- ref[!duplicated(ref$hazard), ]
    own <- cox_basehaz(cx)
    i <- findInterval(ref$time, own$time)
    expect_equal(ifelse(i == 0, 0, own$hazard[pmax(i, 1)]), ref$hazard, tolerance = 1e-10)
  }
})

test_that("aorsf oblique survival forest", {
  s <- rott_split()
  ao <- aorsf::orsf(s$dev, surv_formula, n_tree = 100)
  expect_error(riskRegression::predictRisk(ao, s$val, times = 5))
  check_model(ao, s, 1 - predict(ao, new_data = s$val, pred_type = "surv", pred_horizon = 5))
})

test_that("partykit::ctree survival tree", {
  s <- rott_split()
  ct <- partykit::ctree(surv_formula, data = s$dev)
  expect_error(riskRegression::predictRisk(ct, s$val, times = 5))
  sfs <- predict(ct, newdata = s$val, type = "prob")
  check_model(ct, s, 1 - vapply(sfs, km_at, numeric(1), t = 5))
})

test_that("rstpm2::stpm2 flexible parametric model (needs rstpm2 attached)", {
  s <- rott_split()
  st <- with_rstpm2(rstpm2::stpm2(
    survival::Surv(dtimey, death) ~ age + size + nodes + grade + hormon,
    data = s$dev, df = 3))
  nd <- s$val; nd$dtimey <- 5
  expect_error(nb(s$val, st, surv_outcome, horizon = 5), "run library\\(rstpm2\\)")
  with_rstpm2(
    check_model(st, s, 1 - rstpm2::predict(st, newdata = nd, type = "surv")))
})

test_that("mgcv::gam with the Cox PH family", {
  s <- rott_split()
  g <- mgcv::gam(dtimey ~ s(age) + size + s(nodes, k = 5) + grade + hormon,
                 family = mgcv::cox.ph(), weights = death, data = s$dev)
  nd <- s$val; nd$dtimey <- 5
  check_model(g, s, 1 - predict(g, newdata = nd, type = "response"))
  expect_error(nb(s$val, g, "death"), "Cox PH family\\) is a survival model")
})


# ---- models without a route, and the checks on the model's outcome ---------

test_that("models giving only relative risks or times have no route: clear error", {
  s <- rott_split()
  X <- model.matrix(~ age + nodes + grade + pgr + er + hormon + chemo, s$dev)[, -1]
  yb <- ifelse(s$dev$death == 1, s$dev$dtimey, -s$dev$dtimey)
  xg <- xgboost::xgb.train(params = list(objective = "survival:cox"),
                           data = xgboost::xgb.DMatrix(X, label = yb), nrounds = 5)
  msg <- "Could not get predicted risks from this .* model.*give `predict_function"
  expect_error(nb(s$val, xg, surv_outcome, horizon = 5), msg)
  gb <- suppressWarnings(gbm::gbm(surv_formula, data = s$dev, distribution = "coxph",
                                  n.trees = 50))
  expect_error(nb(s$val, gb, surv_outcome, horizon = 5), msg)
  y <- survival::Surv(s$dev$dtimey, s$dev$death)
  gn <- suppressWarnings(glmnet::glmnet(X, y, family = "cox", lambda = 0.02))
  expect_error(nb(s$val, gn, surv_outcome, horizon = 5), msg)
  set.seed(1)
  cg <- suppressWarnings(glmnet::cv.glmnet(X, y, family = "cox", nfolds = 3))
  expect_error(nb(s$val, cg, surv_outcome, horizon = 5), msg)
  # rpart survival trees: riskRegression returns relative rates, stopped by the check
  rp <- rpart::rpart(surv_formula, data = s$dev)
  expect_error(nb(s$val, rp, surv_outcome, horizon = 5), "outside \\[0, 1\\].*not probabilities")
  # with a predict_function returning their relative risks, the [0, 1] check stops them
  expect_error(nb(s$val, xg, surv_outcome, horizon = 5,
                  predict_function = function(model, newdata, ...)
                    predict(model, as.matrix(newdata[, colnames(X)]))),
               "outside \\[0, 1\\].*not probabilities")
})

test_that("a binary model with a survival outcome works, with a warning", {
  s <- rott_split()
  dev <- s$dev[s$dev$dtimey > 5 | s$dev$death == 1, ]
  dev$dead5 <- as.integer(dev$dtimey <= 5 & dev$death == 1)
  lg <- glm(dead5 ~ age + size + nodes + grade + hormon, family = binomial, data = dev)
  expect_warning(res <- suppress_matching(
    nb(s$val, lg, surv_outcome, horizon = 5, thresholds = surv_z), "^No patient treated"),
    "fitted to a binary outcome, not a survival outcome")
  expect_equal(res$risk, unname(predict(lg, s$val, type = "response")))
  expect_same_as_dcurves_surv(res, s$val, "dtimey", "death", 5)
})

test_that("survival models with a binary outcome are rejected", {
  s <- rott_split()
  models <- list(
    survival::coxph(surv_formula, data = s$dev),
    survival::coxph(surv_formula, data = s$dev, x = TRUE),
    survival::survreg(surv_formula, data = s$dev),
    flexsurv::flexsurvreg(surv_formula, data = s$dev, dist = "weibull"),
    randomForestSRC::rfsrc(rsf_formula, data = s$dev, ntree = 20),
    partykit::ctree(surv_formula, data = s$dev),
    aorsf::orsf(s$dev, surv_formula, n_tree = 5))
  for (m in models)
    expect_error(nb(s$val, m, "death", thresholds = 0.1),
                 "fitted to a survival outcome: use it with a survival outcome")
  # ranger keeps no readable formula: the standard binary routes fail instead
  rg <- ranger::ranger(surv_formula, data = s$dev, num.trees = 20)
  expect_error(nb(s$val, rg, "death", thresholds = 0.1), "Could not get predicted risks")
  # a model whose formula can no longer be found (fitted inside a function):
  # the Surv outcome it keeps is used
  fit_in <- function(d) { f <- surv_formula; flexsurv::flexsurvreg(f, data = d, dist = "weibull") }
  expect_error(nb(s$val, fit_in(s$dev), "death", thresholds = 0.1),
               "fitted to a survival outcome: use it with a survival outcome")
})

test_that("a horizon beyond the model's own follow-up gives a warning", {
  s <- rott_split()
  dev <- s$dev
  cut <- dev$dtimey > 4
  dev$death[cut] <- 0; dev$dtimey[cut] <- 4                 # follow-up ends at 4
  models <- list(
    survival::coxph(surv_formula, data = dev),
    survival::coxph(surv_formula, data = dev, x = TRUE),
    survival::survreg(surv_formula, data = dev),
    flexsurv::flexsurvreg(survival::Surv(dtimey, death) ~ age + nodes, data = dev,
                          dist = "weibull"),
    rms::cph(surv_formula, data = dev, x = TRUE, y = TRUE, surv = TRUE))
  for (m in models)
    expect_warning(suppress_matching(nb(s$val, m, surv_outcome, horizon = 5, thresholds = 0.1),
                                     "^No patient treated"),
                   "beyond the last time in the .* model's own data \\(4\\)")
  # within follow-up: no warning
  expect_no_warning(nb(s$val, models[[1]], surv_outcome, horizon = 3, thresholds = 0.1))
})

test_that("missing covariates are reported for survival models", {
  s <- rott_split()
  v <- s$val; v$nodes <- NULL; v$grade <- NULL
  named <- list(
    survival::coxph(surv_formula, data = s$dev),
    survival::coxph(surv_formula, data = s$dev, x = TRUE),
    survival::survreg(surv_formula, data = s$dev),
    flexsurv::flexsurvreg(surv_formula, data = s$dev, dist = "weibull"),
    rms::cph(surv_formula, data = s$dev, x = TRUE, y = TRUE, surv = TRUE),
    randomForestSRC::rfsrc(rsf_formula, data = s$dev, ntree = 20),
    aorsf::orsf(s$dev, surv_formula, n_tree = 10),
    partykit::ctree(surv_formula, data = s$dev),
    with_rstpm2(rstpm2::stpm2(
      survival::Surv(dtimey, death) ~ age + nodes + grade, data = s$dev, df = 3)),
    mgcv::gam(dtimey ~ age + nodes + grade, family = mgcv::cox.ph(),
              weights = death, data = s$dev))
  for (m in named)
    expect_error(nb(v, m, surv_outcome, horizon = 5, thresholds = 0.1),
                 "missing 2 covariates that the model was developed with: .*(grade|nodes)")
  # ranger: covariates cannot be read generically; ranger's own error is shown
  rg <- ranger::ranger(surv_formula, data = s$dev, num.trees = 20)
  expect_error(nb(v, rg, surv_outcome, horizon = 5, thresholds = 0.1),
               "Could not get predicted risks from this ranger model.*If a covariate is missing")
})

test_that("each model's risk grows with the horizon", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  r2 <- nb(s$val, cx, surv_outcome, horizon = 2, thresholds = 0.1)$risk
  r5 <- nb(s$val, cx, surv_outcome, horizon = 5, thresholds = 0.1)$risk
  expect_true(all(r5 >= r2))
  expect_gt(mean(r5 - r2), 0.05)
})
