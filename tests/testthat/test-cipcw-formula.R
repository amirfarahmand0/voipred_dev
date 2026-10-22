# `cipcw` of nb_val(): the covariates of the censoring as a formula (a Cox
# model fitted by voipred, truncated at 0.02), a fitted censoring model, or
# cipcw_control(); no `censoring` argument.

local_base_formula()

fit_cox <- function(s) survival::coxph(surv_formula, data = s$dev)

# Cox censoring model fitted by hand, as voipred fits it
hand_cox <- function(f, data, weights = NULL)
  do.call(survival::coxph, c(list(formula = f, data = data),
                             if (!is.null(weights)) list(weights = weights)))

test_that("cipcw = ~ covariates is a Cox censoring model truncated at 0.02", {
  s <- rott_split()
  v <- s$val
  cx <- fit_cox(s)
  z <- c(0.1, 0.2, 0.3)
  res <- nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z,
                cipcw = ~ age + nodes + grade)
  expect_equal(res$outcome$censoring, "cipcw")
  ctrl <- res$cipcw$control
  expect_s3_class(ctrl$model, "coxph")
  expect_equal(ctrl$min_prob, 0.02)
  expect_equal(ctrl$formula, ~ age + nodes + grade)
  # the same as the model fitted by hand, with the same truncation
  cens <- hand_cox(survival::Surv(dtimey, death == 0) ~ age + nodes + grade, v)
  by_hand <- nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z,
                    cipcw = cipcw_control(model = cens, min_prob = 0.02))
  expect_equal(res$estimates, by_hand$estimates)
  expect_equal(res$cipcw$weight_event, by_hand$cipcw$weight_event)
  expect_equal(unname(coef(ctrl$model)), unname(coef(cens)))
  # a fitted censoring model given directly, and through cipcw_control()
  expect_equal(nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z,
                      cipcw = cens)$estimates, res$estimates)
  expect_equal(nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z,
                      cipcw = cipcw_control(~ age + nodes + grade))$estimates,
               res$estimates)
  # a two-sided formula is fitted as it is
  expect_equal(nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z,
                      cipcw = survival::Surv(dtimey, death == 0) ~ age + nodes + grade)$estimates,
               res$estimates)
  # the data are not changed
  expect_false(any(grepl("^\\.voipred", names(v))))
})

test_that("the truncation of the formula route: default 0.02, other values, none", {
  # steep censoring in x2: one used weight is above 50 (55) without truncation
  d <- sim_surv(3000, informative = TRUE, seed = 3)
  set.seed(9)
  t_c <- rexp(nrow(d), 0.3 * exp(4 * d$x2))
  d$time <- pmin(d$t_ev, t_c)
  d$status <- ifelse(d$t_ev <= t_c, d$cause, 0L)
  p <- 1 - exp(-0.15 * exp(d$lp1) * 5^1.3)
  out <- ~ survival::Surv(time, status)
  used <- function(r) c(r$cipcw$weight_event[d$status != 0 & d$time <= 5],
                        r$cipcw$weight_horizon[d$time > 5])
  run <- function(cipcw) nb_val(d, p, out, horizon = 5, thresholds = c(0.1, 0.3),
                                cipcw = cipcw)
  expect_no_warning(r_off <- run(cipcw_control(~ x2, min_prob = NULL)))
  expect_gt(max(used(r_off)), 50)
  expect_warning(r02 <- run(~ x2), "^1 censoring probability was truncated at min_prob = 0.02")
  expect_equal(max(used(r02)), 50)
  expect_equal(sum(used(r02) != used(r_off)), 1)
  expect_warning(r05 <- run(cipcw_control(~ x2, min_prob = 0.05)),
                 "censoring probabilities were truncated at min_prob = 0.05")
  expect_equal(max(used(r05)), 20)
  expect_equal(used(r05), pmin(used(r_off), 20))
})

test_that("the formula route corrects informative censoring (closer to the truth than KM)", {
  d <- sim_surv(6000, informative = TRUE, seed = 7)
  p <- 1 - exp(-0.15 * exp(d$lp1) * 5^1.3)
  z <- c(0.2, 0.4, 0.6)
  out <- ~ survival::Surv(time, status)
  truth <- true_nb(d, p, z, 5)
  km <- nb_val(d, p, out, horizon = 5, thresholds = z)
  ci <- suppressWarnings(nb_val(d, p, out, horizon = 5, thresholds = z, cipcw = ~ x2))
  err_km <- max(abs(km$estimates$NB_model - truth))
  err_ci <- max(abs(ci$estimates$NB_model - truth))
  expect_lt(err_ci, err_km)
})

test_that("competing risks: only true censoring is censoring", {
  s <- rott_split()
  v <- s$val
  v$p <- pmin(pmax(v$age / 100, 0.01), 0.99)
  z <- c(0.3, 0.5)
  res <- nb_val(v, "p", cr_outcome, horizon = 5, thresholds = z, cipcw = ~ age + nodes)
  cens <- hand_cox(survival::Surv(cr_time, cr_status == "censored") ~ age + nodes, v)
  expect_equal(res$estimates,
               nb_val(v, "p", cr_outcome, horizon = 5, thresholds = z,
                      cipcw = cens)$estimates)
  # a two-sided formula fitted to the cause of interest is caught
  expect_error(nb_val(v, "p", cr_outcome, horizon = 5, thresholds = z,
                      cipcw = survival::Surv(cr_time, cr_status == "recurrence") ~ age),
               "fitted to the cause of interest")
})

test_that("case weights: the Cox censoring model is fitted with them", {
  s <- rott_split()
  v <- s$val
  cx <- fit_cox(s)
  set.seed(2)
  w <- runif(nrow(v), 0.5, 2)
  z <- c(0.1, 0.3)
  res <- nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z, weights = w,
                cipcw = ~ age + nodes)
  cens <- hand_cox(survival::Surv(dtimey, death == 0) ~ age + nodes, v, weights = w)
  expect_equal(unname(coef(res$cipcw$control$model)), unname(coef(cens)))
  expect_equal(res$estimates,
               nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z, weights = w,
                      cipcw = cens)$estimates)
})

test_that("the formula route refits its Cox model in every draw", {
  s <- rott_split()
  v <- s$val
  cx <- fit_cox(s)
  z <- c(0.1, 0.3)
  res <- nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z, cipcw = ~ age + nodes)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = v)
  by_hand <- nb_val(v, cx, surv_outcome, horizon = 5, thresholds = z, cipcw = cens)
  # the same draws and refits as the model fitted by hand
  set.seed(4); e1 <- evpi(res, n_sim = 15)
  set.seed(4); e2 <- evpi(by_hand, n_sim = 15)
  expect_equal(e1$estimates, e2$estimates)
  expect_equal(e1$p_useful$estimates, e2$p_useful$estimates)
  # the refits do change the weights from draw to draw
  v_info <- voi_info(res)
  refit <- censoring_refitter(v_info)
  set.seed(5)
  w1 <- refit(bb_weights(nrow(v)))
  expect_false(isTRUE(all.equal(w1$weight_horizon, res$cipcw$weight_horizon)))
  set.seed(6); o <- evsi(res, future_n = 500, n_sim = 10)
  expect_s3_class(o, "voipred_evsi")
})

test_that("several models share one censoring model; results combine", {
  s <- rott_split()
  v <- s$val
  cx <- fit_cox(s)
  cx2 <- survival::coxph(survival::Surv(dtimey, death) ~ age, data = s$dev)
  v$p1 <- predict_risk(cx, v, horizon = 5)
  v$p2 <- predict_risk(cx2, v, horizon = 5)
  z <- c(0.1, 0.3)
  both <- nb_val(v, list(Full = "p1", Age = "p2"), surv_outcome, horizon = 5,
                 thresholds = z, cipcw = ~ age + nodes)
  a <- nb_val(v, "p1", surv_outcome, horizon = 5, thresholds = z, cipcw = ~ age + nodes,
              label = "Full")
  b <- nb_val(v, "p2", surv_outcome, horizon = 5, thresholds = z, cipcw = ~ age + nodes,
              label = "Age")
  expect_equal(c(a, b)$estimates, both$estimates)
  set.seed(8); e1 <- evpi(both, n_sim = 10)
  set.seed(8); e2 <- evpi(a, b, n_sim = 10)
  expect_equal(e1$estimates, e2$estimates)
  # Kaplan-Meier and cIPCW cannot be combined
  km <- nb_val(v, "p2", surv_outcome, horizon = 5, thresholds = z, label = "Age")
  expect_error(c(a, km), "different censoring method")
  # different censoring covariates give different weights
  other <- nb_val(v, "p2", surv_outcome, horizon = 5, thresholds = z, cipcw = ~ age,
                  label = "Age")
  expect_error(c(a, other), "cannot be used together")
})

test_that("errors of the formula route", {
  s <- rott_split()
  v <- s$val
  cx <- fit_cox(s)
  run <- function(cipcw, data = v) nb_val(data, cx, surv_outcome, horizon = 5,
                                          thresholds = 0.2, cipcw = cipcw)
  expect_error(run(~ age + nope), "uses a variable that is not in `data`: nope")
  expect_error(run(~ nope1 + nope2), "variables that are not in `data`: nope1, nope2")
  v_na <- v
  v_na$rtime[2] <- NA                  # used by the censoring model only
  expect_error(run(~ age + rtime, v_na), "missing values \\(rtime\\)")
  v_clash <- v
  v_clash$.voipred_time <- 1
  expect_error(run(~ age, v_clash), "has a column named .voipred_time")
  expect_error(run(~ no_such_function(age)), "Could not fit the Cox censoring model")
  expect_error(run("age"), "must be a formula")
  expect_error(run(1:3), "must be a formula")
  # fitted to the event instead of censoring
  expect_error(run(survival::Surv(dtimey, death) ~ age), "fitted to the event")
  # binary outcome
  g <- gusto_split()
  f <- glm(day30 ~ age, family = binomial, data = g$dev)
  expect_error(nb_val(g$val, f, "day30", cipcw = ~ age), "only used with a survival outcome")
})

test_that("printing shows the censoring model", {
  s <- rott_split()
  res <- nb_val(s$val, fit_cox(s), surv_outcome, horizon = 5, thresholds = 0.2,
                cipcw = ~ age + nodes)
  expect_output(print(res), "censoring cIPCW \\(coxph\\)")
  expect_output(print(res$cipcw$control),
                "Cox censoring model ~age \\+ nodes, fitted by voipred, truncation: min_prob = 0.02")
})

test_that("any censoring model works, e.g. a random survival forest with refit_function", {
  skip_if_not_installed("ranger")
  s <- rott_split()
  v <- s$val
  v$censored <- 1 - v$death
  v$p <- pmin(0.9, 0.1 + v$nodes / 40)
  z <- c(0.1, 0.3)
  rf <- ranger::ranger(survival::Surv(dtimey, censored) ~ age + nodes + grade, data = v,
                       num.trees = 50, seed = 1)
  # the formula is read from the call, so the censoring outcome is checked (no warning)
  expect_equal(deparse(model_formula(rf)), "survival::Surv(dtimey, censored) ~ age + nodes + grade")
  expect_no_warning(res <- nb_val(v, "p", surv_outcome, horizon = 5, thresholds = z, cipcw = rf))
  expect_equal(res$outcome$censoring, "cipcw")
  expect_equal(nb_val(v, "p", surv_outcome, horizon = 5, thresholds = z,
                      cipcw = cipcw_control(model = rf))$estimates, res$estimates)
  # a forest fitted to the event is caught
  rf_event <- ranger::ranger(survival::Surv(dtimey, death) ~ age + nodes, data = v,
                             num.trees = 20, seed = 1)
  expect_error(nb_val(v, "p", surv_outcome, horizon = 5, thresholds = z, cipcw = rf_event),
               "fitted to the event, not to censoring")
  # the value-of-information functions need a way to refit it ...
  expect_error(evpi(res, n_sim = 2), "Give `refit_function")
  # ... and use refit_function, once per draw, with the draw's weights
  calls <- 0
  refit <- function(data, weights) {
    calls <<- calls + 1
    expect_equal(nrow(data), length(weights))
    ranger::ranger(survival::Surv(dtimey, censored) ~ age + nodes + grade, data = data,
                   case.weights = weights, num.trees = 50, seed = 1)
  }
  res_r <- nb_val(v, "p", surv_outcome, horizon = 5, thresholds = z,
                  cipcw = cipcw_control(model = rf, refit_function = refit))
  set.seed(1)
  e <- evpi(res_r, n_sim = 3)
  expect_equal(calls, 3)
  expect_s3_class(e, "voipred_evpi")
})

test_that("a parametric censoring model is refitted automatically", {
  skip_if_not_installed("flexsurv")
  s <- rott_split()
  v <- s$val
  v$p <- pmin(0.9, 0.1 + v$nodes / 40)
  fs <- flexsurv::flexsurvreg(survival::Surv(dtimey, death == 0) ~ age + nodes, data = v,
                              dist = "weibull")
  res <- nb_val(v, "p", surv_outcome, horizon = 5, thresholds = c(0.1, 0.3), cipcw = fs)
  set.seed(2)
  expect_s3_class(evpi(res, n_sim = 3), "voipred_evpi")
})

test_that("the formula of a model is read from its call when nothing else gives it", {
  m <- list(call = quote(fit(y ~ x + z, data = d)))
  expect_equal(deparse(call_formula(m)), "y ~ x + z")
  m2 <- list(call = quote(fit(formula = Surv(t, s) ~ x, data = d)))
  expect_equal(deparse(call_formula(m2)), "Surv(t, s) ~ x")
  expect_null(call_formula(list(call = quote(fit(d, "y")))))
  expect_null(call_formula(list()))
  expect_null(call_formula(1:3))
})
