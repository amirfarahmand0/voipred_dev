# The value-of-information functions refit the cIPCW censoring model in every
# draw (cipcw_control(), "Uncertainty").

refit_setup <- function() {
  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(surv_formula, data = s$dev)
  set.seed(5)
  w <- rpois(nrow(v), 1)
  list(s = s, v = v, cx = cx, z = c(0.1, 0.2, 0.3, 0.4), w = w,
       dup = v[rep(seq_len(nrow(v)), w), ],
       cf = survival::Surv(dtimey, death == 0) ~ age + nodes + grade)
}

cipcw_nb <- function(r, data, cens, ...)
  suppressWarnings(nb_val(data, r$cx, surv_outcome, horizon = 5, thresholds = r$z,
                          censoring = "cipcw", cipcw = cipcw_control(model = cens, ...)))

test_that("every known censoring model is refitted with the draw's weights", {
  r <- refit_setup()
  cf <- r$cf
  fits <- list(
    coxph = function(d) survival::coxph(cf, data = d, ties = "breslow"),
    survreg = function(d) survival::survreg(cf, data = d, dist = "weibull"),
    cph = function(d) rms::cph(cf, data = d, x = TRUE, y = TRUE, surv = TRUE,
                               method = "breslow"),
    psm = function(d) rms::psm(cf, data = d, dist = "weibull"),
    flexsurvreg = function(d) flexsurv::flexsurvreg(cf, data = d, dist = "weibull"),
    prodlim = function(d) prodlim::prodlim(prodlim::Hist(dtimey, death == 0) ~ 1, data = d))
  for (nm in names(fits)) {
    res <- cipcw_nb(r, r$v, fits[[nm]](r$v))
    # weights of 1: the refit is the model itself, so the draw is nb_val()
    one <- fixed_draws(res, W = matrix(1, nrow(r$v), 2))
    expect_equal(one$nb[1, , 1], res$estimates$NB_model, tolerance = 1e-12, label = nm)
    expect_equal(one$nb[1, , 2], res$estimates$NB_all, tolerance = 1e-12, label = nm)
    # whole-number weights: the same as nb_val() on the duplicated patients, with the
    # censoring model fitted to them; this fails if the weights were ignored
    ref <- cipcw_nb(r, r$dup, fits[[nm]](r$dup))
    got <- fixed_draws(res, W = cbind(r$w, r$w))
    expect_equal(got$nb[1, , 1], ref$estimates$NB_model, tolerance = 1e-8, label = nm)
    expect_equal(got$nb[2, , 2], ref$estimates$NB_all, tolerance = 1e-8, label = nm)
  }
})

test_that("other censoring models need refit_function", {
  r <- refit_setup()
  cf <- r$cf
  tree <- partykit::ctree(cf, data = r$v)
  res <- cipcw_nb(r, r$v, tree)
  # voipred does not guess how to give it case weights
  set.seed(1)
  expect_error(draws_t(res, n_sim = 3), "does not know how to fit a constparty model")
  expect_error(evpi(res, n_sim = 3), "Give `refit_function")
  expect_error(evsi(res, future_n = 100, n_sim = 3), "Give `refit_function")
  # with refit_function it works: whole-number weights give the same as the
  # model fitted to the duplicated patients
  fn <- function(data, weights)
    partykit::ctree(survival::Surv(dtimey, death == 0) ~ age + nodes + grade, data = data,
                    weights = weights)
  rf <- cipcw_nb(r, r$v, tree, refit_function = fn)
  ref <- cipcw_nb(r, r$dup, partykit::ctree(cf, data = r$dup))
  expect_equal(fixed_draws(rf, W = cbind(r$w, r$w))$nb[1, , 1], ref$estimates$NB_model,
               tolerance = 1e-12)
  set.seed(1)
  expect_no_error(draws_t(rf, n_sim = 3))
})

test_that("refit_function replaces the automatic refit", {
  r <- refit_setup()
  cf <- r$cf
  cens <- survival::coxph(cf, data = r$v)
  a <- cipcw_nb(r, r$v, cens)
  b <- cipcw_nb(r, r$v, cens, refit_function = function(data, weights)
    survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade, data = data,
                    weights = weights))
  set.seed(2)
  W <- matrix(rexp(nrow(r$v) * 5), ncol = 5)
  W[sample(length(W), 2000)] <- 0      # zero weights: left out of the refit
  expect_equal(fixed_draws(a, W = W)$nb, fixed_draws(b, W = W)$nb,
               tolerance = 1e-10)
  # the refit really changes the censoring weights: not the same as nb_val()'s weights
  v <- voi_info(a)
  fixed <- nb_core_cipcw(v$time, v$status, v$k, v$risk, v$thresholds, v$horizon,
                         v$weight_event, v$weight_horizon, W[, 1])
  expect_false(isTRUE(all.equal(fixed_draws(a, W = W)$nb[1, , 1], fixed$NB_model,
                                tolerance = 1e-8)))
  expect_output(print(cipcw_control(model = cens, refit_function = function(data, weights) 1)),
                "with refit_function")
})

test_that("models validated with the same censoring settings share one refit per draw", {
  r <- refit_setup()
  calls <- 0
  ctrl <- cipcw_control(model = survival::coxph(r$cf, data = r$v),
                        refit_function = function(data, weights) {
                          calls <<- calls + 1
                          survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes +
                                            grade, data = data, weights = weights)
                        })
  f2 <- survival::coxph(survival::Surv(dtimey, death) ~ age + nodes, data = r$s$dev)
  a <- nb_val(r$v, r$cx, surv_outcome, horizon = 5, thresholds = r$z, censoring = "cipcw", cipcw = ctrl)
  b <- nb_val(r$v, f2, surv_outcome, horizon = 5, thresholds = r$z, censoring = "cipcw", cipcw = ctrl)
  set.seed(1)
  draws_t(a, b, n_sim = 4)
  expect_equal(calls, 4)
})

test_that("EVPI, P(useful) and EVSI with the refit, also for competing risks", {
  r <- refit_setup()
  res <- cipcw_nb(r, r$v, survival::coxph(r$cf, data = r$v))
  set.seed(3)
  d <- draws_t(res, n_sim = 30)
  expect_true(all(evpi_from_draws(d)$estimates$EVPI >= 0))
  expect_true(all(p_useful_from_draws(d)$estimates$p_useful >= 0 & p_useful_from_draws(d)$estimates$p_useful <= 1))
  set.seed(3)
  o <- evsi(res, future_n = c(500, 5000), n_sim = 20)
  expect_true(all(o$estimates$EVSI <= o$estimates$EVPI + 1e-15))
  set.seed(3)
  ob <- evsi(res, future_n = 500, n_sim = 10)
  expect_true(all(is.finite(ob$estimates$EVSI)))
  # competing risks: only true censoring is the event of the censoring model
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = r$s$dev)
  cens_cr <- survival::coxph(survival::Surv(cr_time, cr_code == 0) ~ age + nodes, data = r$v)
  rc <- nb_val(r$v, csc, cr_outcome, horizon = 5, thresholds = r$z, cause = "death",
               censoring = "cipcw", cipcw = cipcw_control(model = cens_cr))
  dc <- fixed_draws(rc, W = matrix(1, nrow(r$v), 2))
  expect_equal(dc$nb[1, , 1], rc$estimates$NB_model, tolerance = 1e-10)
  set.seed(4)
  expect_true(all(evpi(rc, n_sim = 20)$estimates$EVPI >= 0))
})

test_that("given censoring weights work in nb_val() but not in the value of information", {
  r <- refit_setup()
  cens <- survival::coxph(r$cf, data = r$v)
  m <- cipcw_nb(r, r$v, cens)
  v <- voi_info(m)
  g <- suppressWarnings(nb_val(r$v, r$cx, surv_outcome, horizon = 5, thresholds = r$z,
                               censoring = "cipcw",
                               cipcw = cipcw_control(weight_event = v$weight_event,
                                                 weight_horizon = v$weight_horizon)))
  expect_equal(g$estimates$NB_model, m$estimates$NB_model)
  set.seed(1)
  expect_error(draws_t(g, n_sim = 2), "cannot be re-estimated in each draw")
  expect_error(evpi(g, n_sim = 2), "give the censoring model instead")
  expect_error(evsi(g, future_n = 10, n_sim = 2), "cannot be re-estimated")
  # Kaplan-Meier needs no censoring model
  km <- suppressWarnings(nb_val(r$v, r$cx, surv_outcome, horizon = 5, thresholds = r$z))
  set.seed(1)
  expect_no_error(draws_t(km, n_sim = 2))
})

test_that("refit settings are checked", {
  r <- refit_setup()
  cens <- survival::coxph(r$cf, data = r$v)
  expect_error(cipcw_control(model = cens, refit = TRUE), "There is no `refit` setting")
  expect_error(cipcw_control(weight_event = rep(1, 5), weight_horizon = rep(1, 5),
                             refit_function = function(data, weights) NULL),
               "needs the censoring `model`")
  expect_error(cipcw_control(model = cens, refit_function = "coxph"), "must be a function")
  expect_error(nb_val(r$v, r$cx, surv_outcome, horizon = 5,
                      cipcw = cipcw_control(refit_function = function(data, weights) NULL)),
               "needs the censoring `model`")
  bad <- cipcw_nb(r, r$v, cens, refit_function = function(data, weights) stop("no fit"))
  expect_error(draws_t(bad, n_sim = 2), "Could not refit the censoring model.*no fit")
})
