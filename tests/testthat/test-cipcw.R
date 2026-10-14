# Conditional IPCW net benefit (censoring = "cipcw").

local_base_formula()

fork <- readRDS(test_path("fixtures", "cipcw_fork.rds"))

cipcw_nb <- function(model, data, outcome, ctrl, horizon = 5, thresholds = surv_z, ...)
  nb(model, data, outcome, horizon = horizon, thresholds = thresholds,
     censoring = "cipcw", cipcw = ctrl, ...)

# Muffle the expected "not a survival model" warning of fixed_model().
no_fixed_warning <- function(expr) withCallingHandlers(expr, warning = function(w)
  if (grepl("not a survival model", conditionMessage(w))) invokeRestart("muffleWarning"))

# Direct, one-threshold-at-a-time version of the estimator (the formulas of
# the dcurves cIPCW branch, with risk >= z and risk 0 never treated).
reference_cipcw <- function(time, status, k, p, z, tau, we, wh, w = rep(1, length(time))) {
  ev <- status == k & time <= tau
  comp <- status != 0 & status != k & time <= tau
  free <- time > tau
  prev <- sum(w * ev * we) / sum(w)
  out <- sapply(z, function(t) {
    tr <- p >= t & p > 0
    se <- if (sum(w * ev * we) == 0) 0 else sum(w * (tr & ev) * we) / sum(w * ev * we)
    negw <- w * (free * wh + comp * we)
    sp <- if (sum(negw) == 0) 1 else sum((!tr) * negw) / sum(negw)
    prev * se - (1 - prev) * (1 - sp) * t / (1 - t)
  })
  list(NB_model = out, NB_all = prev - (1 - prev) * z / (1 - z), prevalence = prev)
}

test_that("same as the dcurves cIPCW branch: Cox censoring model", {
  s <- rott_split()
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade,
                          data = s$val)
  res <- no_fixed_warning(cipcw_nb(fixed_model(fork$surv_cox$risk), s$val, surv_outcome,
                                   cipcw_control(model = cens), thresholds = fork$z))
  expect_equal(res$NB_model, fork$surv_cox$NB_model, tolerance = 1e-10)
  expect_equal(res$NB_all, fork$surv_cox$NB_all, tolerance = 1e-10)
  expect_equal(attr(res, "prevalence"), fork$surv_cox$prevalence, tolerance = 1e-10)
  expect_identical(attr(res, "censoring"), "cipcw")
})

test_that("same as the dcurves cIPCW branch: given weights", {
  s <- rott_split()
  f <- fork$surv_weights
  res <- no_fixed_warning(cipcw_nb(fixed_model(f$risk), s$val, surv_outcome,
                                   cipcw_control(weight_event = unname(f$weight_event),
                                                 weight_horizon = unname(f$weight_horizon)),
                                   thresholds = fork$z))
  expect_equal(res$NB_model, f$NB_model, tolerance = 1e-10)
  expect_equal(res$NB_all, f$NB_all, tolerance = 1e-10)
})

test_that("same as the dcurves cIPCW branch: competing risks, Cox censoring model", {
  s <- rott_split()
  cens <- survival::coxph(survival::Surv(cr_time, cr_code == 0) ~ age + nodes + grade,
                          data = s$val)
  res <- no_fixed_warning(cipcw_nb(fixed_model(fork$cr_cox$risk), s$val, cr_outcome,
                                   cipcw_control(model = cens), thresholds = fork$z))
  expect_equal(res$NB_model, fork$cr_cox$NB_model, tolerance = 1e-10)
  expect_equal(res$NB_all, fork$cr_cox$NB_all, tolerance = 1e-10)
  expect_equal(attr(res, "prevalence"), fork$cr_cox$prevalence, tolerance = 1e-10)
})

test_that("the vectorized estimator equals a direct implementation (ties, weights, causes)", {
  set.seed(8)
  n <- 400
  time <- round(rexp(n, 0.25), 1)
  status <- sample(0:2, n, replace = TRUE, prob = c(0.35, 0.45, 0.2))
  p <- round(runif(n), 2); p[1:15] <- 0                 # ties in risk, some risk 0
  we <- runif(n, 1, 4); wh <- runif(n, 1, 3); w <- rpois(n, 1.3)
  z <- c(0, 0.01, seq(0.05, 0.95, 0.05), p[20])         # a threshold equal to a risk
  for (k in 1:2) for (tau in c(2, 4)) {
    a <- nb_core_cipcw(time, status, k, p, z, tau, we, wh, w)
    b <- reference_cipcw(time, status, k, p, z, tau, we, wh, w)
    expect_equal(a$NB_model, b$NB_model, tolerance = 1e-12)
    expect_equal(a$NB_all, b$NB_all, tolerance = 1e-12)
    expect_equal(a$prevalence, b$prevalence, tolerance = 1e-12)
  }
})

test_that("a censoring model gives the same result as its weights given by hand", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade + size,
                          data = s$val, model = TRUE)
  for (ties in c("left", "right")) {
    a <- cipcw_nb(cx, s$val, surv_outcome, cipcw_control(model = cens, ties = ties))
    sf <- survival::survfit(cens, newdata = s$val)          # one curve per patient
    G <- function(t, i) {
      j <- if (ties == "left") findInterval(t, sf$time, left.open = TRUE)
           else findInterval(t, sf$time)
      if (j == 0) 1 else sf$surv[j, i]
    }
    n <- nrow(s$val)
    we <- 1 / vapply(seq_len(n), function(i) G(s$val$dtimey[i], i), 1)
    wh <- 1 / vapply(seq_len(n), function(i) G(5, i), 1)
    b <- cipcw_nb(cx, s$val, surv_outcome, cipcw_control(weight_event = we,
                                                         weight_horizon = wh))
    expect_equal(a$NB_model, b$NB_model, tolerance = 1e-10)
    expect_equal(a$NB_all, b$NB_all, tolerance = 1e-10)
  }
})

test_that("with Kaplan-Meier censoring weights, the risk by the horizon is 1 - KM", {
  # IPCW with the marginal KM of censoring reproduces the KM estimate exactly.
  # (Exact without ties between event and censoring times: at a tie, the
  # reverse KM counts the patient with the event as still at risk of censoring.)
  s <- rott_split()
  set.seed(2); s$val$dtimey <- s$val$dtimey + runif(nrow(s$val), 0, 1e-4)
  cx <- survival::coxph(surv_formula, data = s$dev)
  km_cens <- prodlim::prodlim(prodlim::Hist(dtimey, death == 0) ~ 1, data = s$val)
  for (tau in c(2, 5, 8)) {
    a <- suppressWarnings(cipcw_nb(cx, s$val, surv_outcome,
                                   cipcw_control(model = km_cens), horizon = tau))
    b <- quiet_nb(cx, s$val, surv_outcome, horizon = tau)
    expect_equal(attr(a, "prevalence"), attr(b, "prevalence"), tolerance = 1e-10)
    expect_equal(a$NB_all, b$NB_all, tolerance = 1e-10)
  }
})

test_that("the censoring model is checked against the outcome", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val
  run <- function(m, outcome = surv_outcome, model = cx, ...)
    cipcw_nb(model, v, outcome, cipcw_control(model = m), thresholds = 0.1, ...)
  ev <- survival::coxph(survival::Surv(dtimey, death) ~ age + nodes, data = v)
  expect_error(run(ev), "fitted to the event, not to censoring")
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  c1 <- survival::coxph(survival::Surv(cr_time, cr_code == 1) ~ age + nodes, data = v)
  c2 <- survival::coxph(survival::Surv(cr_time, cr_code == 2) ~ age + nodes, data = v)
  ca <- survival::coxph(survival::Surv(cr_time, cr_code > 0) ~ age + nodes, data = v)
  expect_error(run(c1, cr_outcome, csc), "fitted to the cause of interest \\(recurrence\\)")
  expect_error(run(c2, cr_outcome, csc, cause = "death"),
               "fitted to the cause of interest \\(death\\)")
  expect_error(run(c2, cr_outcome, csc), "fitted to a competing cause \\(death\\)")
  expect_error(run(ca, cr_outcome, csc), "fitted to all events together")
  # the same checks for other model classes that store their outcome
  sr <- survival::survreg(survival::Surv(dtimey, death) ~ age + nodes, data = v)
  expect_error(run(sr), "fitted to the event")
  fs <- flexsurv::flexsurvreg(survival::Surv(dtimey, death) ~ age + nodes, data = v,
                              dist = "weibull")
  expect_error(run(fs), "fitted to the event")
  rf <- randomForestSRC::rfsrc(Surv(dtimey, death) ~ age + nodes, data = v, ntree = 10)
  expect_error(run(rf), "fitted to the event")
  ct <- partykit::ctree(survival::Surv(dtimey, death) ~ age + nodes, data = v)
  expect_error(run(ct), "fitted to the event")
  ao <- aorsf::orsf(v, survival::Surv(dtimey, death) ~ age + nodes, n_tree = 5)
  expect_error(run(ao), "fitted to the event")
  # a censoring model fitted to other data, or one that does not store its outcome
  other <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = s$dev)
  expect_warning(run(other), "fitted to 1500 patients, but `data` has 1482")
  rg <- ranger::ranger(survival::Surv(dtimey, 1 - death) ~ age + nodes, data = v,
                       num.trees = 20)
  expect_warning(run(rg), "Could not verify that the censoring model")
  set.seed(1); v$odd <- rbinom(nrow(v), 1, 0.5)
  od <- survival::coxph(survival::Surv(dtimey, odd) ~ age + nodes, data = v)
  expect_warning(run(od), "its outcome does not match")
  # fitted correctly: no warning
  ok <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = v)
  expect_no_warning(run(ok))
  okc <- survival::coxph(survival::Surv(cr_time, cr_code == 0) ~ age + nodes, data = v)
  expect_no_warning(run(okc, cr_outcome, csc))
})

test_that("the censoring model's covariates must be in the data", {
  s <- rott_split()
  cx <- survival::coxph(survival::Surv(dtimey, death) ~ age + size, data = s$dev)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = s$val)
  v <- s$val; v$nodes <- NULL
  expect_error(cipcw_nb(cx, v, surv_outcome, cipcw_control(model = cens)),
               "missing 1 covariate that the censoring model was fitted with: nodes")
})

test_that("small censoring probabilities: warning, optional truncation, error at 0", {
  d <- sim_surv(3000, informative = TRUE, seed = 3)
  d$x3 <- d$x2
  p <- 1 - exp(-0.15 * exp(d$lp1) * 5^1.3)
  m <- fixed_model(p)
  steep <- fixed_censoring(function(x, t) exp(-0.05 * t * exp(2.5 * x)), vars = "x2")
  run <- function(ctrl) suppress_matching(
    cipcw_nb(m, d, ~ survival::Surv(time, status), ctrl, thresholds = c(0.1, 0.3)),
    "not a survival model|Could not verify")
  expect_warning(run(cipcw_control(model = steep)),
                 "censoring probabilities are below 0.01")
  expect_warning(r_tr <- run(cipcw_control(model = steep, min_prob = 0.05)),
                 "censoring probabilities were truncated at min_prob = 0.05")
  # truncated model weights = the same weights truncated by hand
  g_ev <- ifelse(d$status != 0 & d$time <= 5,
                 pmax(steep$g(d$x2, d$time - 1e-9 * max(d$time)), 0.05), 1)
  g_h <- ifelse(d$time > 5, pmax(steep$g(d$x2, 5), 0.05), 1)
  r_hand <- run(cipcw_control(weight_event = 1 / g_ev, weight_horizon = 1 / g_h))
  expect_equal(r_tr$NB_model, r_hand$NB_model, tolerance = 1e-6)
  # the warning counts only the probabilities the estimator uses
  used <- c(steep$g(d$x2, d$time - 1e-9 * max(d$time))[d$status != 0 & d$time <= 5],
            steep$g(d$x2, 5)[d$time > 5])
  expect_warning(run(cipcw_control(model = steep, min_prob = 0.05)),
                 paste0("^", sum(used < 0.05), " censoring probabilities"))
  # a censoring probability of exactly 0
  zero <- fixed_censoring(function(x, t) ifelse(x > 2, 0, exp(-0.1 * t)), vars = "x2")
  expect_error(run(cipcw_control(model = zero)), "exactly 0, which gives an infinite weight")
  expect_no_error(suppressWarnings(run(cipcw_control(model = zero, min_prob = 0.01))))
  # invalid probabilities from a model
  bad <- fixed_censoring(function(x, t) rep(NA_real_, length(x)), vars = "x2")
  expect_error(run(cipcw_control(model = bad)), "invalid censoring probabilities")
  # given weights: warning above 100, truncation at 1 / min_prob
  we <- rep(2, nrow(d)); wh <- rep(2, nrow(d))
  we[d$status != 0 & d$time <= 5][1:3] <- 500
  expect_warning(run(cipcw_control(weight_event = we, weight_horizon = wh)),
                 "3 cIPCW weights are above 100")
  expect_warning(run(cipcw_control(weight_event = we, weight_horizon = wh, min_prob = 0.02)),
                 "3 cIPCW weights were truncated at 1 / min_prob = 50")
  wu <- we; wu[d$time > 5 & d$status == 0][1] <- 1e4      # unused weights are ignored
  wu[d$time <= 5 & d$status == 0][1] <- 1e4
  expect_warning(run(cipcw_control(weight_event = wu, weight_horizon = wh)),
                 "3 cIPCW weights are above 100")
})

test_that("ties: left uses the censoring survival just before each patient's time", {
  # a censoring jump at exactly an event time: excluded with "left", included with "right"
  d <- data.frame(time = c(1, 2, 2, 3, 4, 6, 7, 8), status = c(1, 0, 1, 0, 1, 0, 0, 0),
                  x = 0)
  cens <- survival::coxph(survival::Surv(time, status == 0) ~ 1, data = d)
  left <- surv_at(cens, d, d$time, "left")
  right <- surv_at(cens, d, d$time, "right")
  bh <- survival::basehaz(cens)
  G <- function(t, left_open) {
    i <- findInterval(t, bh$time, left.open = left_open)
    exp(-ifelse(i == 0, 0, bh$hazard[pmax(i, 1)]))
  }
  expect_equal(left, G(d$time, TRUE))
  expect_equal(right, G(d$time, FALSE))
  expect_true(left[3] > right[3])                  # time 2: a censoring jump at 2
})

test_that("every censoring model class runs and agrees with its own predictions", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val
  cens_formula <- survival::Surv(dtimey, death == 0) ~ age + nodes + grade
  fits <- list(
    coxph = survival::coxph(cens_formula, data = v),
    cph = rms::cph(cens_formula, data = v, x = TRUE, y = TRUE, surv = TRUE),
    survreg = survival::survreg(cens_formula, data = v),
    flexsurv = flexsurv::flexsurvreg(cens_formula, data = v, dist = "weibull"),
    ranger = ranger::ranger(cens_formula, data = v, num.trees = 50, seed = 1),
    rfsrc = randomForestSRC::rfsrc(Surv(dtimey, cens) ~ age + nodes + grade,
                                   data = transform(v, cens = 1 - death), ntree = 50),
    aorsf = aorsf::orsf(v, cens_formula, n_tree = 20),
    ctree = partykit::ctree(cens_formula, data = v),
    stpm2 = with_rstpm2(rstpm2::stpm2(cens_formula, data = v, df = 3)),
    gam = mgcv::gam(dtimey ~ age + nodes + grade, family = mgcv::cox.ph(),
                    weights = 1 - death, data = v))
  with_rstpm2({
    for (nm in names(fits)) {
      res <- suppress_matching(cipcw_nb(cx, v, surv_outcome, cipcw_control(model = fits[[nm]])),
                               "Could not verify|below 0.01|No patient treated|fitted to .* patients")
      expect_true(all(is.finite(res$NB_all)), label = nm)
      expect_true(attr(res, "prevalence") > 0.2 && attr(res, "prevalence") < 0.32, label = nm)
    }
    # the censoring probabilities at the horizon equal each model's own prediction
    nd <- v[1:20, ]
    own <- list(
      survreg = 1 - survival::psurvreg(5, predict(fits$survreg, nd, type = "lp"),
                                       fits$survreg$scale, fits$survreg$dist),
      flexsurv = predict(fits$flexsurv, nd, type = "survival", times = 5)$.pred_survival,
      ranger = 1 - riskRegression::predictRisk(fits$ranger, nd, times = 5),
      rfsrc = 1 - riskRegression::predictRisk(fits$rfsrc, transform(nd, cens = 1 - death),
                                              times = 5),
      aorsf = as.vector(predict(fits$aorsf, new_data = nd, pred_type = "surv", pred_horizon = 5)),
      cph = 1 - riskRegression::predictRisk(fits$cph, nd, times = 5),
      stpm2 = rstpm2::predict(fits$stpm2, newdata = transform(nd, dtimey = 5), type = "surv"),
      gam = predict(fits$gam, newdata = transform(nd, dtimey = 5), type = "response"))
    for (nm in names(own))
      expect_equal(surv_at(fits[[nm]], nd, rep(5, 20), "right"),
                   as.vector(own[[nm]]), tolerance = 1e-6, label = nm)
  })
})

test_that("cIPCW is close to the truth under informative censoring; KM is not", {
  tau <- 5; z <- c(0.1, 0.2, 0.3, 0.4)
  d <- sim_surv(20000, informative = TRUE, seed = 11)
  p <- 1 - exp(-0.15 * exp(d$lp1) * tau^1.3)
  truth <- true_nb(d, p, z, tau)
  cens <- survival::coxph(survival::Surv(time, status == 0) ~ x2, data = d)
  ci <- suppress_matching(cipcw_nb(fixed_model(p), d, ~ survival::Surv(time, status),
                                   cipcw_control(model = cens), horizon = tau, thresholds = z),
                          "not a survival model|below 0.01")
  km <- quiet_nb(fixed_model(p), d, ~ survival::Surv(time, status), horizon = tau,
                 thresholds = z)
  expect_lt(max(abs(ci$NB_model - truth)), 0.01)
  expect_gt(max(abs(km$NB_model - truth)), 0.02)
  # a censoring model without the covariate driving censoring is biased too
  cens0 <- survival::coxph(survival::Surv(time, status == 0) ~ x1, data = d)
  ci0 <- no_fixed_warning(cipcw_nb(fixed_model(p), d, ~ survival::Surv(time, status),
                                   cipcw_control(model = cens0), horizon = tau, thresholds = z))
  expect_gt(max(abs(ci0$NB_model - truth)), 0.02)
})

test_that("cIPCW is close to the truth with competing risks and informative censoring", {
  tau <- 5; z <- c(0.1, 0.2, 0.3)
  d <- sim_surv(20000, informative = TRUE, seed = 12, competing = TRUE)
  p <- plogis(d$lp1 - 0.7)
  truth <- true_nb(d, p, z, tau)
  cens <- survival::coxph(survival::Surv(time, status == 0) ~ x2, data = d)
  ci <- no_fixed_warning(cipcw_nb(fixed_model(p), d, ~ survival::Surv(time, status_f),
                                  cipcw_control(model = cens), horizon = tau, thresholds = z))
  expect_lt(max(abs(ci$NB_model - truth)), 0.01)
  # cause 2
  d2 <- d; d2$cause <- 3L - d2$cause
  truth2 <- true_nb(d2, p, z, tau)
  ci2 <- no_fixed_warning(cipcw_nb(fixed_model(p), d, ~ survival::Surv(time, status_f),
                                   cipcw_control(model = cens), horizon = tau,
                                   thresholds = z, cause = "cause2"))
  expect_lt(max(abs(ci2$NB_model - truth2)), 0.01)
})

test_that("cIPCW with case weights equals duplicated rows", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val[1:500, ]
  set.seed(9); w <- rpois(500, 1)
  we <- runif(500, 1, 3); wh <- runif(500, 1, 2)
  z <- seq(0, 0.6, 0.05)
  a <- cipcw_nb(cx, v, surv_outcome, cipcw_control(weight_event = we, weight_horizon = wh),
                thresholds = z, weights = w)
  i <- rep(seq_len(500), w)
  b <- cipcw_nb(cx, v[i, ], surv_outcome,
                cipcw_control(weight_event = we[i], weight_horizon = wh[i]), thresholds = z)
  expect_equal(a$NB_model, b$NB_model, tolerance = 1e-12)
  expect_equal(attr(a, "prevalence"), attr(b, "prevalence"), tolerance = 1e-12)
})

test_that("a censoring model with case weights equals duplicated rows and dropped rows", {
  # The censoring model is fitted once (unweighted) and kept fixed; case
  # weights then act like duplicating rows (w = 2) or dropping them (w = 0).
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val[1:600, ]
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade,
                          data = v)
  ctrl <- cipcw_control(model = cens)
  z <- seq(0, 0.6, 0.05)
  set.seed(10); w <- rpois(600, 1)
  a <- cipcw_nb(cx, v, surv_outcome, ctrl, thresholds = z, weights = w)
  # duplicated rows (the check of the censoring model cannot match the rows)
  i <- rep(seq_len(600), w)
  b <- suppress_matching(cipcw_nb(cx, v[i, ], surv_outcome, ctrl, thresholds = z),
                         "Could not verify|fitted to .* patients")
  expect_equal(a$NB_model, b$NB_model, tolerance = 1e-12)
  expect_equal(a$NB_all, b$NB_all, tolerance = 1e-12)
  expect_equal(attr(a, "prevalence"), attr(b, "prevalence"), tolerance = 1e-12)
  # the same as giving the model's weights by hand, with the same case weights
  we <- ifelse(v$death == 1 & v$dtimey <= 5,
               1 / surv_at(cens, v, v$dtimey, "left"), 1)
  wh <- 1 / surv_at(cens, v, rep(5, 600), "left")
  c2 <- cipcw_nb(cx, v, surv_outcome, cipcw_control(weight_event = we, weight_horizon = wh),
                 thresholds = z, weights = w)
  expect_equal(a$NB_model, c2$NB_model, tolerance = 1e-12)
  # Dirichlet (non-integer) case weights run and stay finite
  set.seed(11); wd <- rexp(600); wd <- wd / sum(wd) * 600
  d <- cipcw_nb(cx, v, surv_outcome, ctrl, thresholds = z, weights = wd)
  expect_true(all(is.finite(d$NB_model)))
  # competing risks: case weights 0/1 equal dropping the rows with weight 0
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  cens_cr <- survival::coxph(survival::Surv(cr_time, cr_code == 0) ~ age + nodes, data = v)
  keep <- rbinom(600, 1, 0.7)
  e <- cipcw_nb(csc, v, cr_outcome, cipcw_control(model = cens_cr), thresholds = z,
                weights = keep, cause = "death")
  f <- suppress_matching(cipcw_nb(csc, v[keep == 1, ], cr_outcome,
                                  cipcw_control(model = cens_cr), thresholds = z,
                                  cause = "death"), "Could not verify|fitted to .* patients")
  expect_equal(e$NB_model, f$NB_model, tolerance = 1e-12)
  expect_equal(attr(e, "prevalence"), attr(f, "prevalence"), tolerance = 1e-12)
})

test_that("print shows the cIPCW estimator", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age, data = s$val)
  r <- cipcw_nb(cx, s$val, surv_outcome, cipcw_control(model = cens), thresholds = 0.2)
  expect_output(print(r), "survival, horizon = 5, cipcw")
})

test_that("censoring predict_function: same as the model route, per row or as a curve", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade,
                          data = v, model = TRUE)
  # a curve per patient, read by voipred at each patient's own time
  curve <- function(model, newdata, times) {
    sf <- survival::survfit(model, newdata = newdata)
    list(time = sf$time, surv = t(sf$surv))
  }
  for (ties in c("left", "right")) {
    a <- cipcw_nb(cx, v, surv_outcome, cipcw_control(model = cens, ties = ties))
    b <- cipcw_nb(cx, v, surv_outcome, cipcw_control(model = cens, ties = ties,
                                                     predict_function = curve))
    expect_equal(b$NB_model, a$NB_model, tolerance = 1e-10)
    expect_equal(attr(b, "prevalence"), attr(a, "prevalence"), tolerance = 1e-10)
  }
  # one probability per row, at that row's time
  sr <- survival::survreg(survival::Surv(dtimey, death == 0) ~ age + nodes + grade, data = v)
  per_row <- function(model, newdata, times)
    1 - survival::psurvreg(times, predict(model, newdata, type = "lp"), model$scale, model$dist)
  a <- cipcw_nb(cx, v, surv_outcome, cipcw_control(model = sr))
  b <- cipcw_nb(cx, v, surv_outcome, cipcw_control(model = sr, predict_function = per_row))
  expect_equal(b$NB_model, a$NB_model, tolerance = 1e-8)
  expect_output(print(cipcw_control(model = sr, predict_function = per_row)),
                "censoring model \\(survreg, with predict_function\\)")
})

test_that("censoring predict_function makes any censoring model usable", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade,
                          data = v, model = TRUE)
  # a class voipred knows nothing about
  mine <- structure(list(inner = cens), class = "my_censoring_model")
  expect_error(suppress_matching(cipcw_nb(cx, v, surv_outcome, cipcw_control(model = mine)),
                                 "Could not verify"),
               "Could not get survival probabilities from this my_censoring_model model.*`predict_function` in cipcw_control\\(\\)")
  fun <- function(model, newdata, times) {
    sf <- survival::survfit(model$inner, newdata = newdata)
    list(time = sf$time, surv = t(sf$surv))
  }
  expect_warning(b <- cipcw_nb(cx, v, surv_outcome,
                               cipcw_control(model = mine, predict_function = fun)),
                 "Could not verify that the censoring model")
  a <- cipcw_nb(cx, v, surv_outcome, cipcw_control(model = cens))
  expect_equal(b$NB_model, a$NB_model, tolerance = 1e-10)
})

test_that("censoring predict_function: input and output checks", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = v)
  expect_error(cipcw_control(model = cens, predict_function = "f"), "must be a function")
  expect_error(cipcw_control(predict_function = function(...) 1), "needs the censoring `model`")
  run <- function(f) cipcw_nb(cx, v, surv_outcome, cipcw_control(model = cens, predict_function = f))
  expect_error(run(function(model, newdata, times) rep(0.9, 3)), "returned 3 values for")
  expect_error(run(function(model, newdata, times) rep(NA_real_, nrow(newdata))),
               "invalid censoring probabilities")
  expect_error(run(function(model, newdata, times) rep(1.2, nrow(newdata))),
               "invalid censoring probabilities")
  expect_error(run(function(model, newdata, times) cbind(rep(0.9, nrow(newdata)), 0.1)),
               "one probability per row")
  expect_error(run(function(model, newdata, times) list(time = 1:3, surv = matrix(0.9, 2, 3))),
               "curve of the wrong size")
  expect_error(run(function(model, newdata, times) ifelse(times > 4, 0, 0.9)),
               "exactly 0, which gives an infinite weight")
})

test_that("more censoring model classes: psm, flexsurvspline, stratified survreg", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  v <- s$val
  nd <- v[1:20, ]
  cf <- survival::Surv(dtimey, death == 0) ~ age + nodes + grade
  ps <- rms::psm(cf, data = v, dist = "weibull")
  expect_equal(surv_at(ps, nd, rep(5, 20)), as.vector(rms::survest(ps, newdata = nd, times = 5)$surv),
               tolerance = 1e-6)
  fsp <- flexsurv::flexsurvspline(survival::Surv(dtimey, death == 0) ~ age + nodes + grade,
                                  data = v, k = 1)
  expect_equal(surv_at(fsp, nd, rep(5, 20)),
               predict(fsp, nd, type = "survival", times = 5)$.pred_survival, tolerance = 1e-6)
  st <- survival::survreg(survival::Surv(dtimey, death == 0) ~ age + nodes + strata(grade), data = v)
  lp <- predict(st, nd, type = "lp"); sc <- st$scale[paste0("grade=", nd$grade)]
  expect_equal(surv_at(st, nd, rep(5, 20)), unname(exp(-exp((log(5) - lp) / sc))),
               tolerance = 1e-10)
  models <- list(ps, fsp, st)
  for (m in models) {
    res <- suppress_matching(cipcw_nb(cx, v, surv_outcome, cipcw_control(model = m)),
                             "Could not verify|below 0.01|No patient treated|fitted to .* patients")
    expect_true(all(is.finite(res$NB_all)))
    expect_true(attr(res, "prevalence") > 0.2 && attr(res, "prevalence") < 0.32)
  }
})
