# Value of information: draws_t(), evpi() (and relative EVPI), evsi(),
# p_useful() and p_best(). The Bayesian-bootstrap engine (used for survival
# and competing-risks outcomes and case weights) is checked against the code of
# Sadatsafavi et al. (predtools::evpi_val for EVPI; the bootstrap EVSI
# algorithm of the EVSI paper, resplab/papercode/evsiexval) on the binary case
# study, with the same draws; the conjugate posterior of binary outcomes is
# checked in test-bayes.R. And against direct calculations.

test_that("draws with all weights 1 reproduce nb() for every outcome and estimator", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  z <- c(0, 0.02, 0.05, 0.1, 0.3)
  r <- nb(g$val, f, "day30", thresholds = z)
  d <- fixed_draws(r, W = matrix(1, nrow(g$val), 2))
  expect_equal(dim(d$nb), c(2, length(z), 3))
  expect_equal(d$strategies, c("f", "Treat All", "Treat None"))
  expect_equal(d$nb[1, , 1], r$estimates$NB_model, tolerance = 1e-12)
  expect_equal(d$nb[2, , 2], r$estimates$NB_all, tolerance = 1e-12)
  expect_true(all(d$nb[, , 3] == 0))

  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(surv_formula, data = s$dev)
  zs <- c(0.1, 0.2, 0.3, 0.5)
  ones <- matrix(1, nrow(v), 2)
  km <- quiet_nb(v, cx, surv_outcome, horizon = 5, thresholds = zs)
  dk <- fixed_draws(km, W = ones)
  expect_equal(dk$nb[1, , 1], km$estimates$NB_model, tolerance = 1e-12)
  expect_equal(dk$nb[1, , 2], km$estimates$NB_all, tolerance = 1e-12)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade, data = v)
  ci <- nb(v, cx, surv_outcome, horizon = 5, thresholds = zs, censoring = "cipcw",
           cipcw = cipcw_control(model = cens))
  dc <- fixed_draws(ci, W = ones)
  expect_equal(dc$nb[2, , 1], ci$estimates$NB_model, tolerance = 1e-12)
  expect_equal(dc$nb[2, , 2], ci$estimates$NB_all, tolerance = 1e-12)
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  cr <- quiet_nb(v, csc, cr_outcome, horizon = 5, thresholds = zs, cause = "death")
  dr <- suppressWarnings(fixed_draws(cr, W = ones))   # NA where nb() is NA
  expect_equal(dr$nb[1, , 1], cr$estimates$NB_model, tolerance = 1e-12)
  expect_equal(dr$nb[1, , 2], cr$estimates$NB_all, tolerance = 1e-12)
  # integer draw weights equal duplicated rows
  set.seed(3)
  w <- rpois(nrow(v), 1)
  dup <- quiet_nb(v[rep(seq_len(nrow(v)), w), ], cx, surv_outcome, horizon = 5, thresholds = zs)
  expect_equal(fixed_draws(km, W = cbind(w, w))$nb[1, , 1], dup$estimates$NB_model,
               tolerance = 1e-12)
})

test_that("EVPI and P(useful) equal predtools::evpi_val() (Sadatsafavi et al. 2023)", {
  cs <- voi_case()
  ref <- cs$ref
  expect_equal(cs$pi, ref$pi, tolerance = 1e-12)
  res <- nb(cs$val, cs$fit, "day30", thresholds = ref$zs, label = "Model")
  for (m in c("bootstrap", "bayesian_bootstrap")) {        # predtools' method names
    set.seed(42)
    W <- vapply(1:1000, function(i) predtools_bootstrap(500, m == "bayesian_bootstrap"),
                numeric(500))
    d <- fixed_draws(res, W = W)
    e <- evpi_from_draws(d)
    expect_s3_class(e, "voipred_evpi")
    expect_equal(e$estimates$EVPI, ref[[m]]$EVPIv, tolerance = 1e-12, label = m)
    expect_equal(e$estimates$ENB_current, ref[[m]]$ENB_current, tolerance = 1e-12)
    expect_equal(e$estimates$ENB_perfect, ref[[m]]$ENB_perfect, tolerance = 1e-12)
    expect_equal(unname(e$ENB[, "Model"]), ref[[m]]$ENB_model, tolerance = 1e-12)
    expect_equal(unname(e$ENB[, "Treat All"]), ref[[m]]$ENB_all, tolerance = 1e-12)
    # predtools counts a draw where the model ties with a default strategy as
    # useful; voipred does not. Without ties they agree.
    pu <- p_useful_from_draws(d)
    tie <- colSums(d$nb[, , 1] == pmax(d$nb[, , 2], 0)) > 0
    expect_true(any(!tie))
    expect_equal(pu$estimates$p_useful[!tie], ref[[m]]$p_useful[!tie], label = m)
  }
})

test_that("EVSI equals the bootstrap algorithm of Sadatsafavi et al. (2025)", {
  cs <- voi_case()
  fs <- c(500, 1000, 2000, 4000)
  for (z in c(0.017, 0.023)) {
    res <- nb(cs$val, cs$fit, "day30", thresholds = z)
    set.seed(7)
    r <- reference_evsi(cs$val$day30, cs$pi, z, fs, 300, "BBS")
    set.seed(7)
    o <- evsi_bb(res, future_n = fs, n_sim = 300)
    expect_equal(as.vector(o$EVSI), r$EVSI, tolerance = 1e-10, label = paste(z))
    expect_equal(o$EVPI, r$EVPI, tolerance = 1e-10)
  }
  # several thresholds at once give the same as one at a time (same draws)
  res2 <- nb(cs$val, cs$fit, "day30", thresholds = c(0.017, 0.023))
  set.seed(7)
  both <- evsi_bb(res2, future_n = fs, n_sim = 300)
  set.seed(7)
  one <- evsi_bb(nb(cs$val, cs$fit, "day30", thresholds = 0.023), future_n = fs, n_sim = 300)
  expect_equal(both$EVSI[2, ], one$EVSI[1, ], tolerance = 1e-12)
})

test_that("EVSI is between 0 and EVPI, grows with the future sample, and reaches EVPI", {
  cs <- voi_case()
  res <- nb(cs$val, cs$fit, "day30", thresholds = c(0.017, 0.023, 0.3))
  set.seed(11)
  o <- evsi(res, future_n = c(100, 1000, 10000, 1e6), n_sim = 1000)
  expect_true(all(o$estimates$EVSI <= o$estimates$EVPI + 1e-15))
  expect_true(all(o$estimates$EVSI >= -1e-5))
  for (z in c(0.017, 0.023)) {
    e <- o$estimates$EVSI[o$estimates$threshold == z]
    expect_true(all(diff(e) > 0))
    expect_gt(e[4] / o$estimates$EVPI[o$estimates$threshold == z][1], 0.95)
  }
  # the EVPI of evsi() equals evpi() with the same draws (conjugate posterior)
  set.seed(12)
  o2 <- evsi(res, future_n = 100, n_sim = 1000)
  set.seed(12)
  expect_equal(evpi(res, n_sim = 1000)$estimates$EVPI, o2$estimates$EVPI, tolerance = 1e-10)
  # and for the Bayesian bootstrap (same random numbers)
  set.seed(11)
  ob <- evsi_bb(res, future_n = c(100, 1000, 10000, 1e6), n_sim = 1000)
  set.seed(11)
  W <- vapply(1:1000, function(i) {
    w <- stats::rgamma(500, 1, 1)
    for (n in c(100, 1000, 10000, 1e6)) stats::rmultinom(1, n, w / sum(w))
    w
  }, numeric(500))
  expect_equal(evpi_from_draws(fixed_draws(res, W = W))$estimates$EVPI, ob$EVPI, tolerance = 1e-10)
  expect_true(all(ob$EVSI <= ob$EVPI + 1e-15))
})

test_that("P(useful), P(best) and relative EVPI follow their definitions", {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  z <- c(0, 0.02, 0.05, 0.1, 0.2)
  r1 <- nb(g$val, f1, "day30", thresholds = z, label = "Full")
  r2 <- nb(g$val, f2, "day30", thresholds = z, label = "Age")
  set.seed(5)
  d <- draws_t(r1, r2, n_sim = 300)
  expect_equal(d$strategies, c("Full", "Age", "Treat All", "Treat None"))
  nbs <- d$nb
  # P(useful) = P(NB_model > max(NB_all, 0))
  pu <- p_useful_from_draws(d)
  expect_equal(pu$estimates$model, rep(c("Full", "Age"), each = length(z)))
  for (m in 1:2)
    expect_equal(pu$estimates$p_useful[pu$estimates$model == d$models[m]],
                 colMeans(nbs[, , m] > pmax(nbs[, , 3], 0)))
  # at threshold 0 every model treats everyone: it ties with treat all, so it is not useful
  expect_equal(pu$estimates$p_useful[pu$estimates$threshold == 0], c(0, 0))
  # P(best): shares of the draws, ties split, summing to 1
  pb <- p_best_from_draws(d)
  expect_equal(as.vector(tapply(pb$estimates$p_best, pb$estimates$threshold, sum)), rep(1, length(z)))
  best <- apply(nbs, c(1, 2), which.max)
  expect_equal(pb$estimates$p_best[pb$estimates$strategy == "Full" & pb$estimates$threshold == 0.1],
               mean(best[, 4] == 1))
  expect_equal(pb$estimates$p_best[pb$estimates$threshold == 0], c(1, 1, 1, 0) / 3)   # three-way tie at z = 0
  # EVPI and relative EVPI by hand
  e <- evpi_from_draws(d)
  enb <- apply(nbs, c(2, 3), mean)
  perfect <- colMeans(apply(nbs, c(1, 2), max))
  current <- apply(enb, 1, max)
  default <- pmax(enb[, 3], 0)
  expect_equal(e$estimates$EVPI, perfect - current)
  expect_true(all(e$estimates$EVPI >= 0))
  ok <- current > default
  expect_equal(e$estimates$rEVPI[ok], ((perfect - default) / (current - default))[ok])
  expect_true(is.na(e$estimates$rEVPI[1]))                 # z = 0: no gain and no EVPI, 0/0
  # results from nb() directly give the same as from the draws
  set.seed(5)
  expect_equal(evpi(r1, r2, n_sim = 300)$estimates$EVPI, e$estimates$EVPI)
  set.seed(5)
  expect_equal(p_best(r1, r2, n_sim = 300)$estimates$p_best, pb$estimates$p_best)
})

test_that("relative EVPI is Inf when no model is better than the defaults but EVPI > 0", {
  set.seed(8)
  dat <- data.frame(y = rbinom(300, 1, 0.3))
  res <- nb(dat, fixed_model(runif(300)), "y", thresholds = c(0.25, 0.3, 0.35))
  e <- evpi(res, n_sim = 500)
  inf <- e$estimates$ENB_current == pmax(e$ENB[, "Treat All"], 0) & e$estimates$EVPI > 0
  expect_true(any(inf))
  expect_true(all(is.infinite(e$estimates$rEVPI[inf])))
})

test_that("a survival outcome without censoring before the horizon gives the binary results", {
  set.seed(9)
  n <- 400
  dat <- data.frame(x = rnorm(n))
  dat$y <- rbinom(n, 1, plogis(-1 + dat$x))
  dat$time <- ifelse(dat$y == 1, runif(n, 0, 4), 10)
  p <- plogis(-1 + 0.8 * dat$x)
  z <- c(0.1, 0.2, 0.3, 0.5)
  rb <- nb(dat, fixed_model(p), "y", thresholds = z)
  rs <- nb(dat, fixed_model(p), ~ survival::Surv(time, y), horizon = 5, thresholds = z)
  set.seed(1)
  W <- matrix(rexp(n * 200), n)
  db <- fixed_draws(rb, W = W)
  ds <- fixed_draws(rs, W = W)
  expect_equal(ds$nb, db$nb, tolerance = 1e-12)
  expect_equal(evpi_from_draws(ds)$estimates$EVPI, evpi_from_draws(db)$estimates$EVPI, tolerance = 1e-12)
  set.seed(2)
  eb <- evsi_bb(rb, future_n = c(200, 2000), n_sim = 100)   # the bootstrap EVSI of binary
  set.seed(2)
  es <- evsi(rs, future_n = c(200, 2000), n_sim = 100)
  expect_equal(es$draws, "Bayesian bootstrap")
  expect_equal(es$estimates$EVSI, as.vector(eb$EVSI), tolerance = 1e-12)
})

test_that("survival, competing risks and cIPCW: EVPI, EVSI, P(useful), P(best) run and are coherent", {
  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(surv_formula, data = s$dev)
  small <- survival::coxph(survival::Surv(dtimey, death) ~ age + size, data = s$dev)
  zs <- c(0.15, 0.25, 0.35)
  a <- quiet_nb(v, cx, surv_outcome, horizon = 5, thresholds = zs, label = "Cox")
  b <- quiet_nb(v, small, surv_outcome, horizon = 5, thresholds = zs, label = "Small")
  set.seed(4)
  d <- draws_t(a, b, n_sim = 100)
  e <- evpi_from_draws(d)
  expect_true(all(e$estimates$EVPI >= 0))
  expect_equal(nrow(p_useful_from_draws(d)$estimates), 6)
  pb <- p_best_from_draws(d)
  expect_equal(as.vector(tapply(pb$estimates$p_best, pb$estimates$threshold, sum)), rep(1, 3))
  set.seed(4)
  o <- evsi(a, b, future_n = c(500, 5000), n_sim = 60)
  expect_true(all(o$estimates$EVSI <= o$estimates$EVPI + 1e-15))
  # cIPCW (censoring weights fixed) and competing risks
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes + grade, data = v)
  ci <- nb(v, cx, surv_outcome, horizon = 5, thresholds = zs, censoring = "cipcw",
           cipcw = cipcw_control(model = cens))
  set.seed(4)
  expect_true(all(evpi(ci, n_sim = 100)$estimates$EVPI >= 0))
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  cr <- quiet_nb(v, csc, cr_outcome, horizon = 5, thresholds = zs, cause = "death")
  set.seed(4)
  oc <- evsi(cr, future_n = 1000, n_sim = 60)
  expect_true(all(oc$estimates$EVSI <= oc$estimates$EVPI + 1e-15))
})

test_that("thresholds where the net benefit cannot be estimated give NA, with a warning", {
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  z <- seq(0.5, 0.95, by = 0.05)
  r <- quiet_nb(s$val, cx, surv_outcome, horizon = 5, thresholds = z)
  skip_if(!anyNA(r$estimates$NB_model))
  set.seed(1)
  expect_warning(d <- draws_t(r, n_sim = 20), "could not be estimated in some draws")
  e <- evpi_from_draws(d)
  expect_true(all(is.na(e$estimates$EVPI[is.na(r$estimates$NB_model)])))
})

test_that("case weights are part of the draws", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  set.seed(2)
  w0 <- runif(nrow(g$val), 0.5, 2)
  r <- nb(g$val, f, "day30", thresholds = c(0.05, 0.1), weights = w0)
  set.seed(3)
  d <- draws_t(r, n_sim = 50)
  set.seed(3)
  W <- vapply(1:50, function(i) { u <- rgamma(length(w0), 1, 1) * w0; u * length(w0) / sum(u) },
              numeric(length(w0)))
  expect_equal(d$nb, fixed_draws(r, W = W)$nb, tolerance = 1e-12)
})

test_that("inputs are checked", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  r <- nb(g$val, f, "day30", thresholds = c(0.05, 0.1))
  b <- nb(g$val, f, "day30", thresholds = c(0.05, 0.1), label = "b")
  expect_error(evsi(r), "`future_n` must be")
  expect_error(evsi(r, future_n = c(100, 1.5)), "`future_n` must be")
  # no `method`: the draws come from the posterior of nb()'s results
  expect_error(evpi(r, method = "bootstrap"), "Unknown argument `method`")
  expect_error(p_useful(r, method = "ordinary"), "Unknown argument `method`")
  expect_error(evsi(r, future_n = 10, method = "asymptotic"), "Unknown argument `method`")
  # the draws are internal: no nb_draws() and no draws as input
  expect_false("nb_draws" %in% getNamespaceExports("voipred.dev"))
  set.seed(1)
  expect_error(evpi(draws_t(r, n_sim = 10)), "needs results of nb")
  expect_error(evpi(r, nsim = 10), "Unknown argument `nsim`.*evpi\\(\\)")
  expect_error(p_best(r, n_sim = 1), "at least 2")
  expect_error(evpi(r, draw_weights = matrix(1, 10, 3)), "Unknown argument `draw_weights`")
  expect_error(evpi(r, nb(g$dev, f, "day30", thresholds = c(0.05, 0.1), label = "b")),
               "cannot be used together: different data")
  expect_error(p_useful(data.frame(x = 1)), "needs results of nb")
  expect_error(p_useful(r, r), "p_useful\\(\\) takes one result of nb\\(\\)")
  expect_error(p_useful(r, labels = "a"), "Unknown argument `labels`")
  expect_error(evpi(as.data.frame(r)), "needs results of nb")
})

