# nb_informative(): bayesDCA's threshold-varying ("informative") prior.

informative_case <- function() {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  val <- g$val
  val$p1 <- unname(predict(f1, val, type = "response"))
  val$p2 <- unname(predict(f2, val, type = "response"))
  list(val = val)
}

test_that("the prior shapes are bayesDCA's, number for number", {
  ref <- readRDS(test_path("fixtures", "bayesdca_informative_prior.rds"))
  priors <- list(
    nb_informative(),
    nb_informative(),
    nb_informative(),
    nb_informative(sensitivity = c(0.95, 0.5), strength = 10, ignorance = NULL),
    nb_informative(ignorance = c(0.1, 0.9), ignorance_prior = c(1.2, 1.8),
                   prevalence = c(7, 93)))
  for (k in seq_along(ref)) {
    ours <- as.matrix(prior_shapes(priors[[k]], ref[[k]]$thresholds)[, -1])
    expect_equal(unname(ours), unname(ref[[k]]$shapes), tolerance = 1e-12,
                 label = ref[[k]]$label)
  }
  # the setting of the EVPI simulation of Cruz and Korthauer (Table S1)
  s <- prior_shapes(nb_informative(sensitivity = c(0.95, 0.5), strength = 10,
                                   ignorance = NULL), c(0.01, 0.02, 0.05, 0.1))
  expect_equal(s$se_a, c(9.5, 9, 7.5, 5))
  expect_equal(s$se_b, c(0.5, 1, 2.5, 5))
  expect_equal(s$sp_a, s$se_b)
  expect_equal(s$sp_b, s$se_a)
})

test_that("the defaults are bayesDCA's: ends 0.99 / 0.01, strength 5, vague region uniform", {
  s <- prior_shapes(nb_informative(), seq(0, 0.5, by = 0.01))
  size <- s$se_a + s$se_b
  mu <- s$se_a / size
  expect_equal(mu[1], 0.99)
  expect_equal(size[1], 5)
  expect_equal(mu[51], 0.01)
  expect_equal(size[51], 5)
  mid <- s$threshold > 0.125 + 1e-9 & s$threshold <= 0.375
  expect_true(all(abs(s$se_a[mid] - 1) < 1e-12 & abs(s$se_b[mid] - 1) < 1e-12))
  expect_equal(s$prev_a, rep(1, 51))
  expect_equal(s$prev_b, rep(1, 51))
  # the other priors are the same at every threshold
  u <- prior_shapes(nb_uniform(), c(0.1, 0.2))
  expect_equal(unlist(u[1, -1]), unlist(u[2, -1]))
  expect_equal(unname(unlist(u[1, -1])), rep(1, 6))
})

test_that("nb_val(): the posterior means use the prior of each threshold", {
  g <- gusto_split()
  v <- g$val[1:300, ]
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  v$p <- unname(predict(f, v, type = "response"))
  z <- seq(0, 0.4, by = 0.02)
  res <- nb_val(v, "p", "day30", thresholds = z, prior = nb_informative())
  expect_equal(res$prior$type, "informative")
  sh <- prior_shapes(nb_informative(), z)
  y <- v$day30
  tr <- vapply(z, function(t) v$p > 0 & v$p >= t, logical(nrow(v)))
  tp <- colSums(y * tr)
  tn <- colSums((1 - y) * !tr)
  ep <- (1 + sum(y)) / (2 + nrow(v))
  se <- (sh$se_a + tp) / (sh$se_a + sh$se_b + sum(y))
  sp <- (sh$sp_a + tn) / (sh$sp_a + sh$sp_b + sum(1 - y))
  expect_equal(res$estimates$NB_model, ep * se - (1 - ep) * (1 - sp) * z / (1 - z))
  expect_equal(res$estimates$NB_all, ep - (1 - ep) * z / (1 - z))
  # posterior draws (one model, bayesDCA's draws): their means match
  set.seed(1)
  d <- nb_draws(res, 20000)
  expect_lt(max(abs(colMeans(d$nb[, , 1]) - res$estimates$NB_model)), 0.002)
  expect_output(print(res), "prior informative \\(bayesDCA, threshold-varying\\)")
  expect_output(print(nb_informative()),
                "sensitivity mean 0.99 at the lowest threshold to 0.01 at the highest")
})

test_that("several models: each model's sensitivity and specificity has its prior at each threshold", {
  m <- informative_case()
  v <- m$val
  z <- seq(0.01, 0.3, by = 0.01)
  res <- nb_val(v, list(Full = v$p1, Age = v$p2), "day30", thresholds = z,
                prior = nb_informative())
  s <- bayes_setup(res)
  sh <- prior_shapes(nb_informative(), z)
  for (mm in 1:2) {
    expect_equal(as.vector(s$alpha1 %*% s$pos[, mm]), sh$se_a)          # treated, events
    expect_equal(as.vector(s$alpha1 %*% !s$pos[, mm]), sh$se_b)
    expect_equal(as.vector(s$alpha0 %*% !s$pos[, mm]), sh$sp_a)         # untreated, non-events
    expect_equal(as.vector(s$alpha0 %*% s$pos[, mm]), sh$sp_b)
  }
  # each model alone gives the same point estimates
  a <- nb_val(v, list(Full = v$p1), "day30", thresholds = z, prior = nb_informative())
  expect_equal(res$estimates$NB_model[res$estimates$model == "Full"], a$estimates$NB_model)
  # the value-of-information functions use it
  set.seed(2)
  e <- evpi(res, n_sim = 300)
  expect_equal(e$prior$type, "informative")
  expect_s3_class(e$p_best, "voipred_p_best")
  set.seed(3)
  o <- evsi(res, future_n = c(200, 1000), n_sim = 200)
  expect_true(all(o$estimates$EVSI <= o$estimates$EVPI + 1e-12))
  # results with the same prior combine; with another prior they do not
  b <- nb_val(v, list(Age = v$p2), "day30", thresholds = z, prior = nb_informative())
  expect_equal(c(a, b)$estimates, res$estimates)
  u <- nb_val(v, list(Age = v$p2), "day30", thresholds = z)
  expect_error(c(a, u), "different prior")
})

test_that("the fixed priors give the same results as before (one prior row per threshold)", {
  m <- informative_case()
  v <- m$val
  z <- c(0.02, 0.05, 0.1)
  for (pr in list(nb_uniform(), nb_jeffreys(), nb_beta(sensitivity = c(3, 2)))) {
    res <- nb_val(v, list(Full = v$p1, Age = v$p2), "day30", thresholds = z, prior = pr)
    s <- bayes_setup(res)
    expect_equal(nrow(s$alpha1), length(z))
    expect_equal(s$alpha1[1, ], s$alpha1[3, ])
    expect_equal(as.vector(s$alpha1 %*% s$pos[, 1]), rep(pr$sensitivity[1], 3))
  }
})

test_that("nb_informative() checks its inputs and where it can be used", {
  expect_error(nb_informative(sensitivity = 0.9), "two prior means in \\[0, 1\\]")
  expect_error(nb_informative(sensitivity = c(1.2, 0)), "two prior means")
  expect_error(nb_informative(strength = 0), "single positive number")
  expect_error(nb_informative(ignorance = c(0.8, 0.2)), "two increasing fractions")
  expect_error(nb_informative(ignorance_prior = c(0, 0)), "positive total")
  expect_error(nb_informative(prevalence = c(-1, 1)), "non-negative numbers")
  g <- gusto_split()
  v <- g$val[1:300, ]
  v$p <- 0.1 + 0.5 * (v$age - min(v$age)) / diff(range(v$age))
  expect_error(nb_val(v, "p", "day30", thresholds = 0.1, prior = nb_informative()),
               "at least two different thresholds")
  expect_error(nb_val(v, "p", "day30", thresholds = seq(0.2, 0.5, by = 0.05),
                      prior = nb_informative()),
               "vague region of nb_informative\\(\\) starts at 0.125")
  expect_no_error(nb_val(v, "p", "day30", thresholds = seq(0.2, 0.5, by = 0.05),
                         prior = nb_informative(ignorance = NULL)))
  expect_error(nb_val(v, "p", "day30", thresholds = c(0.1, 0.2), prior = nb_informative(),
                      weights = rep(1, nrow(v))), "only be used with the improper prior")
  s <- rott_split()
  s$val$p <- 0.3
  s$val$p[1:10] <- 0.5
  expect_error(nb_val(s$val, "p", surv_outcome, horizon = 5, prior = nb_informative()),
               "other priors are for binary outcomes")
  expect_error(nb_val(v, "p", "day30", prior = list(type = "informative")),
               "nb_informative\\(\\)")
  expect_error(prior_shapes(list(), 0.1), "must be a prior")
})

test_that("the prior shapes are kept in the result; prior_shapes() is internal", {
  g <- gusto_split()
  v <- g$val[1:300, ]
  v$p <- 0.05 + 0.4 * (v$age - min(v$age)) / diff(range(v$age))
  z <- seq(0, 0.4, by = 0.05)
  res <- nb_val(v, "p", "day30", thresholds = z, prior = nb_informative())
  expect_equal(res$prior_shapes, prior_shapes(nb_informative(), z))
  expect_equal(res$prior_shapes$threshold, z)
  # the default: Beta(0, 0) everywhere
  expect_true(all(nb_val(v, "p", "day30", thresholds = z)$prior_shapes[, -1] == 0))
  # survival outcomes: no Beta priors
  s <- rott_split()
  s$val$p <- 0.3
  s$val$p[1:10] <- 0.5
  expect_null(nb_val(s$val, "p", surv_outcome, horizon = 5, thresholds = 0.2)$prior_shapes)
  expect_false("prior_shapes" %in% getNamespaceExports("voipred.dev"))
})

test_that("the EVSI output states the information before and after the future study", {
  g <- gusto_split()
  v <- g$val[1:400, ]
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  v$p <- unname(predict(f, v, type = "response"))
  z <- seq(0.01, 0.2, by = 0.01)
  set.seed(1)
  o <- evsi(nb_val(v, "p", "day30", thresholds = z), future_n = c(500, 2000), n_sim = 50)
  out <- capture.output(print(o))
  expect_match(out[1], "^EVSI: p; binary, n = 400, 50 posterior draws$")
  expect_equal(out[2], paste0("Before the future study (as for EVPI): prior improper Beta(0, 0), ",
                              "updated with the validation data (n = 400)"))
  expect_match(out[3], paste0("^Future study: 500 or 2,000 more patients, drawn from that ",
                              "distribution; then the prior is updated"))
  expect_output(print(summary(o)), "Before the future study \\(as for EVPI\\)")
  set.seed(1)
  oi <- evsi(nb_val(v, "p", "day30", thresholds = z, prior = nb_informative()),
             future_n = 500, n_sim = 50)
  expect_output(print(oi), "prior informative \\(bayesDCA, threshold-varying\\), updated with")
  # survival: Bayesian bootstrap, the data only
  s <- rott_split()
  s$val$p <- 0.2
  s$val$p[1:200] <- 0.4
  set.seed(2)
  os <- evsi(nb_val(s$val, "p", surv_outcome, horizon = 5, thresholds = c(0.1, 0.3)),
             future_n = 300, n_sim = 5)
  out_s <- capture.output(print(os))
  expect_match(out_s[2], "Bayesian bootstrap \\(prior improper Beta\\(0, 0\\): the data only\\)")
  expect_match(out_s[3], "^Future study: 300 more patients, resampled from each draw")
  # EVPI keeps its one-line header
  set.seed(3)
  e <- evpi(nb_val(v, "p", "day30", thresholds = z), n_sim = 50)
  expect_match(capture.output(print(e))[1], "^EVPI: p; binary, n = 400, prior improper Beta\\(0, 0\\), 50 posterior draws$")
})
