# The value of information of the Bayesian net benefit (nb(bayesian = TRUE));
# the point estimates are in test-bayes-nb.R.
# Checked against bayesDCA (giulianonetto/bayesDCA, .dca_binary) and the
# conjugate EVSI() of Sadatsafavi et al. 2025 (resplab/evsiexval): see
# fixtures/make_bayes_reference.R.

bayes_case <- function() {
  cs <- voi_case()
  cs$bref <- readRDS(test_path("fixtures", "bayes_reference.rds"))
  cs
}

test_that("one model: posterior draws equal bayesDCA's", {
  cs <- bayes_case()
  b <- cs$bref
  r <- nb(cs$fit, cs$val, "day30", thresholds = b$zs, bayesian = TRUE)
  set.seed(b$seed)
  d <- draws_t(r, n_sim = b$n_draws)
  expect_equal(d$method, "posterior")
  expect_equal(d$nb[, , 1], b$net_benefit, tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(d$nb[, , 2], b$treat_all, tolerance = 1e-14, ignore_attr = TRUE)
  expect_true(all(d$nb[, , 3] == 0))
  expect_equal(p_useful_from_draws(d)$p_useful, unname(b$p_useful))
})

test_that("EVSI equals the conjugate algorithm of Sadatsafavi et al. (2025)", {
  cs <- bayes_case()
  for (ref in cs$bref$evsi) {
    r <- nb(cs$fit, cs$val, "day30", thresholds = ref$z, bayesian = TRUE)
    set.seed(ref$seed)
    o <- evsi(r, future_n = ref$future_n, n_sim = ref$n_sim)
    expect_equal(attr(o, "method"), "posterior")
    # the expected net benefit of the strategy chosen after the future study is
    # the same; EVSI() subtracts the closed-form posterior mean of the current
    # best strategy, voipred the mean of the same draws
    set.seed(ref$seed)
    d <- draws_t(r, n_sim = ref$n_sim)
    cur_draws <- max(colMeans(d$nb[, 1, ]))
    cur_exact <- max(0, r$NB_model, r$NB_all)
    expect_equal(o$EVSI + cur_draws, ref$EVSI + cur_exact, tolerance = 1e-12)
    expect_equal(o$EVPI[1] + cur_draws, ref$EVPI + cur_exact, tolerance = 1e-12)
    expect_equal(o$EVPI[1], evpi_from_draws(d)$EVPI, tolerance = 1e-12)
  }
})

test_that("several models are drawn jointly, each keeping its own posterior", {
  cs <- bayes_case()
  z <- c(0.02, 0.05, 0.1, 0.2)
  f2 <- glm(day30 ~ age, family = binomial, data = gusto[1:5000, ])
  a <- nb(cs$fit, cs$val, "day30", thresholds = z, bayesian = TRUE, label = "Full")
  b <- nb(f2, cs$val, "day30", thresholds = z, bayesian = TRUE, label = "Age")
  # the Dirichlet prior over the patterns gives each model its Beta prior
  pr <- nb_beta(sensitivity = c(3, 2), specificity = c(1, 4))
  for (pp in list(nb_uniform(), nb_jeffreys(), pr)) {
    s2 <- bayes_setup(lapply(list(a, b, a), function(r) {
      v <- attr(r, "voi"); v$prior <- pp; attr(r, "voi") <- v; r }))
    for (m in 1:3) {
      expect_equal(sum(s2$alpha1[s2$pos[, m]]), pp$sensitivity[1])
      expect_equal(sum(s2$alpha1[!s2$pos[, m]]), pp$sensitivity[2])
      expect_equal(sum(s2$alpha0[!s2$pos[, m]]), pp$specificity[1])
      expect_equal(sum(s2$alpha0[s2$pos[, m]]), pp$specificity[2])
    }
  }
  # so the draws of each model have its posterior mean and spread
  set.seed(5)
  d <- draws_t(a, b, n_sim = 20000)
  expect_equal(colMeans(d$nb[, , 1]), a$NB_model, tolerance = 2e-3)
  expect_equal(colMeans(d$nb[, , 2]), b$NB_model, tolerance = 2e-3)
  set.seed(6)
  da <- draws_t(a, n_sim = 20000)
  expect_equal(apply(d$nb[, , 1], 2, sd), apply(da$nb[, , 1], 2, sd), tolerance = 0.05)
  # and the comparisons are coherent
  e <- evpi_from_draws(d)
  expect_true(all(e$EVPI >= 0))
  pb <- p_best_from_draws(d)
  expect_equal(as.vector(tapply(pb$p_best, pb$threshold, sum)), rep(1, 4))
  set.seed(7)
  o <- evsi(a, b, future_n = c(200, 2000, 20000), n_sim = 3000)
  expect_true(all(o$EVSI <= o$EVPI + 1e-12))
  expect_true(all(o$EVSI >= -1e-6))
  tot <- tapply(o$EVSI, o$future_n, sum)
  expect_true(all(diff(tot) >= -1e-6))
})

test_that("multinomial draws by conditional binomials", {
  set.seed(1)
  prob <- matrix(c(0.1, 0.2, 0.3, 0.4), 5000, 4, byrow = TRUE)
  x <- rmultinom_rows(rep(50, 5000), prob)
  expect_true(all(rowSums(x) == 50))
  expect_equal(colMeans(x) / 50, c(0.1, 0.2, 0.3, 0.4), tolerance = 0.02)
  expect_true(all(rmultinom_rows(c(0, 0), prob[1:2, ]) == 0))
})

test_that("Bayesian results in the value-of-information functions are checked", {
  cs <- bayes_case()
  z <- c(0.05, 0.1)
  r <- nb(cs$fit, cs$val, "day30", thresholds = z, bayesian = TRUE)
  f <- nb(cs$fit, cs$val, "day30", thresholds = z, label = "freq")
  expect_error(evpi(r, method = "ordinary"), "`method` is not used")
  expect_error(evpi(r, method = "asymptotic"), "`method` is not used")
  expect_error(p_useful(r, method = "bayesian"), "`method` is not used")
  expect_error(evsi(r, future_n = 100, method = "ordinary"), "`method` is not used")
  expect_error(draws_t(r, f), "different estimator \\(Bayesian or not\\)\\. ")
  r2 <- nb(cs$fit, cs$val, "day30", thresholds = z, bayesian = TRUE, label = "other",
           prior = nb_jeffreys())
  expect_error(draws_t(r, r2), "different prior")
  expect_error(evpi(r, r2), "different prior")
})
