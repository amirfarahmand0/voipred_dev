# The value of information of a binary outcome: draws from the conjugate
# posterior of nb_val()'s prior; the point estimates are in test-bayes-nb.R.
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
  r <- nb_val(cs$val, cs$fit, "day30", thresholds = b$zs, prior = nb_uniform())
  set.seed(b$seed)
  d <- draws_t(r, n_sim = b$n_draws)
  expect_equal(d$draws, "posterior")
  expect_equal(d$nb[, , 1], b$net_benefit, tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(d$nb[, , 2], b$treat_all, tolerance = 1e-14, ignore_attr = TRUE)
  expect_true(all(d$nb[, , 3] == 0))
  expect_equal(p_useful_from_draws(d)$estimates$p_useful, unname(b$p_useful))
})

test_that("EVSI equals the conjugate algorithm of Sadatsafavi et al. (2025)", {
  cs <- bayes_case()
  for (ref in cs$bref$evsi) {
    r <- nb_val(cs$val, cs$fit, "day30", thresholds = ref$z, prior = nb_uniform())
    set.seed(ref$seed)
    o <- evsi(r, future_n = ref$future_n, n_sim = ref$n_sim)
    expect_equal(o$draws, "posterior")
    # the expected net benefit of the strategy chosen after the future study is
    # the same; EVSI() subtracts the closed-form posterior mean of the current
    # best strategy, voipred the mean of the same draws
    set.seed(ref$seed)
    d <- draws_t(r, n_sim = ref$n_sim)
    cur_draws <- max(colMeans(d$nb[, 1, ]))
    cur_exact <- max(0, r$estimates$NB_model, r$estimates$NB_all)
    expect_equal(o$estimates$EVSI + cur_draws, ref$EVSI + cur_exact, tolerance = 1e-12)
    expect_equal(o$estimates$EVPI[1] + cur_draws, ref$EVPI + cur_exact, tolerance = 1e-12)
    expect_equal(o$estimates$EVPI[1], evpi_from_draws(d)$estimates$EVPI, tolerance = 1e-12)
  }
})

test_that("several models are drawn jointly, each keeping its own posterior", {
  cs <- bayes_case()
  z <- c(0.02, 0.05, 0.1, 0.2)
  f2 <- glm(day30 ~ age, family = binomial, data = cs$fit$data)
  a <- nb_val(cs$val, cs$fit, "day30", thresholds = z, prior = nb_uniform(), label = "Full")
  b <- nb_val(cs$val, f2, "day30", thresholds = z, prior = nb_uniform(), label = "Age")
  # the Dirichlet prior over the patterns gives each model its Beta prior
  pr <- nb_beta(sensitivity = c(3, 2), specificity = c(1, 4))
  for (pp in list(nb_uniform(), nb_jeffreys(), pr)) {
    rs <- lapply(list(a, b, a), function(r) {
      r$prior <- pp
      r
    })
    names(rs) <- c("", "", "Full again")
    s2 <- bayes_setup(do.call(c, rs))
    # the Dirichlet prior has one row per threshold (the same for these priors)
    for (j in seq_len(nrow(s2$alpha1))) for (m in 1:3) {
      expect_equal(sum(s2$alpha1[j, s2$pos[, m]]), pp$sensitivity[1])
      expect_equal(sum(s2$alpha1[j, !s2$pos[, m]]), pp$sensitivity[2])
      expect_equal(sum(s2$alpha0[j, !s2$pos[, m]]), pp$specificity[1])
      expect_equal(sum(s2$alpha0[j, s2$pos[, m]]), pp$specificity[2])
    }
  }
  # so the draws of each model have its posterior mean and spread
  set.seed(5)
  d <- draws_t(a, b, n_sim = 20000)
  expect_equal(colMeans(d$nb[, , 1]), a$estimates$NB_model, tolerance = 2e-3)
  expect_equal(colMeans(d$nb[, , 2]), b$estimates$NB_model, tolerance = 2e-3)
  set.seed(6)
  da <- draws_t(a, n_sim = 20000)
  expect_equal(apply(d$nb[, , 1], 2, sd), apply(da$nb[, , 1], 2, sd), tolerance = 0.05)
  # and the comparisons are coherent
  e <- evpi_from_draws(d)
  expect_true(all(e$estimates$EVPI >= 0))
  pb <- p_best_from_draws(d)
  expect_equal(as.vector(tapply(pb$estimates$p_best, pb$estimates$threshold, sum)), rep(1, 4))
  set.seed(7)
  o <- evsi(a, b, future_n = c(200, 2000, 20000), n_sim = 3000)
  expect_true(all(o$estimates$EVSI <= o$estimates$EVPI + 1e-12))
  expect_true(all(o$estimates$EVSI >= -1e-6))
  tot <- tapply(o$estimates$EVSI, o$estimates$future_n, sum)
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

test_that("results with different priors are not compared; `method` is gone", {
  cs <- bayes_case()
  z <- c(0.05, 0.1)
  r <- nb_val(cs$val, cs$fit, "day30", thresholds = z, prior = nb_uniform())
  f <- nb_val(cs$val, cs$fit, "day30", thresholds = z, label = "improper")
  expect_error(draws_t(r, f), "different prior")
  r2 <- nb_val(cs$val, cs$fit, "day30", thresholds = z, label = "other", prior = nb_jeffreys())
  expect_error(draws_t(r, r2), "different prior")
  expect_error(evpi(r, r2), "different prior")
  expect_error(evpi(r, method = "ordinary"), "Unknown argument `method`")
  expect_error(p_useful(r, method = "bayesian"), "Unknown argument `method`")
  expect_error(evsi(r, future_n = 100, method = "ordinary"), "Unknown argument `method`")
})

test_that("with the improper prior the posterior draws equal the Bayesian bootstrap in distribution", {
  cs <- bayes_case()
  z <- c(0.017, 0.023, 0.05, 0.1)
  f2 <- glm(day30 ~ age, family = binomial, data = cs$fit$data)
  a <- nb_val(cs$val, cs$fit, "day30", thresholds = z, label = "Full")
  b <- nb_val(cs$val, f2, "day30", thresholds = z, label = "Age")
  set.seed(11)
  post <- draws_t(a, b, n_sim = 20000)              # conjugate posterior, improper prior
  expect_equal(post$draws, "posterior")
  set.seed(12)
  boot <- bb_draws_t(a, b, n_sim = 20000)           # Bayesian bootstrap
  for (s in 1:3) {
    se <- apply(boot$nb[, , s], 2, sd)
    expect_lt(max(abs(colMeans(post$nb[, , s]) - colMeans(boot$nb[, , s])) / (se / sqrt(20000))),
              4.5)                                  # same means (Monte Carlo error)
    expect_equal(apply(post$nb[, , s], 2, sd), se, tolerance = 0.03)   # same spread
  }
  # the same correlation between the models
  expect_equal(diag(cor(post$nb[, , 1], post$nb[, , 2])),
               diag(cor(boot$nb[, , 1], boot$nb[, , 2])), tolerance = 0.03)
  # so the same value of information
  ep <- evpi_from_draws(post)$estimates
  eb <- evpi_from_draws(boot)$estimates
  expect_equal(ep$EVPI, eb$EVPI, tolerance = 0.08)
  expect_equal(p_best_from_draws(post)$estimates$p_best,
               p_best_from_draws(boot)$estimates$p_best, tolerance = 0.03)
  # and EVSI: the conjugate algorithm against the bootstrap algorithm
  set.seed(13)
  oc <- evsi(a, future_n = c(500, 2000), n_sim = 20000)$estimates
  set.seed(14)
  ob <- evsi_bb(a, future_n = c(500, 2000), n_sim = 20000)
  expect_equal(oc$EVSI, as.vector(ob$EVSI), tolerance = 0.1)
})
