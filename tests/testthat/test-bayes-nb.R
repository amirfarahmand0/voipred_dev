# The prior of nb(): by default (improper prior, Beta(0, 0)) the data only,
# which gives exactly the usual net benefit; with another prior (binary
# outcome) the net benefit is the posterior mean of the conjugate Beta model
# of bayesDCA.

test_that("with another prior, the net benefit is the posterior mean of the Beta model", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  val <- g$val
  y <- val$day30
  p <- unname(predict(f, val, type = "response"))
  z <- c(0, 0.02, 0.05, 0.2)
  pr <- nb_beta(prevalence = c(2, 30), sensitivity = c(3, 1), specificity = c(1, 2))
  r <- nb(val, p, "day30", thresholds = z, prior = pr)
  n <- length(y); d <- sum(y)
  ep <- (2 + d) / (32 + n)
  tp <- sapply(z, function(t) sum(y == 1 & p > 0 & p >= t))
  tn <- sapply(z, function(t) sum(y == 0 & !(p > 0 & p >= t)))
  se <- (3 + tp) / (4 + d)
  sp <- (1 + tn) / (3 + n - d)
  w <- z / (1 - z)
  expect_equal(r$estimates$NB_model, ep * se - (1 - ep) * (1 - sp) * w)
  expect_equal(r$estimates$NB_all, ep - (1 - ep) * w)
  expect_equal(r$estimates$delta_NB, r$estimates$NB_model - pmax(r$estimates$NB_all, 0))
  expect_equal(r$outcome$prevalence, ep)
  expect_identical(r$prior, pr)
  expect_equal(r$counts$TP, tp)
  expect_equal(r$counts$TN, tn)
  expect_equal(r$counts$TP + r$counts$FN, rep(d, length(z)))
  # uniform and Jeffreys priors, from the risks or from the model
  for (pp in list(nb_uniform(), nb_jeffreys())) {
    a <- pp$prevalence; s1 <- pp$sensitivity; s0 <- pp$specificity
    epp <- (a[1] + d) / (sum(a) + n)
    expect_equal(nb(val, f, "day30", thresholds = z, prior = pp)$estimates$NB_model,
                 epp * (s1[1] + tp) / (sum(s1) + d) -
                   (1 - epp) * (1 - (s0[1] + tn) / (sum(s0) + n - d)) * w)
  }
  # a large sample: close to the frequentist net benefit
  fr <- nb(val, p, "day30", thresholds = z)
  expect_lt(max(abs(nb(val, p, "day30", thresholds = z, prior = nb_uniform())$estimates$NB_model -
                      fr$estimates$NB_model)), 0.005)
  expect_s3_class(plot(r), "ggplot")
  expect_error(plot(r, fr, labels = c("a", "b")), "different prior")
  expect_error(plot(r, b = nb(val, p, "day30", thresholds = z, prior = nb_jeffreys())),
               "different prior")
})

test_that("the improper prior (default) gives exactly the frequentist net benefit", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  val <- g$val
  p <- unname(predict(f, val, type = "response"))
  z <- all_z
  r <- nb(val, p, "day30", thresholds = z)
  expect_identical(r$prior, nb_improper())
  freq <- nb_core(val$day30, p, z)                      # TP/n - FP/n z/(1-z)
  expect_equal(r$estimates$NB_model, freq$NB_model, tolerance = 1e-14)
  expect_equal(r$estimates$NB_all, freq$NB_all, tolerance = 1e-14)
  expect_equal(r$outcome$prevalence, mean(val$day30))
  expect_equal(r$outcome$events, sum(val$day30))
  # with case weights: the weighted frequentist net benefit
  set.seed(3)
  w <- rpois(nrow(val), 1)
  rw <- nb(val, p, "day30", thresholds = z, weights = w)
  fw <- nb_core(val$day30, p, z, w)
  expect_equal(rw$estimates$NB_model, fw$NB_model, tolerance = 1e-14)
  dup <- nb(val[rep(seq_len(nrow(val)), w), ], p[rep(seq_len(nrow(val)), w)], "day30",
            thresholds = z)
  expect_equal(rw$estimates$NB_model, dup$estimates$NB_model, tolerance = 1e-12)
})

test_that("the improper prior handles empty groups (counts of 0) as the frequentist estimate", {
  d <- data.frame(y = c(0, 0, 0, 0), p = c(0.1, 0.3, 0.5, 0.7))
  r <- nb(d, d$p, "y", thresholds = c(0, 0.2, 0.6, 0.9))    # no events at all
  expect_equal(r$estimates$NB_model, nb_core(d$y, d$p, c(0, 0.2, 0.6, 0.9))$NB_model)
  expect_true(all(is.finite(r$estimates$NB_model)))
  d1 <- data.frame(y = c(1, 1, 1), p = c(0.2, 0.5, 0.9))   # only events
  r1 <- nb(d1, d1$p, "y", thresholds = c(0.1, 0.6))
  expect_equal(r1$estimates$NB_model, c(1, 1 / 3))
  # a threshold where the model treats everyone (no false negatives) or no one
  d2 <- data.frame(y = c(1, 0, 1, 0, 0), p = c(0.9, 0.2, 0.4, 0.05, 0.6))
  r2 <- nb(d2, d2$p, "y", thresholds = c(0.01, 0.95))
  expect_equal(r2$estimates$NB_model, nb_core(d2$y, d2$p, c(0.01, 0.95))$NB_model)
  # the posterior draws are finite too
  set.seed(1)
  expect_true(all(is.finite(draws_t(r, n_sim = 50)$nb)))
  expect_true(all(is.finite(draws_t(r1, n_sim = 50)$nb)))
  expect_true(all(is.finite(draws_t(r2, n_sim = 50)$nb)))
  set.seed(1)
  o <- evsi(r2, future_n = c(5, 50), n_sim = 50)
  expect_true(all(is.finite(o$estimates$EVSI)))
})

test_that("priors are checked", {
  g <- gusto_split()
  val <- g$val
  p <- unname(predict(glm(gusto_formula, family = binomial, data = g$dev), val,
                      type = "response"))
  expect_error(nb_beta(prevalence = c(1, -1)), "two non-negative numbers")
  expect_error(nb_beta(sensitivity = 1), "two non-negative numbers")
  expect_s3_class(nb_uniform(), "nb_prior")
  expect_identical(nb_uniform(), nb_beta(c(1, 1), c(1, 1), c(1, 1)))
  expect_identical(nb_improper(), nb_beta(c(0, 0), c(0, 0), c(0, 0)))
  expect_identical(nb_jeffreys(), nb_beta(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5)))
  expect_equal(nb_improper()$type, "improper")
  expect_equal(nb_beta(prevalence = c(7, 93))$type, "beta")
  expect_output(print(nb_improper()), "^Prior: improper \\(Haldane\\), Beta\\(0, 0\\)")
  expect_output(print(nb_uniform()), "^Prior: uniform, Beta\\(1, 1\\)")
  expect_output(print(nb_jeffreys()), "^Prior: Jeffreys, Beta\\(0.5, 0.5\\)")
  expect_output(print(nb_beta(prevalence = c(7, 93))),
                "Prior: Beta: prevalence Beta\\(7, 93\\), sensitivity Beta\\(1, 1\\)")
  expect_error(nb(val, p, "day30", prior = c(1, 1)),
               "nb_improper\\(\\), nb_uniform\\(\\), nb_jeffreys\\(\\) or nb_beta")
  expect_error(nb(val, p, "day30", prior = nb_uniform(), weights = rep(1, nrow(val))),
               "Case weights can only be used with the improper prior")
  expect_error(nb(val, p, "day30", bayesian = TRUE), "unused argument")
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  expect_error(nb(s$val, cx, surv_outcome, horizon = 5, prior = nb_uniform()),
               "the prior is nb_improper\\(\\)")
  expect_no_error(quiet_nb(s$val, cx, surv_outcome, horizon = 5, prior = nb_improper()))
})
