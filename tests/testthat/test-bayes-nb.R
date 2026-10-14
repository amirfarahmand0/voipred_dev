# The Bayesian net benefit (nb(bayesian = TRUE)): posterior means of the
# conjugate Beta model of bayesDCA, binary outcomes only.

test_that("the net benefit is the posterior mean of the Beta model", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  val <- g$val
  y <- val$day30
  p <- unname(predict(f, val, type = "response"))
  z <- c(0, 0.02, 0.05, 0.2)
  pr <- nb_beta(prevalence = c(2, 30), sensitivity = c(3, 1), specificity = c(1, 2))
  r <- nb(p, val, "day30", thresholds = z, bayesian = TRUE, prior = pr)
  n <- length(y); d <- sum(y)
  ep <- (2 + d) / (32 + n)
  tp <- sapply(z, function(t) sum(y == 1 & p > 0 & p >= t))
  tn <- sapply(z, function(t) sum(y == 0 & !(p > 0 & p >= t)))
  se <- (3 + tp) / (4 + d)
  sp <- (1 + tn) / (3 + n - d)
  w <- z / (1 - z)
  expect_equal(r$NB_model, ep * se - (1 - ep) * (1 - sp) * w)
  expect_equal(r$NB_all, ep - (1 - ep) * w)
  expect_equal(r$delta_NB, r$NB_model - pmax(r$NB_all, 0))
  expect_equal(attr(r, "prevalence"), ep)
  expect_true(attr(r, "bayesian"))
  expect_identical(attr(r, "prior"), pr)
  expect_output(print(r), "Bayesian, posterior means")
  # uniform and Jeffreys priors, from the risks or from the model
  for (pp in list(nb_uniform(), nb_jeffreys())) {
    a <- pp$prevalence; s1 <- pp$sensitivity; s0 <- pp$specificity
    epp <- (a[1] + d) / (sum(a) + n)
    expect_equal(nb(f, val, "day30", thresholds = z, bayesian = TRUE, prior = pp)$NB_model,
                 epp * (s1[1] + tp) / (sum(s1) + d) -
                   (1 - epp) * (1 - (s0[1] + tn) / (sum(s0) + n - d)) * w)
  }
  # a large sample: close to the frequentist net benefit
  fr <- nb(p, val, "day30", thresholds = z)
  expect_false(attr(fr, "bayesian"))
  expect_lt(max(abs(nb(p, val, "day30", thresholds = z, bayesian = TRUE)$NB_model -
                      fr$NB_model)), 0.005)
  expect_s3_class(plot_dca(r), "ggplot")
  expect_error(plot_dca(r, fr, labels = c("a", "b")), "different estimator \\(Bayesian or not\\)")
  expect_error(plot_dca(r, b = nb(p, val, "day30", thresholds = z, bayesian = TRUE,
                                  prior = nb_jeffreys())), "different prior")
})

test_that("priors and Bayesian settings are checked", {
  g <- gusto_split()
  val <- g$val
  p <- unname(predict(glm(gusto_formula, family = binomial, data = g$dev), val,
                      type = "response"))
  expect_error(nb_beta(prevalence = c(1, 0)), "two positive numbers")
  expect_error(nb_beta(sensitivity = 1), "two positive numbers")
  expect_s3_class(nb_uniform(), c("prior_uniform", "nb_prior"))
  expect_equal(nb_jeffreys()$sensitivity, c(0.5, 0.5))
  expect_output(print(nb_uniform()), "^Prior: uniform, Beta\\(1, 1\\)")
  expect_output(print(nb_jeffreys()), "^Prior: Jeffreys, Beta\\(0.5, 0.5\\)")
  expect_output(print(nb_beta(prevalence = c(7, 93))),
                "Prior: Beta: prevalence Beta\\(7, 93\\), sensitivity Beta\\(1, 1\\)")
  expect_error(nb(p, val, "day30", bayesian = NA), "TRUE or FALSE")
  expect_error(nb(p, val, "day30", bayesian = TRUE, prior = c(1, 1)),
               "nb_uniform\\(\\), nb_jeffreys\\(\\) or nb_beta")
  expect_error(nb(p, val, "day30", prior = nb_jeffreys()), "only used with `bayesian = TRUE`")
  expect_error(nb(p, val, "day30", bayesian = TRUE, weights = rep(1, nrow(val))),
               "Case weights cannot be used")
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  expect_error(nb(cx, s$val, surv_outcome, horizon = 5, bayesian = TRUE),
               "binary outcomes only")
})
