toy <- data.frame(y = c(1, 0, 1, 0, 0), p = c(0.9, 0.2, 0.4, 0.05, 0.6))

test_that("nb_val() matches a hand calculation", {
  res <- nb_val(toy, fixed_model(toy$p), "y", thresholds = 0.3)
  # treated: risks 0.9 (event), 0.4 (event), 0.6 (no event)
  expect_equal(res$estimates$NB_model, 2 / 5 - 1 / 5 * 0.3 / 0.7)
  expect_equal(res$estimates$NB_all, 2 / 5 - 3 / 5 * 0.3 / 0.7)
  expect_equal(res$estimates$NB_none, 0)
  expect_equal(res$outcome$prevalence, 0.4)
  expect_equal(res$outcome$n, 5)
})

test_that("nb_val() returns one row per threshold with the expected columns", {
  z <- c(0, 0.1, 0.5, 0.95)
  res <- nb_val(toy, fixed_model(toy$p), "y", thresholds = z)
  expect_s3_class(res, "voipred_nb")
  expect_s3_class(res$estimates, "data.frame")
  expect_named(res$estimates, c("threshold", "NB_model", "NB_all", "NB_none", "delta_NB"))
  expect_identical(as.data.frame(res), res$estimates)
  expect_equal(res$estimates$threshold, z)
  expect_true(all(res$estimates$NB_none == 0))
})

test_that("at threshold 0 the model and treat-all both equal the prevalence", {
  res <- nb_val(toy, fixed_model(toy$p), "y", thresholds = 0)
  expect_equal(res$estimates$NB_model, mean(toy$y))
  expect_equal(res$estimates$NB_all, mean(toy$y))
})

test_that("a risk exactly at the threshold is treated (p >= z, as dcurves)", {
  d <- data.frame(y = c(1, 0), p = c(0.2, 0.2))
  res <- nb_val(d, fixed_model(d$p), "y", thresholds = 0.2)
  expect_equal(res$estimates$NB_model, 1 / 2 - 1 / 2 * 0.25)
})

test_that("a risk of exactly 0 is never treated, even at threshold 0 (as dcurves)", {
  d <- data.frame(y = c(1, 1, 0), p = c(0.5, 0, 0))
  res <- nb_val(d, fixed_model(d$p), "y", thresholds = 0)
  expect_equal(res$estimates$NB_model, 1 / 3)   # the event with risk 0 is not treated
  expect_equal(res$estimates$NB_all, 2 / 3)
})

test_that("unsorted thresholds are kept in the order given", {
  z <- c(0.5, 0.1, 0.3)
  res <- nb_val(toy, fixed_model(toy$p), "y", thresholds = z)
  ref <- sapply(z, function(t) nb_val(toy, fixed_model(toy$p), "y", thresholds = t)$estimates$NB_model)
  expect_equal(res$estimates$threshold, z)
  expect_equal(res$estimates$NB_model, ref)
})

test_that("the fast core agrees with a direct loop on random data", {
  set.seed(1)
  n <- 3000
  p <- round(runif(n), 2)          # many ties on the threshold grid
  y <- rbinom(n, 1, p)
  z <- seq(0, 0.99, by = 0.01)
  y[1:20] <- 1; p[1:20] <- 0       # events with risk 0
  slow <- sapply(z, function(t) {
    treat <- p >= t & p > 0
    mean(y * treat) - mean((1 - y) * treat) * t / (1 - t)
  })
  res <- nb_val(data.frame(y = y), fixed_model(p), "y", thresholds = z)
  expect_equal(res$estimates$NB_model, slow, tolerance = 1e-12)
})

test_that("outcome can be a column name, a vector, logical, or a factor", {
  ref <- nb_val(toy, fixed_model(toy$p), "y", thresholds = 0.3)$estimates
  expect_equal(nb_val(toy, fixed_model(toy$p), toy$y, thresholds = 0.3)$estimates$NB_model, ref$NB_model)
  expect_equal(nb_val(toy, fixed_model(toy$p), toy$y == 1, thresholds = 0.3)$estimates$NB_model, ref$NB_model)
  f <- factor(ifelse(toy$y == 1, "dead", "alive"), levels = c("alive", "dead"))
  expect_equal(nb_val(toy, fixed_model(toy$p), f, thresholds = 0.3)$estimates$NB_model, ref$NB_model)
  f01 <- factor(toy$y, levels = c(1, 0))          # level "1" is the event even if first
  expect_equal(nb_val(toy, fixed_model(toy$p), f01, thresholds = 0.3)$estimates$NB_model, ref$NB_model)
})

test_that("weights: all ones equal no weights, integer weights equal duplicated rows", {
  z <- seq(0, 0.9, by = 0.1)
  expect_equal(nb_val(toy, fixed_model(toy$p), "y", z, weights = rep(1, 5))$estimates$NB_model,
               nb_val(toy, fixed_model(toy$p), "y", z)$estimates$NB_model)
  w <- c(2, 0, 1, 3, 1)
  dup <- toy[rep(seq_len(5), w), ]
  expect_equal(nb_val(toy, fixed_model(toy$p), "y", z, weights = w)$estimates$NB_model,
               nb_val(dup, fixed_model(dup$p), "y", z)$estimates$NB_model)
  expect_equal(nb_val(toy, fixed_model(toy$p), "y", z, weights = w)$estimates$NB_all,
               nb_val(dup, fixed_model(dup$p), "y", z)$estimates$NB_all)
})

test_that("bad inputs give clear errors", {
  expect_error(nb_val(toy, fixed_model(toy$p), "nope"), "not in `data`")
  expect_error(nb_val(toy, fixed_model(toy$p), "y", thresholds = 1), "\\[0, 1\\)")
  expect_error(nb_val(toy, fixed_model(toy$p), "y", thresholds = -0.1), "\\[0, 1\\)")
  expect_error(nb_val(toy, fixed_model(toy$p[-1]), "y"), "returned 4 predictions for 5 rows")
  expect_error(nb_val(toy, fixed_model(toy$p), c(1, 2, 0, 0, 1)), "binary")
  expect_error(nb_val(toy, fixed_model(toy$p), c(1, NA, 0, 0, 1)), "missing")
  expect_error(nb_val(toy, fixed_model(c(0.9, 0.2, 1.4, 0.05, 0.6)), "y"), "\\[0, 1\\]")
  expect_error(nb_val(as.list(toy), fixed_model(toy$p), "y"), "data frame")
  expect_error(nb_val(toy, fixed_model(toy$p), "y", weights = c(1, 1, -1, 1, 1)), "non-negative")
})

test_that("the predicted risks can be given directly, without the model", {
  ref <- nb_val(toy, fixed_model(toy$p), "y", thresholds = c(0.1, 0.5))$estimates
  expect_equal(nb_val(toy, toy$p, "y", thresholds = c(0.1, 0.5))$estimates$NB_model, ref$NB_model)
  expect_equal(nb_val(toy, "p", "y", thresholds = c(0.1, 0.5))$estimates$NB_model, ref$NB_model)
  expect_equal(nb_val(predicted_risk = toy$p, data = toy, outcome = "y",
                      thresholds = c(0.1, 0.5))$estimates$NB_model, ref$NB_model)
  expect_error(nb_val(toy, as.integer(toy$y), "y"), "only 0s and 1s")
})

test_that("print shows a header and the table", {
  res <- nb_val(toy, fixed_model(toy$p), "y", thresholds = c(0.1, 0.3))
  expect_output(print(res), "Net benefit of Model: binary, n = 5, 2 events \\(40.0%\\), prior improper Beta\\(0, 0\\)")
  expect_output(print(res), "NB_model")
})

test_that("the first argument is data; the old order gives a clear error", {
  expect_error(nb_val(fixed_model(toy$p), toy, "y"), "The first argument of nb\\(\\) is `data`")
  expect_error(nb_val(toy$p, toy, "y"), "The first argument of nb\\(\\) is `data`")
  res <- nb_val(toy, toy$p, "y", thresholds = 0.3)
  expect_equal(res$call[[1]], as.name("nb"))
  expect_equal(res$label, "p")
})
