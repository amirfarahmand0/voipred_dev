toy <- data.frame(y = c(1, 0, 1, 0, 0), p = c(0.9, 0.2, 0.4, 0.05, 0.6))

test_that("nb() matches a hand calculation", {
  res <- nb(fixed_model(toy$p), toy, "y", thresholds = 0.3)
  # treated: risks 0.9 (event), 0.4 (event), 0.6 (no event)
  expect_equal(res$NB_model, 2 / 5 - 1 / 5 * 0.3 / 0.7)
  expect_equal(res$NB_all, 2 / 5 - 3 / 5 * 0.3 / 0.7)
  expect_equal(res$NB_none, 0)
  expect_equal(attr(res, "prevalence"), 0.4)
  expect_equal(attr(res, "n"), 5)
})

test_that("nb() returns one row per threshold with the expected columns", {
  z <- c(0, 0.1, 0.5, 0.95)
  res <- nb(fixed_model(toy$p), toy, "y", thresholds = z)
  expect_s3_class(res, "data.frame")
  expect_s3_class(res, "voipred_nb")
  expect_named(res, c("threshold", "NB_model", "NB_all", "NB_none"))
  expect_equal(res$threshold, z)
  expect_true(all(res$NB_none == 0))
})

test_that("at threshold 0 the model and treat-all both equal the prevalence", {
  res <- nb(fixed_model(toy$p), toy, "y", thresholds = 0)
  expect_equal(res$NB_model, mean(toy$y))
  expect_equal(res$NB_all, mean(toy$y))
})

test_that("a risk exactly at the threshold is treated (p >= z, as dcurves)", {
  d <- data.frame(y = c(1, 0), p = c(0.2, 0.2))
  res <- nb(fixed_model(d$p), d, "y", thresholds = 0.2)
  expect_equal(res$NB_model, 1 / 2 - 1 / 2 * 0.25)
})

test_that("a risk of exactly 0 is never treated, even at threshold 0 (as dcurves)", {
  d <- data.frame(y = c(1, 1, 0), p = c(0.5, 0, 0))
  res <- nb(fixed_model(d$p), d, "y", thresholds = 0)
  expect_equal(res$NB_model, 1 / 3)   # the event with risk 0 is not treated
  expect_equal(res$NB_all, 2 / 3)
})

test_that("unsorted thresholds are kept in the order given", {
  z <- c(0.5, 0.1, 0.3)
  res <- nb(fixed_model(toy$p), toy, "y", thresholds = z)
  ref <- sapply(z, function(t) nb(fixed_model(toy$p), toy, "y", thresholds = t)$NB_model)
  expect_equal(res$threshold, z)
  expect_equal(res$NB_model, ref)
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
  res <- nb(fixed_model(p), data.frame(y = y), "y", thresholds = z)
  expect_equal(res$NB_model, slow, tolerance = 1e-12)
})

test_that("outcome can be a column name, a vector, logical, or a factor", {
  ref <- nb(fixed_model(toy$p), toy, "y", thresholds = 0.3)
  expect_equal(nb(fixed_model(toy$p), toy, toy$y, thresholds = 0.3)$NB_model, ref$NB_model)
  expect_equal(nb(fixed_model(toy$p), toy, toy$y == 1, thresholds = 0.3)$NB_model, ref$NB_model)
  f <- factor(ifelse(toy$y == 1, "dead", "alive"), levels = c("alive", "dead"))
  expect_equal(nb(fixed_model(toy$p), toy, f, thresholds = 0.3)$NB_model, ref$NB_model)
  f01 <- factor(toy$y, levels = c(1, 0))          # level "1" is the event even if first
  expect_equal(nb(fixed_model(toy$p), toy, f01, thresholds = 0.3)$NB_model, ref$NB_model)
})

test_that("weights: all ones equal no weights, integer weights equal duplicated rows", {
  z <- seq(0, 0.9, by = 0.1)
  expect_equal(nb(fixed_model(toy$p), toy, "y", z, weights = rep(1, 5))$NB_model,
               nb(fixed_model(toy$p), toy, "y", z)$NB_model)
  w <- c(2, 0, 1, 3, 1)
  dup <- toy[rep(seq_len(5), w), ]
  expect_equal(nb(fixed_model(toy$p), toy, "y", z, weights = w)$NB_model,
               nb(fixed_model(dup$p), dup, "y", z)$NB_model)
  expect_equal(nb(fixed_model(toy$p), toy, "y", z, weights = w)$NB_all,
               nb(fixed_model(dup$p), dup, "y", z)$NB_all)
})

test_that("bad inputs give clear errors", {
  expect_error(nb(fixed_model(toy$p), toy, "nope"), "not in `data`")
  expect_error(nb(fixed_model(toy$p), toy, "y", thresholds = 1), "\\[0, 1\\)")
  expect_error(nb(fixed_model(toy$p), toy, "y", thresholds = -0.1), "\\[0, 1\\)")
  expect_error(nb(fixed_model(toy$p[-1]), toy, "y"), "returned 4 predictions for 5 rows")
  expect_error(nb(fixed_model(toy$p), toy, c(1, 2, 0, 0, 1)), "binary")
  expect_error(nb(fixed_model(toy$p), toy, c(1, NA, 0, 0, 1)), "missing")
  expect_error(nb(fixed_model(c(0.9, 0.2, 1.4, 0.05, 0.6)), toy, "y"), "\\[0, 1\\]")
  expect_error(nb(fixed_model(toy$p), as.list(toy), "y"), "data frame")
  expect_error(nb(fixed_model(toy$p), toy, "y", weights = c(1, 1, -1, 1, 1)), "non-negative")
})

test_that("a vector of predicted risks is rejected: nb() needs the model", {
  expect_error(nb(toy$p, toy, "y"), "fitted model, not a vector")
  expect_error(nb(as.integer(toy$y), toy, "y"), "fitted model, not a vector")
})

test_that("print shows a header and the table", {
  res <- nb(fixed_model(toy$p), toy, "y", thresholds = c(0.1, 0.3))
  expect_output(print(res), "Net benefit: n = 5")
  expect_output(print(res), "NB_model")
})
