# Credible intervals of the net benefit: nb_val(interval = ), from the draws of
# the value-of-information functions (posterior for a binary outcome, Bayesian
# bootstrap otherwise), shown by print(), summary(), as.data.frame() and plot().

ci_case <- function(interval = 0.95, n_sim = 400, prior = nb_improper()) {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  v <- g$val[1:500, ]
  v$p1 <- unname(predict(f1, v, type = "response"))
  v$p2 <- unname(predict(f2, v, type = "response"))
  set.seed(7)
  res <- suppressWarnings(nb_val(v, list(Full = "p1", Age = "p2"), "day30",
                                 thresholds = seq(0.01, 0.3, by = 0.01),
                                 interval = interval, n_sim = n_sim, prior = prior))
  list(res = res, v = v)
}

test_that("binary: the quantiles of the posterior draws (the draws of evpi())", {
  cs <- ci_case()
  res <- cs$res
  expect_equal(res$interval, list(level = 0.95, n_sim = 400, draws = "posterior"))
  expect_named(res$estimates, c("threshold", "model", "NB_model", "NB_model_lower",
                                "NB_model_upper", "NB_all", "NB_all_lower", "NB_all_upper",
                                "NB_none", "delta_NB", "delta_NB_lower", "delta_NB_upper",
                                "net_FP_avoided", "net_FP_avoided_lower", "net_FP_avoided_upper"))
  # the same draws as nb_draws() with the same seed
  base <- res
  base$estimates <- base$estimates[, c("threshold", "model", "NB_model", "NB_all", "NB_none",
                                       "delta_NB")]
  set.seed(7)
  d <- nb_draws(base, 400)
  all_d <- d$nb[, , 3]
  full <- res$estimates[res$estimates$model == "Full", ]
  expect_equal(full$NB_model_lower, apply(d$nb[, , 1], 2, quantile, 0.025, names = FALSE))
  expect_equal(full$NB_model_upper, apply(d$nb[, , 1], 2, quantile, 0.975, names = FALSE))
  expect_equal(full$delta_NB_lower,
               apply(d$nb[, , 1] - pmax(all_d, 0), 2, quantile, 0.025, names = FALSE))
  expect_equal(full$NB_all_upper, apply(all_d, 2, quantile, 0.975, names = FALSE))
  # the estimates themselves are unchanged
  r0 <- nb_val(cs$v, list(Full = "p1", Age = "p2"), "day30", thresholds = seq(0.01, 0.3, by = 0.01))
  expect_equal(res$estimates[names(r0$estimates)], r0$estimates)
  expect_null(r0$interval)
  expect_true(all(res$estimates$NB_model_lower <= res$estimates$NB_model_upper))
  # a narrower level gives narrower intervals
  r80 <- ci_case(interval = 0.8)$res
  expect_true(all(r80$estimates$NB_model_upper - r80$estimates$NB_model_lower <=
                    res$estimates$NB_model_upper - res$estimates$NB_model_lower + 1e-12))
})

test_that("the prior is used (e.g. nb_informative())", {
  inf <- ci_case(prior = nb_informative())$res
  imp <- ci_case()$res
  expect_equal(inf$interval$draws, "posterior")
  expect_false(isTRUE(all.equal(inf$estimates$NB_model_lower, imp$estimates$NB_model_lower)))
})

test_that("survival: the Bayesian bootstrap (Kaplan-Meier and cIPCW)", {
  s <- rott_split()
  v <- s$val
  v$p <- pmin(0.9, 0.1 + v$nodes / 40)
  z <- c(0.1, 0.3)
  set.seed(3)
  km <- suppressWarnings(nb_val(v, "p", surv_outcome, horizon = 5, thresholds = z,
                                interval = 0.9, n_sim = 50))
  expect_equal(km$interval$draws, "Bayesian bootstrap")
  expect_true(all(km$estimates$NB_model_lower <= km$estimates$NB_model &
                    km$estimates$NB_model <= km$estimates$NB_model_upper))
  set.seed(3)
  ci <- suppressWarnings(nb_val(v, "p", surv_outcome, horizon = 5, thresholds = z,
                                cipcw = ~ age + nodes, interval = 0.9, n_sim = 20))
  expect_equal(ci$interval$n_sim, 20)
  expect_output(suppressWarnings(print(ci)), "90% credible intervals from 20 Bayesian bootstrap draws")
})

test_that("print(), summary() and as.data.frame() show the intervals", {
  res <- ci_case()$res
  out <- capture.output(print(res, n = 1))
  expect_match(out[1], "; 95% credible intervals from 400 posterior draws$")
  expect_true(any(grepl("95% CrI", out)))
  s <- summary(res, thresholds = c(0.05, 0.1))
  expect_true(all(c("95% CrI", "95% CrI 1") %in% names(s$table)) ||
                sum(names(s$table) == "95% CrI") == 2)
  ci_nb <- s$table[[which(names(s$table) == "95% CrI")[1]]]
  expect_match(ci_nb[1], "^\\(-?[0-9.]+, -?[0-9.]+\\)$")
  # scaled with per
  sp <- summary(res, thresholds = 0.05, per = 1e5)
  e <- res$estimates[res$estimates$model == "Full" & round(res$estimates$threshold, 2) == 0.05, ]
  expect_match(sp$table[[which(names(sp$table) == "95% CrI")[1]]][1],
               paste0("^\\(", round(e$NB_model_lower * 1e5), ", "))
  d <- as.data.frame(res, per = 1e5)
  expect_equal(d$NB_model_lower_per_100000, res$estimates$NB_model_lower * 1e5)
  # without intervals: as before
  r0 <- ci_case(interval = NULL)$res
  expect_false(any(grepl("CrI", capture.output(print(r0)))))
})

test_that("plot(): bands, bands = FALSE, per, and interactive hover", {
  res <- ci_case()$res
  layers <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
  # bands only when asked for (default: none)
  expect_false("GeomRibbon" %in% layers(plot(res)))
  expect_true("GeomRibbon" %in% layers(suppressWarnings(plot(res, bands = TRUE))))
  expect_true("GeomRibbon" %in% layers(suppressWarnings(plot(res, type = "delta", bands = TRUE))))
  expect_error(plot(ci_case(interval = NULL)$res, bands = TRUE),
               "`bands = TRUE` needs credible intervals")
  rib <- function(p) p$layers[[which(layers(p) == "GeomRibbon")]]$data
  b1 <- rib(suppressWarnings(plot(res, bands = TRUE)))
  expect_setequal(unique(as.character(b1$label)), c("Treat All", "Full", "Age"))
  expect_equal(range(rib(suppressWarnings(plot(res, bands = TRUE, per = 1e5)))$upper),
               range(b1$upper) * 1e5)
  expect_error(plot(res, bands = "yes"), "`bands` must be TRUE or FALSE")
  skip_if_not_installed("plotly")
  b <- plotly::plotly_build(suppressWarnings(plot(res, bands = TRUE, interactive = TRUE)))
  fills <- vapply(b$x$data, function(t) identical(t$fill, "toself") || identical(t$fill, "tonexty"),
                  logical(1))
  expect_true(any(fills))
  full <- b$x$data[[which(vapply(b$x$data, function(t) identical(t$name, "Full"), TRUE))[1]]]
  expect_match(full$text[3], "^ \\(-?[0-9.]+ to -?[0-9.]+\\)$")
})

test_that("c() keeps the intervals only if all results have the same", {
  cs <- ci_case()
  v <- cs$v
  z <- seq(0.01, 0.3, by = 0.01)
  set.seed(1)
  a <- suppressWarnings(nb_val(v, "p1", "day30", thresholds = z, label = "A", interval = 0.95,
                               n_sim = 100))
  set.seed(2)
  b <- suppressWarnings(nb_val(v, "p2", "day30", thresholds = z, label = "B", interval = 0.95,
                               n_sim = 100))
  ab <- c(a, b)
  expect_equal(ab$interval$level, 0.95)
  expect_true("NB_model_lower" %in% names(ab$estimates))
  expect_warning(ab2 <- c(a, nb_val(v, "p2", "day30", thresholds = z, label = "B")),
                 "do not have the same credible intervals")
  expect_null(ab2$interval)
  expect_false("NB_model_lower" %in% names(ab2$estimates))
})

test_that("interval and n_sim are checked", {
  g <- gusto_split()
  v <- g$val[1:200, ]
  v$p <- 0.05 + 0.3 * (v$age - min(v$age)) / diff(range(v$age))
  msg <- "`interval` must be the level of the credible intervals"
  expect_error(nb_val(v, "p", "day30", interval = 95), msg)
  expect_error(nb_val(v, "p", "day30", interval = 0), msg)
  expect_error(nb_val(v, "p", "day30", interval = c(0.9, 0.95)), msg)
  expect_error(nb_val(v, "p", "day30", interval = 0.95, n_sim = 1), "at least 2")
})

test_that("every use of credible intervals warns and recommends value of information", {
  res <- ci_case()$res
  msg <- "not a good measure of the uncertainty that matters for choosing a strategy.*evpi\\(\\)"
  g <- gusto_split()
  v <- g$val[1:200, ]
  v$p <- 0.05 + 0.3 * (v$age - min(v$age)) / diff(range(v$age))
  expect_warning(nb_val(v, "p", "day30", thresholds = c(0.05, 0.1), interval = 0.95, n_sim = 20),
                 msg)
  expect_warning(capture.output(print(res)), msg)
  expect_warning(summary(res), msg)
  expect_warning(as.data.frame(res), msg)
  expect_warning(plot(res, bands = TRUE), msg)
  expect_warning(plot(res, type = "delta", bands = TRUE), msg)
  # no intervals shown (the default plot), no warning
  expect_no_warning(plot(res))
  r0 <- ci_case(interval = NULL)$res
  expect_no_warning(capture.output(print(r0)))
  expect_no_warning(summary(r0))
  expect_no_warning(as.data.frame(r0))
  expect_no_warning(plot(r0))
  expect_no_warning(nb_val(v, "p", "day30", thresholds = c(0.05, 0.1)))
})
