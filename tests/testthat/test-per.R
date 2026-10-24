# `per`: the net benefit, EVPI and EVSI per `per` patients in print(), summary(),
# plot() and as.data.frame(); the results themselves are always per patient.

per_case <- function() {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  v <- g$val[1:500, ]
  v$p1 <- unname(predict(f1, v, type = "response"))
  v$p2 <- unname(predict(f2, v, type = "response"))
  res <- nb_val(v, list(Full = "p1", Age = "p2"), "day30", thresholds = seq(0.01, 0.3, by = 0.01))
  set.seed(1)
  e <- evpi(res, n_sim = 200)
  set.seed(2)
  o <- evsi(res, future_n = c(500, 2000), n_sim = 100)
  list(res = res, e = e, o = o)
}

test_that("as.data.frame(): per-patient columns kept, scaled ones added", {
  cs <- per_case()
  K <- 1e5
  d <- as.data.frame(cs$res, per = K)
  expect_equal(d[, 1:6], as.data.frame(cs$res))
  expect_named(d, c("threshold", "model", "NB_model", "NB_all", "NB_none", "delta_NB",
                    "NB_model_per_100000", "NB_all_per_100000", "NB_none_per_100000",
                    "delta_NB_per_100000"))
  expect_equal(d$delta_NB_per_100000, cs$res$estimates$delta_NB * K)
  de <- as.data.frame(cs$e, per = K)
  expect_equal(de$EVPI_per_100000, cs$e$estimates$EVPI * K)
  expect_equal(de$ENB_current_per_100000, cs$e$estimates$ENB_current * K)
  expect_false("rEVPI_per_100000" %in% names(de))              # no unit
  do <- as.data.frame(cs$o, per = K)
  expect_equal(do$EVSI_per_100000, cs$o$estimates$EVSI * K)
  expect_named(as.data.frame(cs$o, per = 1000), c("threshold", "future_n", "EVSI", "EVPI",
                                                  "EVSI_per_1000", "EVPI_per_1000"))
  # per = 1 (default): unchanged; the results themselves are never scaled
  expect_identical(as.data.frame(cs$e), cs$e$estimates)
  expect_identical(as.data.frame(cs$e, per = 1), cs$e$estimates)
  # probabilities have no per
  expect_identical(as.data.frame(cs$e$p_useful), cs$e$p_useful$estimates)
})

test_that("print() and summary(): scaled values and the scale in the header", {
  cs <- per_case()
  K <- 1e5
  out <- capture.output(print(cs$res, n = 1, per = K))
  expect_match(out[1], ", per 100,000 patients$")
  row1 <- strsplit(trimws(out[4]), " +")[[1]]
  e1 <- cs$res$estimates[1, ]
  expect_equal(as.numeric(row1[3]), signif(e1$NB_model * K, 4), tolerance = 1e-3)
  s <- summary(cs$res, per = K)
  expect_equal(s$table$NB_model, (summary(cs$res)$table$NB_model) * K)
  expect_equal(s$table$Se, summary(cs$res)$table$Se)          # no unit
  expect_match(s$header, ", per 100,000 patients$")
  expect_match(s$largest, "^Largest delta_NB: Full [0-9,.]+ at ")
  # EVPI
  se <- summary(cs$e, per = K)
  expect_equal(se$table$EVPI, cs$e$estimates$EVPI * K)
  expect_equal(se$table$rEVPI, cs$e$estimates$rEVPI)
  expect_match(se$header, ", per 100,000 patients$")
  expect_output(print(cs$e, per = K), "per 100,000 patients")
  # EVSI: the scale on the first header line
  out_o <- capture.output(print(cs$o, per = 1000))
  expect_match(out_o[1], ", per 1,000 patients$")
  so <- summary(cs$o, per = K)
  expect_equal(so$table$EVSI, summary(cs$o)$table$EVSI * K)
  expect_equal(so$table[["EVSI / EVPI"]], summary(cs$o)$table[["EVSI / EVPI"]])
  # default: no scale in the header
  expect_false(grepl(", per [0-9]", capture.output(print(cs$e))[1]))
})

test_that("the option voipred.per sets the default of per", {
  cs <- per_case()
  withr::local_options(voipred.per = 1000)
  expect_match(capture.output(print(cs$e))[1], ", per 1,000 patients$")
  expect_true("EVPI_per_1000" %in% names(as.data.frame(cs$e)))
  expect_equal(summary(cs$res)$table$NB_all[1], summary(cs$res, per = 1)$table$NB_all[1] * 1000)
  # a value given in the call wins
  expect_false(grepl(", per [0-9]", capture.output(print(cs$e, per = 1))[1]))
})

test_that("plot(): scaled values, and the scale in the y-axis labels", {
  cs <- per_case()
  K <- 1e5
  ylab <- function(p) p$labels$y
  p <- plot(cs$res, per = K)
  expect_equal(ylab(p), "Net Benefit (per 100,000 patients)")
  b <- ggplot2::ggplot_build(p)$plot$data
  b1 <- ggplot2::ggplot_build(plot(cs$res))$plot$data
  expect_equal(b$value, b1$value * K)
  expect_equal(ylab(plot(cs$res, type = "delta", per = K)), "\u0394NB (per 100,000 patients)")
  expect_equal(ylab(plot(cs$e, per = K)), "EVPI (per 100,000 patients)")
  expect_equal(ylab(plot(cs$o, per = K)), "EVSI (per 100,000 patients)")
  # no unit: unchanged
  expect_equal(ylab(plot(cs$e, type = "relative", per = K)), "Relative EVPI")
  expect_equal(ylab(plot(cs$e, type = "p_best", per = K)), "P(best)")
  # default: the usual labels
  expect_equal(ylab(plot(cs$res)), "Net Benefit")
  expect_equal(ylab(plot(cs$e)), "EVPI")
  skip_if_not_installed("patchwork")
  pp <- plot(cs$res, evpi = cs$e, evsi = cs$o, per = K)
  expect_equal(ylab(pp[[1]]), "Net Benefit (per 100,000 patients)")
  expect_equal(ylab(pp[[2]]), "EVPI / EVSI (per 100,000 patients)")
  bb <- ggplot2::ggplot_build(pp[[2]])$plot$data
  expect_equal(max(bb$value[bb$label == "EVPI"]), max(cs$e$estimates$EVPI) * K)
})

test_that("interactive plots show the scale on the axis", {
  skip_if_not_installed("plotly")
  cs <- per_case()
  p <- plotly::plotly_build(plot(cs$e, per = 1e5, interactive = TRUE))
  title <- p$x$layout$yaxis$title
  if (is.list(title)) title <- title$text
  expect_equal(title, "EVPI (per 100,000 patients)")
})

test_that("per is checked", {
  cs <- per_case()
  msg <- "`per` must be a single positive number"
  expect_error(print(cs$res, per = 0), msg)
  expect_error(summary(cs$e, per = -1), msg)
  expect_error(as.data.frame(cs$o, per = c(10, 100)), msg)
  expect_error(plot(cs$res, per = "a"), msg)
  expect_error(plot(cs$e, per = NA), msg)
})
