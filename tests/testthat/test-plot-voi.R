# Plots of the value-of-information measures, the decision curve with an EVPI /
# EVSI panel, and the interactive (plotly) versions of all plots.

voi_plot_setup <- function() {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  val <- g$val[1:1500, ]
  z <- c(0.01, 0.02, 0.05, 0.1, 0.2)
  r1 <- nb(f1, val, "day30", thresholds = z, label = "Full")
  r2 <- nb(f2, val, "day30", thresholds = z, label = "Age")
  set.seed(1)
  d <- draws_t(r1, r2, n_sim = 100)
  set.seed(2)
  o <- evsi(r1, future_n = c(200, 1000), n_sim = 60)
  list(r1 = r1, r2 = r2, d = d, e = evpi_from_draws(d), e1 = evpi_from_draws(draws_t(r1, n_sim = 100)),
       pu = p_useful_from_draws(d), pb = p_best_from_draws(d), o = o)
}

line_data <- function(p) {
  ld <- ggplot2::layer_data(p, length(p$layers))
  split(ld$y, levels(p$data$label)[ld$group])
}

test_that("plot_evpi(): EVPI or relative EVPI against the threshold", {
  s <- voi_plot_setup()
  p <- plot_evpi(s$e)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$label), "EVPI")
  expect_equal(line_data(p)[["EVPI"]], s$e$EVPI)
  expect_equal(p$labels$y, "EVPI")
  pr <- plot_evpi(s$e, type = "relative")
  expect_equal(ggplot2::layer_data(pr, 1)$yintercept, 1)
  expect_equal(line_data(pr)[["EVPI"]], s$e$rEVPI[is.finite(s$e$rEVPI)])
  # several results, labelled by argument names or `labels`
  p2 <- plot_evpi(s$e, single = s$e1)
  expect_equal(levels(p2$data$label), c("EVPI 1", "single"))
  expect_equal(levels(plot_evpi(s$e, s$e1, labels = c("a", "b"))$data$label), c("a", "b"))
  expect_error(plot_evpi(s$pu), "needs results of evpi")
  expect_error(plot_evpi(s$e, tpye = "relative"), "Unknown argument `tpye`")
})

test_that("plot_evsi(): against the future sample size", {
  s <- voi_plot_setup()
  p <- plot_evsi(s$o)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Future sample size")
  lv <- levels(p$data$label)
  expect_equal(lv, paste0("Threshold ", c(1, 2, 5, 10, 20), "%"))
  first <- p$data[p$data$label == lv[2], ]
  expect_equal(first$x, c(0, 200, 1000))
  expect_equal(first$value, c(0, s$o$EVSI[s$o$threshold == 0.02]))
  ref <- ggplot2::layer_data(p, 1)                      # the dashed EVPI lines
  expect_true(all(ref$linetype == "dashed"))
  expect_error(plot_evsi(s$o, type = "threshold"), "unused argument")   # removed
  expect_error(plot_evsi(s$e), "needs a result of evsi")
})

test_that("plot_p_useful() and plot_p_best(): probabilities, in the colours of plot_dca()", {
  s <- voi_plot_setup()
  p <- plot_p_useful(s$pu)
  expect_equal(levels(p$data$label), c("Full", "Age"))
  expect_equal(line_data(p)[["Age"]], s$pu$p_useful[s$pu$model == "Age"])
  expect_equal(p$coordinates$limits$y, c(0, 1))
  b <- plot_p_best(s$pb)
  expect_equal(levels(b$data$label), c("Treat All", "Treat None", "Full", "Age"))
  expect_equal(line_data(b)[["Treat None"]], s$pb$p_best[s$pb$strategy == "Treat None"])
  col <- function(p, l) unique(ggplot2::layer_data(p, length(p$layers))$colour[
    levels(p$data$label)[ggplot2::layer_data(p, length(p$layers))$group] == l])
  dca <- plot_dca(s$r1, s$r2)
  for (l in c("Full", "Age")) {
    expect_equal(col(p, l), col(dca, l))
    expect_equal(col(b, l), col(dca, l))
  }
  expect_error(plot_p_useful(s$pb), "needs a result of p_useful")
  expect_error(plot_p_best(s$pu), "needs a result of p_best")
})

test_that("plot_dca() with an EVPI / EVSI panel below", {
  skip_if_not_installed("patchwork")
  s <- voi_plot_setup()
  p <- plot_dca(s$r1, evpi = s$e1, evsi = s$o)
  expect_s3_class(p, "patchwork")
  bottom <- p[[2]]
  expect_equal(levels(bottom$data$label), c("EVPI", "EVSI (N = 200)", "EVSI (N = 1,000)"))
  expect_equal(line_data(bottom)[["EVPI"]], s$e1$EVPI)
  # the two panels share the threshold axis
  expect_equal(p[[1]]$coordinates$limits$x, bottom$coordinates$limits$x)
  # EVSI alone: EVPI is taken from evsi()
  q <- plot_dca(s$r1, evsi = s$o)
  expect_equal(line_data(q[[2]])[["EVPI"]], s$o$EVPI[s$o$future_n == 200])
  withr::local_pdf(NULL)
  expect_no_error(print(p))
  expect_error(plot_dca(s$r1, evpi = s$o), "`evpi` must be a result of evpi")
  expect_error(plot_dca(s$r1, evsi = s$e1), "`evsi` must be a result of evsi")
})

test_that("interactive = TRUE gives plotly widgets with the same curves", {
  skip_if_not_installed("plotly")
  s <- voi_plot_setup()
  traces <- function(w) {
    b <- plotly::plotly_build(w)$x$data
    nm <- vapply(b, function(t) if (is.null(t$name)) "" else t$name, character(1))
    stats::setNames(lapply(b, function(t) as.numeric(t$y)), nm)
  }
  w <- plot_dca(s$r1, s$r2, interactive = TRUE)
  expect_s3_class(w, "plotly")
  tr <- traces(w)
  expect_setequal(names(tr), c("Treat All", "Treat None", "Full", "Age"))
  expect_equal(tr[["Full"]], s$r1$NB_model)
  td <- traces(plot_dca(s$r1, s$r2, type = "delta", interactive = TRUE))
  expect_equal(td[["Age"]], s$r2$delta_NB)
  # decision curve with the VoI panel: one shared x axis, two y axes
  wv <- plot_dca(s$r1, evpi = s$e1, evsi = s$o, interactive = TRUE)
  b <- plotly::plotly_build(wv)$x
  expect_true(all(vapply(b$data, function(t) is.null(t$xaxis) || t$xaxis == "x", logical(1))))
  expect_setequal(unique(unlist(lapply(b$data, `[[`, "yaxis"))), c("y", "y2"))
  tv <- traces(wv)
  expect_equal(tv[["EVPI"]], s$e1$EVPI)
  expect_equal(tv[["EVSI (N = 1,000)"]], s$o$EVSI[s$o$future_n == 1000])
  expect_equal(traces(plot_evpi(s$e, interactive = TRUE))[["EVPI"]], s$e$EVPI)
  expect_equal(traces(plot_evsi(s$o, interactive = TRUE))[["Threshold 2%"]],
               c(0, s$o$EVSI[s$o$threshold == 0.02]))
  expect_equal(traces(plot_p_useful(s$pu, interactive = TRUE))[["Full"]],
               s$pu$p_useful[s$pu$model == "Full"])
  expect_equal(traces(plot_p_best(s$pb, interactive = TRUE))[["Treat All"]],
               s$pb$p_best[s$pb$strategy == "Treat All"])
  expect_error(plot_evpi(s$e, interactive = "yes"), "TRUE or FALSE")
})

test_that("static and interactive plots share one look, with the legend on the right", {
  s <- voi_plot_setup()
  p <- plot_dca(s$r1, s$r2)
  expect_equal(p$theme$legend.position, "right")
  expect_null(p$labels$colour)                                   # no legend title
  expect_s3_class(p$theme$panel.grid.minor, "element_blank")
  expect_equal(plain_label(c(0, 6e-4, 0.01, NA, 2000)), c("0", "0.0006", "0.01", NA, "2000"))
  skip_if_not_installed("plotly")
  w <- plotly::plotly_build(plot_dca(s$r1, s$r2, interactive = TRUE))$x$layout
  expect_equal(w$legend$x, 1.02)
  expect_equal(w$legend$yanchor, "middle")
  expect_equal(w$yaxis$tickformat, "~g")
  skip_if_not_installed("patchwork")
  st <- plot_dca(s$r1, evpi = s$e1, evsi = s$o)
  expect_s3_class(st[[1]]$theme$axis.text.x, "element_blank")    # shared x axis
  expect_false(identical(st$patches$layout$guides, "collect"))   # a legend per panel
  skip_if_not(utils::packageVersion("plotly") >= "4.11.0")
  b <- plotly::plotly_build(plot_dca(s$r1, evpi = s$e1, evsi = s$o, interactive = TRUE))$x
  shown <- Filter(function(t) !isFALSE(t$showlegend), b$data)
  leg <- vapply(shown, function(t) if (is.null(t$legend)) "legend" else t$legend, "")
  nm <- vapply(shown, `[[`, "", "name")
  expect_setequal(nm[leg == "legend"], c("Treat All", "Treat None", "Full"))
  expect_setequal(nm[leg == "legend2"], c("EVPI", "EVSI (N = 200)", "EVSI (N = 1,000)"))
  expect_gt(b$layout$legend$y, b$layout$legend2$y)              # each beside its panel
})

test_that("the interactive delta-NB axis title uses codes plotly can show", {
  skip_if_not_installed("plotly")
  s <- voi_plot_setup()
  w <- plotly::plotly_build(plot_dca(s$r1, type = "delta", interactive = TRUE))$x$layout
  ttl <- w$yaxis$title
  if (is.list(ttl)) ttl <- ttl$text
  # plotly decodes numeric codes (&#916;) but not most named ones (&Delta;)
  expect_false(grepl("&[A-Za-z]+;", ttl))
  expect_match(ttl, "&#916;NB", fixed = TRUE)
  expect_match(ttl, "&#8722;", fixed = TRUE)
})

test_that("the value-of-information axis is named after what it shows", {
  s <- voi_plot_setup()
  expect_equal(plot_evpi(s$e)$labels$y, "EVPI")
  expect_equal(plot_evsi(s$o)$labels$y, "EVSI")
  skip_if_not_installed("patchwork")
  expect_equal(plot_dca(s$r1, evpi = s$e1)[[2]]$labels$y, "EVPI")
  expect_equal(plot_dca(s$r1, evpi = s$e1, evsi = s$o)[[2]]$labels$y, "EVPI / EVSI")
  expect_equal(plot_dca(s$r1, evsi = s$o)[[2]]$labels$y, "EVPI / EVSI")
  skip_if_not_installed("plotly")
  w <- plotly::plotly_build(plot_dca(s$r1, evpi = s$e1, interactive = TRUE))$x$layout
  t2 <- w$yaxis2$title
  if (is.list(t2)) t2 <- t2$text
  expect_equal(t2, "EVPI")
})
