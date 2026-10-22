# Optional colours of the curves, and the messages when an optional plotting
# package is missing.

plot_case <- function() {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  z <- seq(0.01, 0.3, by = 0.01)
  res <- nb_val(g$val, list(Full = f1, Age = f2), "day30", thresholds = z)
  set.seed(1)
  e <- evpi(res, n_sim = 100)
  list(res = res, e = e)
}

# The colours a ggplot gives each curve (by label).
gg_colors <- function(p, panel = NULL) {
  b <- ggplot2::ggplot_build(if (is.null(panel)) p else p[[panel]])
  d <- b$data[[length(b$data)]]
  lev <- levels(b$plot$data$label)
  stats::setNames(d$colour[match(seq_along(lev), d$group)], lev)
}

test_that("colors: named by label, partial, or unnamed for the models", {
  cs <- plot_case()
  res <- cs$res
  default <- gg_colors(plot(res))
  # named: only those change
  p <- gg_colors(plot(res, colors = c(Full = "darkred", `Treat All` = "grey50")))
  expect_equal(unname(p[c("Full", "Treat All")]), c("darkred", "grey50"))
  expect_equal(p[c("Treat None", "Age")], default[c("Treat None", "Age")])
  # unnamed: the models in order
  p <- gg_colors(plot(res, colors = c("red", "blue")))
  expect_equal(unname(p[c("Full", "Age")]), c("red", "blue"))
  expect_equal(p[c("Treat All", "Treat None")], default[c("Treat All", "Treat None")])
  # the delta plot has the models only
  p <- gg_colors(plot(res, type = "delta", colors = c(Age = "#123456")))
  expect_equal(unname(p["Age"]), "#123456")
  # hex codes and the same colours for the same model in other plots
  pu <- gg_colors(plot(p_useful(cs$e), colors = c(Full = "darkred")))
  expect_equal(unname(pu["Full"]), "darkred")
  pb <- gg_colors(plot(p_best(cs$e), colors = c(Full = "darkred", `Treat None` = "black")))
  expect_equal(unname(pb[c("Full", "Treat None")]), c("darkred", "black"))
  pb2 <- gg_colors(plot(cs$e, type = "p_best", colors = c("red", "blue")))
  expect_equal(unname(pb2[c("Full", "Age")]), c("red", "blue"))
  ev <- gg_colors(plot(cs$e, colors = "purple"))
  expect_equal(unname(ev), "purple")
})

test_that("colors in the combined plot: the decision curve and the EVPI/EVSI panel", {
  skip_if_not_installed("patchwork")
  cs <- plot_case()
  set.seed(2)
  o <- evsi(cs$res, future_n = c(500, 2000), n_sim = 50)
  p <- plot(cs$res, evpi = cs$e, evsi = o,
            colors = c(Full = "darkred", EVPI = "orange", `EVSI (N = 500)` = "green"))
  top <- gg_colors(p, 1)
  bottom <- gg_colors(p, 2)
  expect_equal(unname(top["Full"]), "darkred")
  expect_equal(unname(bottom[c("EVPI", "EVSI (N = 500)")]), c("orange", "green"))
  # EVSI plot: one curve per threshold
  pe <- gg_colors(plot(o, thresholds = c(0.05, 0.1), colors = c("red", "blue")))
  expect_equal(unname(pe), c("red", "blue"))
})

test_that("colors are checked", {
  cs <- plot_case()
  res <- cs$res
  expect_error(plot(res, colors = "red"),
               "without names must have 2 colours, one per model in order \\(Full, Age\\)")
  expect_error(plot(res, colors = c(Nope = "red")),
               "names \"Nope\", which is not a curve of this plot")
  expect_error(plot(res, colors = c(Full = "notacolour")), "invalid colour: \"notacolour\"")
  expect_error(plot(res, colors = 1:2), "must be a character vector of colours")
  expect_error(plot(res, colors = c(Full = "red", "blue")), "Name every colour")
  expect_error(plot(p_best(cs$e), colors = c("red", "blue", "green")),
               "2 colours, one per model")
})

test_that("interactive plots use the colours", {
  skip_if_not_installed("plotly")
  cs <- plot_case()
  p <- plot(cs$res, colors = c(Full = "darkred"), interactive = TRUE)
  b <- plotly::plotly_build(p)
  cols <- vapply(b$x$data, function(tr) tr$line$color %||% NA_character_, character(1))
  nms <- vapply(b$x$data, function(tr) tr$name %||% NA_character_, character(1))
  expect_equal(unname(cols[nms == "Full"][1]), "darkred")
})

test_that("a missing plotly or patchwork gives a message on how to install it", {
  cs <- plot_case()
  local_mocked_bindings(has_package = function(pkg) FALSE)
  msg <- "need the plotly package. Install it with install.packages\\(\"plotly\"\\)"
  expect_error(plot(cs$res, interactive = TRUE), msg)
  expect_error(plot(cs$e, interactive = TRUE), msg)
  expect_error(plot(p_useful(cs$e), interactive = TRUE), msg)
  expect_error(plot(p_best(cs$e), interactive = TRUE), msg)
  set.seed(3)
  o <- evsi(cs$res, future_n = 500, n_sim = 20)
  expect_error(plot(o, interactive = TRUE), msg)
  # the combined (static) plot needs patchwork
  expect_error(plot(cs$res, evpi = cs$e),
               "needs the patchwork package. Install it with install.packages\\(\"patchwork\"\\)")
  # without them, single static plots still work
  expect_s3_class(plot(cs$res), "ggplot")
})
