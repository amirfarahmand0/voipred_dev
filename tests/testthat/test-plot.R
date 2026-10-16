# delta_NB, labels, and plot() / plot() of nb() results.

toy <- data.frame(y = c(1, 0, 1, 0, 0), p = c(0.9, 0.2, 0.4, 0.05, 0.6))
z <- c(0, 0.1, 0.3, 0.5, 0.7)

curves <- function(p) {
  d <- ggplot2::layer_data(p, i = length(p$layers))   # the geom_line layer
  lab <- levels(p$data$label)
  split(d$y, lab[d$group])
}

test_that("delta_NB is NB_model - max(NB_all, NB_none)", {
  res <- nb(toy, fixed_model(toy$p), "y", thresholds = z)
  expect_equal(res$estimates$delta_NB, res$estimates$NB_model - pmax(res$estimates$NB_all, res$estimates$NB_none))
  # survival and competing risks, with NA where NB_model is NA
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  r <- quiet_nb(s$val, cx, surv_outcome, horizon = 5, thresholds = surv_z)
  expect_equal(r$estimates$delta_NB, r$estimates$NB_model - pmax(r$estimates$NB_all, 0))
  expect_identical(is.na(r$estimates$delta_NB), is.na(r$estimates$NB_model))
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  r2 <- quiet_nb(s$val, csc, cr_outcome, horizon = 5, thresholds = surv_z, cause = "death")
  expect_equal(r2$estimates$delta_NB, r2$estimates$NB_model - pmax(r2$estimates$NB_all, 0))
})

test_that("label: default from the model's name, given, or reserved", {
  fit <- fixed_model(toy$p)
  expect_equal(nb(toy, fit, "y", thresholds = z)$label, "fit")
  expect_equal(nb(toy, fixed_model(toy$p), "y", thresholds = z)$label, "Model")
  expect_equal(nb(toy, fit, "y", thresholds = z, label = "A")$label, "A")
  expect_output(print(nb(toy, fit, "y", thresholds = z, label = "A")),
                "Net benefit of A.*delta_NB")
  expect_error(nb(toy, fit, "y", label = "treat all"), "default strategies")
  expect_error(nb(toy, fit, "y", label = c("a", "b")), "single non-empty string")
  expect_error(nb(toy, fit, "y", label = ""), "single non-empty string")
})

test_that("decision curve: model, Treat All and Treat None", {
  res <- nb(toy, fixed_model(toy$p), "y", thresholds = z, label = "M")
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$label), c("Treat All", "Treat None", "M"))
  cv <- curves(p)
  expect_equal(cv[["M"]], res$estimates$NB_model)
  expect_equal(cv[["Treat All"]], res$estimates$NB_all)
  expect_equal(cv[["Treat None"]], res$estimates$NB_none)
  expect_equal(p$labels$x, "Threshold Probability")
  expect_equal(p$labels$y, "Net Benefit")
  # y-axis from -0.1 * prevalence to the prevalence, as dcurves
  expect_equal(p$coordinates$limits$y, c(-0.04, 0.4))
})

test_that("delta plot: delta_NB of each model, line at 0", {
  res <- nb(toy, fixed_model(toy$p), "y", thresholds = z, label = "M")
  p <- plot(res, type = "delta")
  expect_equal(levels(p$data$label), "M")
  expect_equal(curves(p)[["M"]], res$estimates$delta_NB)
  expect_equal(ggplot2::layer_data(p, 1)$yintercept, 0)
})

test_that("several results on one plot, for both types", {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  zz <- seq(0, 0.5, 0.01)
  r1 <- nb(g$val, f1, "day30", thresholds = zz, label = "Full")
  r2 <- nb(g$val, f2, "day30", thresholds = zz, label = "Age")
  p <- plot(r1, r2)
  expect_equal(levels(p$data$label), c("Treat All", "Treat None", "Full", "Age"))
  expect_equal(curves(p)[["Age"]], r2$estimates$NB_model)
  d <- plot(r1, r2, type = "delta")
  expect_equal(levels(d$data$label), c("Full", "Age"))
  expect_equal(curves(d)[["Full"]], r1$estimates$delta_NB)
  expect_equal(curves(d)[["Age"]], r2$estimates$delta_NB)
  # each model has the same colour in both plots (ggplot2's default hues)
  col_of <- function(p) tapply(ggplot2::layer_data(p, length(p$layers))$colour,
                               levels(p$data$label)[ggplot2::layer_data(p, length(p$layers))$group],
                               unique)
  expect_equal(col_of(p)[c("Full", "Age")], col_of(d)[c("Full", "Age")])
  expect_equal(as.vector(col_of(p)[levels(p$data$label)]),
               c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))
  # names in ... and `labels`
  expect_equal(levels(plot(r1, Simple = r2)$data$label)[3:4], c("Full", "Simple"))
  expect_equal(levels(plot(r1, r2, type = "delta", labels = c("a", "b"))$data$label),
               c("a", "b"))
  # printing draws it
  withr::local_pdf(NULL)
  expect_no_error(print(d))
  # survival
  s <- rott_split()
  c1 <- survival::coxph(surv_formula, data = s$dev)
  c2 <- survival::coxph(survival::Surv(dtimey, death) ~ age + size, data = s$dev)
  s1 <- quiet_nb(s$val, c1, surv_outcome, horizon = 5, thresholds = surv_z, label = "Cox")
  s2 <- quiet_nb(s$val, c2, surv_outcome, horizon = 5, thresholds = surv_z, label = "Small")
  ps <- plot(s1, s2, type = "delta")
  expect_equal(curves(ps)[["Cox"]], s1$estimates$delta_NB[!is.na(s1$estimates$delta_NB)])
  expect_equal(plot(s1)$coordinates$limits$y, c(-0.1, 1) * s1$outcome$prevalence)
})

test_that("results that cannot be drawn together are refused", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  zz <- seq(0, 0.5, 0.01)
  r1 <- nb(g$val, f, "day30", thresholds = zz, label = "A")
  expect_error(plot(r1, nb(g$dev, f, "day30", thresholds = zz, label = "B")),
               "different data")
  expect_error(plot(r1, nb(g$val, f, "day30", thresholds = seq(0, 0.4, 0.01), label = "B")),
               "different thresholds")
  expect_error(plot(r1, nb(g$val[1:4000, ], f, "day30", thresholds = zz, label = "B")),
               "different number of patients")
  s <- rott_split()
  cx <- survival::coxph(surv_formula, data = s$dev)
  a <- quiet_nb(s$val, cx, surv_outcome, horizon = 5, thresholds = surv_z, label = "A")
  b <- quiet_nb(s$val, cx, surv_outcome, horizon = 3, thresholds = surv_z, label = "B")
  expect_error(plot(a, b), "\"A\" and \"B\" cannot be used together: different horizon")
  expect_error(plot(r1, b), "different outcome type")
  # same label twice, unknown argument, other objects
  expect_error(plot(r1, r1), "same label \\(\"A\"\\)")
  expect_error(plot(r1, ty = "delta"), "Unknown argument `ty`")
  expect_error(plot(r1, data.frame(x = 1)), "Only results of nb\\(\\)")
  expect_error(plot(r1, labels = c("a", "b")), "`labels` must be one string")
  expect_error(plot(r1, labels = "Treat None"), "cannot be used as labels")
  expect_error(plot(r1, type = "roc"), "should be one of")
})
