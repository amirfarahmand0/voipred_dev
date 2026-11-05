# Net false positives avoided (net interventions avoided; Vickers and Elkin
# 2006, as in dcurves::net_intervention_avoided()): (NB_model - NB_all) / odds.

fp_case <- function() {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  v <- g$val[1:500, ]
  v$p1 <- unname(predict(f1, v, type = "response"))
  v$p2 <- unname(predict(f2, v, type = "response"))
  res <- nb_val(v, list(Full = "p1", Age = "p2"), "day30", thresholds = c(0, seq(0.01, 0.3, by = 0.01)))
  list(res = res, v = v)
}

test_that("the definition: (NB_model - NB_all) / (t / (1 - t)) = (TN - FN (1 - t) / t) / n", {
  cs <- fp_case()
  e <- cs$res$estimates
  ok <- e$threshold > 0
  odds <- e$threshold / (1 - e$threshold)
  expect_equal(e$net_FP_avoided[ok], ((e$NB_model - e$NB_all) / odds)[ok])
  expect_true(all(is.na(e$net_FP_avoided[!ok])))           # not defined at t = 0
  cn <- cs$res$counts
  cn <- cn[cn$threshold > 0, ]
  direct <- (cn$TN - cn$FN * (1 - cn$threshold) / cn$threshold) / nrow(cs$v)
  m <- merge(e[ok, c("threshold", "model", "net_FP_avoided")], data.frame(threshold = cn$threshold,
             model = cn$model, direct = direct), by = c("threshold", "model"))
  expect_equal(m$net_FP_avoided, m$direct)
  # 0.01 above treat all at 10% = 0.09 net false positives avoided
  expect_equal(net_fp_avoided(0.03, 0.02, 0.1), 0.09)
})

test_that("print(), summary(), as.data.frame() and per", {
  cs <- fp_case()
  res <- cs$res
  out <- capture.output(print(res, n = 1))
  expect_match(out[3], "delta_NB +net_FP_avoided$")
  s <- summary(res, thresholds = c(0.05, 0.1))
  expect_true("net_FP_avoided" %in% names(s$table))
  sp <- summary(res, thresholds = c(0.05, 0.1), per = 100)
  expect_equal(sp$table$net_FP_avoided, s$table$net_FP_avoided * 100)
  d <- as.data.frame(res, per = 100)
  expect_true(all(c("net_FP_avoided", "net_FP_avoided_per_100") %in% names(d)))
})

test_that("plot(type = \"fp_avoided\")", {
  cs <- fp_case()
  p <- plot(cs$res, type = "fp_avoided")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "Net false positives avoided (against treat all)")
  expect_equal(sort(unique(as.character(p$data$label))), c("Age", "Full"))
  expect_equal(plot(cs$res, type = "fp_avoided", per = 100)$labels$y,
               "Net false positives avoided (per 100 patients)")
  b <- ggplot2::ggplot_build(p)$plot$data
  e <- cs$res$estimates
  expect_equal(b$value[b$label == "Full"], e$net_FP_avoided[e$model == "Full" & e$threshold > 0])
})

test_that("credible intervals, subgroups and c()", {
  cs <- fp_case()
  set.seed(1)
  r <- suppressWarnings(nb_val(cs$v, "p1", "day30", thresholds = seq(0.01, 0.3, by = 0.01),
                               interval = 0.95, n_sim = 200))
  e <- r$estimates
  expect_true(all(e$net_FP_avoided_lower <= e$net_FP_avoided &
                    e$net_FP_avoided <= e$net_FP_avoided_upper))
  layers <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_true("GeomRibbon" %in% layers(suppressWarnings(plot(r, type = "fp_avoided", bands = TRUE))))
  sg <- nb_val(cs$v, "p1", "day30", thresholds = c(0.05, 0.1), subgroup = "sex")
  d <- as.data.frame(sg)
  expect_true("net_FP_avoided" %in% names(d))
  expect_equal(d$net_FP_avoided, net_fp_avoided(d$NB_model, d$NB_all, d$threshold))
  a <- nb_val(cs$v, "p1", "day30", thresholds = c(0.05, 0.1), label = "A")
  b <- nb_val(cs$v, "p2", "day30", thresholds = c(0.05, 0.1), label = "B")
  expect_true("net_FP_avoided" %in% names(c(a, b)$estimates))
})

test_that("EVPI and EVSI in net false positives avoided (unit = \"fp\")", {
  cs <- fp_case()
  res <- nb_val(cs$v, list(Full = "p1", Age = "p2"), "day30", thresholds = seq(0.01, 0.3, by = 0.01))
  set.seed(2)
  e <- evpi(res, n_sim = 300)
  odds <- e$estimates$threshold / (1 - e$estimates$threshold)
  s <- summary(e, unit = "fp_avoided")
  expect_equal(s$table$EVPI, e$estimates$EVPI / odds)
  expect_match(s$header, ", in net false positives avoided$")
  expect_equal(summary(e, unit = "fp_avoided", per = 1000)$table$EVPI, e$estimates$EVPI / odds * 1000)
  d <- as.data.frame(e, unit = "fp_avoided")
  expect_equal(d$EVPI_fp_avoided, e$estimates$EVPI / odds)
  expect_equal(d$EVPI, e$estimates$EVPI)                    # per-patient column kept
  expect_false("EVPI_fp_avoided" %in% names(as.data.frame(e)))
  expect_equal(plot(e, unit = "fp_avoided")$labels$y, "EVPI (net false positives avoided)")
  expect_output(print(e, unit = "fp_avoided"), "in net false positives avoided")
  set.seed(3)
  o <- evsi(res, future_n = c(500, 2000), n_sim = 100)
  so <- summary(o, thresholds = 0.05, unit = "fp_avoided")
  expect_equal(so$table$EVSI, summary(o, thresholds = 0.05)$table$EVSI / (0.05 / 0.95))
  expect_true(all(c("EVSI_fp_avoided", "EVPI_fp_avoided") %in% names(as.data.frame(o, unit = "fp_avoided"))))
  expect_equal(plot(o, unit = "fp_avoided")$labels$y, "EVSI (net false positives avoided)")
  expect_error(summary(e, unit = "x"), "`unit` must be \"tp\"")
})
