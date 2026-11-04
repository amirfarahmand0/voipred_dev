# Subgroup utility heterogeneity (Wang, Wang, Van Calster, van Geloven,
# Steyerberg and Wynants): nb_val(subgroup = ), the checks of summary(), the
# K-value, the plots and k_value().

# The motivating example of the paper (Tables 1 and 2): two equal subgroups,
# prevalence 50%, threshold 0.0625; counts per 100 patients.
motivating <- function() {
  mk <- function(tp, fp, fn, tn, g)
    data.frame(y = c(rep(1, tp), rep(0, fp), rep(1, fn), rep(0, tn)),
               p = c(rep(0.9, tp + fp), rep(0.01, fn + tn)), group = g)
  rbind(mk(45, 1, 5, 49, "Group 1"), mk(49, 1, 1, 49, "Group 2"))
}

test_that("the motivating example of the paper is reproduced", {
  d <- motivating()
  res <- nb_val(d, "p", "y", thresholds = 0.0625, label = "Model", subgroup = "group")
  # population NB: 0.4693 (model) and 0.4667 (treat all)
  expect_equal(round(res$estimates$NB_model, 4), 0.4693)
  expect_equal(round(res$estimates$NB_all, 4), 0.4667)
  s <- summary(res, k = c(`Group 1` = 60, `Group 2` = 45))
  t <- s$subgroups$table
  expect_equal(t[["diff_NB_Group 1"]], 0.45 - 0.01 / 15 - (0.5 - 0.5 / 15))
  expect_equal(t[["diff_NB_Group 2"]], 0.49 - 0.01 / 15 - (0.5 - 0.5 / 15))
  expect_equal(t$diff_NB_population, res$estimates$NB_model - res$estimates$NB_all)
  expect_equal(t$signs, "opposite")
  expect_equal(round(t$K, 3), 1.308)
  # population utility favours treat all (183 vs 182.98, summed over the two groups)
  expect_equal(t$diff_utility, (182.98 - 183) / 2, tolerance = 1e-10)
  expect_false(t$consistent)
  expect_match(paste(s$subgroups$lines, collapse = " "),
               "favour different strategies at 0.0625")
  # with equal k, no inconsistency is possible
  s_eq <- summary(res, k = c(`Group 1` = 50, `Group 2` = 50))
  expect_true(s_eq$subgroups$table$consistent)
})

test_that("SCORE2 (the case study): K = 2.156, consistent at 1.867, not at 2.333", {
  kv <- k_value(c(Men = -0.0016, Women = 0.00345))
  expect_equal(round(kv$K, 3), 2.156)
  expect_true(k_value(c(Men = -0.0016, Women = 0.00345), k = c(Men = 0.28, Women = 0.15))$consistent)
  expect_false(k_value(c(Men = -0.0016, Women = 0.00345), k = c(Men = 0.28, Women = 0.12))$consistent)
  expect_output(print(kv), "K-value \\(k_Men / k_Women\\): 2.156")
  expect_true(is.na(k_value(c(A = 0.01, B = 0.02))$K))
  expect_output(print(k_value(c(A = 0.01, B = 0.02))), "same sign")
  # unequal shares
  expect_equal(k_value(c(A = -0.001, B = 0.002), shares = c(A = 3, B = 1))$K, (0.25 * 0.002) / (0.75 * 0.001))
  expect_error(k_value(c(-0.001, 0.002)), "named by subgroup")
  expect_error(k_value(c(A = -0.001, B = 0.002), k = c(A = 1)), "two positive numbers")
})

test_that("the K-value is where the utility-weighted sum changes sign", {
  d <- motivating()
  res <- nb_val(d, "p", "y", thresholds = 0.0625, label = "Model", subgroup = "group")
  K <- summary(res)$subgroups$table$K
  just_below <- summary(res, k = c(`Group 1` = K * 0.999, `Group 2` = 1))$subgroups$table
  just_above <- summary(res, k = c(`Group 1` = K * 1.001, `Group 2` = 1))$subgroups$table
  expect_true(just_below$consistent)
  expect_false(just_above$consistent)
})

test_that("same signs: never inconsistent, whatever k", {
  g <- gusto_split()
  v <- g$val
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  v$p <- unname(predict(f, v, type = "response"))
  res <- nb_val(v, "p", "day30", thresholds = c(0.05, 0.1), label = "Full", subgroup = "sex")
  t <- summary(res)$subgroups$table
  same <- t$signs == "same"
  for (kk in list(c(male = 1, female = 10), c(male = 10, female = 1))) {
    tk <- summary(res, k = kk)$subgroups$table
    expect_true(all(tk$consistent[same]))
  }
  # the subgroup net benefit equals nb_val() on each subgroup
  m <- nb_val(v[v$sex == "male", ], "p", "day30", thresholds = c(0.05, 0.1), label = "Full")
  e <- res$subgroups$estimates
  expect_equal(e$NB_model[e$subgroup == "male"], m$estimates$NB_model)
  expect_equal(e$NB_all[e$subgroup == "male"], m$estimates$NB_all)
  # population delta NB = share-weighted subgroup delta NB (binary: exact)
  expect_equal(t$diff_NB_population, res$estimates$NB_model - res$estimates$NB_all)
})

test_that("more than two subgroups: k as values or ranges, no K-value", {
  set.seed(1)
  n <- 600
  d <- data.frame(y = rbinom(n, 1, 0.3), g = factor(rep(c("a", "b", "c"), each = 200)))
  d$p <- plogis(qlogis(0.3) + 1.2 * (d$y - 0.3) + rnorm(n))
  res <- nb_val(d, "p", "y", thresholds = c(0.2, 0.3, 0.4), subgroup = "g")
  t <- summary(res)$subgroups$table
  expect_null(t$K)
  tv <- summary(res, k = c(a = 1, b = 2, c = 3))$subgroups$table
  ex <- as.matrix(t[, c("diff_NB_a", "diff_NB_b", "diff_NB_c")]) %*% (c(a = 1, b = 2, c = 3) / 3)
  expect_equal(tv$diff_utility, as.vector(ex))
  tr <- summary(res, k = list(a = c(1, 1), b = c(2, 2), c = c(3, 3)))$subgroups$table
  expect_equal(tr$diff_utility_lower, as.vector(ex))
  expect_equal(tr$diff_utility_upper, as.vector(ex))
  wide <- summary(res, k = list(a = c(0.1, 10), b = c(0.1, 10), c = c(0.1, 10)))$subgroups$table
  expect_true(all(wide$diff_utility_lower <= wide$diff_utility_upper))
  expect_error(plot(res, type = "region", threshold = 0.3), "for two subgroups")
})

test_that("comparator and shares", {
  d <- motivating()
  res <- nb_val(d, list(A = "p", B = d$p * 0.9), "y", thresholds = 0.0625, subgroup = "group")
  t0 <- summary(res, comparator = "Treat None")$subgroups$table
  e <- res$subgroups$estimates
  expect_equal(t0[["diff_NB_Group 1"]][t0$model == "A"], e$NB_model[e$model == "A" & e$subgroup == "Group 1"])
  tb <- summary(res, comparator = "B")$subgroups$table
  expect_equal(unique(tb$model), "A")
  td <- summary(res, comparator = "best default")$subgroups$table
  expect_equal(td[["diff_NB_Group 1"]], summary(res)$subgroups$table[["diff_NB_Group 1"]])
  ts <- summary(res, shares = c(`Group 1` = 3, `Group 2` = 1))$subgroups$table
  expect_equal(ts$diff_NB_population,
               0.75 * ts[["diff_NB_Group 1"]] + 0.25 * ts[["diff_NB_Group 2"]])
  expect_error(summary(res, comparator = "Nope"), "`comparator` must be one of")
  expect_error(summary(res, shares = c(x = 1, y = 1)), "named by the subgroups")
  expect_error(summary(res, k = c(`Group 1` = -1, `Group 2` = 1)), "positive values")
  r0 <- nb_val(d, "p", "y", thresholds = 0.0625)
  expect_error(summary(r0, k = c(a = 1)), "are for subgroups")
})

test_that("subgroup can be a column, columns, a factor or a logical vector", {
  d <- motivating()
  d$sex <- rep(c("f", "m"), 100)
  a <- nb_val(d, "p", "y", thresholds = 0.0625, subgroup = "group")
  b <- nb_val(d, "p", "y", thresholds = 0.0625, subgroup = d$group)
  expect_equal(a$subgroups$estimates, b$subgroups$estimates)
  expect_equal(b$subgroups$variable, "group")
  f <- nb_val(d, "p", "y", thresholds = 0.0625,
              subgroup = factor(d$group, levels = c("Group 2", "Group 1")))
  expect_equal(levels(f$subgroups$group), c("Group 2", "Group 1"))
  l <- nb_val(d, "p", "y", thresholds = 0.0625, subgroup = d$group == "Group 1")
  expect_equal(levels(l$subgroups$group), c("FALSE", "TRUE"))
  cr <- nb_val(d, "p", "y", thresholds = 0.0625, subgroup = c("group", "sex"))
  expect_equal(nlevels(cr$subgroups$group), 4)
  expect_match(levels(cr$subgroups$group)[1], "Group 1, f")
  expect_error(nb_val(d, "p", "y", subgroup = d$group[-1]), "values but `data` has 200 rows")
  g2 <- d$group
  g2[3] <- NA
  expect_error(nb_val(d, "p", "y", subgroup = g2), "missing for 1 patients")
  expect_error(nb_val(d, "p", "y", subgroup = seq_len(200) / 200), "numeric with .* different values")
  expect_error(nb_val(d, "p", "y", subgroup = rep("a", 200)), "only one subgroup")
  expect_error(nb_val(d, "p", "y", subgroup = c("group", "nope")), "not in `data`")
  expect_error(nb_val(d, "p", "y", subgroup = "group", interval = 0.95),
               "not available with `subgroup`")
})

test_that("print(), as.data.frame(), plots, c() and survival outcomes", {
  d <- motivating()
  res <- nb_val(d, "p", "y", thresholds = c(0.0625, 0.1), label = "Model", subgroup = "group")
  expect_output(print(res), "# Subgroups \\(group\\): Group 1, Group 2; summary\\(\\) checks them")
  df <- as.data.frame(res)
  expect_named(df, c("threshold", "model", "subgroup", "n", "NB_model", "NB_all", "NB_none",
                     "delta_NB", "net_FP_avoided"))
  expect_equal(levels(df$subgroup), c("All", "Group 1", "Group 2"))
  expect_equal(nrow(df), 2 * 3)
  expect_true("NB_model_per_1000" %in% names(as.data.frame(res, per = 1000)))
  expect_s3_class(plot(res, type = "subgroups"), "ggplot")
  expect_equal(levels(plot(res, type = "subgroups")$data$label), c("Group 1", "Group 2", "Population"))
  p <- plot(res, type = "region", threshold = 0.0625, k = c(`Group 1` = 60, `Group 2` = 45))
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "inside it")
  expect_error(plot(res, type = "region"), "`threshold` must be one of")
  expect_error(plot(nb_val(d, "p", "y", thresholds = 0.1), type = "subgroups"),
               "needs subgroups")
  expect_error(c(res, nb_val(d, "p", "y", thresholds = c(0.0625, 0.1), label = "B")),
               "Results with subgroups cannot be combined")
  # survival: the estimator within each subgroup
  s <- rott_split()
  v <- s$val
  v$p <- pmin(0.9, 0.1 + v$nodes / 40)
  km <- nb_val(v, "p", surv_outcome, horizon = 5, thresholds = c(0.2, 0.4), subgroup = "meno")
  m0 <- suppressWarnings(nb_val(v[v$meno == 0, ], "p", surv_outcome, horizon = 5,
                                thresholds = c(0.2, 0.4)))
  e <- km$subgroups$estimates
  expect_equal(e$NB_model[e$subgroup == "0"], m0$estimates$NB_model)
  ci <- nb_val(v, "p", surv_outcome, horizon = 5, thresholds = c(0.2, 0.4), subgroup = "meno",
               cipcw = ~ age + nodes)
  expect_s3_class(summary(ci), "summary.voipred_nb")
})
