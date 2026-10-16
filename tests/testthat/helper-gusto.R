# Shared GUSTO-I split for the model-based tests: develop on 5,000 patients,
# validate on another 5,000.
gusto_split <- function(seed = 2026) {
  set.seed(seed)
  i <- sample(nrow(gusto))
  list(dev = gusto[i[1:5000], ], val = gusto[i[5001:10000], ])
}

gusto_formula <- day30 ~ age + Killip + sysbp + pulse + pmi + miloc + sex

# All thresholds a user would normally look at.
all_z <- seq(0, 0.99, by = 0.01)

# Net benefit from dcurves, as one data frame with columns all, none, model.
dcurves_nb <- function(risk, data, outcome, thresholds) {
  d <- data.frame(y = data[[outcome]], risk = risk)
  res <- suppressMessages(dcurves::dca(y ~ risk, data = d, thresholds = thresholds))
  tab <- as.data.frame(res$dca)
  pick <- function(v) tab$net_benefit[tab$variable == v][order(tab$threshold[tab$variable == v])]
  data.frame(threshold = sort(unique(tab$threshold)),
             NB_model = pick("risk"), NB_all = pick("all"), NB_none = pick("none"))
}

# nb() result vs dcurves::dca() on the same predicted risks: model, treat all
# and treat none at every threshold.
expect_same_as_dcurves <- function(res, val) {
  ref <- dcurves_nb(res$risk, val, "day30", all_z)
  expect_equal(res$estimates$threshold, ref$threshold)
  expect_equal(res$estimates$NB_model, ref$NB_model, tolerance = 1e-10)
  expect_equal(res$estimates$NB_all, ref$NB_all, tolerance = 1e-10)
  expect_equal(res$estimates$NB_none, ref$NB_none, tolerance = 1e-10)
}

# For models trained on a numeric matrix: the same GUSTO split, with the
# covariates dummy-coded as columns (KillipII, sexfemale, ...).
gusto_matrix_split <- function(seed = 2026) {
  s <- gusto_split(seed)
  to_df <- function(d) {
    x <- stats::model.matrix(gusto_formula, d)[, -1]
    data.frame(x, day30 = d$day30, check.names = FALSE)
  }
  dev <- to_df(s$dev)
  list(dev = dev, val = to_df(s$val),
       x = as.matrix(dev[, setdiff(names(dev), "day30")]), y = dev$day30)
}
