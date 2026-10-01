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
