# Makes voi_reference.rds: reference EVPI and P(useful) from predtools::evpi_val()
# (resplab/predtools, the code of Sadatsafavi et al. 2023), for test-voi.R.
# Not run by the tests; kept to show how the reference values were made.
# Needs predtools (remotes::install_github("resplab/predtools")).
library(predtools)
data(gusto, package = "voipred.dev")
g <- gusto
us <- which(g$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15))
set.seed(123)
val_idx <- us[sample(length(us), 500)]
dev_idx <- setdiff(seq_len(nrow(g)), us)
g$kill <- (as.numeric(g$Killip) > 1) * 1
fit <- glm(day30 ~ age + miloc + pmi + kill + pmin(sysbp, 100) + pulse,
           data = g[dev_idx, ], family = binomial)
val <- g[val_idx, ]
pi <- predict(fit, val, type = "response")
zs <- c(0.017, 0.023, 0.05, 0.1, 0.15, 0.2, 0.3)
stopifnot(!any(pi %in% zs))
out <- list(val_idx = val_idx, zs = zs, pi = unname(pi))
for (m in c("bootstrap", "bayesian_bootstrap")) {
  set.seed(42)
  out[[m]] <- evpi_val(val$day30, pi, method = m, n_sim = 1000, zs = zs)
}
set.seed(1)
out$asymptotic <- evpi_val(val$day30, pi, method = "asymptotic", zs = zs)
saveRDS(out, "voi_reference.rds")
