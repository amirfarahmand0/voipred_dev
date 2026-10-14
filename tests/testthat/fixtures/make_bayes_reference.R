# Makes bayes_reference.rds: reference draws of the Bayesian net benefit from
# bayesDCA (giulianonetto/bayesDCA, R/dca.R, .dca_binary()) and EVSI from the
# conjugate algorithm of Sadatsafavi et al. 2025 (resplab/evsiexval, R/core.R,
# EVSI()), on the GUSTO-I case study of voi_reference.rds, for test-bayes.R.
# Not run by the tests; kept to show how the reference values were made.
# Run from tests/testthat with the two repositories cloned next to it:
#   git clone https://github.com/giulianonetto/bayesDCA
#   git clone https://github.com/resplab/evsiexval
library(magrittr)
library(progress)
source_fun <- function(file, name) {
  ex <- parse(file)
  nm <- vapply(ex, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")))
    as.character(e[[2]]) else "", "")
  eval(ex[[which(nm == name)]], globalenv())
}
source_fun("bayesDCA/R/dca.R", ".dca_binary")
source_fun("evsiexval/R/core.R", "EVSI")

ref0 <- readRDS("fixtures/voi_reference.rds")
data(gusto, package = "voipred.dev")
g <- gusto
g$kill <- (as.numeric(g$Killip) > 1) * 1
us <- which(g$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15))
fit <- glm(day30 ~ age + miloc + pmi + kill + pmin(sysbp, 100) + pulse,
           data = g[setdiff(seq_len(nrow(g)), us), ], family = binomial)
val <- g[ref0$val_idx, ]
pi <- unname(predict(fit, val, type = "response"))
y <- val$day30
zs <- c(0.017, 0.023, 0.05, 0.1, 0.2)
stopifnot(!any(pi %in% zs))

# bayesDCA: uniform priors, one model
tp <- matrix(sapply(zs, function(t) sum(pi[y == 1] >= t)), ncol = 1, dimnames = list(NULL, "m"))
tn <- matrix(sapply(zs, function(t) sum(pi[y == 0] < t)), ncol = 1, dimnames = list(NULL, "m"))
one <- matrix(1, length(zs), 1)
set.seed(2024)
bd <- .dca_binary(length(zs), "m", length(y), sum(y), tp, tn, zs, 1, 1, one, one, one, one,
                  n_draws = 500)
out <- list(zs = zs, seed = 2024, n_draws = 500,
            net_benefit = bd$distributions$net_benefit$m,
            treat_all = bd$distributions$treat_all,
            p_useful = bd$comparisons$p_useful[, 1])

# EVSI(): one threshold at a time, uniform priors; it subtracts the posterior
# mean (closed form) of the current best strategy
fss <- c(200, 1000, 5000)
out$evsi <- lapply(c(0.023, 0.05), function(z) {
  set.seed(77)
  r <- EVSI(NULL, data.frame(Y = y, pi = pi), z, fss, n_sim = 10000)
  list(z = z, future_n = fss, n_sim = 10000, seed = 77, EVPI = r$EVPI, EVSI = r$EVSI)
})
saveRDS(out, "fixtures/bayes_reference.rds")
