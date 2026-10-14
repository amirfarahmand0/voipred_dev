# Reference code of Sadatsafavi et al. for the value-of-information tests.

# The GUSTO-I case study of the papers: model developed outside the US,
# validated on 500 US patients (the indices are saved in the fixture).
voi_case <- function() {
  ref <- readRDS(test_path("fixtures", "voi_reference.rds"))
  g <- gusto
  g$kill <- (as.numeric(g$Killip) > 1) * 1
  us <- which(g$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15))
  fit <- glm(day30 ~ age + miloc + pmi + kill + pmin(sysbp, 100) + pulse,
             data = g[setdiff(seq_len(nrow(g)), us), ], family = binomial)
  val <- g[ref$val_idx, ]
  list(ref = ref, fit = fit, val = val, pi = unname(predict(fit, val, type = "response")))
}

# predtools:::bootstrap() (resplab/predtools): the weights of one draw.
predtools_bootstrap <- function(n, Bayesian = FALSE) {
  if (Bayesian) {
    u <- c(0, sort(stats::runif(n - 1)), 1)
    (u[-1] - u[-length(u)]) * n
  } else {
    as.vector(stats::rmultinom(1, n, rep(1 / n, n)))
  }
}

# EVSI by bootstrap: do_sim("BBS") and do_sim("OBS") of
# resplab/papercode/evsiexval/gusto_evsi_compare.R (one model, one threshold).
reference_evsi <- function(Y, pi, z, future_sample_sizes, n_sim, method = "BBS") {
  n <- length(Y)
  NB1 <- NB2 <- NBbest <- 0
  NBw <- rep(0, length(future_sample_sizes))
  for (i in 1:n_sim) {
    if (method == "BBS") {
      W <- stats::rgamma(n, 1, 1)
      W <- W / sum(W)
    } else {
      W <- stats::rmultinom(1, n, rep(1 / n, n)) / n
    }
    prev <- sum(Y * W)
    se <- sum(Y * (pi >= z) * W) / sum(Y * W)
    sp <- sum((1 - Y) * (pi < z) * W) / sum((1 - Y) * W)
    NB1t <- prev * se - (1 - prev) * (1 - sp) * z / (1 - z)
    NB2t <- prev - (1 - prev) * z / (1 - z)
    NBbest <- NBbest + max(0, NB1t, NB2t)
    NB1 <- NB1 + NB1t
    NB2 <- NB2 + NB2t
    for (j in seq_along(future_sample_sizes)) {
      W_pooled <- stats::rmultinom(1, future_sample_sizes[j], W) + 1
      prev_p <- sum(Y * W_pooled) / sum(W_pooled)
      se_p <- sum(Y * (pi >= z) * W_pooled) / sum(Y * W_pooled)
      sp_p <- sum((1 - Y) * (pi < z) * W_pooled) / sum((1 - Y) * W_pooled)
      NB1_p <- prev_p * se_p - (1 - prev_p) * (1 - sp_p) * z / (1 - z)
      NB2_p <- prev_p - (1 - prev_p) * z / (1 - z)
      NBw[j] <- NBw[j] + c(0, NB1t, NB2t)[which.max(c(0, NB1_p, NB2_p))]
    }
  }
  current <- max(0, NB1 / n_sim, NB2 / n_sim)
  list(EVPI = NBbest / n_sim - current, EVSI = NBw / n_sim - current)
}

# Net benefit draws for given patient weights W (n x draws), for checks against
# direct calculations (the draws of the VoI functions are always random).
fixed_draws <- function(..., W) {
  a <- list(...)
  draws_from(collect_voi(a[[1]], a[-1], NULL, "fixed_draws()", "`W`"), W, "fixed")
}

# Random draws for results of nb(), as the VoI functions make them (internal
# nb_draws()), so several measures can be checked on the same draws.
draws_t <- function(..., n_sim = 1000, method = "bayesian") {
  a <- list(...)
  res <- collect_voi(a[[1]], a[-1], NULL, "draws_t()", "`n_sim` and `method`")
  nb_draws(res, n_sim, method, method_given = !missing(method))
}
