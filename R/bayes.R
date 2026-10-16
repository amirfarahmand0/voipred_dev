# nb() works in a Bayesian framework; by default (improper prior) it uses the
# data only and gives the usual net benefit. A binary outcome has conjugate
# Beta (Dirichlet) posteriors for the prevalence and, at each threshold, the
# sensitivity and specificity of the models (Wynants et al. 2018, Cruz et al. 2024 (bayesDCA),
# Sadatsafavi et al. 2025 (EVSI)). With one model this is bayesDCA's binary
# model; with several, the models are drawn jointly (a Dirichlet over the
# patterns of treat / do not treat across the models), so that their
# correlation on the same patients is kept, and each model's sensitivity and
# specificity still has its Beta prior. Survival and competing-risks outcomes
# use the improper prior, whose posterior is sampled by the Bayesian bootstrap.

#' Priors for the net benefit
#'
#' The prior of [nb()]. By default [nb()] uses no prior information
#' (`nb_improper()`: the data only, the usual net benefit); for a binary
#' outcome, prior information can be added with the other priors. The priors
#' are Beta priors on the prevalence of the outcome and, at every threshold,
#' on the sensitivity and the specificity of the model (the probabilities of being treated among patients with the
#' event, and of not being treated among patients without it). The net benefit
#' at threshold \eqn{z} is
#' \deqn{NB = \pi \, Se - (1 - \pi)(1 - Sp) \frac{z}{1 - z},}
#' and that of treating all \eqn{\pi - (1 - \pi) z / (1 - z)}. The posteriors
#' are Beta distributions (conjugate), so their draws are exact (no MCMC).
#'
#' * `nb_improper()` (default of [nb()]): no prior information. The improper
#'   (Haldane) prior, Beta(0, 0) for all three: the data only. The net benefit
#'   of [nb()] is then the usual net benefit, \eqn{TP/n - FP/n \cdot z/(1-z)}.
#'   It is the only prior for survival and competing-risks outcomes, whose net
#'   benefit is the usual Kaplan-Meier, Aalen-Johansen or cIPCW estimate.
#' * `nb_uniform()`: Beta(1, 1) for all three, as in bayesDCA (binary outcome).
#' * `nb_jeffreys()`: Jeffreys prior, Beta(1/2, 1/2) for all three.
#' * `nb_beta()`: Beta priors with the shapes you give, e.g. an informative
#'   prior on the prevalence from earlier studies. The other priors are
#'   special cases of it.
#'
#' With the improper prior, a posterior is degenerate where a count is 0 (e.g.
#' no patient with the event is left untreated at a threshold): it is then the
#' observed proportion, as in the usual net benefit.
#'
#' @param prevalence,sensitivity,specificity The two shape parameters
#'   \eqn{(a, b)} of each Beta prior, non-negative numbers: prior mean
#'   \eqn{a / (a + b)}, worth \eqn{a + b} patients. Default Beta(1, 1).
#' @return A prior, an object of class `nb_prior`.
#' @references Cruz GNF, Korthauer K. Bayesian decision curve analysis with
#'   bayesDCA. arXiv 2023. Rubin DB. The Bayesian bootstrap. Annals of
#'   Statistics 1981; 9(1): 130-134. Sadatsafavi M, et al. Expected value of
#'   sample information calculations for risk prediction model validation.
#'   Medical Decision Making 2025; 45(3): 232-244.
#' @seealso [nb()]
#' @examples
#' nb_improper()
#' nb_uniform()
#' nb_jeffreys()
#' # an informative prior on the prevalence: mean 0.07, worth 100 patients
#' nb_beta(prevalence = c(7, 93))
#' @name nb_priors
NULL

#' @rdname nb_priors
#' @export
nb_beta <- function(prevalence = c(1, 1), sensitivity = c(1, 1), specificity = c(1, 1)) {
  for (nm in c("prevalence", "sensitivity", "specificity")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 2 || anyNA(v) || any(!is.finite(v)) || any(v < 0))
      stop("`", nm, "` must be two non-negative numbers: the shapes (a, b) of a Beta prior.",
           call. = FALSE)
  }
  ab <- list(prevalence = as.numeric(prevalence), sensitivity = as.numeric(sensitivity),
             specificity = as.numeric(specificity))
  same <- function(v) all(vapply(ab, function(x) identical(x, c(v, v)), logical(1)))
  type <- if (same(0)) "improper" else if (same(1)) "uniform" else if (same(0.5)) "jeffreys"
          else "beta"
  structure(c(list(type = type), ab), class = "nb_prior")
}

#' @rdname nb_priors
#' @export
nb_improper <- function() nb_beta(c(0, 0), c(0, 0), c(0, 0))

#' @rdname nb_priors
#' @export
nb_uniform <- function() nb_beta(c(1, 1), c(1, 1), c(1, 1))

#' @rdname nb_priors
#' @export
nb_jeffreys <- function() nb_beta(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5))

#' @export
print.nb_prior <- function(x, ...) {
  cat("Prior:", prior_text(x, short = FALSE), "\n")
  invisible(x)
}

# The prior in words: long (print) or short (headers of results).
prior_text <- function(prior, short = TRUE) {
  b <- function(v) sprintf("Beta(%g, %g)", v[1], v[2])
  switch(prior$type,
    improper = if (short) "improper Beta(0, 0)" else
      "improper (Haldane), Beta(0, 0) for prevalence, sensitivity and specificity: the data only",
    uniform = if (short) "uniform Beta(1, 1)" else
      "uniform, Beta(1, 1) for prevalence, sensitivity and specificity",
    jeffreys = if (short) "Jeffreys Beta(0.5, 0.5)" else
      "Jeffreys, Beta(0.5, 0.5) for prevalence, sensitivity and specificity",
    sprintf("Beta: prevalence %s, sensitivity %s, specificity %s",
            b(prior$prevalence), b(prior$sensitivity), b(prior$specificity)))
}

# ---- internal ----------------------------------------------------------------

# Treated: risk > 0 and risk >= z (as in the usual net benefit).
treated_at <- function(p, z) p > 0 & p >= z

# Posterior mean net benefit of one model (the point estimate of nb()), with
# the posterior means of the prevalence, sensitivity and specificity. Case
# weights (improper prior only) weight the counts. With the improper prior a
# proportion with an empty denominator (e.g. no events) is multiplied by a
# prevalence (or 1 - prevalence) of 0, so any value gives the same net benefit.
bayes_nb_mean <- function(y, p, thresholds, prior, w = NULL) {
  if (is.null(w)) w <- rep(1, length(y))
  n <- sum(w)
  d <- sum(w * y)
  odds <- thresholds / (1 - thresholds)
  a <- prior$prevalence
  ep <- (a[1] + d) / (sum(a) + n)
  tp <- vapply(thresholds, function(z) sum(w * (y == 1 & treated_at(p, z))), numeric(1))
  tn <- vapply(thresholds, function(z) sum(w * (y == 0 & !treated_at(p, z))), numeric(1))
  se <- (prior$sensitivity[1] + tp) / (sum(prior$sensitivity) + d)
  sp <- (prior$specificity[1] + tn) / (sum(prior$specificity) + n - d)
  se[!is.finite(se)] <- 0
  sp[!is.finite(sp)] <- 1
  list(NB_model = ep * se - (1 - ep) * (1 - sp) * odds,
       NB_all = ep - (1 - ep) * odds, prevalence = ep, sensitivity = se, specificity = sp)
}

# Data and prior of the joint posterior of several models validated on the same
# patients: for each threshold, the counts of the 2^M patterns of
# (treated by model 1, ..., treated by model M) among patients with (c1) and
# without (c0) the event, and the Dirichlet prior of the patterns (alpha1,
# alpha0), chosen so that each model's sensitivity and specificity has its
# Beta prior.
bayes_setup <- function(res) {
  v1 <- voi_info(res[[1]])
  prior <- v1$prior
  y <- v1$y
  z <- v1$thresholds
  M <- length(res)
  K <- 2^M
  pos <- outer(0:(K - 1), 0:(M - 1), function(cell, m) (cell %/% 2^m) %% 2 == 1)  # K x M
  risk <- lapply(res, function(r) r$risk)
  cell <- vapply(z, function(zz) {
    out <- numeric(length(y))
    for (m in seq_len(M)) out <- out + treated_at(risk[[m]], zz) * 2^(m - 1)
    out
  }, numeric(length(y)))                                         # n x Z
  if (!is.matrix(cell)) cell <- matrix(cell, nrow = length(y))
  count <- function(rows) t(apply(cell[rows, , drop = FALSE], 2, function(cc)
    tabulate(cc + 1, nbins = K)))                                 # Z x K
  c1 <- count(y == 1)
  c0 <- count(y == 0)
  dirichlet_prior <- function(ab, positive_is_a) {
    s <- sum(ab)
    if (s == 0) return(rep(0, K))                  # improper prior
    pa <- ab[1] / s
    q <- ifelse(pos, if (positive_is_a) pa else 1 - pa, if (positive_is_a) 1 - pa else pa)
    s * apply(q, 1, prod)
  }
  list(y = y, n = length(y), d = sum(y), z = z, odds = z / (1 - z), M = M, K = K,
       pos = pos, c1 = c1, c0 = c0, prior = prior,
       alpha1 = dirichlet_prior(prior$sensitivity, TRUE),   # treated among events
       alpha0 = dirichlet_prior(prior$specificity, FALSE))  # not treated among non-events
}

# Dirichlet draws (one row each). A row with no data and an improper prior
# (all shapes 0) belongs to an empty group, e.g. the events when there are
# none: its probabilities are multiplied by a prevalence of 0, so they are set
# to 0.
rdirichlet_rows <- function(n, alpha) {
  g <- matrix(stats::rgamma(n * length(alpha), shape = rep(alpha, each = n)), n)
  out <- g / rowSums(g)
  out[!is.finite(out)] <- 0
  out
}

# Draws from the posterior: the prevalence (n_sim), and for each threshold the
# sensitivity and specificity of each model (n_sim x Z x M). With one model
# the draws are made as in bayesDCA (rbeta, in its order).
bayes_params <- function(s, n_sim) {
  a <- s$prior$prevalence
  prev <- stats::rbeta(n_sim, a[1] + s$d, a[2] + s$n - s$d)
  Z <- length(s$z)
  se <- sp <- array(NA_real_, c(n_sim, Z, s$M))
  th1 <- th0 <- NULL
  if (s$M > 1) th1 <- th0 <- array(NA_real_, c(n_sim, Z, s$K))
  for (j in seq_len(Z)) {
    if (s$M == 1) {
      se[, j, 1] <- stats::rbeta(n_sim, s$alpha1[2] + s$c1[j, 2], s$alpha1[1] + s$c1[j, 1])
      sp[, j, 1] <- stats::rbeta(n_sim, s$alpha0[1] + s$c0[j, 1], s$alpha0[2] + s$c0[j, 2])
    } else {
      t1 <- rdirichlet_rows(n_sim, s$alpha1 + s$c1[j, ])
      t0 <- rdirichlet_rows(n_sim, s$alpha0 + s$c0[j, ])
      th1[, j, ] <- t1
      th0[, j, ] <- t0
      se[, j, ] <- t1 %*% s$pos
      sp[, j, ] <- t0 %*% (!s$pos)
    }
  }
  list(prev = prev, se = se, sp = sp, th1 = th1, th0 = th0)
}

# Net benefit of every strategy from prevalence, sensitivity and specificity:
# an array draw x threshold x strategy (models, treat all, treat none).
bayes_nb_array <- function(prev, se, sp, odds) {
  B <- length(prev)
  Z <- length(odds)
  M <- dim(se)[3]
  out <- array(0, c(B, Z, M + 2))
  oz <- matrix(odds, B, Z, byrow = TRUE)
  for (m in seq_len(M))
    out[, , m] <- prev * se[, , m] - (1 - prev) * (1 - sp[, , m]) * oz
  out[, , M + 1] <- prev - (1 - prev) * oz
  out
}

# nb_draws() for a binary outcome without case weights: draws from the posterior.
bayes_draws <- function(res, n_sim) {
  s <- bayes_setup(res)
  th <- bayes_params(s, n_sim)
  bayes_nb_array(th$prev, th$se, th$sp, s$odds)
}

# Multinomial counts (one row per draw) with per-draw sizes and probabilities,
# by conditional binomials.
rmultinom_rows <- function(size, prob) {
  K <- ncol(prob)
  out <- matrix(0, nrow(prob), K)
  left <- size
  rest <- rep(1, nrow(prob))
  for (k in seq_len(K - 1)) {
    pk <- ifelse(rest > 0, pmin(1, pmax(0, prob[, k] / rest)), 0)
    out[, k] <- stats::rbinom(nrow(prob), left, pk)
    left <- left - out[, k]
    rest <- rest - prob[, k]
  }
  out[, K] <- left
  out
}

# EVSI with the conjugate model (Sadatsafavi et al. 2025, EVSI()): for each
# posterior draw of the population, simulate a future sample of N patients
# (events, and treated / not treated among them), update the posterior means
# with it, choose the strategy with the highest posterior mean net benefit,
# and record its true net benefit in the draw.
bayes_evsi <- function(res, future_n, n_sim) {
  s <- bayes_setup(res)
  th <- bayes_params(s, n_sim)
  truth <- bayes_nb_array(th$prev, th$se, th$sp, s$odds)
  Z <- length(s$z)
  current <- apply(apply(truth, c(2, 3), mean), 1, max)
  perfect <- colMeans(matrix(apply(truth, c(1, 2), max), nrow = n_sim))
  a <- s$prior$prevalence
  evsi <- matrix(NA_real_, Z, length(future_n))
  for (i in seq_along(future_n)) {
    N <- future_n[i]
    D <- stats::rbinom(n_sim, N, th$prev)
    ep2 <- (a[1] + s$d + D) / (sum(a) + s$n + N)
    se2 <- sp2 <- array(NA_real_, c(n_sim, Z, s$M))
    for (j in seq_len(Z)) {
      if (s$M == 1) {
        TP <- stats::rbinom(n_sim, D, th$se[, j, 1])
        TN <- stats::rbinom(n_sim, N - D, th$sp[, j, 1])
        se2[, j, 1] <- finite_or((s$alpha1[2] + s$c1[j, 2] + TP) / (sum(s$alpha1) + s$d + D), 0)
        sp2[, j, 1] <- finite_or((s$alpha0[1] + s$c0[j, 1] + TN) /
                                   (sum(s$alpha0) + s$n - s$d + N - D), 1)
      } else {
        X1 <- rmultinom_rows(D, th$th1[, j, ])
        X0 <- rmultinom_rows(N - D, th$th0[, j, ])
        m1 <- finite_or(sweep(X1, 2, s$alpha1 + s$c1[j, ], `+`) / (sum(s$alpha1) + s$d + D), 0)
        m0 <- finite_or(sweep(X0, 2, s$alpha0 + s$c0[j, ], `+`) /
                          (sum(s$alpha0) + s$n - s$d + N - D), 0)
        se2[, j, ] <- m1 %*% s$pos
        sp2[, j, ] <- m0 %*% (!s$pos)
      }
    }
    pooled <- bayes_nb_array(ep2, se2, sp2, s$odds)
    evsi[, i] <- colMeans(winner_value(pooled, truth)) - current
  }
  list(EVSI = evsi, EVPI = perfect - current)
}

# x with its non-finite values (0/0 of an empty group under the improper prior)
# replaced by `value`.
finite_or <- function(x, value) {
  x[!is.finite(x)] <- value
  x
}
