# nb_val() works in a Bayesian framework; by default (improper prior) it uses the
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
#' The prior of [nb_val()]. By default [nb_val()] uses no prior information
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
#' * `nb_improper()` (default of [nb_val()]): no prior information. The improper
#'   (Haldane) prior, Beta(0, 0) for all three: the data only. The net benefit
#'   of [nb_val()] is then the usual net benefit, \eqn{TP/n - FP/n \cdot z/(1-z)}.
#'   It is the only prior for survival and competing-risks outcomes, whose net
#'   benefit is the usual Kaplan-Meier, Aalen-Johansen or cIPCW estimate.
#' * `nb_uniform()`: Beta(1, 1) for all three, as in bayesDCA (binary outcome).
#' * `nb_jeffreys()`: Jeffreys prior, Beta(1/2, 1/2) for all three.
#' * `nb_beta()`: Beta priors with the shapes you give, e.g. an informative
#'   prior on the prevalence from earlier studies. The priors above are
#'   special cases of it.
#' * [nb_informative()]: the informative prior of bayesDCA, whose sensitivity
#'   and specificity priors change with the threshold; for small effective
#'   sample sizes (see its help page).
#'
#' All these priors except [nb_informative()] are the same at every threshold.
#' For a binary outcome, the result of [nb_val()] holds the prior shapes used
#' at each of its thresholds (`prior_shapes`).
#'
#' With the improper prior, a posterior is degenerate where a count is 0 (e.g.
#' no patient with the event is left untreated at a threshold): it is then the
#' observed proportion, as in the usual net benefit.
#'
#' @param prevalence,sensitivity,specificity The two shape parameters
#'   \eqn{(a, b)} of each Beta prior, non-negative numbers: prior mean
#'   \eqn{a / (a + b)}, worth \eqn{a + b} patients. Default Beta(0, 0): no
#'   prior information, so only the parts you give carry prior information
#'   (e.g. `nb_beta(prevalence = c(7, 93))` adds information on the
#'   prevalence only); `nb_beta()` is `nb_improper()`.
#' @return A prior, an object of class `nb_prior`.
#' @references Cruz GNF, Korthauer K. Bayesian decision curve analysis with
#'   bayesDCA. arXiv 2023. Rubin DB. The Bayesian bootstrap. Annals of
#'   Statistics 1981; 9(1): 130-134. Sadatsafavi M, et al. Expected value of
#'   sample information calculations for risk prediction model validation.
#'   Medical Decision Making 2025; 45(3): 232-244.
#' @seealso [nb_val()]
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
nb_beta <- function(prevalence = c(0, 0), sensitivity = c(0, 0), specificity = c(0, 0)) {
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

#' Informative prior of bayesDCA (threshold-varying)
#'
#' The informative prior of bayesDCA (Cruz and Korthauer, 2025), for a binary
#' outcome: Beta priors on the sensitivity and specificity of each model that
#' change with the threshold, and a Beta prior on the prevalence. It reproduces
#' bayesDCA's threshold-varying prior (`dca(threshold_varying_prior = TRUE)`),
#' with the same defaults. Give it to [nb_val()] as `prior = nb_informative()`;
#' the value-of-information functions then use it too.
#'
#' **The logic.** Under the rule "treat if the predicted risk is at or above
#' the threshold", every model treats everyone at a threshold of 0
#' (sensitivity 1, specificity 0: treat all) and no one at a threshold above
#' all predicted risks (sensitivity 0, specificity 1: treat none). Near the
#' ends of the decision curve few patients are untreated (low thresholds) or
#' treated (high thresholds), so the data say little there, and a vague prior
#' lets the net benefit, and the value-of-information measures, swing widely.
#' This prior states, with little weight, that the model behaves like treat
#' all at the lowest threshold and like treat none at the highest, and stays
#' vague in between. Cruz and Korthauer (2025) found that such a prior keeps
#' EVPI decreasing with the sample size in small samples, where it otherwise
#' may not.
#'
#' **The construction** (as in bayesDCA). Each Beta prior is written as a prior
#' mean \eqn{\mu} and a prior strength \eqn{s} (worth \eqn{s} patients):
#' Beta\eqn{(\mu s, (1 - \mu) s)}. With \eqn{t_{min}} and \eqn{t_{max}} the
#' lowest and highest thresholds given to [nb_val()] and the vague
#' ("ignorance") region from `ignorance[1]` \eqn{\times t_{max}} to
#' `ignorance[2]` \eqn{\times t_{max}}, the prior mean of the sensitivity is
#' `sensitivity[1]` at \eqn{t_{min}}, falls linearly to the mean of
#' `ignorance_prior` at the start of the vague region, stays there through it,
#' and falls linearly to `sensitivity[2]` at \eqn{t_{max}}; its strength falls
#' linearly from `strength` to that of `ignorance_prior` and rises back to
#' `strength` in the same way. Without a vague region (`ignorance = NULL`) the
#' mean falls linearly from `sensitivity[1]` to `sensitivity[2]` across all
#' thresholds, with strength `strength`. Prior means are kept within
#' \[0.001, 0.999\], and thresholds are taken within \[1e-9, 0.99\], as in
#' bayesDCA. The prior of the specificity has mean 1 minus that of the
#' sensitivity and the same strength. The priors at different thresholds are
#' independent; each model has the same prior. The prevalence has the Beta
#' prior `prevalence` at every threshold.
#'
#' The prior at a threshold depends on the thresholds given to [nb_val()]
#' (through \eqn{t_{min}}, \eqn{t_{max}} and the vague region), as in
#' bayesDCA. It is for binary outcomes only, without case weights. With
#' reasonably large samples (e.g. 100 events or more) it changes the results
#' little.
#'
#' @param sensitivity Prior mean of the sensitivity at the lowest and at the
#'   highest threshold. Default `c(0.99, 0.01)`.
#' @param strength Prior strength (in patients) of the sensitivity and
#'   specificity priors at the lowest and highest thresholds. Default 5.
#' @param ignorance The vague region, as fractions of the highest threshold:
#'   default `c(0.25, 0.75)`. `NULL`: no vague region.
#' @param ignorance_prior The Beta prior (shapes \eqn{(a, b)}) of the
#'   sensitivity in the vague region. Default `c(1, 1)` (uniform: prior mean
#'   0.5, strength 2).
#' @param prevalence The Beta prior (shapes \eqn{(a, b)}) of the prevalence.
#'   Default `c(1, 1)` (uniform).
#' @return A prior, an object of class `nb_prior`, for [nb_val()].
#' @references Cruz GNF, Korthauer K. Bayesian decision curve analysis with
#'   bayesDCA. Statistics in Medicine 2025.
#' @seealso [nb_priors], [nb_val()]
#' @examples
#' nb_informative()
#' df_dev <- gusto[1:5000, ]
#' df_val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = df_dev)
#' df_val$p <- predict(fit, newdata = df_val, type = "response")
#' res <- nb_val(df_val, "p", "day30", thresholds = seq(0, 0.3, by = 0.01),
#'               prior = nb_informative())
#' res
#' # the prior used at each threshold
#' res$prior_shapes[c(1, 6, 11, 16, 21, 26, 31), ]
#' @export
nb_informative <- function(sensitivity = c(0.99, 0.01), strength = 5,
                           ignorance = c(0.25, 0.75), ignorance_prior = c(1, 1),
                           prevalence = c(1, 1)) {
  if (!is.numeric(sensitivity) || length(sensitivity) != 2 || anyNA(sensitivity) ||
      any(sensitivity < 0 | sensitivity > 1))
    stop("`sensitivity` must be two prior means in [0, 1]: at the lowest and at the ",
         "highest threshold.", call. = FALSE)
  if (!is.numeric(strength) || length(strength) != 1 || is.na(strength) ||
      !is.finite(strength) || strength <= 0)
    stop("`strength` must be a single positive number (a number of patients).",
         call. = FALSE)
  if (!is.null(ignorance) && (!is.numeric(ignorance) || length(ignorance) != 2 ||
      anyNA(ignorance) || any(ignorance < 0 | ignorance > 1) || ignorance[1] > ignorance[2]))
    stop("`ignorance` must be two increasing fractions of the highest threshold, in ",
         "[0, 1] (e.g. c(0.25, 0.75)), or NULL for no vague region.", call. = FALSE)
  for (nm in c("ignorance_prior", "prevalence")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 2 || anyNA(v) || any(!is.finite(v)) || any(v < 0))
      stop("`", nm, "` must be two non-negative numbers: the shapes (a, b) of a Beta ",
           "prior.", call. = FALSE)
  }
  if (sum(ignorance_prior) <= 0)
    stop("`ignorance_prior` must have a positive total (a + b > 0).", call. = FALSE)
  structure(list(type = "informative", sensitivity = as.numeric(sensitivity),
                 strength = as.numeric(strength),
                 ignorance = if (!is.null(ignorance)) as.numeric(ignorance),
                 ignorance_prior = as.numeric(ignorance_prior),
                 prevalence = as.numeric(prevalence)),
            class = "nb_prior")
}

# The shapes (a, b) of the Beta priors of the prevalence and, at each threshold,
# of the sensitivity and specificity, as nb_val() uses them: a data frame with
# one row per threshold (threshold, se_a, se_b, sp_a, sp_b, prev_a, prev_b).
# Stored in the result of nb_val() (binary outcome) as `prior_shapes`.
prior_shapes <- function(prior, thresholds) {
  if (!inherits(prior, "nb_prior"))
    stop("`prior` must be a prior, e.g. nb_informative().", call. = FALSE)
  thresholds <- check_thresholds(thresholds)
  s <- threshold_shapes(prior, thresholds)
  data.frame(threshold = thresholds, se_a = s$se[, 1], se_b = s$se[, 2],
             sp_a = s$sp[, 1], sp_b = s$sp[, 2],
             prev_a = prior$prevalence[1], prev_b = prior$prevalence[2])
}

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
    informative = if (short) "informative (bayesDCA, threshold-varying)" else
      paste0("informative (bayesDCA): threshold-varying Beta priors on sensitivity ",
             "and specificity, sensitivity mean ", prior$sensitivity[1],
             " at the lowest threshold to ", prior$sensitivity[2], " at the highest, ",
             "strength ", prior$strength, " at the ends",
             if (is.null(prior$ignorance)) " (no vague region)"
             else paste0(", ", b(prior$ignorance_prior), " from ",
                         100 * prior$ignorance[1], "% to ", 100 * prior$ignorance[2],
                         "% of the highest threshold"),
             "; prevalence ", b(prior$prevalence)),
    sprintf("Beta: prevalence %s, sensitivity %s, specificity %s",
            b(prior$prevalence), b(prior$sensitivity), b(prior$specificity)))
}

# ---- internal ----------------------------------------------------------------

# Treated: risk > 0 and risk >= z (as in the usual net benefit).
treated_at <- function(p, z) p > 0 & p >= z

# Posterior mean net benefit of one model (the point estimate of nb_val()), with
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
  sh <- threshold_shapes(prior, thresholds)
  se <- (sh$se[, 1] + tp) / (rowSums(sh$se) + d)
  sp <- (sh$sp[, 1] + tn) / (rowSums(sh$sp) + n - d)
  se[!is.finite(se)] <- 0
  sp[!is.finite(sp)] <- 1
  list(NB_model = ep * se - (1 - ep) * (1 - sp) * odds,
       NB_all = ep - (1 - ep) * odds, prevalence = ep, sensitivity = se, specificity = sp)
}

# Data and prior of the joint posterior of several models validated on the same
# patients: for each threshold, the counts of the 2^M patterns of
# (treated by model 1, ..., treated by model M) among patients with (c1) and
# without (c0) the event, and the Dirichlet prior of the patterns at each
# threshold (alpha1, alpha0: threshold x pattern), chosen so that each model's
# sensitivity and specificity has its Beta prior at that threshold.
bayes_setup <- function(res) {
  v1 <- voi_info(res)
  prior <- v1$prior
  y <- v1$y
  z <- v1$thresholds
  M <- ncol(v1$risk)
  K <- 2^M
  pos <- outer(0:(K - 1), 0:(M - 1), function(cell, m) (cell %/% 2^m) %% 2 == 1)  # K x M
  risk <- lapply(seq_len(M), function(m) v1$risk[, m])
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
  sh <- threshold_shapes(prior, z)
  rows <- function(ab, positive_is_a)
    t(vapply(seq_along(z), function(j) dirichlet_prior(ab[j, ], positive_is_a), numeric(K)))
  alpha1 <- rows(sh$se, TRUE)                      # treated among events
  alpha0 <- rows(sh$sp, FALSE)                     # not treated among non-events
  if (K == 1) {
    alpha1 <- matrix(alpha1, ncol = 1)
    alpha0 <- matrix(alpha0, ncol = 1)
  }
  list(y = y, n = length(y), d = sum(y), z = z, odds = z / (1 - z), M = M, K = K,
       pos = pos, c1 = c1, c0 = c0, prior = prior,
       alpha1 = matrix(alpha1, nrow = length(z)), alpha0 = matrix(alpha0, nrow = length(z)))
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
      se[, j, 1] <- stats::rbeta(n_sim, s$alpha1[j, 2] + s$c1[j, 2], s$alpha1[j, 1] + s$c1[j, 1])
      sp[, j, 1] <- stats::rbeta(n_sim, s$alpha0[j, 1] + s$c0[j, 1], s$alpha0[j, 2] + s$c0[j, 2])
    } else {
      t1 <- rdirichlet_rows(n_sim, s$alpha1[j, ] + s$c1[j, ])
      t0 <- rdirichlet_rows(n_sim, s$alpha0[j, ] + s$c0[j, ])
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
  check_joint_size(s, n_sim)
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
  check_joint_size(s, n_sim)
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
        se2[, j, 1] <- finite_or((s$alpha1[j, 2] + s$c1[j, 2] + TP) /
                                   (sum(s$alpha1[j, ]) + s$d + D), 0)
        sp2[, j, 1] <- finite_or((s$alpha0[j, 1] + s$c0[j, 1] + TN) /
                                   (sum(s$alpha0[j, ]) + s$n - s$d + N - D), 1)
      } else {
        X1 <- rmultinom_rows(D, th$th1[, j, ])
        X0 <- rmultinom_rows(N - D, th$th0[, j, ])
        m1 <- finite_or(sweep(X1, 2, s$alpha1[j, ] + s$c1[j, ], `+`) /
                          (sum(s$alpha1[j, ]) + s$d + D), 0)
        m0 <- finite_or(sweep(X0, 2, s$alpha0[j, ] + s$c0[j, ], `+`) /
                          (sum(s$alpha0[j, ]) + s$n - s$d + N - D), 0)
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

# The joint posterior of several models has one probability per pattern of
# treat / do not treat across the models (2^M) at each threshold and in each
# draw, for the patients with and without the event: stop before it needs
# more memory than `limit` bytes.
check_joint_size <- function(s, n_sim, limit = 2e9) {
  if (s$M < 2) return(invisible())
  bytes <- 2 * 8 * n_sim * length(s$z) * s$K
  if (bytes > limit)
    stop("The joint posterior of ", s$M, " models needs about ",
         format(signif(bytes / 1e9, 2)), " GB of memory (2^", s$M, " = ", s$K,
         " patterns of treatment across the models, for each of ", length(s$z),
         " thresholds and ", n_sim, " draws). Use fewer models, thresholds or draws.",
         call. = FALSE)
  invisible()
}

# The Beta shapes of the sensitivity and specificity priors at each threshold:
# list(se, sp), each a matrix with one row per threshold and columns (a, b).
# The same at every threshold, except for nb_informative().
threshold_shapes <- function(prior, thresholds) {
  Z <- length(thresholds)
  if (identical(prior$type, "informative")) return(informative_shapes(prior, thresholds))
  list(se = matrix(prior$sensitivity, Z, 2, byrow = TRUE),
       sp = matrix(prior$specificity, Z, 2, byrow = TRUE))
}

# bayesDCA's threshold-varying prior (.get_threshold_varying_prior_parameters()),
# for one model.
informative_shapes <- function(prior, thresholds) {
  t <- pmax(pmin(thresholds, 0.99), 1e-9)       # as bayesDCA
  min_t <- min(t)
  max_t <- max(t)
  if (max_t <= min_t)
    stop("nb_informative() needs at least two different thresholds: its prior goes ",
         "from the lowest to the highest one.", call. = FALSE)
  hi <- prior$sensitivity[1]
  lo <- prior$sensitivity[2]
  s_end <- prior$strength
  ig_size <- sum(prior$ignorance_prior)
  ig_mean <- prior$ignorance_prior[1] / ig_size
  if (is.null(prior$ignorance)) {
    mu <- hi + (t - min_t) * (lo - hi) / (max_t - min_t)
    size <- rep(s_end, length(t))
  } else {
    cuts <- prior$ignorance * max_t
    if (cuts[1] <= min_t)
      stop("The vague region of nb_informative() starts at ", signif(cuts[1], 3),
           " (", prior$ignorance[1], " of the highest threshold), which is not above ",
           "the lowest threshold (", signif(min_t, 3), "). Give thresholds from lower ",
           "values, another `ignorance`, or ignorance = NULL.", call. = FALSE)
    mu <- size <- numeric(length(t))
    low <- t <= cuts[1]
    mid <- !low & t <= cuts[2]
    top <- !low & !mid
    mu[low] <- hi + (t[low] - min_t) * (ig_mean - hi) / (cuts[1] - min_t)
    size[low] <- s_end + (t[low] - min_t) * (ig_size - s_end) / (cuts[1] - min_t)
    mu[mid] <- ig_mean
    size[mid] <- ig_size
    mu[top] <- ig_mean + (t[top] - cuts[2]) * (lo - ig_mean) / (max_t - cuts[2])
    size[top] <- ig_size + (t[top] - cuts[2]) * (s_end - ig_size) / (max_t - cuts[2])
  }
  mu <- pmax(pmin(mu, 0.999), 0.001)
  list(se = cbind(mu * size, (1 - mu) * size),
       sp = cbind((1 - mu) * size, mu * size))
}
