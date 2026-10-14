# Value of information for model validation: EVPI, relative EVPI, EVSI,
# P(useful) and P(best).
#
# The uncertainty about the net benefit is represented by draws of the
# validation sample: each draw re-weights the patients (Bayesian bootstrap:
# Dirichlet weights; ordinary bootstrap: multinomial counts), and the net
# benefit of every strategy is recomputed with the same estimator as nb()
# (binary, Kaplan-Meier / Aalen-Johansen, or cIPCW with the censoring model
# refitted in each draw). Following Sadatsafavi et al. (2023, 2025).

#' Expected value of perfect information (EVPI) for model validation
#'
#' The expected loss in net benefit due to uncertainty about which strategy
#' (each model, treat all, treat none) is best, at each threshold. With
#' \eqn{NB_s} the net benefit of strategy \eqn{s} and expectations over its
#' uncertainty:
#' \deqn{ENB_{perfect} = E\{\max_s NB_s\}, \qquad
#'       ENB_{current} = \max_s E\{NB_s\},}
#' \deqn{EVPI = ENB_{perfect} - ENB_{current}.}
#' EVPI is in net benefit units (net true positives per patient); it is 0 when
#' the best strategy is certain, and multiplied by the number of patients the
#' decision concerns it gives the population EVPI.
#'
#' The **relative EVPI** compares the gain from perfect information with the
#' gain from using the model(s), both relative to the best default strategy:
#' \deqn{rEVPI = \frac{ENB_{perfect} - \max\{0, E(NB_{all})\}}
#'                   {ENB_{current} - \max\{0, E(NB_{all})\}}.}
#' For example, 1.25 means that perfect information would give 25% more net
#' benefit than using the best model now. It is `Inf` when no model is better
#' than the default strategies in expectation (the denominator is 0) and EVPI
#' is positive, and `NA` when both are 0.
#'
#' @section Uncertainty:
#' The uncertainty about the net benefit is represented by `n_sim` draws, in
#' each of which the net benefit of every model, of treating all and of
#' treating none is recomputed with the same estimator as [nb()] (any outcome
#' type, any number of models). The expected net benefits are the means of the
#' draws (not the point estimates), as in Sadatsafavi et al. (2023), so EVPI
#' is never negative.
#' * `method = "bayesian"` (default): Bayesian bootstrap. The patients'
#'   weights are Dirichlet(1, ..., 1), which gives draws from the posterior
#'   distribution of the net benefit under a non-informative prior (Rubin's
#'   Bayesian bootstrap).
#' * `method = "ordinary"`: ordinary bootstrap. The weights are the numbers of
#'   times each patient is drawn when sampling `n` patients with replacement.
#' * Results of `nb(bayesian = TRUE)` (Bayesian net benefit, binary outcome):
#'   the draws are from the posterior of the prevalence and of each model's
#'   sensitivity and specificity (exact Beta / Dirichlet draws, as in
#'   bayesDCA), and `method` is not used. Several models are drawn jointly (a
#'   Dirichlet over the patterns of treat / do not treat across the models),
#'   so their correlation on the same patients is kept, and each model keeps
#'   its own posterior.
#'
#' With case weights in [nb()] (`weights`), the bootstrap draws are taken from
#' the weighted sample. The predicted risks are kept fixed in every draw. With
#' `censoring = "cipcw"`, the censoring model is refitted in every draw and the
#' censoring weights recomputed from it (see [cipcw_control()]).
#'
#' `evpi()`, [p_useful()] and [p_best()] called with the same seed (and the
#' same results, `n_sim` and `method`) use the same draws.
#'
#' @param x A result of [nb()].
#' @param ... More results of [nb()], for other models evaluated on the same
#'   data, outcome, thresholds (and horizon, cause and censoring). A named
#'   result is labelled with its name.
#' @param n_sim Number of draws. Default 1000.
#' @param method `"bayesian"` (Bayesian bootstrap, default), `"ordinary"`
#'   (ordinary bootstrap) or `"asymptotic"`: see "Uncertainty". `"asymptotic"`
#'   is for a binary outcome and one model only: the net benefits of the model
#'   and of treating all are taken as bivariate normal, with the moments of
#'   Sadatsafavi et al. (2023), and \eqn{E\{\max(0, NB_{model}, NB_{all})\}}
#'   is computed in closed form (it needs the mvtnorm package). Not used (an
#'   error if given) for results of `nb(bayesian = TRUE)`.
#' @param labels Optional names of the models, one per result; default: the
#'   results' labels.
#' @return A data frame of class `voipred_evpi` with one row per threshold
#'   and columns `threshold`, `ENB_current`, `ENB_perfect`, `EVPI` and
#'   `rEVPI`. The expected net benefit of each strategy is in the attribute
#'   `ENB` (a threshold x strategy matrix).
#' @references Sadatsafavi M, Lee TY, Wynants L, Vickers AJ, Gustafson P.
#'   Value-of-information analysis for external validation of risk prediction
#'   models. Medical Decision Making 2023; 43(5): 564-575.
#' @seealso [evsi()], [p_useful()], [p_best()], [plot_evpi()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, "day30", thresholds = c(0.02, 0.05, 0.1, 0.2))
#' set.seed(1)
#' evpi(res, n_sim = 500)
#' evpi(res, n_sim = 500, method = "ordinary")
#' evpi(res, method = "asymptotic")
#' @export
evpi <- function(x, ..., n_sim = 1000, method = c("bayesian", "ordinary", "asymptotic"),
                 labels = NULL) {
  method_given <- !missing(method)
  method <- match.arg(method)
  res <- collect_voi(x, list(...), labels, "evpi()", "`n_sim`, `method` and `labels`")
  if (method == "asymptotic") {
    if (is_bayes(res)) check_bayes_method(TRUE)
    return(evpi_asymptotic(res))
  }
  evpi_from_draws(nb_draws(res, n_sim, method, method_given))
}

#' Probability that a model is useful
#'
#' For each model and threshold, the probability that using the model gives a
#' higher net benefit than both default strategies (treat all and treat
#' none):
#' \deqn{P(useful) = P\{NB_{model} > \max(NB_{all}, 0)\} = P(\Delta NB > 0),}
#' estimated as the share of the draws in which the model's net benefit is
#' strictly higher. A model that ties with a default strategy (e.g. at a
#' threshold where it treats everyone, or no one) is not counted as useful.
#'
#' @inheritSection evpi Uncertainty
#' @inheritParams evpi
#' @param method `"bayesian"` (Bayesian bootstrap, default) or `"ordinary"`
#'   (ordinary bootstrap): see "Uncertainty". Not used (an error if given) for
#'   results of `nb(bayesian = TRUE)`.
#' @return A data frame of class `voipred_p_useful` with columns `threshold`,
#'   `model` and `p_useful`, one row per threshold and model.
#' @seealso [p_best()], [evpi()], [plot_p_useful()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, "day30", thresholds = c(0.02, 0.05, 0.1, 0.2))
#' set.seed(1)
#' p_useful(res, n_sim = 500)
#' @export
p_useful <- function(x, ..., n_sim = 1000, method = c("bayesian", "ordinary"),
                     labels = NULL) {
  method_given <- !missing(method)
  method <- match.arg(method)
  res <- collect_voi(x, list(...), labels, "p_useful()", "`n_sim`, `method` and `labels`")
  p_useful_from_draws(nb_draws(res, n_sim, method, method_given))
}

#' Probability that each strategy is the best
#'
#' For each strategy (each model, treat all and treat none) and threshold,
#' the probability that it has the highest net benefit:
#' \deqn{P(best_s) = P\{NB_s = \max_t NB_t\},}
#' estimated as the share of the draws in which it is the best. Exact ties
#' are shared equally, so the probabilities of all strategies add up to 1 at
#' each threshold.
#'
#' @inheritSection evpi Uncertainty
#' @inheritParams p_useful
#' @return A data frame of class `voipred_p_best` with columns `threshold`,
#'   `strategy` and `p_best`, one row per threshold and strategy.
#' @seealso [p_useful()], [evpi()], [plot_p_best()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' fit2 <- glm(day30 ~ age, family = binomial, data = dev)
#' z <- c(0.02, 0.05, 0.1, 0.2)
#' res1 <- nb(fit1, val, "day30", thresholds = z, label = "Full")
#' res2 <- nb(fit2, val, "day30", thresholds = z, label = "Age only")
#' set.seed(1)
#' p_best(res1, res2, n_sim = 500)
#' @export
p_best <- function(x, ..., n_sim = 1000, method = c("bayesian", "ordinary"),
                   labels = NULL) {
  method_given <- !missing(method)
  method <- match.arg(method)
  res <- collect_voi(x, list(...), labels, "p_best()", "`n_sim`, `method` and `labels`")
  p_best_from_draws(nb_draws(res, n_sim, method, method_given))
}

# ---- draws (internal) ---------------------------------------------------------

# Draws of the net benefit of every strategy for the results `res` (a named
# list from collect_voi()): bootstrap draws ("bayesian" / "ordinary") or, for
# results of nb(bayesian = TRUE), draws from the posterior. A list with `nb`
# (draw x threshold x strategy), `thresholds`, `strategies`, `models`,
# `method` and `n_sim`.
nb_draws <- function(res, n_sim, method = "bayesian", method_given = FALSE) {
  check_n_sim(n_sim)
  if (is_bayes(res)) {
    check_bayes_method(method_given)
    return(finish_draws(res, bayes_draws(res, n_sim), "posterior"))
  }
  v <- attr(res[[1]], "voi")
  n <- length(v$risk)
  W <- vapply(seq_len(n_sim), function(i) draw_weights_one(n, method, v$weights),
              numeric(n))
  if (n == 1) W <- matrix(W, nrow = 1)
  draws_from(res, W, method)
}

# The draws for the results `res` and given patient weights W (n x draws).
draws_from <- function(res, W, method) finish_draws(res, nb_eval(res, W), method)

finish_draws <- function(res, nbs, method) {
  v <- attr(res[[1]], "voi")
  bad <- apply(nbs, 2, anyNA)
  if (any(bad))
    warning("The net benefit could not be estimated in some draws at ",
            if (sum(bad) == 1) "threshold " else "thresholds ",
            format_thresholds(v$thresholds[bad]),
            " (no treated patient followed up to the horizon); the value of ",
            "information is NA there.", call. = FALSE)
  list(nb = nbs, thresholds = v$thresholds,
       strategies = c(names(res), "Treat All", "Treat None"), models = names(res),
       method = method, n_sim = dim(nbs)[1])
}

evpi_from_draws <- function(d) {
  nbs <- d$nb
  S <- dim(nbs)[3]
  enb <- apply(nbs, c(2, 3), mean)
  if (!is.matrix(enb)) enb <- matrix(enb, nrow = dim(nbs)[2])
  colnames(enb) <- d$strategies
  perfect <- colMeans(matrix(apply(nbs, c(1, 2), max), nrow = dim(nbs)[1]))
  current <- apply(enb, 1, max)
  default <- pmax(enb[, S - 1], 0)
  out <- data.frame(threshold = d$thresholds, ENB_current = current,
                    ENB_perfect = perfect, EVPI = perfect - current,
                    rEVPI = relative_evpi(perfect, current, default))
  attr(out, "ENB") <- cbind(threshold = d$thresholds, enb)
  attr(out, "method") <- d$method
  attr(out, "n_sim") <- d$n_sim
  class(out) <- c("voipred_evpi", "data.frame")
  out
}

p_useful_from_draws <- function(d) {
  nbs <- d$nb
  S <- dim(nbs)[3]
  default <- pmax(nbs[, , S - 1, drop = FALSE], 0)[, , 1]
  out <- do.call(rbind, lapply(seq_along(d$models), function(m) {
    better <- nbs[, , m] > default
    if (!is.matrix(better)) better <- matrix(better, nrow = dim(nbs)[1])
    data.frame(threshold = d$thresholds, model = d$models[m],
               p_useful = colMeans(better))
  }))
  rownames(out) <- NULL
  attr(out, "method") <- d$method
  attr(out, "n_sim") <- d$n_sim
  class(out) <- c("voipred_p_useful", "data.frame")
  out
}

p_best_from_draws <- function(d) {
  nbs <- d$nb
  dims <- dim(nbs)
  best <- apply(nbs, c(1, 2), max)                      # draw x threshold
  share <- array(0, dims)
  for (s in seq_len(dims[3])) share[, , s] <- nbs[, , s] >= best - 1e-12
  ties <- apply(share, c(1, 2), sum)
  for (s in seq_len(dims[3])) share[, , s] <- share[, , s] / ties
  out <- do.call(rbind, lapply(seq_len(dims[3]), function(s)
    data.frame(threshold = d$thresholds, strategy = d$strategies[s],
               p_best = colMeans(matrix(share[, , s], nrow = dims[1])))))
  rownames(out) <- NULL
  attr(out, "method") <- d$method
  attr(out, "n_sim") <- d$n_sim
  class(out) <- c("voipred_p_best", "data.frame")
  out
}

#' Expected value of sample information (EVSI) for model validation
#'
#' The expected gain in net benefit from a future validation study of size
#' \eqn{N} (`future_n`), added to the current sample, before deciding which
#' strategy to use. EVSI is between 0 and EVPI, and approaches EVPI as
#' \eqn{N} grows: comparing it with the cost of the study helps decide whether
#' more validation data are worth collecting.
#'
#' It is computed with the bootstrap-based algorithm of Sadatsafavi et al.
#' (2025). For each of `n_sim` draws:
#' 1. draw weights \eqn{W} for the current sample (`method = "bayesian"`:
#'    Bayesian bootstrap, Dirichlet; `"ordinary"`: ordinary bootstrap), which
#'    defines a "true" population,
#'    and compute the true net benefit of every strategy;
#' 2. for each \eqn{N}: draw a future sample of \eqn{N} patients from that
#'    population (multinomial counts with probabilities \eqn{W}), pool it with
#'    the current sample, compute the net benefit of every strategy in the
#'    pooled data, pick the best strategy, and record its *true* net benefit
#'    (if strategies tie exactly in the pooled data, the mean of their true
#'    net benefits, i.e. a random choice among them).
#'
#' For Bayesian results (`nb(bayesian = TRUE)`, binary outcome) the conjugate
#' algorithm of Sadatsafavi et al. (2025) is used instead: each draw from the
#' posterior (prevalence, sensitivity, specificity) is a "true" population;
#' a future sample of \eqn{N} patients is simulated from it (events, and the
#' patients treated among them, binomial / multinomial), the posterior means
#' are updated with it, the strategy with the highest posterior mean net
#' benefit is chosen, and its true net benefit recorded.
#'
#' Then \eqn{ENB_{sample}} is the mean of the recorded true net benefits and
#' \eqn{EVSI = ENB_{sample} - ENB_{current}}, with \eqn{ENB_{current}} as in
#' [evpi()] (from the same draws). The net benefit is computed with the same
#' estimator as [nb()], so this works for binary, survival and
#' competing-risks outcomes (a future patient is a resampled current patient,
#' with the same follow-up), with any number of models. With
#' `censoring = "cipcw"` the censoring model is refitted in every draw and
#' every pooled sample (see [cipcw_control()]).
#'
#' @inheritParams p_useful
#' @param future_n Sizes of the future validation sample: one or more
#'   positive integers.
#' @param method `"bayesian"` (Bayesian bootstrap, default) or `"ordinary"`
#'   (ordinary bootstrap): how the uncertainty about the current sample is
#'   represented. Not used (an error if given) for results of
#'   `nb(bayesian = TRUE)`.
#' @return A data frame of class `voipred_evsi` with columns `threshold`,
#'   `future_n`, `EVSI` and `EVPI` (from the same draws), one row per
#'   threshold and future sample size.
#' @references Sadatsafavi M, et al. Expected value of sample information
#'   calculations for risk prediction model validation. Medical Decision
#'   Making 2025; 45(3): 232-244.
#' @seealso [evpi()], [plot_evsi()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, "day30", thresholds = c(0.02, 0.05, 0.1))
#' set.seed(1)
#' evsi(res, future_n = c(250, 500, 1000, 4000), n_sim = 500)
#' @export
evsi <- function(x, ..., future_n, n_sim = 1000, method = c("bayesian", "ordinary"),
                 labels = NULL) {
  method_given <- !missing(method)
  method <- match.arg(method)
  if (missing(future_n) || !is.numeric(future_n) || !length(future_n) ||
      anyNA(future_n) || any(future_n < 1) || any(future_n != round(future_n)))
    stop("`future_n` must be one or more positive whole numbers: the sizes of the ",
         "future validation sample.", call. = FALSE)
  check_n_sim(n_sim)
  res <- collect_voi(x, list(...), labels, "evsi()",
                     "`future_n`, `n_sim`, `method` and `labels`")
  v <- attr(res[[1]], "voi")
  if (is_bayes(res)) {
    check_bayes_method(method_given)
    r <- bayes_evsi(res, future_n, n_sim)
    Z <- length(v$thresholds)
    out <- data.frame(threshold = rep(v$thresholds, length(future_n)),
                      future_n = rep(future_n, each = Z),
                      EVSI = as.vector(r$EVSI), EVPI = rep(r$EVPI, length(future_n)))
    attr(out, "method") <- "posterior"
    attr(out, "n_sim") <- n_sim
    class(out) <- c("voipred_evsi", "data.frame")
    return(out)
  }
  n <- length(v$risk)
  w0 <- if (is.null(v$weights)) rep(1, n) else v$weights
  Z <- length(v$thresholds)
  J <- length(future_n)
  S <- length(res) + 2
  true_sum <- matrix(0, Z, S)                 # sum over draws of the true NBs
  perfect_sum <- numeric(Z)
  sample_sum <- matrix(0, Z, J)
  na_z <- logical(Z)
  chunk <- 200
  done <- 0
  while (done < n_sim) {
    C <- min(chunk, n_sim - done)
    W <- matrix(0, n, C)
    Wf <- replicate(J, matrix(0, n, C), simplify = FALSE)
    for (b in seq_len(C)) {                    # same random-number order as the reference code
      W[, b] <- evsi_weights_one(n, method, v$weights)
      for (j in seq_len(J))
        Wf[[j]][, b] <- stats::rmultinom(1, future_n[j], W[, b])
    }
    # the true net benefits; with cIPCW the weights are put on the scale of
    # patients (counts, or mean 1), as the censoring-model refit sees them
    truth <- nb_eval(res, if (identical(v$censoring, "cipcw")) W * n else W)  # C x Z x S
    na_z <- na_z | apply(truth, 2, anyNA)
    true_sum <- true_sum + apply(truth, c(2, 3), sum)
    perfect_sum <- perfect_sum + colSums(matrix(apply(truth, c(1, 2), max), nrow = C))
    for (j in seq_len(J)) {
      pooled <- nb_eval(res, Wf[[j]] + w0)
      got <- winner_value(pooled, truth)
      na_z <- na_z | apply(got, 2, anyNA)
      sample_sum[, j] <- sample_sum[, j] + colSums(got)
    }
    done <- done + C
  }
  current <- apply(true_sum / n_sim, 1, max)
  evpi_z <- perfect_sum / n_sim - current
  evsi_zj <- sample_sum / n_sim - current
  evpi_z[na_z] <- NA
  evsi_zj[na_z, ] <- NA
  if (any(na_z))
    warning("The net benefit could not be estimated in some draws at ",
            if (sum(na_z) == 1) "threshold " else "thresholds ",
            format_thresholds(v$thresholds[na_z]),
            " (no treated patient followed up to the horizon); EVSI is NA there.",
            call. = FALSE)
  out <- data.frame(threshold = rep(v$thresholds, J), future_n = rep(future_n, each = Z),
                    EVSI = as.vector(evsi_zj), EVPI = rep(evpi_z, J))
  attr(out, "method") <- method
  attr(out, "n_sim") <- n_sim
  class(out) <- c("voipred_evsi", "data.frame")
  out
}

# ---- internal ----------------------------------------------------------------

# Results of nb() that can be analysed together, named by their labels.
collect_voi <- function(x, more, labels, fun, args) {
  res <- collect_nb(x, more, labels, fun = fun, args = args)
  if (any(vapply(res, function(r) is.null(attr(r, "voi")), logical(1))))
    stop("These results of nb() were made by an older version of voipred; run ",
         "nb() again.", call. = FALSE)
  res
}

# Results of nb(bayesian = TRUE): their uncertainty is their posterior.
is_bayes <- function(res) isTRUE(attr(res[[1]], "bayesian"))

check_bayes_method <- function(method_given) {
  if (method_given)
    stop("These results of nb() are Bayesian (`bayesian = TRUE`): the value of ",
         "information uses draws from their posterior, so `method` is not used.",
         call. = FALSE)
}

check_n_sim <- function(n_sim) {
  if (!is.numeric(n_sim) || length(n_sim) != 1 || is.na(n_sim) || n_sim < 2 ||
      n_sim != round(n_sim))
    stop("`n_sim` must be a whole number of at least 2.", call. = FALSE)
}

# One draw of patient weights, scaled to sum to n. Bayesian bootstrap:
# Dirichlet(1, ..., 1) (times the case weights); ordinary bootstrap:
# multinomial counts (with probabilities proportional to the case weights).
draw_weights_one <- function(n, method, w0 = NULL) {
  if (method == "bayesian") {
    u <- stats::rgamma(n, 1, 1)
    if (!is.null(w0)) u <- u * w0
    u * n / sum(u)
  } else {
    as.vector(stats::rmultinom(1, n, if (is.null(w0)) rep(1 / n, n) else w0))
  }
}

# The weights of the current sample in evsi(), drawn as in the reference code
# of Sadatsafavi et al. (2025): Bayesian bootstrap rgamma(n, 1, 1) normalized;
# ordinary bootstrap rmultinom(1, n, 1/n) / n.
evsi_weights_one <- function(n, method, w0 = NULL) {
  if (method == "bayesian") {
    u <- stats::rgamma(n, 1, 1)
    if (!is.null(w0)) u <- u * w0
    u / sum(u)
  } else {
    as.vector(stats::rmultinom(1, n, if (is.null(w0)) rep(1 / n, n) else w0)) / n
  }
}

# Net benefit of every strategy for each column of weights W (n x B):
# an array draw x threshold x strategy (models, treat all, treat none).
nb_eval <- function(res, W) {
  v1 <- attr(res[[1]], "voi")
  z <- v1$thresholds
  B <- ncol(W)
  Z <- length(z)
  M <- length(res)
  out <- array(NA_real_, c(B, Z, M + 2))
  out[, , M + 2] <- 0
  if (v1$type == "binary") {
    odds <- z / (1 - z)
    gain <- outer(v1$y, rep(1, Z)) - outer(1 - v1$y, odds)      # Y - (1 - Y) z/(1-z)
    model_gain <- lapply(res, function(r) {
      p <- attr(r, "voi")$risk
      gain * (outer(p, z, ">=") & p > 0)
    })
    sw <- colSums(W)
    out[, , M + 1] <- crossprod(W, gain) / sw
    for (m in seq_len(M)) out[, , m] <- crossprod(W, model_gain[[m]]) / sw
  } else {
    refit <- lapply(res, function(r) censoring_refitter(attr(r, "voi"), W))
    # models validated with the same censoring settings share the refit
    share <- vapply(seq_len(M), function(m)
      m > 1 && !is.null(refit[[m]]) && identical(attr(res[[m]], "voi")$cipcw, v1$cipcw),
      logical(1))
    for (b in seq_len(B)) {
      cw1 <- NULL
      for (m in seq_len(M)) {
        cw <- if (is.null(refit[[m]])) NULL else if (share[m]) cw1 else refit[[m]](W[, b])
        if (m == 1) cw1 <- cw
        r <- nb_one(attr(res[[m]], "voi"), W[, b], cw)
        out[b, , m] <- r$NB_model
        if (m == 1) out[b, , M + 1] <- r$NB_all
      }
    }
  }
  out
}

# The true net benefit of the strategy that is best in the pooled data, for
# each draw and threshold (C x Z). Exact ties in the pooled data (possible
# with integer weights) are broken at random, i.e. the true net benefits of
# the tied strategies are averaged.
winner_value <- function(pooled, truth) {
  best <- apply(pooled, c(1, 2), max)
  num <- den <- array(0, dim(best))
  for (s in seq_len(dim(pooled)[3])) {
    tied <- pooled[, , s] >= best - 1e-12
    num <- num + ifelse(tied, truth[, , s], 0)
    den <- den + tied
  }
  num / den
}

nb_one <- function(v, w, cw = NULL) {
  if (v$censoring == "km")
    return(suppressWarnings(nb_core_km(v$time, v$status, v$k, v$risk, v$thresholds,
                                       v$horizon, w)))
  nb_core_cipcw(v$time, v$status, v$k, v$risk, v$thresholds, v$horizon,
                cw$weight_event, cw$weight_horizon, w)
}

# ---- refitting the cIPCW censoring model in each draw ------------------------

# Censoring models whose case weights voipred knows, and the name of the
# fitting function's case-weights argument.
known_weight_args <- c(coxph = "weights", survreg = "weights", cph = "weights",
                       psm = "weights", flexsurvreg = "weights",
                       prodlim = "caseweights")

# For a cIPCW result, a function(w) returning the censoring weights from the
# censoring model refitted with the draw's weights w (NULL for Kaplan-Meier).
# Checks first, for all draws W, that the refit is possible.
censoring_refitter <- function(v, W) {
  if (!identical(v$censoring, "cipcw")) return(NULL)
  ctrl <- v$cipcw
  if (is.null(ctrl$model))
    stop("The cIPCW censoring weights were given directly (cipcw_control(",
         "weight_event =, weight_horizon =)), so they cannot be re-estimated in ",
         "each draw, and keeping them fixed would understate the uncertainty. For ",
         "the value of information, give the censoring model instead: ",
         "cipcw_control(model = ).", call. = FALSE)
  model <- ctrl$model
  cls <- intersect(class(model), names(known_weight_args))
  fit <- if (!is.null(ctrl$refit_function)) {
    ctrl$refit_function
  } else if (length(cls)) {
    function(data, weights) refit_call(model, data, known_weight_args[[cls[1]]], weights)
  } else if (all(apply(W, 2, whole_counts_ok))) {
    function(data, weights)
      refit_call(model, data[rep(seq_len(nrow(data)), whole_counts(weights)), , drop = FALSE])
  } else {
    stop("The censoring model (", class(model)[1], ") must be refitted in each draw ",
         "with the draw's weights, which are not whole numbers (Bayesian bootstrap), ",
         "and voipred does not know how to fit a ", class(model)[1], " model with ",
         "case weights. Give `refit_function = function(data, weights)` in ",
         "cipcw_control(), returning the censoring model fitted with these case ",
         "weights, or use method = \"ordinary\".", call. = FALSE)
  }
  function(w) {
    keep <- w > 0
    d <- v$data[keep, , drop = FALSE]
    m <- tryCatch(fit(d, w[keep]), error = function(e)
      stop("Could not refit the censoring model in a draw (", conditionMessage(e),
           "). Give `refit_function = function(data, weights)` in cipcw_control(), ",
           "returning the censoring model fitted with these case weights.",
           call. = FALSE))
    ctrl$model <- m
    suppressWarnings(cipcw_weights(ctrl, v$data, v$info, v$horizon, v$k))
  }
}

# The model's own call evaluated again with new data (and case weights under
# the argument `wname`), in the environment of the model's formula.
refit_call <- function(model, data, wname = NULL, weights = NULL) {
  cl <- tryCatch(stats::getCall(model), error = function(e) NULL)
  if (is.null(cl)) stop("the model keeps no call to refit it with", call. = FALSE)
  f <- model_formula(model)
  if (is.null(f) && is.list(model)) {        # e.g. flexsurvreg keeps it in all.formulae
    g <- model$all.formulae[[1]]
    if (inherits(g, "formula") && length(g) == 3) f <- g
  }
  env <- if (!is.null(f)) environment(f)
  env <- new.env(parent = if (is.environment(env)) env else globalenv())
  # the formula itself, so a formula kept in a variable needs no lookup
  if (!is.null(f)) {
    nms <- names(cl)
    if ("formula" %in% nms) cl$formula <- f
    else if (length(cl) > 1 && (is.null(nms) || !nzchar(nms[2]))) cl[[2]] <- f
  }
  env$.voipred_data <- data
  cl$data <- quote(.voipred_data)
  if (!is.null(wname)) cl[[wname]] <- weights   # the values, so no name is looked up
  eval(cl, env)
}

# Whole-number weights (up to a common scale, e.g. counts / n): the positive
# weights divided by the smallest one are whole numbers.
whole_counts_ok <- function(w) {
  x <- w[w > 0] / min(w[w > 0])
  all(abs(x - round(x)) < 1e-8 * pmax(1, x))
}

whole_counts <- function(w) {
  if (all(abs(w - round(w)) < 1e-8 * pmax(1, w))) round(w) else round(w / min(w))
}

relative_evpi <- function(perfect, current, default) {
  r <- (perfect - default) / (current - default)
  r[is.nan(r)] <- NA
  r
}

# ---- asymptotic EVPI (binary outcome, one model) -----------------------------

evpi_asymptotic <- function(res) {
  if (length(res) != 1)
    stop("`method = \"asymptotic\"` works for one model only; use `method = ",
         "\"bayesian\"` or `\"ordinary\"` for several models.", call. = FALSE)
  v <- attr(res[[1]], "voi")
  if (v$type != "binary")
    stop("`method = \"asymptotic\"` works for a binary outcome only; use `method = ",
         "\"bayesian\"` or `\"ordinary\"` for survival and competing-risks outcomes.",
         call. = FALSE)
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("The mvtnorm package is needed for `method = \"asymptotic\"`.", call. = FALSE)
  z <- v$thresholds
  nb0 <- nb_core(v$y, v$risk, z, v$weights)
  perfect <- numeric(length(z))
  for (j in seq_along(z)) {
    m <- nb_moments(v$y, v$risk, z[j], v$weights)
    current_j <- max(0, nb0$NB_model[j], nb0$NB_all[j])
    if (!is.finite(m[["rho"]])) { perfect[j] <- current_j; next }
    rho <- min(max(m[["rho"]], -0.999999), 0.999999)
    perfect[j] <- tryCatch(
      mu_max_trunc_bvn(m[["mu1"]], m[["mu2"]], m[["sigma1"]], m[["sigma2"]], rho),
      error = function(e) NA_real_)
  }
  current <- pmax(0, nb0$NB_model, nb0$NB_all)
  default <- pmax(nb0$NB_all, 0)
  out <- data.frame(threshold = z, ENB_current = current, ENB_perfect = perfect,
                    EVPI = perfect - current,
                    rEVPI = relative_evpi(perfect, current, default))
  enb <- cbind(threshold = z, nb0$NB_model, nb0$NB_all, 0)
  colnames(enb) <- c("threshold", names(res), "Treat All", "Treat None")
  attr(out, "ENB") <- enb
  attr(out, "method") <- "asymptotic"
  class(out) <- c("voipred_evpi", "data.frame")
  out
}

# Means, SDs and correlation of the estimated NB of the model and of treating
# all (Sadatsafavi et al. 2023, as in predtools::calc_NB_moments).
nb_moments <- function(y, p, z, w = NULL) {
  n <- length(y)
  a <- as.numeric(p >= z & p > 0)
  if (is.null(w)) w <- rep(1, n)
  rho <- stats::weighted.mean(y, w)
  tpr <- stats::weighted.mean(y * a, w) / rho
  fpr <- stats::weighted.mean(a * (1 - y), w) / (1 - rho)
  tz <- z / (1 - z)
  mu <- c(rho * tpr - tz * (1 - rho) * fpr, rho - tz * (1 - rho))
  s_y <- rho * (1 - rho) / n
  s_tpr <- rho * tpr * (1 - rho * tpr) / n
  s_fpr <- (1 - rho) * fpr * (1 - (1 - rho) * fpr) / n
  s11 <- s_tpr + tz^2 * s_fpr + 2 * tz * rho * (1 - rho) * tpr * fpr / n
  s22 <- (1 + tz)^2 * s_y
  s12 <- s_y * (1 + tz) * (tpr + tz * fpr)
  c(mu1 = mu[1], mu2 = mu[2], sigma1 = sqrt(s11), sigma2 = sqrt(s22),
    rho = s12 / (sqrt(s11) * sqrt(s22)))
}

# E{max(0, X1, X2)} for (X1, X2) bivariate normal (Lee et al.; as in
# predtools::mu_max_trunc_bvn).
mu_max_trunc_bvn <- function(mu1, mu2, sigma1, sigma2, rho,
                             precision = .Machine$double.eps) {
  ind <- function(cond) as.numeric(cond)
  uv <- function(mu1, mu2, sigma1, sigma2, rho) {
    tmp1 <- sigma1 - rho * sigma2
    tmp2 <- (rho * sigma2 * mu1 - sigma1 * mu2) / (sigma1 * sigma2 * sqrt(1 - rho^2))
    mu1 * (ind(tmp1 > precision) + ind(abs(tmp1) < precision) * stats::pnorm(tmp2)) -
      stats::pnorm(tmp2) * (-sigma1 * stats::dnorm(-mu1 / sigma1) +
                              mu1 * stats::pnorm(-mu1 / sigma1))
  }
  vdu <- function(mu1, mu2, sigma1, sigma2, rho) {
    tmp1 <- sigma1 - rho * sigma2
    if (abs(tmp1) < precision) return(0)
    delta <- sign(tmp1)
    alpha <- sigma1 * mu2 - rho * sigma2 * mu1
    beta <- sigma1 * sigma2 * sqrt(1 - rho^2)
    denom <- sqrt(sigma1^2 + sigma2^2 - 2 * rho * sigma1 * sigma2)
    s_rho <- -sigma2 * sqrt(1 - rho^2) / denom
    t1 <- delta * mu1 * (stats::pnorm(delta * (mu2 - mu1) / denom) -
      as.numeric(mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                                  upper = c(delta * (mu2 - mu1) / denom, -delta * alpha / beta),
                                  mean = c(0, 0),
                                  sigma = matrix(c(1, s_rho, s_rho, 1), 2))[[1]]))
    a <- (mu2 - mu1) / tmp1
    b <- (beta / abs(tmp1)) / sigma1
    t <- sqrt(1 + b^2)
    t2 <- -delta * sigma1 / t * stats::dnorm(a / t) *
      (1 - stats::pnorm(-delta * t * alpha / beta + a * b / t))
    t1 + t2
  }
  uv(mu1, mu2, sigma1, sigma2, rho) + uv(mu2, mu1, sigma2, sigma1, rho) -
    vdu(mu1, mu2, sigma1, sigma2, rho) - vdu(mu2, mu1, sigma2, sigma1, rho)
}
