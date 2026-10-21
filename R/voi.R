# Value of information for model validation: P(useful), P(best), EVPI,
# relative EVPI and EVSI (Sadatsafavi et al. 2023, 2025).
#
# The uncertainty about the net benefit is its posterior (see nb_val()): for a
# binary outcome without case weights, draws from the conjugate Beta /
# Dirichlet posterior of the prevalence, sensitivity and specificity (with
# nb_val()'s prior); otherwise (survival, competing risks, case weights; improper
# prior) the Bayesian bootstrap: each draw re-weights the patients with
# Dirichlet(1, ..., 1) weights and the net benefit of every strategy is
# recomputed with the same estimator as nb_val() (Kaplan-Meier / Aalen-Johansen,
# or cIPCW with the censoring model refitted in each draw).

#' Expected value of perfect information (EVPI) for model validation
#'
#' The expected loss in net benefit due to uncertainty about which strategy
#' (each model, treat all, treat none) is best, at each threshold. With
#' \eqn{NB_s} the net benefit of strategy \eqn{s} and expectations over its
#' posterior:
#' \deqn{ENB_{perfect} = E\{\max_s NB_s\}, \qquad
#'       ENB_{current} = \max_s E\{NB_s\},}
#' \deqn{EVPI = ENB_{perfect} - ENB_{current}.}
#' EVPI is in net benefit units (net true positives per patient); it is 0 when
#' the best strategy is certain. It is never negative: both terms are means
#' over the same draws.
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
#' The uncertainty about the net benefit is its posterior (see the prior of
#' [nb_val()]), represented by `n_sim` draws, in each of which the net benefit of
#' every model, of treating all and of treating none is computed:
#' * binary outcome: draws from the conjugate posterior of the prevalence and
#'   of each model's sensitivity and specificity, with the prior of [nb_val()]
#'   (exact Beta / Dirichlet draws, as in bayesDCA). Several models are drawn
#'   jointly (a Dirichlet over the patterns of treat / do not treat across the
#'   models), so their correlation on the same patients is kept, and each model
#'   keeps its own posterior;
#' * survival and competing-risks outcomes, and binary outcomes with case
#'   weights (improper prior): the Bayesian bootstrap. Each
#'   draw gives the patients Dirichlet(1, ..., 1) weights (times the case
#'   weights) and recomputes the net benefit with the same estimator as
#'   [nb_val()]. With cIPCW (`cipcw`), the censoring model is refitted in
#'   every draw and the censoring weights recomputed from it (see
#'   [cipcw_control()]).
#'
#' The predicted risks are kept fixed in every draw. The expected net benefits
#' are the means of the draws. `evpi()` also computes P(useful) and P(best)
#' from its draws, so the three come from the same draws: read them with
#' `p_useful(e)` and `p_best(e)` for a result `e` of `evpi()`. [p_useful()]
#' and [p_best()] given results of [nb_val()] instead draw again; called
#' with the same seed (and the same models and `n_sim`) as `evpi()`, they
#' use the same draws.
#'
#' @param x A result of [nb_val()], for one model or several.
#' @param ... More results of [nb_val()], for other models evaluated on the same
#'   data, outcome, thresholds and prior (and horizon, cause and censoring),
#'   combined as by [c()][c.voipred_nb]. A result with one model given here
#'   with a name is labelled with its name.
#' @param n_sim Number of draws. Default 1000.
#' @param labels Optional names of the models, one per model, in order;
#'   default: the models' labels.
#' @return An object of class `voipred_evpi`, a list with `estimates` (a data
#'   frame with one row per threshold and columns `threshold`, `ENB_current`,
#'   `ENB_perfect`, `EVPI` and `rEVPI`), `ENB` (the expected net benefit of
#'   each strategy at each threshold), `p_useful` and `p_best` (the results of
#'   [p_useful()] and [p_best()] from the same draws, also returned by
#'   `p_useful(e)` and `p_best(e)`), `models`, `outcome`, `censoring`,
#'   `prior`, `draws` (`"posterior"` or `"Bayesian bootstrap"`), `n_sim` and
#'   `call`. It has `print()`, `summary()` (which adds the probability that
#'   the strategy that is best now is the best), `plot()` (EVPI, relative
#'   EVPI, P(useful) or P(best)) and `as.data.frame()` methods; EVPI can also
#'   be drawn under the decision curve with `plot(<nb_val result>, evpi = )`.
#' @references Sadatsafavi M, Lee TY, Wynants L, Vickers AJ, Gustafson P.
#'   Value-of-information analysis for external validation of risk prediction
#'   models. Medical Decision Making 2023; 43(5): 564-575. Rubin DB. The
#'   Bayesian bootstrap. Annals of Statistics 1981; 9(1): 130-134.
#' @seealso [evsi()], [p_useful()], [p_best()], [plot.voipred_evpi()], [summary.voipred_voi()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse,
#'             family = binomial, data = df_dev_binary)
#' df_val_binary$p1 <- predict(fit1, newdata = df_val_binary, type = "response")
#' fit2 <- glm(day30 ~ age, family = binomial, data = df_dev_binary)
#' df_val_binary$p2 <- predict(fit2, newdata = df_val_binary, type = "response")
#' z <- seq(0.01, 0.3, by = 0.01)
#' res1 <- nb_val(data = df_val_binary,
#'                predicted_risk = df_val_binary$p1,
#'                outcome = "day30",
#'                thresholds = z,
#'                label = "Full")
#' res2 <- nb_val(data = df_val_binary,
#'                predicted_risk = df_val_binary$p2,
#'                outcome = "day30",
#'                thresholds = z,
#'                label = "Age only")
#' set.seed(1)
#' e <- evpi(res1, res2)
#' e
#' summary(e)
#' plot(e)
#' plot(e, type = "relative")
#' # P(useful) and P(best), from the same draws
#' p_useful(e)
#' p_best(e)
#' plot(e, type = "p_best")
#' @export
evpi <- function(x, ..., n_sim = 1000, labels = NULL) {
  call <- match.call()
  res <- collect_voi(x, list(...), labels, "evpi()", "`n_sim` and `labels`")
  evpi_from_draws(nb_draws(res, n_sim), res, call)
}

#' Probability that a model is useful
#'
#' For each model, at each threshold, the probability that using the model
#' gives a higher net benefit than both default strategies (treat all and
#' treat none):
#' \deqn{P(useful) = P\{NB_{model} > \max(NB_{all}, 0)\} = P(\Delta NB > 0),}
#' estimated as the share of the posterior draws in which the model's net
#' benefit is strictly higher. A model that ties with a default strategy (e.g.
#' at a threshold where it treats everyone, or no one) is not counted as
#' useful.
#'
#' P(useful) compares each model only with treat all and treat none, so a
#' model's P(useful) does not depend on the other models; to compare the
#' models with each other, use [p_best()].
#'
#' [evpi()] computes P(useful) from its own draws: `p_useful(e)`, with `e` a
#' result of [evpi()], returns it (with no new draws, and no other
#' arguments).
#'
#' @inheritSection evpi Uncertainty
#' @param x A result of [nb_val()], for one model or several, or a result of
#'   [evpi()].
#' @param ... More results of [nb_val()] (see [evpi()]).
#' @param n_sim Number of draws. Default 1000.
#' @param labels Optional names of the models, one per model, in order.
#' @return An object of class `voipred_p_useful`, a list with `estimates` (a
#'   data frame with columns `threshold`, `model` and `p_useful`, one row per
#'   model and threshold), `models` (the labels of the models), `outcome`,
#'   `censoring`, `prior`, `draws`, `n_sim` and `call`, with `print()`,
#'   `summary()`, `plot()` and `as.data.frame()` methods.
#' @seealso [p_best()], [evpi()], [plot.voipred_p_useful()], [summary.voipred_voi()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb_val(data = df_val_binary,
#'               predicted_risk = df_val_binary$p,
#'               outcome = "day30",
#'               thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' pu <- p_useful(res)
#' pu
#' summary(pu)
#' plot(pu)
#' # the same draws as evpi()
#' set.seed(1)
#' e <- evpi(res)
#' all.equal(p_useful(e)$estimates, pu$estimates)
#' @export
p_useful <- function(x, ..., n_sim = 1000, labels = NULL) {
  call <- match.call()
  if (inherits(x, "voipred_evpi"))
    return(from_evpi(x, "p_useful", list(...), missing(n_sim) && is.null(labels)))
  res <- collect_voi(x, list(...), labels, "p_useful()", "`n_sim` and `labels`")
  p_useful_from_draws(nb_draws(res, n_sim), res, call)
}

#' Probability that each strategy is the best
#'
#' For each strategy (each model, treat all and treat none) and threshold,
#' the probability that it has the highest net benefit:
#' \deqn{P(best_s) = P\{NB_s = \max_t NB_t\},}
#' estimated as the share of the posterior draws in which it is the best.
#' Exact ties are shared equally, so the probabilities of all strategies add up
#' to 1 at each threshold.
#'
#' [evpi()] computes P(best) from its own draws: `p_best(e)`, with `e` a
#' result of [evpi()], returns it (with no new draws, and no other
#' arguments).
#'
#' @inheritSection evpi Uncertainty
#' @inheritParams evpi
#' @param x A result of [nb_val()], for one model or several, or a result of
#'   [evpi()].
#' @return An object of class `voipred_p_best`, a list with `estimates` (a
#'   data frame with columns `threshold`, `strategy` and `p_best`, one row per
#'   threshold and strategy), `strategies`, `models`, `outcome`, `censoring`, `prior`,
#'   `draws`, `n_sim` and `call`, with `print()`, `summary()`, `plot()` and
#'   `as.data.frame()` methods.
#' @seealso [p_useful()], [evpi()], [plot.voipred_p_best()], [summary.voipred_voi()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse,
#'             family = binomial, data = df_dev_binary)
#' df_val_binary$p1 <- predict(fit1, newdata = df_val_binary, type = "response")
#' fit2 <- glm(day30 ~ age, family = binomial, data = df_dev_binary)
#' df_val_binary$p2 <- predict(fit2, newdata = df_val_binary, type = "response")
#' z <- seq(0.01, 0.3, by = 0.01)
#' res1 <- nb_val(data = df_val_binary,
#'                predicted_risk = df_val_binary$p1,
#'                outcome = "day30",
#'                thresholds = z,
#'                label = "Full")
#' res2 <- nb_val(data = df_val_binary,
#'                predicted_risk = df_val_binary$p2,
#'                outcome = "day30",
#'                thresholds = z,
#'                label = "Age only")
#' set.seed(1)
#' pb <- p_best(res1, res2)
#' pb
#' summary(pb)
#' plot(pb)
#' @export
p_best <- function(x, ..., n_sim = 1000, labels = NULL) {
  call <- match.call()
  if (inherits(x, "voipred_evpi"))
    return(from_evpi(x, "p_best", list(...), missing(n_sim) && is.null(labels)))
  res <- collect_voi(x, list(...), labels, "p_best()", "`n_sim` and `labels`")
  p_best_from_draws(nb_draws(res, n_sim), res, call)
}

# ---- draws (internal) ---------------------------------------------------------

# Results of nb_val() whose posterior is conjugate: binary outcome, no case weights.
is_conjugate <- function(res) res$outcome$type == "binary" && is.null(res$weights)

draws_name <- function(res) if (is_conjugate(res)) "posterior" else "Bayesian bootstrap"

# Posterior draws of the net benefit of every strategy for the models of the
# result `res` (from collect_voi()): a list with `nb` (draw x threshold x
# strategy), `thresholds`, `strategies`, `models`, `draws` and `n_sim`.
nb_draws <- function(res, n_sim) {
  check_n_sim(n_sim)
  if (is_conjugate(res)) return(finish_draws(res, bayes_draws(res, n_sim), "posterior"))
  v <- voi_info(res)
  n <- nrow(v$risk)
  W <- vapply(seq_len(n_sim), function(i) bb_weights(n, v$weights), numeric(n))
  if (n == 1) W <- matrix(W, nrow = 1)
  draws_from(res, W, "Bayesian bootstrap")
}

# The draws for the models of `res` and given patient weights W (n x draws).
draws_from <- function(res, W, draws) finish_draws(res, nb_eval(res, W), draws)

finish_draws <- function(res, nbs, draws) {
  v <- voi_info(res)
  bad <- apply(nbs, 2, anyNA)
  if (any(bad))
    warning("The net benefit could not be estimated in some draws at ",
            if (sum(bad) == 1) "threshold " else "thresholds ",
            format_thresholds(v$thresholds[bad]),
            " (no treated patient followed up to the horizon); the value of ",
            "information is NA there.", call. = FALSE)
  list(nb = nbs, thresholds = v$thresholds,
       strategies = c(res$models, "Treat All", "Treat None"), models = res$models,
       draws = draws, n_sim = dim(nbs)[1])
}

# The parts every value-of-information result shares.
voi_result <- function(estimates, d, res, call, cls, ...) {
  r1 <- res
  structure(c(list(call = call, estimates = estimates), list(...),
              list(models = d$models, outcome = r1$outcome, prior = r1$prior,
                   censoring = if (!is.null(r1) && r1$outcome$type != "binary")
                     censoring_text(r1),
                   case_weights = !is.null(r1$weights), draws = d$draws, n_sim = d$n_sim)),
            class = cls)
}

evpi_from_draws <- function(d, res = NULL, call = NULL) {
  nbs <- d$nb
  S <- dim(nbs)[3]
  enb <- apply(nbs, c(2, 3), mean)
  if (!is.matrix(enb)) enb <- matrix(enb, nrow = dim(nbs)[2])
  colnames(enb) <- d$strategies
  perfect <- colMeans(matrix(apply(nbs, c(1, 2), max), nrow = dim(nbs)[1]))
  current <- apply(enb, 1, max)
  default <- pmax(enb[, S - 1], 0)
  est <- data.frame(threshold = d$thresholds, ENB_current = current,
                    ENB_perfect = perfect, EVPI = perfect - current,
                    rEVPI = relative_evpi(perfect, current, default))
  ENB <- data.frame(threshold = d$thresholds, enb, check.names = FALSE)
  out <- voi_result(est, d, res, call, "voipred_evpi", ENB = ENB)
  # from the same draws
  out$p_useful <- p_useful_from_draws(d, res, call)
  out$p_best <- p_best_from_draws(d, res, call)
  out
}

# P(useful) or P(best) kept in a result of evpi() (from the same draws).
from_evpi <- function(x, what, dots, no_args) {
  if (length(dots) || !no_args)
    stop(what, "() of a result of evpi() returns the ", what, " computed by evpi() ",
         "from its draws, and takes no other arguments. To compute it again, give ",
         "results of nb_val().", call. = FALSE)
  if (is.null(x[[what]]))
    stop("This result of evpi() has no ", what, ": it was made by an older version ",
         "of voipred; run evpi() again.", call. = FALSE)
  x[[what]]
}

p_useful_from_draws <- function(d, res = NULL, call = NULL) {
  nbs <- d$nb
  S <- dim(nbs)[3]
  default <- pmax(nbs[, , S - 1, drop = FALSE], 0)[, , 1]
  est <- do.call(rbind, lapply(seq_along(d$models), function(m) {
    better <- nbs[, , m] > default
    if (!is.matrix(better)) better <- matrix(better, nrow = dim(nbs)[1])
    data.frame(threshold = d$thresholds, model = d$models[m],
               p_useful = colMeans(better))
  }))
  rownames(est) <- NULL
  voi_result(est, d, res, call, "voipred_p_useful")
}

p_best_from_draws <- function(d, res = NULL, call = NULL) {
  nbs <- d$nb
  dims <- dim(nbs)
  best <- apply(nbs, c(1, 2), max)                      # draw x threshold
  share <- array(0, dims)
  for (s in seq_len(dims[3])) share[, , s] <- nbs[, , s] >= best - 1e-12
  ties <- apply(share, c(1, 2), sum)
  for (s in seq_len(dims[3])) share[, , s] <- share[, , s] / ties
  est <- do.call(rbind, lapply(seq_len(dims[3]), function(s)
    data.frame(threshold = d$thresholds, strategy = d$strategies[s],
               p_best = colMeans(matrix(share[, , s], nrow = dims[1])))))
  rownames(est) <- NULL
  voi_result(est, d, res, call, "voipred_p_best", strategies = d$strategies)
}

#' Expected value of sample information (EVSI) for model validation
#'
#' The expected gain in net benefit from a future validation study of size
#' \eqn{N} (`future_n`), added to the current sample, before deciding which
#' strategy to use. EVSI is between 0 and EVPI, and approaches EVPI as
#' \eqn{N} grows: comparing it with the cost of the study helps decide whether
#' more validation data are worth collecting.
#'
#' For each of `n_sim` posterior draws (see "Uncertainty"), which define a
#' "true" population with the true net benefit of every strategy, and for each
#' \eqn{N}, a future sample of \eqn{N} patients is drawn from that population
#' and combined with the current data, the strategy with the highest net
#' benefit in the combined data is chosen, and its *true* net benefit
#' recorded (Sadatsafavi et al. 2025):
#' * binary outcome: the conjugate algorithm. The future sample is simulated
#'   from the draw (events, and the patients treated among them, binomial /
#'   multinomial), and the strategy with the highest posterior mean net benefit
#'   after updating the posterior with it is chosen;
#' * survival and competing-risks outcomes (and case weights): the bootstrap
#'   algorithm. The draw's Bayesian-bootstrap weights \eqn{W} define the
#'   population; the future sample is \eqn{N} patients drawn from it
#'   (multinomial counts with probabilities \eqn{W}, a future patient being a
#'   current patient, with the same follow-up), pooled with the current sample,
#'   and the net benefit recomputed with the same estimator as [nb_val()]. With
#'   cIPCW (`cipcw`) the censoring model is refitted in every draw and
#'   every pooled sample (see [cipcw_control()]).
#'
#' If strategies tie exactly in the combined data, the mean of their true net
#' benefits is recorded (a random choice among them). Then \eqn{ENB_{sample}}
#' is the mean of the recorded true net benefits and
#' \eqn{EVSI = ENB_{sample} - ENB_{current}}, with \eqn{ENB_{current}} as in
#' [evpi()] (from the same draws). The true EVSI is never negative, but its
#' estimate can be slightly negative where EVPI is close to 0 (Monte Carlo
#' error; more draws make it smaller): read such values as 0.
#'
#' @inheritSection evpi Uncertainty
#' @inheritParams evpi
#' @param future_n Sizes of the future validation sample: one or more
#'   positive integers.
#' @return An object of class `voipred_evsi`, a list with `estimates` (a data
#'   frame with columns `threshold`, `future_n`, `EVSI` and `EVPI`, one row per
#'   threshold and future sample size), `future_n`, `models`, `outcome`,
#'   `censoring`, `prior`, `draws`, `n_sim` and `call`. It has `print()`, `summary()`,
#'   `plot()` and `as.data.frame()` methods; EVSI can also be drawn under the
#'   decision curve with `plot(<nb_val result>, evsi = )`.
#' @references Sadatsafavi M, et al. Expected value of sample information
#'   calculations for risk prediction model validation. Medical Decision
#'   Making 2025; 45(3): 232-244.
#' @seealso [evpi()], [plot.voipred_evsi()], [summary.voipred_voi()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb_val(data = df_val_binary,
#'               predicted_risk = df_val_binary$p,
#'               outcome = "day30",
#'               thresholds = c(0.02, 0.05, 0.1))
#' set.seed(1)
#' o <- evsi(res, future_n = c(250, 500, 1000, 4000))
#' o
#' summary(o)
#' plot(o)
#' @export
evsi <- function(x, ..., future_n, n_sim = 1000, labels = NULL) {
  call <- match.call()
  if (missing(future_n) || !is.numeric(future_n) || !length(future_n) ||
      anyNA(future_n) || any(future_n < 1) || any(future_n != round(future_n)))
    stop("`future_n` must be one or more positive whole numbers: the sizes of the ",
         "future validation sample.", call. = FALSE)
  check_n_sim(n_sim)
  res <- collect_voi(x, list(...), labels, "evsi()", "`future_n`, `n_sim` and `labels`")
  r <- if (is_conjugate(res)) bayes_evsi(res, future_n, n_sim)
       else evsi_bootstrap(res, future_n, n_sim)
  z <- nb_thresholds(res)
  Z <- length(z)
  est <- data.frame(threshold = rep(z, length(future_n)),
                    future_n = rep(future_n, each = Z),
                    EVSI = as.vector(r$EVSI), EVPI = rep(r$EVPI, length(future_n)))
  d <- list(models = res$models, draws = draws_name(res), n_sim = n_sim)
  voi_result(est, d, res, call, "voipred_evsi", future_n = future_n)
}

# EVSI with the bootstrap algorithm of Sadatsafavi et al. (2025), Bayesian
# bootstrap version: a Z x J matrix of EVSI and the EVPI of the same draws.
evsi_bootstrap <- function(res, future_n, n_sim) {
  v <- voi_info(res)
  n <- nrow(v$risk)
  w0 <- if (is.null(v$weights)) rep(1, n) else v$weights
  Z <- length(v$thresholds)
  J <- length(future_n)
  S <- length(res$models) + 2
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
      u <- stats::rgamma(n, 1, 1)
      if (!is.null(v$weights)) u <- u * v$weights
      W[, b] <- u / sum(u)
      for (j in seq_len(J))
        Wf[[j]][, b] <- stats::rmultinom(1, future_n[j], W[, b])
    }
    # the true net benefits; with cIPCW the weights are put on the scale of
    # patients (mean 1), as the censoring-model refit sees them
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
  list(EVSI = evsi_zj, EVPI = evpi_z)
}

# ---- internal ----------------------------------------------------------------

# Results of nb_val() that can be analysed together, merged into one.
collect_voi <- function(x, more, labels, fun, args) {
  collect_nb(x, more, labels, fun = fun, args = args)
}

check_n_sim <- function(n_sim) {
  if (!is.numeric(n_sim) || length(n_sim) != 1 || is.na(n_sim) || n_sim < 2 ||
      n_sim != round(n_sim))
    stop("`n_sim` must be a whole number of at least 2.", call. = FALSE)
}

# One draw of Bayesian-bootstrap patient weights, Dirichlet(1, ..., 1) (times
# the case weights), scaled to sum to n.
bb_weights <- function(n, w0 = NULL) {
  u <- stats::rgamma(n, 1, 1)
  if (!is.null(w0)) u <- u * w0
  u * n / sum(u)
}

# Net benefit of every strategy for each column of weights W (n x B):
# an array draw x threshold x strategy (models, treat all, treat none).
nb_eval <- function(res, W) {
  v <- voi_info(res)
  z <- v$thresholds
  B <- ncol(W)
  Z <- length(z)
  M <- ncol(v$risk)
  out <- array(NA_real_, c(B, Z, M + 2))
  out[, , M + 2] <- 0
  if (v$type == "binary") {
    odds <- z / (1 - z)
    gain <- outer(v$y, rep(1, Z)) - outer(1 - v$y, odds)      # Y - (1 - Y) z/(1-z)
    sw <- colSums(W)
    out[, , M + 1] <- crossprod(W, gain) / sw
    for (m in seq_len(M)) {
      p <- v$risk[, m]
      out[, , m] <- crossprod(W, gain * (outer(p, z, ">=") & p > 0)) / sw
    }
  } else {
    # one censoring model, refitted once per draw, for all the models
    refit <- censoring_refitter(v, W)
    for (b in seq_len(B)) {
      cw <- if (is.null(refit)) NULL else refit(W[, b])
      for (m in seq_len(M)) {
        r <- nb_one(v, v$risk[, m], W[, b], cw)
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

# Net benefit of the model with risks p, with patient weights w (and the
# censoring weights cw of cIPCW).
nb_one <- function(v, p, w, cw = NULL) {
  if (v$censoring == "km")
    return(suppressWarnings(nb_core_km(v$time, v$status, v$k, p, v$thresholds,
                                       v$horizon, w)))
  nb_core_cipcw(v$time, v$status, v$k, p, v$thresholds, v$horizon,
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
censoring_refitter <- function(v, W = NULL) {
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
  } else {
    stop("The censoring model (", class(model)[1], ") must be refitted in each draw ",
         "with the draw's weights (Bayesian bootstrap), and voipred does not know how ",
         "to fit a ", class(model)[1], " model with case weights. Give ",
         "`refit_function = function(data, weights)` in cipcw_control(), returning ",
         "the censoring model fitted with these case weights.", call. = FALSE)
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

relative_evpi <- function(perfect, current, default) {
  r <- (perfect - default) / (current - default)
  r[is.nan(r)] <- NA
  r
}

