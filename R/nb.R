#' Net benefit of a prediction model in a validation sample
#'
#' Computes the net benefit of using the model (treat if predicted risk is at
#' or above the threshold), of treating everyone, and of treating no one, at
#' each risk threshold. This is the core quantity every other voipred
#' validation metric (ΔNB, P(useful), P(best), EVPI, EVSI) is built on.
#'
#' With outcome \eqn{Y}, predicted risk \eqn{\hat\pi} and threshold \eqn{z},
#' \deqn{NB_{model}(z) = \frac{1}{n}\sum_i Y_i T_i - \frac{1}{n}\sum_i (1-Y_i) T_i \frac{z}{1-z},\quad T_i = 1(\hat\pi_i \ge z)}
#' \deqn{NB_{all}(z) = \bar Y - (1-\bar Y)\frac{z}{1-z}, \qquad NB_{none}(z) = 0.}
#' With `weights`, the means are weighted means. Following dcurves, a patient
#' whose predicted risk is exactly 0 is never treated, even at \eqn{z = 0}.
#'
#' @param model A fitted model that [predict_risk()] understands (`glm`,
#'   `ranger`, `randomForest`, bayespred). A vector of predicted risks is not
#'   accepted.
#' @param data A data frame: the validation sample.
#' @param outcome The binary outcome: a column name in `data`, or a vector with
#'   one value per row. The event is 1: accepts 0/1, logical (TRUE = 1), or a
#'   two-level factor (level `"1"`, or the second level if the levels are not
#'   0/1).
#' @param thresholds Risk thresholds in \[0, 1).
#' @param weights Optional non-negative observation weights.
#' @param bayespred_method Only for bayespred models: `"pm"` (default, the
#'   posterior mean of risk) or `"pe"` (the point estimate). See
#'   [predict_risk()].
#' @return A data frame of class `voipred_nb` with one row per threshold and
#'   columns `threshold`, `NB_model`, `NB_all`, `NB_none`. Attributes hold
#'   `n`, `prevalence`, and the predicted risks (`risk`).
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:10000, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, outcome = "day30", thresholds = seq(0, 0.3, by = 0.01))
#' head(res)
#' @export
nb <- function(model, data, outcome, thresholds = seq(0, 0.99, by = 0.01),
               weights = NULL, bayespred_method = "pm") {
  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  y <- get_outcome(outcome, data)
  thresholds <- check_thresholds(thresholds)
  p <- predict_risk(model, data, bayespred_method = bayespred_method)
  if (length(p) != length(y))
    stop("Got ", length(p), " predicted risks for ", length(y), " outcomes.",
         call. = FALSE)
  w <- check_weights(weights, length(y))

  res <- nb_core(y, p, thresholds, w)
  out <- data.frame(threshold = thresholds,
                    NB_model = res$NB_model,
                    NB_all = res$NB_all,
                    NB_none = 0)
  attr(out, "n") <- length(y)
  attr(out, "prevalence") <- res$prevalence
  attr(out, "risk") <- p
  class(out) <- c("voipred_nb", "data.frame")
  out
}

#' @export
print.voipred_nb <- function(x, ...) {
  cat("Net benefit: n = ", attr(x, "n"), ", prevalence = ",
      format(attr(x, "prevalence"), digits = 3), ", ", nrow(x),
      " thresholds\n", sep = "")
  y <- x
  attr(y, "risk") <- NULL
  print(as.data.frame(unclass(y))[, c("threshold", "NB_model", "NB_all", "NB_none")],
        row.names = FALSE, ...)
  invisible(x)
}


# ---- internal ---------------------------------------------------------------

# Net benefit for all thresholds at once, O(n log n). Weights let the same core
# serve the (Bayesian) bootstrap later. Treat if p >= z; a risk of exactly 0 is
# never treated, even at z = 0. Both rules are the ones dcurves uses.
nb_core <- function(y, p, thresholds, w = NULL) {
  if (is.null(w)) w <- rep(1, length(y))
  sw <- sum(w)
  p[p == 0] <- -1
  o <- order(p)
  ps <- p[o]
  cum_ev <- c(0, cumsum((w * y)[o]))
  cum_ne <- c(0, cumsum((w * (1 - y))[o]))
  tot_ev <- cum_ev[length(cum_ev)]
  tot_ne <- cum_ne[length(cum_ne)]
  # number of sorted risks strictly below each threshold
  k <- findInterval(thresholds, ps, left.open = TRUE)
  tp <- (tot_ev - cum_ev[k + 1]) / sw
  fp <- (tot_ne - cum_ne[k + 1]) / sw
  odds <- thresholds / (1 - thresholds)
  prev <- tot_ev / sw
  list(NB_model = tp - fp * odds,
       NB_all = prev - (1 - prev) * odds,
       prevalence = prev)
}

get_outcome <- function(outcome, data) {
  if (is.character(outcome) && length(outcome) == 1) {
    if (!outcome %in% names(data))
      stop("Column '", outcome, "' is not in `data`.", call. = FALSE)
    y <- data[[outcome]]
  } else {
    y <- outcome
  }
  if (length(y) != nrow(data))
    stop("`outcome` has length ", length(y), " but `data` has ", nrow(data),
         " rows.", call. = FALSE)
  if (anyNA(y)) stop("`outcome` has missing values.", call. = FALSE)
  if (is.factor(y)) {
    if (nlevels(y) != 2)
      stop("A factor outcome must have exactly two levels.", call. = FALSE)
    ev <- if ("1" %in% levels(y)) "1" else levels(y)[2]
    y <- as.integer(y == ev)
  } else if (is.logical(y)) {
    y <- as.integer(y)
  } else {
    y <- as.numeric(unclass(y))
    if (!all(y %in% c(0, 1)))
      stop("`outcome` must be binary (0/1, logical, or a two-level factor).",
           call. = FALSE)
  }
  y
}

check_thresholds <- function(thresholds) {
  if (!is.numeric(thresholds) || length(thresholds) == 0 || anyNA(thresholds))
    stop("`thresholds` must be a non-empty numeric vector.", call. = FALSE)
  if (any(thresholds < 0 | thresholds >= 1))
    stop("`thresholds` must lie in [0, 1).", call. = FALSE)
  thresholds
}

check_weights <- function(weights, n) {
  if (is.null(weights)) return(NULL)
  if (!is.numeric(weights) || length(weights) != n || anyNA(weights) ||
      any(weights < 0) || sum(weights) <= 0)
    stop("`weights` must be ", n, " non-negative numbers with a positive sum.",
         call. = FALSE)
  weights
}
