#' Settings for conditional IPCW net benefit
#'
#' `cipcw_control()` says how the censoring weights are obtained when
#' [nb()] is called with `censoring = "cipcw"`. Conditional inverse
#' probability of censoring weighting (cIPCW) allows censoring to depend on
#' covariates: it assumes censoring is independent of the event time *given*
#' the covariates of the censoring model.
#'
#' There are two ways to give the censoring information:
#'
#' **1. A fitted censoring model** (`model =`). Fit any survival model to the
#' *censoring* times in the validation data, i.e. with the censored patients
#' as the "events": `Surv(time, status == 0)`. With competing risks, only true
#' censoring is censoring: patients with a competing event are not censored
#' (`Surv(time, status == "censored")`, not `status != cause`). voipred
#' predicts from the model the probability of remaining uncensored just before
#' each patient's own time and at the horizon, in the same way as the risks
#' of [predict_risk()]: through `riskRegression::predictRisk()` (e.g.
#' `coxph` with `x = TRUE`, `rms::cph`, `survreg`, `flexsurvreg`, `ranger`,
#' `randomForestSRC`, a Kaplan-Meier `prodlim` fit), or voipred's own code for
#' models riskRegression cannot handle (`coxph` without `x = TRUE`, `aorsf`,
#' `partykit` survival trees, `rstpm2::stpm2`, `mgcv::gam` with `cox.ph()`).
#' For any other model, give `predict_function` (or the weights, option 2).
#' Where the model's
#' outcome can be read from its formula, voipred checks it against `data`: it
#' stops if the model was fitted to the event (for competing risks: to the
#' cause of interest, to another cause or to all events together), and warns
#' if the check is not possible or if the model was fitted to a different
#' number of patients than `data` has.
#'
#' **2. Weights** (`weight_event =`, `weight_horizon =`): one value per row of
#' `data`, the inverse probabilities of remaining uncensored just before each
#' patient's own time and at the horizon, i.e. \eqn{1/\hat G(T_i- \mid X_i)}
#' and \eqn{1/\hat G(\tau \mid X_i)}. Only the weights of patients with an
#' event by the horizon (`weight_event`; including competing events) and of
#' patients event-free at the horizon (`weight_horizon`) are used. They must be
#' **unstabilized**: stabilized weights bias this estimator.
#'
#' **Flexible censoring models.** A flexible model (e.g. a forest) fitted and
#' evaluated on the same patients overfits their censoring times, and its
#' weights are biased. Use cross-fitting: split `data` into folds, fit the
#' censoring model on all folds but one, predict the censoring probabilities
#' for the held-out fold, and give the resulting weights through option 2.
#'
#' **Truncation** is off by default. With `min_prob`, censoring probabilities
#' below it (option 1) or weights above `1 / min_prob` (option 2) are
#' truncated, with a warning saying how many of the used ones were. Without
#' it, voipred warns when some used probabilities are below 0.01 (weights
#' above 100), a sign of near-violations of positivity, and stops if a
#' probability is exactly 0.
#'
#' **Uncertainty.** The censoring model or weights are taken as fixed: in
#' resampling (e.g. bootstrap draws for the value of information) the same
#' weights are used in every draw, which ignores the uncertainty in the
#' censoring model.
#'
#' @param model Option 1: a survival model fitted to the censoring times of
#'   the validation data `data` passed to [nb()] (same patients), with
#'   censoring as the event, e.g. `coxph(Surv(time, status == 0) ~ x,
#'   data = val)`. Give either `model` or the two weights, not both.
#' @param weight_event,weight_horizon Option 2: two numeric vectors, one
#'   positive finite value per row of `data` (in its row order), given
#'   together: `weight_event` is the inverse probability of remaining
#'   uncensored just before each patient's own time, \eqn{1/\hat G(T_i- \mid
#'   X_i)}, used for patients with an event (of any cause) by the horizon;
#'   `weight_horizon` is the inverse probability of remaining uncensored at
#'   the horizon, \eqn{1/\hat G(\tau \mid X_i)}, used for patients event-free
#'   at the horizon. Values for the other patients are ignored (any positive
#'   number will do). Unstabilized weights only.
#' @param min_prob Optional truncation, for both options: censoring
#'   probabilities below `min_prob` (option 1), or weights above
#'   `1 / min_prob` (option 2), are set to that limit. A single number
#'   strictly between 0 and 1. Default `NULL`: no truncation.
#' @param predict_function Option 1, optional: a `function(model, newdata,
#'   times)` returning, for each row of `newdata`, the probability of
#'   **remaining uncensored** at that row's time in `times` (the censoring
#'   survival \eqn{G(t \mid X)}, not a risk). It may instead return a curve
#'   per row, `list(time = <grid>, surv = <n x grid matrix>)`, which voipred
#'   reads at each row's own time. Use it for a censoring model voipred cannot
#'   predict from by itself, or to predict in your own way; it then replaces
#'   voipred's own prediction. Requires `model`.
#' @param ties Option 1 only, and only matters for a censoring model whose
#'   survival curve has steps (e.g. Cox, Kaplan-Meier, forests, trees):
#'   `"left"` (default) uses the probability of remaining uncensored just
#'   before each patient's own time, so a censoring step at exactly that time
#'   is not included (the usual convention, as for Kaplan-Meier: an event
#'   tied with a censoring comes first); `"right"` includes it. For smooth
#'   (parametric) censoring models both give practically the same result.
#' @return An object of class `cipcw_control`.
#' @examples
#' \dontrun{
#' library(survival)
#' cens_fit <- coxph(Surv(time, status == 0) ~ age + nodes, data = val)
#' nb(fit, val, ~ Surv(time, status), horizon = 5,
#'    censoring = "cipcw", cipcw = cipcw_control(model = cens_fit))
#' }
#' @export
cipcw_control <- function(model = NULL, weight_event = NULL,
                          weight_horizon = NULL, min_prob = NULL,
                          ties = c("left", "right"), predict_function = NULL) {
  ties <- match.arg(ties)
  has_w <- !is.null(weight_event) || !is.null(weight_horizon)
  if (!is.null(model) && has_w)
    stop("Give either a censoring `model` or the weights, not both.",
         call. = FALSE)
  if (!is.null(predict_function)) {
    if (!is.function(predict_function))
      stop("`predict_function` must be a function(model, newdata, times).",
           call. = FALSE)
    if (is.null(model))
      stop("`predict_function` needs the censoring `model` it predicts from.",
           call. = FALSE)
  }
  if (xor(is.null(weight_event), is.null(weight_horizon)))
    stop("`weight_event` and `weight_horizon` must be given together.",
         call. = FALSE)
  for (nm in c("weight_event", "weight_horizon")) {
    v <- get(nm)
    if (!is.null(v) && (!is.numeric(v) || anyNA(v) || any(!is.finite(v)) ||
                        any(v <= 0)))
      stop("`", nm, "` must be positive, finite numbers without missing values.",
           call. = FALSE)
  }
  if (has_w && length(weight_event) != length(weight_horizon))
    stop("`weight_event` and `weight_horizon` must have the same length.",
         call. = FALSE)
  if (!is.null(model) && (is.numeric(model) || is.data.frame(model)))
    stop("`model` must be a fitted censoring model. For weights, use ",
         "`weight_event` and `weight_horizon`.", call. = FALSE)
  if (!is.null(min_prob) && (!is.numeric(min_prob) || length(min_prob) != 1 ||
                             is.na(min_prob) || min_prob <= 0 || min_prob >= 1))
    stop("`min_prob` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  structure(list(model = model, predict_function = predict_function,
                 weight_event = weight_event, weight_horizon = weight_horizon,
                 min_prob = min_prob, ties = ties),
            class = "cipcw_control")
}

#' @export
print.cipcw_control <- function(x, ...) {
  how <- if (!is.null(x$model)) paste0("censoring model (", class(x$model)[1],
                                       if (!is.null(x$predict_function))
                                         ", with predict_function", ")")
         else if (!is.null(x$weight_event)) "given weights" else "not set"
  cat("cIPCW settings: ", how, ", truncation: ",
      if (is.null(x$min_prob)) "off" else paste0("min_prob = ", x$min_prob),
      ", ties = ", x$ties, "\n", sep = "")
  invisible(x)
}

is_default_cipcw <- function(cipcw) {
  is.null(cipcw$model) && is.null(cipcw$weight_event) &&
    is.null(cipcw$min_prob) && identical(cipcw$ties, "left")
}

# The two cIPCW weight vectors for `data`. Only the weights that the estimator
# uses (event-by-horizon rows, event-free rows) are checked and reported.
cipcw_weights <- function(cipcw, data, info, horizon, k) {
  n <- nrow(data)
  time <- info$time
  needs_event <- info$status != 0 & time <= horizon
  needs_horizon <- time > horizon

  if (!is.null(cipcw$weight_event)) {
    if (length(cipcw$weight_event) != n)
      stop("The cIPCW weights have length ", length(cipcw$weight_event),
           " but `data` has ", n, " rows.", call. = FALSE)
    we <- cipcw$weight_event
    wh <- cipcw$weight_horizon
    used <- c(we[needs_event], wh[needs_horizon])
    if (!is.null(cipcw$min_prob)) {
      cap <- 1 / cipcw$min_prob
      n_trunc <- sum(used > cap)
      if (n_trunc > 0)
        warning(n_trunc, " cIPCW weight", if (n_trunc > 1) "s were" else " was",
                " truncated at 1 / min_prob = ", signif(cap, 4), ".", call. = FALSE)
      we <- pmin(we, cap); wh <- pmin(wh, cap)
    } else if (any(used > 100)) {
      warning(sum(used > 100), " cIPCW weight", if (sum(used > 100) > 1) "s are" else " is",
              " above 100 (censoring probability below 0.01), a sign of near-",
              "violations of positivity. Consider truncation with `min_prob`.",
              call. = FALSE)
    }
    return(list(weight_event = we, weight_horizon = wh))
  }

  model <- cipcw$model
  fun <- cipcw$predict_function
  # with predict_function the model may be anything, so its covariates are not
  # checked; whether it was fitted to censoring is still checked when possible
  if (is.null(fun)) check_covariates(model, data, what = "censoring model")
  check_censoring_model(model, data, info, k)
  g_at <- function(rows, times) {
    nd <- data[rows, , drop = FALSE]
    if (is.null(fun)) return(surv_at(model, nd, times, cipcw$ties))
    censoring_from_function(fun, model, nd, times, cipcw$ties)
  }
  g_event <- rep(1, n); g_horizon <- rep(1, n)
  if (any(needs_event))
    g_event[needs_event] <- g_at(which(needs_event), time[needs_event])
  if (any(needs_horizon))
    g_horizon[needs_horizon] <- g_at(which(needs_horizon), rep(horizon, sum(needs_horizon)))
  used <- c(g_event[needs_event], g_horizon[needs_horizon])
  if (anyNA(used) || any(used < 0 | used > 1))
    stop("The censoring model returned invalid censoring probabilities ",
         "(missing, or outside [0, 1]).", call. = FALSE)
  if (!is.null(cipcw$min_prob)) {
    n_trunc <- sum(used < cipcw$min_prob)
    if (n_trunc > 0)
      warning(n_trunc, " censoring probabilit", if (n_trunc > 1) "ies were" else "y was",
              " truncated at min_prob = ", cipcw$min_prob, ". This signals ",
              "near-violations of positivity.", call. = FALSE)
    g_event <- pmax(g_event, cipcw$min_prob)
    g_horizon <- pmax(g_horizon, cipcw$min_prob)
  } else {
    if (any(used == 0))
      stop(sum(used == 0), " censoring probabilit", if (sum(used == 0) > 1) "ies are" else "y is",
           " exactly 0, which gives an infinite weight (positivity is violated). ",
           "Set `min_prob` to truncate them.", call. = FALSE)
    if (any(used < 0.01))
      warning(sum(used < 0.01), " censoring probabilit",
              if (sum(used < 0.01) > 1) "ies are" else "y is",
              " below 0.01, a sign of near-violations of positivity. ",
              "Consider truncation with `min_prob`.", call. = FALSE)
  }
  list(weight_event = 1 / g_event, weight_horizon = 1 / g_horizon)
}

# Censoring probabilities from a user's predict_function(model, newdata, times):
# either one probability per row (at that row's time), or a curve per row,
# list(time = <grid>, surv = <n x grid matrix>), read at each row's own time.
censoring_from_function <- function(fun, model, newdata, times, ties) {
  n <- nrow(newdata)
  eps <- if (ties == "left") 1e-9 * max(1, max(times)) else 0
  out <- fun(model, newdata, times - eps)
  if (is.list(out) && !is.data.frame(out) && "time" %in% names(out)) {
    surv <- as.matrix(out$surv)
    if (nrow(surv) != n || ncol(surv) != length(out$time))
      stop("The censoring `predict_function` returned a curve of the wrong ",
           "size: `surv` must be ", n, " rows by length(time) columns.", call. = FALSE)
    i <- if (ties == "left") findInterval(times, out$time, left.open = TRUE)
         else findInterval(times, out$time)
    return(ifelse(i == 0, 1, surv[cbind(seq_len(n), pmax(i, 1))]))
  }
  if (is.data.frame(out) || is.matrix(out)) {
    if (ncol(as.matrix(out)) != 1)
      stop("The censoring `predict_function` must return one probability per ",
           "row (or a curve list(time, surv)).", call. = FALSE)
    out <- as.matrix(out)[, 1]
  }
  if (!is.numeric(out) || length(out) != n)
    stop("The censoring `predict_function` returned ", length(out),
         " values for ", n, " rows; it must return one probability of remaining ",
         "uncensored per row.", call. = FALSE)
  as.vector(unname(out))
}

# Stop if a censoring model was clearly fitted to the event rather than to
# censoring; warn if this cannot be verified. Generic: the model's outcome is
# read from its formula (evaluated in `data`), and its size with insight.
check_censoring_model <- function(model, data, info, k) {
  y <- model_lhs(model, data)
  if (!inherits(y, c("Surv", "Hist"))) {
    warning("Could not verify that the censoring model was fitted to the ",
            "censoring times (its outcome could not be read from its formula). ",
            "It must be fitted with the censored patients as events, e.g. ",
            "Surv(time, status == 0).", call. = FALSE)
    return(invisible(FALSE))
  }
  y <- as.integer(outcome_status(y) != 0)
  st <- info$status
  refit <- " Refit it with only censoring as the event, e.g. Surv(time, status == 0)."
  if (!identical(y, as.integer(st == 0))) {
    if (info$type == "survival" && identical(y, as.integer(st == 1)))
      stop("The censoring model appears to have been fitted to the event, ",
           "not to censoring.", refit, call. = FALSE)
    if (info$type == "competing") {
      if (identical(y, as.integer(st == k)))
        stop("The censoring model appears to have been fitted to the cause of ",
             "interest (", info$states[k], "), not to censoring.", refit, call. = FALSE)
      if (identical(y, as.integer(st != 0)))
        stop("The censoring model appears to have been fitted to all events ",
             "together, not to censoring.", refit, call. = FALSE)
      for (j in setdiff(seq_along(info$states), k))
        if (identical(y, as.integer(st == j)))
          stop("The censoring model appears to have been fitted to a competing ",
               "cause (", info$states[j], "), not to censoring.", refit, call. = FALSE)
    }
    warning("Could not verify that the censoring model was fitted to the ",
            "censoring times of `data`: its outcome does not match. It must be ",
            "fitted with the censored patients as events, e.g. ",
            "Surv(time, status == 0).", call. = FALSE)
    return(invisible(FALSE))
  }
  n_fit <- tryCatch(suppressWarnings(insight::n_obs(model)), error = function(e) NULL)
  if (is.numeric(n_fit) && length(n_fit) == 1 && n_fit != nrow(data))
    warning("The censoring model was fitted to ", n_fit, " patients, but `data` ",
            "has ", nrow(data), ". It should be fitted to the censoring times of ",
            "`data` (the validation sample).", call. = FALSE)
  invisible(TRUE)
}
