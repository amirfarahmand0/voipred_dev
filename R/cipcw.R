#' Settings for conditional IPCW net benefit (optional)
#'
#' Conditional inverse probability of censoring weighting (cIPCW) allows
#' censoring to depend on covariates: it assumes censoring is independent of
#' the event time *given* the covariates of the censoring model. It is used
#' when [nb_val()] is given `cipcw`, and most analyses need nothing more than
#' the covariates, `cipcw = ~ age + nodes` (option 1 below, with the default
#' truncation). `cipcw_control()` is only needed for the other options: your
#' own censoring model, given weights, another truncation, `ties`, or your own
#' predict and refit functions.
#'
#' There are three ways to give the censoring information:
#'
#' **1. The covariates of the censoring** (`formula =`, or simply
#' `cipcw = ~ age + nodes` in [nb_val()]). voipred fits a Cox model of the
#' censoring times of the validation data on them,
#' `coxph(Surv(time, censored) ~ age + nodes)`, with only true censoring as
#' the event (with competing risks, patients with a competing event are not
#' censored), weighted by the case weights of [nb_val()] if any, and refits it
#' in every draw of the value-of-information functions. A two-sided formula
#' (e.g. `Surv(time, status == 0) ~ age + nodes`) is fitted as it is.
#'
#' **2. A fitted censoring model** (`model =`, or `cipcw = <model>` in
#' [nb_val()]). Fit any survival model to the
#' *censoring* times in the validation data, i.e. with the censored patients
#' as the "events": `Surv(time, status == 0)`. With competing risks, only true
#' censoring is censoring: patients with a competing event are not censored
#' (`Surv(time, status == "censored")`, not `status != cause`). voipred
#' predicts from the model the probability of remaining uncensored just before
#' each patient's own time and at the horizon, in the same way as it predicts
#' risks from a model in [nb_val()]: through `riskRegression::predictRisk()` (e.g.
#' `coxph` with `x = TRUE`, `rms::cph`, `survreg`, `flexsurvreg`, `ranger`,
#' `randomForestSRC`, a Kaplan-Meier `prodlim` fit), or voipred's own code for
#' models riskRegression cannot handle (`coxph` without `x = TRUE`, `aorsf`,
#' `partykit` survival trees, `rstpm2::stpm2`, `mgcv::gam` with `cox.ph()`).
#' For any other model, give `predict_function` (or the weights, option 3).
#' Where the model's
#' outcome can be read from its formula, voipred checks it against `data`: it
#' stops if the model was fitted to the event (for competing risks: to the
#' cause of interest, to another cause or to all events together), and warns
#' if the check is not possible or if the model was fitted to a different
#' number of patients than `data` has.
#'
#' **3. Weights** (`weight_event =`, `weight_horizon =`): one value per row of
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
#' for the held-out fold, and give the resulting weights through option 3 (in
#' [nb_val()] only: see "Uncertainty").
#'
#' **Truncation.** By default, censoring probabilities below `min_prob = 0.02`
#' (options 1 and 2), or weights above `1 / min_prob = 50` (option 3), are
#' truncated, with a warning saying how many of the used ones were. With
#' `min_prob = NULL` there is no truncation: voipred then warns when some used
#' probabilities are below 0.01 (weights above 100), a sign of near-violations
#' of positivity, and stops if a probability is exactly 0.
#'
#' **Uncertainty (value of information).** The censoring model is estimated
#' from the validation data, so its uncertainty is part of the uncertainty
#' about the net benefit. The value-of-information functions ([evpi()],
#' [evsi()]) therefore refit the censoring
#' model in every draw of the Bayesian bootstrap, with that draw's weights as
#' case weights, and recompute the censoring weights from it. The refit is
#' done:
#' * for option 1, automatically (voipred's own Cox model);
#' * with `refit_function`, if given;
#' * otherwise automatically for the censoring models whose case weights
#'   voipred knows: `survival::coxph`, `survival::survreg`, `rms::cph`,
#'   `rms::psm`, `flexsurv::flexsurvreg` (also `flexsurvspline`) and
#'   `prodlim::prodlim`: the model's own call, with the draw's data and the
#'   draw's weights as case weights;
#' * otherwise it stops and asks for `refit_function`: other models take
#'   weights under other names or with another meaning (e.g. `ranger`'s
#'   `case.weights` are sampling probabilities), so voipred does not guess.
#'
#' Patients with a weight of 0 in a draw are left out of the refit. Given
#' weights (option 3) cannot be refitted, so they can be used in [nb_val()] but not
#' in the value-of-information functions. `nb_val()` itself always uses the model
#' as given.
#'
#' @param formula Option 1: the covariates of the censoring, a one-sided
#'   formula such as `~ age + nodes` (variables of `data`), for a Cox model
#'   fitted by voipred. Give one of `formula`, `model` or the two weights.
#' @param model Option 2: a survival model fitted to the censoring times of
#'   the validation data `data` passed to [nb_val()] (same patients), with
#'   censoring as the event, e.g. `coxph(Surv(time, status == 0) ~ x,
#'   data = val)`.
#' @param weight_event,weight_horizon Option 3: two numeric vectors, one
#'   positive finite value per row of `data` (in its row order), given
#'   together: `weight_event` is the inverse probability of remaining
#'   uncensored just before each patient's own time, \eqn{1/\hat G(T_i- \mid
#'   X_i)}, used for patients with an event (of any cause) by the horizon;
#'   `weight_horizon` is the inverse probability of remaining uncensored at
#'   the horizon, \eqn{1/\hat G(\tau \mid X_i)}, used for patients event-free
#'   at the horizon. Values for the other patients are ignored (any positive
#'   number will do). Unstabilized weights only.
#' @param min_prob Truncation, for every option: censoring probabilities
#'   below `min_prob` (options 1 and 2), or weights above `1 / min_prob`
#'   (option 3), are set to that limit. A single number strictly between 0
#'   and 1; default 0.02 (weights at most 50). `NULL`: no truncation.
#' @param predict_function Option 2, optional: a `function(model, newdata,
#'   times)` returning, for each row of `newdata`, the probability of
#'   **remaining uncensored** at that row's time in `times` (the censoring
#'   survival \eqn{G(t \mid X)}, not a risk). It may instead return a curve
#'   per row, `list(time = <grid>, surv = <n x grid matrix>)`, which voipred
#'   reads at each row's own time. Use it for a censoring model voipred cannot
#'   predict from by itself, or to predict in your own way; it then replaces
#'   voipred's own prediction. Requires `model`.
#' @param ties Options 1 and 2, and only matters for a censoring model whose
#'   survival curve has steps (e.g. Cox, Kaplan-Meier, forests, trees):
#'   `"left"` (default) uses the probability of remaining uncensored just
#'   before each patient's own time, so a censoring step at exactly that time
#'   is not included (the usual convention, as for Kaplan-Meier: an event
#'   tied with a censoring comes first); `"right"` includes it. For smooth
#'   (parametric) censoring models both give practically the same result.
#' @param refit_function Option 2, optional: a `function(data, weights)`
#'   returning the censoring model fitted to `data` (the validation data given
#'   to [nb_val()], restricted to the patients with a positive weight in a draw)
#'   with case weights `weights` (the Bayesian-bootstrap weights of the draw,
#'   positive and averaging 1, times the case weights of [nb_val()] if any; in
#'   [evsi()] also whole numbers for the pooled samples), e.g.
#'   `function(data, weights)
#'   coxph(Surv(time, status == 0) ~ x, data = data, weights = weights)`. Used
#'   by the value-of-information functions to refit the censoring model in
#'   every draw (see "Uncertainty"). Needed for a censoring model voipred
#'   cannot refit with case weights by itself.
#' @return An object of class `cipcw_control`.
#' @examples
#' library(survival)
#' rott <- rotterdam
#' rott$years <- rott$dtime / 365.25
#' df_dev_survival <- rott[1:1500, ]
#' df_val_survival <- rott[1501:2982, ]
#' cox <- coxph(Surv(years, death) ~ age + size + nodes + grade,
#'              data = df_dev_survival)
#' # the predicted risks of death by 5 years
#' sf <- survfit(cox, newdata = df_val_survival)
#' df_val_survival$p5 <- 1 - summary(sf, times = 5)$surv[1, ]
#'
#' # usually only the covariates of the censoring are needed (no cipcw_control())
#' nb_val(data = df_val_survival,
#'        predicted_risk = df_val_survival$p5,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3),
#'        cipcw = ~ age + nodes + grade)
#'
#' # the same covariates, truncated at 0.05 instead of 0.02
#' nb_val(data = df_val_survival,
#'        predicted_risk = df_val_survival$p5,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3),
#'        cipcw = cipcw_control(~ age + nodes + grade, min_prob = 0.05))
#'
#' # your own censoring model: fitted on the validation data, censoring as the
#' # event, here without truncation
#' cens_fit <- coxph(Surv(years, death == 0) ~ age + nodes + grade,
#'                   data = df_val_survival)
#' nb_val(data = df_val_survival,
#'        predicted_risk = df_val_survival$p5,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3),
#'        cipcw = cipcw_control(model = cens_fit, min_prob = NULL))
#'
#' # a flexible censoring model: a random survival forest (ranger), fitted on
#' # the validation data with the censored patients as the events
#' df_val_survival$censored <- 1 - df_val_survival$death
#' cens_forest <- ranger::ranger(Surv(years, censored) ~ age + size + nodes + grade,
#'                               data = df_val_survival, num.trees = 200, seed = 1)
#' res_forest <- nb_val(data = df_val_survival,
#'                      predicted_risk = df_val_survival$p5,
#'                      outcome = ~ Surv(years, death),
#'                      horizon = 5,
#'                      thresholds = c(0.1, 0.2, 0.3),
#'                      cipcw = cipcw_control(
#'                        model = cens_forest,
#'                        # how to refit it in each draw of the value-of-information
#'                        # functions, with the draw's case weights
#'                        refit_function = function(data, weights)
#'                          ranger::ranger(Surv(years, censored) ~ age + size + nodes + grade,
#'                                         data = data, case.weights = weights,
#'                                         num.trees = 200, seed = 1)))
#' res_forest
#' \donttest{
#' # the forest is refitted in every draw, so this takes a while
#' set.seed(1)
#' evpi(res_forest, n_sim = 20)
#' }
#' @export
cipcw_control <- function(formula = NULL, model = NULL, weight_event = NULL,
                          weight_horizon = NULL, min_prob = 0.02,
                          ties = c("left", "right"), predict_function = NULL,
                          refit_function = NULL) {
  ties <- match.arg(ties)
  if (!is.null(formula) && !inherits(formula, "formula")) {
    # a fitted censoring model given first, as in cipcw_control(cens_fit)
    if (is.null(model) && is.object(formula) && !is.data.frame(formula)) {
      model <- formula
      formula <- NULL
    } else {
      stop("`formula` must be a formula of the covariates of the censoring model, ",
           "e.g. ~ age + nodes.", call. = FALSE)
    }
  }
  if (!is.null(formula)) {
    if (!is.null(model) || !is.null(weight_event) || !is.null(weight_horizon))
      stop("Give one of `formula`, a censoring `model` or the weights, not several.",
           call. = FALSE)
    if (!is.null(predict_function) || !is.null(refit_function))
      stop("`predict_function` and `refit_function` are for your own censoring ",
           "`model`. With `formula`, voipred fits, predicts from and refits the Cox ",
           "censoring model itself.", call. = FALSE)
  }
  if (is.logical(refit_function))
    stop("There is no `refit` setting: the value-of-information functions always ",
         "refit the censoring model in each draw. `refit_function` is an optional ",
         "function(data, weights) returning the refitted censoring model.", call. = FALSE)
  if (!is.null(refit_function) && !is.function(refit_function))
    stop("`refit_function` must be a function(data, weights) returning the ",
         "refitted censoring model.", call. = FALSE)
  if (!is.null(refit_function) && is.null(model))
    stop("`refit_function` needs the censoring `model` it refits; given weights ",
         "cannot be refitted.", call. = FALSE)
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
    stop("`min_prob` must be a single number strictly between 0 and 1, or NULL ",
         "for no truncation.", call. = FALSE)
  if (is.null(formula) && is.null(model) && !has_w)
    stop("cipcw_control() needs the covariates of the censoring model (`formula`, ",
         "e.g. ~ age + nodes), a fitted censoring `model`, or the weights ",
         "(`weight_event` and `weight_horizon`).", call. = FALSE)
  structure(list(formula = formula, model = model, predict_function = predict_function,
                 weight_event = weight_event, weight_horizon = weight_horizon,
                 min_prob = min_prob, ties = ties, refit_function = refit_function),
            class = "cipcw_control")
}

#' @export
print.cipcw_control <- function(x, ...) {
  how <- if (!is.null(x$formula))
           paste0("Cox censoring model ", paste(deparse(x$formula), collapse = " "),
                  ", fitted by voipred")
         else if (!is.null(x$model)) paste0("censoring model (", class(x$model)[1],
                                       if (!is.null(x$predict_function))
                                         ", with predict_function",
                                       if (!is.null(x$refit_function))
                                         ", with refit_function", ")")
         else if (!is.null(x$weight_event)) "given weights" else "not set"
  cat("cIPCW settings: ", how, ", truncation: ",
      if (is.null(x$min_prob)) "off" else paste0("min_prob = ", x$min_prob),
      ", ties = ", x$ties, "\n", sep = "")
  invisible(x)
}

# `cipcw` of nb_val(): NULL (Kaplan-Meier), a formula of the censoring
# covariates, a fitted censoring model, or cipcw_control().
as_cipcw <- function(x) {
  if (is.null(x) || inherits(x, "cipcw_control")) return(x)
  if (inherits(x, "formula")) return(cipcw_control(formula = x))
  if (is.object(x) && !is.data.frame(x)) return(cipcw_control(model = x))
  stop("`cipcw` must be a formula of the covariates of the censoring model ",
       "(e.g. cipcw = ~ age + nodes), a fitted censoring model, or cipcw_control() ",
       "for more options. Leave it out (NULL) for Kaplan-Meier.", call. = FALSE)
}

# The Cox censoring model of cipcw = ~ covariates, fitted to `data` with
# censoring (only censoring; competing events are not censored) as the event,
# and with the case weights `w` if any. Returns the model and the data it was
# fitted to (`data` with the censoring outcome added, for the refits). A
# formula with a left-hand side is used as it is.
fit_censoring_cox <- function(formula, data, info, w = NULL) {
  rhs <- formula[[length(formula)]]
  vars <- all.vars(rhs)
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop("The censoring formula (`cipcw`) uses ", if (length(miss) == 1)
           "a variable that is" else "variables that are", " not in `data`: ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  na <- vars[vapply(vars, function(v) anyNA(data[[v]]), logical(1))]
  if (length(na))
    stop("The covariates of the censoring model have missing values (",
         paste(na, collapse = ", "), "), but every patient needs a censoring ",
         "weight. Impute them, or leave them out of `cipcw`.", call. = FALSE)
  d <- data
  f <- formula
  if (length(formula) == 2) {
    if (any(c(".voipred_time", ".voipred_censored") %in% names(d)))
      stop("`data` has a column named .voipred_time or .voipred_censored, which ",
           "voipred uses for the censoring model; rename it.", call. = FALSE)
    d$.voipred_time <- info$time
    d$.voipred_censored <- as.integer(info$status == 0)
    f <- stats::as.formula(call("~", quote(survival::Surv(.voipred_time, .voipred_censored)),
                                rhs), env = environment(formula))
  }
  env <- new.env(parent = if (is.environment(environment(formula))) environment(formula)
                          else globalenv())
  env$.voipred_data <- d
  # without x = TRUE: voipred's own prediction from coxph, faster in the refits
  cl <- as.call(c(list(quote(survival::coxph), formula = f, data = quote(.voipred_data)),
                  if (!is.null(w)) list(weights = w)))
  m <- tryCatch(eval(cl, env), error = function(e)
    stop("Could not fit the Cox censoring model ", paste(deparse(formula), collapse = " "),
         " (", conditionMessage(e), "). Give your own censoring model with ",
         "cipcw = cipcw_control(model = ).", call. = FALSE))
  list(model = m, data = d)
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
