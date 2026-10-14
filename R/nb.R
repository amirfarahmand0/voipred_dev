#' Net benefit of a prediction model in a validation sample
#'
#' Computes the net benefit of using the model (treat if predicted risk is at
#' or above the threshold), of treating everyone, and of treating no one, at
#' each risk threshold, for a binary, survival or competing-risks outcome.
#' This is the core quantity every other voipred validation metric (ΔNB,
#' P(useful), P(best), EVPI, EVSI) is built on.
#'
#' **Binary outcome.** With outcome \eqn{Y}, predicted risk \eqn{\hat\pi} and
#' threshold \eqn{z},
#' \deqn{NB_{model}(z) = \frac{1}{n}\sum_i Y_i T_i - \frac{1}{n}\sum_i (1-Y_i) T_i \frac{z}{1-z},\quad T_i = 1(\hat\pi_i \ge z)}
#' \deqn{NB_{all}(z) = \bar Y - (1-\bar Y)\frac{z}{1-z}, \qquad NB_{none}(z) = 0.}
#'
#' **Survival outcome** (`Surv(time, status)`, with a `horizon`). The event is
#' "the event by the horizon". With \eqn{q} the share treated at threshold
#' \eqn{z} and \eqn{F_+} the risk of the event by the horizon among them,
#' \eqn{NB_{model} = q F_+ - q (1 - F_+) z/(1-z)}. Two estimators:
#' * `censoring = "km"` (default): \eqn{F_+} and the overall risk are
#'   Kaplan-Meier estimates (Aalen-Johansen for competing risks), as in
#'   dcurves. If none of the treated patients is followed up to the horizon,
#'   the net benefit at that threshold is `NA`, with a warning.
#' * `censoring = "cipcw"`: conditional inverse probability of censoring
#'   weighting, which allows censoring to depend on covariates. Patients with
#'   the event by the horizon are weighted by \eqn{1/\hat G(T_i- \mid X_i)},
#'   patients event-free at the horizon by \eqn{1/\hat G(\tau \mid X_i)}, and
#'   patients censored before the horizon get no weight. The risk by the
#'   horizon is \eqn{\sum_{events} w_i / n}; sensitivity and specificity are
#'   weighted proportions, and \eqn{NB = \pi Se - (1-\pi)(1-Sp) z/(1-z)}. See
#'   [cipcw_control()].
#'
#' **Competing risks** (`Surv(time, status)` with a factor status whose first
#' level is censoring). The event is `cause` (default: the first event level)
#' by the horizon; a patient with a competing event by the horizon counts as
#' not having the event. The model must predict the cumulative incidence of
#' that cause; for competing risks only true censoring is censoring.
#'
#' A patient is treated when the predicted risk is at or above the threshold,
#' and a predicted risk of exactly 0 is never treated, as in dcurves.
#'
#' **Bayesian net benefit** (binary outcome, `bayesian = TRUE`). The prevalence
#' \eqn{\pi} and, at each threshold, the sensitivity \eqn{Se} and specificity
#' \eqn{Sp} of the model get Beta priors ([nb_uniform()] by default), updated
#' with the validation data to Beta posteriors (as in bayesDCA). The net
#' benefit columns are then the posterior means,
#' \eqn{E[\pi] E[Se] - (1 - E[\pi])(1 - E[Sp]) z/(1-z)} and
#' \eqn{E[\pi] - (1 - E[\pi]) z/(1-z)}.
#'
#' **Predicted risks.** Validation needs only the predicted risks of the
#' patients in `data`: give them with `predicted_risk`, a numeric vector with
#' one risk per row of `data`, e.g. `predicted_risk = val$p`. The name of the
#' column of `data` that holds them works too: `predicted_risk = "p"` is the
#' same as `predicted_risk = val$p` when `data = val`. The model
#' object is not needed, and the risks can come from anywhere (an R model,
#' published coefficients, other software). The outcome alone decides the type
#' of analysis (binary, survival or competing risks).
#'
#' Optionally, voipred can predict the risks itself from the fitted model,
#' given as `model` (or in place of `predicted_risk`, e.g. `nb(fit, val,
#' "y")`). It predicts them once, and every later step (the
#' value-of-information functions, plots) reuses them. Then
#' voipred keeps no list of models; there are two routes:
#' * **Automatically**, from the model alone: voipred calls [predict_risk()],
#'   which uses a method for the model's class if there is one (voipred ships
#'   some, and the author of a modelling package can add one), otherwise the
#'   standard routes, `predict()` with `type = "response"` or `"prob"` and
#'   `riskRegression::predictRisk()`. This covers the common models, each
#'   tested against its package's own predictions: e.g. `glm`, `gam`,
#'   random forests (`ranger`, `randomForest`), boosting (`gbm`, `xgboost`,
#'   `lightgbm`, `catboost`), `glmnet`, `svm`, bayespred; survival: `coxph`,
#'   `rms::cph`/`psm`, `survreg`, `flexsurvreg`, `rstpm2::stpm2`,
#'   bayespred.dev Bayesian Cox and parametric models (posterior mean of risk),
#'   survival forests (`ranger`, `randomForestSRC`, `aorsf`), survival trees;
#'   competing risks: `riskRegression::CSC` and `FGR`, `cmprsk::crr`,
#'   competing-risks `randomForestSRC`, survival's multi-state `coxph`. See
#'   [predict_risk()] for the full list. Any other model whose `predict()`
#'   returns probabilities in the usual way, or that has a
#'   `riskRegression::predictRisk()` method, works too.
#' * **With `predict_function`**, a function you write that returns the risks.
#'   It works for any model, and is the way to go when voipred stops with
#'   "Could not get predicted risks from this ... model". If given, it is
#'   always used.
#'
#' Whichever way they are given, the risks are checked: one per row, no `NA`,
#' values in \[0, 1\] (log-odds, relative risks or predicted times are
#' rejected as scores), and not only 0s and 1s (rejected as class labels).
#' With a model, where the model's own outcome can be read from its formula,
#' voipred also stops or
#' warns if it does not match `outcome` (e.g. a survival model with a binary
#' outcome, or a binary model with a survival outcome), and warns if
#' `horizon` is beyond the model's own follow-up.
#'
#' What voipred cannot check is the meaning of valid probabilities. With
#' `predict_function`, make sure it returns the risk of the event (not the
#' survival probability), by `horizon` and for `cause`.
#'
#' @param predicted_risk The predicted risks of the patients in `data`: a
#'   numeric vector, one probability per row of `data` (in its row order),
#'   e.g. `val$p`. Alternatively, the name of the column of `data` that holds
#'   them, e.g. `"p"`. For survival outcomes,
#'   the risk of the event by `horizon` (competing risks: of `cause` by
#'   `horizon`). This is all the validation needs. A fitted model can also be
#'   given here instead; it is then used as `model`.
#' @param data A data frame: the validation sample, one row per patient, with
#'   the outcome columns when `outcome` names them (and, with a `model`, every
#'   covariate it uses; with `censoring = "cipcw"`, every covariate of the
#'   censoring model).
#' @param outcome The observed outcome, one per row of `data`:
#'   * binary: a column name (`"y"`), a vector, or a one-sided formula
#'     (`~ y`); 0/1, logical (`TRUE` is the event) or a two-level factor (the
#'     event is level `"1"` if present, otherwise the second level);
#'   * survival: `~ Surv(time, status)` (or a `Surv` object, or the name of a
#'     `Surv` column), with status 1 = event and 0 = censored;
#'   * competing risks: `~ Surv(time, status)` where `status` is a **factor**
#'     whose first level is censoring and the other levels are the event
#'     types, e.g. `factor(code, 0:2, c("censored", "relapse", "death"))`.
#'     Numeric codes 0/1/2 are not accepted (Surv() would turn 2 into `NA`).
#'
#'   Formulas are evaluated in `data`. Only right-censored outcomes are
#'   supported; missing or negative times and missing statuses are errors.
#' @param thresholds Risk thresholds in \[0, 1), in any order (kept as
#'   given). At threshold `z` a patient is treated when the predicted risk is
#'   `>= z`; a risk of exactly 0 is never treated.
#' @param horizon Survival and competing-risks outcomes only (required
#'   there, an error for a binary outcome): the time by which the risk is
#'   predicted and the event counted, a single positive number in the units
#'   of the survival times. It must be within the follow-up of at least one
#'   patient.
#' @param cause Competing-risks outcomes only: the event type of interest,
#'   one of the event levels of the outcome (default: the first event level).
#'   The net benefit is for "this cause by the horizon"; patients with
#'   another cause by then count as not having the event. The cause is passed
#'   to the model by its label, or by its position among the event levels if
#'   the model numbers its causes; so a model that numbers its causes must
#'   number them in the order of the outcome's event levels.
#' @param censoring Survival and competing-risks outcomes only: how censoring
#'   is handled. `"km"` (default): Kaplan-Meier / Aalen-Johansen, which
#'   assumes censoring is independent of the outcome. `"cipcw"`: conditional
#'   inverse probability of censoring weights, which allows censoring to
#'   depend on covariates; it needs `cipcw`.
#' @param cipcw With `censoring = "cipcw"` only: a censoring model (with,
#'   optionally, your own `predict_function` for it) or the censoring
#'   weights, and optional truncation, made with [cipcw_control()]. Giving settings with `censoring = "km"` is an error.
#' @param weights Optional case weights, one non-negative number per row of
#'   `data` (with a positive sum), e.g. sampling weights or the weights of a
#'   (Bayesian) bootstrap draw. They weight each patient in the net benefit
#'   (an integer weight is the same as repeating the row; 0 drops the
#'   patient). They are not passed to the prediction model or the censoring
#'   model.
#' @param model Optional, instead of `predicted_risk`: a fitted prediction
#'   model, from which voipred predicts the risks of the patients in `data`,
#'   once. See "Predicted risks" for the models that work automatically and
#'   for `predict_function`. Give either `predicted_risk` or `model`.
#' @param predict_function Optional, with `model` only: a
#'   `function(model, newdata, horizon, cause)`
#'   that returns the predicted risks, used instead of the automatic routes
#'   (it always takes precedence). voipred calls it once with
#'   `newdata = data`; `horizon` and `cause` are those given to `nb()`
#'   (`NULL` when not used: `horizon` for a binary outcome, `cause` unless
#'   competing risks; `cause` is the event level's label). It must return,
#'   in the row order of `newdata`, the probability of the event (binary),
#'   of the event by `horizon` (survival), or of `cause` by `horizon`
#'   (competing risks): either a vector, one column per outcome class (the
#'   event column is taken: `"1"` or `".pred_1"`, otherwise the second of
#'   two), or a curve per patient, `list(time = <grid>, surv = <n x grid
#'   matrix>)` or `list(time = <grid>, cif = <n x grid matrix>)`, read at the
#'   last grid time `<= horizon`. A function that ignores some arguments can
#'   be written `function(model, newdata, ...)`.
#' @param bayesian Binary outcomes only: `TRUE` for the Bayesian net benefit
#'   (see "Bayesian net benefit"): the net benefit columns are posterior means,
#'   and the value-of-information functions ([evpi()], [evsi()], [p_useful()],
#'   [p_best()]) draw from the posterior instead of resampling the patients.
#'   Default `FALSE`. Case weights cannot be used with it.
#' @param prior With `bayesian = TRUE`: the prior, [nb_uniform()] (default,
#'   Beta(1, 1) for the prevalence, sensitivity and specificity, as in
#'   bayesDCA), [nb_jeffreys()] or your own Beta priors with [nb_beta()].
#' @param label Name of the model in printed output and plots, a single
#'   string. Default: the name of the object given as `predicted_risk` or
#'   `model` (e.g. `p`, the column name, or `fit`), or `"Model"` for an
#'   expression. `"Treat All"` and `"Treat None"` are reserved.
#' @return A data frame of class `voipred_nb` with one row per threshold
#'   and columns `threshold`, `NB_model` (using the model; `NA` where it
#'   cannot be estimated, see Details), `NB_all` (treat everyone), `NB_none`
#'   (treat no one, always 0) and `delta_NB`, the difference in net benefit
#'   between the model and the best default strategy,
#'   \eqn{NB_{model} - \max(NB_{all}, NB_{none})}: positive where the model
#'   is useful. Plot the results with [plot_dca()]. Attributes: `label`, `n` (rows of `data`),
#'   `prevalence` (binary: the weighted event rate; survival: the estimated
#'   risk of the event, or of `cause`, by the horizon), `outcome_type`
#'   (`"binary"`, `"survival"` or `"competing"`), `risk` (the predicted
#'   risks, one per row of `data`), `bayesian` (and with `bayesian = TRUE` the
#'   `prior`; `prevalence` is then its posterior mean) and, for survival
#'   outcomes, `horizon`, `cause` (the event level used) and `censoring`.
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:10000, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#'
#' # The predicted risks are all that is needed
#' val$p <- predict(fit, val, type = "response")
#' res <- nb(predicted_risk = val$p, data = val, outcome = "day30",
#'           thresholds = seq(0, 0.3, by = 0.01))
#' head(res)
#'
#' # The same, with the name of the column of `data` that holds the risks
#' nb(predicted_risk = "p", data = val, outcome = "day30", thresholds = c(0.05, 0.1))
#'
#' # Or let voipred predict them from the model
#' nb(fit, val, "day30", thresholds = c(0.05, 0.1))
#'
#' # Bayesian net benefit (posterior means), as in bayesDCA
#' nb(predicted_risk = val$p, data = val, outcome = "day30",
#'    thresholds = c(0.05, 0.1), bayesian = TRUE)
#' nb(predicted_risk = val$p, data = val, outcome = "day30",
#'    thresholds = c(0.05, 0.1), bayesian = TRUE, prior = nb_jeffreys())
#'
#' # Any model, with your own prediction function
#' nb(model = fit, data = val, outcome = "day30", thresholds = c(0.05, 0.1),
#'    predict_function = function(model, newdata, ...)
#'      predict(model, newdata, type = "response"))
#'
#' # Survival: death by 5 years (Rotterdam breast cancer data)
#' library(survival)
#' rott <- rotterdam
#' rott$years <- rott$dtime / 365.25
#' cox <- coxph(Surv(years, death) ~ age + size + nodes + grade, data = rott[1:1500, ],
#'              x = TRUE)
#' val_s <- rott[1501:2982, ]
#' nb(cox, val_s, ~ Surv(years, death), horizon = 5, thresholds = c(0.1, 0.2, 0.3))
#'
#' # The same with conditional IPCW: censoring modelled on the validation data
#' cens <- coxph(Surv(years, death == 0) ~ age + nodes + grade, data = val_s)
#' nb(cox, val_s, ~ Surv(years, death), horizon = 5, thresholds = c(0.1, 0.2, 0.3),
#'    censoring = "cipcw", cipcw = cipcw_control(model = cens))
#' @export
nb <- function(predicted_risk = NULL, data, outcome,
               thresholds = seq(0, 0.99, by = 0.01), horizon = NULL, cause = NULL,
               censoring = c("km", "cipcw"), cipcw = cipcw_control(), weights = NULL,
               model = NULL, predict_function = NULL, bayesian = FALSE,
               prior = nb_uniform(), label = NULL) {
  # a fitted model given in place of the risks (e.g. nb(fit, val, "y")) is the model
  risk_expr <- substitute(predicted_risk)
  model_expr <- substitute(model)
  if (!is.null(predicted_risk) && !is_risk_input(predicted_risk)) {
    if (!is.null(model))
      stop("Give either `predicted_risk` or `model`, not both.", call. = FALSE)
    model <- predicted_risk
    model_expr <- risk_expr
    predicted_risk <- NULL
  }
  if (is.null(label))
    label <- if (!is.null(predicted_risk)) risk_label(risk_expr, predicted_risk)
             else if (is.name(model_expr)) as.character(model_expr)
  label <- check_label(label, NULL)
  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  censoring <- match.arg(censoring)
  if (!inherits(cipcw, "cipcw_control"))
    stop("`cipcw` must be made with cipcw_control().", call. = FALSE)
  info <- parse_outcome(outcome, data)
  k <- check_survival_args(info, horizon, cause, censoring, cipcw)
  info$k <- k
  thresholds <- check_thresholds(thresholds)
  n <- nrow(data)
  w <- check_weights(weights, n)
  check_bayesian(bayesian, prior, info, weights)

  # the risks: given, or predicted once from the model
  p <- if (!is.null(predicted_risk)) {
    if (!is.null(model))
      stop("Give either `predicted_risk` or `model`, not both.", call. = FALSE)
    if (!is.null(predict_function))
      stop("`predict_function` predicts from `model`; it is not used with ",
           "`predicted_risk`.", call. = FALSE)
    given_risk(predicted_risk, data)
  } else {
    if (is.null(model))
      stop("Give the predicted risks with `predicted_risk` (or a fitted `model`).",
           call. = FALSE)
    get_risk(model, data, info, horizon, predict_function)
  }
  # what the VoI functions need to recompute the net benefit with other weights
  voi <- list(type = info$type, thresholds = thresholds, risk = p, weights = w)
  if (info$type == "binary") {
    res <- if (bayesian) bayes_nb_mean(info$y, p, thresholds, prior)
           else nb_core(info$y, p, thresholds, w)
    voi$y <- info$y
    if (bayesian) voi$prior <- prior
  } else {
    voi[c("time", "status", "k", "horizon", "censoring")] <-
      list(info$time, info$status, k, horizon, censoring)
    if (censoring == "km") {
      res <- nb_core_km(info$time, info$status, k, p, thresholds, horizon, w)
    } else {
      cw <- cipcw_weights(cipcw, data, info, horizon, k)
      res <- nb_core_cipcw(info$time, info$status, k, p, thresholds, horizon,
                           cw$weight_event, cw$weight_horizon, w)
      voi[c("weight_event", "weight_horizon")] <- list(cw$weight_event, cw$weight_horizon)
      # what refitting the censoring model in each draw needs
      voi[c("cipcw", "data", "info")] <- list(cipcw, data, info)
    }
  }

  out <- data.frame(threshold = thresholds,
                    NB_model = res$NB_model,
                    NB_all = res$NB_all,
                    NB_none = 0)
  out$delta_NB <- out$NB_model - pmax(out$NB_all, 0)
  attr(out, "n") <- n
  attr(out, "prevalence") <- res$prevalence
  attr(out, "outcome_type") <- info$type
  attr(out, "bayesian") <- bayesian
  if (bayesian) attr(out, "prior") <- prior
  if (info$type != "binary") {
    attr(out, "horizon") <- horizon
    attr(out, "cause") <- info$states[k]
    attr(out, "censoring") <- censoring
  }
  attr(out, "risk") <- p
  attr(out, "label") <- label
  attr(out, "voi") <- voi
  class(out) <- c("voipred_nb", "data.frame")
  out
}

#' @export
print.voipred_nb <- function(x, ...) {
  type <- attr(x, "outcome_type")
  if (is.null(type)) type <- "binary"
  if (type == "binary") {
    cat("Net benefit of ", attr(x, "label"),
        if (isTRUE(attr(x, "bayesian"))) " (Bayesian, posterior means)",
        ": n = ", attr(x, "n"), ", prevalence = ",
        format(attr(x, "prevalence"), digits = 3), ", ", nrow(x),
        " thresholds\n", sep = "")
  } else {
    cat("Net benefit of ", attr(x, "label"), " (", type, if (type == "competing")
          paste0(", cause = ", attr(x, "cause")),
        ", horizon = ", attr(x, "horizon"), ", ", attr(x, "censoring"),
        "): n = ", attr(x, "n"), ", risk by horizon = ",
        format(attr(x, "prevalence"), digits = 3), ", ", nrow(x),
        " thresholds\n", sep = "")
  }
  y <- x
  attr(y, "risk") <- NULL
  print(as.data.frame(unclass(y))[, c("threshold", "NB_model", "NB_all", "NB_none", "delta_NB")],
        row.names = FALSE, ...)
  invisible(x)
}


# ---- internal: predicted risks ----------------------------------------------

# predict_function if given, otherwise predict_risk() (a method for the class, or
# the standard routes). Always read and checked the same way.
get_risk <- function(model, data, info, horizon, predict_function) {
  check_model_object(model)
  cause <- if (info$type == "competing") info$states[info$k] else NULL
  if (!is.null(predict_function)) {
    if (!is.function(predict_function))
      stop("`predict_function` must be a function(model, newdata, horizon, cause).",
           call. = FALSE)
    out <- predict_function(model, data, horizon = horizon, cause = cause)
    return(finish_risk(out, data, horizon))
  }
  check_outcome_match(model, data, info, horizon)
  out <- if (info$type == "competing")
    predict_risk(model, data, horizon = horizon, cause = cause, cause_index = info$k)
  else predict_risk(model, data, horizon = horizon)
  finish_risk(out, data, horizon)
}


# ---- internal: binary -------------------------------------------------------

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

check_thresholds <- function(thresholds) {
  if (!is.numeric(thresholds) || length(thresholds) == 0 || anyNA(thresholds))
    stop("`thresholds` must be a non-empty numeric vector.", call. = FALSE)
  if (any(thresholds < 0 | thresholds >= 1))
    stop("`thresholds` must lie in [0, 1).", call. = FALSE)
  thresholds
}

check_label <- function(label, model_expr) {
  if (is.null(label))
    label <- if (is.name(model_expr)) as.character(model_expr) else "Model"
  if (!is.character(label) || length(label) != 1 || is.na(label) || !nzchar(label))
    stop("`label` must be a single non-empty string.", call. = FALSE)
  if (tolower(label) %in% c("treat all", "treat none"))
    stop("`label` cannot be \"", label, "\": \"Treat All\" and \"Treat None\" ",
         "are the default strategies.", call. = FALSE)
  label
}

check_weights <- function(weights, n) {
  if (is.null(weights)) return(NULL)
  if (!is.numeric(weights) || length(weights) != n || anyNA(weights) ||
      any(weights < 0) || sum(weights) <= 0)
    stop("`weights` must be ", n, " non-negative numbers with a positive sum.",
         call. = FALSE)
  weights
}


# ---- internal: outcome ------------------------------------------------------
#
# Outcome parsing: binary, survival (right-censored) or competing risks.
#
# `outcome` can be a column name, a vector / Surv object with one value per
# row, or a one-sided formula evaluated in `data` (e.g. ~ Surv(time, status)).
# A Surv object with 0/1 status is a survival outcome; a Surv object built
# from a factor status (first level = censored) is a competing-risks outcome.

parse_outcome <- function(outcome, data) {
  n <- nrow(data)
  if (inherits(outcome, "formula")) {
    if (length(outcome) != 2)
      stop("An `outcome` formula must be one-sided, e.g. ~ Surv(time, status).",
           call. = FALSE)
    y <- tryCatch(eval(outcome[[2]], data, environment(outcome)),
                  error = function(e)
                    stop("Could not evaluate the `outcome` formula in `data`: ",
                         conditionMessage(e), call. = FALSE))
  } else if (is.character(outcome) && length(outcome) == 1) {
    if (!outcome %in% names(data))
      stop("Column '", outcome, "' is not in `data`.", call. = FALSE)
    y <- data[[outcome]]
  } else {
    y <- outcome
  }
  if (inherits(y, "Surv")) return(parse_surv(y, n))
  list(type = "binary", y = binary_outcome(y, n))
}

binary_outcome <- function(y, n) {
  if (length(y) != n)
    stop("`outcome` has length ", length(y), " but `data` has ", n,
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
      stop("`outcome` must be binary (0/1, logical, or a two-level factor) ",
           "or a survival outcome made with Surv().", call. = FALSE)
  }
  y
}

parse_surv <- function(y, n) {
  type <- attr(y, "type")
  if (!type %in% c("right", "mright"))
    stop("Only right-censored survival outcomes are supported: ",
         "Surv(time, status). This outcome is of type '", type, "'.",
         call. = FALSE)
  if (nrow(y) != n)
    stop("`outcome` has ", nrow(y), " rows but `data` has ", n, " rows.",
         call. = FALSE)
  time <- as.numeric(y[, "time"])
  status <- as.integer(y[, "status"])
  if (anyNA(time) || anyNA(status))
    stop("The survival outcome has missing times or statuses. (For competing ",
         "risks, give the status as a factor whose first level is censoring, ",
         "e.g. Surv(time, factor(status, 0:2)); Surv() turns other numeric ",
         "codes into NA.)", call. = FALSE)
  if (any(time < 0))
    stop("The survival outcome has negative times.", call. = FALSE)
  if (type == "right") {
    return(list(type = "survival", time = time, status = status,
                states = "event"))
  }
  states <- attr(y, "states")       # event levels; status 0 = censored
  list(type = "competing", time = time, status = status, states = states)
}

# Check horizon / cause / censoring arguments against the outcome type and
# return the index of the cause of interest (1 for a standard survival outcome).
check_survival_args <- function(info, horizon, cause, censoring, cipcw) {
  if (info$type == "binary") {
    if (!is.null(horizon))
      stop("`horizon` is only used with a survival outcome made with Surv().",
           call. = FALSE)
    if (!is.null(cause))
      stop("`cause` is only used with a competing-risks outcome.", call. = FALSE)
    if (censoring != "km" || !is_default_cipcw(cipcw))
      stop("`censoring` and `cipcw` are only used with a survival outcome.",
           call. = FALSE)
    return(invisible(NULL))
  }
  if (is.null(horizon))
    stop("A survival outcome needs `horizon`: the time at which the risk is ",
         "predicted (in the same units as the survival times).", call. = FALSE)
  if (!is.numeric(horizon) || length(horizon) != 1 || is.na(horizon) ||
      !is.finite(horizon) || horizon <= 0)
    stop("`horizon` must be a single positive number.", call. = FALSE)
  if (horizon > max(info$time))
    stop("`horizon` (", horizon, ") is beyond the follow-up of every patient ",
         "in `data` (longest follow-up: ", signif(max(info$time), 4), ").",
         call. = FALSE)
  if (info$type == "survival") {
    if (!is.null(cause))
      stop("`cause` is only used with a competing-risks outcome, i.e. ",
           "Surv(time, status) with a factor status.", call. = FALSE)
    k <- 1L
  } else {
    if (is.null(cause)) {
      k <- 1L
    } else {
      if (length(cause) != 1 || !as.character(cause) %in% info$states)
        stop("`cause` must be one of the event types of the outcome: ",
             paste(info$states, collapse = ", "), ".", call. = FALSE)
      k <- match(as.character(cause), info$states)
    }
  }
  if (censoring == "km" && !is_default_cipcw(cipcw))
    stop("`cipcw` settings were given but `censoring = \"km\"`. Use ",
         "`censoring = \"cipcw\"` to apply them, or remove `cipcw`.",
         call. = FALSE)
  if (censoring == "cipcw" && is.null(cipcw$model) && is.null(cipcw$weight_event))
    stop("`censoring = \"cipcw\"` needs a censoring model or weights: ",
         "cipcw = cipcw_control(model = ...) or ",
         "cipcw_control(weight_event = ..., weight_horizon = ...).", call. = FALSE)
  k
}


# ---- internal: survival estimators ------------------------------------------
#
# Net benefit for survival and competing-risks outcomes at a horizon.
#
# Both estimators treat a patient when the predicted risk is >= the threshold;
# a predicted risk of exactly 0 is never treated (the same rules as binary and
# as dcurves). `k` is the index of the cause of interest (1 for standard
# survival; status codes: 0 = censored, 1..K = event types).

# ---- Kaplan-Meier / Aalen-Johansen (dcurves' estimator) -------------------

nb_core_km <- function(time, status, k, p, thresholds, horizon, w = NULL) {
  if (is.null(w)) w <- rep(1, length(time))
  odds <- thresholds / (1 - thresholds)
  # compiled loop over the thresholds (src/km_nb.cpp): the same arithmetic,
  # in the same order, as km_cuminc() applied to the treated patients
  o <- order(time)
  r <- km_nb_sorted(as.numeric(time[o]), as.integer(status[o]), as.integer(k),
                    as.numeric(p[o]), as.numeric(w[o]), as.numeric(thresholds),
                    as.numeric(horizon))
  undefined <- r$undefined
  if (any(undefined))
    warning("No patient treated at ", if (sum(undefined) == 1) "threshold " else
              "thresholds ", format_thresholds(thresholds[undefined]),
            " is followed up to the horizon (", horizon, "), so the risk among ",
            "them cannot be estimated: NB_model is NA there.", call. = FALSE)
  list(NB_model = r$NB_model,
       NB_all = r$risk_all - (1 - r$risk_all) * odds,
       prevalence = r$risk_all)
}

# Kaplan-Meier / Aalen-Johansen cumulative incidence of cause k by tau
# (the R version of the compiled code used by nb_core_km()).
km_cuminc <- function(time, status, k, w, tau) {
  keep <- w > 0
  time <- time[keep]; status <- status[keep]; w <- w[keep]
  o <- order(time)
  time <- time[o]; status <- status[o]; w <- w[o]
  ut <- unique(time)
  g <- match(time, ut)
  d_all <- rowsum(w * (status != 0), g, reorder = FALSE)[, 1]
  d_k   <- rowsum(w * (status == k), g, reorder = FALSE)[, 1]
  at_g  <- rowsum(w, g, reorder = FALSE)[, 1]
  n_risk <- rev(cumsum(rev(at_g)))
  s <- cumprod(1 - d_all / n_risk)                 # all-cause survival after each time
  s_before <- c(1, s[-length(s)])                  # survival just before each time
  inc <- cumsum(s_before * d_k / n_risk)           # cumulative incidence of cause k
  i <- findInterval(tau, ut)
  if (i == 0) 0 else unname(inc[i])
}

format_thresholds <- function(z) {
  if (length(z) <= 6) return(paste(signif(z, 3), collapse = ", "))
  paste0(signif(min(z), 3), " to ", signif(max(z), 3), " (", length(z), " thresholds)")
}

# ---- Conditional IPCW (as in the dcurves cIPCW branch) --------------------
#
# Patients with the event of interest by the horizon get weight_event =
# 1 / G(T_i- | X_i); event-free at the horizon get weight_horizon =
# 1 / G(tau | X_i); a competing event by the horizon counts as a negative with
# weight_event. Patients censored before the horizon get no weight.
# prevalence = sum(event weights) / n; Se and Sp are weighted proportions.

nb_core_cipcw <- function(time, status, k, p, thresholds, horizon,
                          weight_event, weight_horizon, w = NULL) {
  if (is.null(w)) w <- rep(1, length(time))
  event_k <- status == k & time <= horizon
  competing <- status != 0 & status != k & time <= horizon
  event_free <- time > horizon
  we <- ifelse(event_k | competing, weight_event, 0)
  wh <- ifelse(event_free, weight_horizon, 0)
  pos_w <- w * event_k * we                        # weighted events of interest
  neg_w <- w * (event_free * wh + competing * we)  # weighted negatives
  den_se <- sum(pos_w)
  den_sp <- sum(neg_w)
  prev <- den_se / sum(w)

  # weighted sums over the treated (risk >= z, risk 0 never treated), all z at once
  p2 <- p
  p2[p2 == 0] <- -1
  o <- order(p2)
  ps <- p2[o]
  cum_pos <- c(0, cumsum(pos_w[o]))
  cum_neg <- c(0, cumsum(neg_w[o]))
  m <- findInterval(thresholds, ps, left.open = TRUE)
  tp_w <- cum_pos[length(cum_pos)] - cum_pos[m + 1]  # treated events
  fp_w <- cum_neg[length(cum_neg)] - cum_neg[m + 1]  # treated negatives
  se <- if (den_se == 0) rep(0, length(thresholds)) else tp_w / den_se
  sp <- if (den_sp == 0) rep(1, length(thresholds)) else 1 - fp_w / den_sp

  odds <- thresholds / (1 - thresholds)
  list(NB_model = prev * se - (1 - prev) * (1 - sp) * odds,
       NB_all = prev - (1 - prev) * odds,
       prevalence = prev)
}

# ---- internal: predicted risks given directly -------------------------------

# Risks (a vector, or a column name) rather than a fitted model.
is_risk_input <- function(x) is.atomic(x)

# Predicted risks given with `predicted_risk`: a vector, or the name of a
# column of `data`.
given_risk <- function(risk, data) {
  n <- nrow(data)
  if (is.character(risk) && length(risk) == 1) {
    if (!risk %in% names(data))
      stop("`predicted_risk` names a column (\"", risk, "\") that is not in `data`.",
           call. = FALSE)
    risk <- data[[risk]]
  }
  if (!is.numeric(risk) || (is.matrix(risk) && ncol(risk) != 1))
    stop("`predicted_risk` must be numeric: one predicted probability per row of ",
         "`data`, or the name of the column that holds them.", call. = FALSE)
  risk <- as.vector(risk)
  if (length(risk) != n)
    stop("`predicted_risk` has ", length(risk), " values but `data` has ", n, " rows.",
         call. = FALSE)
  if (anyNA(risk))
    stop(sum(is.na(risk)), " of the risks in `predicted_risk` are NA.", call. = FALSE)
  if (any(risk < 0 | risk > 1))
    stop("`predicted_risk` has values outside [0, 1] (from ", signif(min(risk), 3),
         " to ", signif(max(risk), 3), "), so these are not probabilities (e.g. ",
         "log-odds must be turned into probabilities first).", call. = FALSE)
  if (all(risk == 0 | risk == 1))
    stop("`predicted_risk` has only 0s and 1s. These look like predicted classes, ",
         "not probabilities.", call. = FALSE)
  unname(risk)
}

# Default label for given risks: the column name, or the name of the vector
# (`p`, or `p` in `val$p`).
risk_label <- function(expr, risk) {
  if (is.character(risk) && length(risk) == 1) return(risk)
  if (is.name(expr)) return(as.character(expr))
  if (is.call(expr) && identical(expr[[1]], as.name("$"))) return(as.character(expr[[3]]))
  "Model"
}

check_bayesian <- function(bayesian, prior, info, weights) {
  if (!is.logical(bayesian) || length(bayesian) != 1 || is.na(bayesian))
    stop("`bayesian` must be TRUE or FALSE.", call. = FALSE)
  if (!inherits(prior, "nb_prior"))
    stop("`prior` must be nb_uniform(), nb_jeffreys() or nb_beta().", call. = FALSE)
  if (!bayesian) {
    if (!identical(prior, nb_uniform()))
      stop("`prior` is only used with `bayesian = TRUE`.", call. = FALSE)
    return(invisible())
  }
  if (info$type != "binary")
    stop("The Bayesian net benefit (`bayesian = TRUE`) is available for binary ",
         "outcomes only.", call. = FALSE)
  if (!is.null(weights))
    stop("Case weights cannot be used with `bayesian = TRUE`.", call. = FALSE)
  invisible()
}
