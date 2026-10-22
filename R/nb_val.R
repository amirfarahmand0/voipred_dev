#' Net benefit of prediction models in a validation sample
#'
#' Computes the net benefit of using a model (treat if predicted risk is at
#' or above the threshold), of treating everyone, and of treating no one, at
#' each risk threshold, for a binary, survival or competing-risks outcome,
#' for one model or several models validated on the same patients.
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
#' * Kaplan-Meier (default, `cipcw` not given): \eqn{F_+} and the overall risk are
#'   Kaplan-Meier estimates (Aalen-Johansen for competing risks), as in
#'   dcurves. If none of the treated patients is followed up to the horizon,
#'   the net benefit at that threshold is `NA`, with a warning.
#' * cIPCW, when `cipcw` is given: conditional inverse probability of censoring
#'   weighting, which allows censoring to depend on covariates. Patients with
#'   the event by the horizon are weighted by \eqn{1/\hat G(T_i- \mid X_i)},
#'   patients event-free at the horizon by \eqn{1/\hat G(\tau \mid X_i)}, and
#'   patients censored before the horizon get no weight. The risk by the
#'   horizon is \eqn{\sum_{events} w_i / n}; sensitivity and specificity are
#'   weighted proportions, and \eqn{NB = \pi Se - (1-\pi)(1-Sp) z/(1-z)}.
#'   Usually `cipcw` is just the covariates of the censoring, e.g.
#'   `cipcw = ~ age + nodes`: voipred fits a Cox model of the censoring times
#'   of `data` on them and truncates the probabilities of remaining uncensored
#'   at 0.02 (weights at 50). See [cipcw_control()] for your own censoring
#'   model, given weights, or other truncation.
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
#' **Prior.** `nb_val()` works in a Bayesian framework, which the
#' value-of-information functions ([p_useful()], [p_best()], [evpi()],
#' [evsi()]) use to represent the uncertainty of the net benefit. By default
#' it uses no prior information: with the improper prior [nb_improper()]
#' (Beta(0, 0)) only the data are used, and the net benefit is exactly the
#' usual one, \eqn{TP/n - FP/n \cdot z/(1-z)}. For a binary outcome, prior
#' information can be added with another prior ([nb_uniform()],
#' [nb_jeffreys()], your own with [nb_beta()], or bayesDCA's threshold-varying
#' informative prior [nb_informative()]; see [nb_priors]): the
#' prevalence \eqn{\pi} and, at each threshold, the sensitivity \eqn{Se}
#' and specificity \eqn{Sp} of the model then have Beta priors, updated with
#' the validation data (as in bayesDCA), and the net benefit columns are the
#' posterior means, \eqn{E[\pi] E[Se] - (1 - E[\pi])(1 - E[Sp]) z/(1-z)} and
#' \eqn{E[\pi] - (1 - E[\pi]) z/(1-z)}. Survival and competing-risks
#' outcomes use only the improper prior: their net benefit is the usual
#' Kaplan-Meier, Aalen-Johansen or cIPCW estimate.
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
#' **Several models.** To compare models, give them in one call as a list,
#' each element one model: a vector of risks, the name of a column of `data`,
#' or a fitted model, mixed freely, e.g.
#' `predicted_risk = list(Full = fit1, Age = "p_age", Old = val$p_old)`. All
#' are predicted for the same patients (the rows of `data`), and every other
#' argument (outcome, thresholds, horizon, cause, `cipcw`, prior, weights)
#' applies to all of them. (`c()` cannot hold several models: it joins vectors
#' into one and takes fitted models apart, so voipred asks for `list()`.) The
#' result is the same as computing each model separately and combining the
#' results with [c()][c.voipred_nb], and the value-of-information functions
#' draw the models jointly, so their correlation on the same patients is kept.
#' With cIPCW, one censoring model (and one set of censoring weights) is used
#' for all the models. To compare a different set of models, call `nb_val()`
#' with those models.
#'
#' **Labels.** Each model is labelled, in printed output and plots, by the
#' first of:
#' 1. `label`, if given (one per model);
#' 2. its name in the list, e.g. `list(Full = ...)`;
#' 3. the column name, when it is given as one, e.g. `"p_glm"`;
#' 4. the name in the call: of the object (`p_glm`, `fit_glm`), or after `$`
#'    (`val$p_glm` gives `p_glm`);
#' 5. otherwise `"Model"` (one model) or `"Model 1"`, `"Model 2"`, ... (by
#'    position).
#'
#' The labels must be different, and cannot be `"Treat All"` or
#' `"Treat None"` (the default strategies); otherwise `nb_val()` stops and asks
#' for names.
#'
#' **Predicting the risks from a model (optional).** Instead of the risks,
#' the fitted model can be given as `predicted_risk` (e.g.
#' `nb_val(val, fit, "y")`, or as an element of a list). voipred then predicts the risks of the patients
#' in `data` from it, once, and every later step (the value-of-information
#' functions, plots) reuses them: the probability of the event (binary), of
#' the event by `horizon` (survival), or of `cause` by `horizon` (competing
#' risks). There are two ways:
#' * **By default**, from the model alone. voipred keeps no list of models: it
#'   uses `predict()` with `type = "response"` or `type = "prob"`, or
#'   `riskRegression::predictRisk()` (at `horizon`, and for `cause`), and its
#'   own prediction code only for models these get wrong. Checked in the
#'   tests against each package's own predictions:
#'   * binary: `glm` (logit, probit, cloglog), `logistf`, `mgcv::gam`,
#'     `earth`, `rpart`, `partykit::ctree`, `randomForest`, `ranger`
#'     (`probability = TRUE`), `gbm`, `xgboost`, `lightgbm`, `catboost`,
#'     `glmnet`, `cv.glmnet`, `e1071::svm` (`probability = TRUE`),
#'     `e1071::naiveBayes`, `naivebayes`, `nnet`, `MASS::lda`, `MASS::qda`,
#'     `kknn`, and bayespred `bpmfit` and `bpm` (posterior mean of risk);
#'   * survival: `coxph`, `rms::cph`, `survreg`, `rms::psm`, `flexsurvreg`,
#'     `flexsurvspline`, `rstpm2::stpm2`, `mgcv::gam` with `cox.ph()`,
#'     `ranger`, `randomForestSRC`, `aorsf`, `partykit::ctree`, `prodlim`, and
#'     bayespred.dev Bayesian Cox and parametric survival models (posterior
#'     mean of risk);
#'   * competing risks: `riskRegression::CSC`, `riskRegression::FGR` (fitted
#'     for the requested cause), `cmprsk::crr`, competing-risks
#'     `randomForestSRC`, and survival's multi-state `coxph`.
#'
#'   Any other model whose `predict()` returns probabilities in the usual way,
#'   or that has a `riskRegression::predictRisk()` method, works too. Models
#'   that only give relative risks or times (e.g. `gbm` and `glmnet` Cox
#'   models, `rpart` survival trees) have no route to probabilities. Before
#'   predicting, the covariates the model uses are compared with the columns
#'   of `data`, and missing ones are named in an error.
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
#' @param data A data frame: the validation sample, one row per patient, with
#'   the outcome columns when `outcome` names them (and, for a fitted model, every
#'   covariate it uses; with cIPCW, every covariate of the censoring model).
#' @param predicted_risk The predicted risks of the patients in `data`: a
#'   numeric vector, one probability per row of `data` (in its row order),
#'   e.g. `val$p`. Alternatively, the name of the column of `data` that holds
#'   them, e.g. `"p"`. For survival outcomes,
#'   the risk of the event by `horizon` (competing risks: of `cause` by
#'   `horizon`). This is all the validation needs. A fitted model can also be
#'   given instead: voipred predicts the risks from it (see "Predicting the
#'   risks from a model"). For several models: a list, each element a vector
#'   of risks, a column name or a fitted model, e.g.
#'   `list(Full = fit1, Age = "p_age")` (see "Several models").
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
#' @param cipcw Survival and competing-risks outcomes only: how censoring is
#'   handled. `NULL` (default): Kaplan-Meier / Aalen-Johansen, which assumes
#'   censoring is independent of the outcome. To let censoring depend on
#'   covariates (conditional IPCW), give:
#'   * the covariates of the censoring, as a one-sided formula, e.g.
#'     `cipcw = ~ age + nodes + grade`: voipred fits
#'     `coxph(Surv(time, censored) ~ age + nodes + grade)` to the patients of
#'     `data` (with competing risks, only true censoring is censoring; with
#'     `weights`, weighted by them), predicts from it the probabilities of
#'     remaining uncensored, truncated at 0.02, and refits it in every draw of
#'     the value-of-information functions;
#'   * or your own censoring model fitted to the censoring times of `data`,
#'     e.g. `cipcw = coxph(Surv(time, status == 0) ~ age, data = val)`;
#'   * or [cipcw_control()], for more options: another truncation
#'     (`cipcw_control(~ age + nodes, min_prob = 0.05)`, or `min_prob = NULL`
#'     for none), given weights, `ties`, or your own predict and refit
#'     functions for your censoring model.
#' @param weights Optional case weights, one non-negative number per row of
#'   `data` (with a positive sum), e.g. sampling weights. They weight each
#'   patient in the net benefit (an integer weight is the same as repeating the
#'   row; 0 drops the patient), and the value-of-information functions then
#'   use the Bayesian bootstrap of the weighted sample. The Cox censoring
#'   model of `cipcw = ~ covariates` is fitted with them; they are not passed
#'   to the prediction model, nor to a censoring model you fitted yourself.
#'   Only with the improper prior.
#' @param prior The prior (see "Prior" and [nb_priors]). Default
#'   [nb_improper()]: no prior information, the data only, which gives the
#'   usual net benefit. For a binary outcome, prior information can be added
#'   with [nb_uniform()], [nb_jeffreys()], your own Beta priors with
#'   [nb_beta()], or [nb_informative()] (bayesDCA's informative prior, whose
#'   sensitivity and specificity priors change with the threshold; for small
#'   effective sample sizes). Survival and competing-risks outcomes: only
#'   [nb_improper()].
#' @param predict_function Optional, for fitted models only: a
#'   `function(model, newdata, horizon, cause)`
#'   that returns the predicted risks, used instead of the automatic routes
#'   (it always takes precedence). voipred calls it once with
#'   `newdata = data`; `horizon` and `cause` are those given to `nb_val()`
#'   (`NULL` when not used: `horizon` for a binary outcome, `cause` unless
#'   competing risks; `cause` is the event level's label). It must return,
#'   in the row order of `newdata`, the probability of the event (binary),
#'   of the event by `horizon` (survival), or of `cause` by `horizon`
#'   (competing risks): either a vector, one column per outcome class (the
#'   event column is taken: `"1"` or `".pred_1"`, otherwise the second of
#'   two), or a curve per patient, `list(time = <grid>, surv = <n x grid
#'   matrix>)` or `list(time = <grid>, cif = <n x grid matrix>)`, read at the
#'   last grid time `<= horizon`. A function that ignores some arguments can
#'   be written `function(model, newdata, ...)`. With several fitted models,
#'   one function is used for each of them; or give a list of functions named
#'   by the labels of the models they are for, e.g.
#'   `list(A = predict_a, B = predict_b)` (the other models use the automatic
#'   routes). If voipred has no way to predict from a fitted model and no
#'   function is given for it, `nb_val()` stops before predicting, naming the
#'   models.
#' @param label Optional names of the models in printed output and plots: a
#'   single string, or one per model, in order. Default: see "Labels".
#' @return An object of class `voipred_nb`, a list with:
#'   * `estimates`: a data frame with one row per model and threshold (by
#'     model, then threshold) and columns `threshold`, `model` (its label),
#'     `NB_model` (using the model; `NA` where it cannot be
#'     estimated, see Details), `NB_all` (treat everyone), `NB_none` (treat no
#'     one, always 0) and `delta_NB`, the difference in net benefit between the
#'     model and the best default strategy,
#'     \eqn{NB_{model} - \max(NB_{all}, NB_{none})}: positive where the model
#'     is useful;
#'   * `models` (the labels of the models), `call` and `prior`;
#'   * `outcome`: the outcome type (`"binary"`, `"survival"` or
#'     `"competing"`), the number of patients `n`, the `prevalence` (binary:
#'     the posterior mean event rate; survival: the estimated risk of the
#'     event, or of `cause`, by the horizon), `events` (binary), and for
#'     survival outcomes `horizon`, `cause` and `censoring` (`"km"` or
#'     `"cipcw"`);
#'   * `risk`: the predicted risks, a matrix with one row per row of `data`
#'     and one column per model (named by its label);
#'   * `counts` (binary): the numbers of true and false positives and
#'     negatives of each model at each threshold;
#'   * `observed`, `weights` and `cipcw`: the outcome, case weights and
#'     censoring settings, used by the value-of-information functions.
#'
#'   Like a `glm` fit, it is a list whose parts are read with `$` (e.g.
#'   `res$estimates`), and it has methods: `print()` shows what was evaluated
#'   and the first rows of `estimates` (several models: their net benefits
#'   side by side); [summary()][summary.voipred_nb] gives
#'   the best strategy, sensitivity and specificity at each threshold;
#'   [plot()][plot.voipred_nb] draws the decision curve or the
#'   \eqn{\Delta NB} curve of every model; `as.data.frame()` returns
#'   `estimates`; [c()][c.voipred_nb] combines results.
#' @seealso [summary.voipred_nb()], [plot.voipred_nb()], [c.voipred_nb()], [nb_priors],
#'   [cipcw_control()], and the value-of-information functions [p_useful()],
#'   [p_best()], [evpi()] and [evsi()].
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:10000, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#'
#' # Binary outcome: the predicted risks are all that is needed
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb_val(data = df_val_binary,
#'               predicted_risk = df_val_binary$p,
#'               outcome = "day30",
#'               thresholds = seq(0, 0.3, by = 0.01))
#' res
#' summary(res)
#' plot(res)
#'
#' # The same, with the name of the column of `data` that holds the risks
#' nb_val(data = df_val_binary,
#'        predicted_risk = "p",
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1))
#'
#' # Several models: their risks in one call
#' fit_age <- glm(day30 ~ age, family = binomial, data = df_dev_binary)
#' df_val_binary$p_age <- predict(fit_age, newdata = df_val_binary, type = "response")
#' both <- nb_val(data = df_val_binary,
#'                predicted_risk = list(Full = df_val_binary$p, `Age only` = df_val_binary$p_age),
#'                outcome = "day30",
#'                thresholds = seq(0, 0.3, by = 0.01))
#' both
#' summary(both, thresholds = c(0.05, 0.1))
#' plot(both)
#' # the same with column names, with the fitted models, or a mix
#' nb_val(data = df_val_binary,
#'        predicted_risk = list(Full = "p", `Age only` = "p_age"),
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1))
#' nb_val(data = df_val_binary,
#'        predicted_risk = list(Full = fit, `Age only` = fit_age),
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1))
#' nb_val(data = df_val_binary,
#'        predicted_risk = list(Full = fit, `Age only` = "p_age"),
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1))
#'
#' # Optionally, the model instead of the risks: by default, or with your own
#' # prediction function
#' nb_val(data = df_val_binary,
#'        predicted_risk = fit,
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1))
#' nb_val(data = df_val_binary,
#'        predicted_risk = fit,
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1),
#'        predict_function = function(model, newdata, ...)
#'      predict(model, newdata, type = "response"))
#'
#' # Optional prior information (binary outcome); the default uses the data only
#' nb_val(data = df_val_binary,
#'        predicted_risk = df_val_binary$p,
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1),
#'        prior = nb_uniform())
#' nb_val(data = df_val_binary,
#'        predicted_risk = df_val_binary$p,
#'        outcome = "day30",
#'        thresholds = c(0.05, 0.1),
#'        prior = nb_beta(prevalence = c(7, 93)))
#'
#' # Survival: death by 5 years (Rotterdam breast cancer data)
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
#' nb_val(data = df_val_survival,
#'        predicted_risk = df_val_survival$p5,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3))
#'
#' # Optionally, the model instead of the risks: by default, or with your own
#' # prediction function (here the same risks as above)
#' nb_val(data = df_val_survival,
#'        predicted_risk = cox,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3))
#' nb_val(data = df_val_survival,
#'        predicted_risk = cox,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3),
#'        predict_function = function(model, newdata, horizon, ...)
#'      1 - summary(survfit(model, newdata = newdata), times = horizon)$surv[1, ])
#'
#' # The same with conditional IPCW: censoring depends on age, nodes and grade
#' # (a Cox model of the censoring times, fitted by voipred)
#' nb_val(data = df_val_survival,
#'        predicted_risk = df_val_survival$p5,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3),
#'        cipcw = ~ age + nodes + grade)
#'
#' # or your own censoring model, fitted to the validation data
#' cens <- coxph(Surv(years, death == 0) ~ age + nodes + grade,
#'               data = df_val_survival)
#' nb_val(data = df_val_survival,
#'        predicted_risk = df_val_survival$p5,
#'        outcome = ~ Surv(years, death),
#'        horizon = 5,
#'        thresholds = c(0.1, 0.2, 0.3),
#'        cipcw = cens)
#' @export
nb_val <- function(data, predicted_risk = NULL, outcome,
                   thresholds = seq(0, 0.99, by = 0.01), horizon = NULL, cause = NULL,
                   cipcw = NULL, weights = NULL,
                   prior = nb_improper(), predict_function = NULL,
                   label = NULL) {
  call <- match.call()
  if (!is.data.frame(data)) {
    if (is.data.frame(predicted_risk))
      stop("`data` must be a data frame. The first argument of nb_val() is `data`: ",
           "nb_val(data, predicted_risk, outcome, ...).", call. = FALSE)
    stop("`data` must be a data frame.", call. = FALSE)
  }
  # the models: predicted risks or fitted models, one or several
  inputs <- model_inputs(predicted_risk, substitute(predicted_risk), data)
  labels <- model_labels(inputs, label)
  M <- length(inputs)
  # censoring: Kaplan-Meier, or cIPCW when `cipcw` is given
  cipcw <- as_cipcw(cipcw)
  censoring <- if (is.null(cipcw)) "km" else "cipcw"
  info <- parse_outcome(outcome, data)
  k <- check_survival_args(info, horizon, cause, cipcw)
  info$k <- k
  thresholds <- check_thresholds(thresholds)
  n <- nrow(data)
  w <- check_weights(weights, n)
  check_prior(prior, info, w)
  pf <- check_predict_function(predict_function, inputs, labels)
  check_prediction_routes(inputs, labels, pf, info)

  # the risks (one column per model): given, or predicted once from each model
  P <- matrix(NA_real_, n, M, dimnames = list(NULL, labels))
  for (m in seq_len(M)) {
    who <- if (M > 1) labels[m]
    P[, m] <- if (inputs[[m]]$is_model)
      with_model_name(get_risk(inputs[[m]]$value, data, info, horizon, pf[[labels[m]]]), who)
    else given_risk(inputs[[m]]$value, data, who)
  }

  # the net benefit of each model; that of treating all is the same for all
  outcome_info <- list(type = info$type, n = n)
  cw <- NULL
  counts <- NULL
  nb_model <- vector("list", M)
  if (info$type == "binary") {
    for (m in seq_len(M)) {
      res <- bayes_nb_mean(info$y, P[, m], thresholds, prior, w)
      nb_model[[m]] <- res$NB_model
    }
    observed <- list(y = info$y)
    outcome_info$events <- if (is.null(w)) sum(info$y) else NA_real_
    counts <- do.call(rbind, lapply(seq_len(M), function(m)
      binary_counts(info$y, P[, m], thresholds, w, labels[m])))
  } else {
    observed <- list(time = info$time, status = info$status)
    outcome_info[c("horizon", "cause", "k", "states", "censoring")] <-
      list(horizon, info$states[k], k, info$states, censoring)
    if (censoring == "km") {
      for (m in seq_len(M)) {
        res <- nb_core_km(info$time, info$status, k, P[, m], thresholds, horizon, w,
                          label = if (M > 1) labels[m])
        nb_model[[m]] <- res$NB_model
      }
    } else {
      # one censoring model, and one set of censoring weights, for every model;
      # with a formula, a Cox model fitted here to the censoring times of `data`
      cdata <- data
      if (!is.null(cipcw$formula)) {
        fc <- fit_censoring_cox(cipcw$formula, data, info, w)
        cipcw$model <- fc$model
        cdata <- fc$data
      }
      cws <- cipcw_weights(cipcw, cdata, info, horizon, k)
      for (m in seq_len(M)) {
        res <- nb_core_cipcw(info$time, info$status, k, P[, m], thresholds, horizon,
                             cws$weight_event, cws$weight_horizon, w)
        nb_model[[m]] <- res$NB_model
      }
      # what refitting the censoring model in each draw of the value of
      # information needs
      cw <- list(control = cipcw, weight_event = cws$weight_event,
                 weight_horizon = cws$weight_horizon, data = cdata, info = info)
    }
  }
  outcome_info$prevalence <- res$prevalence

  est <- do.call(rbind, lapply(seq_len(M), function(m)
    data.frame(threshold = thresholds, model = labels[m], NB_model = nb_model[[m]],
               NB_all = res$NB_all, NB_none = 0)))
  est$delta_NB <- est$NB_model - pmax(est$NB_all, 0)
  rownames(est) <- NULL
  structure(list(call = call, models = labels, estimates = est, prior = prior,
                 outcome = outcome_info, risk = P, observed = observed, weights = w,
                 counts = counts, cipcw = cw),
            class = "voipred_nb")
}

# Counts (weighted with case weights) of treated and untreated patients with
# and without the event at each threshold, for the model `label`.
binary_counts <- function(y, p, thresholds, w = NULL, label = "Model") {
  if (is.null(w)) w <- rep(1, length(y))
  tr <- vapply(thresholds, function(z) treated_at(p, z), logical(length(y)))
  if (!is.matrix(tr)) tr <- matrix(tr, nrow = length(y))
  data.frame(threshold = thresholds, model = label,
             TP = colSums(w * y * tr), FP = colSums(w * (1 - y) * tr),
             FN = colSums(w * y * !tr), TN = colSums(w * (1 - y) * !tr))
}

# The thresholds and the net benefit of treating all of a result of nb_val()
# (the same for every model).
nb_thresholds <- function(x) x$estimates$threshold[x$estimates$model == x$models[1]]
nb_all <- function(x) x$estimates$NB_all[x$estimates$model == x$models[1]]

# What the value-of-information functions need from a result of nb_val(), as
# one list (internal). `risk` has one column per model.
voi_info <- function(r) {
  o <- r$outcome
  v <- list(type = o$type, thresholds = nb_thresholds(r), risk = r$risk,
            weights = r$weights, prior = r$prior)
  if (o$type == "binary") {
    v$y <- r$observed$y
  } else {
    v[c("time", "status", "k", "horizon", "censoring")] <-
      list(r$observed$time, r$observed$status, o$k, o$horizon, o$censoring)
    if (identical(o$censoring, "cipcw"))
      v[c("weight_event", "weight_horizon", "cipcw", "data", "info")] <-
        list(r$cipcw$weight_event, r$cipcw$weight_horizon, r$cipcw$control,
             r$cipcw$data, r$cipcw$info)
  }
  v
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
check_survival_args <- function(info, horizon, cause, cipcw) {
  if (info$type == "binary") {
    if (!is.null(horizon))
      stop("`horizon` is only used with a survival outcome made with Surv().",
           call. = FALSE)
    if (!is.null(cause))
      stop("`cause` is only used with a competing-risks outcome.", call. = FALSE)
    if (!is.null(cipcw))
      stop("`cipcw` is only used with a survival outcome.", call. = FALSE)
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

nb_core_km <- function(time, status, k, p, thresholds, horizon, w = NULL,
                       label = NULL) {
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
    warning(if (!is.null(label)) paste0("\"", label, "\": "),
            "No patient treated at ", if (sum(undefined) == 1) "threshold " else
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
       prevalence = prev, Se = se, Sp = sp)
}

# ---- internal: the models (predicted risks or fitted models) -----------------

# Risks (a vector, or a column name) rather than a fitted model.
is_risk_input <- function(x) is.atomic(x)

# A list made with list(), not a fitted model (most models are lists with a
# class).
is_plain_list <- function(x) is.list(x) && !is.object(x)

# The models given with `predicted_risk`, one element each:
# list(value, is_model, label), with label NULL when the call gives none. One
# model: risks (a vector or a column name) or a fitted model; several: a list
# whose elements are any of these.
model_inputs <- function(predicted_risk, risk_expr, data) {
  x <- predicted_risk
  if (is.null(x))
    stop("Give the predicted risks with `predicted_risk`: a vector of risks, the ",
         "name of the column of `data` that holds them, or a fitted model (several ",
         "models: a list of these).", call. = FALSE)
  if (is.data.frame(x) || (is.matrix(x) && ncol(x) > 1))
    stop("To give several models, use a list: predicted_risk = list(A = \"p1\", ",
         "B = \"p2\"), with column names, vectors of risks or fitted models.",
         call. = FALSE)
  if (is.character(x) && length(x) > 1) {
    # risks stored as text (one per row) rather than several column names
    if (length(x) == nrow(data) && !any(x %in% names(data)))
      stop("`predicted_risk` must be numeric: one predicted probability per row of ",
           "`data`, or the name of the column that holds them.", call. = FALSE)
    nms <- if (is.null(names(x)) || !all(nzchar(names(x)))) LETTERS[seq_along(x)]
           else names(x)
    stop("To give several models, use list() rather than c(): predicted_risk = list(",
         paste0(nms, " = \"", x, "\"", collapse = ", "), ").", call. = FALSE)
  }
  if (is_plain_list(x)) return(list_inputs(x, risk_expr))
  if (is_risk_input(x)) return(list(model_input(x, FALSE, risk_label(risk_expr, x))))
  # a fitted model
  list(model_input(x, TRUE, risk_label(risk_expr, NULL)))
}

model_input <- function(value, is_model, label) list(value = value, is_model = is_model,
                                                     label = label)

# A list of models (fitted models, column names or vectors of risks),
# labelled by its names, or else by the expressions in the call
# (list(val$p1, fit2) gives "p1" and "fit2").
list_inputs <- function(x, expr) {
  if (!length(x)) stop("`predicted_risk` is an empty list.", call. = FALSE)
  nms <- names(x)
  if (is.null(nms)) nms <- rep("", length(x))
  args <- if (is.call(expr) && identical(expr[[1]], as.name("list")))
    as.list(expr)[-1]
  if (length(args) != length(x)) args <- NULL
  lapply(seq_along(x), function(i) {
    v <- x[[i]]
    if (is.null(v) || is.data.frame(v) || is_plain_list(v) || (is.matrix(v) && ncol(v) > 1) ||
        (is.character(v) && length(v) != 1))
      stop("Each element of the list `predicted_risk` must be one model: a vector of ",
           "risks, the name of a column of `data`, or a fitted model",
           if (nzchar(nms[i])) paste0(" (element \"", nms[i], "\" is not)")
           else paste0(" (element ", i, " is not)"), ".", call. = FALSE)
    is_model <- !is_risk_input(v)
    lab <- if (nzchar(nms[i])) nms[i]
           else if (!is.null(args)) risk_label(args[[i]], if (is_model) NULL else v)
           else if (is.character(v)) v
    model_input(v, is_model, lab)
  })
}

# Label of risks given directly: the column name, or the name of the vector
# (`p`, or `p` in `val$p`); NULL if there is none.
risk_label <- function(expr, risk) {
  if (is.character(risk) && length(risk) == 1) return(risk)
  if (is.name(expr)) return(as.character(expr))
  if (is.call(expr) && identical(expr[[1]], as.name("$"))) return(as.character(expr[[3]]))
  NULL
}

# The labels of the models: `label` if given, otherwise those found in the
# call, otherwise "Model" (one model) or "Model 1", "Model 2", ...
model_labels <- function(inputs, label) {
  M <- length(inputs)
  if (!is.null(label)) {
    if (!is.character(label) || length(label) != M || anyNA(label) || !all(nzchar(label)))
      stop(if (M == 1) "`label` must be a single non-empty string."
           else paste0("`label` must be ", M, " non-empty strings, one per model."),
           call. = FALSE)
    lab <- unname(label)
  } else {
    lab <- vapply(inputs, function(i) if (is.null(i$label)) NA_character_ else i$label,
                  character(1))
    miss <- is.na(lab)
    lab[miss] <- if (M == 1) "Model" else paste("Model", which(miss))
  }
  check_model_labels(lab, "give each a different one with `label`, or with names")
  lab
}

check_model_labels <- function(lab, how) {
  bad <- tolower(lab) %in% c("treat all", "treat none")
  if (any(bad))
    stop("A model cannot be labelled \"", lab[bad][1], "\": \"Treat All\" and ",
         "\"Treat None\" are the default strategies.", call. = FALSE)
  if (anyDuplicated(lab))
    stop("Two models have the same label (\"", lab[anyDuplicated(lab)], "\"); ", how, ".",
         call. = FALSE)
  invisible(lab)
}

# `predict_function` as a list of functions named by the labels of the models
# they predict from: one function for all the models, or one per model.
check_predict_function <- function(pf, inputs, labels) {
  if (is.null(pf)) return(list())
  is_model <- vapply(inputs, function(i) i$is_model, logical(1))
  if (!any(is_model))
    stop("`predict_function` predicts risks from fitted models, but `predicted_risk` ",
         "has no fitted model.", call. = FALSE)
  if (is.function(pf)) return(stats::setNames(rep(list(pf), sum(is_model)), labels[is_model]))
  if (is_plain_list(pf) && length(pf) && !is.null(names(pf)) && all(nzchar(names(pf))) &&
      all(vapply(pf, is.function, logical(1)))) {
    bad <- setdiff(names(pf), labels[is_model])
    if (length(bad))
      stop("`predict_function` has a function for \"", bad[1], "\", which is not one of ",
           "the fitted models: ", paste0("\"", labels[is_model], "\"", collapse = ", "), ".",
           call. = FALSE)
    return(pf)
  }
  stop("`predict_function` must be a function(model, newdata, horizon, cause), or a ",
       "list of such functions named by the labels of the models.", call. = FALSE)
}

# Every fitted model needs a way to predict its risks: voipred's own route for
# its class, riskRegression::predictRisk(), predict() (binary outcomes), or
# `predict_function`. Checked before predicting, for all the models at once.
check_prediction_routes <- function(inputs, labels, pf, info) {
  miss <- vapply(seq_along(inputs), function(m)
    inputs[[m]]$is_model && is.null(pf[[labels[m]]]) &&
      !has_route(inputs[[m]]$value, info), logical(1))
  if (!any(miss)) return(invisible())
  cls <- vapply(inputs[miss], function(i) class(i$value)[1], character(1))
  several <- length(inputs) > 1
  many <- sum(miss) > 1
  what <- if (several) paste0("\"", labels[miss], "\" (", cls, ")", collapse = ", ")
          else paste0("this ", cls, " model")
  stop("voipred cannot predict risks from ", what, ": there is no built-in way to ",
       "predict from ", if (many) "these models" else "it", " (see \"Predicting the ",
       "risks from a model\" in ?nb_val). If ", if (many) "they come" else "it comes",
       " from a package, load that package first. Otherwise give `predict_function`, ",
       "a function(model, newdata, horizon, cause) returning the risk of the event",
       if (info$type != "binary") " by the horizon", ", e.g. ",
       if (several) paste0("predict_function = list(", r_name(labels[miss][1]),
                           " = my_predict)")
       else "predict_function = my_predict", ".", call. = FALSE)
}

# A label as it is written in R code (with backticks when it is not a
# syntactic name).
r_name <- function(x) if (identical(make.names(x), x)) x else paste0("`", x, "`")

# Whether voipred knows how to predict risks from `model` for this outcome.
# Each generic is looked up from where it is defined (methods defined by the
# user in the global environment are found too).
has_route <- function(model, info) {
  cls <- class(model)
  found <- function(generic, envir)
    any(vapply(cls, function(k) !is.null(utils::getS3method(generic, k, optional = TRUE,
                                                            envir = envir)),
               logical(1)))
  if (found("predict_risk", topenv(environment(has_route))) ||
      found("predictRisk", asNamespace("riskRegression"))) return(TRUE)
  info$type == "binary" && found("predict", globalenv())
}

# Errors while predicting from one of several models name that model.
with_model_name <- function(expr, label) {
  if (is.null(label)) return(expr)
  tryCatch(expr, error = function(e)
    stop("\"", label, "\": ", conditionMessage(e), call. = FALSE))
}

# Predicted risks given with `predicted_risk`: a vector, or the name of a
# column of `data`. `label` names the model in the messages (several models).
given_risk <- function(risk, data, label = NULL) {
  arg <- if (is.null(label)) "`predicted_risk`"
         else paste0("`predicted_risk` (\"", label, "\")")
  n <- nrow(data)
  if (is.character(risk) && length(risk) == 1) {
    if (!risk %in% names(data))
      stop(arg, " names a column (\"", risk, "\") that is not in `data`.", call. = FALSE)
    risk <- data[[risk]]
  }
  if (!is.numeric(risk) || (is.matrix(risk) && ncol(risk) != 1))
    stop(arg, " must be numeric: one predicted probability per row of ",
         "`data`, or the name of the column that holds them.", call. = FALSE)
  risk <- as.vector(risk)
  if (length(risk) != n) {
    k <- length(risk) / n
    stop(arg, " has ", length(risk), " values but `data` has ", n, " rows.",
         if (k >= 2 && k == round(k))
           paste0(" To give ", k, " models, use list() rather than c(), which joins ",
                  "the risks into one vector: predicted_risk = list(A = p1, B = p2)."),
         call. = FALSE)
  }
  if (anyNA(risk))
    stop(sum(is.na(risk)), " of the risks in ", arg, " are NA.", call. = FALSE)
  if (any(risk < 0 | risk > 1))
    stop(arg, " has values outside [0, 1] (from ", signif(min(risk), 3),
         " to ", signif(max(risk), 3), "), so these are not probabilities (e.g. ",
         "log-odds must be turned into probabilities first).", call. = FALSE)
  if (all(risk == 0 | risk == 1))
    stop(arg, " has only 0s and 1s. These look like predicted classes, ",
         "not probabilities.", call. = FALSE)
  unname(risk)
}

check_prior <- function(prior, info, weights) {
  if (!inherits(prior, "nb_prior"))
    stop("`prior` must be nb_improper(), nb_uniform(), nb_jeffreys(), nb_beta() or ",
         "nb_informative().", call. = FALSE)
  if (prior$type == "improper") return(invisible())
  if (info$type != "binary")
    stop("For survival and competing-risks outcomes the prior is nb_improper() ",
         "(the data only); other priors are for binary outcomes.", call. = FALSE)
  if (!is.null(weights))
    stop("Case weights can only be used with the improper prior, nb_improper().",
         call. = FALSE)
  invisible()
}
