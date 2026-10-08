#' Predicted risks from a fitted model
#'
#' Turns a fitted model into one predicted risk per row of `newdata`: the
#' probability of the event (binary outcome), of the event by `horizon`
#' (survival outcome), or of the event `cause` by `horizon` (competing risks).
#'
#' voipred does not keep a list of supported models. Risks are obtained in
#' this order:
#'
#' 1. **A `predict_risk()` method for the model's class**, if one exists. This
#'    is an ordinary S3 generic: anyone (a user, or the author of a modelling
#'    package) can add a method, e.g.
#'    `predict_risk.mymodel <- function(model, newdata, horizon = NULL, cause = NULL, ...)`,
#'    and voipred then works with that model. `horizon = NULL` means a
#'    binary outcome, `cause = NULL` a non-competing outcome; the method may
#'    return any of the output shapes described below. voipred itself only
#'    ships methods for models that the standard routes below cannot handle.
#' 2. **The standard routes** (the default method):
#'    * binary outcome: `predict(model, newdata, type = "response")`, then
#'      `predict(model, newdata, type = "prob")`, then
#'      `riskRegression::predictRisk(model, newdata)`;
#'    * survival outcome:
#'      `riskRegression::predictRisk(model, newdata, times = horizon)`;
#'    * competing risks: the same with `cause = cause` (the event level's
#'      label), and if that fails `cause = cause_index` (its position among
#'      the event levels), for models that number their causes.
#'
#'    If a route fails, it is tried again with only the model's covariate
#'    columns of `newdata` (some `predict()` methods trip over other columns).
#'
#' In [nb()], a `predict_function` takes precedence over both, and works for any
#' model.
#'
#' Whatever the route, the output is read in a generic way: a vector of
#' probabilities; one column per outcome class (the event column is taken:
#' column `"1"` or `".pred_1"`, otherwise the second of two); or a survival
#' curve per patient, `list(time = <grid>, surv = <n x grid matrix>)` (or
#' `cif =` for a cumulative incidence), read at the horizon. It must then be
#' probabilities: values outside \[0, 1\] (e.g. log-odds, relative risks,
#' predicted times) and predictions that are only 0s and 1s (class labels) are
#' rejected, as are `NA`s and a wrong number of predictions.
#'
#' **Works through the standard routes** (checked in the tests against each
#' package's own predictions): binary: `glm`, `mgcv::gam`, `logistf`, `rpart`,
#' `partykit::ctree`, `randomForest`, `ranger` (`probability = TRUE`), `gbm`,
#' `xgboost()`, `naivebayes::naive_bayes`, `earth`, `kknn` (attached),
#' bayespred `bpmfit` and `bpm`; survival: `coxph` (fitted with `x = TRUE`),
#' `rms::cph`, `survreg`, `rms::psm`, `flexsurvreg`, `flexsurvspline`,
#' `ranger`, `randomForestSRC`, `prodlim`; competing risks:
#' `riskRegression::CSC`, `riskRegression::FGR`, competing-risks
#' `randomForestSRC`.
#'
#' **Methods shipped by voipred**, only where the standard routes fail or
#' give a wrong answer: binary: `e1071::svm`, `e1071::naiveBayes`, `nnet`,
#' `MASS::lda`/`qda`, `glmnet`/`cv.glmnet`, `xgb.train()` boosters (and the
#' level order of `xgboost()`), `lightgbm`, `catboost`; survival: `coxph`
#' without `x = TRUE` or with case weights, stratified `survreg`,
#' `flexsurvreg` when its formula was passed as a variable, `aorsf`,
#' `partykit::ctree` survival trees, `rstpm2::stpm2` (attach rstpm2),
#' `mgcv::gam` with `cox.ph()`, bayespred.dev Bayesian Cox and parametric
#' survival models (`bpmfit` survival fits, and `bpm` objects packaged at the
#' horizon with `posterior(fit, time = )`; the risk is the posterior mean of
#' risk); competing risks: `cmprsk::crr`, survival's
#' multi-state `coxph` (`coxph(Surv(time, factor_status) ~ x, id = id)`:
#' cause-specific Cox models combined by Aalen-Johansen, as
#' `riskRegression::CSC`), and a check that an `FGR` model was fitted for the
#' requested cause. Models trained on a
#' numeric matrix need `newdata` columns with the same names as the features.
#'
#' Models that only give relative risks or times (e.g. `xgboost` survival
#' objectives, `gbm` and `glmnet` Cox models, `rpart` survival trees) have no
#' route to probabilities: voipred stops and asks for `predict_function`, or the
#' \[0, 1\] check stops their output.
#'
#' Before predicting, the covariates the model uses (read with
#' `insight::find_predictors()`, or from a formula the model keeps) are
#' compared with the columns of `newdata`; missing ones are listed in an
#' error. When they cannot be read (e.g. `ranger`), the model's own error is
#' shown instead.
#'
#' @param model A fitted model. A vector of risks is not accepted.
#' @param newdata A data frame with the covariates the model needs; one risk
#'   is returned per row, in the same order.
#' @param horizon For survival and competing-risks outcomes: the time by
#'   which the risk is predicted. `NULL` for a binary outcome.
#' @param cause For competing risks: the cause of interest (the outcome's
#'   event level). `NULL` otherwise.
#' @param ... Passed to methods. voipred also passes `cause_index`, the
#'   position of `cause` among the outcome's event levels.
#' @return A numeric vector of length `nrow(newdata)` with values in \[0, 1\],
#'   checked as described above. Calling `predict_risk()` directly gives the
#'   risks [nb()] uses (without `predict_function`).
#' @examples
#' fit <- glm(day30 ~ age + Killip, family = binomial, data = gusto[1:2000, ])
#' p <- predict_risk(fit, gusto[2001:3000, ])
#' summary(p)
#'
#' library(survival)
#' cox <- coxph(Surv(dtime, death) ~ age + nodes, data = rotterdam[1:1500, ], x = TRUE)
#' summary(predict_risk(cox, rotterdam[1501:2000, ], horizon = 1826))
#' @export
predict_risk <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  UseMethod("predict_risk")
}

#' @rdname predict_risk
#' @export
predict_risk.default <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  check_model_object(model)
  check_covariates(model, newdata)
  out <- standard_route(model, newdata, horizon, cause, ...)
  # some predict() methods trip over columns the model does not use (e.g. a
  # partial outcome, or missing values in unused columns): retry without them
  vars <- model_covariates(model)
  if (!usable(out) && !is.null(vars) && length(vars) < ncol(newdata))
    out2 <- standard_route(model, newdata[, intersect(vars, names(newdata)), drop = FALSE],
                           horizon, cause, ...)
  else out2 <- NULL
  if (usable(out2)) out <- out2
  if (!usable(out))
    stop("Could not get predicted risks from this ", class(model)[1], " model",
         if (inherits(out, "error")) paste0(" (", trimws(conditionMessage(out)), ")")
         else if (!is.null(out)) " (it returned only NA; is `horizon` beyond the model's follow-up?)",
         ". If a covariate is missing from the data, add it. Otherwise give ",
         "`predict_function = function(model, newdata, horizon, cause)` to nb(), ",
         "returning the risk of the event", if (!is.null(horizon)) " by the horizon",
         ", or add a predict_risk() method for this class (see ?predict_risk).",
         call. = FALSE)
  finish_risk(out, newdata, horizon)
}

standard_route <- function(model, newdata, horizon, cause, ...) {
  if (is.null(horizon)) {
    out <- try_route(stats::predict(model, newdata = newdata, type = "response"))
    if (!usable(out)) out <- try_route(stats::predict(model, newdata = newdata, type = "prob"))
    if (!usable(out)) out <- try_route(riskRegression::predictRisk(model, newdata = newdata))
  } else if (is.null(cause)) {
    out <- try_route(riskRegression::predictRisk(model, newdata = newdata, times = horizon))
  } else {
    # the cause as labelled in the outcome; if the model codes causes by
    # number, its position among the outcome's event levels
    out <- try_route(riskRegression::predictRisk(model, newdata = newdata,
                                                 times = horizon, cause = cause))
    k <- list(...)$cause_index
    if (!usable(out) && !is.null(k))
      out <- try_route(riskRegression::predictRisk(model, newdata = newdata,
                                                   times = horizon, cause = k))
  }
  out
}

try_route <- function(expr) {
  tryCatch(suppressWarnings(expr), error = function(e) e)
}

# A route's output is usable if it is numeric (or a data frame / list of
# shapes as_risk() understands) and not entirely missing.
usable <- function(out) {
  if (is.null(out) || inherits(out, "error")) return(FALSE)
  if (is.list(out) && !is.data.frame(out)) return(all(c("time") %in% names(out)))
  if (is.data.frame(out)) out <- as.matrix(out)
  is.numeric(out) && !all(is.na(out))
}

finish_risk <- function(out, newdata, horizon = NULL) {
  check_risk(as_risk(out, horizon), newdata)
}

check_model_object <- function(model) {
  if (is.numeric(model) || is.logical(model))
    stop("`model` must be a fitted model, not a vector of predicted risks. ",
         "voipred needs the model itself.", call. = FALSE)
}


# ---- reading and checking predictions (generic) -----------------------------

# Common output shapes -> one risk per row.
as_risk <- function(out, horizon = NULL) {
  if (is.list(out) && !is.data.frame(out) && "time" %in% names(out)) {
    if (is.null(horizon))
      stop("A survival curve was returned, but there is no `horizon`.", call. = FALSE)
    i <- findInterval(horizon, out$time)                 # last grid time <= horizon
    if (!is.null(out$cif)) {
      cif <- as.matrix(out$cif)
      return(if (i == 0) rep(0, nrow(cif)) else unname(cif[, i]))
    }
    surv <- as.matrix(out$surv)
    return(if (i == 0) rep(0, nrow(surv)) else unname(1 - surv[, i]))
  }
  if (is.data.frame(out) || is.matrix(out)) return(event_column(as.matrix(out)))
  if (is.factor(out)) return(out)      # classes: rejected by check_risk()
  as.vector(unname(out))
}

# One column per outcome class: take the event's ("1", ".pred_1"), otherwise
# the second of two. A single column is the risk itself.
event_column <- function(pr) {
  if (ncol(pr) == 1) return(as.vector(unname(pr[, 1])))
  cn <- colnames(pr)
  hit <- which(cn %in% c("1", ".pred_1", "TRUE", "X1"))
  if (length(hit) > 0) return(as.vector(unname(pr[, hit[1]])))
  if (ncol(pr) != 2)
    stop("The model predicts ", ncol(pr), " outcome classes; the outcome ",
         "must be binary.", call. = FALSE)
  as.vector(unname(pr[, 2]))
}

check_risk <- function(p, newdata = NULL) {
  if (is.factor(p) || is.character(p))
    stop("The model returned classes, not probabilities. Make it return ",
         "probabilities (or give `predict_function`).", call. = FALSE)
  if (is.logical(p)) p <- as.integer(p)
  if (!is.numeric(p))
    stop("Predicted risks must be numeric.", call. = FALSE)
  if (!is.null(newdata) && length(p) != nrow(newdata))
    stop("The model returned ", length(p), " predictions for ", nrow(newdata),
         " rows. Check for missing values in the covariates.", call. = FALSE)
  if (anyNA(p))
    stop(sum(is.na(p)), " predicted risks are NA. Check for missing values ",
         "in the covariates.", call. = FALSE)
  if (any(p < 0 | p > 1))
    stop("The model returned values outside [0, 1] (from ",
         signif(min(p), 3), " to ", signif(max(p), 3), "), so these are not ",
         "probabilities. They may be scores from the model (e.g. log-odds, ",
         "relative risks or predicted times) that must be turned into ",
         "probabilities first.", call. = FALSE)
  if (all(p == 0 | p == 1))
    stop("The model returned only 0s and 1s. These look like predicted ",
         "classes (labels), not probabilities. Make the model return ",
         "probabilities (e.g. fit ranger or svm with `probability = TRUE`).",
         call. = FALSE)
  as.vector(p)
}


# ---- covariates and the model's own outcome (generic, via insight) ----------

# Names of the covariates a model uses, or NULL if they cannot be read.
model_covariates <- function(model) {
  v <- tryCatch(suppressWarnings(suppressMessages(
    insight::find_predictors(model, flatten = TRUE, verbose = FALSE))),
    error = function(e) NULL)
  if (length(v) == 0) {
    # any formula or terms object kept inside the model (e.g. sub-models)
    v <- unique(unlist(lapply(stored_formulas(model), function(f)
      tryCatch(all.vars(stats::delete.response(stats::terms(f))), error = function(e) NULL))))
  }
  if (length(v) == 0) NULL else unique(setdiff(v, "."))
}

# Formulas / terms stored in a model object (top level and one level down),
# found without knowing the class.
stored_formulas <- function(model) {
  out <- list()
  visit <- function(x, depth) {
    if (inherits(x, "formula")) { out[[length(out) + 1]] <<- x; return(invisible()) }
    if (depth < 4 && is.list(x) && !is.data.frame(x))
      for (el in x) visit(el, depth + 1)
  }
  if (is.list(model)) visit(unclass(model), 1)
  out
}

check_covariates <- function(model, newdata, what = "model", needed = NULL) {
  if (!is.data.frame(newdata))
    stop("`data` must be a data frame.", call. = FALSE)
  if (is.null(needed)) needed <- model_covariates(model)
  if (is.null(needed)) return(invisible(TRUE))
  missing <- setdiff(needed, names(newdata))
  if (length(missing) > 0)
    stop("The validation data is missing ", length(missing),
         if (length(missing) == 1) " covariate" else " covariates",
         " that the ", what, " was ",
         if (what == "model") "developed" else "fitted", " with: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  invisible(TRUE)
}

# The model's formula, if it can be read.
model_formula <- function(model) {
  f <- tryCatch(suppressWarnings(insight::find_formula(model, verbose = FALSE)$conditional),
                error = function(e) NULL)
  if (is.null(f))
    f <- tryCatch(suppressWarnings(stats::formula(model)), error = function(e) NULL)
  if (!inherits(f, "formula") || length(f) < 3) NULL else f
}

# The kind of outcome the model was fitted to: from its formula evaluated in
# `data`, otherwise from the response insight reads from the model.
model_outcome_kind <- function(model, data) {
  kind <- outcome_kind(model_lhs(model, data))
  if (kind != "unknown") return(kind)
  kind <- outcome_kind(stored_surv(model))          # e.g. the Surv / Hist it keeps
  if (kind != "unknown") return(kind)
  y <- tryCatch(suppressWarnings(suppressMessages(insight::get_response(model, verbose = FALSE))),
                error = function(e) NULL)
  if (is.null(y) || is.data.frame(y)) return("unknown")
  if (outcome_kind(y) == "binary") "binary" else "unknown"
}

# The left-hand side of the model's formula evaluated in `data` (e.g. the
# Surv() object), or NULL.
model_lhs <- function(model, data) {
  f <- model_formula(model)
  if (is.null(f)) return(NULL)
  env <- new.env(parent = if (is.null(environment(f))) globalenv() else environment(f))
  env$Surv <- survival::Surv
  if (requireNamespace("prodlim", quietly = TRUE)) env$Hist <- prodlim::Hist
  tryCatch(suppressWarnings(eval(f[[2]], data, env)), error = function(e) NULL)
}

# "binary", "survival", "competing" or "unknown", from the model's own outcome.
outcome_kind <- function(y) {
  if (is.null(y)) return("unknown")
  if (inherits(y, "Surv")) {
    tp <- attr(y, "type")
    # Surv(time, status) with numeric codes 0/1/2 (a competing-risks coding)
    # turns the codes above 1 into NA
    if (tp == "right" && anyNA(y[, "status"]) && !anyNA(y[, "time"])) return("competing")
    return(if (tp == "mright") "competing" else if (tp == "right") "survival" else "unknown")
  }
  if (inherits(y, "Hist")) {
    return(if (identical(attr(y, "model"), "competing.risks")) "competing" else "survival")
  }
  if (is.logical(y) || (is.factor(y) && nlevels(y) == 2) ||
      (is.numeric(y) && all(y[!is.na(y)] %in% c(0, 1))))
    return("binary")
  "unknown"
}

# Status (0 = censored) and times of a Surv / Hist object.
outcome_status <- function(y) {
  if (inherits(y, "Hist")) {
    if ("event" %in% colnames(y))                 # competing risks: event = cause
      return(ifelse(y[, "status"] == 0, 0, as.numeric(y[, "event"])))
    return(as.numeric(y[, "status"]))
  }
  as.numeric(y[, "status"])
}

outcome_time <- function(y) {
  as.numeric(y[, "time"])
}

# Generic consistency checks between the model's own outcome and `outcome`.
check_outcome_match <- function(model, data, info, horizon) {
  kind <- model_outcome_kind(model, data)
  if (kind == "unknown") return(invisible())
  cls <- class(model)[1]
  if (info$type == "binary") {
    if (kind %in% c("survival", "competing"))
      stop("This ", cls, " model was fitted to a survival outcome: use it with a ",
           "survival outcome made with Surv() and give `horizon`.", call. = FALSE)
    return(invisible())
  }
  if (kind == "binary")
    warning("The model (", cls, ") was fitted to a binary outcome, not a ",
            "survival outcome. Its predicted probabilities are taken as the ",
            "risk of the event by the horizon (", horizon, "), which can be ",
            "misleading if it was not built for this horizon and outcome.",
            call. = FALSE)
  if (info$type == "survival" && kind == "competing")
    stop("The model (", cls, ") is a competing-risks model, but the outcome ",
         "has a single event type. Give a competing-risks outcome: ",
         "Surv(time, status) with a factor status.", call. = FALSE)
  if (info$type == "competing" && kind == "survival")
    warning("The model (", cls, ") was fitted to a standard survival outcome: ",
            "its risk is used as the risk of '", info$states[info$k], "', but ",
            "it is not the cumulative incidence of that cause when there are ",
            "competing events, and can overstate the risk. Use a ",
            "competing-risks model.", call. = FALSE)
  if (kind %in% c("survival", "competing")) {
    tmax <- model_max_time(model)
    if (!is.null(tmax) && is.finite(tmax) && horizon > tmax)
      warning("`horizon` (", horizon, ") is beyond the last time in the ", cls,
              " model's own data (", signif(tmax, 4), "); its predictions there ",
              "are extrapolated or carry the last estimate forward.", call. = FALSE)
  }
  invisible()
}

# Longest follow-up in the model's own (training) data, if it can be read.
model_max_time <- function(model) {
  d <- tryCatch(suppressWarnings(suppressMessages(insight::get_data(model, verbose = FALSE))),
                error = function(e) NULL)
  y <- if (is.data.frame(d)) model_lhs(model, d) else NULL
  if (!inherits(y, c("Surv", "Hist"))) y <- stored_surv(model)
  if (!inherits(y, c("Surv", "Hist"))) return(NULL)
  max(outcome_time(y), na.rm = TRUE)
}

# A Surv / Hist object kept inside a model object (e.g. its outcome), found
# without knowing the class.
stored_surv <- function(model) {
  found <- NULL
  visit <- function(x, depth) {
    if (!is.null(found)) return(invisible())
    if (inherits(x, c("Surv", "Hist"))) { found <<- x; return(invisible()) }
    if (depth < 3 && is.list(x)) for (el in x) visit(el, depth + 1)
  }
  if (is.list(model)) visit(unclass(model), 1)
  found
}


# ---- methods for models the standard routes cannot handle -------------------
#
# Binary-outcome models. Each needs its own predict() call.

#' @rdname predict_risk
#' @export
predict_risk.svm <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is.null(horizon)) return(NextMethod())
  if (!isTRUE(model$compprob))
    stop("This svm was fitted without `probability = TRUE`, so it only ",
         "predicts classes. Refit it with `probability = TRUE` to get ",
         "predicted risks.", call. = FALSE)
  check_covariates(model, newdata)
  # predict.svm drops rows with a missing value in any column, used or not
  vars <- model_covariates(model)
  nd <- if (is.null(vars)) newdata else newdata[, vars, drop = FALSE]
  pr <- stats::predict(model, newdata = nd, probability = TRUE)
  finish_risk(attr(pr, "probabilities"), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.naiveBayes <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is.null(horizon)) return(NextMethod())
  vars <- names(model$tables)
  check_covariates(model, newdata, needed = vars)
  pr <- stats::predict(model, newdata = newdata[, vars, drop = FALSE], type = "raw")
  finish_risk(pr, newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.nnet <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is.null(horizon)) return(NextMethod())
  check_covariates(model, newdata)
  finish_risk(stats::predict(model, newdata = newdata, type = "raw"), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.lda <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is.null(horizon)) return(NextMethod())
  check_covariates(model, newdata)
  finish_risk(stats::predict(model, newdata = newdata)$posterior, newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.qda <- predict_risk.lda

#' @rdname predict_risk
#' @export
predict_risk.glmnet <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!inherits(model, "lognet")) return(NextMethod())
  if (length(model$lambda) != 1)
    stop("This glmnet model has ", length(model$lambda), " values of lambda. ",
         "Fit it with a single lambda, or use cv.glmnet().", call. = FALSE)
  x <- feature_matrix(model, newdata, rownames(model$beta))
  finish_risk(stats::predict(model, newx = x, type = "response"), newdata)
}

#' @rdname predict_risk
#' @method predict_risk cv.glmnet
#' @export
predict_risk.cv.glmnet <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!inherits(model$glmnet.fit, "lognet")) return(NextMethod())
  x <- feature_matrix(model, newdata, rownames(model$glmnet.fit$beta))
  finish_risk(stats::predict(model, newx = x, type = "response"), newdata)  # lambda.1se
}

#' @rdname predict_risk
#' @method predict_risk xgb.Booster
#' @export
predict_risk.xgb.Booster <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is.null(horizon)) return(NextMethod())
  if (inherits(model, "xgboost")) {
    # xgboost() (version 3+) predicts from a data frame through the standard
    # route; the prediction is P(second outcome level), so turn it round if the
    # levels were given as c(1, 0)
    p <- predict_risk.default(model, newdata)
    lev <- attr(model, "metadata")$y_levels
    if (length(lev) == 2 && lev[1] == "1") p <- 1 - p
    return(p)
  }
  feats <- tryCatch(stats::variable.names(model), error = function(e) NULL)
  if (is.null(feats)) feats <- model$feature_names              # xgboost < 2
  finish_risk(stats::predict(model, feature_matrix(model, newdata, feats)), newdata)
}

#' @rdname predict_risk
#' @method predict_risk lgb.Booster
#' @export
predict_risk.lgb.Booster <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is.null(horizon)) return(NextMethod())
  js <- model$dump_model()
  m <- regmatches(js, regexpr('"feature_names":\\[[^]]*\\]', js))
  feats <- if (length(m) == 0) NULL else
    gsub('^\\s*"|"\\s*$', "", strsplit(sub('^"feature_names":\\[(.*)\\]$', "\\1", m), ",")[[1]])
  finish_risk(stats::predict(model, feature_matrix(model, newdata, feats)), newdata)
}

#' @rdname predict_risk
#' @method predict_risk catboost.Model
#' @export
predict_risk.catboost.Model <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is.null(horizon)) return(NextMethod())
  if (!requireNamespace("catboost", quietly = TRUE))
    stop("The catboost package is needed to predict from this model.", call. = FALSE)
  feats <- rownames(model$feature_importances)
  check_covariates(model, newdata, needed = feats)
  pool <- catboost::catboost.load_pool(newdata[, feats, drop = FALSE])
  finish_risk(catboost::catboost.predict(model, pool, prediction_type = "Probability"),
              newdata)
}

# For models trained on a numeric matrix: the feature columns of `newdata`,
# in the model's order.
feature_matrix <- function(model, newdata, feats) {
  if (is.null(feats))
    stop("Could not read the feature names from this ", class(model)[1],
         " model. Train it on a matrix with column names.", call. = FALSE)
  check_covariates(model, newdata, needed = feats)
  x <- newdata[, feats, drop = FALSE]
  not_num <- feats[!vapply(x, is.numeric, logical(1))]
  if (length(not_num) > 0)
    stop("This ", class(model)[1], " model was trained on a numeric matrix, ",
         "but these columns of the validation data are not numeric: ",
         paste(not_num, collapse = ", "), ". Code them the same way as in ",
         "the development data (e.g. as dummy columns).", call. = FALSE)
  as.matrix(x)
}

# Survival models. Each gives S(t | x); the risk by the horizon is 1 - S.

#' @rdname predict_risk
#' @export
predict_risk.coxph <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (is.null(horizon)) return(NextMethod())
  # with the design matrix stored (x = TRUE, or rms::cph), riskRegression
  # works, except for case weights and beyond the last event time
  if (!is.null(model[["x"]])) {
    r <- try_route(predict_risk.default(model, newdata, horizon = horizon, cause = cause, ...))
    if (!inherits(r, "error")) return(r)
  }
  check_covariates(model, newdata)
  s <- surv_at(model, newdata, rep(horizon, nrow(newdata)))
  finish_risk(1 - s, newdata)
}

#' @rdname predict_risk
#' @method predict_risk ObliqueForestSurvival
#' @export
predict_risk.ObliqueForestSurvival <- function(model, newdata, horizon = NULL,
                                               cause = NULL, ...) {
  if (is.null(horizon)) return(NextMethod())
  check_covariates(model, newdata)
  finish_risk(1 - surv_at(model, newdata, rep(horizon, nrow(newdata))), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.constparty <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (is.null(horizon) || !inherits(model$fitted[["(response)"]], "Surv"))
    return(NextMethod())
  check_covariates(model, newdata)
  finish_risk(1 - surv_at(model, newdata, rep(horizon, nrow(newdata))), newdata)
}

#' @rdname predict_risk
#' @method predict_risk stpm2
#' @export
predict_risk.stpm2 <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (is.null(horizon)) return(NextMethod())
  check_covariates(model, newdata, needed = stpm2_covariates(model))
  finish_risk(1 - surv_at(model, newdata, rep(horizon, nrow(newdata))), newdata)
}

#' @rdname predict_risk
#' @method predict_risk pstpm2
#' @export
predict_risk.pstpm2 <- predict_risk.stpm2

#' @rdname predict_risk
#' @export
predict_risk.gam <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is_gam_cox(model)) return(NextMethod())
  if (is.null(horizon))
    stop("This gam model (Cox PH family) is a survival model: use it with a ",
         "survival outcome made with Surv() and give `horizon`.", call. = FALSE)
  check_covariates(model, newdata)
  finish_risk(1 - surv_at(model, newdata, rep(horizon, nrow(newdata))), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.flexsurvreg <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (is.null(horizon)) return(NextMethod())
  # riskRegression fails when the formula was passed as a variable
  r <- try_route(predict_risk.default(model, newdata, horizon = horizon, ...))
  if (!inherits(r, "error")) return(r)
  check_covariates(model, newdata)
  finish_risk(1 - surv_at(model, newdata, rep(horizon, nrow(newdata))), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.survreg <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  # riskRegression mismatches the scales of a stratified survreg
  if (is.null(horizon) || length(model$scale) == 1) return(NextMethod())
  check_covariates(model, newdata)
  finish_risk(1 - surv_at(model, newdata, rep(horizon, nrow(newdata))), newdata)
}

# bayespred.dev Bayesian survival models (Cox and parametric). A `bpmfit`
# survival fit gives risks only once a horizon is chosen: its posterior at the
# horizon is a `bpm` whose predict() is the posterior mean of the risk, the same
# convention as for bayespred's binary models.

#' @rdname predict_risk
#' @export
predict_risk.bpmfit <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (!is_bp_survival(model)) return(NextMethod())          # binary: standard route
  if (is.null(horizon))
    stop("This bpmfit model is a survival model: use it with a survival outcome ",
         "made with Surv() and give `horizon`.", call. = FALSE)
  check_covariates(model, newdata)
  post <- bp_dev("posterior")(model, time = horizon)
  finish_risk(bp_dev_predict(post, newdata), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.bpm <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (inherits(model, "bpm_lp"))
    stop("This bpm object has no baseline hazard, so it gives no risks. Package ",
         "the fit at the horizon: posterior(fit, time = ", if (is.null(horizon))
           "<horizon>" else horizon, ").", call. = FALSE)
  if (is.null(model$time)) return(NextMethod())             # binary: standard route
  if (is.null(horizon))
    stop("This bpm object is a survival model packaged at time = ", model$time,
         ": use it with a survival outcome made with Surv() and horizon = ",
         model$time, ".", call. = FALSE)
  if (abs(horizon - model$time) > sqrt(.Machine$double.eps) * max(1, horizon))
    stop("This bpm object gives risks at time = ", model$time, ", but `horizon` ",
         "is ", horizon, ". Use horizon = ", model$time, ", or package the fit at ",
         "the horizon: posterior(fit, time = ", horizon, ").", call. = FALSE)
  check_covariates(model, newdata)
  finish_risk(bp_dev_predict(model, newdata), newdata)
}

is_bp_survival <- function(model) {
  !is.null(model$cox) || !is.null(model$survreg)
}

# bayespred.dev's own functions (bayespred also defines bpmfit / bpm classes,
# so the methods are taken from bayespred.dev's namespace explicitly).
bp_dev <- function(name) {
  if (!requireNamespace("bayespred.dev", quietly = TRUE))
    stop("The bayespred.dev package is needed for this survival model.", call. = FALSE)
  getExportedValue("bayespred.dev", name)
}

bp_dev_predict <- function(post, newdata) {
  bp_dev("posterior")                                  # checks it is installed
  pred <- utils::getFromNamespace("predict.bpm", "bayespred.dev")
  as.vector(pred(post, newdata = newdata, type = "response"))
}

# Competing-risks models.

#' @rdname predict_risk
#' @export
predict_risk.crr <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (is.null(horizon) || is.null(cause))
    stop("This crr model is a competing-risks model: use it with a ",
         "competing-risks outcome, Surv(time, status) with a factor status, ",
         "and give `horizon`.", call. = FALSE)
  failcode <- model$call$failcode
  if (is.null(failcode)) failcode <- 1
  k <- list(...)$cause_index
  label <- if (is.null(cause)) "" else as.character(cause)
  ok <- identical(as.character(failcode), label) ||
    (!is.null(k) && identical(as.numeric(failcode), as.numeric(k)))
  if (!is.null(cause) && !ok)
    stop("This crr model was fitted for failcode ", failcode, ", but the net ",
         "benefit is for cause '", label, "'.", call. = FALSE)
  x <- feature_matrix(model, newdata, names(model$coef))
  pc <- stats::predict(model, cov1 = x)
  tt <- pc[, 1]
  i <- findInterval(horizon, tt)
  r <- if (i == 0) rep(0, nrow(x)) else pc[i, -1]
  finish_risk(as.vector(r), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.coxphms <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  # survival's multi-state Cox model, coxph(Surv(time, factor_status) ~ x, id = id).
  # riskRegression cannot predict from it; survival::survfit() gives the state
  # probabilities. stype = 1, ctype = 1: Aalen-Johansen with Breslow-type
  # hazard increments, the same estimator as riskRegression::CSC.
  if (is.null(horizon) || is.null(cause))
    stop("This coxph model is a multi-state (competing-risks) model: use it ",
         "with a competing-risks outcome, Surv(time, status) with a factor ",
         "status, and give `horizon`.", call. = FALSE)
  check_covariates(model, newdata)
  sf <- tryCatch(survival::survfit(model, newdata = newdata, ctype = 1, stype = 1),
                 error = function(e)
                   stop("Could not compute the cumulative incidences from this ",
                        "multi-state coxph model: ", conditionMessage(e), ". Fit it ",
                        "with `model = TRUE` (or keep its data available).", call. = FALSE))
  states <- sf$states
  events <- states[-1]                       # the first state is the initial one
  k <- list(...)$cause_index
  j <- if (as.character(cause) %in% events) match(as.character(cause), states)
       else if (!is.null(k) && k <= length(events)) k + 1
       else stop("Could not match the cause '", cause, "' to the states of the ",
                 "multi-state coxph model (", paste(events, collapse = ", "), ").",
                 call. = FALSE)
  ps <- summary(sf, times = horizon, extend = TRUE)$pstate
  ps <- array(ps, dim = c(1, nrow(newdata), length(states)))
  finish_risk(ps[1, , j], newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.FGR <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  # riskRegression returns the fitted cause's risk whatever cause is asked for,
  # so check the cause here before using the standard route.
  if (!is.null(cause)) {
    states <- as.character(attr(model$response, "states"))
    k <- list(...)$cause_index
    wanted <- if (as.character(cause) %in% states) as.character(cause)
              else if (!is.null(k) && k <= length(states)) states[k] else as.character(cause)
    if (!identical(wanted, as.character(model$cause)))
      stop("This FGR model was fitted for cause '", model$cause, "', but the ",
           "net benefit is for cause '", wanted, "'. Fit FGR with cause = '",
           wanted, "'.", call. = FALSE)
  }
  utils::capture.output(                        # predict.FGR prints covariate names
    r <- predict_risk.default(model, newdata, horizon = horizon, cause = cause, ...))
  r
}


# ---- survival probabilities S(t_i | x_i) at per-row times -------------------
#
# Used for the risk by the horizon (models above) and for the censoring
# probabilities of cIPCW. side = "left" gives the value just before t_i.

surv_at <- function(model, newdata, times, side = "right") UseMethod("surv_at")

# Standard route: riskRegression::predictRisk on a grid of the needed times.
#' @exportS3Method
surv_at.default <- function(model, newdata, times, side = "right") {
  eps <- if (side == "left") 1e-9 * max(1, max(times)) else 0
  out <- numeric(nrow(newdata))
  for (ch in split(seq_len(nrow(newdata)), ceiling(seq_len(nrow(newdata)) / 250))) {
    tg <- sort(unique(times[ch] - eps))
    r <- try_route(riskRegression::predictRisk(model, newdata = newdata[ch, , drop = FALSE],
                                               times = tg))
    if (inherits(r, "error") && !is.null(vars <- model_covariates(model)))
      r <- try_route(riskRegression::predictRisk(
        model, newdata = newdata[ch, intersect(vars, names(newdata)), drop = FALSE], times = tg))
    if (inherits(r, "error"))
      stop("Could not get survival probabilities from this ", class(model)[1],
           " model (", trimws(conditionMessage(r)), ").", call. = FALSE)
    n <- length(ch)
    if (length(r) == length(tg) && n > 1)            # a model without covariates
      r <- matrix(rep(as.vector(r), each = n), nrow = n)
    if (length(r) != n * length(tg))
      stop("Unexpected output from riskRegression::predictRisk for this ",
           class(model)[1], " model.", call. = FALSE)
    r <- matrix(r, nrow = n)
    out[ch] <- 1 - r[cbind(seq_len(n), match(times[ch] - eps, tg))]
  }
  out
}

# coxph without the stored design matrix: S(t | x) = exp(-H0(t) exp(lp)) from
# what the model object keeps (outcome, linear predictors, weights).
#' @exportS3Method
surv_at.coxph <- function(model, newdata, times, side = "right") {
  if (!is.null(model[["x"]])) {
    r <- try_route(surv_at.default(model, newdata, times, side))
    if (!inherits(r, "error") && !anyNA(r)) return(r)
  }
  if (cox_simple(model)) {
    bh <- cox_basehaz(model)
    lp <- stats::predict(model, newdata = newdata, type = "lp")
    i <- if (side == "left") findInterval(times, bh$time, left.open = TRUE)
         else findInterval(times, bh$time)
    h <- ifelse(i == 0, 0, bh$hazard[pmax(i, 1)])
    return(unname(exp(-h * exp(lp))))
  }
  # stratified / time-dependent: survfit, which needs the model's data
  sf <- tryCatch(survival::survfit(model, newdata = newdata),
                 error = function(e)
                   stop("Could not compute survival from this coxph model: ",
                        conditionMessage(e), ". Fit it with `x = TRUE` (or ",
                        "`model = TRUE`).", call. = FALSE))
  vapply(seq_len(nrow(newdata)), function(r) {
    sr <- sf[r]
    i <- findInterval(times[r], sr$time, left.open = side == "left")
    if (i == 0) 1 else sr$surv[i]
  }, numeric(1))
}

# An unstratified Cox model without time-dependent terms that keeps its outcome.
cox_simple <- function(model) {
  !is.null(model$y) && is.null(model$strata) &&
    is.null(attr(stats::terms(model), "specials")$strata) &&
    is.null(attr(stats::terms(model), "specials")$tt) &&
    !is.null(model$linear.predictors) && ncol(model$y) == 2
}

# Centered baseline cumulative hazard (Efron or Breslow, as survival::basehaz).
cox_basehaz <- function(model) {
  y <- model$y
  tt <- y[, 1]; ev <- y[, 2]
  w <- if (is.null(model$weights)) rep(1, length(tt)) else model$weights
  r <- w * exp(model$linear.predictors)
  ut <- sort(unique(tt[ev == 1]))
  efron <- identical(model$method, "efron")
  o <- order(tt, decreasing = TRUE)
  risk_sum <- cumsum(r[o])
  at_risk <- function(t) risk_sum[sum(tt[o] >= t)]
  inc <- vapply(ut, function(t) {
    dd <- ev == 1 & tt == t
    R <- at_risk(t)
    if (!efron) return(sum(w[dd]) / R)
    D <- sum(r[dd]); m <- sum(dd)
    sum((sum(w[dd]) / m) / (R - (0:(m - 1)) / m * D))
  }, numeric(1))
  data.frame(hazard = cumsum(inc), time = ut)
}

#' @exportS3Method
surv_at.survreg <- function(model, newdata, times, side = "right") {
  if (length(model$scale) == 1) return(NextMethod())
  # one scale per stratum: S(t) = 1 - F((log t - lp) / scale_stratum)
  tt <- stats::terms(model)
  sp <- attr(tt, "specials")$strata
  st_call <- attr(tt, "variables")[[1 + sp[1]]]
  env <- new.env(parent = environment(tt))
  env$strata <- survival::strata
  lab <- as.character(eval(st_call, newdata, env))
  sc <- model$scale[lab]
  if (anyNA(sc))
    stop("Some strata of the validation data are not among the survreg ",
         "model's strata.", call. = FALSE)
  eps <- if (side == "left") 1e-9 * max(1, max(times)) else 0
  lp <- stats::predict(model, newdata = newdata, type = "lp")
  unname(1 - survival::psurvreg(times - eps, mean = lp, scale = sc,
                                distribution = model$dist))
}

#' @exportS3Method
surv_at.flexsurvreg <- function(model, newdata, times, side = "right") {
  r <- try_route(surv_at.default(model, newdata, times, side))
  if (!inherits(r, "error")) return(r)
  eps <- if (side == "left") 1e-9 * max(1, max(times)) else 0
  out <- numeric(nrow(newdata))
  for (tt in unique(times)) {
    rows <- which(times == tt)
    out[rows] <- summary(model, newdata = newdata[rows, , drop = FALSE], type = "survival",
                         t = tt - eps, ci = FALSE, tidy = TRUE)$est
  }
  out
}

#' @exportS3Method
surv_at.bpmfit <- function(model, newdata, times, side = "right") {
  # bayespred.dev survival fit as a censoring model: its posterior at each
  # needed time (posterior mean of the survival probability)
  if (!is_bp_survival(model)) return(NextMethod())
  eps <- if (side == "left") 1e-9 * max(1, max(times)) else 0
  first_event <- if (!is.null(model$cox$time))
    min(model$cox$time[model$cox$status == 1]) else -Inf
  out <- numeric(nrow(newdata))
  for (tt in unique(times)) {
    rows <- which(times == tt)
    t_use <- tt - eps
    out[rows] <- if (t_use < first_event) 1 else
      1 - bp_dev_predict(bp_dev("posterior")(model, time = t_use),
                         newdata[rows, , drop = FALSE])
  }
  out
}

#' @exportS3Method
surv_at.ObliqueForestSurvival <- function(model, newdata, times, side = "right") {
  eps <- if (side == "left") 1e-9 * max(1, max(times)) else 0
  out <- numeric(nrow(newdata))
  for (ch in split(seq_len(nrow(newdata)), ceiling(seq_len(nrow(newdata)) / 250))) {
    tg <- sort(unique(times[ch] - eps))
    s <- matrix(stats::predict(model, new_data = newdata[ch, , drop = FALSE],
                               pred_type = "surv", pred_horizon = tg), nrow = length(ch))
    out[ch] <- s[cbind(seq_along(ch), match(times[ch] - eps, tg))]
  }
  out
}

# partykit survival trees: Kaplan-Meier of the training data in each node.
#' @exportS3Method
surv_at.constparty <- function(model, newdata, times, side = "right") {
  resp <- model$fitted[["(response)"]]
  if (!inherits(resp, "Surv")) return(NextMethod())
  node_train <- model$fitted[["(fitted)"]]
  wt <- model$fitted[["(weights)"]]
  if (is.null(wt)) wt <- rep(1, length(node_train))
  node_new <- stats::predict(model, newdata = newdata, type = "node")
  out <- numeric(length(node_new))
  for (nd in unique(node_new)) {
    sel <- node_train == nd
    sf <- survival::survfit(resp[sel] ~ 1, weights = wt[sel])
    rows <- node_new == nd
    i <- if (side == "left") findInterval(times[rows], sf$time, left.open = TRUE)
         else findInterval(times[rows], sf$time)
    out[rows] <- ifelse(i == 0, 1, sf$surv[pmax(i, 1)])
  }
  out
}

# rstpm2 flexible parametric models (predict needs rstpm2 attached).
#' @exportS3Method
surv_at.stpm2 <- function(model, newdata, times, side = "right") {
  if (!"package:rstpm2" %in% search())
    stop("rstpm2 must be attached to predict from this ", class(model)[1],
         " model: run library(rstpm2) first.", call. = FALSE)
  eps <- if (side == "left") 1e-9 * max(1, max(times)) else 0
  nd <- newdata
  nd[[model@timeVar]] <- times - eps
  pred <- getExportedValue("rstpm2", "predict")       # an S4 generic in rstpm2
  as.vector(pred(model, newdata = nd, type = "surv"))
}
#' @exportS3Method
surv_at.pstpm2 <- surv_at.stpm2

stpm2_covariates <- function(model) {
  setdiff(all.vars(model@call.formula[[3]]), model@timeVar)
}

# mgcv's Cox PH family: predict(type = "response") is the survival at the
# time in the response column of newdata.
#' @exportS3Method
surv_at.gam <- function(model, newdata, times, side = "right") {
  if (!is_gam_cox(model)) return(NextMethod())
  eps <- if (side == "left") 1e-9 * max(1, max(times)) else 0
  nd <- newdata
  nd[[all.vars(stats::formula(model)[[2]])[1]]] <- times - eps
  as.vector(stats::predict(model, newdata = nd, type = "response"))
}

is_gam_cox <- function(model) {
  inherits(model, "gam") && identical(model$family$family, "Cox PH")
}
