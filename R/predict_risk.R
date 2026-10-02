#' Predicted risks from a fitted model
#'
#' Turns a fitted model into a numeric vector of predicted risks
#' P(Y = 1 | x), one per row of `newdata`.
#'
#' Any model whose `predict()` returns probabilities with
#' `type = "response"` or `type = "prob"` works through the default method
#' (e.g. `glm`, `mgcv::gam`, `logistf`, `rpart`, `partykit::ctree`, `gbm`,
#' `kknn` (attach it with `library(kknn)`), `earth` with a binomial glm). Methods are provided for models that return
#' probabilities in another way:
#'
#' * `ranger` (fit with `probability = TRUE`), `randomForest`;
#' * bayespred (`bpmfit`, `bpm`);
#' * `e1071::svm` (fit with `probability = TRUE`), `e1071::naiveBayes`,
#'   `naivebayes::naive_bayes`, `nnet`, `MASS::lda`, `MASS::qda`;
#' * models fitted to a numeric matrix: `glmnet` (one lambda) and
#'   `cv.glmnet` (at `lambda.1se`, glmnet's default), `xgboost` (objective
#'   `binary:logistic`), `lightgbm` (objective `binary`) and `catboost`
#'   (loss `Logloss`). For these, `newdata` must contain columns with the same
#'   names as the features the model was trained on (e.g. the dummy columns).
#'
#' A vector of predicted risks is not accepted: voipred needs the fitted model
#' itself. The event is always Y = 1: when a model predicts one probability per
#' outcome level, the risk is the probability of level `"1"` (or of the second
#' level if the levels are not 0/1).
#'
#' Before predicting, the covariates the model was fitted with are compared
#' with the columns of `newdata`. If any are missing, an error lists them.
#' The result must be probabilities: values outside \[0, 1\] (e.g. log-odds or
#' other scores) and predictions that are only 0s and 1s (classes) are
#' rejected.
#'
#' @param model A fitted model.
#' @param newdata A data frame with the covariates the model needs.
#' @param bayespred_method Only for bayespred models: how the risk is
#'   predicted. `"pm"` (default) is the posterior mean of risk; `"pe"` uses the
#'   point estimate of the coefficients; `"pm_mackay"` is MacKay's
#'   approximation to the posterior mean.
#' @param ... Not used; for methods.
#' @return A numeric vector of length `nrow(newdata)` with values in \[0, 1\].
#' @examples
#' fit <- glm(day30 ~ age + Killip, family = binomial, data = gusto[1:2000, ])
#' p <- predict_risk(fit, gusto[2001:3000, ])
#' summary(p)
#' @export
predict_risk <- function(model, newdata, ...) {
  UseMethod("predict_risk")
}

#' @rdname predict_risk
#' @export
predict_risk.glm <- function(model, newdata, ...) {
  fam <- stats::family(model)$family
  if (!fam %in% c("binomial", "quasibinomial"))
    stop("The glm must have a binomial family; this one is '", fam, "'.",
         call. = FALSE)
  check_covariates(model, newdata)
  p <- stats::predict(model, newdata = newdata, type = "response")
  check_risk(as.vector(p), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.ranger <- function(model, newdata, ...) {
  if (identical(model$treetype, "Classification"))
    stop("This ranger forest only predicts classes. Refit it with ",
         "`probability = TRUE` to get predicted risks.", call. = FALSE)
  check_covariates(model, newdata)
  pr <- stats::predict(model, data = newdata)$predictions
  if (is.matrix(pr)) pr <- event_column(pr)
  check_risk(as.vector(pr), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.randomForest <- function(model, newdata, ...) {
  check_covariates(model, newdata)
  if (identical(model$type, "classification")) {
    pr <- stats::predict(model, newdata = newdata, type = "prob")
    pr <- event_column(pr)
  } else {
    pr <- stats::predict(model, newdata = newdata)
  }
  check_risk(as.vector(pr), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.bpmfit <- function(model, newdata, bayespred_method = "pm", ...) {
  bayespred_method <- match.arg(bayespred_method, c("pm", "pe", "pm_mackay"))
  check_covariates(model, newdata)
  p <- stats::predict(model, newdata = newdata, type = "response",
                      method = bayespred_method)
  check_risk(as.vector(p), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.bpm <- predict_risk.bpmfit

#' @rdname predict_risk
#' @export
predict_risk.svm <- function(model, newdata, ...) {
  if (!isTRUE(model$compprob))
    stop("This svm was fitted without `probability = TRUE`, so it only ",
         "predicts classes. Refit it with `probability = TRUE` to get ",
         "predicted risks.", call. = FALSE)
  check_covariates(model, newdata)
  # predict.svm drops rows with a missing value in any column, used or not
  pr <- stats::predict(model, newdata = model_columns(model, newdata),
                       probability = TRUE)
  check_risk(as.vector(event_column(attr(pr, "probabilities"))), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.naiveBayes <- function(model, newdata, ...) {
  check_covariates(model, newdata)
  pr <- stats::predict(model, newdata = model_columns(model, newdata), type = "raw")
  check_risk(as.vector(event_column(pr)), newdata)
}

#' @rdname predict_risk
#' @method predict_risk naive_bayes
#' @export
predict_risk.naive_bayes <- function(model, newdata, ...) {
  check_covariates(model, newdata)
  pr <- stats::predict(model, newdata = model_columns(model, newdata), type = "prob")
  check_risk(as.vector(event_column(pr)), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.nnet <- function(model, newdata, ...) {
  check_covariates(model, newdata)
  pr <- stats::predict(model, newdata = newdata, type = "raw")
  if (is.matrix(pr)) pr <- event_column(pr)
  check_risk(as.vector(pr), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.lda <- function(model, newdata, ...) {
  check_covariates(model, newdata)
  pr <- stats::predict(model, newdata = newdata)$posterior
  check_risk(as.vector(event_column(pr)), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.qda <- predict_risk.lda

#' @rdname predict_risk
#' @export
predict_risk.glmnet <- function(model, newdata, ...) {
  if (!inherits(model, "lognet"))
    stop("The glmnet model must be fitted with family = \"binomial\".",
         call. = FALSE)
  if (length(model$lambda) != 1)
    stop("This glmnet model has ", length(model$lambda), " values of lambda. ",
         "Fit it with a single lambda, or use cv.glmnet().", call. = FALSE)
  x <- feature_matrix(model, newdata)
  p <- stats::predict(model, newx = x, type = "response")
  check_risk(as.vector(p), newdata)
}

#' @rdname predict_risk
#' @method predict_risk cv.glmnet
#' @export
predict_risk.cv.glmnet <- function(model, newdata, ...) {
  if (!inherits(model$glmnet.fit, "lognet"))
    stop("The glmnet model must be fitted with family = \"binomial\".",
         call. = FALSE)
  x <- feature_matrix(model, newdata)
  p <- stats::predict(model, newx = x, type = "response")   # s = "lambda.1se"
  check_risk(as.vector(p), newdata)
}

#' @rdname predict_risk
#' @method predict_risk xgb.Booster
#' @export
predict_risk.xgb.Booster <- function(model, newdata, ...) {
  lev <- attr(model, "metadata")$y_levels
  if (inherits(model, "xgboost")) {
    # xgboost() (version 3+) takes a data frame, factors included, and
    # predicts the probability of the second outcome level
    check_covariates(model, newdata)
    p <- stats::predict(model, model_columns(model, newdata))
    if (length(lev) == 2 && lev[1] == "1") p <- 1 - p
  } else {
    p <- stats::predict(model, feature_matrix(model, newdata))
  }
  check_risk(as.vector(p), newdata)
}

#' @rdname predict_risk
#' @method predict_risk lgb.Booster
#' @export
predict_risk.lgb.Booster <- function(model, newdata, ...) {
  x <- feature_matrix(model, newdata)
  check_risk(as.vector(stats::predict(model, x)), newdata)
}

#' @rdname predict_risk
#' @method predict_risk catboost.Model
#' @export
predict_risk.catboost.Model <- function(model, newdata, ...) {
  if (!requireNamespace("catboost", quietly = TRUE))
    stop("The catboost package is needed to predict from this model.",
         call. = FALSE)
  check_covariates(model, newdata)
  pool <- catboost::catboost.load_pool(newdata[, model_covariates(model), drop = FALSE])
  p <- catboost::catboost.predict(model, pool, prediction_type = "Probability")
  check_risk(as.vector(p), newdata)
}

#' @rdname predict_risk
#' @export
predict_risk.default <- function(model, newdata, ...) {
  if (is.numeric(model) || is.logical(model))
    stop("`model` must be a fitted model, not a vector of predicted risks. ",
         "voipred needs the model itself (e.g. a glm, ranger, randomForest ",
         "or bayespred fit).", call. = FALSE)
  check_covariates(model, newdata)
  p <- tryCatch(stats::predict(model, newdata = newdata, type = "response"),
                error = function(e) NULL)
  if (is.null(p) || !is.numeric(p)) {
    p <- tryCatch(stats::predict(model, newdata = newdata, type = "prob"),
                  error = function(e) NULL)
  }
  if (is.matrix(p) || is.data.frame(p)) p <- event_column(as.matrix(p))
  if (is.null(p) || !is.numeric(p))
    stop("Could not get predicted risks from an object of class '",
         class(model)[1], "'. See ?predict_risk for the supported models; ",
         "others work if their predict() returns probabilities with ",
         "type = \"response\" or \"prob\".", call. = FALSE)
  check_risk(as.vector(p), newdata)
}


# ---- internal helpers -------------------------------------------------------

# Names of the covariates a fitted model needs, or NULL if they can't be read.
# The outcome is excluded: it is not needed to predict.
model_covariates <- function(model) {
  if (inherits(model, "bpmfit")) return(model_covariates(model$posterior))
  if (inherits(model, "ranger")) return(model$forest$independent.variable.names)
  if (inherits(model, c("naiveBayes", "naive_bayes"))) return(names(model$tables))
  if (inherits(model, "cv.glmnet")) return(rownames(model$glmnet.fit$beta))
  if (inherits(model, "glmnet")) return(rownames(model$beta))
  if (inherits(model, "xgb.Booster")) return(xgb_feature_names(model))
  if (inherits(model, "lgb.Booster")) return(lgb_feature_names(model))
  if (inherits(model, "catboost.Model")) return(rownames(model$feature_importances))
  vars <- tryCatch(all.vars(stats::delete.response(stats::terms(model))),
                   error = function(e) NULL)
  if (is.null(vars) && inherits(model, "randomForest"))
    vars <- rownames(model$importance)   # x/y interface: no terms stored
  vars
}

xgb_feature_names <- function(model) {
  fn <- tryCatch(model$feature_names, error = function(e) NULL)       # xgboost < 2
  if (is.null(fn))
    fn <- tryCatch(stats::variable.names(model), error = function(e) NULL)
  fn
}

lgb_feature_names <- function(model) {
  js <- model$dump_model()
  m <- regmatches(js, regexpr('"feature_names":\\[[^]]*\\]', js))
  if (length(m) == 0) return(NULL)
  inner <- sub('^"feature_names":\\[(.*)\\]$', "\\1", m)
  gsub('^\\s*"|"\\s*$', "", strsplit(inner, ",")[[1]])
}

check_covariates <- function(model, newdata) {
  if (!is.data.frame(newdata))
    stop("`data` must be a data frame.", call. = FALSE)
  needed <- model_covariates(model)
  if (is.null(needed)) return(invisible(TRUE))
  missing <- setdiff(needed, names(newdata))
  if (length(missing) > 0)
    stop("The validation data is missing ", length(missing),
         if (length(missing) == 1) " covariate" else " covariates",
         " that the model was developed with: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  invisible(TRUE)
}

# Only the columns the model was fitted with (some predict() methods drop rows
# with missing values in unused columns, or warn about extra columns).
model_columns <- function(model, newdata) {
  vars <- model_covariates(model)
  if (is.null(vars)) newdata else newdata[, vars, drop = FALSE]
}

# For models fitted to a numeric matrix: take the model's feature columns from
# the validation data, in the model's order.
feature_matrix <- function(model, newdata) {
  check_covariates(model, newdata)
  feats <- model_covariates(model)
  if (is.null(feats))
    stop("Could not read the feature names from this ", class(model)[1],
         " model. Train it on a matrix with column names.", call. = FALSE)
  x <- newdata[, feats, drop = FALSE]
  not_num <- feats[!vapply(x, is.numeric, logical(1))]
  if (length(not_num) > 0)
    stop("This ", class(model)[1], " model was trained on a numeric matrix, ",
         "but these columns of the validation data are not numeric: ",
         paste(not_num, collapse = ", "), ". Code them the same way as in ",
         "the development data (e.g. as dummy columns).", call. = FALSE)
  as.matrix(x)
}

# From a matrix with one probability column per outcome level, take the risk
# of the event: column "1" if present, otherwise the second of two columns.
# A single column is taken as the risk itself.
event_column <- function(pr) {
  if (ncol(pr) == 1) return(pr[, 1])
  if (!is.null(colnames(pr)) && "1" %in% colnames(pr)) return(pr[, "1"])
  if (ncol(pr) != 2)
    stop("The model predicts ", ncol(pr), " outcome classes; the outcome ",
         "must be binary.", call. = FALSE)
  pr[, 2]
}

check_risk <- function(p, newdata = NULL) {
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
         "probabilities. They may be wrong, or scores from the model (e.g. ",
         "log-odds or decision values) that must be turned into probabilities ",
         "first.", call. = FALSE)
  if (all(p == 0 | p == 1))
    stop("The model returned only 0s and 1s. These look like predicted ",
         "classes, not probabilities. Make the model return probabilities ",
         "(e.g. fit ranger or svm with `probability = TRUE`).", call. = FALSE)
  p
}
