#' Predicted risks from a fitted model
#'
#' Turns a fitted model into a numeric vector of predicted risks
#' P(Y = 1 | x), one per row of `newdata`.
#'
#' Methods are provided for `glm`, `ranger`, `randomForest`, and bayespred
#' (`bpmfit` and `bpm`) objects. The default method tries
#' `predict(model, newdata, type = "response")` and then
#' `predict(model, newdata, type = "prob")`. A vector of predicted risks is not
#' accepted: voipred needs the fitted model itself.
#'
#' The event is always Y = 1. Random forests fitted to a factor outcome
#' predict one probability per outcome level; the risk is the probability of
#' level `"1"` (or of the second level if the levels are not 0/1).
#'
#' Before predicting, the covariates the model was fitted with are compared
#' with the columns of `newdata`. If any are missing, an error lists them.
#' This check is skipped when the covariate names cannot be read from the
#' model.
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
         class(model)[1], "'. Supported models: glm, ranger, randomForest, ",
         "bayespred (bpmfit, bpm), or any model whose predict() method ",
         "returns probabilities with type = \"response\" or \"prob\".",
         call. = FALSE)
  check_risk(as.vector(p), newdata)
}


# ---- internal helpers -------------------------------------------------------

# Names of the covariates a fitted model needs, or NULL if they can't be read.
# The outcome is excluded: it is not needed to predict.
model_covariates <- function(model) {
  if (inherits(model, "bpmfit")) return(model_covariates(model$posterior))
  if (inherits(model, "ranger")) return(model$forest$independent.variable.names)
  vars <- tryCatch(all.vars(stats::delete.response(stats::terms(model))),
                   error = function(e) NULL)
  if (is.null(vars) && inherits(model, "randomForest"))
    vars <- rownames(model$importance)   # x/y interface: no terms stored
  vars
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

# From a matrix with one probability column per outcome level, take the risk
# of the event: column "1" if present, otherwise the second of two columns.
event_column <- function(pr) {
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
    stop("Predicted risks must lie in [0, 1].", call. = FALSE)
  if (all(p == 0 | p == 1))
    stop("The model returned only 0s and 1s. These look like predicted ",
         "classes, not probabilities. Make the model return probabilities ",
         "(e.g. fit ranger or svm with `probability = TRUE`).", call. = FALSE)
  p
}
