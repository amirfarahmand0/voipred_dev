# A test-only "model" whose predict() returns risks we choose. It lets the
# net benefit tests use hand-picked risks while nb() still receives a model,
# not a vector. It goes through predict_risk.default like any unknown model.
fixed_model <- function(p) structure(list(p = p), class = "fixed_risk")

predict.fixed_risk <- function(object, newdata, type = "response", ...) object$p

registerS3method("predict", "fixed_risk", predict.fixed_risk)

# For survival outcomes the standard route (riskRegression::predictRisk) does
# not know this test class, so it gets an (internal) predict_risk() method.
predict_risk.fixed_risk <- function(model, newdata, horizon = NULL, cause = NULL, ...) {
  if (is.null(horizon)) return(NextMethod())
  model$p
}
registerS3method("predict_risk", "fixed_risk", predict_risk.fixed_risk,
                 envir = asNamespace("voipred.dev"))
