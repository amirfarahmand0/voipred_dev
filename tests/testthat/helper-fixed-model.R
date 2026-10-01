# A test-only "model" whose predict() returns risks we choose. It lets the
# net benefit tests use hand-picked risks while nb() still receives a model,
# not a vector. It goes through predict_risk.default like any unknown model.
fixed_model <- function(p) structure(list(p = p), class = "fixed_risk")

predict.fixed_risk <- function(object, newdata, type = "response", ...) object$p

registerS3method("predict", "fixed_risk", predict.fixed_risk)
