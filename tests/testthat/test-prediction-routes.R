# How risks are obtained from a model: predict_function > voipred's own
# (internal) method for the class > the standard routes; output shapes read
# generically; the same checks for every model. Uses small test-only classes.

toy <- function(cls, ...) structure(list(...), class = cls)
reg_predict <- function(cls, f) registerS3method("predict", cls, f)
reg_method <- function(cls, f) registerS3method("predict_risk", cls, f,
                                                envir = asNamespace("voipred.dev"))
reg_predictrisk <- function(cls, f) registerS3method("predictRisk", cls, f,
                                                     envir = asNamespace("riskRegression"))

d <- data.frame(y = c(1, 0, 1, 0, 0, 1), x = 1:6)
p0 <- c(0.8, 0.2, 0.6, 0.3, 0.1, 0.7)

test_that("predict_function takes precedence over voipred's own method and the standard route", {
  reg_predict("toyA", function(object, newdata, type = "response", ...) rep(0.2, nrow(newdata)))
  m <- toy("toyA")
  expect_equal(nb(d, m, "y", thresholds = 0.1)$risk, rep(0.2, 6))      # standard
  reg_method("toyA", function(model, newdata, ...) rep(0.3, nrow(newdata)))
  expect_equal(nb(d, m, "y", thresholds = 0.1)$risk, rep(0.3, 6))      # method
  expect_equal(nb(d, m, "y", thresholds = 0.1,
                       predict_function = function(model, newdata, ...) p0)$risk, p0)
})

test_that("predict_function gets horizon and cause; it is checked like any model", {
  seen <- NULL
  f <- function(model, newdata, horizon, cause) { seen <<- list(horizon, cause); p0 }
  nb(d, toy("toyB"), "y", thresholds = 0.1, predict_function = f)
  expect_null(seen[[1]]); expect_null(seen[[2]])
  ds <- data.frame(t = c(1, 2, 6, 7, 8, 9), st = factor(c(1, 2, 0, 1, 0, 2), 0:2,
                                                        c("cens", "a", "b")))
  nb(ds, toy("toyB"), ~ survival::Surv(t, st), horizon = 5, cause = "b", thresholds = 0.1,
     predict_function = f)
  expect_equal(seen, list(5, "b"))
  m <- toy("toyB")
  expect_error(nb(d, m, "y", predict_function = "p"), "must be a function")
  expect_error(nb(d, m, "y", predict_function = function(...) p0[-1]),
               "returned 5 predictions for 6 rows")
  expect_error(nb(d, m, "y", predict_function = function(...) replace(p0, 2, NA)), "1 predicted risks are NA")
  expect_error(nb(d, m, "y", predict_function = function(...) qlogis(p0)),
               "outside \\[0, 1\\].*not probabilities")
  expect_error(nb(d, m, "y", predict_function = function(...) c(1, 0, 1, 0, 0, 1)),
               "only 0s and 1s.*classes \\(labels\\)")
  expect_error(nb(d, m, "y", predict_function = function(...) p0 > 0.5),
               "only 0s and 1s.*classes \\(labels\\)")
  expect_error(nb(d, m, "y", predict_function = function(...) factor(p0 > 0.5)),
               "returned classes, not probabilities")
  expect_error(nb(d, m, "y", predict_function = function(...) as.character(p0)),
               "returned classes, not probabilities")
  expect_error(nb(d, p0, "y", predict_function = function(...) p0), "not used with `predicted_risk`")
})

test_that("output shapes are read generically", {
  m <- toy("toyC")
  r <- function(out, ...) nb(d, m, "y", thresholds = 0.1,
                                  predict_function = function(...) out, ...)$risk
  two <- cbind(`0` = 1 - p0, `1` = p0)
  expect_equal(r(two), p0)                                     # column "1"
  expect_equal(r(two[, 2:1]), p0)                              # wherever it is
  expect_equal(r(data.frame(.pred_0 = 1 - p0, .pred_1 = p0)), p0)  # tidymodels
  expect_equal(r(unname(two)), p0)                             # second of two
  expect_equal(r(matrix(p0, ncol = 1)), p0)                    # one column
  expect_equal(r(setNames(p0, letters[1:6])), p0)              # names dropped
  expect_error(r(cbind(a = p0 / 2, b = p0 / 2, c = 1 - p0)), "predicts 3 outcome classes")
  expect_error(r(list(time = 1:3, surv = matrix(0.5, 6, 3))),
               "survival curve was returned, but there is no `horizon`")
  # survival curve on a grid: the value at the last grid time <= horizon
  ds <- data.frame(t = c(1, 2, 6, 7, 8, 9), st = c(1, 1, 0, 1, 0, 1))
  surv <- outer(1 - p0, c(1, 3, 6), function(a, tt) a^(tt / 3))
  res <- nb(ds, m, ~ survival::Surv(t, st), horizon = 5, thresholds = 0.1,
            predict_function = function(...) list(time = c(1, 3, 6), surv = surv))
  expect_equal(res$risk, 1 - surv[, 2])
  # cumulative incidence curve (competing risks)
  dc <- data.frame(t = c(1, 2, 6, 7, 8, 9),
                   st = factor(c(1, 2, 0, 1, 0, 2), 0:2, c("cens", "a", "b")))
  cif <- outer(p0, c(1, 3, 6), function(a, tt) a * tt / 6)
  res <- nb(dc, m, ~ survival::Surv(t, st), horizon = 5, thresholds = 0.1,
            predict_function = function(...) list(time = c(1, 3, 6), cif = cif))
  expect_equal(res$risk, cif[, 2])
})

test_that("standard binary route: type = 'response', then 'prob', then riskRegression", {
  reg_predict("toyD", function(object, newdata, type = "response", ...) {
    if (type == "response") stop("no response type")
    cbind(`0` = 1 - p0, `1` = p0)
  })
  expect_equal(predict_risk(toy("toyD"), d), p0)
  reg_predict("toyE", function(object, newdata, type = "response", ...)
    factor(p0 > 0.5))                         # classes for every type
  reg_predictrisk("toyE", function(object, newdata, ...) p0)
  expect_equal(predict_risk(toy("toyE"), d), p0)
})

test_that("standard survival route: riskRegression::predictRisk, cause by label or position", {
  ds <- data.frame(t = c(1, 2, 6, 7, 8, 9),
                   st = factor(c(1, 2, 0, 1, 0, 2), 0:2, c("cens", "a", "b")))
  # a model whose predictRisk() only knows numeric causes (like rfsrc)
  reg_predictrisk("toyF", function(object, newdata, times, cause, ...) {
    if (!is.numeric(cause)) stop("cause must be a number")
    rep(c(0.2, 0.6)[cause], nrow(newdata)) + seq(0, 0.05, length.out = nrow(newdata))
  })
  r_b <- nb(ds, toy("toyF"), ~ survival::Surv(t, st), horizon = 5, cause = "b",
                 thresholds = 0.1)$risk
  expect_equal(r_b, 0.6 + seq(0, 0.05, length.out = 6))
  # a model whose predictRisk() knows the labels (like CSC fitted on a factor)
  reg_predictrisk("toyG", function(object, newdata, times, cause, ...) {
    rep(c(a = 0.3, b = 0.5)[[cause]], nrow(newdata)) + seq(0, 0.05, length.out = nrow(newdata))
  })
  r_a <- nb(ds, toy("toyG"), ~ survival::Surv(t, st), horizon = 5, cause = "a",
                 thresholds = 0.1)$risk
  expect_equal(r_a, 0.3 + seq(0, 0.05, length.out = 6))
})

test_that("a model with no usable route gets one clear error with the ways out", {
  d2 <- data.frame(t = c(1, 2, 6, 7, 8, 9), s = c(1, 1, 0, 1, 0, 1))
  expect_error(nb(d, toy("toyH"), "y"),
               "Could not get predicted risks from this toyH model.*give `predict_function")
  expect_error(nb(d, toy("toyH"), "y"), "returning the risk of the event\\.$")
  expect_error(nb(d2, toy("toyH"), ~ survival::Surv(t, s), horizon = 5),
               "Could not get predicted risks from this toyH model.*by the horizon")
})

test_that("a new survival model works with predict_function", {
  s <- rott_split()
  # a model from a hypothetical new package: its predict() returns a survival
  # curve on a time grid
  fit_new <- function(formula, data)
    structure(list(cox = survival::coxph(formula, data, x = TRUE)), class = "newsurv")
  reg_predict("newsurv", function(object, newdata, ...) {
    sf <- survival::survfit(object$cox, newdata = newdata)
    list(times = sf$time, survival = t(sf$surv))
  })
  m <- fit_new(surv_formula, s$dev)
  expect_error(nb(s$val, m, surv_outcome, horizon = 5), "Could not get predicted risks")
  ref <- predict_risk(m$cox, s$val, horizon = 5)
  res <- nb(s$val, m, surv_outcome, horizon = 5, thresholds = c(0.1, 0.3),
            predict_function = function(model, newdata, ...) {
              pr <- predict(model, newdata); list(time = pr$times, surv = pr$survival) })
  expect_equal(res$risk, ref, tolerance = 1e-12)
  # the same as giving the risks
  expect_equal(nb(s$val, ref, surv_outcome, horizon = 5, thresholds = c(0.1, 0.3))$estimates,
               res$estimates)
})

test_that("predict_risk() is internal", {
  ns <- readLines(system.file("NAMESPACE", package = "voipred.dev"))
  expect_false(any(grepl("^export\\(predict_risk\\)", ns)))
  expect_true(any(grepl("^S3method\\(predict_risk,coxph\\)", ns)))   # its methods are registered
})
