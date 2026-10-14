# The predicted risks given directly with `predicted_risk` (no model object),
# and a model's risks predicted once.

test_that("binary: risks given directly give the same results as the model", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  val <- g$val
  z <- c(0.02, 0.05, 0.1, 0.2)
  p <- unname(predict(f, val, type = "response"))
  m <- nb(f, val, "day30", thresholds = z)                         # model given first
  expect_equal(nb(model = f, data = val, outcome = "day30", thresholds = z)$NB_model,
               m$NB_model)                                         # or as `model`
  r <- nb(predicted_risk = p, data = val, outcome = "day30", thresholds = z)
  expect_equal(as.data.frame(r)[, -1], as.data.frame(m)[, -1])
  expect_equal(attr(r, "risk"), attr(m, "risk"))
  expect_equal(attr(r, "label"), "p")                            # the vector's name
  expect_equal(attr(m, "label"), "f")                            # the model's name
  val$pred <- p
  expect_equal(attr(nb("pred", val, "day30", thresholds = z), "label"), "pred")
  expect_equal(attr(nb(val$pred, val, "day30", thresholds = z), "label"), "pred")
  expect_equal(attr(nb(p, val, "day30", thresholds = z, label = "Full"), "label"), "Full")
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  expect_s3_class(plot_dca(r, Age = nb(f2, val, "day30", thresholds = z)), "ggplot")
})

test_that("survival and competing risks, with Kaplan-Meier and cIPCW", {
  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(surv_formula, data = s$dev)
  z <- c(0.1, 0.2, 0.3)
  p <- predict_risk(cx, v, horizon = 5)
  m <- quiet_nb(cx, v, surv_outcome, horizon = 5, thresholds = z)
  r <- quiet_nb(p, v, surv_outcome, horizon = 5, thresholds = z)
  expect_equal(r$NB_model, m$NB_model)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = v)
  mc <- nb(cx, v, surv_outcome, horizon = 5, thresholds = z, censoring = "cipcw",
           cipcw = cipcw_control(model = cens))
  rc <- nb(p, v, surv_outcome, horizon = 5, thresholds = z, censoring = "cipcw",
           cipcw = cipcw_control(model = cens))
  expect_equal(rc$NB_model, mc$NB_model)
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  pc <- predict_risk(csc, v, horizon = 5, cause = "death", cause_index = 2)
  mcr <- quiet_nb(csc, v, cr_outcome, horizon = 5, thresholds = z, cause = "death")
  rcr <- quiet_nb(pc, v, cr_outcome, horizon = 5, thresholds = z, cause = "death")
  expect_equal(rcr$NB_model, mcr$NB_model)
})

test_that("a model's risks are predicted once", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  calls <- 0
  pf <- function(model, newdata, ...) {
    calls <<- calls + 1
    predict(model, newdata, type = "response")
  }
  r <- nb(model = f, data = g$val, outcome = "day30", thresholds = c(0.05, 0.1),
          predict_function = pf)
  plot_dca(r)
  plot_dca(r, type = "delta")
  expect_equal(calls, 1)
})

test_that("`predicted_risk` is checked", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  val <- g$val
  p <- unname(predict(f, val, type = "response"))
  expect_error(nb(data = val, outcome = "day30"), "Give the predicted risks with `predicted_risk`")
  expect_error(nb(p, val, "day30", model = f), "either `predicted_risk` or `model`, not both")
  expect_error(nb(f, val, "day30", model = f), "either `predicted_risk` or `model`, not both")
  expect_error(nb(p, val, "day30", predict_function = function(...) p),
               "not used with `predicted_risk`")
  expect_error(nb(p[-1], val, "day30"), "values but `data` has")
  expect_error(nb(replace(p, 3, NA), val, "day30"), "1 of the risks")
  expect_error(nb(qlogis(p), val, "day30"), "outside \\[0, 1\\]")
  expect_error(nb(as.numeric(p > 0.1), val, "day30"), "only 0s and 1s")
  expect_error(nb("nope", val, "day30"), "not in `data`")
  expect_error(nb(as.character(p), val, "day30"), "must be numeric")
  expect_error(nb(model = p, data = val, outcome = "day30"), "give them with `predicted_risk =`")
})
