# Predicted risks given directly with `predicted_risk` (no model object), and the risks
# of a model predicted once and reused by every later step.

test_that("binary: risks given directly give the same results as the model", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  val <- g$val
  z <- c(0.02, 0.05, 0.1, 0.2)
  p <- unname(predict(f, val, type = "response"))
  m <- nb(val, f, "day30", thresholds = z)
  r <- nb(predicted_risk = p, data = val, outcome = "day30", thresholds = z)
  expect_equal(as.data.frame(r)[, -1], as.data.frame(m)[, -1])
  expect_equal(r$risk, m$risk)
  expect_equal(r$label, "p")                       # the vector's name
  val$pred <- p
  expect_equal(nb(predicted_risk = "pred", data = val, outcome = "day30", thresholds = z)$label, "pred")                        # a column of data
  expect_equal(nb(predicted_risk = val$pred, data = val, outcome = "day30", thresholds = z)$label, "pred")
  expect_equal(nb(predicted_risk = p, data = val, outcome = "day30", thresholds = z,
                       label = "Full")$label, "Full")
  # the same draws, so the same value of information
  set.seed(1); e_m <- evpi(m, n_sim = 200)
  set.seed(1); e_r <- evpi(r, n_sim = 200)
  expect_equal(e_r$estimates$EVPI, e_m$estimates$EVPI)
  set.seed(2); o_m <- evsi(m, future_n = c(200, 2000), n_sim = 50)
  set.seed(2); o_r <- evsi(r, future_n = c(200, 2000), n_sim = 50)
  expect_equal(o_r$estimates$EVSI, o_m$estimates$EVSI)
  # another prior, and models mixed with given risks
  expect_equal(nb(predicted_risk = p, data = val, outcome = "day30", thresholds = z,
                  prior = nb_uniform())$estimates$NB_model,
               nb(val, f, "day30", thresholds = z, prior = nb_uniform())$estimates$NB_model)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  set.seed(3)
  expect_equal(nrow(p_best(r, Age = nb(val, f2, "day30", thresholds = z), n_sim = 50)$estimates), 16)
  expect_s3_class(plot(r), "ggplot")
})

test_that("survival and competing risks, with Kaplan-Meier and cIPCW", {
  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(surv_formula, data = s$dev)
  z <- c(0.1, 0.2, 0.3)
  p <- predict_risk(cx, v, horizon = 5)
  m <- quiet_nb(v, cx, surv_outcome, horizon = 5, thresholds = z)
  r <- quiet_nb(predicted_risk = p, data = v, outcome = surv_outcome, horizon = 5, thresholds = z)
  expect_equal(r$estimates$NB_model, m$estimates$NB_model)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = v)
  mc <- nb(v, cx, surv_outcome, horizon = 5, thresholds = z, censoring = "cipcw",
           cipcw = cipcw_control(model = cens))
  rc <- nb(predicted_risk = p, data = v, outcome = surv_outcome, horizon = 5, thresholds = z,
           censoring = "cipcw", cipcw = cipcw_control(model = cens))
  expect_equal(rc$estimates$NB_model, mc$estimates$NB_model)
  # the censoring model is still refitted in each draw: same results
  set.seed(4); e_m <- evpi(mc, n_sim = 20)
  set.seed(4); e_r <- evpi(rc, n_sim = 20)
  expect_equal(e_r$estimates$EVPI, e_m$estimates$EVPI)
  # competing risks: the risk of the cause by the horizon
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = s$dev)
  pc <- predict_risk(csc, v, horizon = 5, cause = "death", cause_index = 2)
  mcr <- quiet_nb(v, csc, cr_outcome, horizon = 5, thresholds = z, cause = "death")
  rcr <- quiet_nb(predicted_risk = pc, data = v, outcome = cr_outcome, horizon = 5, thresholds = z,
                  cause = "death")
  expect_equal(rcr$estimates$NB_model, mcr$estimates$NB_model)
})

test_that("a model's risks are predicted once and reused by every later step", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  calls <- 0
  pf <- function(model, newdata, ...) {
    calls <<- calls + 1
    predict(model, newdata, type = "response")
  }
  r <- nb(g$val, f, "day30", thresholds = c(0.05, 0.1), predict_function = pf)
  expect_equal(calls, 1)
  set.seed(1)
  e <- evpi(r, n_sim = 50)
  p_useful(r, n_sim = 50)
  p_best(r, n_sim = 50)
  o <- evsi(r, future_n = c(100, 1000), n_sim = 20)
  plot(r)
  plot(e)
  rb <- nb(g$val, f, "day30", thresholds = c(0.05, 0.1), predict_function = pf,
           prior = nb_uniform())
  evpi(rb, n_sim = 50)
  evsi(rb, future_n = 100, n_sim = 50)
  expect_equal(calls, 2)                  # once per call of nb(), never again
})

test_that("`predicted_risk` is checked", {
  g <- gusto_split()
  f <- glm(gusto_formula, family = binomial, data = g$dev)
  val <- g$val
  p <- unname(predict(f, val, type = "response"))
  expect_error(nb(data = val, outcome = "day30"), "Give the predicted risks with `predicted_risk`")
  expect_error(nb(val, p, "day30", model = f), "either `predicted_risk` or `model`, not both")
  expect_error(nb(predicted_risk = p, data = val, outcome = "day30",
                  predict_function = function(...) p), "not used with `predicted_risk`")
  expect_error(nb(predicted_risk = p[-1], data = val, outcome = "day30"), "values but `data` has")
  expect_error(nb(predicted_risk = replace(p, 3, NA), data = val, outcome = "day30"), "1 of the risks")
  expect_error(nb(predicted_risk = qlogis(p), data = val, outcome = "day30"), "outside \\[0, 1\\]")
  expect_error(nb(predicted_risk = as.numeric(p > 0.1), data = val, outcome = "day30"), "only 0s and 1s")
  expect_error(nb(predicted_risk = "nope", data = val, outcome = "day30"), "not in `data`")
  expect_error(nb(predicted_risk = as.character(p), data = val, outcome = "day30"), "must be numeric")
  expect_error(nb(model = p, data = val, outcome = "day30"), "give them with `predicted_risk =`")
})
