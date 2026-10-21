# Several models in one call of nb_val(), combining results with c(), and
# P(useful) and P(best) computed by evpi() from its draws.

multi_case <- function() {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  val <- g$val
  val$p1 <- unname(predict(f1, val, type = "response"))
  val$p2 <- unname(predict(f2, val, type = "response"))
  list(f1 = f1, f2 = f2, val = val, z = c(0.02, 0.05, 0.1, 0.2))
}

test_that("every way of giving several models gives the same result", {
  m <- multi_case()
  v <- m$val
  z <- m$z
  ref <- nb_val(v, list(Full = v$p1, Age = v$p2), "day30", thresholds = z)
  expect_s3_class(ref, "voipred_nb")
  expect_equal(ref$models, c("Full", "Age"))
  expect_equal(dim(ref$risk), c(nrow(v), 2))
  expect_equal(colnames(ref$risk), c("Full", "Age"))
  expect_named(ref$estimates, c("threshold", "model", "NB_model", "NB_all", "NB_none",
                                "delta_NB"))
  expect_equal(ref$estimates$model, rep(c("Full", "Age"), each = length(z)))
  same <- function(r) {
    expect_equal(r$estimates, ref$estimates)
    expect_equal(r$risk, ref$risk)
    expect_equal(r$counts, ref$counts)
  }
  same(nb_val(v, c(Full = "p1", Age = "p2"), "day30", thresholds = z))
  same(nb_val(v, data.frame(Full = v$p1, Age = v$p2), "day30", thresholds = z))
  same(nb_val(v, cbind(Full = v$p1, Age = v$p2), "day30", thresholds = z))
  same(nb_val(v, model = list(Full = m$f1, Age = m$f2), outcome = "day30", thresholds = z))
  same(nb_val(v, list(Full = m$f1, Age = v$p2), "day30", thresholds = z))
  same(nb_val(v, list(Full = "p1", Age = "p2"), "day30", thresholds = z))
  same(nb_val(v, c("p1", "p2"), "day30", thresholds = z, label = c("Full", "Age")))
  # each model as computed alone
  a <- nb_val(v, v$p1, "day30", thresholds = z)
  b <- nb_val(v, v$p2, "day30", thresholds = z)
  expect_equal(ref$estimates$NB_model[1:4], a$estimates$NB_model)
  expect_equal(ref$estimates$NB_model[5:8], b$estimates$NB_model)
  expect_equal(ref$estimates$NB_all[5:8], a$estimates$NB_all)
})

test_that("labels: names, expressions, columns, defaults and `label`", {
  m <- multi_case()
  v <- m$val
  z <- 0.1
  expect_equal(nb_val(v, list(v$p1, m$f2), "day30", thresholds = z)$models, c("p1", "f2"))
  expect_equal(nb_val(v, c("p1", "p2"), "day30", thresholds = z)$models, c("p1", "p2"))
  expect_equal(nb_val(v, v[c("p1", "p2")], "day30", thresholds = z)$models, c("p1", "p2"))
  expect_equal(nb_val(v, cbind(v$p1, v$p2), "day30", thresholds = z)$models,
               c("Model 1", "Model 2"))
  expect_equal(nb_val(v, list(v$p1 + 0, A = v$p2), "day30", thresholds = z)$models,
               c("Model 1", "A"))
  # one column of a data frame is one model, labelled by the column name
  one <- nb_val(v, v["p1"], "day30", thresholds = z)
  expect_equal(one$models, "p1")
  expect_equal(nb_val(v, c(Full = "p1"), "day30", thresholds = z)$models, "Full")
  expect_error(nb_val(v, list(A = v$p1, A = v$p2), "day30", thresholds = z),
               "same label")
  expect_error(nb_val(v, list(`Treat all` = v$p1, B = v$p2), "day30", thresholds = z),
               "Treat None")
  expect_error(nb_val(v, list(A = v$p1, B = v$p2), "day30", thresholds = z, label = "x"),
               "2 non-empty strings")
  expect_error(nb_val(v, v$p1, "day30", thresholds = z, label = c("a", "b")),
               "single non-empty string")
})

test_that("several models: errors name the model", {
  m <- multi_case()
  v <- m$val
  z <- 0.1
  expect_error(nb_val(v, list(A = v$p1, B = qlogis(v$p2)), "day30", thresholds = z),
               "`predicted_risk` \\(\"B\"\\) has values outside \\[0, 1\\]")
  expect_error(nb_val(v, list(A = v$p1, B = v$p2[-1]), "day30", thresholds = z),
               "\\(\"B\"\\) has .* values but `data` has")
  expect_error(nb_val(v, c("p1", "nope"), "day30", thresholds = z),
               "names a column \\(\"nope\"\\) that is not in `data`")
  expect_error(nb_val(v, model = list(A = m$f1, B = v$p2), outcome = "day30", thresholds = z),
               "\"B\": `model` must be a fitted model")
  expect_error(nb_val(v, list(A = v$p1, B = list(v$p2)), "day30", thresholds = z),
               "Each element of `predicted_risk`")
  expect_error(nb_val(v, list(), "day30", thresholds = z), "empty list")
  expect_error(nb_val(v, data.frame(), "day30", thresholds = z), "no columns")
  # predict_function: one for all the models, or one per model by label
  calls <- character()
  pf <- function(model, newdata, ...) {
    calls <<- c(calls, "pf")
    predict(model, newdata, type = "response")
  }
  r <- nb_val(v, model = list(A = m$f1, B = m$f2), outcome = "day30", thresholds = z,
              predict_function = pf)
  expect_equal(length(calls), 2)
  r2 <- nb_val(v, model = list(A = m$f1, B = m$f2), outcome = "day30", thresholds = z,
               predict_function = list(B = pf))
  expect_equal(length(calls), 3)
  expect_equal(r2$risk, r$risk)
  expect_error(nb_val(v, model = list(A = m$f1, B = m$f2), outcome = "day30",
                      thresholds = z, predict_function = list(C = pf)),
               "not one of the models")
  expect_error(nb_val(v, list(A = v$p1, B = v$p2), "day30", thresholds = z,
                      predict_function = pf), "not used with `predicted_risk`")
})

test_that("c() combines results, and plot() and the VoI functions combine theirs", {
  m <- multi_case()
  v <- m$val
  z <- m$z
  ref <- nb_val(v, list(Full = v$p1, Age = v$p2), "day30", thresholds = z)
  a <- nb_val(v, v$p1, "day30", thresholds = z, label = "Full")
  b <- nb_val(v, v$p2, "day30", thresholds = z)
  ab <- c(a, Age = b)
  expect_equal(ab$estimates, ref$estimates)
  expect_equal(ab$risk, ref$risk)
  expect_equal(ab$counts, ref$counts)
  expect_equal(c(ab)$models, ab$models)
  expect_equal(c(ref, Third = b)$models, c("Full", "Age", "Third"))
  expect_error(c(Both = ref, b), "relabels a result with one model")
  expect_error(c(a, a), "same label")
  expect_error(c(a, 1), "results of nb_val\\(\\) only")
  other <- nb_val(v, v$p2, "day30", thresholds = c(0.02, 0.05))
  expect_error(c(a, other), "different thresholds")
  # the same draws and results whichever way the models are given
  set.seed(1); e1 <- evpi(ref, n_sim = 200)
  set.seed(1); e2 <- evpi(a, Age = b, n_sim = 200)
  expect_equal(e1$estimates, e2$estimates)
  expect_equal(e1$models, c("Full", "Age"))
  set.seed(1); e3 <- evpi(ref, n_sim = 200, labels = c("x", "y"))
  expect_equal(colnames(e3$ENB)[2:3], c("x", "y"))
  expect_equal(levels(plot(ref)$data$label), c("Treat All", "Treat None", "Full", "Age"))
  expect_equal(levels(plot(a, Age = b)$data$label), c("Treat All", "Treat None", "Full", "Age"))
  expect_equal(levels(plot(ref, type = "delta", labels = c("x", "y"))$data$label),
               c("x", "y"))
  expect_error(plot(ref, labels = "x"), "2 strings, one per model")
})

test_that("evpi() computes P(useful) and P(best) from the same draws", {
  m <- multi_case()
  v <- m$val
  res <- nb_val(v, list(Full = v$p1, Age = v$p2), "day30", thresholds = m$z)
  set.seed(3); e <- evpi(res, n_sim = 300)
  set.seed(3); pu <- p_useful(res, n_sim = 300)
  set.seed(3); pb <- p_best(res, n_sim = 300)
  expect_s3_class(p_useful(e), "voipred_p_useful")
  expect_s3_class(p_best(e), "voipred_p_best")
  expect_equal(p_useful(e)$estimates, pu$estimates)
  expect_equal(p_best(e)$estimates, pb$estimates)
  expect_identical(p_useful(e), e$p_useful)
  expect_equal(unique(pu$estimates$model), c("Full", "Age"))
  expect_error(p_useful(e, n_sim = 10), "takes no other arguments")
  expect_error(p_best(e, res), "takes no other arguments")
  # summary() adds P(best) of the strategy that is best now
  s <- summary(e)
  expect_true("P(best)" %in% names(s$table))
  pbw <- wide(pb$estimates, "strategy", "p_best", pb$strategies)
  expect_equal(s$table[["P(best)"]],
               as.matrix(pbw[, -1])[cbind(seq_len(nrow(pbw)),
                                          match(s$table[["best now"]], pb$strategies))])
  expect_output(print(e), "From the same draws")
  expect_s3_class(plot(e, type = "p_useful"), "ggplot")
  expect_equal(levels(plot(e, type = "p_best")$data$label),
               c("Treat All", "Treat None", "Full", "Age"))
  expect_error(plot(e, e, type = "p_best"), "one result of evpi")
  # P(useful) of several models: plot and summary
  expect_equal(levels(plot(pu)$data$label), c("Full", "Age"))
  expect_equal(levels(plot(pu, labels = c("a", "b"))$data$label), c("a", "b"))
  expect_error(plot(pu, labels = "a"), "one per curve")
  at <- summary(pu, thresholds = c(0.05, 0.1))$at
  expect_named(at, c("threshold", "Full", "Age"))
})

test_that("several models with survival outcomes, Kaplan-Meier and cIPCW", {
  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(surv_formula, data = s$dev)
  cx2 <- survival::coxph(survival::Surv(dtimey, death) ~ age, data = s$dev)
  z <- c(0.1, 0.2, 0.3)
  v$p1 <- predict_risk(cx, v, horizon = 5)
  v$p2 <- predict_risk(cx2, v, horizon = 5)
  km <- quiet_nb(v, c(Full = "p1", Age = "p2"), surv_outcome, horizon = 5, thresholds = z)
  a <- quiet_nb(v, "p1", surv_outcome, horizon = 5, thresholds = z)
  b <- quiet_nb(v, "p2", surv_outcome, horizon = 5, thresholds = z)
  expect_equal(km$estimates$NB_model, c(a$estimates$NB_model, b$estimates$NB_model))
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = v)
  ci <- nb_val(v, c(Full = "p1", Age = "p2"), surv_outcome, horizon = 5, thresholds = z,
               censoring = "cipcw", cipcw = cipcw_control(model = cens))
  ca <- nb_val(v, "p1", surv_outcome, horizon = 5, thresholds = z,
               censoring = "cipcw", cipcw = cipcw_control(model = cens), label = "Full")
  cb <- nb_val(v, "p2", surv_outcome, horizon = 5, thresholds = z,
               censoring = "cipcw", cipcw = cipcw_control(model = cens), label = "Age")
  expect_equal(ci$estimates, c(ca, cb)$estimates)
  # the summary has the sensitivity and specificity of each model
  sm <- summary(ci)
  expect_equal(sm$table$Se[sm$table$model == "Full"], sens_spec_cipcw(ca)$Se)
  smk <- summary(km)
  expect_equal(smk$table$Sp[smk$table$model == "Age"], sens_spec_km(b)$Sp)
  # the censoring model is refitted once per draw for both models: the same
  # draws as the results computed separately
  set.seed(5); e1 <- evpi(ci, n_sim = 10)
  set.seed(5); e2 <- evpi(ca, cb, n_sim = 10)
  expect_equal(e1$estimates, e2$estimates)
  expect_equal(p_best(e1)$estimates, p_best(e2)$estimates)
  # results with different censoring models cannot be combined
  cens2 <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age, data = v)
  cc <- nb_val(v, "p2", surv_outcome, horizon = 5, thresholds = z,
               censoring = "cipcw", cipcw = cipcw_control(model = cens2), label = "Age")
  expect_error(c(ca, cc), "cannot be used together")
})

test_that("the joint posterior of many models stops before using too much memory", {
  m <- multi_case()
  v <- m$val
  risks <- stats::setNames(lapply(1:12, function(i) pmin(v$p1 * (0.9 + i / 50), 1)),
                           paste0("m", 1:12))
  res <- nb_val(v, risks, "day30", thresholds = seq(0.01, 0.5, by = 0.01))
  expect_error(evpi(res, n_sim = 1000), "needs about .* GB of memory")
  expect_silent(check_joint_size(list(M = 3, K = 8, z = 1:50), 1000))
})
