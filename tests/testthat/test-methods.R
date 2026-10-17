# The S3 classes of the results: their components and their print(),
# summary() and as.data.frame() methods (plot() is in test-plot*.R).

methods_setup <- function() {
  g <- gusto_split()
  f1 <- glm(gusto_formula, family = binomial, data = g$dev)
  f2 <- glm(day30 ~ age, family = binomial, data = g$dev)
  val <- g$val[1:500, ]
  z <- seq(0.01, 0.3, by = 0.01)
  r1 <- nb(val, f1, "day30", thresholds = z, label = "Full")
  r2 <- nb(val, f2, "day30", thresholds = z, label = "Age")
  set.seed(1)
  list(val = val, z = z, r1 = r1, r2 = r2,
       e = evpi(r1, r2, n_sim = 200), pu = p_useful(r1, n_sim = 200),
       pb = p_best(r1, r2, n_sim = 200), o = evsi(r1, r2, future_n = c(250, 1000), n_sim = 200))
}

test_that("nb() returns a list with its parts, not attributes", {
  s <- methods_setup()
  r <- s$r1
  expect_s3_class(r, "voipred_nb")
  expect_type(r, "list")
  expect_named(r, c("call", "label", "estimates", "prior", "outcome", "risk", "observed",
                    "weights", "counts", "cipcw"))
  expect_null(attributes(r$estimates)$voi)
  expect_equal(r$outcome$type, "binary")
  expect_equal(r$outcome$n, 500)
  expect_equal(r$outcome$events, sum(s$val$day30))
  expect_equal(r$outcome$prevalence, mean(s$val$day30))
  expect_length(r$risk, 500)
  expect_equal(r$counts$threshold, s$z)
  expect_equal(as.data.frame(r), r$estimates)
  # survival: horizon, cause and censoring; no counts
  sv <- rott_split()
  cx <- survival::coxph(surv_formula, data = sv$dev)
  k <- quiet_nb(sv$val, cx, surv_outcome, horizon = 5, thresholds = c(0.1, 0.3))
  expect_equal(k$outcome[c("type", "horizon", "censoring")],
               list(type = "survival", horizon = 5, censoring = "km"))
  expect_null(k$counts)
  expect_null(k$cipcw)
})

test_that("print() and summary() of nb() are short and say what the result is", {
  s <- methods_setup()
  out <- capture.output(print(s$r1))
  expect_match(out[1], "^Net benefit of Full: binary, n = 500, [0-9]+ events \\([0-9.]+%\\); prior improper Beta\\(0, 0\\)$")
  expect_length(out, 1 + 1 + 1 + 10 + 1)                     # header, blank, column names, 10 rows, footer
  expect_match(out[length(out)], "^# 20 more rows; use summary\\(\\) or as.data.frame\\(\\)$")
  expect_length(capture.output(print(s$r1, n = 3)), 1 + 1 + 1 + 3 + 1)
  sm <- summary(s$r1)
  expect_s3_class(sm, "summary.voipred_nb")
  so <- capture.output(sm)
  expect_match(so[2], "^Best strategy: ")
  expect_match(so[3], "^Largest delta_NB: [0-9.]+ at 0\\.[0-9]+$")
  expect_named(sm$table, c("threshold", "Se", "Sp", "NB_model", "NB_all", "delta_NB", "best"))
  # sensitivity and specificity are the observed proportions with the improper prior
  y <- s$val$day30
  tr <- s$r1$risk >= 0.1 & s$r1$risk > 0
  i <- which(abs(s$z - 0.1) < 1e-9)
  expect_equal(sm$table$Se[i], sum(tr & y == 1) / sum(y == 1))
  expect_equal(sm$table$Sp[i], sum(!tr & y == 0) / sum(y == 0))
  # the best strategy follows delta_NB and NB_all
  e <- s$r1$estimates
  expect_equal(sm$table$best, ifelse(e$delta_NB > 1e-12, "model",
                                     ifelse(e$NB_all > 1e-12, "treat all", "treat none")))
  # chosen thresholds: all of them shown
  ch <- capture.output(summary(s$r1, thresholds = c(0.05, 0.2)))
  expect_false(any(grepl("more rows", ch)))
  expect_equal(nrow(summary(s$r1, thresholds = c(0.05, 0.2))$table), 2)
  expect_error(summary(s$r1, thresholds = 0.055), "not among the thresholds")
  # survival: the share treated instead of sensitivity and specificity
  sv <- rott_split()
  cx <- survival::coxph(surv_formula, data = sv$dev)
  cens <- survival::coxph(survival::Surv(dtimey, death == 0) ~ age + nodes, data = sv$val)
  k <- nb(sv$val, cx, surv_outcome, horizon = 5, thresholds = c(0.1, 0.3),
          censoring = "cipcw", cipcw = cipcw_control(model = cens), label = "Cox")
  expect_match(capture.output(print(k))[1],
               "^Net benefit of Cox: survival, horizon 5, n = 1482, risk by horizon [0-9.]+%; censoring cIPCW \\(coxph\\); prior improper")
  st <- summary(k)$table
  expect_named(st, c("threshold", "Se", "Sp", "NB_model", "NB_all", "delta_NB", "best"))
  # Se and Sp are those implied by the net benefit and the risk by the horizon
  F <- k$outcome$prevalence
  odds <- c(0.1, 0.3) / c(0.9, 0.7)
  expect_equal(F * st$Se - (1 - F) * (1 - st$Sp) * odds, st$NB_model, tolerance = 1e-12)
  # with cIPCW they are the estimator's own IPCW-weighted proportions
  v <- voi_info(k)
  we <- ifelse(v$status == 1 & v$time <= 5, v$weight_event, 0)
  wh <- ifelse(v$time > 5, v$weight_horizon, 0)
  tr <- outer(k$risk, c(0.1, 0.3), ">=")
  expect_equal(st$Se, unname(colSums(tr * we) / sum(we)), tolerance = 1e-12)
  expect_equal(st$Sp, unname(1 - colSums(tr * wh) / sum(wh)), tolerance = 1e-12)
  csc <- riskRegression::CSC(prodlim::Hist(cr_time, cr_code) ~ age + nodes, data = sv$dev)
  cr <- quiet_nb(sv$val, csc, cr_outcome, horizon = 5, thresholds = 0.2, cause = "death")
  expect_match(capture.output(print(cr))[1], "competing risks \\(cause death\\), horizon 5")
  # case weights
  w <- rep(1:2, 250)
  rw <- nb(s$val, s$val$day30 * 0.5 + 0.1, "day30", thresholds = 0.2, weights = w, label = "W")
  expect_match(capture.output(print(rw))[1], "n = 500 \\(case weights\\), weighted prevalence")
})

test_that("the value-of-information results are lists with print(), summary() and as.data.frame()", {
  s <- methods_setup()
  for (x in list(s$e, s$pu, s$pb, s$o)) {
    models <- if (inherits(x, "voipred_p_useful")) "Full" else c("Full", "Age")
    expect_type(x, "list")
    expect_true(all(c("call", "estimates", "models", "outcome", "prior", "draws", "n_sim") %in%
                      names(x)))
    expect_equal(x$models, models)
    expect_equal(x$draws, "posterior")
    expect_equal(x$n_sim, 200)
    expect_identical(as.data.frame(x), x$estimates)
    out <- capture.output(print(x))
    expect_match(out[1], paste0(paste(models, collapse = ", "), "; binary, n = 500, prior ",
                                "improper Beta\\(0, 0\\); 200 posterior draws$"))
    expect_lte(length(out), 15)
    expect_output(print(summary(x)), paste(models, collapse = ", "))
  }
  expect_s3_class(s$e, "voipred_evpi")
  expect_s3_class(s$o, "voipred_evsi")
  expect_s3_class(s$pu, "voipred_p_useful")
  expect_s3_class(s$pb, "voipred_p_best")
  expect_named(s$e$ENB, c("threshold", "Full", "Age", "Treat All", "Treat None"))
  expect_equal(s$o$future_n, c(250, 1000))
  expect_equal(s$pb$strategies, c("Full", "Age", "Treat All", "Treat None"))
  # print(): wide tables
  expect_match(capture.output(print(s$pu))[3], "threshold +Full$")
  expect_match(capture.output(print(s$pb))[3], "threshold +Full +Age +Treat All +Treat None")
  expect_match(capture.output(print(s$o))[3], "threshold +EVPI +N=250 +N=1000")
  # summary(): the main findings
  se <- summary(s$e)
  expect_s3_class(se, "summary.voipred_evpi")
  expect_match(se$lines, "^EVPI > 0 at .*; largest [0-9.]+ at 0\\.[0-9]+$")
  enb <- as.matrix(s$e$ENB[, -1])
  expect_equal(se$table[["best now"]], colnames(enb)[apply(enb, 1, which.max)])
  so <- summary(s$o)
  expect_equal(nrow(so$table), 3)                       # the 3 thresholds with the largest EVPI
  top <- s$o$estimates[s$o$estimates$future_n == 250, ]
  expect_setequal(so$table$threshold, top$threshold[order(-top$EVPI)][1:3])
  expect_match(so$table[["N=1000"]][1], "^[0-9-]+%$|^-$")
  expect_match(so$lines, "with the largest EVPI")
  expect_equal(nrow(summary(s$o, thresholds = 0.05)$table), 1)
  expect_match(summary(s$o, thresholds = 0.05)$lines, "at the chosen thresholds")
  sp <- summary(s$pu)
  expect_named(sp$table, c("model", "P > 0.95", "0.5 to 0.95", "P < 0.5"))
  expect_equal(sp$table$model, "Full")
  sb <- summary(s$pb)
  w <- as.matrix(wide(s$pb$estimates, "strategy", "p_best", s$pb$strategies)[, -1])
  expect_equal(sb$table[["P(best)"]], apply(w, 1, max))
  # survival: Bayesian bootstrap draws
  sv <- rott_split()
  cx <- survival::coxph(surv_formula, data = sv$dev)
  k <- quiet_nb(sv$val, cx, surv_outcome, horizon = 5, thresholds = c(0.2, 0.3), label = "Cox")
  set.seed(1)
  ek <- evpi(k, n_sim = 20)
  expect_equal(ek$draws, "Bayesian bootstrap")
  expect_match(capture.output(print(ek))[1],
               "^EVPI: Cox; survival, horizon 5, n = 1482; 20 Bayesian bootstrap draws$")
})

test_that("ranges of thresholds are written compactly", {
  z <- seq(0.01, 0.1, by = 0.01)
  expect_equal(segments_of(z, z <= 0.03 | z >= 0.09), "0.01-0.03, 0.09-0.10")
  expect_equal(segments_of(z, z == 0.05), "0.05")
  expect_equal(segments_of(z, rep(FALSE, 10)), "none")
  expect_equal(segments_text(z, rep(c("model", "treat none"), c(6, 4))),
               "model 0.01-0.06, treat none 0.07-0.10")
})

test_that("summary() has the same columns for every outcome; without censoring Se and Sp equal the binary ones", {
  # KM (dcurves) Se and Sp: from the Kaplan-Meier risk overall and among the treated
  g <- gusto_split()
  v <- g$val[1:2000, ]
  p <- unname(predict(glm(gusto_formula, family = binomial, data = g$dev), v, type = "response"))
  z <- c(0.02, 0.05, 0.1, 0.2)
  b <- summary(nb(data = v, predicted_risk = p, outcome = "day30", thresholds = z))$table
  # the same outcome as survival data with everyone followed to the horizon
  v$t <- ifelse(v$day30 == 1, 10, 40)
  k <- summary(nb(data = v, predicted_risk = p, outcome = ~ survival::Surv(t, day30),
                  horizon = 30, thresholds = z))$table
  expect_named(k, names(b))
  expect_equal(k$Se, b$Se, tolerance = 1e-12)
  expect_equal(k$Sp, b$Sp, tolerance = 1e-12)
  expect_equal(k$NB_model, b$NB_model, tolerance = 1e-12)
})
