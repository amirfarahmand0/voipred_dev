# Shared data and helpers for the survival tests.
#
# Rotterdam breast cancer data (survival package): death within 5 years
# (standard survival), and recurrence vs death without recurrence (competing
# risks). Developed on 1,500 patients, validated on the other 1,482.

rott_split <- function(seed = 1) {
  d <- survival::rotterdam
  d$dtimey <- d$dtime / 365.25
  d$cr_time <- ifelse(d$recur == 1, d$rtime, d$dtime) / 365.25
  d$cr_code <- ifelse(d$recur == 1, 1L, ifelse(d$death == 1, 2L, 0L))
  d$cr_status <- factor(d$cr_code, 0:2, c("censored", "recurrence", "death"))
  set.seed(seed)
  i <- sample(nrow(d))
  list(dev = d[i[1:1500], ], val = d[i[1501:2982], ])
}

strata <- survival::strata      # for stratified models fitted in the tests

surv_formula <- survival::Surv(dtimey, death) ~ age + meno + size + grade + nodes +
  pgr + er + hormon + chemo
# randomForestSRC needs a plain Surv() on the left-hand side.
rsf_formula <- Surv(dtimey, death) ~ age + meno + size + grade + nodes + pgr + er +
  hormon + chemo
surv_outcome <- ~ survival::Surv(dtimey, death)
cr_outcome <- ~ survival::Surv(cr_time, cr_status)
surv_z <- seq(0, 0.99, by = 0.01)

# nb_val() for a survival outcome vs dcurves::dca() (KM / Aalen-Johansen) on the
# same predicted risks. `levels` reorders a competing-risks status so that the
# cause of interest comes first, as dcurves always uses the first event level.
expect_same_as_dcurves_surv <- function(res, data, time, status, horizon,
                                        levels = NULL) {
  d <- data.frame(tt = data[[time]], ss = data[[status]], risk = res$risk[, 1])
  if (!is.null(levels)) d$ss <- factor(d$ss, levels = levels)
  ref <- suppressMessages(suppressWarnings(
    dcurves::dca(survival::Surv(tt, ss) ~ risk, data = d, time = horizon,
                 thresholds = res$estimates$threshold)))$dca
  m <- ref$net_benefit[ref$variable == "risk"]
  a <- ref$net_benefit[ref$variable == "all"]
  # dcurves quirk: NA when no treated patient has an event of interest by the
  # horizon (the risk among the treated is 0). nb_val() then gives -q z / (1 - z).
  ev <- as.integer(d$ss) == if (is.factor(d$ss)) 2L else 1L
  quirk <- is.na(m) & !is.na(res$estimates$NB_model)
  for (j in which(quirk)) {
    z <- res$estimates$threshold[j]; tr <- d$risk >= z & d$risk > 0
    expect_false(any(tr & ev & d$tt <= horizon))
    expect_equal(res$estimates$NB_model[j], -mean(tr) * z / (1 - z), tolerance = 1e-12)
  }
  expect_false(any(is.na(res$estimates$NB_model) & !is.na(m)))
  expect_equal(res$estimates$NB_model[!is.na(m)], m[!is.na(m)], tolerance = 1e-10)
  expect_equal(res$estimates$NB_all, a, tolerance = 1e-10)
}

quiet_nb <- function(...) suppressWarnings(nb_val(...))

# Simulated survival data with known truth. Event times are Weibull with a
# hazard depending on x1 and x2; censoring is exponential, either independent
# of the covariates or depending strongly on x2 (informative unless x2 is
# adjusted for). `true_nb()` uses the uncensored event times.
sim_surv <- function(n, informative = TRUE, seed = 1, competing = FALSE) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  lp1 <- 0.8 * x1 + 0.8 * x2
  t1 <- (-log(runif(n)) / (0.15 * exp(lp1)))^(1 / 1.3)
  if (competing) {
    t2 <- (-log(runif(n)) / (0.10 * exp(0.6 * x2 - 0.5 * x1)))^(1 / 1.1)
    t_ev <- pmin(t1, t2); cause <- ifelse(t1 <= t2, 1L, 2L)
  } else {
    t_ev <- t1; cause <- rep(1L, n)
  }
  rate_c <- if (informative) 0.12 * exp(1.2 * x2) else rep(0.12, n)
  t_c <- rexp(n, rate_c)
  time <- pmin(t_ev, t_c)
  status <- ifelse(t_ev <= t_c, cause, 0L)
  d <- data.frame(x1, x2, time, status, t_ev, cause, lp1)
  if (competing) d$status_f <- factor(status, 0:2, c("censored", "cause1", "cause2"))
  d
}

# True net benefit at horizon tau for risks p, from the uncensored data.
true_nb <- function(d, p, z, tau) {
  y <- as.integer(d$t_ev <= tau & d$cause == 1)
  sapply(z, function(t) { tr <- p >= t & p > 0
    mean(y * tr) - mean((1 - y) * tr) * t / (1 - t) })
}

# A test-only censoring "model": riskRegression::predictRisk returns
# 1 - G(t) from a fixed function, so censoring probabilities can be set exactly.
fixed_censoring <- function(g_fun, vars = "x2") structure(list(g = g_fun, vars = vars),
                                                            class = "fixed_censoring")
predictRisk.fixed_censoring <- function(object, newdata, times, ...) {
  1 - outer(seq_len(nrow(newdata)), times, function(i, t) object$g(newdata[i, object$vars], t))
}
if (requireNamespace("riskRegression", quietly = TRUE))
  registerS3method("predictRisk", "fixed_censoring", predictRisk.fixed_censoring,
                   envir = asNamespace("riskRegression"))

# Evaluate `expr`, muffling only the warnings whose message matches `pattern`.
suppress_matching <- function(expr, pattern) withCallingHandlers(expr, warning = function(w)
  if (grepl(pattern, conditionMessage(w))) invokeRestart("muffleWarning"))

# randomForestSRC and aorsf parse formulas with as.character(); formula.tools
# (loaded by logistf in the binary-model tests) overrides that method for
# formulas and breaks them. Use the base behaviour for the rest of the calling
# test file, then restore.
local_base_formula <- function(env = parent.frame()) {
  if (!"formula.tools" %in% loadedNamespaces()) return(invisible())
  registerS3method("as.character", "formula", function(x, ...) as.character(unclass(x)))
  withr::defer(registerS3method("as.character", "formula",
                                getFromNamespace("as.character.formula", "formula.tools")),
               envir = env)
}

# rstpm2 only fits and predicts when attached; attaching it also attaches
# survival and splines. Detach whatever it added afterwards.
with_rstpm2 <- function(code) {
  before <- search()
  suppressPackageStartupMessages(attachNamespace("rstpm2"))
  on.exit(for (p in setdiff(search(), before))
    detach(p, character.only = TRUE, unload = FALSE), add = TRUE)
  force(code)
}
