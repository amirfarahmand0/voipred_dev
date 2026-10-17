# print(), summary() and as.data.frame() methods of the results of nb(),
# p_useful(), p_best(), evpi() and evsi(). (The plot() methods are in plot.R
# and plot_voi.R.)

# ---- nb() --------------------------------------------------------------------

#' Summarise the net benefit of a model
#'
#' `print()` shows the net benefit at the first thresholds; `summary()` adds
#' the best strategy at each threshold, where the model is useful, and the
#' model's sensitivity and specificity;
#' `as.data.frame()` returns the estimates at every threshold.
#'
#' @param object,x A result of [nb()].
#' @param thresholds Optional thresholds to show in the table (among those of
#'   the result); default: all, of which the first 10 are printed.
#' @param n Number of thresholds printed. Default 10.
#' @param row.names,optional Not used (as in [as.data.frame()]).
#' @param ... Not used.
#' @return `summary()`: an object of class `summary.voipred_nb`, printed as
#'   shown. `as.data.frame()`: the `estimates` data frame (columns `threshold`,
#'   `NB_model`, `NB_all`, `NB_none`, `delta_NB`).
#' @details The best strategy at a threshold is the model where
#'   `delta_NB > 0`, otherwise treat all if its net benefit is positive,
#'   otherwise treat none. For a binary outcome, sensitivity and
#'   specificity are posterior means (the observed proportions with the default
#'   improper prior). For survival and competing-risks outcomes they are those
#'   of the event by the horizon, from the same estimates as the net benefit:
#'   with `censoring = "km"`, from the Kaplan-Meier (Aalen-Johansen) risk by
#'   the horizon \eqn{F}, the share treated \eqn{q} and the risk among the
#'   treated \eqn{F_+}, \eqn{Se = q F_+ / F} and
#'   \eqn{Sp = 1 - q (1 - F_+) / (1 - F)}; with `censoring = "cipcw"`, the
#'   IPCW-weighted proportions of the cIPCW estimator.
#' @seealso [nb()], [plot.voipred_nb()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb(data = df_val_binary,
#'           predicted_risk = df_val_binary$p,
#'           outcome = "day30",
#'           thresholds = seq(0.01, 0.3, by = 0.01))
#' res
#' summary(res)
#' summary(res, thresholds = c(0.05, 0.1, 0.2))
#' head(as.data.frame(res))
#' @name summary.voipred_nb
NULL

#' @rdname summary.voipred_nb
#' @export
print.voipred_nb <- function(x, n = 10, ...) {
  cat(nb_header(x), "\n\n", sep = "")
  print_rows(x$estimates, n, "use summary() or as.data.frame()", digits = 4)
  invisible(x)
}

#' @rdname summary.voipred_nb
#' @export
summary.voipred_nb <- function(object, thresholds = NULL, ...) {
  e <- object$estimates
  best <- best_strategy(e)
  o <- order(e$threshold)
  tab <- data.frame(threshold = e$threshold)
  ss <- if (object$outcome$type == "binary") sens_spec(object)
        else if (object$outcome$censoring == "cipcw") sens_spec_cipcw(object)
        else sens_spec_km(object)
  tab$Se <- ss$Se
  tab$Sp <- ss$Sp
  tab$NB_model <- e$NB_model
  tab$NB_all <- e$NB_all
  tab$delta_NB <- e$delta_NB
  tab$best <- best
  rows <- select_rows(e$threshold, thresholds)
  ok <- is.finite(e$delta_NB)
  largest <- if (any(ok & e$delta_NB > 1e-12)) {
    j <- which(ok)[which.max(e$delta_NB[ok])]
    sprintf("Largest delta_NB: %s at %s", fmt_num(e$delta_NB[j]), fmt_zs(e$threshold)[j])
  } else "The model is not better than both default strategies at any threshold"
  structure(list(header = nb_header(object),
                 best = paste("Best strategy:", segments_text(e$threshold[o], best[o])),
                 largest = largest, table = tab[rows, , drop = FALSE],
                 all_rows = is.null(thresholds)),
            class = "summary.voipred_nb")
}

#' @export
print.summary.voipred_nb <- function(x, n = 10, ...) {
  cat(x$header, "\n", x$best, "\n", x$largest, "\n\n", sep = "")
  print_rows(x$table, if (x$all_rows) n else Inf,
             "summary(x, thresholds = ) to choose", digits = 3)
  invisible(x)
}

#' @rdname summary.voipred_nb
#' @export
as.data.frame.voipred_nb <- function(x, row.names = NULL, optional = FALSE, ...)
  x$estimates

# ---- value of information ----------------------------------------------------

#' Summarise value-of-information results
#'
#' `print()` shows the results at the first thresholds; `summary()` gives the
#' main findings:
#' * [p_useful()]: the thresholds where P(useful) is above
#'   0.95, between 0.5 and 0.95, and below 0.5;
#' * [p_best()]: the most likely best strategy at each threshold;
#' * [evpi()]: where EVPI is positive and where it is largest, and the
#'   strategy that is best now (highest expected net benefit) with its
#'   expected net benefit;
#' * [evsi()]: the share of EVPI that a future study of each size recovers
#'   (EVSI / EVPI), at the thresholds with the largest EVPI.
#'
#' `as.data.frame()` returns the full results (the `estimates`, one row per
#' threshold, and model, strategy or future sample size).
#'
#' @param object,x A result of [p_useful()], [p_best()], [evpi()] or [evsi()].
#' @param thresholds Optional thresholds to show (among those of the result).
#'   Default: all (the first 10 are printed); for `evsi()`, the 3 with the
#'   largest EVPI.
#' @param n Number of thresholds printed. Default 10.
#' @param row.names,optional Not used (as in [as.data.frame()]).
#' @param ... Not used.
#' @return `summary()`: an object of class `summary.<class>`, printed as shown.
#'   `as.data.frame()`: the `estimates` data frame.
#' @seealso [p_useful()], [p_best()], [evpi()], [evsi()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse,
#'             family = binomial, data = df_dev_binary)
#' df_val_binary$p1 <- predict(fit1, newdata = df_val_binary, type = "response")
#' fit2 <- glm(day30 ~ age, family = binomial, data = df_dev_binary)
#' df_val_binary$p2 <- predict(fit2, newdata = df_val_binary, type = "response")
#' z <- seq(0.01, 0.3, by = 0.01)
#' res1 <- nb(data = df_val_binary,
#'            predicted_risk = df_val_binary$p1,
#'            outcome = "day30",
#'            thresholds = z,
#'            label = "Full")
#' res2 <- nb(data = df_val_binary,
#'            predicted_risk = df_val_binary$p2,
#'            outcome = "day30",
#'            thresholds = z,
#'            label = "Age only")
#' set.seed(1)
#' e <- evpi(res1, res2)
#' e
#' summary(e)
#' set.seed(1)
#' summary(p_useful(res1))
#' set.seed(1)
#' summary(evsi(res1, res2, future_n = c(500, 1000, 2000)))
#' @name summary.voipred_voi
NULL

#' @rdname summary.voipred_voi
#' @export
print.voipred_evpi <- function(x, n = 10, ...) {
  cat(voi_header("EVPI", x), "\n\n", sep = "")
  print_rows(x$estimates[, c("threshold", "EVPI", "rEVPI")], n,
             "use summary() or as.data.frame()", digits = 4)
  invisible(x)
}

#' @rdname summary.voipred_voi
#' @export
summary.voipred_evpi <- function(object, thresholds = NULL, ...) {
  e <- object$estimates
  enb <- as.matrix(object$ENB[, -1, drop = FALSE])
  best_now <- colnames(enb)[apply(enb, 1, which.max)]
  pos <- is.finite(e$EVPI) & e$EVPI > 1e-12
  line <- if (any(pos)) {
    j <- which(pos)[which.max(e$EVPI[pos])]
    o <- order(e$threshold)
    sprintf("EVPI > 0 at %s; largest %s at %s",
            segments_of(e$threshold[o], pos[o]), fmt_num(e$EVPI[j]),
            fmt_zs(e$threshold)[j])
  } else "EVPI is 0 at every threshold: the best strategy is certain"
  tab <- data.frame(threshold = e$threshold, EVPI = e$EVPI, rEVPI = e$rEVPI,
                    best_now = best_now, ENB_best_now = apply(enb, 1, max),
                    check.names = FALSE)
  names(tab)[4:5] <- c("best now", "ENB best now")
  rows <- select_rows(e$threshold, thresholds)
  structure(list(header = voi_header("EVPI", object), lines = line,
                 table = tab[rows, , drop = FALSE], all_rows = is.null(thresholds)),
            class = "summary.voipred_evpi")
}

#' @rdname summary.voipred_voi
#' @export
print.voipred_evsi <- function(x, n = 10, ...) {
  cat(voi_header("EVSI", x), "\n\n", sep = "")
  print_rows(evsi_wide(x, function(e) e$EVSI), n, "use summary() or as.data.frame()",
             digits = 4)
  invisible(x)
}

#' @rdname summary.voipred_voi
#' @export
summary.voipred_evsi <- function(object, thresholds = NULL, ...) {
  zs <- pick_thresholds(object, thresholds, n_default = 3, max_all = 3)
  share <- evsi_wide(object, function(e)
    ifelse(is.finite(e$EVPI) & e$EVPI > 1e-12,
           paste0(round(100 * e$EVSI / e$EVPI), "%"), "-"))
  share <- share[vapply(share$threshold, function(t) any(near(zs, t)), logical(1)), ,
                 drop = FALSE]
  structure(list(header = voi_header("EVSI", object),
                 lines = if (is.null(thresholds))
                   "EVSI / EVPI at the thresholds with the largest EVPI:"
                 else "EVSI / EVPI at the chosen thresholds:",
                 table = share, all_rows = FALSE, given = !is.null(thresholds)),
            class = "summary.voipred_evsi")
}

#' @rdname summary.voipred_voi
#' @export
print.voipred_p_useful <- function(x, n = 10, ...) {
  cat(voi_header("P(useful)", x), "\n\n", sep = "")
  print_rows(wide(x$estimates, "model", "p_useful", x$models), n,
             "use summary() or as.data.frame()", digits = 3, fixed = TRUE)
  invisible(x)
}

#' @rdname summary.voipred_voi
#' @export
summary.voipred_p_useful <- function(object, ...) {
  e <- object$estimates
  bands <- t(vapply(object$models, function(m) {
    d <- e[e$model == m, ]
    d <- d[order(d$threshold), ]
    p <- d$p_useful
    c(segments_of(d$threshold, p > 0.95, "-"),
      segments_of(d$threshold, p >= 0.5 & p <= 0.95, "-"),
      segments_of(d$threshold, p < 0.5, "-"))
  }, character(3)))
  tab <- data.frame(model = object$models, bands, check.names = FALSE)
  names(tab) <- c("model", "P > 0.95", "0.5 to 0.95", "P < 0.5")
  structure(list(header = voi_header("P(useful)", object),
                 lines = "Thresholds where the probability that the model is useful is:",
                 table = tab, all_rows = FALSE),
            class = "summary.voipred_p_useful")
}

#' @rdname summary.voipred_voi
#' @export
print.voipred_p_best <- function(x, n = 10, ...) {
  cat(voi_header("P(best)", x), "\n\n", sep = "")
  print_rows(wide(x$estimates, "strategy", "p_best", x$strategies), n,
             "use summary() or as.data.frame()", digits = 3, fixed = TRUE)
  invisible(x)
}

#' @rdname summary.voipred_voi
#' @export
summary.voipred_p_best <- function(object, thresholds = NULL, ...) {
  w <- wide(object$estimates, "strategy", "p_best", object$strategies)
  pm <- as.matrix(w[, -1, drop = FALSE])
  j <- apply(pm, 1, which.max)
  tab <- data.frame(threshold = w$threshold,
                    most = colnames(pm)[j],
                    p = pm[cbind(seq_len(nrow(pm)), j)], check.names = FALSE)
  names(tab)[2:3] <- c("most likely best", "P(best)")
  rows <- select_rows(w$threshold, thresholds)
  structure(list(header = voi_header("P(best)", object), lines = NULL,
                 table = tab[rows, , drop = FALSE], all_rows = is.null(thresholds)),
            class = "summary.voipred_p_best")
}

#' @export
print.summary.voipred_evpi <- function(x, n = 10, ...) print_voi_summary(x, n, 4)

#' @export
print.summary.voipred_evsi <- function(x, n = 10, ...) print_voi_summary(x, n, 4)

#' @export
print.summary.voipred_p_useful <- function(x, n = 10, ...) print_voi_summary(x, n, 3)

#' @export
print.summary.voipred_p_best <- function(x, n = 10, ...) print_voi_summary(x, n, 3, TRUE)

#' @rdname summary.voipred_voi
#' @export
as.data.frame.voipred_evpi <- function(x, row.names = NULL, optional = FALSE, ...)
  x$estimates

#' @rdname summary.voipred_voi
#' @export
as.data.frame.voipred_evsi <- function(x, row.names = NULL, optional = FALSE, ...)
  x$estimates

#' @rdname summary.voipred_voi
#' @export
as.data.frame.voipred_p_useful <- function(x, row.names = NULL, optional = FALSE, ...)
  x$estimates

#' @rdname summary.voipred_voi
#' @export
as.data.frame.voipred_p_best <- function(x, row.names = NULL, optional = FALSE, ...)
  x$estimates

# ---- internal ----------------------------------------------------------------

# The first line of print() and summary() of nb().
nb_header <- function(x) {
  o <- x$outcome
  paste0("Net benefit of ", x$label, ": ", outcome_text(o, x$weights),
         if (o$type == "binary") {
           if (is.null(x$weights))
             sprintf(", %s events (%s)", fmt_count(o$events), fmt_pct(o$prevalence))
           else sprintf(", weighted prevalence %s", fmt_pct(o$prevalence))
         } else sprintf(", risk by horizon %s; censoring %s", fmt_pct(o$prevalence),
                        censoring_text(x)),
         "; prior ", prior_text(x$prior))
}

outcome_text <- function(o, weights = NULL) {
  w <- if (!is.null(weights)) " (case weights)" else ""
  switch(o$type,
    binary = paste0("binary, n = ", o$n, w),
    survival = paste0("survival, horizon ", format(o$horizon), ", n = ", o$n, w),
    competing = paste0("competing risks (cause ", o$cause, "), horizon ",
                       format(o$horizon), ", n = ", o$n, w))
}

censoring_text <- function(x) {
  if (!identical(x$outcome$censoring, "cipcw")) return("KM")
  m <- x$cipcw$control$model
  if (is.null(m)) "cIPCW (given weights)" else paste0("cIPCW (", class(m)[1], ")")
}

# The first line of print() and summary() of the value-of-information results.
voi_header <- function(what, x) {
  paste0(what, ": ", paste(x$models, collapse = ", "), "; ",
         outcome_text(x$outcome, if (isTRUE(x$case_weights)) TRUE),
         if (x$outcome$type == "binary") paste0(", prior ", prior_text(x$prior)), "; ",
         x$n_sim, if (x$draws == "posterior") " posterior draws" else
           " Bayesian bootstrap draws")
}

print_voi_summary <- function(x, n, digits, fixed = FALSE) {
  cat(x$header, "\n", sep = "")
  if (length(x$lines)) cat(x$lines, sep = "\n")
  cat("\n")
  print_rows(x$table, if (x$all_rows) n else Inf, "summary(x, thresholds = ) to choose",
             digits = digits, fixed = fixed)
  invisible(x)
}

# Print the first n rows of a table, and how many more there are. Numbers
# with `digits` significant digits (`fixed`: with `digits` decimals, for
# probabilities), never in scientific notation; thresholds with a common
# number of decimals.
print_rows <- function(tab, n, hint, digits = 4, fixed = FALSE) {
  m <- nrow(tab)
  shown <- if (is.finite(n) && m > n) tab[seq_len(n), , drop = FALSE] else tab
  for (k in names(shown)) {
    v <- shown[[k]]
    if (!is.numeric(v)) next
    shown[[k]] <- if (k == "threshold") fmt_zs(v) else fmt_values(v, digits, fixed)
  }
  print(shown, row.names = FALSE, right = TRUE)
  if (nrow(shown) < m)
    cat("# ", m - nrow(shown), " more ", if (m - nrow(shown) == 1) "row" else "rows",
        "; ", hint, "\n", sep = "")
}

# A long table (threshold, key, value) as one row per threshold and one column
# per key.
wide <- function(e, key, value, keys) {
  z <- unique(e$threshold)
  out <- data.frame(threshold = z)
  for (k in keys) out[[k]] <- e[[value]][e[[key]] == k][match(z, e$threshold[e[[key]] == k])]
  out
}

# EVSI results with one column per future sample size, after EVPI.
evsi_wide <- function(x, value) {
  e <- x$estimates
  z <- unique(e$threshold)
  out <- data.frame(threshold = z, EVPI = e$EVPI[match(z, e$threshold)])
  for (N in x$future_n) {
    d <- e[e$future_n == N, ]
    out[[paste0("N=", N)]] <- value(d)[match(z, d$threshold)]
  }
  out
}

# Rows of the thresholds asked for (all if NULL).
select_rows <- function(z, thresholds) {
  if (is.null(thresholds)) return(seq_along(z))
  if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds))
    stop("`thresholds` must be one or more of the thresholds of the result.",
         call. = FALSE)
  rows <- vapply(thresholds, function(t) {
    i <- which(near(z, t))
    if (!length(i))
      stop("Threshold ", format(t), " is not among the thresholds of the result.",
           call. = FALSE)
    i[1]
  }, integer(1))
  rows
}

# Best strategy at each threshold of the estimates of nb().
best_strategy <- function(e) {
  ifelse(is.na(e$delta_NB), NA_character_,
         ifelse(e$delta_NB > 1e-12, "model",
                ifelse(e$NB_all > 1e-12, "treat all", "treat none")))
}

# Posterior mean sensitivity and specificity at each threshold (binary).
# KM (dcurves): sensitivity and specificity implied by the Kaplan-Meier net
# benefit, the share treated q and the risk by the horizon F:
# with F+ the risk among the treated, NB = q F+ - q (1 - F+) z/(1-z), so
# q F+ = (NB + q z/(1-z)) / (1 + z/(1-z)); Se = q F+ / F and
# Sp = 1 - q (1 - F+) / (1 - F). Then F Se - (1 - F)(1 - Sp) z/(1-z) = NB.
sens_spec_km <- function(x) {
  e <- x$estimates
  F <- x$outcome$prevalence
  w <- if (is.null(x$weights)) rep(1, length(x$risk)) else x$weights
  q <- vapply(e$threshold, function(z) sum(w * treated_at(x$risk, z)) / sum(w),
              numeric(1))                                  # share treated
  odds <- e$threshold / (1 - e$threshold)
  qf <- (e$NB_model + q * odds) / (1 + odds)              # q F+
  se <- if (F > 0) qf / F else rep(NA_real_, length(q))
  sp <- if (F < 1) 1 - (q - qf) / (1 - F) else rep(NA_real_, length(q))
  clamp <- function(v) ifelse(is.finite(v), pmin(pmax(v, 0), 1), NA)
  list(Se = clamp(se), Sp = clamp(sp))
}

# cIPCW: the estimator's own Se and Sp, IPCW-weighted proportions (the net
# benefit is prevalence * Se - (1 - prevalence) * (1 - Sp) * z/(1-z)).
sens_spec_cipcw <- function(x) {
  v <- voi_info(x)
  r <- nb_core_cipcw(v$time, v$status, v$k, v$risk, v$thresholds, v$horizon,
                     v$weight_event, v$weight_horizon, v$weights)
  list(Se = r$Se, Sp = r$Sp)
}

sens_spec <- function(x) {
  cn <- x$counts
  p <- x$prior
  se <- (p$sensitivity[1] + cn$TP) / (sum(p$sensitivity) + cn$TP + cn$FN)
  sp <- (p$specificity[1] + cn$TN) / (sum(p$specificity) + cn$TN + cn$FP)
  list(Se = ifelse(is.finite(se), se, NA), Sp = ifelse(is.finite(sp), sp, NA))
}

# "a 0.01-0.26, b 0.27-0.30": runs of the same value over sorted thresholds.
segments_text <- function(z, value) {
  zl <- fmt_zs(z)
  r <- rle(ifelse(is.na(value), "NA", value))
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  paste(vapply(seq_along(r$values), function(i)
    paste(r$values[i], range_text(zl[starts[i]], zl[ends[i]])), ""), collapse = ", ")
}

# The thresholds (sorted) where `flag` is TRUE, as ranges: "0.02-0.10, 0.15".
segments_of <- function(z, flag, none = "none") {
  flag[is.na(flag)] <- FALSE
  if (!any(flag)) return(none)
  zl <- fmt_zs(z)
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  keep <- which(r$values)
  paste(vapply(keep, function(i) range_text(zl[starts[i]], zl[ends[i]]), ""), collapse = ", ")
}

range_text <- function(a, b) if (a == b) a else paste0(a, "-", b)

# Thresholds written with the same number of decimals ("0.30", not "0.3").
fmt_zs <- function(z) format(signif(z, 4), trim = TRUE, scientific = FALSE)

# One number of decimals per column: `digits` significant digits for its
# largest value (`fixed`: `digits` decimals).
fmt_values <- function(v, digits, fixed = FALSE) {
  dec <- digits
  if (!fixed) {
    m <- suppressWarnings(max(abs(v[is.finite(v)])))
    dec <- if (!is.finite(m) || m == 0) digits else
      min(10, max(0, digits - floor(log10(m)) - 1))
  }
  out <- formatC(v, format = "f", digits = dec)
  out[is.na(v)] <- "NA"
  out[is.infinite(v)] <- ifelse(v[is.infinite(v)] > 0, "Inf", "-Inf")
  trimws(out)
}

fmt_num <- function(v) format(signif(v, 3), trim = TRUE, scientific = FALSE)
fmt_pct <- function(p) paste0(format(round(100 * p, 1), nsmall = 1, trim = TRUE), "%")
fmt_count <- function(k) format(round(k), trim = TRUE)
