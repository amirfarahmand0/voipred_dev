#' Plot the expected value of perfect information
#'
#' The plot of results of [evpi()]: EVPI (or relative EVPI) against the
#' threshold, for one or more results, e.g. for different validation samples
#' or sets of models; or, for one result, P(useful) or P(best) from the same
#' draws (as [plot.voipred_p_useful()] and [plot.voipred_p_best()] draw
#' them). To draw EVPI under the decision curve, use
#' `plot(<nb_val result>, evpi = )` ([plot.voipred_nb()]).
#'
#' @param x A result of [evpi()].
#' @param ... More results of [evpi()] to draw on the same plot. A result given
#'   here with a name is labelled with its name.
#' @param type `"evpi"` (default), `"relative"` (relative EVPI, with a dashed
#'   line at 1; infinite values are left out), `"p_useful"` (P(useful) of each
#'   model) or `"p_best"` (P(best) of each strategy); the last two for one
#'   result.
#' @param labels Optional names of the curves, one per result, in order
#'   (`type = "p_useful"`: one per model).
#'   Default: `"EVPI"` for one result, the argument names (or `"EVPI 1"`,
#'   `"EVPI 2"`, ...) for several.
#' @param colors Optional colours of the curves, instead of the defaults:
#'   named by curve label, or unnamed, one per curve in order (with
#'   `type = "p_useful"` or `"p_best"`, as in [plot.voipred_p_useful()] and
#'   [plot.voipred_p_best()]).
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [evpi()], [plot.voipred_nb()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb_val(data = df_val_binary,
#'               predicted_risk = df_val_binary$p,
#'               outcome = "day30",
#'               thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' e <- evpi(res, n_sim = 500)
#' plot(e)
#' plot(e, type = "relative")
#' plot(e, type = "p_useful")
#' plot(e, type = "p_best")
#' @export
plot.voipred_evpi <- function(x, ..., type = c("evpi", "relative", "p_useful", "p_best"),
                              labels = NULL, colors = NULL, interactive = FALSE) {
  type <- match.arg(type)
  check_interactive(interactive)
  if (type %in% c("p_useful", "p_best")) {
    # P(useful) or P(best) from the draws of evpi()
    if (length(list(...)))
      stop("plot(type = \"", type, "\") draws one result of evpi().", call. = FALSE)
    if (type == "p_best" && !is.null(labels))
      stop("`labels` is not used with type = \"p_best\".", call. = FALSE)
    obj <- from_evpi(x, type, list(), TRUE)
    return(if (type == "p_useful")
             plot(obj, labels = labels, colors = colors, interactive = interactive)
           else plot(obj, colors = colors, interactive = interactive))
  }
  res <- collect_results(x, list(...), labels, "voipred_evpi", "evpi()", "plot()",
                         "`type`, `labels`, `colors` and `interactive`", "EVPI")
  col <- if (type == "evpi") "EVPI" else "rEVPI"
  d <- do.call(rbind, Map(function(r, l) curve_data(r$estimates$threshold, l,
                                                     r$estimates[[col]]),
                          res, names(res)))
  d$label <- factor(d$label, levels = names(res))
  sp <- new_spec(d, stats::setNames(series_colors(length(res)), names(res)),
                 "Threshold Probability",
                 if (type == "evpi") "EVPI" else "Relative EVPI",
                 hline = if (type == "relative") 1, fmt = if (type == "evpi") ".5f" else ".3f")
  render_plot(recolor(list(sp), colors, names(res), "result"), interactive)
}

#' Plot the expected value of sample information
#'
#' The plot of a result of [evsi()]: EVSI against the future sample size, one
#' curve per threshold, starting at 0 for no new data, with that threshold's
#' EVPI as a dashed line (the limit as the future sample grows), as in
#' Sadatsafavi et al. (2025). To draw EVSI under the decision curve, use
#' `plot(<nb_val result>, evsi = )` ([plot.voipred_nb()]).
#'
#' @param x A result of [evsi()].
#' @param ... Not used.
#' @param thresholds The thresholds to draw, among those of `x`. Default: all
#'   of them if there are at most 5, otherwise the 4 with the largest EVPI.
#' @param colors Optional colours of the curves (one per threshold), instead
#'   of the defaults: named by curve label (e.g. `"Threshold 10%"`), or
#'   unnamed, one per curve in order.
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [evsi()], [plot.voipred_nb()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb_val(data = df_val_binary,
#'               predicted_risk = df_val_binary$p,
#'               outcome = "day30",
#'               thresholds = c(0.02, 0.05, 0.1))
#' set.seed(1)
#' o <- evsi(res, future_n = c(250, 500, 1000, 2000), n_sim = 500)
#' plot(o)
#' plot(o, thresholds = 0.05)
#' @export
plot.voipred_evsi <- function(x, ..., thresholds = NULL, colors = NULL, interactive = FALSE) {
  check_interactive(interactive)
  check_no_dots(list(...), "plot()", "`thresholds`, `colors` and `interactive`")
  o <- x$estimates
  zs <- pick_thresholds(x, thresholds, n_default = 4, max_all = 5)
  lev <- paste0("Threshold ", format(100 * zs, trim = TRUE), "%")
  ns <- sort(unique(o$future_n))
  d <- do.call(rbind, lapply(seq_along(zs), function(i) {
    e <- o[near(o$threshold, zs[i]), ]
    e <- e[order(e$future_n), ]
    curve_data(c(0, e$future_n), lev[i], c(0, e$EVSI))
  }))
  ref <- do.call(rbind, lapply(seq_along(zs), function(i) {
    evpi <- o$EVPI[near(o$threshold, zs[i])][1]
    curve_data(c(0, max(ns)), lev[i], rep(evpi, 2))
  }))
  d$label <- factor(d$label, levels = lev)
  ref$label <- factor(ref$label, levels = lev)
  sp <- new_spec(d, stats::setNames(series_colors(length(lev)), lev),
                 "Future sample size", "EVSI", ref = ref, percent = FALSE, fmt = ".5f")
  render_plot(recolor(list(sp), colors, lev, "threshold"), interactive)
}

#' Plot the probability that a model is useful
#'
#' The plot of results of [p_useful()]: P(useful) against the threshold, one
#' curve per model of each result, in the colours of [plot.voipred_nb()].
#'
#' @param x A result of [p_useful()].
#' @param ... More results of [p_useful()] to draw on the same plot. A result
#'   with one model given here with a name is labelled with its name.
#' @param labels Optional names of the curves, one per model, in order.
#'   Default: the models' labels (or the argument names).
#' @param colors Optional colours of the curves, instead of the defaults:
#'   named by label, e.g. `c(Full = "darkred")`, or unnamed, one per curve in
#'   order.
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [p_useful()], [plot.voipred_p_best()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb_val(data = df_val_binary,
#'               predicted_risk = df_val_binary$p,
#'               outcome = "day30",
#'               thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' pu <- p_useful(res, n_sim = 500)
#' plot(pu)
#'
#' # two models
#' fit_age <- glm(day30 ~ age, family = binomial, data = df_dev_binary)
#' df_val_binary$p_age <- predict(fit_age, newdata = df_val_binary, type = "response")
#' both <- nb_val(data = df_val_binary,
#'                predicted_risk = list(Full = df_val_binary$p, `Age only` = df_val_binary$p_age),
#'                outcome = "day30",
#'                thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' plot(p_useful(both, n_sim = 500))
#' @export
plot.voipred_p_useful <- function(x, ..., labels = NULL, colors = NULL, interactive = FALSE) {
  check_interactive(interactive)
  dots <- list(...)
  collect_results(x, dots, NULL, "voipred_p_useful", "p_useful()", "plot()",
                  "`labels`, `colors` and `interactive`", "P(useful)")
  items <- c(list(x), dots)
  nms <- c("", if (is.null(names(dots))) rep("", length(dots)) else names(dots))
  # one curve per model of each result, labelled by the model, or by the
  # argument name for a result with one model
  curves <- list()
  lev <- character()
  for (i in seq_along(items)) {
    r <- items[[i]]
    if (nzchar(nms[i]) && length(r$models) > 1)
      stop("A name relabels a result with one model, but the result named \"", nms[i],
           "\" has ", length(r$models), "; use `labels`.", call. = FALSE)
    for (m in r$models) {
      curves[[length(curves) + 1]] <- r$estimates[r$estimates$model == m, ]
      lev <- c(lev, if (nzchar(nms[i])) nms[i] else m)
    }
  }
  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != length(lev) || anyNA(labels))
      labels_error(length(lev), "curve")
    lev <- labels
  }
  if (anyDuplicated(lev))
    stop("Two curves have the same label (\"", lev[anyDuplicated(lev)], "\"); give ",
         "each a different one with `labels`.", call. = FALSE)
  d <- do.call(rbind, Map(function(e, l) curve_data(e$threshold, l, e$p_useful),
                          curves, lev))
  d$label <- factor(d$label, levels = lev)
  # the models keep the colours they have in the decision curve
  cols <- default_hues(length(lev) + 2)[-(1:2)]
  sp <- new_spec(d, stats::setNames(cols, lev), "Threshold Probability", "P(useful)",
                 ylim = c(0, 1), fmt = ".3f")
  render_plot(recolor(list(sp), colors, lev, "curve"), interactive)
}

#' Plot the probability that each strategy is the best
#'
#' The plot of a result of [p_best()]: P(best) against the threshold, one curve
#' per strategy (each model, treat all and treat none), in the colours of
#' [plot.voipred_nb()].
#'
#' @param x A result of [p_best()].
#' @param ... Not used.
#' @param colors Optional colours of the curves, instead of the defaults:
#'   named by strategy, e.g. `c(Full = "darkred", "Treat All" = "grey50")`, or
#'   unnamed, one per model in order (Treat All and Treat None keep theirs).
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [p_best()], [plot.voipred_p_useful()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse,
#'            family = binomial, data = df_dev_binary)
#' df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")
#' res <- nb_val(data = df_val_binary,
#'               predicted_risk = df_val_binary$p,
#'               outcome = "day30",
#'               thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' plot(p_best(res, n_sim = 500))
#' @export
plot.voipred_p_best <- function(x, ..., colors = NULL, interactive = FALSE) {
  check_interactive(interactive)
  check_no_dots(list(...), "plot()", "`colors` and `interactive`")
  e <- x$estimates
  lev <- c("Treat All", "Treat None", x$models)
  d <- curve_data(e$threshold, e$strategy, e$p_best)
  d$label <- factor(d$label, levels = lev)
  sp <- new_spec(d, stats::setNames(default_hues(length(lev)), lev),
                 "Threshold Probability", "P(best)", ylim = c(0, 1), fmt = ".3f")
  render_plot(recolor(list(sp), colors, x$models, "model"), interactive)
}

# ---- internal ----------------------------------------------------------------

series_colors <- function(n) if (n == 1) "#000000" else default_hues(n)

near <- function(a, b) abs(a - b) < 1e-9

check_no_dots <- function(dots, fun, args) {
  if (length(dots))
    stop("Unknown argument to ", fun, ". Its arguments are ", args, ".", call. = FALSE)
}

labels_error <- function(n, what = "result")
  stop("`labels` must be ", if (n == 1) "one string" else paste(n, "strings"),
       ", one per ", what, ".", call. = FALSE)

# The thresholds of an evsi() result to show: the ones asked for, or by default
# all (at most `max_all`) or the `n_default` with the largest EVPI.
pick_thresholds <- function(x, thresholds, n_default, max_all) {
  o <- x$estimates
  zs <- unique(o$threshold)
  if (!is.null(thresholds)) {
    if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds))
      stop("`thresholds` must be one or more of the thresholds of the result.",
           call. = FALSE)
    ok <- vapply(thresholds, function(t) any(near(zs, t)), logical(1))
    if (!all(ok))
      stop("Threshold ", format(thresholds[!ok][1]), " is not among the thresholds ",
           "of the result.", call. = FALSE)
    return(zs[vapply(thresholds, function(t) which(near(zs, t))[1], integer(1))])
  }
  if (length(zs) <= max_all) return(zs)
  evpi <- o$EVPI[match(zs, o$threshold)]
  top <- order(-ifelse(is.na(evpi), -Inf, evpi))[seq_len(n_default)]
  sort(zs[top])
}

# One or more results of class `cls`, named by `labels`, the argument names,
# or `single` / "<single> 1", "<single> 2", ...
collect_results <- function(x, more, labels, cls, made_by, fun, args, single) {
  nms <- names(more)
  if (is.null(nms)) nms <- rep("", length(more))
  bad <- !vapply(more, inherits, logical(1), cls)
  if (any(bad & nzchar(nms)))
    stop("Unknown argument `", nms[bad & nzchar(nms)][1], "`. The arguments of ",
         fun, " are ", args, "; other arguments must be results of ", made_by, ".",
         call. = FALSE)
  if (!inherits(x, cls) || any(bad))
    stop(fun, " needs results of ", made_by, ".", call. = FALSE)
  res <- c(list(x), more)
  lab <- if (length(res) == 1) single else paste(single, seq_along(res))
  lab[-1][nzchar(nms)] <- nms[nzchar(nms)]
  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != length(res) || anyNA(labels))
      labels_error(length(res))
    lab <- labels
  }
  if (anyDuplicated(lab))
    stop("Two results have the same label (\"", lab[anyDuplicated(lab)], "\"); give ",
         "each a different one with `labels`.", call. = FALSE)
  names(res) <- lab
  res
}
