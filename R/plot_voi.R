#' Plot the expected value of perfect information
#'
#' Plots EVPI (or relative EVPI) against the threshold, from one or more
#' results of [evpi()], e.g. for different validation samples or sets of
#' models. To draw EVPI under the decision curve, use
#' `plot_dca(..., evpi = )`.
#'
#' @param x A result of [evpi()].
#' @param ... More results of [evpi()] to draw on the same plot. A named
#'   result is labelled with its name.
#' @param type `"evpi"` (default) or `"relative"` (relative EVPI, with a dashed
#'   line at 1; infinite values are left out).
#' @param labels Optional names of the curves, one per result. Default:
#'   `"EVPI"` for one result, the argument names (or `"EVPI 1"`, `"EVPI 2"`,
#'   ...) for several.
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [evpi()], [plot_dca()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, "day30", thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' e <- evpi(res, n_sim = 200)
#' plot_evpi(e)
#' plot_evpi(e, type = "relative")
#' @export
plot_evpi <- function(x, ..., type = c("evpi", "relative"), labels = NULL,
                      interactive = FALSE) {
  type <- match.arg(type)
  check_interactive(interactive)
  res <- collect_results(x, list(...), labels, "voipred_evpi", "evpi()", "plot_evpi()",
                         "`type`, `labels` and `interactive`", "EVPI")
  col <- if (type == "evpi") "EVPI" else "rEVPI"
  d <- do.call(rbind, Map(function(r, l) curve_data(r$threshold, l, r[[col]]),
                          res, names(res)))
  d$label <- factor(d$label, levels = names(res))
  sp <- new_spec(d, stats::setNames(series_colors(length(res)), names(res)),
                 "Threshold Probability",
                 if (type == "evpi") "EVPI" else "Relative EVPI",
                 hline = if (type == "relative") 1, fmt = if (type == "evpi") ".5f" else ".3f")
  render_plot(list(sp), interactive)
}

#' Plot the expected value of sample information
#'
#' Plots EVSI from [evsi()] against the future sample size, one curve per
#' threshold, starting at 0 for no new data, with that threshold's EVPI as a
#' dashed line (the limit as the future sample grows), as in Sadatsafavi et
#' al. (2025). To draw EVSI under the decision curve, use
#' `plot_dca(..., evsi = )`.
#'
#' @param x A result of [evsi()].
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [evsi()], [plot_dca()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, "day30", thresholds = c(0.02, 0.05, 0.1))
#' set.seed(1)
#' o <- evsi(res, future_n = c(250, 500, 1000, 2000), n_sim = 200)
#' plot_evsi(o)
#' @export
plot_evsi <- function(x, interactive = FALSE) {
  check_interactive(interactive)
  if (!inherits(x, "voipred_evsi"))
    stop("plot_evsi() needs a result of evsi().", call. = FALSE)
  zs <- unique(x$threshold)
  lev <- paste0("Threshold ", format(100 * zs, trim = TRUE), "%")
  ns <- sort(unique(x$future_n))
  d <- do.call(rbind, lapply(seq_along(zs), function(i) {
    e <- x[x$threshold == zs[i], ]
    e <- e[order(e$future_n), ]
    curve_data(c(0, e$future_n), lev[i], c(0, e$EVSI))
  }))
  ref <- do.call(rbind, lapply(seq_along(zs), function(i) {
    evpi <- x$EVPI[x$threshold == zs[i]][1]
    curve_data(c(0, max(ns)), lev[i], rep(evpi, 2))
  }))
  d$label <- factor(d$label, levels = lev)
  ref$label <- factor(ref$label, levels = lev)
  sp <- new_spec(d, stats::setNames(series_colors(length(lev)), lev),
                 "Future sample size", "EVSI", ref = ref, percent = FALSE, fmt = ".5f")
  render_plot(list(sp), interactive)
}

#' Plot the probability that each model is useful
#'
#' Plots P(useful), from [p_useful()], against the threshold, one curve per
#' model.
#'
#' @param x A result of [p_useful()].
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [p_useful()], [plot_p_best()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, "day30", thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' plot_p_useful(p_useful(res, n_sim = 200))
#' @export
plot_p_useful <- function(x, interactive = FALSE) {
  check_interactive(interactive)
  if (!inherits(x, "voipred_p_useful"))
    stop("plot_p_useful() needs a result of p_useful().", call. = FALSE)
  lev <- unique(x$model)
  d <- curve_data(x$threshold, x$model, x$p_useful)
  d$label <- factor(d$label, levels = lev)
  # the models keep the colours they have in plot_dca()
  cols <- default_hues(length(lev) + 2)[-(1:2)]
  sp <- new_spec(d, stats::setNames(cols, lev), "Threshold Probability", "P(useful)",
                 ylim = c(0, 1), fmt = ".3f")
  render_plot(list(sp), interactive)
}

#' Plot the probability that each strategy is the best
#'
#' Plots P(best), from [p_best()], against the threshold, one curve per
#' strategy (each model, treat all and treat none), in the colours of
#' [plot_dca()].
#'
#' @param x A result of [p_best()].
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, or a plotly widget with `interactive = TRUE`.
#' @seealso [p_best()], [plot_p_useful()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:5500, ]
#' fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' res <- nb(fit, val, "day30", thresholds = seq(0.01, 0.3, by = 0.01))
#' set.seed(1)
#' plot_p_best(p_best(res, n_sim = 200))
#' @export
plot_p_best <- function(x, interactive = FALSE) {
  check_interactive(interactive)
  if (!inherits(x, "voipred_p_best"))
    stop("plot_p_best() needs a result of p_best().", call. = FALSE)
  strat <- unique(x$strategy)
  models <- setdiff(strat, c("Treat All", "Treat None"))
  lev <- c("Treat All", "Treat None", models)
  d <- curve_data(x$threshold, x$strategy, x$p_best)
  d$label <- factor(d$label, levels = lev)
  sp <- new_spec(d, stats::setNames(default_hues(length(lev)), lev),
                 "Threshold Probability", "P(best)", ylim = c(0, 1), fmt = ".3f")
  render_plot(list(sp), interactive)
}

# ---- internal ----------------------------------------------------------------

series_colors <- function(n) if (n == 1) "#000000" else default_hues(n)

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
      stop("`labels` must be ", length(res), " strings, one per result.", call. = FALSE)
    lab <- labels
  }
  if (anyDuplicated(lab))
    stop("Two results have the same label (\"", lab[anyDuplicated(lab)], "\"); give ",
         "each a different one with `labels`.", call. = FALSE)
  names(res) <- lab
  res
}
