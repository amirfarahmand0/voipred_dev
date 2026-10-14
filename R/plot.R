#' Plot decision curves
#'
#' Plots the results of [nb()]: the decision curve (`type = "nb"`), or the
#' difference in net benefit between the model and the best default strategy
#' (`type = "delta"`). Give several results to draw their curves together.
#' Optionally, the expected value of perfect and sample information (from
#' [evpi()] and [evsi()]) is drawn in a panel below, on the same threshold
#' axis.
#'
#' * `type = "nb"`: the net benefit of each model, of treating everyone
#'   ("Treat All") and of treating no one ("Treat None"), against the
#'   threshold. The y-axis is shown from \eqn{-0.1\pi} to \eqn{\pi}, with
#'   \eqn{\pi} the prevalence (risk by the horizon for survival outcomes), as
#'   in dcurves.
#' * `type = "delta"`: \eqn{\Delta NB = NB_{model} - \max(NB_{all},
#'   NB_{none})} (column `delta_NB` of [nb()]) of each model, with a dashed
#'   line at 0. Above 0 the model is better than both default strategies.
#'
#' Results drawn together must come from the same data, outcome, horizon,
#' cause, censoring method and thresholds (voipred checks this). Thresholds
#' where the net benefit is `NA` are left out. Each curve is named by its
#' `label` (see [nb()]), unless `labels` is given; a model has the same colour
#' in both types of plot.
#'
#' The static plot is a ggplot object (a patchwork object with the
#' value-of-information panel): it is drawn when printed, and can be changed
#' with the usual `+` (e.g. `+ ggplot2::coord_cartesian(xlim = c(0, 0.5))`) or
#' saved with `ggplot2::ggsave()`. With `interactive = TRUE` it is a plotly
#' widget: hovering shows the values of every curve at a threshold, the panels
#' zoom together, and clicking a legend entry hides or shows a curve; save it
#' with `htmlwidgets::saveWidget()`.
#'
#' @param x A result of [nb()].
#' @param ... More results of [nb()] to draw on the same plot. A named result
#'   is labelled with its name.
#' @param type `"nb"` (decision curve, default) or `"delta"` (\eqn{\Delta NB}).
#' @param labels Optional names of the curves, one per result, in order;
#'   default: the results' labels.
#' @param evpi,evsi Optional results of [evpi()] and [evsi()], drawn in a panel
#'   below the decision curve: EVPI, and EVSI for each future sample size,
#'   against the threshold. The static version needs the patchwork package.
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, a patchwork object (with `evpi` or `evsi`), or a
#'   plotly widget (`interactive = TRUE`).
#' @seealso [plot_evpi()], [plot_evsi()], [plot_p_useful()], [plot_p_best()]
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:10000, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' fit2 <- glm(day30 ~ age, family = binomial, data = dev)
#' z <- seq(0, 0.5, by = 0.01)
#' res1 <- nb(fit1, val, "day30", thresholds = z, label = "Full model")
#' res2 <- nb(fit2, val, "day30", thresholds = z, label = "Age only")
#' plot_dca(res1)
#' plot_dca(res1, res2)
#' plot_dca(res1, res2, type = "delta")
#' \donttest{
#' set.seed(1)
#' e <- evpi(res1, n_sim = 200)
#' if (requireNamespace("patchwork", quietly = TRUE)) plot_dca(res1, evpi = e)
#' if (requireNamespace("plotly", quietly = TRUE)) plot_dca(res1, evpi = e, interactive = TRUE)
#' }
#' @export
plot_dca <- function(x, ..., type = c("nb", "delta"), labels = NULL, evpi = NULL,
                     evsi = NULL, interactive = FALSE) {
  type <- match.arg(type)
  check_interactive(interactive)
  res <- collect_nb(x, list(...), labels,
                    args = "`type`, `labels`, `evpi`, `evsi` and `interactive`")
  top <- dca_spec(res, type)
  if (is.null(evpi) && is.null(evsi)) return(render_plot(list(top), interactive))
  bottom <- voi_panel_spec(evpi, evsi)
  render_plot(list(top, bottom), interactive, heights = c(0.62, 0.38))
}

# ---- plot specifications -----------------------------------------------------
# A plot is built from one or more "specs": the curves (columns x, value,
# label), dashed reference curves (`ref`), axis labels (`ylab` for ggplot2,
# `ylab_text` for plotly), colours, y limits, an optional dashed line at
# `hline`, and the number format of the hover text. render_plot() draws them
# with ggplot2 (stacked with patchwork when there are several) or plotly.

new_spec <- function(data, colors, xlab, ylab, ylab_text = ylab, ylim = NULL,
                     hline = NULL, ref = NULL, percent = TRUE, fmt = ".4f") {
  data <- data[is.finite(data$value), , drop = FALSE]
  if (!is.null(ref)) ref <- ref[is.finite(ref$value), , drop = FALSE]
  list(data = data, colors = colors, xlab = xlab, ylab = ylab, ylab_text = ylab_text,
       ylim = ylim, hline = hline, ref = ref, percent = percent, fmt = fmt)
}

dca_spec <- function(res, type) {
  r1 <- res[[1]]
  lab <- names(res)
  # ggplot2's default colours for Treat All, Treat None and the models, so a
  # model has the same colour in both types of plot
  cols <- default_hues(length(lab) + 2)
  if (type == "nb") {
    d <- rbind(curve_data(r1$threshold, "Treat All", r1$NB_all),
               curve_data(r1$threshold, "Treat None", r1$NB_none),
               do.call(rbind, Map(function(r, l) curve_data(r$threshold, l, r$NB_model),
                                  res, lab)))
    d$label <- factor(d$label, levels = c("Treat All", "Treat None", lab))
    prev <- attr(r1, "prevalence")
    new_spec(d, stats::setNames(cols, levels(d$label)), "Threshold Probability",
             "Net Benefit", ylim = c(-0.1 * prev, prev))
  } else {
    d <- do.call(rbind, Map(function(r, l) curve_data(r$threshold, l, r$delta_NB),
                            res, lab))
    d$label <- factor(d$label, levels = lab)
    new_spec(d, stats::setNames(cols[-(1:2)], lab), "Threshold Probability",
             quote(Delta * NB ~ "(" * NB[model] - "max(" * NB[treat ~ all] * "," ~
                     NB[treat ~ none] * "))"),
             ylab_text = "&#916;NB (NB<sub>model</sub> &#8722; max(NB<sub>treat all</sub>, NB<sub>treat none</sub>))",
             hline = 0)
  }
}

# EVPI and EVSI (for each future sample size) against the threshold.
# The y-axis says what is drawn: "EVPI" (EVPI only) or "EVPI / EVSI" (both).
voi_panel_spec <- function(evpi, evsi) {
  if (!is.null(evpi) && !inherits(evpi, "voipred_evpi"))
    stop("`evpi` must be a result of evpi().", call. = FALSE)
  if (!is.null(evsi) && !inherits(evsi, "voipred_evsi"))
    stop("`evsi` must be a result of evsi().", call. = FALSE)
  d <- NULL
  if (!is.null(evpi)) d <- curve_data(evpi$threshold, "EVPI", evpi$EVPI)
  else {
    first <- evsi[evsi$future_n == evsi$future_n[1], ]
    d <- curve_data(first$threshold, "EVPI", first$EVPI)
  }
  lev <- "EVPI"
  if (!is.null(evsi)) {
    ns <- unique(evsi$future_n)
    lev <- c(lev, paste0("EVSI (N = ", format(ns, big.mark = ",", trim = TRUE), ")"))
    for (i in seq_along(ns)) {
      e <- evsi[evsi$future_n == ns[i], ]
      d <- rbind(d, curve_data(e$threshold, lev[i + 1], e$EVSI))
    }
  }
  d$label <- factor(d$label, levels = lev)
  ylab <- if (is.null(evsi)) "EVPI" else "EVPI / EVSI"
  new_spec(d, stats::setNames(c("#000000", default_hues(length(lev) - 1)), lev),
           "Threshold Probability", ylab, ylab_text = ylab, fmt = ".5f")
}

# ---- rendering -----------------------------------------------------------------

# One look for the static (ggplot2) and the interactive (plotly) plots: white
# background, light grid, no frame, 2 px lines, grey text, and the legend on
# the right without a title (as in dcurves).
plot_style <- list(text = "#444444", grid = "#EEEEEE", ref = "#555555", font_size = 12)

render_plot <- function(specs, interactive, heights = NULL) {
  if (interactive) return(render_plotly(specs, heights))
  if (length(specs) == 1) return(render_gg(specs[[1]]))
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("The patchwork package is needed to stack the value-of-information panel ",
         "under the decision curve; install it, or use `interactive = TRUE`.",
         call. = FALSE)
  xr <- range(unlist(lapply(specs, function(sp) sp$data$x)))
  plots <- lapply(seq_along(specs), function(i) {
    p <- render_gg(specs[[i]], xlim = xr)
    # the panels share the x axis: its labels only under the bottom one
    if (i < length(specs))
      p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    p
  })
  # each panel keeps its own legend, on its right
  Reduce(`/`, plots) + patchwork::plot_layout(heights = heights)
}

render_gg <- function(sp, xlim = NULL) {
  st <- plot_style
  p <- ggplot(sp$data, aes(x = .data$x, y = .data$value, color = .data$label))
  if (!is.null(sp$hline))
    p <- p + geom_hline(yintercept = sp$hline, linetype = "dashed", colour = st$ref,
                        linewidth = 0.4)
  if (!is.null(sp$ref) && nrow(sp$ref))
    p <- p + geom_line(data = sp$ref, linetype = "dashed", linewidth = 0.6)
  p <- p + geom_line(linewidth = 0.75) +
    scale_color_manual(values = sp$colors) +
    labs(x = sp$xlab, y = sp$ylab, color = NULL)
  if (!is.null(sp$ylim) || !is.null(xlim))
    p <- p + coord_cartesian(xlim = xlim, ylim = sp$ylim)
  p <- p + scale_x_continuous(labels = if (sp$percent) percent_label else plain_label) +
    scale_y_continuous(labels = plain_label)
  p + theme_minimal(base_size = st$font_size) +
    theme(panel.grid.major = element_line(colour = st$grid, linewidth = 0.4),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = st$text),
          axis.title = element_text(colour = st$text),
          legend.text = element_text(colour = st$text),
          legend.position = "right",
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
}

render_plotly <- function(specs, heights = NULL) {
  st <- plot_style
  K <- length(specs)
  # a legend per panel needs plotly.js 2.23 or later (R plotly 4.11.0); with an
  # older plotly the panels share one legend
  own <- K > 1 && utils::packageVersion("plotly") >= "4.11.0"
  panels <- lapply(seq_len(K), function(i)
    plotly_panel(specs[[i]], i, legend = if (own) legend_id(i)))
  margin <- 0.04
  if (is.null(heights)) heights <- rep(1 / K, K)
  p <- if (K == 1) panels[[1]] else
    plotly::subplot(panels, nrows = K, shareX = TRUE, titleY = TRUE,
                    heights = heights, margin = margin)
  # each legend on the right of its panel, vertically centred on it, as in dcurves
  legends <- lapply(seq_len(K), function(i) {
    top <- 1 - sum(heights[seq_len(i - 1)]) - if (i > 1) margin else 0
    bottom <- 1 - sum(heights[seq_len(i)]) + if (i < K) margin else 0
    list(orientation = "v", x = 1.02, xanchor = "left",
         y = if (own || K == 1) (top + bottom) / 2 else 0.5, yanchor = "middle")
  })
  lay <- if (own) stats::setNames(legends, vapply(seq_len(K), legend_id, ""))
         else list(legend = legends[[1]])
  sp1 <- specs[[1]]
  p <- do.call(plotly::layout, c(list(p), lay))
  plotly::config(
    plotly::layout(p, hovermode = "x unified",
                   font = list(size = st$font_size, color = st$text),
                   paper_bgcolor = "white", plot_bgcolor = "white",
                   xaxis = list(title = sp1$xlab, tickformat = if (sp1$percent) ".0%" else "~g",
                                gridcolor = st$grid, zeroline = FALSE,
                                showspikes = TRUE, spikemode = "across", spikedash = "dot",
                                spikecolor = "#888888", spikethickness = 1)),
    displaylogo = FALSE)
}

legend_id <- function(i) if (i == 1) "legend" else paste0("legend", i)

plotly_panel <- function(sp, i, legend = NULL) {
  p <- plotly::plot_ly()
  # the legend this panel's curves go to (several panels), or the default one
  add <- function(p, ...) {
    if (is.null(legend)) plotly::add_lines(p, ...) else plotly::add_lines(p, ..., legend = legend)
  }
  hover <- paste0("%{y:", sp$fmt, "}")
  xr <- range(sp$data$x)
  if (!is.null(sp$hline))
    p <- plotly::add_lines(p, x = xr, y = rep(sp$hline, 2), showlegend = FALSE,
                           hoverinfo = "skip", line = list(color = plot_style$ref, dash = "dash"))
  for (l in levels(sp$data$label)) {
    d <- sp$data[sp$data$label == l, ]
    if (!nrow(d)) next
    p <- add(p, x = d$x, y = d$value, name = l, legendgroup = paste0(i, l),
                           line = list(color = unname(sp$colors[l])),
                           hovertemplate = hover)
    if (!is.null(sp$ref)) {
      r <- sp$ref[sp$ref$label == l, ]
      if (nrow(r))
        p <- add(p, x = r$x, y = r$value, name = paste(l, "(EVPI)"),
                               legendgroup = paste0(i, l), showlegend = FALSE,
                               line = list(color = unname(sp$colors[l]), dash = "dash"),
                               hovertemplate = hover)
    }
  }
  yaxis <- list(title = sp$ylab_text, zeroline = FALSE, gridcolor = plot_style$grid,
                tickformat = "~g")                  # plain numbers, as plain_label()
  if (!is.null(sp$ylim)) yaxis$range <- sp$ylim
  if (!sp$percent) yaxis$rangemode <- "tozero"
  plotly::layout(p, yaxis = yaxis,
                 xaxis = list(title = sp$xlab, gridcolor = plot_style$grid, zeroline = FALSE))
}

check_interactive <- function(interactive) {
  if (!is.logical(interactive) || length(interactive) != 1 || is.na(interactive))
    stop("`interactive` must be TRUE or FALSE.", call. = FALSE)
  if (interactive && !requireNamespace("plotly", quietly = TRUE))
    stop("The plotly package is needed for `interactive = TRUE`.", call. = FALSE)
}

# ---- internal ----------------------------------------------------------------

curve_data <- function(x, label, value)
  data.frame(x = x, label = label, value = value)

default_hues <- function(n)
  grDevices::hcl(h = seq(15, 375, length.out = n + 1), l = 65, c = 100)[seq_len(n)]

# Axis numbers written out (0.0006, not 6e-04), as plotly's "~g" format does.
plain_label <- function(x)
  vapply(x, function(v) if (is.na(v)) NA_character_ else
    format(signif(v, 6), scientific = FALSE, trim = TRUE, drop0trailing = TRUE),
    character(1))

percent_label <- function(x) ifelse(is.na(x), NA, paste0(format(100 * x, trim = TRUE), "%"))

# The results to plot, named by their labels, after checking that they can be
# drawn together.
collect_nb <- function(object, more, labels, fun = "plot_dca()",
                       args = "`type` and `labels`") {
  nms <- names(more)
  if (is.null(nms)) nms <- rep("", length(more))
  bad <- !vapply(more, inherits, logical(1), "voipred_nb")
  if (any(bad & nzchar(nms)))
    stop("Unknown argument `", nms[bad & nzchar(nms)][1], "`. The arguments of ",
         fun, " are ", args, "; other arguments must be results of nb().",
         call. = FALSE)
  if (!inherits(object, "voipred_nb") || any(bad))
    stop(if (fun == "plot_dca()") "Only results of nb() can be plotted together."
         else paste0(fun, " needs results of nb()."), call. = FALSE)
  res <- c(list(object), more)
  lab <- vapply(res, function(r) {
    l <- attr(r, "label")
    if (is.null(l)) "Model" else l
  }, character(1))
  lab[-1][nzchar(nms)] <- nms[nzchar(nms)]
  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != length(res) || anyNA(labels))
      stop("`labels` must be ", length(res), " strings, one per result.", call. = FALSE)
    lab <- labels
  }
  if (any(tolower(lab) %in% c("treat all", "treat none")))
    stop("\"Treat All\" and \"Treat None\" cannot be used as labels.", call. = FALSE)
  if (anyDuplicated(lab))
    stop("The results have the same label (\"", lab[anyDuplicated(lab)], "\"); ",
         "give each a different one with `labels` or nb(label = ).", call. = FALSE)
  for (i in seq_along(res)[-1]) check_comparable(res[[1]], res[[i]], lab[1], lab[i])
  names(res) <- lab
  res
}

check_comparable <- function(a, b, la, lb) {
  same_attr <- function(k) identical(attr(a, k), attr(b, k))
  what <- c(
    "outcome type" = !same_attr("outcome_type"),
    "horizon" = !same_attr("horizon"),
    "cause" = !same_attr("cause"),
    "censoring method" = !same_attr("censoring"),
    "estimator (Bayesian or not)" = !identical(isTRUE(attr(a, "bayesian")),
                                               isTRUE(attr(b, "bayesian"))),
    "prior" = isTRUE(attr(a, "bayesian")) && isTRUE(attr(b, "bayesian")) &&
      !same_attr("prior"),
    "number of patients" = !same_attr("n"),
    "thresholds" = !isTRUE(all.equal(a$threshold, b$threshold)))
  if (!any(what))
    what["data (net benefit of treating all)"] <-
      !isTRUE(all.equal(a$NB_all, b$NB_all, tolerance = 1e-10))
  if (any(what))
    stop("\"", la, "\" and \"", lb, "\" cannot be used together: different ",
         paste(names(what)[what], collapse = ", "), ". Compute them with nb() on ",
         "the same data, outcome, horizon, cause, censoring and thresholds.",
         call. = FALSE)
}
