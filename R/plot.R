#' Plot decision curves
#'
#' The plot of results of [nb_val()]: the decision curve (`type = "nb"`), or the
#' difference in net benefit between each model and the best default strategy
#' (`type = "delta"`). A result of [nb_val()] with several models draws all of
#' them; several results are drawn together, combined as by
#' [c()][c.voipred_nb].
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
#'   NB_{none})} (column `delta_NB` of the estimates) of each model, with a
#'   dashed line at 0. Above 0 the model is better than both default
#'   strategies.
#'
#' Results drawn together must come from the same data, outcome, horizon,
#' cause, censoring method, prior and thresholds (voipred checks this).
#' Thresholds where the net benefit is `NA` are left out. Each curve is named
#' by its `label` (see [nb_val()]), unless `labels` is given; a model has the same
#' colour in both types of plot.
#'
#' The panel below shows EVPI (from `evpi`) and EVSI for each future sample
#' size (from `evsi`) against the threshold; without `evpi`, the EVPI of the
#' `evsi` result is drawn. (EVSI against the threshold is a voipred addition:
#' Sadatsafavi et al. (2025) plot EVSI against the future sample size, as
#' [plot.voipred_evsi()] does.)
#'
#' The static plot is a ggplot object (a patchwork object with the
#' value-of-information panel): it is drawn when printed, and can be changed
#' with the usual `+` (e.g. `+ ggplot2::coord_cartesian(xlim = c(0, 0.5))`) or
#' saved with `ggplot2::ggsave()`. With `interactive = TRUE` it is a plotly
#' widget: hovering shows the values of every curve at a threshold, the panels
#' zoom together, and clicking a legend entry hides or shows a curve; save it
#' with `htmlwidgets::saveWidget()`.
#'
#' @param x A result of [nb_val()], for one model or several.
#' @param ... More results of [nb_val()] to draw on the same plot. A result with
#'   one model given here with a name is labelled with its name.
#' @param type `"nb"` (decision curve, default) or `"delta"` (\eqn{\Delta NB}).
#' @param labels Optional names of the curves, one per model, in order;
#'   default: the models' labels.
#' @param evpi,evsi Optional results of [evpi()] and [evsi()], drawn in a panel
#'   below the decision curve: EVPI, and EVSI for each future sample size,
#'   against the threshold. The static version needs the patchwork package.
#' @param bands If the results have credible intervals (`nb_val(interval = )`),
#'   `TRUE` (default) draws them as shaded bands around each model's curve
#'   (and Treat All in the decision curve), in the curves' colours; `FALSE`
#'   draws the curves only. Interactive plots show the interval when
#'   hovering.
#' @param colors Optional colours of the curves, instead of the defaults:
#'   named by label, e.g. `c(Full = "darkred", "Treat All" = "grey50")` (any
#'   curve, including Treat All, Treat None and, with `evpi`/`evsi`, the EVPI
#'   and EVSI curves; the others keep their colours), or unnamed, one per
#'   model in order. Colour names or hex codes. By default, a model has the
#'   same colour in every plot of voipred.
#' @param per Scale of the net benefit (and of EVPI and EVSI in the panel
#'   below): the number of patients it is expressed per. Default 1 (per
#'   patient); e.g. `per = 100000` draws the net number of true positives per
#'   100,000 patients, and the y axes say "(per 100,000 patients)". The
#'   default can be set with `options(voipred.per = )`.
#' @param interactive `TRUE` for an interactive plot (needs the plotly
#'   package). Default `FALSE`.
#' @return A ggplot object, a patchwork object (with `evpi` or `evsi`), or a
#'   plotly widget (`interactive = TRUE`).
#' @seealso [nb_val()], [plot.voipred_evpi()], [plot.voipred_evsi()],
#'   [plot.voipred_p_useful()], [plot.voipred_p_best()]
#' @examples
#' df_dev_binary <- gusto[1:5000, ]
#' df_val_binary <- gusto[5001:10000, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse,
#'             family = binomial, data = df_dev_binary)
#' df_val_binary$p1 <- predict(fit1, newdata = df_val_binary, type = "response")
#' fit2 <- glm(day30 ~ age, family = binomial, data = df_dev_binary)
#' df_val_binary$p2 <- predict(fit2, newdata = df_val_binary, type = "response")
#' z <- seq(0, 0.5, by = 0.01)
#' res1 <- nb_val(data = df_val_binary,
#'                predicted_risk = df_val_binary$p1,
#'                outcome = "day30",
#'                thresholds = z,
#'                label = "Full model")
#' res2 <- nb_val(data = df_val_binary,
#'                predicted_risk = df_val_binary$p2,
#'                outcome = "day30",
#'                thresholds = z,
#'                label = "Age only")
#' plot(res1)
#' plot(res1, res2)
#' plot(res1, res2, type = "delta")
#' \donttest{
#' set.seed(1)
#' e <- evpi(res1, n_sim = 200)
#' if (requireNamespace("patchwork", quietly = TRUE)) plot(res1, evpi = e)
#' if (requireNamespace("plotly", quietly = TRUE)) plot(res1, evpi = e, interactive = TRUE)
#' }
#' @export
plot.voipred_nb <- function(x, ..., type = c("nb", "delta"), labels = NULL, evpi = NULL,
                            evsi = NULL, colors = NULL, per = getOption("voipred.per", 1),
                            bands = TRUE, interactive = FALSE) {
  type <- match.arg(type)
  check_interactive(interactive)
  per <- check_per(per)
  if (!is.logical(bands) || length(bands) != 1 || is.na(bands))
    stop("`bands` must be TRUE or FALSE.", call. = FALSE)
  res <- collect_nb(x, list(...), labels, fun = "plot()",
                    args = "`type`, `labels`, `evpi`, `evsi`, `colors`, `per`, `bands` and `interactive`")
  specs <- list(scale_spec(dca_spec(res, type, bands), per,
                           if (type == "nb") "Net Benefit" else "\u0394NB"))
  if (!is.null(evpi) || !is.null(evsi)) {
    bottom <- voi_panel_spec(evpi, evsi)
    specs[[2]] <- scale_spec(bottom, per, bottom$ylab_text)
  }
  specs <- recolor(specs, colors, res$models, "model")
  render_plot(specs, interactive, heights = if (length(specs) > 1) c(0.62, 0.38))
}

# ---- plot specifications -----------------------------------------------------
# A plot is built from one or more "specs": the curves (columns x, value,
# label), dashed reference curves (`ref`), axis labels (`ylab` for ggplot2,
# `ylab_text` for plotly), colours, y limits, an optional dashed line at
# `hline`, and the number format of the hover text. render_plot() draws them
# with ggplot2 (stacked with patchwork when there are several) or plotly.

new_spec <- function(data, colors, xlab, ylab, ylab_text = ylab, ylim = NULL,
                     hline = NULL, ref = NULL, percent = TRUE, fmt = ".4f", band = NULL) {
  data <- data[is.finite(data$value), , drop = FALSE]
  if (!is.null(ref)) ref <- ref[is.finite(ref$value), , drop = FALSE]
  if (!is.null(band)) {
    band <- band[is.finite(band$lower) & is.finite(band$upper), , drop = FALSE]
    band$label <- factor(band$label, levels = levels(data$label))
  }
  list(data = data, colors = colors, xlab = xlab, ylab = ylab, ylab_text = ylab_text,
       ylim = ylim, hline = hline, ref = ref, percent = percent, fmt = fmt, band = band)
}

dca_spec <- function(x, type, bands = TRUE) {
  e <- x$estimates
  lab <- x$models
  z <- nb_thresholds(x)
  # credible intervals (if computed): a band around each model's curve, and
  # around Treat All in the decision curve
  band <- NULL
  if (bands && !is.null(x$interval)) {
    k <- if (type == "nb") "NB_model" else "delta_NB"
    band <- do.call(rbind, lapply(lab, function(l) {
      d <- e[e$model == l, ]
      data.frame(x = d$threshold, label = l, lower = d[[paste0(k, "_lower")]],
                 upper = d[[paste0(k, "_upper")]])
    }))
    if (type == "nb") {
      d <- e[e$model == lab[1], ]
      band <- rbind(data.frame(x = d$threshold, label = "Treat All", lower = d$NB_all_lower,
                               upper = d$NB_all_upper), band)
    }
  }
  # ggplot2's default colours for Treat All, Treat None and the models, so a
  # model has the same colour in both types of plot
  cols <- default_hues(length(lab) + 2)
  models <- do.call(rbind, lapply(lab, function(l) {
    d <- e[e$model == l, ]
    curve_data(d$threshold, l, if (type == "nb") d$NB_model else d$delta_NB)
  }))
  if (type == "nb") {
    d <- rbind(curve_data(z, "Treat All", nb_all(x)),
               curve_data(z, "Treat None", rep(0, length(z))),
               models)
    d$label <- factor(d$label, levels = c("Treat All", "Treat None", lab))
    prev <- x$outcome$prevalence
    new_spec(d, stats::setNames(cols, levels(d$label)), "Threshold Probability",
             "Net Benefit", ylim = c(-0.1 * prev, prev), band = band)
  } else {
    d <- models
    d$label <- factor(d$label, levels = lab)
    new_spec(d, stats::setNames(cols[-(1:2)], lab), "Threshold Probability",
             quote(Delta * NB ~ "(" * NB[model] - "max(" * NB[treat ~ all] * "," ~
                     NB[treat ~ none] * "))"),
             ylab_text = "&#916;NB (NB<sub>model</sub> &#8722; max(NB<sub>treat all</sub>, NB<sub>treat none</sub>))",
             hline = 0, band = band)
  }
}

# A plot spec with its values (curves, reference lines, limits) in net
# benefit per `per` patients, and the scale in the axis label `label`.
scale_spec <- function(sp, per, label) {
  if (per == 1) return(sp)
  sp$data$value <- sp$data$value * per
  if (!is.null(sp$ref)) sp$ref$value <- sp$ref$value * per
  if (!is.null(sp$ylim)) sp$ylim <- sp$ylim * per
  if (!is.null(sp$hline)) sp$hline <- sp$hline * per
  if (!is.null(sp$band)) {
    sp$band$lower <- sp$band$lower * per
    sp$band$upper <- sp$band$upper * per
  }
  sp$ylab <- sp$ylab_text <- paste0(label, " (", per_text(per), ")")
  sp$fmt <- ".4~r"                                  # 4 significant digits in the hover
  sp
}

# EVPI and EVSI (for each future sample size) against the threshold.
# The y-axis says what is drawn: "EVPI" (EVPI only) or "EVPI / EVSI" (both).
voi_panel_spec <- function(evpi, evsi) {
  if (!is.null(evpi) && !inherits(evpi, "voipred_evpi"))
    stop("`evpi` must be a result of evpi().", call. = FALSE)
  if (!is.null(evsi) && !inherits(evsi, "voipred_evsi"))
    stop("`evsi` must be a result of evsi().", call. = FALSE)
  d <- NULL
  if (!is.null(evpi)) d <- curve_data(evpi$estimates$threshold, "EVPI", evpi$estimates$EVPI)
  else {
    o <- evsi$estimates
    first <- o[o$future_n == o$future_n[1], ]
    d <- curve_data(first$threshold, "EVPI", first$EVPI)
  }
  lev <- "EVPI"
  if (!is.null(evsi)) {
    o <- evsi$estimates
    ns <- unique(o$future_n)
    lev <- c(lev, paste0("EVSI (N = ", format(ns, big.mark = ",", trim = TRUE), ")"))
    for (i in seq_along(ns)) {
      e <- o[o$future_n == ns[i], ]
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
  if (!has_package("patchwork"))
    stop("Drawing the value-of-information panel under the decision curve needs the ",
         "patchwork package. Install it with install.packages(\"patchwork\"), or use ",
         "`interactive = TRUE`.", call. = FALSE)
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
  # credible intervals: shaded bands in the colours of the curves
  if (!is.null(sp$band) && nrow(sp$band))
    p <- p + geom_ribbon(data = sp$band, aes(x = .data$x, ymin = .data$lower,
                                             ymax = .data$upper, fill = .data$label),
                         inherit.aes = FALSE, alpha = 0.15, show.legend = FALSE) +
      scale_fill_manual(values = sp$colors, guide = "none")
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

# A bound of a credible interval in the hover: 4 significant digits.
fmt_hover <- function(v) trimws(formatC(v, format = "fg", digits = 4))

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
    b <- if (!is.null(sp$band)) sp$band[sp$band$label == l, ]
    tmpl <- hover
    txt <- NULL
    if (!is.null(b) && nrow(b)) {
      # the band, behind the curve, and the interval in the hover
      rgb <- grDevices::col2rgb(sp$colors[l])
      fill <- sprintf("rgba(%d, %d, %d, 0.15)", rgb[1], rgb[2], rgb[3])
      p <- plotly::add_ribbons(p, x = b$x, ymin = b$lower, ymax = b$upper,
                               fillcolor = fill, line = list(width = 0),
                               legendgroup = paste0(i, l), showlegend = FALSE,
                               hoverinfo = "skip", name = paste(l, "interval"))
      j <- match(d$x, b$x)
      txt <- ifelse(is.na(j), "", sprintf(" (%s to %s)", fmt_hover(b$lower[j]),
                                          fmt_hover(b$upper[j])))
      tmpl <- paste0(hover, "%{text}")
    }
    p <- add(p, x = d$x, y = d$value, name = l, legendgroup = paste0(i, l),
                           line = list(color = unname(sp$colors[l])),
                           text = txt, hovertemplate = tmpl)
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
  if (interactive && !has_package("plotly"))
    stop("Interactive plots (`interactive = TRUE`) need the plotly package. Install it ",
         "with install.packages(\"plotly\").", call. = FALSE)
}

# Whether an optional package is installed (a function, so tests can mock it).
has_package <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# The plot specs with the colours of `colors` instead of the defaults: named
# by curve label (any curve of any panel), or unnamed for the curves
# `ordered` (e.g. the models), in order. The other curves keep their colours.
recolor <- function(specs, colors, ordered, what = "model") {
  if (is.null(colors)) return(specs)
  defaults <- unlist(lapply(specs, function(sp) sp$colors))
  new <- plot_colors(defaults, colors, ordered, what)
  lapply(specs, function(sp) {
    sp$colors[] <- new[names(sp$colors)]
    sp
  })
}

plot_colors <- function(cols, colors, ordered, what = "model") {
  if (!is.character(colors) || !length(colors) || anyNA(colors))
    stop("`colors` must be a character vector of colours, named by label, e.g. ",
         "colors = c(Full = \"darkred\"), or in order, e.g. ",
         "colors = c(\"darkred\", \"steelblue\").", call. = FALSE)
  ok <- vapply(colors, function(k)
    !inherits(tryCatch(grDevices::col2rgb(k), error = function(e) e), "error"), logical(1))
  if (!all(ok))
    stop("`colors` has an invalid colour: \"", colors[!ok][1], "\". Use colour names ",
         "(e.g. \"darkred\") or hex codes (e.g. \"#8B0000\").", call. = FALSE)
  nms <- names(colors)
  if (is.null(nms) || !any(nzchar(nms))) {
    if (length(colors) != length(ordered))
      stop("`colors` without names must have ", length(ordered), " colour",
           if (length(ordered) > 1) "s", ", one per ", what, " in order (",
           paste(ordered, collapse = ", "), "); or name them by label, e.g. ",
           "colors = c(", r_name(ordered[1]), " = \"darkred\").", call. = FALSE)
    cols[ordered] <- unname(colors)
    return(cols)
  }
  if (!all(nzchar(nms)))
    stop("Name every colour in `colors`, or none.", call. = FALSE)
  unknown <- setdiff(nms, names(cols))
  if (length(unknown))
    stop("`colors` names \"", unknown[1], "\", which is not a curve of this plot: ",
         paste0("\"", names(cols), "\"", collapse = ", "), ".", call. = FALSE)
  cols[nms] <- unname(colors)
  cols
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

# The results of nb_val() given to plot() or a value-of-information function,
# checked and merged into one result with all their models.
collect_nb <- function(object, more, labels, fun = "plot()",
                       args = "`type` and `labels`") {
  nms <- names(more)
  if (is.null(nms)) nms <- rep("", length(more))
  bad <- !vapply(more, inherits, logical(1), "voipred_nb")
  if (any(bad & nzchar(nms)))
    stop("Unknown argument `", nms[bad & nzchar(nms)][1], "`. The arguments of ",
         fun, " are ", args, "; other arguments must be results of nb_val().",
         call. = FALSE)
  if (!inherits(object, "voipred_nb") || any(bad))
    stop(if (fun == "plot()") "Only results of nb_val() can be plotted together."
         else paste0(fun, " needs results of nb_val()."), call. = FALSE)
  res <- c(list(object), more)
  names(res) <- c("", nms)
  merge_nb(res, labels)
}

#' Combine results of nb_val()
#'
#' Combines results of [nb_val()] for different models, computed on the same
#' data, outcome, thresholds and prior (and horizon, cause and censoring),
#' into one result with all their models, as if they had been given together
#' to [nb_val()]. A result given with a name (e.g. `c(Full = res1, res2)`) is
#' labelled with it (results with one model only). `plot()` and the
#' value-of-information functions combine the results given to them in the
#' same way.
#'
#' @param ... Results of [nb_val()].
#' @return A result of [nb_val()] with the models of all the results, in
#'   order.
#' @seealso [nb_val()]
#' @examples
#' df_val_binary <- gusto[5001:5500, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial,
#'             data = gusto[1:5000, ])
#' fit2 <- glm(day30 ~ age, family = binomial, data = gusto[1:5000, ])
#' z <- seq(0.01, 0.3, by = 0.01)
#' res1 <- nb_val(df_val_binary, fit1, "day30", thresholds = z)
#' res2 <- nb_val(df_val_binary, fit2, "day30", thresholds = z)
#' both <- c(Full = res1, `Age only` = res2)
#' both
#' @export
c.voipred_nb <- function(...) {
  res <- list(...)
  if (!all(vapply(res, inherits, logical(1), "voipred_nb")))
    stop("c() combines results of nb_val() only.", call. = FALSE)
  merge_nb(res)
}

# One result with the models of the results `res` (a list, whose names, if
# any, relabel results with one model), relabelled with `labels` if given.
merge_nb <- function(res, labels = NULL) {
  nms <- names(res)
  if (is.null(nms)) nms <- rep("", length(res))
  groups <- lapply(seq_along(res), function(i) {
    r <- res[[i]]
    if (!nzchar(nms[i])) return(r$models)
    if (length(r$models) > 1)
      stop("A name relabels a result with one model, but the result named \"", nms[i],
           "\" has ", length(r$models), " (", paste(r$models, collapse = ", "),
           "). Relabel its models with `labels` or nb_val(label = ).", call. = FALSE)
    nms[i]
  })
  sizes <- lengths(groups)
  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != sum(sizes) || anyNA(labels))
      labels_error(sum(sizes), "model")
    groups <- split(unname(labels), rep(seq_along(res), sizes))
  }
  lab <- unlist(groups, use.names = FALSE)
  check_model_labels(lab, "give each a different one with `labels` or nb_val(label = )")
  for (i in seq_along(res)[-1])
    check_comparable(res[[1]], res[[i]], groups[[1]], groups[[i]])
  # credible intervals are kept only if every result has the same ones
  ints <- lapply(res, function(r) r$interval)
  if (!all(vapply(ints, identical, logical(1), ints[[1]]))) {
    warning("The results do not have the same credible intervals (`interval` and ",
            "`n_sim` of nb_val()), so the combined result has none.", call. = FALSE)
    res <- lapply(res, function(r) {
      r$estimates <- r$estimates[, c("threshold", "model", "NB_model", "NB_all", "NB_none",
                                     "delta_NB")]
      r$interval <- NULL
      r
    })
  }
  relabel <- function(df, r, new) {
    if (is.null(df)) return(NULL)
    df$model <- new[match(df$model, r$models)]
    df
  }
  out <- res[[1]]
  out$models <- lab
  out$estimates <- do.call(rbind, Map(function(r, g) relabel(r$estimates, r, g), res, groups))
  rownames(out$estimates) <- NULL
  if (!is.null(out$counts)) {
    out$counts <- do.call(rbind, Map(function(r, g) relabel(r$counts, r, g), res, groups))
    rownames(out$counts) <- NULL
  }
  out$risk <- do.call(cbind, lapply(res, function(r) r$risk))
  colnames(out$risk) <- lab
  out
}

check_comparable <- function(a, b, la, lb) {
  same <- function(k) identical(a$outcome[[k]], b$outcome[[k]])
  what <- c(
    "outcome type" = !same("type"),
    "horizon" = !same("horizon"),
    "cause" = !same("cause"),
    "censoring method" = !same("censoring"),
    "prior" = !identical(a$prior, b$prior),
    "case weights" = !identical(a$weights, b$weights),
    "number of patients" = !same("n"),
    "thresholds" = !isTRUE(all.equal(nb_thresholds(a), nb_thresholds(b))))
  if (!any(what)) {
    what["data (net benefit of treating all)"] <-
      !isTRUE(all.equal(nb_all(a), nb_all(b), tolerance = 1e-10))
    if (identical(a$outcome$censoring, "cipcw"))
      what["censoring weights"] <-
        !isTRUE(all.equal(a$cipcw[c("weight_event", "weight_horizon")],
                          b$cipcw[c("weight_event", "weight_horizon")]))
  }
  if (any(what))
    stop("\"", paste(la, collapse = ", "), "\" and \"", paste(lb, collapse = ", "),
         "\" cannot be used together: different ",
         paste(names(what)[what], collapse = ", "), ". Compute them with nb_val() on ",
         "the same data, outcome, horizon, cause, censoring, prior and thresholds, ",
         "or together in one call of nb_val().", call. = FALSE)
}
