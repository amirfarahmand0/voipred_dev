#' Plot decision curves
#'
#' Plots the results of [nb()]: the decision curve (`type = "nb"`), or the
#' difference in net benefit between the model and the best default strategy
#' (`type = "delta"`). Give several results to draw their curves together.
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
#' `label` (see [nb()]), unless `labels` is given.
#'
#' `autoplot()` returns the plot as a ggplot object, which you can change with
#' the usual `+` (e.g. `+ ggplot2::coord_cartesian(xlim = c(0, 0.5))`).
#' `plot()` draws it and returns it invisibly.
#'
#' @param object,x A result of [nb()].
#' @param y,... More results of [nb()] to draw on the same plot. A named
#'   result is labelled with its name.
#' @param type `"nb"` (decision curve, default) or `"delta"` (\eqn{\Delta NB}).
#' @param labels Optional names of the curves, one per result, in order;
#'   default: the results' labels.
#' @return A ggplot object (invisibly for `plot()`).
#' @examples
#' dev <- gusto[1:5000, ]
#' val <- gusto[5001:10000, ]
#' fit1 <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
#' fit2 <- glm(day30 ~ age, family = binomial, data = dev)
#' z <- seq(0, 0.5, by = 0.01)
#' res1 <- nb(fit1, val, "day30", thresholds = z, label = "Full model")
#' res2 <- nb(fit2, val, "day30", thresholds = z, label = "Age only")
#' plot(res1)
#' plot(res1, res2)
#' plot(res1, res2, type = "delta")
#' @name plot.voipred_nb
NULL

#' @rdname plot.voipred_nb
#' @export
autoplot.voipred_nb <- function(object, ..., type = c("nb", "delta"), labels = NULL) {
  type <- match.arg(type)
  res <- collect_nb(object, list(...), labels)
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
    p <- ggplot(d[!is.na(d$value), ],
                aes(x = .data$threshold, y = .data$value, color = .data$label)) +
      geom_line() +
      scale_color_manual(values = cols) +
      coord_cartesian(ylim = c(-0.1 * prev, prev)) +
      labs(x = "Threshold Probability", y = "Net Benefit", color = "")
  } else {
    d <- do.call(rbind, Map(function(r, l) curve_data(r$threshold, l, r$delta_NB),
                            res, lab))
    d$label <- factor(d$label, levels = lab)
    p <- ggplot(d[!is.na(d$value), ],
                aes(x = .data$threshold, y = .data$value, color = .data$label)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_line() +
      scale_color_manual(values = cols[-(1:2)]) +
      labs(x = "Threshold Probability",
           y = quote(Delta * NB ~ "(" * NB[model] - "max(" * NB[treat ~ all] * "," ~
                       NB[treat ~ none] * "))"),
           color = "")
  }
  p + scale_x_continuous(labels = percent_label) + theme_bw()
}

#' @rdname plot.voipred_nb
#' @export
plot.voipred_nb <- function(x, y, ..., type = c("nb", "delta"), labels = NULL) {
  more <- if (missing(y)) list(...) else c(list(y), list(...))
  p <- do.call(autoplot.voipred_nb, c(list(x), more, list(type = type, labels = labels)))
  print(p)
  invisible(p)
}

# ---- internal ----------------------------------------------------------------

curve_data <- function(threshold, label, value)
  data.frame(threshold = threshold, label = label, value = value)

default_hues <- function(n)
  grDevices::hcl(h = seq(15, 375, length.out = n + 1), l = 65, c = 100)[seq_len(n)]

percent_label <- function(x) ifelse(is.na(x), NA, paste0(format(100 * x, trim = TRUE), "%"))

# The results to plot, named by their labels, after checking that they can be
# drawn together.
collect_nb <- function(object, more, labels) {
  nms <- names(more)
  if (is.null(nms)) nms <- rep("", length(more))
  bad <- !vapply(more, inherits, logical(1), "voipred_nb")
  if (any(bad & nzchar(nms)))
    stop("Unknown argument `", nms[bad & nzchar(nms)][1], "`. The arguments are ",
         "`type` and `labels`; other arguments must be results of nb().", call. = FALSE)
  if (!inherits(object, "voipred_nb") || any(bad))
    stop("Only results of nb() can be plotted together.", call. = FALSE)
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
    "number of patients" = !same_attr("n"),
    "thresholds" = !isTRUE(all.equal(a$threshold, b$threshold)))
  if (!any(what))
    what["data (net benefit of treating all)"] <-
      !isTRUE(all.equal(a$NB_all, b$NB_all, tolerance = 1e-10))
  if (any(what))
    stop("\"", la, "\" and \"", lb, "\" cannot be plotted together: different ",
         paste(names(what)[what], collapse = ", "), ". Compute them with nb() on ",
         "the same data, outcome, horizon, cause, censoring and thresholds.",
         call. = FALSE)
}
