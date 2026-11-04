# Subgroup analysis of net benefit under subgroup utility heterogeneity, as in
# Wang, Wang, Van Calster, van Geloven, Steyerberg and Wynants (Population-level
# decision curve analysis may mislead the evaluation of prediction model
# usefulness under subgroup utility heterogeneity): decision curve analysis
# within subgroups (Step 1), whether the subgroup-specific differences in net
# benefit against a comparator (diff_NB; the Delta NB_i of the paper) have the same
# sign (Step 2), and the inconsistency conditions (Step 3): the inconsistency
# region and the K-value for two subgroups, and the utility-weighted sum for
# any number of subgroups.

# ---- the subgroups (nb_val(subgroup = )) -------------------------------------

# The subgroup of each patient: a factor (length n) and the name of the
# variable. `subgroup` is a column name of `data`, several column names (their
# combinations), or a factor or vector with one value per row of `data`.
parse_subgroup <- function(subgroup, data, expr) {
  n <- nrow(data)
  if (is.character(subgroup) && length(subgroup) >= 1 && length(subgroup) < n &&
      all(subgroup %in% names(data))) {
    name <- paste(subgroup, collapse = ", ")
    vals <- lapply(subgroup, function(v) check_subgroup_values(data[[v]], v))
    g <- if (length(vals) == 1) vals[[1]]
         else interaction(vals, sep = ", ", lex.order = TRUE)
  } else if (is.character(subgroup) && length(subgroup) != n &&
             (length(subgroup) <= 5 || any(subgroup %in% names(data)))) {
    stop("`subgroup` names columns that are not in `data`: ",
         paste0("\"", setdiff(subgroup, names(data)), "\"", collapse = ", "), ".",
         call. = FALSE)
  } else {
    if (length(subgroup) != n)
      stop("`subgroup` has ", length(subgroup), " values but `data` has ", n, " rows: ",
           "give one subgroup per patient, or the name of a column of `data`.",
           call. = FALSE)
    name <- if (is.name(expr)) as.character(expr)
            else if (is.call(expr) && identical(expr[[1]], as.name("$")))
              as.character(expr[[3]]) else "subgroup"
    g <- check_subgroup_values(subgroup, name)
  }
  used <- levels(droplevels(g))
  if (length(used) < length(levels(g)))
    message("Subgroups with no patients are left out: ",
            paste(setdiff(levels(g), used), collapse = ", "), ".")
  g <- droplevels(g)
  if (nlevels(g) < 2)
    stop("`subgroup` has only one subgroup (", levels(g), "): it needs at least two.",
         call. = FALSE)
  list(group = g, variable = name)
}

check_subgroup_values <- function(x, name) {
  if (anyNA(x))
    stop("`subgroup` (", name, ") is missing for ", sum(is.na(x)), " patients: every ",
         "patient needs a subgroup.", call. = FALSE)
  if (is.factor(x)) return(x)
  if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
  if (is.numeric(x) && length(unique(x)) > 10)
    stop("`subgroup` (", name, ") is numeric with ", length(unique(x)), " different ",
         "values: give categories instead, e.g. with cut().", call. = FALSE)
  if (!is.numeric(x) && !is.character(x))
    stop("`subgroup` (", name, ") must be a factor, character, logical or a few ",
         "numeric values.", call. = FALSE)
  factor(x, levels = sort(unique(x)))
}

# The net benefit of each model, treat all and treat none within each subgroup,
# with the estimator of nb_val() (with cIPCW, the censoring weights of the whole
# sample): a data frame with one row per threshold, model and subgroup.
subgroup_nb <- function(g, info, P, thresholds, prior, w, censoring, k, horizon, cws,
                        labels) {
  out <- list()
  for (lev in levels(g)) {
    i <- which(g == lev)
    wi <- if (is.null(w)) NULL else w[i]
    size <- if (is.null(w)) length(i) else sum(w[i])
    for (m in seq_along(labels)) {
      r <- if (info$type == "binary") {
        bayes_nb_mean(info$y[i], P[i, m], thresholds, prior, wi)
      } else if (censoring == "km") {
        suppressWarnings(nb_core_km(info$time[i], info$status[i], k, P[i, m], thresholds,
                                    horizon, wi))
      } else {
        nb_core_cipcw(info$time[i], info$status[i], k, P[i, m], thresholds, horizon,
                      cws$weight_event[i], cws$weight_horizon[i], wi)
      }
      out[[length(out) + 1]] <- data.frame(threshold = thresholds, model = labels[m],
                                           subgroup = lev, n = size,
                                           NB_model = r$NB_model, NB_all = r$NB_all,
                                           NB_none = 0)
    }
  }
  d <- do.call(rbind, out)
  d$subgroup <- factor(d$subgroup, levels = levels(g))
  rownames(d) <- NULL
  d
}

# ---- the checks of the framework (summary(), plot()) -------------------------

# The comparator of the subgroup analysis: "Treat All" (default), "Treat None",
# "best default" (treat all or treat none, whichever has the higher net
# benefit in the whole sample at the threshold, the same for every subgroup),
# or the label of another model.
check_comparator <- function(comparator, models) {
  ok <- c("Treat All", "Treat None", "best default", models)
  if (!is.character(comparator) || length(comparator) != 1 || !comparator %in% ok)
    stop("`comparator` must be one of ", paste0("\"", ok, "\"", collapse = ", "), ".",
         call. = FALSE)
  comparator
}

# The subgroup shares: from the sample (default), or given, named by subgroup.
subgroup_shares <- function(sg, shares) {
  lev <- levels(sg$estimates$subgroup)
  n <- tapply(sg$estimates$n, sg$estimates$subgroup, `[`, 1)[lev]
  if (is.null(shares)) return(n / sum(n))
  if (!is.numeric(shares) || is.null(names(shares)) || anyNA(shares) || any(shares < 0) ||
      !setequal(names(shares), lev) || sum(shares) <= 0)
    stop("`shares` must be non-negative numbers named by the subgroups (",
         paste(lev, collapse = ", "), ").", call. = FALSE)
  s <- shares[lev]
  s / sum(s)
}

# `k`, the incremental utility of a true positive relative to a false negative
# (a - c) in each subgroup, e.g. the relative risk reduction of the treatment:
# values (a named vector) or ranges (a named list of c(lower, upper)).
check_k <- function(k, lev) {
  if (is.null(k)) return(NULL)
  bad <- function()
    stop("`k` must be positive values named by the subgroups (", paste(lev, collapse = ", "),
         "), e.g. k = c(", lev[1], " = 0.28, ", lev[2], " = 0.15), or ranges, e.g. ",
         "k = list(", lev[1], " = c(0.2, 0.3), ...).", call. = FALSE)
  if (is.null(names(k)) || !setequal(names(k), lev)) bad()
  if (is.list(k)) {
    if (!all(vapply(k, function(v) is.numeric(v) && length(v) == 2 && !anyNA(v) &&
                                   all(v > 0) && v[1] <= v[2], logical(1)))) bad()
    r <- do.call(rbind, k[lev])
    return(list(lower = r[, 1], upper = r[, 2], range = TRUE))
  }
  if (!is.numeric(k) || anyNA(k) || any(k <= 0)) bad()
  list(lower = k[lev], upper = k[lev], range = FALSE)
}

# The checks of the framework at each threshold and for each model: the
# subgroup-specific differences in net benefit against the comparator (diff_NB),
# the population difference (their share-weighted sum), whether they have the same sign (Step 2), the
# K-value (two subgroups with opposite signs) and, with k, whether population
# NB and population utility favour different strategies (Step 3).
subgroup_check <- function(x, comparator = "Treat All", k = NULL, shares = NULL) {
  sg <- x$subgroups
  comparator <- check_comparator(comparator, x$models)
  lev <- levels(sg$estimates$subgroup)
  G <- length(lev)
  s <- subgroup_shares(sg, shares)
  kk <- check_k(k, lev)
  e <- sg$estimates
  z <- nb_thresholds(x)
  models <- setdiff(x$models, if (comparator %in% x$models) comparator)
  pop_all <- nb_all(x)
  rows <- list()
  for (m in models) {
    d <- matrix(NA_real_, length(z), G, dimnames = list(NULL, lev))
    for (j in seq_len(G)) {
      ej <- e[e$model == m & e$subgroup == lev[j], ]
      comp <- switch(comparator,
                     "Treat All" = ej$NB_all,
                     "Treat None" = rep(0, length(z)),
                     "best default" = ifelse(pop_all > 0, ej$NB_all, 0),
                     e$NB_model[e$model == comparator & e$subgroup == lev[j]])
      d[, j] <- ej$NB_model - comp
    }
    pop <- as.vector(d %*% s)
    conflict <- apply(d, 1, function(v) any(v > 1e-12, na.rm = TRUE) &&
                                          any(v < -1e-12, na.rm = TRUE))
    out <- data.frame(threshold = z, model = m, check.names = FALSE)
    for (j in seq_len(G)) out[[paste0("diff_NB_", lev[j])]] <- d[, j]
    out$diff_NB_population <- pop
    out$signs <- ifelse(apply(d, 1, anyNA), NA, ifelse(conflict, "opposite", "same"))
    if (G == 2) {
      K <- -(s[2] * d[, 2]) / (s[1] * d[, 1])
      out$K <- ifelse(conflict, K, NA_real_)
    }
    if (!is.null(kk)) {
      if (!kk$range) {
        u <- as.vector(d %*% (s * kk$lower))
        out$diff_utility <- u
        out$consistent <- ifelse(is.na(u) | is.na(pop), NA,
                                 !(sign(u) * sign(pop) < 0))
      } else {
        lo <- as.vector(pmin(d * rep(kk$lower, each = nrow(d)),
                             d * rep(kk$upper, each = nrow(d))) %*% s)
        hi <- as.vector(pmax(d * rep(kk$lower, each = nrow(d)),
                             d * rep(kk$upper, each = nrow(d))) %*% s)
        out$diff_utility_lower <- lo
        out$diff_utility_upper <- hi
        out$consistent <- ifelse(is.na(pop), NA,
                                 !((pop > 0 & lo < 0) | (pop < 0 & hi > 0)))
      }
    }
    rows[[m]] <- out
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  list(table = tab, comparator = comparator, shares = s, k = kk, levels = lev,
       variable = sg$variable)
}

# The lines of summary() about the subgroups.
subgroup_lines <- function(chk) {
  s <- chk$shares
  lev <- chk$levels
  lines <- paste0("Subgroups (", chk$variable, "): ",
                  paste0(lev, " (", round(100 * s), "%)", collapse = ", "),
                  "; difference in net benefit (diff_NB) against ", chk$comparator)
  tab <- chk$table
  for (m in unique(tab$model)) {
    t <- tab[tab$model == m, ]
    o <- order(t$threshold)
    opp <- !is.na(t$signs) & t$signs == "opposite"
    who <- if (length(unique(tab$model)) > 1) paste0(m, ": ") else ""
    if (!any(opp)) {
      lines <- c(lines, paste0(who, "the subgroup differences in net benefit have the same sign at every ",
                               "threshold, so differences in a - c between subgroups ",
                               "cannot change the population conclusion"))
      next
    }
    lines <- c(lines, paste0(who, "opposite signs at ", segments_of(t$threshold[o], opp[o])))
    if (length(lev) == 2) {
      K <- t$K[opp]
      lines <- c(lines, paste0(who, "K-value (k_", lev[1], " / k_", lev[2], "): ",
                               if (length(K) == 1) fmt_num(K)
                               else paste(fmt_num(min(K)), "to", fmt_num(max(K)))))
    }
    if (!is.null(chk$k)) {
      inc <- !is.na(t$consistent) & !t$consistent
      lines <- c(lines, paste0(who, if (!any(inc))
        "with k, population NB and population utility favour the same strategy at every threshold"
        else paste0("with k, population NB and population utility ",
                    if (chk$k$range) "may favour" else "favour",
                    " different strategies at ", segments_of(t$threshold[o], inc[o]))))
    }
  }
  lines
}

#' The K-value from subgroup-specific net benefits
#'
#' The K-value of Wang et al. for two subgroups, from their differences in net
#' benefit (e.g. published values), without patient data: the critical ratio
#' \eqn{k_1/k_2} of the incremental utility of a true positive relative to a
#' false negative (\eqn{a - c}) at which population-level net benefit and
#' population utility start to favour different strategies,
#' \deqn{K = -\frac{n_2 \Delta NB_2}{n_1 \Delta NB_1}.}
#' It is defined when the two \eqn{\Delta NB} have opposite signs. When
#' \eqn{k_1 > k_2}, population net benefit and utility favour different
#' strategies if \eqn{k_1/k_2 > K}; when \eqn{k_1 < k_2}, if \eqn{k_1/k_2 < K}.
#' With patient data, `summary()` of a result of [nb_val()] with `subgroup`
#' gives the same check at every threshold.
#'
#' @param diff_nb The difference in net benefit between the model and the
#'   comparator (the same strategy in both subgroups) in each of the two subgroups, named by subgroup, e.g.
#'   `c(Men = -0.0016, Women = 0.00345)`.
#' @param shares Optional relative sizes of the two subgroups, named like `diff_nb`.
#'   Default: equal.
#' @param k Optional \eqn{a - c} of each subgroup (e.g. the relative risk
#'   reductions of the treatment), named like `diff_nb`.
#' @return An object of class `voipred_k_value`, a list with `diff_nb`, `shares`,
#'   `K` (NA when the signs are the same), `k` and `consistent` (with `k`:
#'   whether population net benefit and population utility favour the same
#'   strategy), printed as a short report.
#' @references Wang J, Wang KZ, Van Calster B, van Geloven N, Steyerberg E,
#'   Wynants L. Population-level decision curve analysis may mislead the
#'   evaluation of prediction model usefulness under subgroup utility
#'   heterogeneity.
#' @examples
#' # SCORE2 at a 3% threshold (the case study of Wang et al.)
#' k_value(c(Men = -0.0016, Women = 0.00345))
#' k_value(c(Men = -0.0016, Women = 0.00345), k = c(Men = 0.28, Women = 0.15))
#' k_value(c(Men = -0.0016, Women = 0.00345), k = c(Men = 0.28, Women = 0.12))
#' @export
k_value <- function(diff_nb, shares = NULL, k = NULL) {
  if (!is.numeric(diff_nb) || length(diff_nb) != 2 || anyNA(diff_nb) || is.null(names(diff_nb)) ||
      any(!nzchar(names(diff_nb))))
    stop("`diff_nb` must be the differences in net benefit of two subgroups, named by subgroup, e.g. ",
         "c(Men = -0.0016, Women = 0.00345).", call. = FALSE)
  lev <- names(diff_nb)
  s <- if (is.null(shares)) c(0.5, 0.5) else {
    if (!is.numeric(shares) || !setequal(names(shares), lev) || any(shares <= 0))
      stop("`shares` must be two positive numbers named like `diff_nb`.", call. = FALSE)
    shares[lev] / sum(shares)
  }
  names(s) <- lev
  opposite <- sign(diff_nb[1]) * sign(diff_nb[2]) < 0
  K <- if (opposite) unname(-(s[2] * diff_nb[2]) / (s[1] * diff_nb[1])) else NA_real_
  kk <- NULL
  consistent <- NA
  if (!is.null(k)) {
    if (!is.numeric(k) || !setequal(names(k), lev) || any(k <= 0))
      stop("`k` must be two positive numbers named like `diff_nb`.", call. = FALSE)
    kk <- k[lev]
    pop <- sum(s * diff_nb)
    u <- sum(s * kk * diff_nb)
    consistent <- !(sign(u) * sign(pop) < 0)
  }
  structure(list(diff_nb = diff_nb, shares = s, K = K, k = kk, consistent = consistent),
            class = "voipred_k_value")
}

#' @export
print.voipred_k_value <- function(x, ...) {
  lev <- names(x$diff_nb)
  cat("Subgroup difference in net benefit (diff_NB): ", paste(lev, fmt_num(x$diff_nb), collapse = ", "),
      if (is.na(x$K)) " (same sign)" else " (opposite signs)", "\n", sep = "")
  cat("Shares: ", paste0(lev, " ", round(100 * x$shares), "%", collapse = ", "), "\n", sep = "")
  if (is.na(x$K)) {
    cat("The signs are the same: differences in a - c between the subgroups cannot\n",
        "change the population conclusion; there is no K-value.\n", sep = "")
    return(invisible(x))
  }
  cat("K-value (k_", lev[1], " / k_", lev[2], "): ", fmt_num(x$K), "\n", sep = "")
  cat("Population NB and population utility favour different strategies if k_", lev[1],
      " / k_", lev[2], " is beyond ", fmt_num(x$K), " (above it when k_", lev[1], " > k_",
      lev[2], ", below it when k_", lev[1], " < k_", lev[2], ").\n", sep = "")
  if (!is.null(x$k)) {
    r <- unname(x$k[1] / x$k[2])
    cat("With k: k_", lev[1], " / k_", lev[2], " = ", fmt_num(r), ", so population NB and ",
        "population utility ", if (x$consistent) "favour the same strategy"
        else "favour different strategies", ".\n", sep = "")
  }
  invisible(x)
}

# ---- plots -------------------------------------------------------------------

# The subgroup-specific differences in net benefit against the comparator
# (diff_NB) across thresholds, one curve per subgroup (and model), with the
# population difference.
subgroup_curve_spec <- function(x, comparator, shares) {
  chk <- subgroup_check(x, comparator = comparator, shares = shares)
  tab <- chk$table
  several <- length(unique(tab$model)) > 1
  curves <- list()
  for (m in unique(tab$model)) {
    t <- tab[tab$model == m, ]
    for (lev in c(chk$levels, "population")) {
      col <- paste0("diff_NB_", lev)
      lab <- if (lev == "population") "Population" else lev
      if (several) lab <- paste0(m, ": ", lab)
      curves[[length(curves) + 1]] <- curve_data(t$threshold, lab, t[[col]])
    }
  }
  d <- do.call(rbind, curves)
  d$label <- factor(d$label, levels = unique(d$label))
  lev <- levels(d$label)
  new_spec(d, stats::setNames(default_hues(length(lev)), lev), "Threshold Probability",
           paste0("Difference in NB against ", chk$comparator), hline = 0)
}

# The inconsistency region at one threshold (two subgroups), as in Figures 2
# and 5 of Wang et al.: the boundary lines y = -(n1/n2) x and, with k,
# y = -(n1 k1)/(n2 k2) x, the region between them shaded, the observed point
# (diff_NB_1, diff_NB_2), and the line through it, whose slope gives the
# K-value.
subgroup_region_plot <- function(x, threshold, model, comparator, k, shares) {
  chk <- subgroup_check(x, comparator = comparator, k = k, shares = shares)
  lev <- chk$levels
  if (length(lev) != 2)
    stop("The inconsistency region is for two subgroups; with ", length(lev),
         " subgroups use summary(res, k = ) (Wang et al., multiple subgroups).",
         call. = FALSE)
  if (!is.null(chk$k) && chk$k$range)
    stop("The region plot needs values of `k`, not ranges.", call. = FALSE)
  tab <- chk$table
  if (is.null(model)) model <- tab$model[1]
  if (!model %in% tab$model)
    stop("`model` must be one of ", paste0("\"", unique(tab$model), "\"", collapse = ", "),
         ".", call. = FALSE)
  if (missing(threshold) || is.null(threshold) || !is.numeric(threshold) ||
      length(threshold) != 1)
    stop("`threshold` must be one of the thresholds of the result, e.g. threshold = 0.05.",
         call. = FALSE)
  t <- tab[tab$model == model, ]
  j <- which(abs(t$threshold - threshold) < 1e-9)
  if (!length(j))
    stop("`threshold` (", threshold, ") is not one of the thresholds of the result.",
         call. = FALSE)
  px <- t[[paste0("diff_NB_", lev[1])]][j]
  py <- t[[paste0("diff_NB_", lev[2])]][j]
  s <- chk$shares
  s1 <- -s[1] / s[2]
  lim <- 1.6 * max(abs(c(px, py)), 1e-6)
  xs <- c(-lim, lim)
  lines <- data.frame(x = xs, y = s1 * xs, line = "y = -(n1/n2) x")
  shade <- NULL
  if (!is.null(chk$k)) {
    s2 <- -(s[1] * chk$k$lower[1]) / (s[2] * chk$k$lower[2])
    lines <- rbind(lines, data.frame(x = xs, y = s2 * xs, line = "y = -(n1 k1)/(n2 k2) x"))
    shade <- rbind(data.frame(x = c(0, lim, lim), y = c(0, s1 * lim, s2 * lim), part = 1),
                   data.frame(x = c(0, -lim, -lim), y = c(0, -s1 * lim, -s2 * lim), part = 2))
  }
  K <- t$K[j]
  if (!is.na(K))
    lines <- rbind(lines, data.frame(x = xs, y = (py / px) * xs,
                                     line = paste0("observed (K = ", fmt_num(K), ")")))
  lines$line <- factor(lines$line, levels = unique(lines$line))
  cols <- c("black", "steelblue", "firebrick")[seq_along(levels(lines$line))]
  p <- ggplot(lines, aes(x = .data$x, y = .data$y))
  if (!is.null(shade))
    p <- p + ggplot2::geom_polygon(data = shade, aes(group = .data$part), fill = "grey80",
                                   alpha = 0.6, inherit.aes = TRUE)
  p + geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    geom_line(aes(colour = .data$line), linetype = "dashed", linewidth = 0.7) +
    ggplot2::geom_point(data = data.frame(x = px, y = py), size = 3) +
    scale_color_manual(values = cols, name = NULL) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = paste0("Difference in NB in ", lev[1]), y = paste0("Difference in NB in ", lev[2]),
         title = paste0(model, " at threshold ", threshold, ", against ", chk$comparator),
         subtitle = if (!is.null(chk$k))
           paste0("Grey: inconsistency region (k_", lev[1], " / k_", lev[2], " = ",
                  fmt_num(chk$k$lower[1] / chk$k$lower[2]), "); the point is ",
                  if (isTRUE(t$consistent[j])) "outside it" else "inside it")
           else if (!is.na(K)) "Opposite signs; give k to draw the inconsistency region"
           else "Same signs: no inconsistency is possible") +
    theme_minimal() + theme(legend.position = "bottom")
}
