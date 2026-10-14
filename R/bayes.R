# Bayesian net benefit for a binary outcome (nb(bayesian = TRUE)): conjugate
# Beta posteriors for the prevalence and, at each threshold, the sensitivity
# and specificity of the model, as in bayesDCA (Cruz et al.); the net benefit
# is reported as its posterior mean.

#' Priors for the Bayesian net benefit
#'
#' Priors for [nb()] with `bayesian = TRUE`, written like bayespred's priors
#' (`flat()`, `jeffreys()`, `log_f()`). The Bayesian net benefit puts Beta
#' priors on the prevalence of the outcome and, at every threshold, on the
#' sensitivity and the specificity of the model (the probabilities of being
#' treated among patients with the event, and of not being treated among
#' patients without it). The net benefit at threshold \eqn{z} is
#' \deqn{NB = \pi \, Se - (1 - \pi)(1 - Sp) \frac{z}{1 - z},}
#' and that of treating all \eqn{\pi - (1 - \pi) z / (1 - z)}. The posteriors
#' are Beta distributions (conjugate), so no MCMC is needed.
#'
#' * `nb_uniform()` (default): Beta(1, 1) for all three, as in bayesDCA.
#' * `nb_jeffreys()`: Jeffreys prior, Beta(1/2, 1/2) for all three.
#' * `nb_beta()`: your own Beta priors, e.g. an informative prior on the
#'   prevalence from earlier studies.
#'
#' @param prevalence,sensitivity,specificity The two shape parameters
#'   \eqn{(a, b)} of each Beta prior, positive numbers: prior mean
#'   \eqn{a / (a + b)}, worth \eqn{a + b} patients. Default Beta(1, 1).
#' @return A prior specification object of class
#'   `c("prior_<type>", "nb_prior")`.
#' @references Cruz GNF, Korthauer K. Bayesian decision curve analysis with
#'   bayesDCA. arXiv 2023.
#' @seealso [nb()]
#' @examples
#' nb_uniform()
#' nb_jeffreys()
#' # an informative prior on the prevalence: mean 0.07, worth 100 patients
#' nb_beta(prevalence = c(7, 93))
#' @name nb_priors
NULL

#' @rdname nb_priors
#' @export
nb_uniform <- function() new_nb_prior("uniform", c(1, 1), c(1, 1), c(1, 1))

#' @rdname nb_priors
#' @export
nb_jeffreys <- function() new_nb_prior("jeffreys", c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5))

#' @rdname nb_priors
#' @export
nb_beta <- function(prevalence = c(1, 1), sensitivity = c(1, 1), specificity = c(1, 1)) {
  for (nm in c("prevalence", "sensitivity", "specificity")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 2 || anyNA(v) || any(!is.finite(v)) || any(v <= 0))
      stop("`", nm, "` must be two positive numbers: the shapes (a, b) of a Beta prior.",
           call. = FALSE)
  }
  new_nb_prior("beta", as.numeric(prevalence), as.numeric(sensitivity),
               as.numeric(specificity))
}

new_nb_prior <- function(type, prevalence, sensitivity, specificity) {
  structure(list(type = type, prevalence = prevalence, sensitivity = sensitivity,
                 specificity = specificity),
            class = c(paste0("prior_", type), "nb_prior"))
}

#' @export
print.nb_prior <- function(x, ...) {
  b <- function(v) sprintf("Beta(%g, %g)", v[1], v[2])
  msg <- switch(x$type,
    uniform  = "uniform, Beta(1, 1) for prevalence, sensitivity and specificity",
    jeffreys = "Jeffreys, Beta(0.5, 0.5) for prevalence, sensitivity and specificity",
    sprintf("Beta: prevalence %s, sensitivity %s, specificity %s",
            b(x$prevalence), b(x$sensitivity), b(x$specificity)))
  cat("Prior:", msg, "\n")
  invisible(x)
}

# ---- internal ----------------------------------------------------------------

# Treated: risk > 0 and risk >= z (as in the frequentist net benefit).
treated_at <- function(p, z) p > 0 & p >= z

# Posterior mean net benefit of one model (the point estimate of nb()): with
# independent Beta posteriors, E[NB] = E[pi] E[Se] - (1 - E[pi])(1 - E[Sp]) z/(1-z).
bayes_nb_mean <- function(y, p, thresholds, prior) {
  n <- length(y)
  d <- sum(y)
  odds <- thresholds / (1 - thresholds)
  a <- prior$prevalence
  ep <- (a[1] + d) / (sum(a) + n)
  tp <- vapply(thresholds, function(z) sum(y == 1 & treated_at(p, z)), numeric(1))
  tn <- vapply(thresholds, function(z) sum(y == 0 & !treated_at(p, z)), numeric(1))
  se <- (prior$sensitivity[1] + tp) / (sum(prior$sensitivity) + d)
  sp <- (prior$specificity[1] + tn) / (sum(prior$specificity) + n - d)
  list(NB_model = ep * se - (1 - ep) * (1 - sp) * odds,
       NB_all = ep - (1 - ep) * odds, prevalence = ep)
}
