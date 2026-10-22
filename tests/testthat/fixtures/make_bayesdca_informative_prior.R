# Reference shapes of bayesDCA's threshold-varying ("informative") prior, for
# test-informative.R. Run from the package root with a copy of the bayesDCA
# source (https://github.com/giulianonetto/bayesdca, commit e5dd009) in
# /tmp/bayesdca. Only its .get_threshold_varying_prior_parameters() is used
# (bayesDCA itself, which needs rstan, does not have to be installed).
src <- parse("/tmp/bayesdca/R/utils.R")
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                   as.character(e[[2]]) == ".get_threshold_varying_prior_parameters") eval(e)
# bayesDCA's dca() clamps the thresholds to [1e-9, 0.99] before building the prior
ref_shapes <- function(th, ...) {
  th2 <- unique(pmax(pmin(th, 0.99), 1e-9))
  pr <- .get_threshold_varying_prior_parameters(thresholds = th2, n_strategies = 1, ...)
  cbind(pr$Se1[, 1], pr$Se2[, 1], pr$Sp1[, 1], pr$Sp2[, 1], pr$p1, pr$p2)
}
cases <- list(
  list("default grid 0-0.5", seq(0, 0.5, length = 51), list()),
  list("grid 0.01-0.3", seq(0.01, 0.3, by = 0.01), list()),
  list("grid 0-0.99", seq(0, 0.99, by = 0.01), list()),
  list("no vague region, strength 10", seq(0.01, 0.1, by = 0.01),
       list(ignorance_region_cutpoints = NULL, max_sens_prior_mean = 0.95,
            min_sens_prior_mean = 0.5, max_sens_prior_sample_size = 10)),
  list("custom region and prevalence", seq(0, 0.4, by = 0.02),
       list(ignorance_region_cutpoints = c(0.1, 0.9) * 0.4, ignorance_region_mean = 0.4,
            ignorance_region_sample_size = 3, prev_prior_mean = 0.07,
            prev_prior_sample_size = 100)))
ref <- lapply(cases, function(cs) list(label = cs[[1]], thresholds = cs[[2]],
                                       shapes = do.call(ref_shapes, c(list(cs[[2]]), cs[[3]]))))
saveRDS(ref, "tests/testthat/fixtures/bayesdca_informative_prior.rds")
