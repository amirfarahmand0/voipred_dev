# voipred.dev

Value of information (VoI) analysis for clinical prediction models
(development version of voipred).

## Installation

```r
# install.packages("remotes")
remotes::install_github("amirfarahmand0/voipred_dev", upgrade = "never")
library(voipred.dev)
```

`upgrade = "never"` installs only the missing dependencies and does not stop
to offer updates of packages you already have. voipred.dev needs
riskRegression 2026.3.11 or later: if yours is older, update it first with
`install.packages("riskRegression")`. The
[guide](vignettes/voipred_nb_guide.Rmd) has an installation script that
also installs the packages its examples use.

## Current scope

This version supports **model validation**: given the predicted risks of a
validation sample and a range of risk thresholds, `nb()` returns the net
benefit of the model, of treating everyone and of treating no one at each
threshold, and `delta_NB`, the gain of the model over the best default
strategy, for a **binary**, **survival** or **competing-risks** outcome.
The result is an object with `print()`, `summary()` and `plot()` methods:
`plot()` draws the decision curve and the ΔNB curve, for one model or
several together.

- Survival and competing-risks outcomes are evaluated at a time horizon, with
  Kaplan-Meier / Aalen-Johansen estimates or with conditional inverse
  probability of censoring weights (cIPCW), which allow censoring to depend on
  covariates.
- Only the predicted risks are needed (`predicted_risk`, a vector such as
  `val$p`, or the name of a column of the data such as `"p"`), from any model
  or software. Optionally, the fitted model can be given instead, and
  voipred.dev predicts the risks from it once: the common binary, survival
  and competing-risks models in R work by default, including bayespred and
  bayespred.dev Bayesian models, and any other model with a
  `predict_function` that returns its predicted risks (see `?nb`).
- The predictions must be probabilities: scores and class labels are rejected.
- **Bayesian net benefit**: the net benefit is the posterior mean of the
  conjugate Beta model of bayesDCA. The default prior, `nb_improper()`
  (Beta(0, 0)), gives exactly the usual (frequentist) net benefit; for a
  binary outcome the prior can also be `nb_uniform()`, `nb_jeffreys()` or
  `nb_beta()`.

**Value of information** for validation (Sadatsafavi et al. 2023, 2025), from
the same `nb()` results, for every outcome type and any number of models:

- `evpi()`: expected value of perfect information and relative EVPI;
- `evsi()`: expected value of sample information of a future validation study
  of a given size;
- `p_useful()` and `p_best()`: the probability that a model is better than
  treating all and none, and that each strategy is the best;
- the uncertainty is represented with `n_sim` draws from the posterior of
  the net benefit: for a binary outcome from the conjugate Beta model (and
  `evsi()` uses the conjugate algorithm of Sadatsafavi et al. 2025); for
  survival and competing risks, and for case weights, with the Bayesian
  bootstrap, which is the posterior under the improper prior;
- with cIPCW, the censoring model is refitted in every draw, so its
  uncertainty is included: automatically for `coxph`, `survreg`, `cph`,
  `psm`, `flexsurvreg` and `prodlim`, with `refit_function` in
  `cipcw_control()` for other models.

**Results**: every result (`nb()`, `evpi()`, `evsi()`, `p_useful()`,
`p_best()`) has `print()` (a short table), `summary()` (the main findings),
`as.data.frame()` and `plot()`; `plot(res, evpi = , evsi = )` draws the
decision curve with EVPI and EVSI in a panel below. Every plot takes
`interactive = TRUE` for an interactive (plotly) version with a shared hover
across the curves.

```r
library(voipred.dev)
fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = gusto[1:5000, ])
val <- gusto[5001:10000, ]
val$p <- predict(fit, val, type = "response")

res <- nb(val, predicted_risk = val$p, outcome = "day30",
          thresholds = seq(0.01, 0.3, by = 0.01))
# or by column name: nb(val, "p", "day30", ...)
res
summary(res)
plot(res)
plot(res, type = "delta")

e <- evpi(res)
p_useful(res)
p_best(res)
o <- evsi(res, future_n = c(1000, 2000, 4000))
summary(o)
plot(o)
plot(res, evpi = e, evsi = o)
plot(res, evpi = e, evsi = o, interactive = TRUE)

nb(val, fit, "day30")                 # optional: the model instead of the risks

# another prior (binary outcome)
res_u <- nb(val, val$p, "day30", thresholds = seq(0.01, 0.3, by = 0.01),
            prior = nb_uniform())
evpi(res_u)
```

## Documentation

- [Net benefit guide](vignettes/voipred_nb_guide.Rmd): binary, survival and
  competing-risks examples, predicted risks or a fitted model, Bayesian net
  benefit, every option explained, plots, and the list of supported models.
- [Tests and messages](reports/voipred_nb_tests_and_messages.pdf): every
  test (checked against `dcurves`, each modelling package's own predictions,
  predtools, the EVSI paper's code and bayesDCA) and every error and warning
  message.
- The help pages give the full details: `?nb`, `?nb_priors`,
  `?cipcw_control`, `?evpi`, `?evsi`, `?p_useful`,
  `?p_best`, `?summary.voipred_nb`, `?summary.voipred_voi`,
  `?plot.voipred_nb`, `?plot.voipred_evpi`, `?plot.voipred_evsi`,
  `?plot.voipred_p_useful`, `?plot.voipred_p_best`.

## Next version

Subgroup analysis of the net benefit (subgroup utility heterogeneity).
