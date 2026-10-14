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
`plot_dca()` draws the decision curve and the ΔNB curve, for one model or
several together.

- Survival and competing-risks outcomes are evaluated at a time horizon, with
  Kaplan-Meier / Aalen-Johansen estimates or with conditional inverse
  probability of censoring weights (cIPCW), which allow censoring to depend on
  covariates.
- Only the predicted risks are needed (`predicted_risk`, a vector such as
  `val$p`, or the name of a column of the data such as `"p"`), from any model
  or software. Optionally, voipred.dev predicts
  them once from a fitted model: the common binary, survival and
  competing-risks models in R work automatically, including bayespred and
  bayespred.dev Bayesian models, and any other model with a function that
  returns its predicted risks.
- The predictions must be probabilities: scores and class labels are rejected.
- **Bayesian net benefit** (binary outcome): `nb(..., bayesian = TRUE)` gives
  the posterior means of the conjugate Beta model of bayesDCA, with the prior
  `nb_uniform()` (default), `nb_jeffreys()` or `nb_beta()`.

**Value of information** for validation (Sadatsafavi et al. 2023, 2025), from
the same `nb()` results, for every outcome type and any number of models:

- `evpi()`: expected value of perfect information and relative EVPI;
- `evsi()`: expected value of sample information of a future validation study
  of a given size;
- `p_useful()` and `p_best()`: the probability that a model is better than
  treating all and none, and that each strategy is the best;
- the uncertainty is represented with `method = "bayesian"` (Bayesian
  bootstrap, default) or `"ordinary"` (ordinary bootstrap), with `n_sim`
  draws; `evpi()` also has the asymptotic method for a binary outcome and one
  model. For Bayesian net benefit results the draws come from the posterior,
  and `evsi()` uses the conjugate algorithm of Sadatsafavi et al. (2025);
- with cIPCW, the censoring model is refitted in every draw, so its
  uncertainty is included: automatically for `coxph`, `survreg`, `cph`,
  `psm`, `flexsurvreg` and `prodlim`, with `refit_function` in
  `cipcw_control()` for other models.

**Plots**: `plot_dca()`, `plot_evpi()`, `plot_evsi()`, `plot_p_useful()`,
`plot_p_best()`, and `plot_dca(res, evpi = , evsi = )` for the decision curve
with EVPI and EVSI in a panel below. Every plot function takes
`interactive = TRUE` for an interactive (plotly) version with a shared hover
across the curves.

```r
library(voipred.dev)
fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = gusto[1:5000, ])
val <- gusto[5001:10000, ]
val$p <- predict(fit, val, type = "response")

res <- nb(predicted_risk = val$p, data = val, outcome = "day30",
          thresholds = seq(0.01, 0.3, by = 0.01))
# or by column name: nb(predicted_risk = "p", data = val, ...)
plot_dca(res)
plot_dca(res, type = "delta")

e <- evpi(res)                        # Bayesian bootstrap
evpi(res, method = "ordinary")        # ordinary bootstrap
p_useful(res)
p_best(res)
o <- evsi(res, future_n = c(1000, 2000, 4000))
plot_evsi(o)
plot_dca(res, evpi = e, evsi = o)
plot_dca(res, evpi = e, evsi = o, interactive = TRUE)

nb(model = fit, data = val, outcome = "day30")       # or let voipred.dev predict from the model

# Bayesian net benefit, with draws from its posterior
res_b <- nb(predicted_risk = val$p, data = val, outcome = "day30",
            thresholds = seq(0.01, 0.3, by = 0.01), bayesian = TRUE)
evpi(res_b)
```

## Documentation

- [Net benefit guide](vignettes/voipred_nb_guide.Rmd): binary, survival and
  competing-risks examples, predicted risks or a fitted model, Bayesian net
  benefit, every option explained, plots, and the list of supported models.
- [Tests and messages](reports/voipred_nb_tests_and_messages.pdf): every
  test (checked against `dcurves`, each modelling package's own predictions,
  predtools, the EVSI paper's code and bayesDCA) and every error and warning
  message.
- The help pages give the full details: `?nb`, `?plot_dca`, `?cipcw_control`,
  `?predict_risk`, `?evpi`, `?evsi`, `?p_useful`, `?p_best`,
  `?plot_evpi`, `?plot_evsi`, `?plot_p_useful`, `?plot_p_best`, `?nb_priors`.

## Next version

Subgroup analysis of the net benefit (subgroup utility heterogeneity).
