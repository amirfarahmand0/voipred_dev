# voipred.dev

Value of information (VoI) analysis for clinical prediction models
(development version of voipred).

## Installation

```r
# install.packages("remotes")
remotes::install_github("amirfarahmand0/voipred_dev")
library(voipred.dev)
```

## Current scope

This version supports **model validation**: given the predicted risks of a
validation sample and a range of risk thresholds, `nb()` returns the net benefit of the
model, of treating everyone and of treating no one at each threshold, and
`delta_NB`, the gain of the model over the best default strategy, for a
**binary**, **survival** or **competing-risks** outcome. `plot_dca()` draws the
decision curve and the ΔNB curve, for one model or several together.

- Survival and competing-risks outcomes are evaluated at a time horizon, with
  Kaplan-Meier / Aalen-Johansen estimates or with conditional inverse
  probability of censoring weights (cIPCW), which allow censoring to depend on
  covariates.
- Only the predicted risks are needed (`predicted_risk`, a vector or a column
  of the data), from any model or software. Optionally, voipred.dev predicts
  them once from a fitted model: the common binary, survival and
  competing-risks models in R work automatically, including bayespred and
  bayespred.dev Bayesian models, and any other model with a function that
  returns its predicted risks.
- The predictions must be probabilities: scores and class labels are rejected.
- **Bayesian net benefit** (binary outcome): `nb(..., bayesian = TRUE)` gives
  the posterior means of the conjugate Beta model of bayesDCA, with the prior
  `nb_uniform()` (default), `nb_jeffreys()` or `nb_beta()`.

```r
library(voipred.dev)
fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = gusto[1:5000, ])
val <- gusto[5001:10000, ]
val$p <- predict(fit, val, type = "response")

res <- nb(predicted_risk = "p", data = val, outcome = "day30",
          thresholds = seq(0, 0.3, by = 0.01))
plot_dca(res)
plot_dca(res, type = "delta")

nb(model = fit, data = val, outcome = "day30")       # or let voipred.dev predict from the model
nb("p", val, "day30", bayesian = TRUE)               # Bayesian net benefit (posterior means)
```

## Documentation

- [Net benefit guide](vignettes/voipred_nb_guide.Rmd): binary, survival and
  competing-risks examples, every option explained, plots, and the list of
  supported models ([rendered HTML](reports/voipred_nb_guide.html)).
- [Tests and messages](reports/voipred_nb_tests_and_messages.pdf): every
  test (checked against `dcurves` and each modelling package's own
  predictions) and every error and warning message.
- The help pages give the full details: `?nb`, `?plot_dca`, `?cipcw_control`,
  `?predict_risk`, `?nb_priors`.

## Next version

The VoI and uncertainty measures: P(useful), P(best), EVPI and EVSI.
