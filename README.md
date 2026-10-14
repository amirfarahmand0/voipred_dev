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

This version supports **model validation**: given a fitted model, a validation
data set and a range of risk thresholds, `nb()` returns the net benefit of the
model, of treating everyone and of treating no one at each threshold, and
`delta_NB`, the gain of the model over the best default strategy, for a
**binary**, **survival** or **competing-risks** outcome. `plot_dca()` draws the
decision curve and the ΔNB curve, for one model or several together.

- Survival and competing-risks outcomes are evaluated at a time horizon, with
  Kaplan-Meier / Aalen-Johansen estimates or with conditional inverse
  probability of censoring weights (cIPCW), which allow censoring to depend on
  covariates.
- Predicted risks are obtained from the model itself. The common binary,
  survival and competing-risks models in R work automatically, including
  bayespred and bayespred.dev Bayesian models; any other model, from any
  software, can be used by giving a function that returns its predicted
  risks.
- The predictions must be probabilities: scores and class labels are rejected,
  and the validation data must contain every covariate the model uses.

```r
library(voipred.dev)
fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = gusto[1:5000, ])
res <- nb(fit, gusto[5001:10000, ], outcome = "day30", thresholds = seq(0, 0.3, by = 0.01))
plot_dca(res)
plot_dca(res, type = "delta")
```

## Documentation

- [Net benefit guide](vignettes/voipred_nb_guide.Rmd): binary, survival and
  competing-risks examples, every option explained, plots, and the list of
  supported models ([rendered HTML](reports/voipred_nb_guide.html)).
- [Tests and messages](reports/voipred_nb_tests_and_messages.pdf): every
  test (checked against `dcurves` and each modelling package's own
  predictions) and every error and warning message.
- The help pages give the full details: `?nb`, `?plot_dca`, `?cipcw_control`,
  `?predict_risk`.

## Next version

The VoI and uncertainty measures: P(useful), P(best), EVPI and EVSI.
