# voipred

Value of information (VoI) analysis for clinical prediction models.

## Current scope

This version supports **model validation**: given a fitted model, a validation
data set and a range of risk thresholds, `nb()` returns the net benefit of the
model, of treating everyone and of treating no one at each threshold, for a
**binary**, **survival** or **competing-risks** outcome.

- Survival and competing-risks outcomes are evaluated at a time horizon, with
  Kaplan-Meier / Aalen-Johansen estimates or with conditional inverse
  probability of censoring weights (cIPCW), which allow censoring to depend on
  covariates.
- Predicted risks are obtained from the model itself. The common binary,
  survival and competing-risks models work automatically; any other model can
  be used by giving a function that returns its predicted risks.
- The predictions must be probabilities: scores and class labels are rejected,
  and the validation data must contain every covariate the model uses.

Results are checked against `dcurves` and against each modelling package's own
predictions in the test suite.

Model development and the uncertainty and VoI measures (ΔNB, P(useful),
P(best), EVPI, EVSI) will be added in later versions.
