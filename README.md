# voipred

Value of information (VoI) analysis for clinical prediction models.

## Current scope

This version supports **model validation** for a **binary outcome** only.

Given a fitted model, a validation data set and a range of risk thresholds,
`nb()` returns the net benefit of the model, of treating everyone and of
treating no one at each threshold. Supported models are logistic regression
(`glm`), random forests (`ranger`, `randomForest`) and Bayesian logistic
regression from bayespred. The validation data must contain every covariate
the model was developed with.

Results are checked against `dcurves` in the test suite.

Model development, survival outcomes and the uncertainty and VoI measures
(ΔNB, P(useful), P(best), EVPI, EVSI) will be added in later versions.


