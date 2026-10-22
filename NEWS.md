# voipred.dev 0.0.2.9000

## Breaking changes

* `nb()` is renamed `nb_val()`.
* `nb_val()` has one input argument, `predicted_risk`: a vector of risks, a
  column name or a fitted model, or, for several models, a `list()` of these.
  The `model` argument is removed (give the fitted model as `predicted_risk`),
  and so are data frames and matrices of risks; `c()` of several models is
  refused with a message pointing to `list()`.
* The result of `nb_val()` holds one or more models: `models` (the labels)
  replaces `label`, `risk` is a matrix with one column per model, and
  `estimates` and `counts` have a `model` column.
* The `censoring` argument of `nb_val()` is removed: cIPCW is used when
  `cipcw` is given, Kaplan-Meier otherwise.
* cIPCW truncates censoring probabilities at `min_prob = 0.02` (weights at
  most 50) by default; `min_prob = NULL` turns truncation off (the previous
  default).
* `nb_beta()`: the parts not given are Beta(0, 0) (no prior information)
  instead of Beta(1, 1); `nb_beta()` is `nb_improper()`.

## New features

* Several models in one call of `nb_val()`, as a list that can mix fitted
  models, column names and vectors of risks. `print()`, `summary()`, `plot()`
  and the value-of-information functions compare all of them; the models are
  drawn jointly. `c()` combines results computed separately.
* Labels follow one procedure: `label`, the names in the list, the column
  name, the name in the call, then `"Model 1"`, `"Model 2"`, ...
* `predict_function` can be one function for every fitted model, or a list
  named by model. `nb_val()` stops before predicting if voipred has no way to
  predict from a model and no function is given for it, naming the models.
* `cipcw = ~ covariates`: voipred fits a Cox model of the censoring times on
  those covariates (refitted in every draw of the value-of-information
  functions). `cipcw` also takes a fitted censoring model, and
  `cipcw_control()` is optional, for the other settings.
* `evpi()` also computes P(useful) and P(best) from its draws: `p_useful(e)`,
  `p_best(e)`, `plot(e, type = "p_useful")` and `plot(e, type = "p_best")`;
  `summary(e)` gives the probability that the strategy that is best now is
  the best. `p_useful()` takes several models.
* `nb_informative()`: the informative prior of bayesDCA, whose priors on the
  sensitivity and specificity change with the threshold (binary outcomes),
  as `nb_val(prior = nb_informative())`, with bayesDCA's defaults; its shapes
  match bayesDCA's exactly. `prior_shapes()` gives the prior shapes of any
  prior at each threshold.
* Plots take optional colours, by label or in order: `colors =`.
* Clear messages when plotly (interactive plots) or patchwork (stacked
  panels) is not installed, with how to install them.
* The joint posterior of several models stops with a message before it would
  need more than about 2 GB of memory.
