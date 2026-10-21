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
validation sample and a range of risk thresholds, `nb_val()` returns the net
benefit of the model, of treating everyone and of treating no one at each
threshold, and `delta_NB`, the gain of the model over the best default
strategy, for a **binary**, **survival** or **competing-risks** outcome.
The result is an object with `print()`, `summary()` and `plot()` methods:
`plot()` draws the decision curve and the ΔNB curve, for one model or
several together.

- **Several models** are validated in one call, on the same patients: a
  named list of risks (`predicted_risk = list(Full = val$p1, Age = val$p2)`),
  several column names (`c(Full = "p1", Age = "p2")`), a data frame or
  matrix with one column per model, or fitted models
  (`model = list(Full = fit1, Age = fit2)`). Results computed separately can
  be combined with `c(res1, res2)`.

- Survival and competing-risks outcomes are evaluated at a time horizon, with
  Kaplan-Meier / Aalen-Johansen estimates or with conditional inverse
  probability of censoring weights (cIPCW), which allow censoring to depend on
  covariates: give them as a formula, `cipcw = ~ age + nodes`, and voipred
  fits a Cox model of the censoring times (with weights truncated at 50);
  `cipcw_control()` is only for other options.
- Only the predicted risks are needed (`predicted_risk`, a vector such as
  `val$p`, or the name of a column of the data such as `"p"`), from any model
  or software. Optionally, the fitted model can be given instead, and
  voipred.dev predicts the risks from it once: the common binary, survival
  and competing-risks models in R work by default, including bayespred and
  bayespred.dev Bayesian models, and any other model with a
  `predict_function` that returns its predicted risks (see `?nb_val`).
- The predictions must be probabilities: scores and class labels are rejected.
- **Prior information (optional)**: `nb_val()` works in a Bayesian framework.
  By default (`nb_improper()`, Beta(0, 0)) it uses the data only, with no
  prior information, and gives the usual net benefit. For a binary outcome,
  prior information can be added with `nb_uniform()`, `nb_jeffreys()` or
  your own Beta priors with `nb_beta()` (the conjugate Beta model of
  bayesDCA).

**Value of information** for validation (Sadatsafavi et al. 2023, 2025), from
the same `nb_val()` results, for every outcome type and any number of models:

- `evpi()`: expected value of perfect information and relative EVPI;
- `evsi()`: expected value of sample information of a future validation study
  of a given size;
- `p_useful()`: the probability that each model is better than treating all
  and treating none;
- `p_best()`: the probability that each strategy (each model, treat all,
  treat none) is the best;
- `evpi()` also computes P(useful) and P(best) from its draws: `p_useful(e)`
  and `p_best(e)` return them, and `plot(e, type = "p_best")` draws them;
- the uncertainty is represented with `n_sim` draws from the posterior of
  the net benefit: for a binary outcome from the conjugate Beta model (and
  `evsi()` uses the conjugate algorithm of Sadatsafavi et al. 2025); for
  survival and competing risks, and for case weights, with the Bayesian
  bootstrap;
- with cIPCW, the censoring model is refitted in every draw, so its
  uncertainty is included: automatically for `cipcw = ~ covariates` and for
  your own `coxph`, `survreg`, `cph`, `psm`, `flexsurvreg` and `prodlim`
  censoring models, with `refit_function` in `cipcw_control()` for others.

**Results**: every result (`nb_val()`, `evpi()`, `evsi()`, `p_useful()`,
`p_best()`) has `print()` (a short table), `summary()` (the main findings),
`as.data.frame()` and `plot()`; `plot(res, evpi = , evsi = )` draws the
decision curve with EVPI and EVSI in a panel below. Every plot takes
`interactive = TRUE` for an interactive (plotly) version with a shared hover
across the curves.

```r
library(voipred.dev)
df_dev_binary <- gusto[1:5000, ]
df_val_binary <- gusto[5001:10000, ]
fit <- glm(day30 ~ age + Killip + sysbp + pulse,
           family = binomial, data = df_dev_binary)
df_val_binary$p <- predict(fit, newdata = df_val_binary, type = "response")

res <- nb_val(data = df_val_binary,
              predicted_risk = df_val_binary$p,      # or the column name: "p"
              outcome = "day30",
              thresholds = seq(0.01, 0.3, by = 0.01),
              label = "Logistic")
print(res)
summary(res)
plot(res)
plot(res, type = "delta")

e <- evpi(res)
p_useful(e)            # from the same draws as EVPI
p_best(e)
o <- evsi(res, future_n = c(1000, 2000, 4000))
summary(o)
plot(o)
plot(res, evpi = e, evsi = o)
plot(res, evpi = e, evsi = o, interactive = TRUE)

# several models in one call
fit_age <- glm(day30 ~ age, family = binomial, data = df_dev_binary)
df_val_binary$p_age <- predict(fit_age, newdata = df_val_binary, type = "response")
both <- nb_val(data = df_val_binary,
               predicted_risk = c(Logistic = "p", `Age only` = "p_age"),
               outcome = "day30",
               thresholds = seq(0.01, 0.3, by = 0.01))
plot(both)
e2 <- evpi(both)
summary(e2)
plot(e2, type = "p_best")

# optional: the model instead of the risks
nb_val(data = df_val_binary,
       model = fit,
       outcome = "day30")

# optional prior information (binary outcome); the default uses the data only
res_u <- nb_val(data = df_val_binary,
                predicted_risk = df_val_binary$p,
                outcome = "day30",
                thresholds = seq(0.01, 0.3, by = 0.01),
                prior = nb_uniform())
evpi(res_u)
```

## Documentation

- [Net benefit guide](vignettes/voipred_nb_guide.Rmd)
  ([rendered](reports/voipred_nb_guide.html)): binary, survival and
  competing-risks examples, predicted risks or a fitted model, the prior
  (the data only by default, or prior information), every option explained,
  plots, and the list of supported models.
- [Value of information guide](vignettes/voipred_voi_guide.Rmd)
  ([rendered](reports/voipred_voi_guide.html)): `p_useful()`, `p_best()`,
  `evpi()` and `evsi()` for binary, survival and competing-risks outcomes,
  with their `print()`, `summary()` and plots, EVPI and EVSI below the
  decision curve, and how to extract values.
- [Tests and messages](reports/voipred_nb_tests_and_messages.pdf): every
  test (checked against `dcurves`, each modelling package's own predictions,
  predtools, the EVSI paper's code and bayesDCA) and every error and warning
  message.
- The help pages give the full details: `?nb_val`, `?nb_priors`,
  `?cipcw_control`, `?evpi`, `?evsi`, `?p_useful`,
  `?p_best`, `?summary.voipred_nb`, `?summary.voipred_voi`,
  `?plot.voipred_nb`, `?plot.voipred_evpi`, `?plot.voipred_evsi`,
  `?plot.voipred_p_useful`, `?plot.voipred_p_best`.

## Next version

Subgroup analysis of the net benefit (subgroup utility heterogeneity).
