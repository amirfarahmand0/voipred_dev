# voipred.dev

Decision curve analysis and value of information analysis for the validation
of clinical prediction models (development version of voipred).

## Tutorial

**For how to use voipred, go to our tutorial:**
**<https://amirfarahmand0.github.io/voipred/vignettes/voipred.html>**

It walks through a complete analysis, for binary and survival outcomes: net
benefit and decision curves, comparing models, prior information, value of
information analysis (EVPI, P(useful), P(best), EVSI), credible intervals,
informative censoring (cIPCW), competing risks, subgroups, and reporting. The
[package website](https://amirfarahmand0.github.io/voipred/) also has a
Reference page for every function.

## Installation

```r
# install.packages("remotes")
remotes::install_github("amirfarahmand0/voipred_dev", upgrade = "never")
library(voipred.dev)
```

Install (or update) in a new R session (RStudio: Session > Restart R): R
cannot replace a package that is loaded, which gives the warning "package
'voipred.dev' is in use and will not be installed". Running it again updates
voipred.dev when GitHub has a newer version. voipred.dev needs riskRegression
2026.3.11 or later, and a compiler for its C++ code (Rtools on Windows, the
Xcode command line tools on macOS).

## What it does

For a validation sample, from the predicted risks of one model or several
(from any model or software, or the fitted models themselves) and the observed
outcomes (binary, survival or competing risks):

- **Net benefit and decision curves** (`nb_val()`): the net benefit of each
  model, of treating everyone and of treating no one, at each threshold; also
  as net false positives avoided compared with treating everyone. For
  survival outcomes, with Kaplan-Meier or conditional inverse probability of
  censoring weights (cIPCW) for informative censoring.
- **Value of information analysis**: `evpi()` (EVPI, and the probabilities that
  each model is useful and that each strategy is the best) and `evsi()` (the
  value of a further validation study).
- **A Bayesian framework**: the usual net benefit by default, with optional
  prior information, including the informative prior of bayesDCA
  (`nb_informative()`).
- **Subgroups**: whether population-level net benefit can mislead when the
  benefit of treatment differs between subgroups (Wang et al., 2026).
- **Reporting**: `print()`, `summary()`, `plot()` and `as.data.frame()` for
  every result, per 100,000 patients (`per`) if wanted, and optional credible
  intervals (value of information analysis is recommended instead).

```r
library(voipred.dev)
dev <- gusto[1:5000, ]
val <- gusto[5001:5500, ]
fit <- glm(day30 ~ age + Killip + sysbp + pulse, family = binomial, data = dev)
val$p <- predict(fit, newdata = val, type = "response")

res <- nb_val(data = val, predicted_risk = "p", outcome = "day30",
              thresholds = seq(0.01, 0.3, by = 0.01), label = "Logistic")
summary(res)
plot(res)

e <- evpi(res)
summary(e)
plot(res, evpi = e)
```

## Guides for collaborators

Detailed guides that go through every option and check (the tutorial above is
the main documentation):

- [Net benefit guide](vignettes/voipred_nb_guide.Rmd)
  ([rendered](reports/voipred_nb_guide.html))
- [Value of information guide](vignettes/voipred_voi_guide.Rmd)
  ([rendered](reports/voipred_voi_guide.html))
- [Tests and messages](reports/voipred_nb_tests_and_messages.pdf) (October 2026
  version): the tests (checked against dcurves, the modelling packages' own
  predictions, predtools, the EVSI paper's code and bayesDCA) and the error and
  warning messages.

**The guides install and update automatically.** When a guide is knitted, its
first chunk installs the packages it uses that are missing, and installs or
updates voipred.dev from GitHub, before anything is loaded (RStudio's Knit
button runs it in a new R session). `install_github()` asks GitHub for the
latest version, which GitHub allows 60 times an hour per internet address
without a login; on a shared network, a GitHub token raises this
(`usethis::create_github_token()`, then `gitcreds::gitcreds_set()`). If GitHub
cannot be reached, the guides use the installed version.

## Authors

Amirhosein Farahmand (maintainer), Mohsen Sadatsafavi, Laure Wynants.
