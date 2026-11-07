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
