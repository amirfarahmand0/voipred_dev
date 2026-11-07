# voipred

voipred (the R package `voipred.dev`) evaluates clinical prediction models in
a validation sample by their **clinical utility**, and quantifies how much the
remaining uncertainty matters for the decision of whether to use them.

- **Net benefit and decision curves** (`nb_val()`), for binary, survival and
  competing-risk outcomes, for one model or several compared on the same
  patients; with Kaplan-Meier or conditional inverse probability of
  censoring weights (cIPCW) for informative censoring.
- **Value of information analysis**: the expected value of perfect
  information (EVPI), the probability that each model is useful and that
  each strategy is the best (`evpi()`), and the expected value of sample
  information of a further validation study (`evsi()`).
- Results as net true positives (net benefit) or as net false positives
  avoided compared with treating everyone.
- **Subgroups**: whether population-level net benefit can mislead when the
  benefit of treatment differs between subgroups (Wang et al., 2026).
- A **Bayesian framework**, which by default gives the usual net benefit,
  with optional prior information, including the informative prior of
  bayesDCA (`nb_informative()`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("amirfarahmand0/voipred_dev")
```

## Get started

The article [Introduction to voipred](vignettes/voipred.qmd) walks through a
complete analysis: net benefit and decision curves, comparing models, value
of information, prior information, survival and competing-risk outcomes, and
reporting. The Reference (in the menu) documents every function.

## Authors

- **Amirhosein Farahmand**, author and maintainer
- **Mohsen Sadatsafavi**, author
- **Laure Wynants**, author

To cite voipred, see [Citation](CITATION.html).
