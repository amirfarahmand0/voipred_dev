#' GUSTO-I trial data
#'
#' Patients with acute myocardial infarction from the GUSTO-I trial, used for
#' examples and tests. The outcome is 30-day mortality.
#'
#' @format A data frame with 40,830 rows and 29 variables. The outcome is
#'   `day30` (30-day mortality, 0/1). Predictors used in the examples are `age`
#'   (years), `sex`, `Killip` (Killip class I to IV), `sysbp` (systolic blood
#'   pressure, mmHg), `pulse` (heart rate), `pmi` (previous myocardial
#'   infarction) and `miloc` (location of the infarction).
#' @source GUSTO Investigators. An international randomized trial comparing
#'   four thrombolytic strategies for acute myocardial infarction. N Engl J Med
#'   1993;329:673-682.
"gusto"
