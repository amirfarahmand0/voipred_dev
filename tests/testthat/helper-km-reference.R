# The R implementation of the Kaplan-Meier / Aalen-Johansen net benefit that
# the compiled code in src/km_nb.cpp replaces (kept to test that they agree).
nb_core_km_reference <- function(time, status, k, p, thresholds, horizon, w = NULL) {
  if (is.null(w)) w <- rep(1, length(time))
  sw <- sum(w)
  odds <- thresholds / (1 - thresholds)
  risk_all <- km_cuminc(time, status, k, w, horizon)
  nb_model <- numeric(length(thresholds))
  for (j in seq_along(thresholds)) {
    tr <- p > 0 & p >= thresholds[j] & w > 0
    q <- sum(w[tr]) / sw
    if (q == 0) { nb_model[j] <- 0; next }
    if (max(time[tr]) < horizon) { nb_model[j] <- NA_real_; next }
    f <- km_cuminc(time[tr], status[tr], k, w[tr], horizon)
    nb_model[j] <- q * f - q * (1 - f) * odds[j]
  }
  list(NB_model = nb_model, NB_all = risk_all - (1 - risk_all) * odds, prevalence = risk_all)
}
