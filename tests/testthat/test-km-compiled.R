# The compiled Kaplan-Meier / Aalen-Johansen net benefit (src/km_nb.cpp) gives
# the same results as the R implementation it replaces.

test_that("compiled KM / Aalen-Johansen net benefit equals the R implementation", {
  s <- rott_split()
  v <- s$val
  cx <- survival::coxph(surv_formula, data = s$dev)
  p <- predict_risk(cx, v, horizon = 5)
  z <- c(0, 0.01, seq(0.05, 0.95, by = 0.05), 0.99)
  tt <- v$cr_time
  st <- v$cr_code
  set.seed(1)
  for (b in 1:40) {
    w <- switch(b %% 4 + 1, rep(1, nrow(v)), rexp(nrow(v)), rpois(nrow(v), 1),
                rbinom(nrow(v), 1, 0.2) * rexp(nrow(v)))
    for (k in 1:2) {
      a <- suppressWarnings(nb_core_km(tt, st, k, p, z, 5, w))
      r <- nb_core_km_reference(tt, st, k, p, z, 5, w)
      expect_equal(a, r, tolerance = 1e-14)
      expect_identical(is.na(a$NB_model), is.na(r$NB_model))
    }
  }
  # survival (one cause), ties in the times, a risk of exactly 0, and a horizon
  # before any event
  tt2 <- round(v$dtimey, 1)
  p2 <- replace(p, 1:20, 0)
  a <- suppressWarnings(nb_core_km(tt2, v$death, 1, p2, z, 5, NULL))
  expect_equal(a, nb_core_km_reference(tt2, v$death, 1, p2, z, 5, NULL), tolerance = 1e-14)
  early <- nb_core_km(tt2, v$death, 1, p2, z, 0.01, NULL)
  expect_equal(early, nb_core_km_reference(tt2, v$death, 1, p2, z, 0.01, NULL), tolerance = 1e-14)
})
