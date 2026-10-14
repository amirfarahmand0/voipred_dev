#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Work buffers, allocated once per call of km_nb_sorted() and reused.
struct Buffers {
  std::vector<double> t, at, dall, dk, n_risk;
  explicit Buffers(int n) : t(n), at(n), dall(n), dk(n), n_risk(n) {}
};

// Kaplan-Meier (Aalen-Johansen for competing risks) cumulative incidence of
// cause k by tau among the patients with keep[i], for patients sorted by time.
// Same arithmetic, in the same order, as the R function km_cuminc().
static double cuminc_sorted(const double* time, const int* status, int k, const double* w,
                            const std::vector<char>& keep, double tau, int n, Buffers& b) {
  int G = 0, i = 0;
  while (i < n) {
    if (!keep[i]) { i++; continue; }
    double t = time[i], at = 0, d_all = 0, d_k = 0;
    int j = i;
    while (j < n && time[j] == t) {
      if (keep[j]) {
        at += w[j];
        d_all += w[j] * (status[j] != 0);
        d_k += w[j] * (status[j] == k);
      }
      j++;
    }
    b.t[G] = t; b.at[G] = at; b.dall[G] = d_all; b.dk[G] = d_k;
    G++;
    i = j;
  }
  if (G == 0) return 0.0;
  // R's cumsum() and cumprod() accumulate in long double and store doubles
  long double cum = 0;
  for (int g = G - 1; g >= 0; g--) { cum += b.at[g]; b.n_risk[g] = (double) cum; }
  double s_before = 1.0, inc = 0.0;
  long double prod = 1.0, sum = 0.0;
  for (int g = 0; g < G && b.t[g] <= tau; g++) {
    double term = s_before * b.dk[g] / b.n_risk[g];
    sum += term;
    inc = (double) sum;
    prod *= (1 - b.dall[g] / b.n_risk[g]);
    s_before = (double) prod;
  }
  return inc;
}

// Net benefit of the model (and the risk by the horizon) at every threshold,
// for patients sorted by time: treated = risk > 0 & risk >= z & w > 0;
// NB = q F - q (1 - F) z / (1 - z), 0 when nobody is treated, NA when no
// treated patient is followed up to the horizon.
// [[Rcpp::export]]
List km_nb_sorted(NumericVector time, IntegerVector status, int k, NumericVector risk,
                  NumericVector w, NumericVector thresholds, double horizon) {
  int n = time.size(), Z = thresholds.size();
  long double sw_l = 0;               // R's sum() accumulates in long double
  for (int i = 0; i < n; i++) sw_l += w[i];
  double sw = (double) sw_l;
  const double* tp = time.begin();
  const int* sp = status.begin();
  const double* wp = w.begin();
  Buffers buf(n);
  std::vector<char> keep(n);
  for (int i = 0; i < n; i++) keep[i] = w[i] > 0;
  double risk_all = cuminc_sorted(tp, sp, k, wp, keep, horizon, n, buf);
  NumericVector nb_model(Z);
  LogicalVector undefined(Z);
  for (int jz = 0; jz < Z; jz++) {
    double z = thresholds[jz], tmax = R_NegInf;
    long double qw_l = 0;
    for (int i = 0; i < n; i++) {
      keep[i] = risk[i] > 0 && risk[i] >= z && w[i] > 0;
      if (keep[i]) { qw_l += w[i]; if (time[i] > tmax) tmax = time[i]; }
    }
    double qw = (double) qw_l;
    if (qw == 0) { nb_model[jz] = 0; continue; }
    if (tmax < horizon) { nb_model[jz] = NA_REAL; undefined[jz] = true; continue; }
    double q = qw / sw;
    double f = cuminc_sorted(tp, sp, k, wp, keep, horizon, n, buf);
    nb_model[jz] = q * f - q * (1 - f) * (z / (1 - z));
  }
  return List::create(_["NB_model"] = nb_model, _["risk_all"] = risk_all,
                      _["undefined"] = undefined);
}
