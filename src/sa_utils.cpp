//---------------------Simulated annealing utility functions -------------------

#include "sa_utils.h"
#include "objective_utils.h"
#include "update_candidate.h"

#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

//---------------------- Local helper functions --------------------------------

namespace {

  static inline double clamp01(double x, double eps = 1e-15) {
    if (x < eps) return eps;
    if (x > 1.0 - eps) return 1.0 - eps;
    return x;
  }

  static double log_sum_exp_negE_over_T(const std::vector<double>& E, double invT) {
    if (E.empty()) return -std::numeric_limits<double>::infinity();
    double m = -std::numeric_limits<double>::infinity();
    for (double e : E) {
      const double v = -e * invT;
      if (v > m) m = v;
    }
    double s = 0.0;
    for (double e : E) s += std::exp(-e * invT - m);
    return m + std::log(s);
  }

}

//--------------------- Main functions -----------------------------------------

InitTempResult estimate_initial_temperature_benameur_cpp(
    std::vector<double> X_seed,
    const std::vector<double>& W,
    const std::vector<double>& U,
    const std::vector<double>& C,
    const std::vector<double>& E,
    const std::vector<double>& O,
    const std::vector<double>& SD,
    const std::vector<double>& SxH,
    const std::vector<int>&    D,
    int    max_disp_thres,
    int    disp_boundary,
    int    n_h,
    int    n_s,
    int    dim_x,
    int    dim_y,
    double alpha_scaled,
    double beta_scaled,
    double gamma_scaled,
    double base_prob_X0,
    double step_proportion,
    double step_probability,
    int    num_samples,
    double chi0,                // target acceptance
    double p,                   // damping exponent in Ben-Ameur’s update (≥1)
    double tol_logchi,          // tolerance on |log χ̂ − log χ0|
    int    max_iters,           // cap on fixed-point iterations
    double T1,                  // optional initial guess; ≤0 => use heuristic
    int    max_tries_factor,    // limit on total neighbour proposals during
    int    rng_seed,
    bool   verbose
) {
  if (!(chi0 > 0.0 && chi0 < 1.0)) throw std::invalid_argument("chi0 must be in (0,1)");
  if (p < 1.0) throw std::invalid_argument("p must be >= 1");
  if (num_samples <= 0) throw std::invalid_argument("num_samples must be > 0");

  // Set random seeds
  std::mt19937 rng;
  if (rng_seed >= 0) {
    rng.seed(static_cast<uint32_t>(rng_seed));
  } else {
    std::random_device rd; rng.seed(rd());
  }
  std::uniform_real_distribution<> runif(0.0, 1.0);

  // Vectors for storing uphill transitions
  std::vector<double> Emins, Emaxs, dpos;
  Emins.reserve(num_samples);
  Emaxs.reserve(num_samples);
  dpos.reserve(num_samples);

  // Initial
  double E_base = compute_H(X_seed, C, E, O, SD, SxH, D,
                            alpha_scaled, beta_scaled, gamma_scaled,
                            n_h, dim_x, dim_y, max_disp_thres, disp_boundary).H;

  const int max_tries = std::max(num_samples * max_tries_factor, 100);
  int tries = 0;

  // Tune initial temp
  while ((int)Emins.size() < num_samples && tries++ < max_tries) {
    std::vector<double> cand = update_candidate(W, U, X_seed,
                                                step_proportion, step_probability,
                                                n_h, dim_x, dim_y);
    if (cand == X_seed) continue;

    double E_cand = compute_H(cand, C, E, O, SD, SxH, D,
                              alpha_scaled, beta_scaled, gamma_scaled,
                              n_h, dim_x, dim_y, max_disp_thres, disp_boundary).H;

    const double dE = E_cand - E_base;
    if (dE > 0.0) {
      Emins.push_back(E_base);
      Emaxs.push_back(E_cand);
      dpos.push_back(dE);
    }

    // random walk to diversify bases
    if (runif(rng) < 0.5) { X_seed.swap(cand); E_base = E_cand; }
  }

  InitTempResult res{};
  res.diag.tries = tries;

  if (Emins.empty()) {
    // No uphill moves observed -> tiny T0 and diagnostics
    res.T0 = 1e-6;
    res.diag.samples_used  = 0;
    res.diag.chi0          = chi0;
    res.diag.chi_hat_final = 0.0;
    res.diag.iters         = 0;
    res.diag.dE_mean       = std::numeric_limits<double>::quiet_NaN();
    res.diag.dE_median     = std::numeric_limits<double>::quiet_NaN();
    if (verbose) std::cerr << "[init_temp] No uphill transitions found; returning T0=1e-6\n";
    return res;
  }

  // stats
  const double dE_mean = std::accumulate(dpos.begin(), dpos.end(), 0.0) / (double)dpos.size();
  std::nth_element(dpos.begin(), dpos.begin() + dpos.size()/2, dpos.end());
  const double dE_median = dpos[dpos.size()/2];

  // Initial temp guess
  double T = T1;
  if (!(T > 0.0)) {
    T = dE_median / (-std::log(chi0));
    if (!std::isfinite(T) || T <= 0.0) T = 1.0;
  }

  // chi_hat(T)
  auto chi_hat = [&](double Tcur) {
    const double invT = 1.0 / Tcur;
    const double lnum = log_sum_exp_negE_over_T(Emaxs, invT);
    const double lden = log_sum_exp_negE_over_T(Emins, invT);
    const double val  = std::exp(lnum - lden);
    return clamp01(val);
  };

  // Ben-Ameur iterations
  const double log_chi0 = std::log(clamp01(chi0));
  int iters = 0;
  for (; iters < max_iters; ++iters) {
    const double ch  = chi_hat(T);
    const double lch = std::log(clamp01(ch));

    if (std::fabs(lch - log_chi0) <= tol_logchi) break;

    double ratio = (log_chi0 != 0.0) ? (lch / log_chi0) : 1.0; // both negative
    if (!(ratio > 0.0)) ratio = 1.0;
    T *= std::pow(ratio, 1.0 / p);

    if (!std::isfinite(T) || T <= 1e-16) T = 1e-6;
    if (T > 1e16) T = 1e16;
  }

  res.T0 = T;
  res.diag.iters         = iters;
  res.diag.samples_used  = (int)Emins.size();
  res.diag.chi0          = chi0;
  res.diag.chi_hat_final = chi_hat(T);
  res.diag.dE_mean       = dE_mean;
  res.diag.dE_median     = dE_median;

  if (verbose) {
    std::cout << "[init_temp] T0=" << res.T0
              << " | chi_hat(T0)=" << res.diag.chi_hat_final
              << " | samples=" << res.diag.samples_used
              << " | iters=" << res.diag.iters
              << " | dE_med=" << res.diag.dE_median
              << " | dE_mean=" << res.diag.dE_mean
              << std::endl;
  }
  return res;
}
