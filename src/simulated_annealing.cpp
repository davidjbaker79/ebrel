// -------------------------------Simulated annealing----------------------------

#include "simulated_annealing.h"
#include "dispersal_utils.h"
#include "update_candidate.h"
#include "objective_utils.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
#include <stdexcept>

/* Simulated annealing algorithm with:
  - option for Lam et al style tuning schedule [Lam, J. and Delosme, J.M., 1988.
    An efficient simulated annealing schedule: derivation. New Haven, CT: Yale
    Electrical Engineering Department, 8816.]
  - or geometric cooling schedule with a option to control rate of cooling
*/

SAResult simulated_annealing(
    // --------- Input data and parameters
    std::vector<double> X0,
    const std::vector<double>& W,
    const std::vector<double>& U,
    const std::vector<double>& C,
    const std::vector<double>& E,
    const std::vector<double>& O,
    const std::vector<double>& SD,
    const std::vector<double>& SxH,
    const std::vector<int>& D,
    int max_disp_thres,
    int disp_boundary,
    int n_h,
    int n_s,
    int dim_x,
    int dim_y,
    // ---------- Objective function related
    double alpha_scaled = 1,
    double beta_scaled = 25,
    double gamma_scaled = 100,
    // ---------- Simulated Annealing
    double step_proportion = 0.05,
    double step_probability = 0.05,
    int n_iterations = 10000,
    double temp = 2000,
    double cooling_rate_c = 1,        // constant for tuning cooling rate
    // ---------- LAM: user controls (default = disabled) ----------
    bool   lam_enabled = false,       // turn Lam-style online control on/off
    double lam_target_mid = 0.44,     // target uphill acceptance during early/mid run
    double lam_target_final = 0.05,   // target uphill acceptance at the end
    double lam_hold_frac = 0.60,      // fraction of run to hold lam_target_mid before decaying
    double lam_p = 2.0,               // damping exponent in Ben-Ameur correction (>=1)
    // ---------- Early stopping
    int min_iterations = 1000,        // require at least this many iterations
    int acceptance_window = 200,      // window length for acceptance rate
    double acceptance_thres = 0.001,  // "low" acceptance threshold
    int iter_no_improve = 1000,       // consecutive iters with no meaningful improvement
    double improve_eps = 1e-6,        // relative improvement needed to reset patience
    // --------- Verbose set to false
    bool verbose = false
) {

  // Printing control
  std::ios_base::sync_with_stdio(false);
  std::cout.tie(nullptr);

  // RNG
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<> runif(0.0, 1.0);

  // Histories
  std::vector<double> H_history;
  std::vector<double> F_history;
  std::vector<double> F1_history;
  std::vector<double> F2_history;

  // Initial evaluation of X0 (only need H here)
  double  best_score = compute_H(X0, C, E, O, SD, SxH, D,
                                alpha_scaled, beta_scaled, gamma_scaled,
                                n_h, dim_x, dim_y, max_disp_thres, disp_boundary).H;

  // Current & best states
  std::vector<double> best = X0;
  std::vector<double> curr = X0;
  double curr_eval = best_score;
  std::vector<double> g_best(n_s, 0.0);

  // Ensure stopping parameters are in a 'sensible' range
  const int    min_iters    = std::max(0,  min_iterations);
  const int    win          = std::max(1,  acceptance_window); // avoid % 0
  const double acc_min      = std::clamp(acceptance_thres, 0.0, 1.0); // keep in [0,1]
  const int    patience     = std::max(1,  iter_no_improve);
  const double improve_thr  = std::max(0.0, improve_eps);

  // Diagnostics to check if it is freezing early
  std::vector<double> acc_rate_trace; // per window
  int    early_stop_iter = -1;
  int    attempted_total = 0;
  int    accepted_total  = 0;

  // For tracking combined stopping criterion
  int    attempted_in_win = 0; // tracks whether solution was attempted
  int    accepted_in_win  = 0; // tracks whether solution was accepted
  int    no_improve       = 0; // tracks improvements
  double last_best        = best_score;   // last H(x) that counted as an improvement

  // ---------- LAM: helpers & counters ----------
  auto clamp01 = [](double x, double eps = 1e-12) {
    if (x < eps) return eps;
    if (x > 1.0 - eps) return 1.0 - eps;
    return x;
  };
  auto target_acceptance = [&](double progress01) {
    progress01 = std::min(1.0, std::max(0.0, progress01));
    if (progress01 <= lam_hold_frac) return lam_target_mid;
    // log-linear decay from mid to final across the tail
    double u = (progress01 - lam_hold_frac) / (1.0 - lam_hold_frac);
    return std::exp((1.0 - u) * std::log(lam_target_mid) + u * std::log(lam_target_final));
  };
  int uphill_attempted_in_win = 0;
  int uphill_accepted_in_win  = 0;
  double temp_curr = temp; // current temperature (used if lam_enabled)


  // Main loop over iterations
  for (int z = 0; z < n_iterations; ++z) {

    // Propose candidate using distance weighted approach
    std::vector<double> candidate = update_candidate(W, U, curr,
                                                     step_proportion,
                                                     step_probability,
                                                     n_h, dim_x, dim_y);

    // If no difference, skip so no-op proposals don't affect acceptance stats
    // Doesn't seem valuable to count candidates with no change as an iteration
    if (candidate == curr) {
      continue;
    }

    attempted_in_win++; // iterate forward # of real propose candidates attempted
    attempted_total++;  // count proposals

    // Calculate scores for candidate
    HResult scores = compute_H(candidate, C, E, O, SD, SxH, D,
                               alpha_scaled, beta_scaled, gamma_scaled,
                               n_h, dim_x, dim_y, max_disp_thres, disp_boundary);
    double candidate_eval = scores.H; // score for candidate

    // Trace histories (per proposal)
    H_history.push_back(candidate_eval);
    F_history.push_back(scores.Fx);
    F1_history.push_back(scores.F1);
    F2_history.push_back(scores.F2);

    // Update best if improved
    if (candidate_eval < best_score) {
      best = candidate; // if candidate H(x) < current best then update best to candidate
      best_score = candidate_eval; // update scores
      g_best = scores.g; // track g
      // std::cout  << "Iteration = " << z
      //            << " | Current Best H(x) = " << best_score
      //            << " | Current Best F(x) = " << scores.Fx << std::endl;
    }

    // ------------- Metropolis acceptance -------------------------------

    // Energy difference and uphill/downhill/flat classification
    const double delta = candidate_eval - curr_eval;  // >0 = uphill (worse), <0 = downhill (better)
    const double epsE  = 1e-12 * (1.0 + std::max(std::abs(curr_eval), std::abs(candidate_eval)));

    const bool downhill = (delta < -epsE);
    const bool uphill   = (delta >  epsE);
    const bool flat     = !downhill && !uphill;  // treat tiny diffs as ties to avoid noise effects

    // LAM: count uphill attempts (only when meaningfully uphill)
    if (uphill) uphill_attempted_in_win++;

    // Temperature approach
    const double t = (lam_enabled)
      ? temp_curr
    : (temp / (1.0 + cooling_rate_c * z));  // default harmonic schedule

    // Metropolis acceptance in log-space to avoid under/overflow issues
    bool accepted = false;

    if (downhill || flat) {
      // Always accept downhill; accept flats to avoid noise effects
      accepted = true;
    } else if (t > 0.0 && std::isfinite(t)) {
      double u = runif(rng);
      if (u <= 0.0) u = std::numeric_limits<double>::min(); // avoid log(0)
      const double logu = std::log(u);     // (-inf, 0]
      const double logA = -delta / t;      // log acceptance ratio
      accepted = (logu < std::min(0.0, logA));
    } else {
      accepted = false;
    }

    // ---- Apply acceptance, update counters ----
    if (accepted) {
      curr = candidate; // update to candidate
      curr_eval = candidate_eval; // canditate H(x)
      accepted_in_win++; // track acceptance window
      accepted_total++; // total accepted in window
      if (uphill) uphill_accepted_in_win++; // uphill (i.e. search mode)
    }

    // ---- Combined stopping criterion components ----

    // Using a joint stopping criterion that tracks both the relative magnitude
    // of improvement and acceptance rate.

    // -- No meaningful improvement tracker (relative to magnitude)
    {
      double denom = std::max(1.0, std::abs(last_best));
      double rel_delta = std::abs(best_score - last_best) / denom;
      if (rel_delta >= improve_thr) {
        no_improve = 0;
        last_best = best_score;
      } else {
        no_improve++;
      }
    }

    // -- Low acceptance over a window and Lam cooling schedule controller
    if ((z + 1) % win == 0) {

      // Acceptance rate controller
      double acc_rate = (attempted_in_win > 0)
      ? static_cast<double>(accepted_in_win) / attempted_in_win
      : 0.0;

      // record diagnostics trace
      acc_rate_trace.push_back(acc_rate); // acceptance trace

      // Early stopping
      if ((z + 1) >= min_iters && acc_rate < acc_min && no_improve >= patience) {
        std::cout << "Stopping early at iteration " << (z + 1)
                  << " due to combined criterion: low acceptance AND no improvement.\n";
        early_stop_iter = (z + 1);
        break;
      }

      // -- LAM adaptive cooling schedule (e.g.update temp_curr)
      // lowers the temperature at every step and keeps the system in
      // quasi-equilibrium at all times, using quasi-equilibrium criterion (i.e.
      // he controller will push t down if uphill acceptance exceeds the target,
      // and up if it’s below the target based on the target acceptance curve)
      if (lam_enabled) {
        // Estimate uphill acceptance in this window
        double hat_up = (uphill_attempted_in_win > 0)
        ? static_cast<double>(uphill_accepted_in_win) / uphill_attempted_in_win
        : ((attempted_in_win > 0) ? static_cast<double>(accepted_in_win) / attempted_in_win : 0.0);

        double progress = static_cast<double>(z + 1) / static_cast<double>(n_iterations);
        double a_star   = target_acceptance(progress);

        double num = std::log(clamp01(hat_up));
        double den = std::log(clamp01(a_star));
        double ratio = (den != 0.0) ? (num / den) : 1.0; // both logs negative in (0,1)
        if (!(ratio > 0.0) || !std::isfinite(ratio)) ratio = 1.0;

        temp_curr *= std::pow(ratio, 1.0 / lam_p);
        if (!std::isfinite(temp_curr) || temp_curr < 1e-9) temp_curr = 1e-9;
        if (temp_curr > 1e16) temp_curr = 1e16;
      }

      // reset window counters for stopping criteria
      attempted_in_win = 0;
      accepted_in_win  = 0;

      // LAM: reset uphill counters
      uphill_attempted_in_win = 0;
      uphill_accepted_in_win  = 0;

    }

    // Print options
    if (verbose && ((z + 1) % win == 0)) {
      std::cout << "[SA] it=" << (z + 1)
                << " T=" << (lam_enabled ? temp_curr : (temp / (1.0 + cooling_rate_c * z)))
                << " Best H(x) = " << best_score
                << " Best F(x) = " << scores.Fx
                << '\n';
    }

  }

  // Result of SA
  SAResult out;
  out.X_best   = best;
  out.H_best   = best_score;
  out.H_trace  = H_history;
  out.F_trace  = F_history;
  out.F1_trace = F1_history;
  out.F2_trace = F2_history;
  out.g_best = g_best;

  // Diagnostics (NEW)
  out.diag.acceptance_rate_trace = std::move(acc_rate_trace);
  out.diag.early_stop_iter       = early_stop_iter;
  out.diag.attempted_total       = attempted_total;
  out.diag.accepted_total        = accepted_total;

  // concise summary (robust)
  const int proposals_from_trace = static_cast<int>(H_history.size());
  const int proposals_print = (attempted_total > 0) ? attempted_total : proposals_from_trace;

  double overall_acc = (proposals_print > 0)
    ? static_cast<double>(accepted_total) / proposals_print
  : std::numeric_limits<double>::quiet_NaN();

  std::cout << "[SA] proposals=" << proposals_print
            << " | accepted=" << accepted_total
            << " | overall_acc=" << overall_acc
            << (early_stop_iter > 0 ? " | early_stop_at=" + std::to_string(early_stop_iter) : "")
            << std::endl;

  return out;
}
