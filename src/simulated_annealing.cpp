//------------------------------ Simulated annealing ---------------------------

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
#include <cstdint>

/* Simulated annealing algorithm with:
 - option for Lam et al style tuning schedule [Lam, J. and Delosme, J.M., 1988.
 An efficient simulated annealing schedule: derivation. New Haven, CT: Yale
 Electrical Engineering Department, 8816.]
 - or harmonic/geometric cooling schedule with a option to control rate of cooling
 */
SAResult simulated_annealing(
    // --------- Input data and parameters
    std::vector<double> X0,
    const std::vector<double>& W,
    const std::vector<double>& U,
    const std::vector<double>& C,
    const std::vector<double>& O,
    const std::vector<double>& SxH,
    const std::vector<int>& D,
    const std::vector<int>& E_h_of_cell,
    const std::vector<std::vector<std::size_t>>& Etiles_per_h,
    const std::vector<int>&    cell_r,
    const std::vector<int>&    cell_c,
    const std::vector<SpeciesDispData>& species_info,
    // ---------- Dimensions
    int n_h,
    int n_s,
    int dim_x,
    int dim_y,
    int universal_disp_thres,
    int max_disp_steps,
    // ---------- Cap for size of ROI (in grid cells)
    int roi_cap,
    // ---------- For avoiding extra work on sea cells
    const std::vector<uint8_t>& LM,
    const std::vector<int>& row_first_land,
    const std::vector<int>& row_last_land,
    const std::vector<int>& col_first_land,
    const std::vector<int>& col_last_land,
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

  // Initial evaluation of X0
  HResult init_scores = compute_H(X0, C, O, SxH, D,
                                  alpha_scaled, beta_scaled, gamma_scaled,
                                  n_h, n_s,
                                  dim_x, dim_y,
                                  universal_disp_thres, max_disp_steps, roi_cap,
                                  LM,
                                  row_first_land, row_last_land,
                                  col_first_land, col_last_land,
                                  E_h_of_cell,
                                  Etiles_per_h,
                                  cell_r, cell_c,
                                  species_info);

  double best_score = init_scores.H;

  // Current & best states
  std::vector<double> best = X0;
  std::vector<double> curr = X0;
  double curr_eval = best_score;
  std::vector<double> g_best = init_scores.g;

  // Ensure stopping parameters are in a sensible range
  const int    min_iters    = std::max(0, min_iterations);
  const int    win          = std::max(1, acceptance_window);
  const double acc_min      = std::clamp(acceptance_thres, 0.0, 1.0);
  const int    patience     = std::max(1, iter_no_improve);
  const double improve_thr  = std::max(0.0, improve_eps);

  // Diagnostics
  std::vector<double> acc_rate_trace; // per window
  int    early_stop_iter = -1;
  int    attempted_total = 0;
  int    accepted_total  = 0;

  // For stopping criterion
  int    attempted_in_win = 0;
  int    accepted_in_win  = 0;
  int    no_improve       = 0;
  double last_best        = best_score;

  // ---------- LAM: helpers & counters ----------
  auto clamp01 = [](double x, double eps = 1e-12) {
    if (x < eps) return eps;
    if (x > 1.0 - eps) return 1.0 - eps;
    return x;
  };
  auto target_acceptance = [&](double progress01) {
    progress01 = std::min(1.0, std::max(0.0, progress01));
    if (progress01 <= lam_hold_frac) return lam_target_mid;
    double u = (progress01 - lam_hold_frac) / (1.0 - lam_hold_frac);
    return std::exp((1.0 - u) * std::log(lam_target_mid) + u * std::log(lam_target_final));
  };
  int uphill_attempted_in_win = 0;
  int uphill_accepted_in_win  = 0;
  double temp_curr = temp; // for LAM

  // Main loop
  for (int z = 0; z < n_iterations; ++z) {

    // Propose candidate
    std::vector<double> candidate = update_candidate(W, U, curr,
                                                     step_proportion,
                                                     step_probability,
                                                     n_h, dim_x, dim_y);

    if (candidate == curr) {
      continue; // skip no-op proposals
    }

    attempted_in_win++;
    attempted_total++;

    // Evaluate candidate with new params
    HResult scores = compute_H(candidate, C, O, SxH, D,
                               alpha_scaled, beta_scaled, gamma_scaled,
                               n_h, n_s,
                               dim_x, dim_y,
                               universal_disp_thres, max_disp_steps, roi_cap,
                               LM,
                               row_first_land, row_last_land,
                               col_first_land, col_last_land,
                               E_h_of_cell,
                               Etiles_per_h,
                               cell_r, cell_c,
                               species_info);
    double candidate_eval = scores.H;

    // Trace histories
    H_history.push_back(candidate_eval);
    F_history.push_back(scores.Fx);
    F1_history.push_back(scores.F1);
    F2_history.push_back(scores.F2);

    // Update best
    if (candidate_eval < best_score) {
      best = candidate;
      best_score = candidate_eval;
      g_best = scores.g;
    }

    // ------------- Metropolis acceptance -------------------------------
    const double delta = candidate_eval - curr_eval;
    const double epsE  = 1e-12 * (1.0 + std::max(std::abs(curr_eval), std::abs(candidate_eval)));

    const bool downhill = (delta < -epsE);
    const bool uphill   = (delta >  epsE);
    const bool flat     = !downhill && !uphill;

    if (uphill) uphill_attempted_in_win++;

    const double t = (lam_enabled)
      ? temp_curr
    : (temp / (1.0 + cooling_rate_c * z));  // harmonic default

    bool accepted = false;
    if (downhill || flat) {
      accepted = true;
    } else if (t > 0.0 && std::isfinite(t)) {
      double u = runif(rng);
      if (u <= 0.0) u = std::numeric_limits<double>::min();
      const double logu = std::log(u);
      const double logA = -delta / t;
      accepted = (logu < std::min(0.0, logA));
    }

    if (accepted) {
      curr = candidate;
      curr_eval = candidate_eval;
      accepted_in_win++;
      accepted_total++;
      if (uphill) uphill_accepted_in_win++;
    }

    // ---- Stopping logic ----
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

    if ((z + 1) % win == 0) {

      double acc_rate = (attempted_in_win > 0)
      ? static_cast<double>(accepted_in_win) / attempted_in_win
      : 0.0;

      acc_rate_trace.push_back(acc_rate);

      if ((z + 1) >= min_iters && acc_rate < acc_min && no_improve >= patience) {
        std::cout << "Stopping early at iteration " << (z + 1)
                  << " due to combined criterion: low acceptance AND no improvement.\n";
        early_stop_iter = (z + 1);
        break;
      }

      if (lam_enabled) {
        double hat_up = (uphill_attempted_in_win > 0)
        ? static_cast<double>(uphill_accepted_in_win) / uphill_attempted_in_win
        : ((attempted_in_win > 0) ? static_cast<double>(accepted_in_win) / attempted_in_win : 0.0);

        double progress = static_cast<double>(z + 1) / static_cast<double>(n_iterations);
        double a_star   = target_acceptance(progress);

        double num = std::log(clamp01(hat_up));
        double den = std::log(clamp01(a_star));
        double ratio = (den != 0.0) ? (num / den) : 1.0;
        if (!(ratio > 0.0) || !std::isfinite(ratio)) ratio = 1.0;

        temp_curr *= std::pow(ratio, 1.0 / lam_p);
        if (!std::isfinite(temp_curr) || temp_curr < 1e-9) temp_curr = 1e-9;
        if (temp_curr > 1e16) temp_curr = 1e16;
      }

      // reset window counters
      attempted_in_win = 0;
      accepted_in_win  = 0;
      uphill_attempted_in_win = 0;
      uphill_accepted_in_win  = 0;
    }

    if (verbose && ((z + 1) % win == 0)) {
      std::cout << "[SA] it=" << (z + 1)
                << " T=" << (lam_enabled ? temp_curr : (temp / (1.0 + cooling_rate_c * z)))
                << " Best H(x) = " << best_score
                << " Best F(x) = " << scores.Fx
                << '\n' << std::flush;
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
  out.g_best   = g_best;

  // Diagnostics
  out.diag.acceptance_rate_trace = std::move(acc_rate_trace);
  out.diag.early_stop_iter       = early_stop_iter;
  out.diag.attempted_total       = attempted_total;
  out.diag.accepted_total        = accepted_total;

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
