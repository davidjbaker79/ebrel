// -------------------------------Simulated annealing----------------------------

#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include <vector>
#include <stdexcept>


// For SA diagnostics
struct SADiagnostics {
  std::vector<double> acceptance_rate_trace; // acceptance per 'acceptance_window'
  int                 early_stop_iter = -1;  // iteration (1-based) is stopped early, -1 otherwise
  int                 attempted_total = 0;   // total proposals attempted (non no-op)
  int                 accepted_total  = 0;   // total accepted
};

// Result container for simulated annealing (Rcpp-free)
struct SAResult {
  std::vector<double> X_best;   // flattened [n_cells * n_h]
  double H_best;                // best objective
  std::vector<double> H_trace;  // candidate H evaluations
  std::vector<double> F_trace;  // optional Fx trace from compute_H
  std::vector<double> F1_trace; // optional F1 trace
  std::vector<double> F2_trace; // optional F2 trace
  std::vector<double> g_best; // Short-fall in targets
  SADiagnostics diag;
};

// Simulated annealing
SAResult simulated_annealing(
    std::vector<double> X0,           // size: n_cells * n_h
    const std::vector<double>& W,     // size: n_cells * n_h (or your chosen layout)
    const std::vector<double>& U,     // size: n_cells * n_h
    const std::vector<double>& C,     // per compute_H
    const std::vector<double>& E,     // per compute_H
    const std::vector<double>& O,     // per compute_H
    const std::vector<double>& SD,    // per compute_H
    const std::vector<double>& SxH,   // flattened matrix for compute_H
    const std::vector<int>& D,        // dispersal thresholds, etc.
    int max_disp_thres,               // threshold after which universal dispersal is assumed
    int disp_boundary,                // limit for dispersal as D[sp] x disp_boundary
    int n_h,                          // number of habitats/classes per cell
    int n_s,                          // number of species/features
    int dim_x,                        // X dimension (e.g. longitude or easting)
    int dim_y,                        // Y dimension (e.g. latitude or northing)
    double alpha_scaled,              // scaling on targets
    double beta_scaled,               // scaling on spatial aggregation
    double gamma_scaled,              // scaling on costs
    double step_proportion,           // proportion of eligible cells to update
    double step_probability,          // probability of assigning any habitat
    int n_iterations,                 // max number of iterations for SA
    double temp,                      // initial temperature
    double cooling_rate_c,            // control parameter for cooling rate (used only if lam is false)
    bool   lam_enabled,               // turn Lam-style online control on/off
    double lam_target_mid,            // target uphill acceptance during early/mid run
    double lam_target_final,          // target uphill acceptance at the end
    double lam_hold_frac,             // fraction of run to hold lam_target_mid before decaying
    double lam_p,                     // damping exponent in Ben-Ameur correction (>=1)
    int min_iterations,               // require at least this many iterations
    int acceptance_window,            // window length for acceptance rate
    double acceptance_thres,          // "low" acceptance threshold
    int iter_no_improve,              // consecutive iters with no meaningful improvement
    double improve_eps,               // relative improvement needed to reset patience
    bool verbose
);

#endif // SIMULATED_ANNEALING_H
