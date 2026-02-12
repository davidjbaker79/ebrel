// -------------------------------Simulated annealing----------------------------

#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <string>

#include "offsets.h"
#include "species_plan.h"

// ------------------------------- Main functions ------------------------------

// For SA diagnostics
struct SADiagnostics {
  std::vector<double> acceptance_rate_trace; // acceptance per 'acceptance_window'
  int                 early_stop_iter = -1;  // iteration (1-based) is stopped early, -1 otherwise
  int                 attempted_total = 0;   // total proposals attempted (non no-op)
  int                 accepted_total  = 0;   // total accepted
  long long           iter_ms_total   = 0;   // total iterations time in ms
  int                 iter_count      = 0;   // number of iterations
};

// Result container for simulated annealing (Rcpp-free)
struct SAResult {
  std::vector<int8_t> X_best;   // flattened [n_cells]
  double H_best;                // best objective
  std::vector<double> H_trace;  // candidate H evaluations
  std::vector<double> F_trace;  // optional Fx trace from compute_H
  std::vector<double> F1_trace; // optional F1 trace
  std::vector<double> F2_trace; // optional F2 trace
  std::vector<double> g_best;   // Short-fall in targets
  SADiagnostics diag;
};

// Simulated annealing
SAResult simulated_annealing(
    const std::vector<int8_t>& X0,              // size: n_cells
    const std::vector<double>& W,               // size: n_cells * n_h
    const std::vector<uint8_t>& U,              // size: n_cells * n_h
    const std::vector<double>& C,               // per compute_H
    const std::vector<double>& O,               // per compute_H
    const std::vector<double>& SxH,             // flattened matrix for compute_H
    const std::vector<int>& D,                  // dispersal thresholds, etc.
    const std::vector<int8_t>& E,               // -1,0,1,...,n_h-1 encoded E
    const std::vector<std::vector<std::size_t>>& Etiles_per_h, // lookup Etiles
    const std::vector<int>& cell_r,             // rows indexes
    const std::vector<int>& cell_c,             // cols indexes
    const RowRunsCache& rowruns_cache,          // Caches rowruns
    const SpeciesPlan& species_plan,            // Species dispersal info e.g. ROI
    int n_h,                                    // number of habitats/classes per cell
    int n_s,                                    // number of species/features
    int dim_x,                                  // X dimension (e.g. longitude or easting)
    int dim_y,                                  // Y dimension (e.g. latitude or northing)
    int universal_disp_thres,                   // threshold after which universal dispersal is assumed
    int max_disp_steps,                         // limit for dispersal as D[sp] x universal_disp_thres
    int roi_cap,                                // ROI cap in cells (upper limit for D[sp] x universal_disp_thres)
    const std::vector<uint8_t>& LM,             // Land mask
    const std::vector<int16_t>& row_first_land, // Vector of first land index per row
    const std::vector<int16_t>& row_last_land,  // vector of last land index per row
    const std::vector<int16_t>& col_first_land, // vector of first land index per col
    const std::vector<int16_t>& col_last_land,  // vector of last land index per col
    double alpha_scaled = 1,                    // scaling on targets
    double beta_scaled = 25,                    // scaling on spatial aggregation
    double gamma_scaled = 100,                  // scaling on cost
    double step_proportion = 0.05,              // proportion of eligible cells to update
    double step_probability = 0.05,             // probability of assigning any habitat
    int n_iterations = 10000,                   // max number of iterations for SA
    double temp = 5000,                         // initial temperature
    double cooling_rate_c = 1.0,                // control parameter for cooling rate (used only if lam is false)
    bool   lam_enabled = false,                 // turn Lam-style online control on/off
    double lam_target_mid = 0.44,               // target uphill acceptance during early/mid run
    double lam_target_final = 0.05,             // target uphill acceptance at the end
    double lam_hold_frac = 0.60,                // fraction of run to hold lam_target_mid before decaying
    double lam_p = 2.0,                         // damping exponent in Ben-Ameur correction (>=1)
    int min_iterations = 1000,                  // require at least this many iterations
    int acceptance_window = 200,                // window length for acceptance rate
    double acceptance_thres = 0.001,            // "low" acceptance threshold
    int iter_no_improve = 1000,                 // consecutive iters with no meaningful improvement
    double improve_eps = 1e-6,                  // relative improvement needed to reset patience
    int  write_every = 0,                       // 0 = disabled
    std::string trace_file = "",                // empty = disabled
    bool verbose = false
);

#endif // SIMULATED_ANNEALING_H
