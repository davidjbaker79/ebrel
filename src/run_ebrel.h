//-------------------------- Run Ebrel model -----------------------------------

#ifndef RUN_EBREL_H
#define RUN_EBREL_H

#include <vector>
#include <stdexcept>
#include <limits>

//---------------------------- Main function -----------------------------------

// User supplied options for running ebrel model
struct RunEBRELOptions {
  int max_disp_thres = 50;
  int disp_boundary = 30;
  double sigma            = 0.05;
  const std::vector<double>* X0 = nullptr;   // X0 supply optional
  double base_prob_X0     = 0.85;
  double step_proportion  = 0.05;
  double step_probability = 0.15;
  int n_iterations        = 10000;
  double temp             = 2000;
  double cooling_rate_c   = 1;
  bool   lam_enabled      = false;  // turn Lam-style online control on/off
  double lam_target_mid   = 0.44;   // target uphill acceptance during early/mid run
  double lam_target_final = 0.05;   // target uphill acceptance at the end
  double lam_hold_frac    = 0.60;   // fraction of run to hold lam_target_mid before decaying
  double lam_p            = 2.0;     // damping exponent in Ben-Ameur correction (>=1)
  int min_iterations      = 1000;  // require at least this many iterations
  int acceptance_window   = 1000;  // window length for acceptance rate
  double acceptance_thres = 0.01;  // "low" acceptance threshold
  int iter_no_improve     = 1000;  // consecutive iters with no meaningful improvement
  double improve_eps      = 1e-6;  // relative improvement needed to reset patience
  int rng_seed            = -1;
  bool verbose            = false;
};

// Input data structure
struct RunEBRELInput {
  int dim_x = 0, dim_y = 0, n_h = 0, n_s = 0, max_disp_thres = 0, disp_boundary = 0;
  double alpha = 1.0, beta = 25.0, gamma = 100.0;

  std::vector<double> U, C, E, O, SD, SxH;
  std::vector<int>    D;
  std::vector<double> X0; // optional
};

// Results structure
struct RunEBRELResult {
  std::vector<double> X_best;
  double H_best = 0.0;
  std::vector<double> H_trace, F_trace, F1_trace, F2_trace, g_best;
  int iterations_run = 0;
  // --- Diagnostics from SA ---
  std::vector<double> acc_rate_trace;  // acceptance per window
  int    early_stop_iter = -1;         // iteration where combined stop fired, or -1
  int    proposals = 0;                // total evaluated proposals (attempted_total)
  int    accepted  = 0;                // total accepted proposals
  double overall_acc = std::numeric_limits<double>::quiet_NaN(); // accepted/proposals
};

// Declare the wrapper
RunEBRELResult run_ebrel(const RunEBRELInput& in, const RunEBRELOptions& opt);

#endif // RUN_EBREL_H
