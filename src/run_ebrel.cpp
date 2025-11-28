//-------------------------- Run Ebrel model -----------------------------------

#include "run_ebrel.h"
#include "simulated_annealing.h"
#include "objective_utils.h"
#include "dispersal_utils.h"
#include "generate_x_zero.h"
#include "optimisation_utils.h"

#include <numeric>
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdint>

//---------------------- Local helper functions --------------------------------

namespace {

  inline void ensure(bool ok, const char* msg) {
    if (!ok) throw std::invalid_argument(msg);
  }

  inline size_t n_cells(int dim_x, int dim_y) {
    return static_cast<size_t>(dim_x) * static_cast<size_t>(dim_y);
  }

  inline void validate_shapes(const RunEBRELInput& in) {
    const size_t cells = n_cells(in.dim_x, in.dim_y);
    const size_t Hsz   = cells * static_cast<size_t>(in.n_h);
    const size_t Ssz   = cells * static_cast<size_t>(in.n_s);
    const size_t SxHsz = static_cast<size_t>(in.n_h) * static_cast<size_t>(in.n_s);

    ensure(in.dim_x > 0 && in.dim_y > 0, "dim_x and dim_y must be positive");
    ensure(in.n_h   > 0, "n_h must be positive");
    ensure(in.n_s   > 0, "n_s must be positive");

    ensure(in.U.size()   == Hsz,   "U size mismatch (expected dim_x*dim_y*n_h)");
    ensure(in.C.size()   == Hsz,   "C size mismatch (expected dim_x*dim_y*n_h)");
    ensure(in.E.size()   == Hsz,   "E size mismatch (expected dim_x*dim_y*n_h)");
    ensure(in.O.size()   == static_cast<size_t>(in.n_s), "O size mismatch (expected n_s)");
    ensure(in.SD.size()  == Ssz,   "SD size mismatch (expected dim_x*dim_y*n_s)");
    ensure(in.SxH.size() == SxHsz, "SxH size mismatch (expected n_h*n_s)");
    ensure(in.D.size()   == static_cast<size_t>(in.n_s), "D size mismatch (expected n_s)");

    if (!in.X0.empty()) {
      ensure(in.X0.size() == Hsz, "X0 size mismatch (expected dim_x*dim_y*n_h)");
    }

  }

}

//--------------------- Main functions -----------------------------------------

RunEBRELResult run_ebrel(const RunEBRELInput& in, const RunEBRELOptions& opt) {
  validate_shapes(in);

  // --- X0: provided or generated ---
  std::vector<double> X0;
  if (!in.X0.empty()) {
    X0 = in.X0; // user-supplied
  } else {
    X0 = generate_X0_A(in.U, in.n_h, in.dim_x, in.dim_y,
                         opt.base_prob_X0, opt.rng_seed);
  };

  // --- Initial objectives (for scaling) ---
  const double F1 = compute_F1(X0, in.C, in.n_h, in.dim_x, in.dim_y);
  const double F2 = compute_F2(X0, in.E, in.n_h, in.dim_x, in.dim_y);

  // Use the same G variant as SA so scaling matches the run
  std::vector<double> G_init;
  G_init = compute_G(
    X0,
    in.E,
    in.SD,
    in.SxH,
    in.D,
    in.n_h,
    in.n_s,
    in.dim_x,
    in.dim_y,
    in.universal_disp_thres,
    in.max_disp_steps,
    in.roi_cap,
    in.LM,
    in.row_first_land,
    in.row_last_land,
    in.col_first_land,
    in.col_last_land
  );
  const double G_sum = std::accumulate(G_init.begin(), G_init.end(), 0.0);

  // Safe rescaling
  const double eps = 1e-12;
  const double alpha_scaled = in.alpha;
  const double beta_scaled  = (std::abs(F2) > eps) ? (in.beta  * F1 / F2)  : in.beta;
  const double gamma_scaled = (std::abs(G_sum) > eps) ? (in.gamma * F1 / G_sum) : in.gamma;

  // Distance weights
  std::vector<double> W = compute_distance_weights(in.E, in.U, in.n_h, in.dim_x, in.dim_y, opt.sigma);
  {
    const size_t Hsz = n_cells(in.dim_x, in.dim_y) * static_cast<size_t>(in.n_h);
    ensure(W.size() == Hsz, "W size mismatch after compute_distance_weights");
  }

  if (opt.verbose) {
    std::cout << "[run_ebrel] Starting SA with "
              << in.dim_x << "x" << in.dim_y
              << ", n_h=" << in.n_h << ", n_s=" << in.n_s
              << ", iterations=" << opt.n_iterations << std::endl;
  }

  SAResult sa = simulated_annealing(
    std::move(X0),
    W,
    in.U,
    in.C,
    in.E,
    in.O,
    in.SD,
    in.SxH,
    in.D,
    in.n_h,
    in.n_s,
    in.dim_x,
    in.dim_y,
    in.universal_disp_thres,
    in.max_disp_steps,
    in.roi_cap,
    in.LM,
    in.row_first_land,
    in.row_last_land,
    in.col_first_land,
    in.col_last_land,
    alpha_scaled,
    beta_scaled,
    gamma_scaled,
    opt.step_proportion,
    opt.step_probability,
    opt.n_iterations,
    opt.temp,
    opt.cooling_rate_c,
    opt.lam_enabled,
    opt.lam_target_mid,
    opt.lam_target_final,
    opt.lam_hold_frac,
    opt.lam_p,
    opt.min_iterations,
    opt.acceptance_window,
    opt.acceptance_thres,
    opt.iter_no_improve,
    opt.improve_eps,
    opt.verbose
  );

  RunEBRELResult out;
  out.X_best         = std::move(sa.X_best);
  out.H_best         = sa.H_best;
  out.H_trace        = std::move(sa.H_trace);
  out.F_trace        = std::move(sa.F_trace);
  out.F1_trace       = std::move(sa.F1_trace);
  out.F2_trace       = std::move(sa.F2_trace);
  out.iterations_run = static_cast<int>(out.H_trace.size());
  out.g_best         = std::move(sa.g_best);

  // --- NEW: diagnostics ---
  out.acc_rate_trace   = std::move(sa.diag.acceptance_rate_trace);
  out.early_stop_iter  = sa.diag.early_stop_iter;

  // Prefer SA’s attempted_total; fall back to trace length if zero
  out.proposals = (sa.diag.attempted_total > 0)
    ? sa.diag.attempted_total
  : static_cast<int>(out.H_trace.size());
  out.accepted  = sa.diag.accepted_total;
  out.overall_acc = (out.proposals > 0)
    ? static_cast<double>(out.accepted) / out.proposals
  : std::numeric_limits<double>::quiet_NaN();

  if (opt.verbose) {
    std::cout << "[run_ebrel] Done. Best H = " << out.H_best
              << " (cand evals recorded: " << out.iterations_run << ")"
              << std::endl;
  }
  return out;
}
