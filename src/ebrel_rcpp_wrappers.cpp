//----------------------- R wrappers for cpp functions -------------------------

#include <Rcpp.h>
#include "ebrel_setup.h"
#include "generate_x_zero.h"
#include "objective_utils.h"
#include "dispersal_utils.h"
#include "optimisation_utils.h"
#include "update_candidate.h"
#include "simulated_annealing.h"
#include "sa_utils.h"
#include "run_ebrel.h"
#include "ebrel_builder.h"

#include <algorithm>
#include <cctype>
#include <numeric>
#include <cmath>
#include <limits>
#include <vector>
#include <cstdint>
#include <queue>
#include <string>
#include <cstring>

using namespace Rcpp;

//-------------------------- Local helpers -------------------------------------

namespace {

// ---- helpers for safe option reads ----
inline int opt_int(const Rcpp::List& opt, const char* name, int def) {
  return opt.containsElementNamed(name) ? Rcpp::as<int>(opt[name]) : def;
}

inline double opt_double(const Rcpp::List& opt, const char* name, double def) {
  if (!opt.containsElementNamed(name)) return def;
  Rcpp::NumericVector v = opt[name];
  if (v.size() == 0 || Rcpp::NumericVector::is_na(v[0])) return def;
  return static_cast<double>(v[0]);
}

inline bool opt_bool(const Rcpp::List& opt, const char* name, bool def) {
  return opt.containsElementNamed(name) ? Rcpp::as<bool>(opt[name]) : def;
}

}

//------------------------------- Main functions -------------------------------


// ---- Tuning helper ----
struct TempTuneOptions {
  bool   enabled = false;
  int    num_samples = 400;
  double chi0 = 0.8;
  double p = 2.0;
  double tol_logchi = 1e-3;
  int    max_iters = 50;
  double T1 = -1.0;              // <=0 => heuristic
  int    max_tries_factor = 50;
  int    rng_seed = -1;
  bool   verbose = false;
};

// ---- Parse temperature tuning options
inline TempTuneOptions parse_temp_tune_opts(const Rcpp::List& opt) {
  TempTuneOptions t;

  if (opt.containsElementNamed("tune_temp"))
    t.enabled = Rcpp::as<bool>(opt["tune_temp"]);
  if (!t.enabled) return t;

  if (opt.containsElementNamed("num_samples")) t.num_samples = Rcpp::as<int>(opt["num_samples"]);
  if (opt.containsElementNamed("chi0"))        t.chi0 = Rcpp::as<double>(opt["chi0"]);
  if (opt.containsElementNamed("p"))           t.p = Rcpp::as<double>(opt["p"]);
  if (opt.containsElementNamed("tol_logchi"))  t.tol_logchi = Rcpp::as<double>(opt["tol_logchi"]);
  if (opt.containsElementNamed("max_iters"))   t.max_iters = Rcpp::as<int>(opt["max_iters"]);
  if (opt.containsElementNamed("max_tries_factor"))
    t.max_tries_factor = Rcpp::as<int>(opt["max_tries_factor"]);

  // T1 can be NA in R; treat NA as "not provided"
  if (opt.containsElementNamed("T1")) {
    Rcpp::NumericVector T1r = opt["T1"];
    if (T1r.size() > 0 && Rcpp::NumericVector::is_na(T1r[0]) == false) {
      t.T1 = static_cast<double>(T1r[0]);
    }
  }

  // Prefer user supplied seed
  if (opt.containsElementNamed("seed")) t.rng_seed = Rcpp::as<int>(opt["seed"]);
  else if (opt.containsElementNamed("rng_seed")) t.rng_seed = Rcpp::as<int>(opt["rng_seed"]);

  if (opt.containsElementNamed("verbose")) t.verbose = Rcpp::as<bool>(opt["verbose"]);

  return t;
}

// [[Rcpp::export]]
Rcpp::List run_ebrel_cpp(
    const std::vector<int>& E_int,
    const std::vector<double>& C,
    const std::vector<double>& SD,
    const std::vector<int>& D,
    const std::vector<double>& SxH,
    const std::vector<double>& O,
    const std::vector<int>& LM_int,
    int dim_x,
    int dim_y,
    int n_h,
    int n_s,
    double sentinel,
    Rcpp::List opt
) {
  // ---- cast compact types ----
  std::vector<int8_t>  E(E_int.begin(), E_int.end());
  std::vector<uint8_t> LM(LM_int.begin(), LM_int.end());

  // ---- For build ebrel options ----
  const int    universal_disp_thres = opt_int(opt, "universal_disp_thres", 20);
  const int    max_disp_steps       = opt_int(opt, "max_disp_steps", 10);
  const int    roi_cap              = opt_int(opt, "roi_cap", 100);
  const bool   precompute_W         = opt_bool(opt, "precompute_W", true);

  // ---- Sigma_in: NA / missing => compute default in C++ ----
  const double sigma_in = opt_double(opt, "sigma_in", -1);

  // --- build fully-initialized input + defaults ---
  auto built = build_ebrel(
    std::move(E),
    std::vector<double>(C.begin(), C.end()),
    std::vector<double>(SD.begin(), SD.end()),
    std::vector<int>(D.begin(), D.end()),
    std::vector<double>(SxH.begin(), SxH.end()),
    std::vector<double>(O.begin(), O.end()),
    std::move(LM),
    dim_x, dim_y, n_h, n_s,
    sentinel,
    sigma_in,
    universal_disp_thres,
    max_disp_steps,
    roi_cap,
    precompute_W
  );

  // ---- Parse temperature tuning options ----
  TempTuneOptions tune = parse_temp_tune_opts(opt);

  // ---- For generating X0, if not supplied ----
  if (opt.containsElementNamed("base_prob_X0"))
    built.opt.base_prob_X0 = Rcpp::as<double>(opt["base_prob_X0"]);

  // ---- Choose / generate ONE X0 ----
  std::vector<int8_t> X0_seed;

  if (opt.containsElementNamed("X0")) {
    Rcpp::IntegerVector X0r = opt["X0"];
    X0_seed.assign(X0r.begin(), X0r.end());
  } else {
    // Choose a seed deterministically.
    int seed_for_X0 = opt.containsElementNamed("rng_seed")
    ? Rcpp::as<int>(opt["rng_seed"])
      : tune.rng_seed;

    X0_seed = generate_X0_A(
      built.in.U, built.in.n_h, built.in.dim_x, built.in.dim_y,
      built.opt.base_prob_X0,
      seed_for_X0
    );
  }
  built.in.X0 = X0_seed;

  // ---- Optional tuning of initial temperature ----
  Rcpp::List tune_diag; // empty by default
  if (tune.enabled) {

    // ---- Compute scaling consistent with run_ebrel ----
    const double F1 = compute_F1(built.in.X0, built.in.C, built.in.n_h, built.in.dim_x, built.in.dim_y);
    const double F2 = compute_F2(built.in.X0, built.in.Etiles_per_h, built.in.n_h, built.in.dim_x, built.in.dim_y);

    std::vector<double> G_init = compute_G(
      built.in.X0,
      built.in.LM,
      built.in.row_first_land, built.in.row_last_land,
      built.in.col_first_land, built.in.col_last_land,
      built.in.cell_r, built.in.cell_c,
      built.in.E,
      built.in.dim_x, built.in.dim_y,
      built.in.universal_disp_thres, built.in.max_disp_steps, built.in.roi_cap,
      built.in.rowruns_cache,
      built.in.species_plan
    );
    const double G_sum = std::accumulate(G_init.begin(), G_init.end(), 0.0);

    const double eps = 1e-12;
    const double alpha_scaled = built.in.alpha;
    const double beta_scaled  = (std::abs(F2) > eps)   ? (built.in.beta  * F1 / F2)   : built.in.beta;
    const double gamma_scaled = (std::abs(G_sum) > eps)? (built.in.gamma * F1 / G_sum): built.in.gamma;

    // ---- estimate T0 ----
    InitTempResult t0 = estimate_initial_temperature_benameur_cpp(
      built.in.X0,
      built.W,
      built.in.U,
      built.in.C,
      built.in.O,
      built.in.SxH,
      built.in.D,
      built.in.E,
      built.in.Etiles_per_h,
      built.in.cell_r,
      built.in.cell_c,
      built.in.rowruns_cache,
      built.in.species_plan,
      built.in.universal_disp_thres,
      built.in.max_disp_steps,
      built.in.roi_cap,
      built.in.LM,
      built.in.row_first_land,
      built.in.row_last_land,
      built.in.col_first_land,
      built.in.col_last_land,
      built.in.n_h,
      built.in.n_s,
      built.in.dim_x,
      built.in.dim_y,
      alpha_scaled,
      beta_scaled,
      gamma_scaled,
      built.opt.step_proportion,
      built.opt.step_probability,
      tune.num_samples,
      tune.chi0,
      tune.p,
      tune.tol_logchi,
      tune.max_iters,
      tune.T1,
      tune.max_tries_factor,
      tune.rng_seed,
      tune.verbose
    );

    built.opt.temp = t0.T0;

    tune_diag = Rcpp::List::create(
      Rcpp::Named("enabled")        = true,
      Rcpp::Named("T0")             = t0.T0,
      Rcpp::Named("chi_hat_final")  = t0.diag.chi_hat_final,
      Rcpp::Named("chi0")           = t0.diag.chi0,
      Rcpp::Named("iters")          = t0.diag.iters,
      Rcpp::Named("samples_used")   = t0.diag.samples_used,
      Rcpp::Named("tries")          = t0.diag.tries,
      Rcpp::Named("dE_median")      = t0.diag.dE_median,
      Rcpp::Named("dE_mean")        = t0.diag.dE_mean
    );
  } else {
    tune_diag = Rcpp::List::create(Rcpp::Named("enabled") = false);
  }

  // ---- Apply SA/run options from R to built.opt (override defaults) ----
  // (You can factor this into a helper.)
  if (opt.containsElementNamed("base_prob_X0"))     built.opt.base_prob_X0     = Rcpp::as<double>(opt["base_prob_X0"]);
  if (opt.containsElementNamed("step_proportion"))  built.opt.step_proportion  = Rcpp::as<double>(opt["step_proportion"]);
  if (opt.containsElementNamed("step_probability")) built.opt.step_probability = Rcpp::as<double>(opt["step_probability"]);
  if (opt.containsElementNamed("n_iterations"))     built.opt.n_iterations     = Rcpp::as<int>(opt["n_iterations"]);
  if (opt.containsElementNamed("temp"))             built.opt.temp             = Rcpp::as<double>(opt["temp"]);
  if (opt.containsElementNamed("cooling_rate_c"))   built.opt.cooling_rate_c   = Rcpp::as<double>(opt["cooling_rate_c"]);
  if (opt.containsElementNamed("lam_enabled"))      built.opt.lam_enabled      = Rcpp::as<bool>(opt["lam_enabled"]);
  if (opt.containsElementNamed("lam_target_mid"))   built.opt.lam_target_mid   = Rcpp::as<double>(opt["lam_target_mid"]);
  if (opt.containsElementNamed("lam_target_final")) built.opt.lam_target_final = Rcpp::as<double>(opt["lam_target_final"]);
  if (opt.containsElementNamed("lam_hold_frac"))    built.opt.lam_hold_frac    = Rcpp::as<double>(opt["lam_hold_frac"]);
  if (opt.containsElementNamed("lam_p"))            built.opt.lam_p            = Rcpp::as<double>(opt["lam_p"]);
  if (opt.containsElementNamed("min_iterations"))   built.opt.min_iterations   = Rcpp::as<int>(opt["min_iterations"]);
  if (opt.containsElementNamed("acceptance_window"))built.opt.acceptance_window= Rcpp::as<int>(opt["acceptance_window"]);
  if (opt.containsElementNamed("acceptance_thres")) built.opt.acceptance_thres = Rcpp::as<double>(opt["acceptance_thres"]);
  if (opt.containsElementNamed("iter_no_improve"))  built.opt.iter_no_improve  = Rcpp::as<int>(opt["iter_no_improve"]);
  if (opt.containsElementNamed("improve_eps"))      built.opt.improve_eps      = Rcpp::as<double>(opt["improve_eps"]);
  if (opt.containsElementNamed("rng_seed"))         built.opt.rng_seed         = Rcpp::as<int>(opt["rng_seed"]);
  if (opt.containsElementNamed("write_every"))      built.opt.write_every      = Rcpp::as<int>(opt["write_every"]);
  if (opt.containsElementNamed("trace_file"))       built.opt.trace_file       = Rcpp::as<std::string>(opt["trace_file"]);
  if (opt.containsElementNamed("verbose"))          built.opt.verbose          = Rcpp::as<bool>(opt["verbose"]);

  // ---- User supplied alpha, beta, gamma ----
  if (opt.containsElementNamed("alpha")) built.in.alpha = Rcpp::as<double>(opt["alpha"]);
  if (opt.containsElementNamed("beta"))  built.in.beta  = Rcpp::as<double>(opt["beta"]);
  if (opt.containsElementNamed("gamma")) built.in.gamma = Rcpp::as<double>(opt["gamma"]);

  // ---- Run Ebrel----
  RunEBRELResult res = run_ebrel(built.in, built.opt);

  // ---- Convert to integer vector ----
  Rcpp::IntegerVector X_best_iv(res.X_best.size());
  for (std::size_t i = 0; i < res.X_best.size(); ++i) {
    X_best_iv[i] = static_cast<int>(res.X_best[i]);
  }

  // ---- X0
  Rcpp::IntegerVector X0(X0_seed.size());
  for (std::size_t i = 0; i < X0_seed.size(); ++i) {
    X0[i] = static_cast<int>(X0_seed[i]);
  }

  const int cells = built.in.dim_x * built.in.dim_y;
  const int nh    = built.in.n_h;

  // built.in.U is std::vector<uint8_t> of length cells*nh
  Rcpp::IntegerVector U_iv(built.in.U.size());
  for (std::size_t i = 0; i < built.in.U.size(); ++i) {
    U_iv[i] = static_cast<int>(built.in.U[i]);  // 0/1 (or 0..255 if not binary)
  }
  U_iv.attr("dim") = Rcpp::IntegerVector::create(cells, nh);

  // ---- Return ----
  return Rcpp::List::create(
    Rcpp::Named("X0") = X0,
    Rcpp::Named("U") = U_iv,
    Rcpp::Named("X_best") = X_best_iv,
    Rcpp::Named("H_best") = res.H_best,
    Rcpp::Named("iterations_run") = res.iterations_run,
    Rcpp::Named("accepted") = res.accepted,
    Rcpp::Named("proposals") = res.proposals,
    Rcpp::Named("overall_acc") = res.overall_acc,
    Rcpp::Named("early_stop_iter") = res.early_stop_iter,
    Rcpp::Named("H_trace") = res.H_trace,
    Rcpp::Named("F_trace") = res.F_trace,
    Rcpp::Named("F1_trace") = res.F1_trace,
    Rcpp::Named("F2_trace") = res.F2_trace,
    Rcpp::Named("g_best") = res.g_best,
    Rcpp::Named("acc_rate_trace") = res.acc_rate_trace,
    Rcpp::Named("iter_ms_total") = res.iter_ms_total,
    Rcpp::Named("iter_count") = res.iter_count
  );

}

// [[Rcpp::export]]
Rcpp::IntegerVector generate_X0_A_R(const std::vector<uint8_t>& U,
                                    int n_h,
                                    int dim_x,
                                    int dim_y,
                                    double base_prob,
                                    int seed) {
  std::vector<int8_t> x0 = generate_X0_A(U, n_h, dim_x, dim_y, base_prob, seed);

  Rcpp::IntegerVector out(x0.size());
  for (std::size_t i = 0; i < x0.size(); ++i) out[i] = static_cast<int>(x0[i]);

  return out;
}
