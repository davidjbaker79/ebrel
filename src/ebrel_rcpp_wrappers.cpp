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

#include <algorithm>
#include <cctype>
#include <numeric>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List create_ebrel_class_object_R(std::vector<double> E,
                                       std::vector<double> C,
                                       std::vector<double> SD,
                                       std::vector<int> D,
                                       std::vector<double> SxH,
                                       std::vector<double> O,
                                       int dim_x,
                                       int dim_y,
                                       int n_h,
                                       int n_s,
                                       double sentinel = 1e10) {
  // Convert Rcpp inputs to STL vectors
  std::vector<double> E_v(E.begin(), E.end());
  std::vector<double> C_v(C.begin(), C.end());
  std::vector<double> SD_v(SD.begin(), SD.end());
  std::vector<int> D_v(D.begin(), D.end());
  std::vector<double> SxH_v(SxH.begin(), SxH.end());
  std::vector<double> O_v(O.begin(), O.end());

  std::vector<double> U_v;

  // Call core C++ function
  create_ebrel_class_object(E_v, C_v, SD_v, D_v, SxH_v, O_v,
                            dim_x, dim_y, n_h, n_s, sentinel, U_v);

  // Return as R list
  return List::create(
    _["E"] = E,
    _["C"] = C,
    _["U"] = NumericVector(U_v.begin(), U_v.end()),
    _["SD"] = SD,
    _["D"] = D,
    _["SxH"] = SxH,
    _["O"] = O,
    _["dim1"] = dim_x,
    _["dim2"] = dim_y,
    _["nh"] = n_h,
    _["ns"] = n_s
  );
}

// [[Rcpp::export]]
Rcpp::List run_ebrel_R(
    Rcpp::List ebrel_obj,
    Rcpp::Nullable<Rcpp::NumericVector> X0 = R_NilValue,
    double base_prob_X0     = 0.85,
    double sigma            = 0.05,
    int max_disp_thres      = 50,
    int disp_boundary       = 30,
    double alpha            = 1.0,
    double beta             = 25.0,
    double gamma            = 100.0,
    double step_proportion  = 0.05,
    double step_probability = 0.05,
    int    n_iterations     = 10000,
    double temp             = 1000.0,
    double cooling_rate_c   = 1.0,
    bool   lam_enabled      = false,
    double lam_target_mid   = 0.44,
    double lam_target_final = 0.05,
    double lam_hold_frac    = 0.60,
    double lam_p            = 2.0,
    int min_iterations      = 1000,
    int acceptance_window   = 1000,
    double acceptance_thres = 0.01,
    int iter_no_improve     = 1000,
    double improve_eps      = 1e-6,
    Rcpp::Nullable<Rcpp::IntegerVector> seed = R_NilValue,
    bool   verbose          = false

) {
  try {
    Rcpp::NumericVector E   = ebrel_obj["E"];
    Rcpp::NumericVector C   = ebrel_obj["C"];
    Rcpp::NumericVector U   = ebrel_obj["U"];
    Rcpp::NumericVector SD  = ebrel_obj["SD"];
    Rcpp::IntegerVector D   = ebrel_obj["D"];
    Rcpp::NumericVector SxH = ebrel_obj["SxH"];
    Rcpp::NumericVector O   = ebrel_obj["O"];   // O = targets

    int dim_x = Rcpp::as<int>(ebrel_obj["dim1"]);
    int dim_y = Rcpp::as<int>(ebrel_obj["dim2"]);
    int n_h   = Rcpp::as<int>(ebrel_obj["nh"]);
    int n_s   = Rcpp::as<int>(ebrel_obj["ns"]);

    // O length should be n_s
    if (O.size() != n_s) Rcpp::stop("`O` length must equal ns.");

    RunEBRELInput in;
    in.dim_x = dim_x; in.dim_y = dim_y;
    in.n_h   = n_h;   in.n_s   = n_s;

    in.alpha = alpha; in.beta  = beta; in.gamma = gamma;

    in.U   = Rcpp::as<std::vector<double>>(U);
    in.C   = Rcpp::as<std::vector<double>>(C);
    in.E   = Rcpp::as<std::vector<double>>(E);
    in.O   = Rcpp::as<std::vector<double>>(O);
    in.SD  = Rcpp::as<std::vector<double>>(SD);
    in.SxH = Rcpp::as<std::vector<double>>(SxH);
    in.D   = Rcpp::as<std::vector<int>>(D);
    in.max_disp_thres = max_disp_thres;
    in.disp_boundary = disp_boundary;

    // X0 - Null (and therefore naive) or user supplied
    const std::size_t expected_len =
      static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y) * static_cast<std::size_t>(n_h);
    if (X0.isNotNull()) {
      Rcpp::NumericVector X0_r = X0.get();
      if (static_cast<std::size_t>(X0_r.size()) != expected_len) {
        Rcpp::stop("`X0` length mismatch: got %d, expected %d",
                   X0_r.size(), static_cast<int>(expected_len));
      }
      in.X0 = Rcpp::as<std::vector<double>>(X0_r);
    }

    RunEBRELOptions opt;
    opt.max_disp_thres   = max_disp_thres;
    opt.step_proportion  = step_proportion;
    opt.step_probability = step_probability;
    opt.n_iterations     = n_iterations;
    opt.temp             = temp;
    opt.cooling_rate_c   = cooling_rate_c;
    opt.lam_enabled      = lam_enabled;
    opt.lam_target_mid   = lam_target_mid;
    opt.lam_target_final = lam_target_final;
    opt.lam_hold_frac    = lam_hold_frac;
    opt.lam_p            = lam_p;
    opt.min_iterations   = min_iterations;
    opt.acceptance_window= acceptance_window;
    opt.acceptance_thres = acceptance_thres;
    opt.iter_no_improve  = iter_no_improve;
    opt.improve_eps      = improve_eps;
    opt.base_prob_X0     = base_prob_X0;
    opt.sigma            = sigma;
    opt.rng_seed         = seed.isNotNull() ? Rcpp::as<int>(seed.get()) : -1;
    opt.verbose          = verbose;


    RunEBRELResult res = run_ebrel(in, opt);

    return Rcpp::List::create(
      Rcpp::Named("X_best")         = Rcpp::NumericVector(res.X_best.begin(), res.X_best.end()),
      Rcpp::Named("H_best")         = res.H_best,
      Rcpp::Named("g_best")         = Rcpp::NumericVector(res.g_best.begin(),  res.g_best.end()),
      Rcpp::Named("H_trace")        = Rcpp::NumericVector(res.H_trace.begin(),  res.H_trace.end()),
      Rcpp::Named("F_trace")        = Rcpp::NumericVector(res.F_trace.begin(),  res.F_trace.end()),
      Rcpp::Named("F1_trace")       = Rcpp::NumericVector(res.F1_trace.begin(), res.F1_trace.end()),
      Rcpp::Named("F2_trace")       = Rcpp::NumericVector(res.F2_trace.begin(), res.F2_trace.end()),
      Rcpp::Named("iterations_run") = res.iterations_run,
      // --- Diagnostics ---
      Rcpp::Named("acc_rate_trace")   = Rcpp::NumericVector(res.acc_rate_trace.begin(),res.acc_rate_trace.end()),
      Rcpp::Named("early_stop_iter")  = res.early_stop_iter,
      Rcpp::Named("proposals")        = res.proposals,
      Rcpp::Named("accepted")         = res.accepted,
      Rcpp::Named("overall_acc")      = res.overall_acc
    );
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("Unknown error in run_ebrel_R");
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector generate_X0_A_R(const std::vector<double>& U,
                                    int n_h,
                                    int dim_x,
                                    int dim_y,
                                    double base_prob,
                                    int seed) {
  std::vector<double> result = generate_X0_A(U, n_h, dim_x, dim_y, base_prob, seed);
  return NumericVector(result.begin(), result.end());
}

// Function to calculate initial temperature based on method proposed in
// Ben-Ameur, Walid. "Computing the initial temperature of simulated annealing."
// Computational Optimization and Applications 29, no. 3 (2004): 369-385.
// [[Rcpp::export]]
Rcpp::List estimate_initial_temp_R(
    Rcpp::List ebrel_obj,
    Rcpp::Nullable<Rcpp::NumericVector> X0 = R_NilValue,
    double base_prob_X0 = 0.85,
    int    max_disp_thres   = 50,
    int    disp_boundary    = 30,
    double alpha            = 1.0,
    double beta             = 25.0,
    double gamma            = 100.0,
    double step_proportion  = 0.05,
    double step_probability = 0.05,
    int    num_samples      = 400,
    double chi0             = 0.8,
    double p                = 2.0,
    double tol_logchi       = 1e-3,
    int    max_iters        = 50,
    Rcpp::Nullable<Rcpp::NumericVector> T1 = R_NilValue,
    int    max_tries_factor = 50,
    double sigma            = 0.05,
    Rcpp::Nullable<Rcpp::IntegerVector> seed = R_NilValue,
    bool   verbose          = false
) {
  try {
    Rcpp::NumericVector E_r   = ebrel_obj["E"];
    Rcpp::NumericVector C_r   = ebrel_obj["C"];
    Rcpp::NumericVector U_r   = ebrel_obj["U"];
    Rcpp::NumericVector SD_r  = ebrel_obj["SD"];
    Rcpp::IntegerVector D_r   = ebrel_obj["D"];
    Rcpp::NumericVector SxH_r = ebrel_obj["SxH"];
    Rcpp::NumericVector O_r   = ebrel_obj["O"];

    int dim_x = Rcpp::as<int>(ebrel_obj["dim1"]);
    int dim_y = Rcpp::as<int>(ebrel_obj["dim2"]);
    int n_h   = Rcpp::as<int>(ebrel_obj["nh"]);
    int n_s   = Rcpp::as<int>(ebrel_obj["ns"]);

    if (O_r.size() != n_s) Rcpp::stop("`O` length must equal ns.");

    // Cast to std::vector
    std::vector<double> E   = Rcpp::as<std::vector<double>>(E_r);
    std::vector<double> C   = Rcpp::as<std::vector<double>>(C_r);
    std::vector<double> U   = Rcpp::as<std::vector<double>>(U_r);
    std::vector<double> SD  = Rcpp::as<std::vector<double>>(SD_r);
    std::vector<int>    D   = Rcpp::as<std::vector<int>>(D_r);
    std::vector<double> SxH = Rcpp::as<std::vector<double>>(SxH_r);
    std::vector<double> O   = Rcpp::as<std::vector<double>>(O_r);

    // lengths & rng
    const std::size_t expected_len =
      static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y) * static_cast<std::size_t>(n_h);
    if (U.size() != expected_len)
      Rcpp::stop("`U` size mismatch: expected %d, got %d.",
                 static_cast<int>(expected_len), static_cast<int>(U.size()));
    if (!(base_prob_X0 >= 0.0 && base_prob_X0 <= 1.0))
      Rcpp::stop("`base_prob_X0` must be in [0,1].");
    int rng_seed = seed.isNotNull() ? Rcpp::as<int>(seed.get()) : -1;

    // X_seed: PROVIDED or GENERATED
    std::vector<double> X_seed;
    if (X0.isNotNull()) {
      Rcpp::NumericVector Xs = X0.get();
      if (static_cast<std::size_t>(Xs.size()) != expected_len)
        Rcpp::stop("`Xseed` length mismatch: got %d, expected %d",
                   Xs.size(), static_cast<int>(expected_len));
      X_seed = Rcpp::as<std::vector<double>>(Xs);
    } else {
      X_seed = generate_X0_A(U, n_h, dim_x, dim_y, base_prob_X0, rng_seed);
      if (X_seed.size() != expected_len)
        Rcpp::stop("Internal error: generated X_seed has wrong length.");
    }

    // initial objective function evaluations and scaling as per run_ebrel function
    const double F1 = compute_F1(X_seed, C, n_h, dim_x, dim_y);
    const double F2 = compute_F2(X_seed, E, n_h, dim_x, dim_y);
    std::vector<double> G = compute_G(X_seed, E, SD, SxH, D,
                                      n_h, dim_x, dim_y,
                                      max_disp_thres, disp_boundary);
    const double G_sum = std::accumulate(G.begin(), G.end(), 0.0);
    const double eps   = 1e-12;
    const double alpha_scaled = alpha;
    const double beta_scaled  = (std::abs(F2) > eps) ? (beta  * F1 / F2)  : beta;
    const double gamma_scaled = (std::abs(G_sum) > eps) ? (gamma * F1 / G_sum) : gamma;

    std::vector<double> W = compute_distance_weights(E, U, n_h, dim_x, dim_y, sigma);
    if (W.size() != expected_len) Rcpp::stop("W size mismatch after compute_distance_weights");

    // T1: extract scalar from Nullable
    double T1_val = -1.0;
    if (T1.isNotNull()) {
      Rcpp::NumericVector t = T1.get();
      if (t.size() > 0) T1_val = t[0];
    }

    // IMPORTANT: pass X_seed (std::vector<double>), not X0 (Nullable)
    InitTempResult out = estimate_initial_temperature_benameur_cpp(
      X_seed, W, U, C, E, O, SD, SxH, D,
      max_disp_thres, disp_boundary,
      n_h, n_s, dim_x, dim_y,
      alpha_scaled, beta_scaled, gamma_scaled, base_prob_X0,
      step_proportion, step_probability,
      num_samples, chi0, p, tol_logchi,
      max_iters, T1_val, max_tries_factor,
      rng_seed, verbose
    );

    return Rcpp::List::create(
      Rcpp::Named("T0")            = out.T0,
      Rcpp::Named("chi_hat_final") = out.diag.chi_hat_final,
      Rcpp::Named("chi0")          = out.diag.chi0,
      Rcpp::Named("iters")         = out.diag.iters,
      Rcpp::Named("samples_used")  = out.diag.samples_used,
      Rcpp::Named("tries")         = out.diag.tries,
      Rcpp::Named("dE_median")     = out.diag.dE_median,
      Rcpp::Named("dE_mean")       = out.diag.dE_mean
    );
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("Unknown error in estimate_initial_temp_R");
  }
}
