#' @keywords internal
#' @useDynLib ebrel, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

# Add here any cpp wrapper functions that are to be exposed in R. (I want to
# keep this minimal).

#' Create an ebrel data object
#'
#' Initialises and returns an `ebrel` object used by `run_ebrel_R()`.
#'
#' @details
#' Let \code{n_cells <- dim_x * dim_y}. Vectors are flattened arrays with the
#' indicated lengths:
#' \itemize{
#'   \item \code{E}: cell × habitat (length \code{n_cells * n_h}); stacked by habitat within cell,
#'         i.e., index \code{cell + n_cells * (h - 1)}.
#'   \item \code{C}: cell × habitat costs (length \code{n_cells * n_h}); same layout as \code{E}.
#'   \item \code{SD}: cell × species baseline distribution (length \code{n_cells * n_s}).
#'   \item \code{SxH}: species × habitat association scores (length \code{n_s * n_h});
#'         column-major flattening of an \code{n_s × n_h} matrix, i.e., index
#'         \code{sp + n_s * (h - 1)}.
#' }
#' Cells with no data/unavailable should be marked with \code{sentinel}.
#'
#' @param E Numeric vector of length \code{n_cells * n_h}; encodes current habitat states.
#' @param C Numeric vector of length \code{n_cells * n_h}; encodes per-cell, per-habitat costs.
#' @param SD Numeric vector of length \code{n_cells * n_s}; encodes current species/feature distribution.
#' @param D Numeric vector of length \code{n_s}; species/feature dispersal potential (e.g., max distance).
#' @param SxH Numeric vector of length \code{n_s * n_h}; species–habitat association strengths.
#' @param O Numeric vector of length \code{n_s}; species/feature targets (e.g., proportional increase in occupancy).
#' @param dim_x Integer; number of grid cells along the x/longitude/easting dimension.
#' @param dim_y Integer; number of grid cells along the y/latitude/northing dimension.
#' @param n_h Integer; number of habitats considered in spatial prioritisation.
#' @param n_s Integer; number of species/features considered in spatial prioritisation.
#' @param sentinel Numeric; sentinel value for missing/unavailable cells.
#'
#' @return An list with the data needed to run ebrel.
#'
#' @seealso [run_ebrel_R()]
#'
#' @examples
#' \dontrun{
#' n_x <- 100; n_y <- 80
#' obj <- create_ebrel_class_object_R(E, C, SD, D, SxH, O, n_x, n_y, 5, 12)
#' }
#' @export
create_ebrel_class_object_R <- function(E, C, SD, D, SxH, O,  dim_x, dim_y, n_h, n_s, sentinel = 1e10) {
  .Call("_ebrel_create_ebrel_class_object_R", PACKAGE = "ebrel", E, C, SD, D, SxH, O, dim_x, dim_y, n_h, n_s, sentinel)
}

#' Run the ebrel optimization
#'
#' Runs simulated-annealing optimisation on an `ebrel` object.
#'
#' @param ebrel_obj External pointer created by `create_ebrel_class_object_R()`.
#' @param X0 Optional initial configuration.
#' @param base_prob_X0 Numeric; probability of assigning **no** habitat in a creation step.
#' @param sigma Numeric; controls distance weighting when selecting candidates.
#' @param max_disp_thres Numeric threshold after which universal dispersal is assumed.
#' @param disp_boundary Numeric factor setting the dispersal limit as \code{D[sp] * disp_boundary}.
#' @param alpha Numeric weight for target attainment.
#' @param beta Numeric weight for spatial aggregation.
#' @param gamma Numeric weight for costs.
#' @param temp Numeric initial temperature for simulated annealing.
#' @param step_proportion Numeric in (0, 1]; proportion of eligible cells to update per step.
#' @param step_probability Numeric in [0, 1]; probability of assigning any habitat to a cell.
#' @param n_iterations Integer > 0; maximum number of SA iterations.
#' @param cooling_rate_c Numeric > 0; cooling-rate control (used only when \code{lam_enabled = FALSE}).
#' @param lam_enabled Logical; enable Lam-style adaptive cooling rate control.
#' @param lam_target_mid Numeric target uphill acceptance rate for early/mid run.
#' @param lam_target_final Numeric target uphill acceptance rate at the end.
#' @param lam_hold_frac Numeric in [0, 1]; fraction of the run to hold \code{lam_target_mid} before decaying.
#' @param lam_p Numeric >= 1; damping exponent in Ben-Ameur correction.
#' @param min_iterations Integer >= 0; require at least this many iterations.
#' @param acceptance_window Integer > 0; window length for acceptance-rate estimation.
#' @param acceptance_thres Numeric; "low" acceptance threshold.
#' @param iter_no_improve Integer >= 0; consecutive iterations with no meaningful improvement.
#' @param improve_eps Numeric > 0; relative improvement needed to reset patience.
#' @param seed Integer; RNG seed (use -1 for non-deterministic/auto seed).
#' @param verbose Logical; print progress.
#'
#' @return A list. Replace with the exact structure you return, e.g.:
#' \itemize{
#'   \item \code{X_best}: best configuration found.
#'   \item \code{objective_trace}: numeric vector of objective values over iterations.
#'   \item \code{acceptance_rate}: numeric; overall or rolling acceptance rate.
#'   \item \code{iters}: integer; iterations executed.
#' }
#'
#' @examples
#' \dontrun{
#' res <- run_ebrel_R(obj, n_iterations = 2000, verbose = TRUE)
#' }
#' @export
run_ebrel_R <- function(ebrel_obj,
                        X0 = NULL,
                        max_disp_thres = 50,
                        disp_boundary = 30,
                        alpha = 1,
                        beta = 25,
                        gamma = 100,
                        sigma = 0.05,
                        base_prob_X0 = 0.85,
                        temp = 1000,
                        step_proportion = 0.05,
                        step_probability = 0.15,
                        n_iterations = 10000,
                        cooling_rate_c = 1,
                        lam_enabled = FALSE,
                        lam_target_mid = 0.44,
                        lam_target_final = 0.05,
                        lam_hold_frac = 0.60,
                        lam_p = 2.0,
                        min_iterations = 1000,
                        acceptance_window = 1000,
                        acceptance_thres = 0.01,
                        iter_no_improve = 1000,
                        improve_eps = 1e-6,
                        seed = -1,
                        verbose = FALSE
                        ) {
  .Call("_ebrel_run_ebrel_R", PACKAGE = "ebrel",
        ebrel_obj,
        X0,
        base_prob_X0,
        sigma,
        max_disp_thres,
        disp_boundary,
        alpha,
        beta,
        gamma,
        step_proportion,
        step_probability,
        n_iterations,
        temp,
        cooling_rate_c,
        lam_enabled,
        lam_target_mid,
        lam_target_final,
        lam_hold_frac,
        lam_p,
        min_iterations,
        acceptance_window,
        acceptance_thres,
        iter_no_improve,
        improve_eps,
        seed,
        verbose
        )
}


#' Estimate initial temperature for simulated annealing
#'
#' Chooses an initial SA temperature \eqn{T_1} so that the *uphill* move
#' acceptance rate is close to a target \code{chi0}. The routine samples random
#' neighbour proposals from \code{Xseed} (or an internal seed if \code{NULL}),
#' estimates the empirical uphill acceptance at the current temperature, and
#' updates \eqn{T} via a Ben-Ameur–style fixed-point iteration until
#' \eqn{|log(chi_hat) - log(chi0)| <= tol_logchi} or \code{max_iters} is reached.
#'
#' @param ebrel_obj External pointer created by \code{create_ebrel_class_object_R()}.
#' @param X0 Optional initial configuration (type/shape as expected by \code{run_ebrel_R()}).
#' @param max_disp_thres Numeric; distance (in cells) after which universal dispersal is assumed.
#' @param disp_boundary Numeric; dispersal limit factor applied as \code{D[sp] * disp_boundary}.
#' @param alpha Numeric; weight for target attainment in the objective.
#' @param beta Numeric; weight for spatial aggregation in the objective.
#' @param gamma Numeric; weight for costs in the objective.
#' @param step_proportion Numeric in (0, 1]; proportion of eligible cells considered per proposal.
#' @param step_probability Numeric in [0, 1]; probability of assigning any habitat during a proposal.
#' @param num_samples Integer > 0; number of random neighbour proposals used to estimate acceptance.
#' @param chi0 Numeric in (0, 1); target uphill acceptance rate (e.g., 0.8).
#' @param p Numeric >= 1; damping exponent in Ben-Ameur’s update.
#' @param tol_logchi Numeric > 0; tolerance for convergence on \code{log-acceptance}.
#' @param max_iters Integer >= 1; maximum number of fixed-point iterations.
#' @param T1_in Optional numeric; initial guess for \eqn{T_1}. If \code{NULL} or \eqn{\le 0}, a heuristic is used.
#' @param max_tries_factor Integer >= 1; cap on total neighbour proposals (as a multiple of \code{num_samples}) to avoid stalls.
#' @param sigma Numeric; distance-weighting parameter for candidate selection.
#' @param seed Optional integer; RNG seed (use \code{NULL} for non-deterministic seeding).
#' @param verbose Logical; print progress messages.
#'
#' @return
#' A numeric scalar or a list containing the chosen initial temperature and diagnostics.
#' Replace with the exact structure returned by the C++ backend, e.g.:
#' \itemize{
#'   \item \code{T1}: estimated initial temperature.
#'   \item \code{chi_hat}: empirical uphill acceptance at \code{T1}.
#'   \item \code{iters}: number of fixed-point iterations performed.
#'   \item \code{samples}: number of proposals evaluated.
#' }
#'
#' @seealso \code{\link{run_ebrel_R}}, \code{\link{create_ebrel_class_object_R}}
#'
#' @references
#' Ben-Ameur, W. (2004). Computing the initial temperature of simulated annealing.
#' \emph{Computational Optimization and Applications}, 29, 369–385.
#'
#' @examples
#' \dontrun{
#' # Assume `obj` was created by create_ebrel_class_object_R(...)
#' T0 <- estimate_initial_temp_R(
#'   ebrel_obj = obj,
#'   num_samples = 400,
#'   chi0 = 0.8,
#'   step_proportion = 0.05,
#'   step_probability = 0.05,
#'   verbose = TRUE
#' )
#' }
#' @export
estimate_initial_temp_R <- function(ebrel_obj,
                                    X0 = NULL,
                                    base_prob_X0 = 0.85,
                                    max_disp_thres = 50,
                                    disp_boundary = 30,
                                    alpha = 1.0,
                                    beta = 25.0,
                                    gamma = 100.0,
                                    step_proportion = 0.05,
                                    step_probability = 0.05,
                                    num_samples = 400,
                                    chi0 = 0.8, # target acceptance
                                    p = 2.0, # damping exponent in Ben-Ameur’s update (≥1)
                                    tol_logchi = 1e-3, # tolerance on |log χ̂ − log χ0|
                                    max_iters = 50, # cap on fixed-point iterations
                                    T1 = NULL, # optional initial guess; ≤0 => use heuristic
                                    max_tries_factor = 50, # limit on total neighbour proposals during
                                    sigma = 0.05,
                                    seed = NULL,
                                    verbose = FALSE) {
  .Call("_ebrel_estimate_initial_temp_R",
        PACKAGE = "ebrel",
        ebrel_obj,
        X0,
        base_prob_X0,
        max_disp_thres,
        disp_boundary,
        alpha,
        beta,
        gamma,
        step_proportion,
        step_probability,
        num_samples,
        chi0,
        p,
        tol_logchi,
        max_iters,
        T1,
        max_tries_factor,
        sigma,
        seed,
        verbose
  )
}
