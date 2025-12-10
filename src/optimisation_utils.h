//--------------------- Compute distance weights -------------------------------

#ifndef OPTIMISATION_UTILS_H
#define OPTIMISATION_UTILS_H

#include <vector>
#include <chrono>

/**
 * Compute distance–decay weights W for each (cell, habitat).
 *
 * Layout:
 * - Grid has dim_y rows and dim_x cols; n_cells = dim_x * dim_y.
 * - Inputs E and U are flattened in **habitat-major** order:
 *     index(h, cell) = h * n_cells + cell
 *   i.e., n_h contiguous blocks of length n_cells (one block per habitat).
 * - Output W uses the same habitat-major layout and size (n_h * n_cells).
 *
 * Distances:
 * - Multi-source BFS from all seed cells where E==1 for a given habitat.
 * - 8-connected neighbourhood with **unit-cost** steps (using Chebyshev).
 *
 * Weights:
 * - For each cell: w = exp(-sigma * d), where d is the integer hop distance.
 * - Cells are masked to zero weight if U==1 (unavailable) or E==1 (already present).
 * - For each habitat h, weights W[h,*] are normalized to sum to 1 if any positive mass exists.
 *
 */
std::vector<double> compute_distance_weights(
    const std::vector<double>& E,  // size: n_h * n_cells (habitat-major)
    const std::vector<double>& U,  // size: n_h * n_cells (habitat-major)
    int n_h,
    int dim_x,
    int dim_y,
    double sigma
);

inline double ms_since(const std::chrono::steady_clock::time_point& t0) {
  using namespace std::chrono;
  return duration_cast<duration<double, std::milli>>(steady_clock::now() - t0).count();
}

#endif // OPTIMISATION_UTILS_H
