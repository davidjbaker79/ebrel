//--------------------- Compute distance weights -------------------------------

#ifndef OPTIMISATION_UTILS_H
#define OPTIMISATION_UTILS_H

#include <vector>
#include <chrono>
#include <stdexcept>

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
    const std::vector<int8_t>& E,  // size: n_h * n_cells (habitat-major)
    const std::vector<uint8_t>& U,  // size: n_h * n_cells (habitat-major)
    int n_h,
    int dim_x,
    int dim_y,
    double sigma
);

// For benchmarking
inline double ms_since(const std::chrono::steady_clock::time_point& t0) {
  using namespace std::chrono;
  return duration_cast<duration<double, std::milli>>(steady_clock::now() - t0).count();
}


// Convert dense habitat-major X (h*cells + g) to per-cell habitat id
// Rule: argmax per cell; if max <= eps → -1
inline std::vector<int> dense_X_to_cell_habitat(const std::vector<double>& X,
                                                int n_h,
                                                std::size_t cells,
                                                double eps = 0.0)
{
  const std::size_t expected =
    static_cast<std::size_t>(n_h) * cells;

  if (X.size() != expected) {
    throw std::runtime_error(
        "`X` must have length n_h * (dim_x * dim_y).");
  }

  std::vector<int> X_h_of_cell(cells, -1);
  const double* Xptr = X.data();

  for (std::size_t g = 0; g < cells; ++g) {
    int best_h = -1;
    double best_v = eps;

    for (int h = 0; h < n_h; ++h) {
      const double v = Xptr[static_cast<std::size_t>(h) * cells + g];
      if (v > best_v) {
        best_v = v;
        best_h = h;
      }
    }

    X_h_of_cell[g] = best_h; // -1 if none
  }

  return X_h_of_cell;
}

#endif // OPTIMISATION_UTILS_H
