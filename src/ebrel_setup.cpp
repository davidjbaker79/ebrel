//--------------------------- Ebrel setup functions ----------------------------

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstddef>

// ------------------------- Local helper functions ----------------------------

namespace {

  // Just to ensure that unavailable areas aren't misclassified due to small numerical difference
  inline bool approx_equal(double x, double y, double rel = 1e-12, double abs_tol = 1e-12) {
    const double diff  = std::abs(x - y);
    const double scale = std::max(std::abs(x), std::abs(y));
    return diff <= std::max(abs_tol, rel * scale);
  }

  // Calculates number of cells based on dims
  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  }

}

//---------------------- Main functions ----------------------------------------

// Create binary mask: 1 if cost ≈ sentinel, else 0
std::vector<double> get_unavailable_mask(const std::vector<double>& C,
                                         int n_h,
                                         int dim_x,
                                         int dim_y,
                                         double sentinel) {
  const std::size_t cells    = n_cells(dim_x, dim_y);
  const std::size_t expected = cells * static_cast<std::size_t>(n_h);
  if (C.size() != expected) {
    throw std::runtime_error("C size mismatch in get_unavailable_mask: got " +
                             std::to_string(C.size()) + ", expected " + std::to_string(expected));
  }

  // Compare using approx_equal (avoids some issues with small differences in sentinel)
  std::vector<double> U(expected);
  for (std::size_t i = 0; i < expected; ++i) {
    U[i] = approx_equal(C[i], sentinel) ? 1.0 : 0.0;
  }
  return U;
}

// Ensure each cell in E has at most one habitat set to 1 (habitat-major)
// i.e. in this binary model a grid cell can have at most one existing habitat type.
// Will expand to consider fractional cover, but that has implications for connectivity so requires thought.
void validate_existing_habitat(const std::vector<double>& E,
                               int dim_x,
                               int dim_y,
                               int n_h) {
  const std::size_t cells    = n_cells(dim_x, dim_y);
  const std::size_t expected = cells * static_cast<std::size_t>(n_h);
  if (E.size() != expected) {
    throw std::runtime_error("E size mismatch in validate_existing_habitat: got " +
                             std::to_string(E.size()) + ", expected " + std::to_string(expected));
  }

  for (std::size_t cell = 0; cell < cells; ++cell) {
    int count = 0;
    for (int h = 0; h < n_h; ++h) {
      const std::size_t idx = static_cast<std::size_t>(h) * cells + cell;
      if (E[idx] == 1.0 && ++count > 1) {
        throw std::runtime_error(
            "Each cell in E must have at most one habitat with value 1 (found >1 at cell " +
              std::to_string(cell) + ")."
        );
      }
    }
  }
}
