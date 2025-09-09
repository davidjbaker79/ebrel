//------------------------- X0 functions ---------------------------------------

#include "generate_x_zero.h"

#include <vector>
#include <random>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <cctype>

//-------------------------- Local helper functions ----------------------------

namespace {

  // Number of cells
  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  }

  inline bool in_unit_interval(double p) {
    return std::isfinite(p) && p >= 0.0 && p <= 1.0;
  }

}

//----------------------------- Main functions ---------------------------------

// Naive initialize of habitat configuration matrix X0
std::vector<double> generate_X0_A(const std::vector<double>& U,
                                  int n_h,
                                  int dim_x,
                                  int dim_y,
                                  double base_prob,
                                  int rng_seed) {
  // --- validation ---
  if (dim_x <= 0 || dim_y <= 0 || n_h <= 0) {
    throw std::invalid_argument("dim_x, dim_y, and n_h must be positive.");
  }
  if (!in_unit_interval(base_prob)) {
    throw std::invalid_argument("base_prob must be in [0, 1].");
  }
  const std::size_t cells   = n_cells(dim_x, dim_y);
  const std::size_t expected = cells * static_cast<std::size_t>(n_h);
  if (U.size() != expected) {
    throw std::invalid_argument("U size mismatch: got " + std::to_string(U.size()) +
                                ", expected " + std::to_string(expected));
  }

  // Empty case
  if (base_prob >= 1.0) {
    return std::vector<double>(expected, 0.0);
  }

  std::vector<double> X0(expected, 0.0);

  // RNG: deterministic if rng_seed >= 0, random_device otherwise (for checking)
  std::mt19937 rng;
  if (rng_seed >= 0) {
    rng.seed(static_cast<std::mt19937::result_type>(rng_seed));
  } else {
    std::random_device rd;
    rng.seed(rd());
  }
  std::uniform_real_distribution<double> unif01(0.0, 1.0);
  std::uniform_int_distribution<int> habitat_dist(0, n_h - 1);

  // For each tile, decide empty vs assign a random habitat (if available)
  // This version tries to fill the tile with one of the available habitats
  // rather than not filling if the habitat is unavailable
  for (std::size_t tile = 0; tile < cells; ++tile) {
    const double u = unif01(rng);
    if (u < base_prob) {
      continue;
    }
    int h0 = habitat_dist(rng);
    for (int k = 0; k < n_h; ++k) {
      int h = (h0 + k) % n_h; // probe all habitats in a rotated order
      std::size_t idx = static_cast<std::size_t>(h) * cells + tile;
      if (U[idx] == 0.0) { X0[idx] = 1.0; break; }
    }
  }

  return X0;
}
