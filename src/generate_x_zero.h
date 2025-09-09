//------------------------- X0 functions ---------------------------------------

#ifndef GENERATE_X_ZERO_H
#define GENERATE_X_ZERO_H

#include <vector>

// Generate X0 (initial configuration) with an "empty" option:
//  - With probability base_prob, leave the tile empty (all 0s across habitats).
//  - Otherwise choose a habitat uniformly at random and set that [tile, habitat] to 1,
//    but only if U[idx] == 0 (available). If unavailable, leaves empty. If the chosen
//    habitat is not available then continue trying until filled.
// rng_seed < 0 => non-deterministic; >=0 => deterministic

// Random X0 generation (naive) - this is entirely random
std::vector<double> generate_X0_A(const std::vector<double>& U,
                                  int n_h,
                                  int dim_x,
                                  int dim_y,
                                  double base_prob,
                                  int rng_seed = -1);

#endif  // GENERATE_X_ZERO_H
