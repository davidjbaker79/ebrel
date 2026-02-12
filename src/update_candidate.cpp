//--------------------------- Update candidate ---------------------------------

#include "update_candidate.h"
#include <random>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <cmath>

//-------------------------- Local helpers -------------------------------------

namespace {

  // Helper function for getting cell index for dim_x * dim_y * n_h vector
  static inline int idx(int cell, int h, int n_cells) {
    return h * n_cells + cell;
  }

}

//-------------------------- Main functions ------------------------------------

// Update candidate solution based on W (distance decay weighting)
std::vector<int8_t> update_candidate(
    const std::vector<double>& W,            // [n_h * n_cells], habitat-major
    const std::vector<uint8_t>& U,           // [n_h * n_cells], 1=unavailable
    const std::vector<int8_t>& candidate_in, // [n_cells], -1 or habitat index
    double step_proportion,                  // fraction of eligible cells to consider
    double step_probability,                 // mass for "pick a habitat"
    int n_h,
    int dim_x,
    int dim_y,
    int rng_seed = -1
) {
  const int n_cells = dim_x * dim_y;

  if ((int)candidate_in.size() != n_cells) {
    throw std::runtime_error("candidate_in must have length dim_x*dim_y");
  }
  if ((int)W.size() != n_h * n_cells || (int)U.size() != n_h * n_cells) {
    throw std::runtime_error("W and U must have length n_h*dim_x*dim_y");
  }
  if (step_proportion <= 0.0 || step_probability < 0.0) return candidate_in;
  if (step_probability > 1.0) throw std::runtime_error("step_probability must be <= 1");


  // Make a mutable copy to edit and return
  std::vector<int8_t> candidate = candidate_in;

  // Collect eligible cells (any habitat available)
  std::vector<int> eligible;
  eligible.reserve(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int h = 0; h < n_h; ++h) {
      if (U[idx(i, h, n_cells)] == 0) {
        eligible.push_back(i);
        break;
      }
    }
  }
  if (eligible.empty()) return candidate;

  // How many cells to try
  int step_size = static_cast<int>(std::round(step_proportion * static_cast<double>(eligible.size())));
  if (step_size <= 0) return candidate;

  // Per-cell probs p_i = sum_h W[h,i]
  std::vector<double> cell_probs(n_cells, 0.0);
  for (int i = 0; i < n_cells; ++i) {
    double s = 0.0;
    for (int h = 0; h < n_h; ++h) {
      s += W[idx(i, h, n_cells)];
    }
    cell_probs[i] = s;
  }

  // Restrict to eligible & normalize
  std::vector<double> eligible_probs;
  eligible_probs.reserve(eligible.size());
  for (int i : eligible) eligible_probs.push_back(cell_probs[i]);

  double total = std::accumulate(eligible_probs.begin(), eligible_probs.end(), 0.0);
  if (total <= 0.0) throw std::runtime_error("update_candidate: All values in W are zero.");
  for (double& p : eligible_probs) p /= total;

  // Random seeds
  std::mt19937 gen;
  if (rng_seed >= 0) gen.seed(static_cast<uint32_t>(rng_seed));
  else gen.seed(std::random_device{}());

  // Weights and probs
  std::discrete_distribution<> cell_dist(eligible_probs.begin(), eligible_probs.end());
  std::vector<double> habitat_weights(n_h);
  std::vector<double> outcome_probs(n_h + 1); // [0]=no habitat, [1..n_h]=habitats

  // Update habitat (with replacement until filled)
  for (int s = 0; s < step_size; ++s) {
    const int i = eligible[cell_dist(gen)];

    // Masked habitat weights for this cell
    double sum_w = 0.0;
    for (int h = 0; h < n_h; ++h) {
      const double w = (U[idx(i, h, n_cells)] == 0) ? W[idx(i, h, n_cells)] : 0.0;
      habitat_weights[h] = w;
      sum_w += w;
    }

    // Outcome: [no habitat] + per-habitat
    outcome_probs[0] = 1.0 - step_probability;
    if (sum_w > 0.0) {
      const double scale = step_probability / sum_w;
      for (int h = 0; h < n_h; ++h) {
        outcome_probs[h + 1] = habitat_weights[h] * scale;
      }
    } else {
      std::fill(outcome_probs.begin() + 1, outcome_probs.end(), 0.0);
      outcome_probs[0] = 1.0;
    }

    std::discrete_distribution<> outcome_dist(outcome_probs.begin(), outcome_probs.end());
    const int outcome = outcome_dist(gen); // 0 = no habitat, 1..n_h = chosen habitat

    if (outcome == 0) candidate[i] = static_cast<int8_t>(-1);
    else candidate[i] = static_cast<int8_t>(outcome - 1); // 0..n_h-1
  }
  return candidate;
}
