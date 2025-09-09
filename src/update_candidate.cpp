//--------------------------- Update candidate ---------------------------------

#include "update_candidate.h"
#include <random>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <vector>

//-------------------------- Local helpers -------------------------------------

namespace {

  // Helper function for getting cell index for dim_x * dim_y * n_h vector
  static inline int idx(int cell, int h, int n_cells) {
    return h * n_cells + cell;
  }

}

//-------------------------- Main functions ------------------------------------


// Update candidate solution based on W (distance decay weighting)
std::vector<double> update_candidate(
    const std::vector<double>& W,            // [n_h * n_cells], habitat-major
    const std::vector<double>& U,            // [n_h * n_cells], 1=unavailable
    const std::vector<double>& candidate_in, // [n_h * n_cells], one-hot
    double step_proportion,                  // fraction of eligible cells to consider
    double step_probability,                 // mass for "pick a habitat"
    int n_h,
    int dim_x,
    int dim_y
) {
  const int n_cells = dim_x * dim_y;

  // Make a mutable copy to edit and return
  std::vector<double> candidate = candidate_in;

  // Collect eligible cells (any habitat available)
  std::vector<int> eligible;
  eligible.reserve(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int h = 0; h < n_h; ++h) {
      if (U[idx(i, h, n_cells)] == 0.0) { eligible.push_back(i); break; }
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
    for (int h = 0; h < n_h; ++h) s += W[idx(i, h, n_cells)];
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
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> cell_dist(eligible_probs.begin(), eligible_probs.end());

  std::vector<double> habitat_weights(n_h);
  std::vector<double> outcome_probs(n_h + 1); // [0]=no habitat, [1..n_h]=habitats

  // Update habitat (with replacement until filled)
  for (int s = 0; s < step_size; ++s) {
    const int i = eligible[cell_dist(gen)];

    // Masked habitat weights for this cell
    double sum_w = 0.0;
    for (int h = 0; h < n_h; ++h) {
      double w = (U[idx(i, h, n_cells)] == 0.0) ? W[idx(i, h, n_cells)] : 0.0;
      habitat_weights[h] = w;
      sum_w += w;
    }

    // Outcome: [no habitat] + per-habitat
    outcome_probs[0] = 1.0 - step_probability;
    if (sum_w > 0.0) {
      const double scale = step_probability / sum_w;
      for (int h = 0; h < n_h; ++h) outcome_probs[h + 1] = habitat_weights[h] * scale;
    } else {
      std::fill(outcome_probs.begin() + 1, outcome_probs.end(), 0.0);
      outcome_probs[0] = 1.0;
    }

    std::discrete_distribution<> outcome_dist(outcome_probs.begin(), outcome_probs.end());
    const int outcome = outcome_dist(gen); // 0 = no habitat, 1..n_h = chosen habitat

    // Clear existing assignment, then apply outcome
    for (int h = 0; h < n_h; ++h) candidate[idx(i, h, n_cells)] = 0.0;
    if (outcome > 0) candidate[idx(i, outcome - 1, n_cells)] = 1.0;
  }

  return candidate;
}
