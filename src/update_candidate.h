//--------------------------- Update candidate ---------------------------------

#ifndef UPDATE_CANDIDATE_H
#define UPDATE_CANDIDATE_H

#include <vector>
#include <cstdint>

/**
 * Update candidate solution for simulated annealing
 *
 * All matrices are flattened row-major arrays of size [n_cells * n_h].
 * - W: per-(cell, habitat) weights, typically column-normalized.
 * - U: availability mask (1 = unavailable, 0 = available).
 * - candidate_in: current one-hot configuration (values 0/1).
 *
 * Returns a new flattened candidate vector after applying updates.
 */
std::vector<int8_t> update_candidate(
    const std::vector<double>& W,             // [n_cells * n_h]
    const std::vector<uint8_t>& U,             // [n_cells * n_h], 1 = unavailable
    const std::vector<int8_t>& candidate_in,  // [n_cells * n_h], one-hot
    double step_proportion,                   // proportion of eligible cells to update
    double step_probability,                  // probability of assigning any habitat
    int n_h,
    int dim_x,
    int dim_y,
    int rng_seed);

#endif // UPDATE_CANDIDATE_H
