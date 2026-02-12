//-------------------- Pre-compute species plan of modelling -------------------

#include <cmath>
#include <algorithm>
#include <vector>

#include "species_plan.h"

//-------------------------- Local helper functions ----------------------------

namespace {

// Get number of cells
inline std::size_t n_cells(int dim_x, int dim_y) {
  return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
}

}

//----------------------------- Main functions ---------------------------------

// ---- Function to build species' plan (e.g., seed, roi etc) ----
SpeciesPlan build_species_plan(
    const std::vector<double>& SD,
    const std::vector<double>& SxH,   // row-major sp*n_h + h
    const std::vector<int>&    D,     // disp_raw per species
    const std::vector<double>& O,     // targets
    int n_h, int n_s,
    int dim_x, int dim_y,
    int universal_disp_thres,
    int max_disp_steps,
    int roi_cap,
    const std::vector<int16_t>& row_first_land,
    const std::vector<int16_t>& row_last_land,
    const std::vector<int16_t>& col_first_land,
    const std::vector<int16_t>& col_last_land,
    const std::vector<int8_t>& E,
    const std::vector<int>& cell_r,
    const std::vector<int>& cell_c
) {
  const std::size_t cells = n_cells(dim_x, dim_y);

  if (SD.size() != static_cast<std::size_t>(n_s) * cells)
    throw std::runtime_error("SD size mismatch in build_species_plan_soa");
  if (SxH.size() != static_cast<std::size_t>(n_s) * static_cast<std::size_t>(n_h))
    throw std::runtime_error("SxH size mismatch in build_species_plan_soa");
  if (D.size() != static_cast<std::size_t>(n_s))
    throw std::runtime_error("D size mismatch in build_species_plan_soa");

  SpeciesPlan P;
  P.n_s = n_s;
  P.n_h = n_h;

  P.active.assign(n_s, 0u);
  P.fast_path.assign(n_s, 0u);
  P.disp.assign(n_s, 0);
  P.r0.assign(n_s, 0); P.r1.assign(n_s, -1);
  P.c0.assign(n_s, 0); P.c1.assign(n_s, -1);
  P.W.assign(n_s, 0);  P.H.assign(n_s, 0);
  P.O_n.assign(n_s, 0u);

  P.suitable.assign(static_cast<std::size_t>(n_s) * static_cast<std::size_t>(n_h), 0u);

  P.seed_start.assign(n_s, 0u);
  P.seed_len.assign(n_s, 0u);
  P.seed_pool.clear();
  P.seed_pool.reserve(1024u * static_cast<std::size_t>(n_s)); // heuristic

  for (int sp = 0; sp < n_s; ++sp) {

    // ---- Suitable flags for this species (contiguous in P.suitable)
    const std::size_t base_h = static_cast<std::size_t>(sp) * static_cast<std::size_t>(n_h);
    for (int h = 0; h < n_h; ++h) {
      P.suitable[base_h + static_cast<std::size_t>(h)] =
        (SxH[base_h + static_cast<std::size_t>(h)] > 0.0) ? 1u : 0u;
    }

    const int disp_raw = D[sp];
    const int disp = std::min(disp_raw, universal_disp_thres);
    P.disp[sp] = static_cast<int16_t>(disp);
    P.fast_path[sp] = (disp_raw >= universal_disp_thres) ? 1u : 0u;

    // ---- Seeds + bbox (but store seeds into global flat pool)
    const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;

    int min_r = dim_y, max_r = -1, min_c = dim_x, max_c = -1;

    const uint32_t start = static_cast<uint32_t>(P.seed_pool.size());
    P.seed_start[sp] = start;

    for (std::size_t k = 0; k < cells; ++k) {
      if (SD[base_sd + k] != 1.0) continue;

      const int ehx = E[k];
      if (ehx == -1) continue;
      if (!P.suitable[base_h + static_cast<std::size_t>(ehx)]) continue;

      P.seed_pool.push_back(static_cast<SeedT>(k));

      const int r = cell_r[k];
      const int c = cell_c[k];
      if (r < min_r) min_r = r;
      if (r > max_r) max_r = r;
      if (c < min_c) min_c = c;
      if (c > max_c) max_c = c;
    }

    // ---- Number of seeds
    const uint32_t len = static_cast<uint32_t>(P.seed_pool.size()) - start;
    P.seed_len[sp] = len;

    // ---- Targets in number of cells to gain
    const double o = O[static_cast<std::size_t>(sp)]; // in [0,1]
    if (len == 0u) {
      P.O_n[sp] = 0u;   // If no seeds -> species inactive anyway, but keep target 0
    } else {
      // target = ceil(len * o), with tiny epsilon to avoid FP bump-ups
      const double x = static_cast<double>(len) * o;
      uint32_t target = static_cast<uint32_t>(std::ceil(x - 1e-12));

      // Enforce minimum 1
      if (target < 1u) target = 1u;

      // Enforce maximum len (for safety clamp)
      if (target > len) target = len;

      P.O_n[sp] = target;
    }

    // ---- Mark inactive species
    if (len == 0u) {
      // inactive species
      P.active[sp] = 0u;
      continue;
    }

    // ---- ROI expand
    const int R_raw = (max_disp_steps > 0)
      ? (disp_raw * max_disp_steps)
      : std::max(dim_x, dim_y);
    const int R = std::min(R_raw, roi_cap);

    int r0 = std::max(0,         min_r - R);
    int r1 = std::min(dim_y - 1, max_r + R);
    int c0 = std::max(0,         min_c - R);
    int c1 = std::min(dim_x - 1, max_c + R);

    // tighten vertically using row spans
    while (r0 <= r1 && (row_last_land[r0] < c0 || row_first_land[r0] > c1)) ++r0;
    while (r1 >= r0 && (row_last_land[r1] < c0 || row_first_land[r1] > c1)) --r1;
    if (r0 > r1) { P.active[sp] = 0u; continue; }

    // tighten horizontally using col spans
    while (c0 <= c1 && (col_last_land[c0] < r0 || col_first_land[c0] > r1)) ++c0;
    while (c1 >= c0 && (col_last_land[c1] < r0 || col_first_land[c1] > r1)) --c1;
    if (c0 > c1) { P.active[sp] = 0u; continue; }

    const int W = c1 - c0 + 1;
    const int H = r1 - r0 + 1;
    if (W <= 0 || H <= 0) { P.active[sp] = 0u; continue; }

    P.r0[sp] = r0; P.r1[sp] = r1;
    P.c0[sp] = c0; P.c1[sp] = c1;
    P.W[sp]  = W;  P.H[sp]  = H;
    P.active[sp] = 1u;
  }

  return P;
}
