//------------------------------ Dispersal-aware G function --------------------

#include "dispersal_utils.h"

#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <limits> //  for std::numeric_limits
#include <cmath>   // for std::abs

#ifdef _OPENMP
  #include <omp.h>
#endif

//-------------------------- Local helpers -------------------------------------

namespace {

  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  }

  // Chebyshev (chessboard) neighbourhood offsets up to distance d
  // - used to better handle diagonal movements
  std::vector<std::pair<int,int>> make_offsets(int d) {
    std::vector<std::pair<int,int>> off;
    off.reserve((2*d + 1) * (2*d + 1) - 1);
    for (int dr = -d; dr <= d; ++dr) {
      for (int dc = -d; dc <= d; ++dc) {
        if (dr == 0 && dc == 0) continue;
        if (std::max(std::abs(dr), std::abs(dc)) <= d) {
          off.emplace_back(dr, dc);
        }
      }
    }
    return off;
  }

}

//------------------------- G Function -----------------------------------------

// Compute G(x): per-species count of reachable newly created suitable habitat (X),
// using E (existing) and suitable-X as transit.
// Compute G(x): per-species count of reachable newly created suitable habitat (X),
// using E (existing) and suitable-X as transit.
std::vector<double> compute_G(const std::vector<double>& X,
                              const std::vector<double>& E,
                              const std::vector<double>& SD,
                              const std::vector<double>& SxH, // row-major: sp * n_h + h.
                              const std::vector<int>& D,
                              int n_h,
                              int n_s,
                              int dim_x,
                              int dim_y,
                              int universal_disp_thres,
                              int max_disp_steps,
                              int roi_cap,
                              // Data to avoid wasted work over sea (important for UK analysis)
                              const std::vector<uint8_t>& LM, // Land mask
                              const std::vector<int>& row_first_land,
                              const std::vector<int>& row_last_land,
                              const std::vector<int>& col_first_land,
                              const std::vector<int>& col_last_land
                              )
{
  // ---- Number of grid cells ----
  const std::size_t cells = n_cells(dim_x, dim_y);

  // ---- Vector of storing results per species ----
  std::vector<double> result(n_s, 0.0);

  // ---- Precompute row/col coords [TESTED in codePad & WORKS] ----
  std::vector<int> cell_r(cells), cell_c(cells);
  for (std::size_t k = 0; k < cells; ++k) {
    cell_r[k] = static_cast<int>(k / static_cast<std::size_t>(dim_x));
    cell_c[k] = static_cast<int>(k % static_cast<std::size_t>(dim_x));
  }

  // ---- Cache offsets for all distances [TESTED in codePad & WORKS] ----
  std::vector<std::vector<std::pair<int,int>>> offsets_cache(universal_disp_thres + 1);
  for (int d = 0; d <= universal_disp_thres; ++d) {
    offsets_cache[d] = make_offsets(d);
  }

  // ---- One-hot habitat index per cell for X and E (-1 = none) [TESTED in codePad & WORKS] ----
  std::vector<int> X_h_of_cell(cells, -1), E_h_of_cell(cells, -1);
  for (int h = 0; h < n_h; ++h) {
    const std::size_t base = static_cast<std::size_t>(h) * cells;
    for (std::size_t k = 0; k < cells; ++k) {
      if (X[base + k] == 1.0) X_h_of_cell[k] = h;
      if (E[base + k] == 1.0) E_h_of_cell[k] = h;
    }
  }


  // Steps distance in BFS layers.
  using step_t = uint16_t;
  const step_t STEP_MAX = std::numeric_limits<step_t>::max();

#ifdef _OPENMP

  //======================= PARALLEL IMPLEMENTATION ============================
  // Uses a parallel region with thread-local scratch buffers.

#pragma omp parallel
{
  // Thread-local scratch arrays - unsized and adjusted by species ROI
  std::vector<uint8_t> transit;   // ROI-sized per species
  std::vector<uint8_t> habitat;   // ROI-sized per species
  std::vector<step_t>  steps;     // ROI-sized per species

  std::vector<std::size_t> frontier;
  std::vector<std::size_t> next;

  // ---- Parallelise over species with dynamic scheduling ----
#pragma omp for schedule(dynamic)
  for (int sp = 0; sp < n_s; ++sp) {

    // ---- Per-species suitability flags ----
    std::vector<uint8_t> suitable_h_flag(n_h, 0u);
    const std::size_t base_h = static_cast<std::size_t>(sp) * static_cast<std::size_t>(n_h);
    for (int h = 0; h < n_h; ++h) {
      suitable_h_flag[h] = (SxH[base_h + static_cast<std::size_t>(h)] > 0.0) ? 1u : 0u;
    }

    // ---- Species dispersal (bounded by max) ----
    const int disp_raw = D[sp];
    const int disp = std::min(disp_raw, universal_disp_thres);

    // ---- Find seeds and their bounding box ----
    const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;
    std::vector<std::size_t> seed_idxs;
    seed_idxs.reserve(1024);

    // This creates a tight bbox around distribution seeds
    int min_r = dim_y, max_r = -1, min_c = dim_x, max_c = -1;
    for (std::size_t k = 0; k < cells; ++k) {
      if (SD[base_sd + k] != 1.0) continue;
      const int ehx = E_h_of_cell[k];
      if (ehx == -1 || !suitable_h_flag[ehx]) continue;
      seed_idxs.push_back(k);
      const int r = cell_r[k], c = cell_c[k];
      if (r < min_r) min_r = r;
      if (r > max_r) max_r = r;
      if (c < min_c) min_c = c;
      if (c > max_c) max_c = c;
    }
    if (seed_idxs.empty()) {
      result[sp] = 0.0; // if no seeds found then zero and skip
      continue;
    }

    // ---- ROI ()expand seed bbox by R = disp_raw * max_disp_steps ----
    // Here, the ROI is expanded
    // [later, if the universal dispersal threshold is tripped all the suitable
    // cells are considered accessible]
    const int R_raw = (max_disp_steps > 0) ? (disp_raw * max_disp_steps) : std::max(dim_x, dim_y); //defaults to universal
    const int R = std::min(R_raw, roi_cap); // Sets upper limit on ROI expansion outwards

    int r0 = std::max(0,         min_r - R); // create ROI boundaries
    int r1 = std::min(dim_y - 1, max_r + R);
    int c0 = std::max(0,         min_c - R);
    int c1 = std::min(dim_x - 1, max_c + R);

    // Tighten vertically using row spans vs. [c0,c1] to clip sea
    while (r0 <= r1 && (row_last_land[r0] < c0 || row_first_land[r0] > c1)) ++r0;
    while (r1 >= r0 && (row_last_land[r1] < c0 || row_first_land[r1] > c1)) --r1;
    if (r0 > r1) { result[sp] = 0.0; continue; }

    // Tighten horizontally using column spans vs. [r0,r1] to clip sea
    while (c0 <= c1 && (col_last_land[c0] < r0 || col_first_land[c0] > r1)) ++c0;
    while (c1 >= c0 && (col_last_land[c1] < r0 || col_first_land[c1] > r1)) --c1;
    if (c0 > c1) { result[sp] = 0.0; continue; }

    const int W = c1 - c0 + 1;
    const int H = r1 - r0 + 1;

    // ---- Adjust the size of the masks by species ROI ----
    const std::size_t Nroi = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);

    if (Nroi > transit.capacity()) {
      transit.reserve(Nroi);
      habitat.reserve(Nroi);
      steps.reserve(Nroi);
    }
    transit.resize(Nroi);
    habitat.resize(Nroi);
    steps.resize(Nroi);

    frontier.clear(); next.clear();
    frontier.reserve(Nroi); next.reserve(Nroi);

    // ---- Fast-path (effectively global dispersal, skipping BFS) ----
    if (disp_raw >= universal_disp_thres) {
      double total = 0.0;
      for (int r = r0; r <= r1; ++r) {
        std::size_t g = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(c0);
        for (int c = c0; c <= c1; ++c, ++g) {
          if (!LM[g]) continue; // skip sea cells
          const int hx = X_h_of_cell[g];
          if (hx != -1 && suitable_h_flag[hx]) {
            total += 1.0;
          }
        }
      }
      result[sp] = total;
      continue; // [skip BFS]
    }

    // ---- Build transit/habitat masks within ROI (with sea skipping) ----
    std::size_t roiH = 0;
    {
      for (int r = r0; r <= r1; ++r) {
        // Row-local starting indices
        const std::size_t rowOff = static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W);

        // If the row has no land overlapping [c0,c1], just zero and continue
        if (row_last_land[r] < c0 || row_first_land[r] > c1) {
          // zero the whole row slice once
          std::fill_n(&habitat[rowOff], W, 0u);
          std::fill_n(&transit[rowOff], W, 0u);
          std::fill_n(&steps[rowOff],   W, STEP_MAX);
          continue;
        }

        // Pre-zero entire row slice (only set non-zeros in the land span)
        std::fill_n(&habitat[rowOff], W, 0u);
        std::fill_n(&transit[rowOff], W, 0u);
        std::fill_n(&steps[rowOff],   W, STEP_MAX);

        // Work only where this row actually has land within the ROI
        const int cc0 = std::max(c0, row_first_land[r]);
        const int cc1 = std::min(c1, row_last_land[r]);

        // Global index for the first land cell in this row's ROI
        std::size_t g = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
          + static_cast<std::size_t>(cc0);

        // Row-local index for the first land cell
        std::size_t idx = rowOff + static_cast<std::size_t>(cc0 - c0);

        for (int c = cc0; c <= cc1; ++c, ++idx, ++g) {

          uint8_t hmask = 0u, tmask = 0u; // Masks

          const int hx = X_h_of_cell[g];
          if (hx != -1 && suitable_h_flag[hx]) {
            hmask = 1u;
            tmask = 1u; // X counts as transit too
          }
          if (!tmask) {
            const int he = E_h_of_cell[g];
            if (he != -1 && suitable_h_flag[he]) {
              tmask = 1u; // Suitable E as transit only
            }
          }

          habitat[idx] = hmask;
          transit[idx] = tmask;

          // Early out if no suitable habitat
          roiH += hmask;
        }
      }
    }

    // ---- Early out (no new suitable X in ROI) ----
    if (roiH == 0) {
      result[sp] = 0.0;
      continue;
    }

    // ---- Seed frontier (ROI-local, only on transit) ----
    frontier.clear();
    {
      for (std::size_t k : seed_idxs) {
        const int r = cell_r[k];
        const int c = cell_c[k];
        if (r < r0 || r > r1 || c < c0 || c > c1) continue;

        // Local index
        const std::size_t idx =
          static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W)
          + static_cast<std::size_t>(c - c0);

        // global index
        const std::size_t global_k =
        static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
          + static_cast<std::size_t>(c);

        if (!LM[global_k]) continue;        // Check just in case
        if (!transit[idx]) continue;        // Must be transit
        if (steps[idx] == 0) continue;      // Avoid duplicate seeds

        steps[idx] = 0;                     // Distance 0 at seeds
        frontier.push_back(idx);
      }
    }

    if (frontier.empty()) {
      result[sp] = 0.0;
      continue;
    }

    // ---- Frontier-by-layers BFS (depth-limited) ----
    // Count reachable habitat on-the-fly
    const auto& offsets = offsets_cache[disp];
    double g_sp = 0.0;

    for (step_t d = 0;
         !frontier.empty() &&
           (max_disp_steps <= 0 || d < static_cast<step_t>(max_disp_steps));
         ++d)
    {
      next.clear();

      for (std::size_t i : frontier) {

        const std::size_t row = i / static_cast<std::size_t>(W);
        const std::size_t col = i - row * static_cast<std::size_t>(W);
        const int r = r0 + static_cast<int>(row);
        const int c = c0 + static_cast<int>(col);

        for (const auto& off : offsets) {
          const int nr = r + off.first;
          const int nc = c + off.second;

          // Bounds check
          if (static_cast<unsigned>(nr - r0) >= static_cast<unsigned>(H) ||
              static_cast<unsigned>(nc - c0) >= static_cast<unsigned>(W)) {
            continue;
          }

          const std::size_t j =
            static_cast<std::size_t>(nr - r0) * static_cast<std::size_t>(W)
            + static_cast<std::size_t>(nc - c0);

          if (!transit[j]) continue;

          const step_t next_step = static_cast<step_t>(d + 1);
          if (steps[j] <= next_step) continue; // Already reached earlier / good

          steps[j] = next_step;

          // Count new suitable X habitat as soon as it becomes reachable.
          // Seeds were at step 0; only increment on first reach with next_step>=1.
          if (habitat[j]) {
            g_sp += 1.0;
          }

          next.push_back(j);
        }
      }

      frontier.swap(next);
    }

    result[sp] = g_sp;
  } // for sp
} // parallel region

#else  // !_OPENMP

//======================= SERIAL IMPLEMENTATION ==============================

// ---- Single set of scratch arrays - unsized as will adjust by species ROI ----
std::vector<uint8_t> transit;   // ROI-sized per species
std::vector<uint8_t> habitat;   // ROI-sized per species
std::vector<step_t>  steps;     // ROI-sized per species

std::vector<std::size_t> frontier;
std::vector<std::size_t> next;

for (int sp = 0; sp < n_s; ++sp) {

  // ---- Per-species suitability flags ----
  std::vector<uint8_t> suitable_h_flag(n_h, 0u);
  const std::size_t base_h = static_cast<std::size_t>(sp) * static_cast<std::size_t>(n_h);
  for (int h = 0; h < n_h; ++h) {
    suitable_h_flag[h] = (SxH[base_h + static_cast<std::size_t>(h)] > 0.0) ? 1u : 0u;
  }

  // ---- Species dispersal ----
  const int disp_raw = D[sp];
  const int disp = std::min(disp_raw, universal_disp_thres);

  // ---- ROI (expand seed bbox by R = disp_raw * max_disp_steps) ----
  const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;
  std::vector<std::size_t> seed_idxs;
  seed_idxs.reserve(1024);

  // This creates a tight bbox around distribution seeds
  int min_r = dim_y, max_r = -1, min_c = dim_x, max_c = -1;
  for (std::size_t k = 0; k < cells; ++k) {
    if (SD[base_sd + k] != 1.0) continue;
    const int ehx = E_h_of_cell[k];
    if (ehx == -1 || !suitable_h_flag[ehx]) continue;
    seed_idxs.push_back(k);
    const int r = cell_r[k], c = cell_c[k];
    if (r < min_r) min_r = r;
    if (r > max_r) max_r = r;
    if (c < min_c) min_c = c;
    if (c > max_c) max_c = c;
  }
  if (seed_idxs.empty()) {
    result[sp] = 0.0;
    continue;
  }

  // ---- ROI (expand seed bbox by R = disp_raw * max_disp_steps) ----
  const int R_raw = (max_disp_steps > 0) ? (disp_raw * max_disp_steps) : std::max(dim_x, dim_y);
  const int R = std::min(R_raw, roi_cap);

  int r0 = std::max(0,         min_r - R);
  int r1 = std::min(dim_y - 1, max_r + R);
  int c0 = std::max(0,         min_c - R);
  int c1 = std::min(dim_x - 1, max_c + R);

  // Tighten vertically using row spans vs. [c0,c1] to clip sea
  while (r0 <= r1 && (row_last_land[r0] < c0 || row_first_land[r0] > c1)) ++r0;
  while (r1 >= r0 && (row_last_land[r1] < c0 || row_first_land[r1] > c1)) --r1;
  if (r0 > r1) { result[sp] = 0.0; continue; }

  // Tighten horizontally using column spans vs. [r0,r1]  to clip sea
  while (c0 <= c1 && (col_last_land[c0] < r0 || col_first_land[c0] > r1)) ++c0;
  while (c1 >= c0 && (col_last_land[c1] < r0 || col_first_land[c1] > r1)) --c1;
  if (c0 > c1) { result[sp] = 0.0; continue; }

  const int W = c1 - c0 + 1;
  const int H = r1 - r0 + 1;

  // ---- Adjust the size of the masks by species ROI ----
  const std::size_t Nroi = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);

  if (Nroi > transit.capacity()) {
    transit.reserve(Nroi);
    habitat.reserve(Nroi);
    steps.reserve(Nroi);
  }
  transit.resize(Nroi);
  habitat.resize(Nroi);
  steps.resize(Nroi);

  frontier.clear(); next.clear();
  frontier.reserve(Nroi); next.reserve(Nroi);

  // ---- Fast-path for large dispersal ----
  if (disp_raw >= universal_disp_thres) {
    double total = 0.0;
    for (int r = r0; r <= r1; ++r) {
      std::size_t g = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
      + static_cast<std::size_t>(c0);
      for (int c = c0; c <= c1; ++c, ++g) {
        if (!LM[g]) continue; // skip sea cells
        const int hx = X_h_of_cell[g];
        if (hx != -1 && suitable_h_flag[hx]) {
          total += 1.0;
        }
      }
    }
    result[sp] = total;
    continue;
  }

  // ---- Build transit/habitat masks within ROI (with sea skipping) ----
  std::size_t roiH = 0;
  {
    for (int r = r0; r <= r1; ++r) {
      // Row-local starting indices
      const std::size_t rowOff = static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W);

      // If the row has no land overlapping [c0,c1], just zero and continue
      if (row_last_land[r] < c0 || row_first_land[r] > c1) {
        // zero the whole row slice once
        std::fill_n(&habitat[rowOff], W, 0u);
        std::fill_n(&transit[rowOff], W, 0u);
        std::fill_n(&steps[rowOff],   W, STEP_MAX);
        continue;
      }

      // Pre-zero entire row slice (only set non-zeros in the land span)
      std::fill_n(&habitat[rowOff], W, 0u);
      std::fill_n(&transit[rowOff], W, 0u);
      std::fill_n(&steps[rowOff],   W, STEP_MAX);

      // Work only where this row actually has land within the ROI
      const int cc0 = std::max(c0, row_first_land[r]);
      const int cc1 = std::min(c1, row_last_land[r]);

      // Global index for the first land cell in this row's ROI
      std::size_t g = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(cc0);

      // Row-local index for the first land cell
      std::size_t idx = rowOff + static_cast<std::size_t>(cc0 - c0);

      for (int c = cc0; c <= cc1; ++c, ++idx, ++g) {

        uint8_t hmask = 0u, tmask = 0u; // Masks

        const int hx = X_h_of_cell[g];
        if (hx != -1 && suitable_h_flag[hx]) {
          hmask = 1u;
          tmask = 1u; // X counts as transit too
        }
        if (!tmask) {
          const int he = E_h_of_cell[g];
          if (he != -1 && suitable_h_flag[he]) {
            tmask = 1u; // Suitable E as transit only
          }
        }

        habitat[idx] = hmask;
        transit[idx] = tmask;

        // Early out if no suitable habitat
        roiH += hmask;
      }
    }
  }

  // Early out if no suitable habitat
  if (roiH == 0) {
    result[sp] = 0.0;
    continue;
  }

  // ---- Seed frontier ----
  frontier.clear();
  {
    for (std::size_t k : seed_idxs) {
      const int r = cell_r[k];
      const int c = cell_c[k];
      if (r < r0 || r > r1 || c < c0 || c > c1) continue;

      // Local index
      const std::size_t idx =
        static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W)
        + static_cast<std::size_t>(c - c0);

      // global index
      const std::size_t global_k =
      static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(c);

      if (!LM[global_k]) continue;     // Check just in case
      if (!transit[idx]) continue;     // must be transit
      if (steps[idx] == 0) continue;   // avoid duplicate seeds

      steps[idx] = 0;
      frontier.push_back(idx);
    }
  }

  if (frontier.empty()) {
    result[sp] = 0.0;
    continue;
  }

  // ---- BFS with on-the-fly counting ----
  const auto& offsets = offsets_cache[disp];
  double g_sp = 0.0;

  for (step_t d = 0;
       !frontier.empty() &&
         (max_disp_steps <= 0 || d < static_cast<step_t>(max_disp_steps));
       ++d)
  {
    next.clear();

    for (std::size_t i : frontier) {

      const std::size_t row = i / static_cast<std::size_t>(W);
      const std::size_t col = i - row * static_cast<std::size_t>(W);
      const int r = r0 + static_cast<int>(row);
      const int c = c0 + static_cast<int>(col);

      for (const auto& off : offsets) {
        const int nr = r + off.first;
        const int nc = c + off.second;

        if (static_cast<unsigned>(nr - r0) >= static_cast<unsigned>(H) ||
            static_cast<unsigned>(nc - c0) >= static_cast<unsigned>(W)) {
          continue;
        }

        const std::size_t j =
          static_cast<std::size_t>(nr - r0) * static_cast<std::size_t>(W)
          + static_cast<std::size_t>(nc - c0);

        if (!transit[j]) continue;

        const step_t next_step = static_cast<step_t>(d + 1);
        if (steps[j] <= next_step) continue;

        steps[j] = next_step;

        if (habitat[j]) { g_sp += 1.0; }
        next.push_back(j);

      }
    }
    frontier.swap(next);
  }
  result[sp] = g_sp;
}

#endif // _OPENMP

return result;
}
