//------------------------------ Dispersal-aware G function --------------------

#include "dispersal_utils.h"

#include <vector>
#include <queue>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <limits> //  for std::numeric_limits

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
                              const std::vector<double>& SxH, // nh x ns
                              const std::vector<int>& D,
                              int n_h, int dim_x, int dim_y,
                              int max_disp_thres,
                              int disp_boundary)
{
  const std::size_t cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  const int n_s = static_cast<int>(D.size());

  std::vector<double> result(n_s, 0.0);

  // Precompute row/col coords
  std::vector<int> cell_r(cells), cell_c(cells);
  for (std::size_t k = 0; k < cells; ++k) {
    cell_r[k] = static_cast<int>(k / static_cast<std::size_t>(dim_x));
    cell_c[k] = static_cast<int>(k % static_cast<std::size_t>(dim_x));
  }

  // Cache offsets for all distances
  std::vector<std::vector<std::pair<int,int>>> offsets_cache(max_disp_thres + 1);
  for (int d = 0; d <= max_disp_thres; ++d) {
    offsets_cache[d] = make_offsets(d);
  }

  // One-hot habitat index per cell for X and E (-1 = none)
  std::vector<int> X_h_of_cell(cells, -1), E_h_of_cell(cells, -1);
  for (int h = 0; h < n_h; ++h) {
    const std::size_t base = static_cast<std::size_t>(h) * cells;
    for (std::size_t k = 0; k < cells; ++k) {
      if (X[base + k] == 1.0) X_h_of_cell[k] = h;
      if (E[base + k] == 1.0) E_h_of_cell[k] = h;
    }
  }

#pragma omp parallel for schedule(dynamic) if (n_s > 1)
  for (int sp = 0; sp < n_s; ++sp) {

    // Per-species suitability flags
    std::vector<uint8_t> suitable_h_flag(n_h, 0u);
    const std::size_t base_h = static_cast<std::size_t>(sp) * static_cast<std::size_t>(n_h);
    for (int h = 0; h < n_h; ++h) {
      suitable_h_flag[h] = (SxH[base_h + static_cast<std::size_t>(h)] > 0.0) ? 1u : 0u;
    }

    // Species dispersal (bounded by max)
    const int disp_raw = D[sp];
    const int disp = std::min(disp_raw, max_disp_thres);

    // === Unbounded species fast-path: take precedence over everything ===
    if (disp_raw >= max_disp_thres) {
      // Count all X-cells whose habitat type is suitable for this species.
      double total = 0.0;
      for (std::size_t k = 0; k < cells; ++k) {
        const int hx = X_h_of_cell[k];
        if (hx != -1 && suitable_h_flag[hx]) total += 1.0;
      }
      result[sp] = total;
      continue;
    }

    // Buffers (kept full-sized for simplicity; work is ROI-limited below)
    std::vector<uint8_t> transit(cells, 0u);
    std::vector<uint8_t> habitat(cells, 0u);
    std::vector<uint8_t> reachable(cells, 0u);
    std::vector<uint16_t> steps(cells, std::numeric_limits<uint16_t>::max());
    std::queue<std::size_t> q;

    // --- Find seeds and their bounding box
    const std::size_t base_sd = static_cast<std::size_t>(sp) * cells;
    std::vector<std::size_t> seed_idxs;
    seed_idxs.reserve(1024);

    int min_r = dim_y, max_r = -1, min_c = dim_x, max_c = -1;
    for (std::size_t k = 0; k < cells; ++k) {
      if (SD[base_sd + k] != 1.0) continue;
      const int hx = X_h_of_cell[k];
      if (hx != -1 && !suitable_h_flag[hx]) continue; // don't seed if X is unsuitable
      seed_idxs.push_back(k);
      const int r = cell_r[k], c = cell_c[k];
      if (r < min_r) { min_r = r; }
      if (r > max_r) { max_r = r; }
      if (c < min_c) { min_c = c; }
      if (c > max_c) { max_c = c; }
    }
    if (seed_idxs.empty()) { result[sp] = 0.0; continue; }

    // --- ROI: expand seed bbox by R = disp * disp_boundary
    int R = (disp_boundary > 0) ? (disp * disp_boundary) : std::max(dim_x, dim_y);
    int r0 = std::max(0, min_r - R), r1 = std::min(dim_y - 1, max_r + R);
    int c0 = std::max(0, min_c - R), c1 = std::min(dim_x - 1, max_c + R);

    // --- Build transit/habitat masks only within ROI
    for (int r = r0; r <= r1; ++r) {
      for (int c = c0; c <= c1; ++c) {
        const std::size_t k = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(c);
        uint8_t hmask = 0u, tmask = 0u;
        const int hx = X_h_of_cell[k];
        if (hx != -1 && suitable_h_flag[hx]) { hmask = 1u; tmask = 1u; }
        if (!tmask) {
          const int he = E_h_of_cell[k];
          if (he != -1 && suitable_h_flag[he]) tmask = 1u;
        }
        habitat[k] = hmask;
        transit[k] = tmask;
      }
    }

    // --- Seed queue (only from seed_idxs)
    for (auto k : seed_idxs) {
      // seed may lie just outside ROI in degenerate cases; guard anyway
      const int r = cell_r[k], c = cell_c[k];
      if (r < r0 || r > r1 || c < c0 || c > c1) continue;
      reachable[k] = 1u;
      steps[k] = 0;
      q.push(k);
    }
    if (q.empty()) { result[sp] = 0.0; continue; }

    // --- Depth-limited BFS
    const auto& offsets = offsets_cache[disp];
    while (!q.empty()) {
      const std::size_t curr = q.front(); q.pop();
      const int r = cell_r[curr], c = cell_c[curr];

      // stop expanding beyond boundary (if set)
      if (disp_boundary > 0 && steps[curr] >= static_cast<uint16_t>(disp_boundary))
        continue;

      for (const auto& [dr, dc] : offsets) {
        const int nr = r + dr, nc = c + dc;
        if (nr < r0 || nr > r1 || nc < c0 || nc > c1) continue; // ROI guard
        const std::size_t nb = static_cast<std::size_t>(nr) * static_cast<std::size_t>(dim_x)
          + static_cast<std::size_t>(nc);
        if (!transit[nb]) continue;

        const uint16_t next_step = static_cast<uint16_t>(steps[curr] + 1u);
        if (steps[nb] <= next_step) continue; // already discovered with fewer/equal steps

        steps[nb] = next_step;
        reachable[nb] = 1u;
        q.push(nb);
      }
    }

    // --- Sum reachable habitat within ROI
    double g_sp = 0.0;
    for (int r = r0; r <= r1; ++r) {
      for (int c = c0; c <= c1; ++c) {
        const std::size_t k = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(c);
        if (reachable[k] && habitat[k]) g_sp += 1.0;
      }
    }
    result[sp] = g_sp;
  }

  return result;
}
