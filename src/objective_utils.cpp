
//-------------------------- Ebrel objective functions -------------------------

#include "objective_utils.h"
#include "dispersal_utils.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>  // std::runtime_error, std::invalid_argument
#include <string>     // std::to_strings
//#include <iostream>

//-------------------------- Local helper functions ----------------------------
namespace {

  // Get number of cells
  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  }

  inline std::size_t suitable_idx(const SpeciesPlan& P, int sp, int h) {
    return static_cast<std::size_t>(sp) * static_cast<std::size_t>(P.n_h) + static_cast<std::size_t>(h);
  }

  inline bool is_suitable(const SpeciesPlan& P, int sp, int h) {
    return P.suitable[suitable_idx(P, sp, h)] != 0u;
  }

  // SeedT is presumably uint32_t or similar; adapt if needed.
  inline const SeedT* seed_begin(const SpeciesPlan& P, int sp) {
    return P.seed_pool.data() + P.seed_start[sp];
  }
  inline const SeedT* seed_end(const SpeciesPlan& P, int sp) {
    return P.seed_pool.data() + P.seed_start[sp] + P.seed_len[sp];
  }

  // For Euclidean distance
  inline int row_of(std::size_t tile, int dim_x) {
    return static_cast<int>(tile / static_cast<std::size_t>(dim_x));
  }

  // For Euclidean distance
  inline int col_of(std::size_t tile, int dim_x) {
    return static_cast<int>(tile % static_cast<std::size_t>(dim_x));
  }

  // Compute discounted seeds
  inline int compute_m_i_from_plan(
      int sp,
      const SpeciesPlan& species_plan,
      const std::vector<int8_t>& X
  ) {
    int m_i = 0;

    const SeedT* b = seed_begin(species_plan, sp);
    const SeedT* e = seed_end(species_plan, sp);

    for (const SeedT* it = b; it != e; ++it) {
      const std::size_t k = static_cast<std::size_t>(*it);

      const int hx = X[k];
      if (hx < 0) {
        ++m_i; // not converted -> seed remains
      } else {
        if (is_suitable(species_plan, sp, hx)) ++m_i; // converted but still suitable -> keep
        // else discounted
      }
    }
    return m_i;
  }


  // Sum of Euclidean distances for unordered pairs within one set (i<j).
  // - can multiply by 2 to obtain the ordered-pair sum to match JAE
  double sum_pairwise_self_unordered(const std::vector<std::size_t>& tiles, int dim_x) {
    const std::size_t n = tiles.size();
    if (n < 2) return 0.0;
    double s = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
      const int ri = row_of(tiles[i], dim_x), ci = col_of(tiles[i], dim_x);
      for (std::size_t j = i + 1; j < n; ++j) {
        const int rj = row_of(tiles[j], dim_x), cj = col_of(tiles[j], dim_x);
        const double dr = static_cast<double>(ri - rj);
        const double dc = static_cast<double>(ci - cj);
        s += std::sqrt(dr*dr + dc*dc);
      }
    }
    return s;
  }

  // Sum of Euclidean distances across two sets (e.g. E × X0).
  double sum_pairwise_cross(const std::vector<std::size_t>& A,
                            const std::vector<std::size_t>& B,
                            int dim_x) {
    if (A.empty() || B.empty()) return 0.0;
    double s = 0.0;
    for (std::size_t ia = 0; ia < A.size(); ++ia) {
      const int ra = row_of(A[ia], dim_x), ca = col_of(A[ia], dim_x);
      for (std::size_t ib = 0; ib < B.size(); ++ib) {
        const int rb = row_of(B[ib], dim_x), cb = col_of(B[ib], dim_x);
        const double dr = static_cast<double>(ra - rb);
        const double dc = static_cast<double>(ca - cb);
        s += std::sqrt(dr*dr + dc*dc);
      }
    }
    return s;
  }

}

// -------------------- Objective functions ------------------------------------

// F1: sum_{h,cell} X * C  (habitat-major layout: idx = h*cells + tile)
double compute_F1(const std::vector<int8_t>& X,
                  const std::vector<double>& C,
                  int n_h, int dim_x, int dim_y) {

  const int cells = n_cells(dim_x, dim_y);
  double f1 = 0.0;

  for (int cell = 0; cell < cells; ++cell) {
    const int h = X[cell];
    if (h >= 0) {
      const std::size_t idx =
        static_cast<std::size_t>(h) * static_cast<std::size_t>(cells)
      + static_cast<std::size_t>(cell);
      f1 += C[idx];
    }
  }

  return f1;
}

// F2: total = w_xx * (X–X sum once) + w_ex * (E–X sum once)
//  - original formulation used ordered pairwise sum, but this weights xx more than xe
double compute_F2(const std::vector<int8_t>& X,
                  const std::vector<std::vector<std::size_t>>& Etiles_per_h,
                  int n_h, int dim_x, int dim_y) {

  const int cells = n_cells(dim_x, dim_y);

  // Build X tiles per habitat in one pass over cells
  std::vector<std::vector<std::size_t>> Xtiles_per_h(static_cast<std::size_t>(n_h));
  // Optional: rough reserve heuristic (tune or drop)
  for (int h = 0; h < n_h; ++h) Xtiles_per_h[static_cast<std::size_t>(h)].reserve(cells / 16);

  for (int tile = 0; tile < cells; ++tile) {
    const int h = X[static_cast<std::size_t>(tile)];
    if (h >= 0) {
      Xtiles_per_h[static_cast<std::size_t>(h)].push_back(static_cast<std::size_t>(tile));
    }
  }

  double total = 0.0;

  for (int h = 0; h < n_h; ++h) {
    const auto& Etiles = Etiles_per_h[static_cast<std::size_t>(h)];
    const auto& Xtiles = Xtiles_per_h[static_cast<std::size_t>(h)];

    if (!Xtiles.empty()) {
      total += sum_pairwise_self_unordered(Xtiles, dim_x);     // X–X
      if (!Etiles.empty()) total += sum_pairwise_cross(Etiles, Xtiles, dim_x); // E–X
    }
  }

  return total;
}

// ---- Compute H
HResult compute_H(const std::vector<int8_t>& X,
                  const std::vector<double>& C,
                  const std::vector<double>& O,
                  const std::vector<double>& SxH,
                  const std::vector<int>& D,
                  double alpha_scaled,
                  double beta_scaled,
                  double gamma_scaled,
                  int n_h,
                  int n_s,
                  int dim_x,
                  int dim_y,
                  int universal_disp_thres,
                  int max_disp_steps,
                  int roi_cap,
                  const std::vector<uint8_t>& LM,
                  const std::vector<int16_t>& row_first_land,
                  const std::vector<int16_t>& row_last_land,
                  const std::vector<int16_t>& col_first_land,
                  const std::vector<int16_t>& col_last_land,
                  const std::vector<int8_t>& E,
                  const std::vector<std::vector<std::size_t>>& Etiles_per_h,
                  const std::vector<int>& cell_r,
                  const std::vector<int>& cell_c,
                  const RowRunsCache& rowruns_cache,
                  const SpeciesPlan& species_plan)
{

  // ---- F1 and F2 -----
  const double f1 = compute_F1(X, C, n_h, dim_x, dim_y);
  const double f2 = compute_F2(X, Etiles_per_h, n_h, dim_x, dim_y);

  // ---- Compute G ----
  std::vector<double> G = compute_G(
    X,
    LM,
    row_first_land, row_last_land,
    col_first_land, col_last_land,
    cell_r, cell_c,
    E,
    dim_x, dim_y,
    universal_disp_thres, max_disp_steps, roi_cap,
    rowruns_cache,
    species_plan);

  // ---- Shortfall per species relative to target o_i * m_i ----
  std::vector<double> g(static_cast<std::size_t>(n_s), 0.0);

  for (int sp = 0; sp < n_s; ++sp) {

    if (!species_plan.active[sp]) {
      g[static_cast<std::size_t>(sp)] = 0.0;
      continue;
    }

    // This needs to be target in n cells based on original distribution

    const int m_i = compute_m_i_from_plan(sp, species_plan, X);

    const double target = species_plan.O_n[static_cast<std::size_t>(sp)];
    const double m0   = static_cast<double>(species_plan.seed_len[static_cast<std::size_t>(sp)]);
    const double kept = static_cast<double>(m_i);
    const double lost = m0 - kept;

    // net gain = new reachable X - lost original seeds
    const double Gi = G[static_cast<std::size_t>(sp)] - lost;

    // Proportional target shortfall
    g[static_cast<std::size_t>(sp)] =
      (target > 0.0) ? std::max(0.0, (target - Gi) / target) : 0.0;

  }

  const double gx_val = std::accumulate(g.begin(), g.end(), 0.0);
  const double Fx_val = alpha_scaled * f1 + beta_scaled * f2;
  const double H_val  = Fx_val + gamma_scaled * gx_val;

  return { H_val, Fx_val, gx_val, f1, f2, g };
}
