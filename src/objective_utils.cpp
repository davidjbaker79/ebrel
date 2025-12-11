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

  // For Euclidean distance
  inline int row_of(std::size_t tile, int dim_x) {
    return static_cast<int>(tile / static_cast<std::size_t>(dim_x));
  }

  // For Euclidean distance
  inline int col_of(std::size_t tile, int dim_x) {
    return static_cast<int>(tile % static_cast<std::size_t>(dim_x));
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
double compute_F1(const std::vector<double>& X,
                  const std::vector<double>& C,
                  int n_h, int dim_x, int dim_y) {

  const int cells = n_cells(dim_x, dim_y);
  double f1 = 0.0;

  for (int h = 0; h < n_h; ++h) {
    const std::size_t base = static_cast<std::size_t>(h) * static_cast<std::size_t>(cells);
    for (int tile = 0; tile < cells; ++tile) {
      std::size_t idx = base + static_cast<std::size_t>(tile);
      const double x = X[idx];
      if (x != 0.0) f1 += x * C[idx];  // skip zero in X as these will dominate
    }
  }
  return f1;
}

// F2: total = w_xx * (X–X sum once) + w_ex * (E–X sum once)
//  - original formulation used ordered pairwise sum, but this weights xx more than xe
double compute_F2(const std::vector<double>& X,
                  const std::vector<std::vector<std::size_t>>& Etiles_per_h,
                  int n_h, int dim_x, int dim_y) {

  const int cells = n_cells(dim_x, dim_y);
  double total = 0.0;

  // For each habitat
  for (int h = 0; h < n_h; ++h) {
    const std::size_t base = static_cast<std::size_t>(h)* static_cast<std::size_t>(cells);

    // Precomputed E tiles for this habitat
    const std::vector<std::size_t>& Etiles = Etiles_per_h[static_cast<std::size_t>(h)];

    // Collect X tiles for this habitat under the current configuration
    std::vector<std::size_t> Xtiles;
    Xtiles.reserve(cells / 8); // or reserve(Etiles.size()) if you like
    for (int tile = 0; tile < cells; ++tile) {
      std::size_t idx = base + static_cast<std::size_t>(tile);
      if (X[idx] == 1.0) {
        Xtiles.push_back(tile);
      }
    }

    const double xx_once = sum_pairwise_self_unordered(Xtiles, dim_x);  // X–X
    const double ex_once = sum_pairwise_cross(Etiles, Xtiles, dim_x); // E–X

    total += xx_once + ex_once;
  }

  return total;
}


// // F(X) = alpha * F1 + beta * F2
// double compute_F(const std::vector<double>& X,
//                  const std::vector<double>& C,
//                  const std::vector<double>& E,
//                  double alpha, double beta,
//                  int n_h, int dim_x, int dim_y) {
//   const double f1 = compute_F1(X, C, n_h, dim_x, dim_y);
//   const double f2 = compute_F2(X, E, n_h, dim_x, dim_y);
//   return alpha * f1 + beta * f2;
// }

// ---- Compute H
HResult compute_H(const std::vector<double>& X,
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
                  const std::vector<int>& row_first_land,
                  const std::vector<int>& row_last_land,
                  const std::vector<int>& col_first_land,
                  const std::vector<int>& col_last_land,
                  const std::vector<int>& E_h_of_cell,
                  const std::vector<std::vector<std::size_t>>& Etiles_per_h,
                  const std::vector<int>& cell_r,
                  const std::vector<int>& cell_c,
                  const std::vector<SpeciesDispData>& species_info)
{

  const int cells = n_cells(dim_x, dim_y);

  // ---- F1 and F2 -----
  const double f1 = compute_F1(X, C, n_h, dim_x, dim_y);
  const double f2 = compute_F2(X, Etiles_per_h, n_h, dim_x, dim_y);

  // ---- Compute G ----
  std::vector<double> G;
  G = compute_G(X, n_h, n_s,
                dim_x, dim_y,
                universal_disp_thres, max_disp_steps, roi_cap,
                LM,
                row_first_land, row_last_land,
                col_first_land, col_last_land,
                cell_r, cell_c,
                E_h_of_cell,
                species_info);

  // ---- Shortfall per species relative to target o_i * m_i ----
  std::vector<double> g(n_s, 0.0);

  // ---- Create a habitat mask for species so that cells that are
  // converted and no longer suitable are discounted.
  for (int sp = 0; sp < n_s; ++sp) {

    // Per-species suitability flags ----
    const SpeciesDispData& S = species_info[sp];

    // Check whether the species can go anywhere
    if (!S.active) {
      // No seeds / empty ROI -> no shortfall contribution
      g[sp] = 0.0;
      continue;
    }

    // X mask by species ----
    std::vector<double> sp_sd_x_mask(static_cast<std::size_t>(cells), 1u);
    for (int h = 0; h < n_h; ++h) {
      if (!S.suitable_h_flag[h]) {
        // Cells in habitat h that are "on" in X become unsuitable
        const std::size_t base_h = static_cast<std::size_t>(h) * static_cast<std::size_t>(cells);
        for (std::size_t k = 0; k < static_cast<std::size_t>(cells); ++k) {
          if (X[base_h + k] > 0.0) {
            sp_sd_x_mask[k] = 0u;
          }
        }
      }
    }

    // m_i: initial occupied tiles for species i (from SD)
    int m_i = 0;
    for (std::size_t idx : S.seed_idxs) {
      if (sp_sd_x_mask[idx] == 1u) {
        ++m_i;
      }
    }

    const double target = O[static_cast<std::size_t>(sp)] * static_cast<double>(m_i);
    const double Gi = G[static_cast<std::size_t>(sp)];
    // Absolute targets
    //g[static_cast<std::size_t>(sp)] = std::max(0.0, target - Gi);
    // Proportional targets
    g[static_cast<std::size_t>(sp)] = (target > 0.0) ? std::max(0.0, (target - Gi) / target): 0.0;

    // std::cout << "sp=" << sp
    //           << " m_i=" << m_i
    //           << " target=" << target
    //           << " G=" << Gi
    //           << " shortfall=" << std::max(0.0, target - Gi)
    //           << std::endl;

  }

  const double gx_val = std::accumulate(g.begin(), g.end(), 0.0);
  const double Fx_val = alpha_scaled * f1 + beta_scaled * f2;
  const double H_val  = Fx_val + gamma_scaled * gx_val;

  return { H_val, Fx_val, gx_val, f1, f2, g };
}

