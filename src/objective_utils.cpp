//-------------------------- Ebrel objective functions -------------------------

#include "objective_utils.h"
#include "dispersal_utils.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstddef>

// ------------------------- Local helper functions ----------------------------
namespace {

  // Get number of cells
  inline std::size_t n_cells(int dim_x, int dim_y) {
    return static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  }

  inline int row_of(std::size_t tile, int dim_x) {
    return static_cast<int>(tile / static_cast<std::size_t>(dim_x));
  }

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
  const std::size_t cells = n_cells(dim_x, dim_y);
  double f1 = 0.0;

  for (int h = 0; h < n_h; ++h) {
    const std::size_t base = static_cast<std::size_t>(h) * cells;
    for (std::size_t tile = 0; tile < cells; ++tile) {
      const std::size_t idx = base + tile;
      const double x = X[idx];
      if (x != 0.0) f1 += x * C[idx];  // skip zero in X as these will dominate
    }
  }
  return f1;
}

// F2: total = w_xx * (X–X sum once) + w_ex * (E–X sum once)
//  - original formulation used ordered pairwise sum, but this weights xx more than xe
double compute_F2(const std::vector<double>& X,
                  const std::vector<double>& E,
                  int n_h, int dim_x, int dim_y) {
  const std::size_t cells = n_cells(dim_x, dim_y);
  double total = 0.0;

  for (int h = 0; h < n_h; ++h) {
    const std::size_t base = static_cast<std::size_t>(h) * cells;

    // Collect tile indices (integers) instead of allocating 2D coords
    std::vector<std::size_t> Xtiles;
    std::vector<std::size_t> Etiles;
    Xtiles.reserve(cells / 8); // heuristic; OK to drop
    Etiles.reserve(cells / 8);

    for (std::size_t tile = 0; tile < cells; ++tile) {
      const std::size_t idx = base + tile;
      if (X[idx] == 1.0) Xtiles.push_back(tile);
      if (E[idx] == 1.0) Etiles.push_back(tile);
    }

    const double xx_once = sum_pairwise_self_unordered(Xtiles, dim_x); // unordered (i<j) only
    const double ex_once = sum_pairwise_cross(Etiles, Xtiles, dim_x);  // E×X once

    // Calculate total (xx_once using unordered pairs; 2 * xx_once to get ordered)
    total += 1 * xx_once + 1 * ex_once;
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

// H(X) = alpha*f1 + beta*f2 + gamma*sum_i
// - also returns all the other objective function statistics
HResult compute_H(const std::vector<double>& X,
                  const std::vector<double>& C,
                  const std::vector<double>& E,
                  const std::vector<double>& O,
                  const std::vector<double>& SD,
                  const std::vector<double>& SxH,
                  const std::vector<int>& D,
                  double alpha_scaled,
                  double beta_scaled,
                  double gamma_scaled,
                  int n_h, int dim_x, int dim_y,
                  int max_disp_thres,
                  int disp_boundary) {
  const double f1 = compute_F1(X, C, n_h, dim_x, dim_y);
  const double f2 = compute_F2(X, E, n_h, dim_x, dim_y);

  // G computed using new dispersal BFS model
  std::vector<double> G = compute_G(X, E, SD, SxH, D,
                                    n_h, dim_x, dim_y,
                                    max_disp_thres, disp_boundary);

  const int n_species = static_cast<int>(G.size());
  const std::size_t cells = n_cells(dim_x, dim_y);

  std::vector<double> g(n_species, 0.0);  // shortfall per species
  for (int i = 0; i < n_species; ++i) {
    int m_i = 0;
    const std::size_t base = static_cast<std::size_t>(i) * cells;
    for (std::size_t k = 0; k < cells; ++k) {
      if (SD[base + k] == 1.0) ++m_i;
    }
    const double target = O[static_cast<std::size_t>(i)] * static_cast<double>(m_i);
    g[static_cast<std::size_t>(i)] = std::max(0.0, target - G[static_cast<std::size_t>(i)]);
  }

  const double gx_val = std::accumulate(g.begin(), g.end(), 0.0);
  const double Fx_val = alpha_scaled * f1 + beta_scaled * f2;
  const double H_val  = Fx_val + gamma_scaled * gx_val;

  return { H_val, Fx_val, gx_val, f1, f2, g };
}
