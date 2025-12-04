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

/* Pre-compute the species-specific dispersal information needed within
 * compute_G that remains constant across iterations */
std::vector<SpeciesDispData> precompute_species_data(
    const std::vector<double>& SD,
    const std::vector<double>& SxH,
    const std::vector<int>&    D,
    int n_h, int n_s,
    int dim_x, int dim_y,
    int universal_disp_thres,
    int max_disp_steps,
    int roi_cap,
    const std::vector<int>& row_first_land,
    const std::vector<int>& row_last_land,
    const std::vector<int>& col_first_land,
    const std::vector<int>& col_last_land,
    const std::vector<int>& E_h_of_cell,
    const std::vector<int>& cell_r,
    const std::vector<int>& cell_c
) {

  const int cells = n_cells(dim_x, dim_y);
  std::vector<SpeciesDispData> out(n_s);

  for (int sp = 0; sp < n_s; ++sp) {

    SpeciesDispData S;
    S.sp = sp;
    S.active = false;

    // Suitable_h_flag
    S.suitable_h_flag.assign(n_h, 0u);
    const std::size_t base_h = static_cast<std::size_t>(sp) * n_h;
    for (int h = 0; h < n_h; ++h) {
      if (SxH[base_h + h] > 0.0) S.suitable_h_flag[h] = 1u;
    }

    S.disp_raw = D[sp];
    S.disp = std::min(S.disp_raw, universal_disp_thres);

    // --- Seeds + bbox
    const std::size_t base_sd = static_cast<std::size_t>(sp) * static_cast<std::size_t>(cells);
    int min_r = dim_y, max_r = -1, min_c = dim_x, max_c = -1;

    S.seed_idxs.clear();
    for (int k = 0; k < cells; ++k) {

      std::size_t kk = static_cast<std::size_t>(k);

      if (SD[base_sd + kk] != 1.0) continue;

      // This condition could be set - not for now
      //const int ehx = E_h_of_cell[k]; // precomputed once
      //if (ehx == -1 || !S.suitable_h_flag[ehx]) continue;

      S.seed_idxs.push_back(kk);

      const int r = cell_r[k];
      const int c = cell_c[k];
      if (r < min_r) min_r = r;
      if (r > max_r) max_r = r;
      if (c < min_c) min_c = c;
      if (c > max_c) max_c = c;

    }

    if (S.seed_idxs.empty()) {
      out[sp] = std::move(S);
      continue; // inactive species
    }

    // ROI expansion (depends only on disp_raw, max_disp_steps, roi_cap)
    const int R_raw = (max_disp_steps > 0)
      ? (S.disp_raw * max_disp_steps)
      : std::max(dim_x, dim_y);
    const int R = std::min(R_raw, roi_cap);

    int r0 = std::max(0,         min_r - R);
    int r1 = std::min(dim_y - 1, max_r + R);
    int c0 = std::max(0,         min_c - R);
    int c1 = std::min(dim_x - 1, max_c + R);

    // Clip ROI using land spans – exactly as in your function
    while (r0 <= r1 && (row_last_land[r0] < c0 || row_first_land[r0] > c1)) ++r0;
    while (r1 >= r0 && (row_last_land[r1] < c0 || row_first_land[r1] > c1)) --r1;
    if (r0 > r1) {
      out[sp] = std::move(S);
      continue;
    }

    while (c0 <= c1 && (col_last_land[c0] < r0 || col_first_land[c0] > r1)) ++c0;
    while (c1 >= c0 && (col_last_land[c1] < r0 || col_first_land[c1] > r1)) --c1;
    if (c0 > c1) {
      out[sp] = std::move(S);
      continue;
    }

    S.r0 = r0; S.r1 = r1; S.c0 = c0; S.c1 = c1;

    S.use_fast_path = (S.disp_raw >= universal_disp_thres);
    S.active = true;

    out[sp] = std::move(S);
  }
  return out;
}

/* Compute G(x): per-species count of reachable newly created suitable habitat (X),
 * using E (existing) and suitable-X as transit */
std::vector<double> compute_G(
    const std::vector<double>& X,
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
    const std::vector<int>& cell_r,
    const std::vector<int>& cell_c,
    const std::vector<int>& E_h_of_cell,
    const std::vector<SpeciesDispData>& species_info)
{
  const int cells = n_cells(dim_x, dim_y);
  std::vector<double> result(n_s, 0.0);

  // offsets cache (you *can* precompute this once and pass in too if desired)
  std::vector<std::vector<std::pair<int,int>>> offsets_cache(universal_disp_thres + 1);
  for (int d = 0; d <= universal_disp_thres; ++d) {
    offsets_cache[d] = make_offsets(d);
  }

  // One-hot for X, which changes each iteration
  std::vector<int> X_h_of_cell(cells, -1);
  for (int h = 0; h < n_h; ++h) {
    const std::size_t base = static_cast<std::size_t>(h) * static_cast<std::size_t>(cells);
    for (int k = 0; k < cells; ++k) {
      std::size_t kk = static_cast<std::size_t>(k);
      if (X[base + kk] == 1.0) X_h_of_cell[k] = h;
    }
  }

  using step_t = uint16_t;
  const step_t STEP_MAX = std::numeric_limits<step_t>::max();

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
  // Thread-local scratch
  std::vector<uint8_t> transit;
  std::vector<uint8_t> habitat;
  std::vector<step_t>  steps;
  std::vector<std::size_t> frontier;
  std::vector<std::size_t> next;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (int sp = 0; sp < n_s; ++sp) {

    const SpeciesDispData& S = species_info[sp];
    if (!S.active) {
      result[sp] = 0.0;
      continue;
    }

    //const int disp_raw = S.disp_raw;
    const int disp = S.disp;

    const int r0 = S.r0;
    const int r1 = S.r1;
    const int c0 = S.c0;
    const int c1 = S.c1;

    const int W = c1 - c0 + 1;
    const int H = r1 - r0 + 1;
    const std::size_t Nroi =
      static_cast<std::size_t>(W) * static_cast<std::size_t>(H);

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

    const auto& suitable_h_flag = S.suitable_h_flag;

    // fast path: effectively global dispersal for this species
    if (S.use_fast_path) {
      double total = 0.0;
      for (int r = r0; r <= r1; ++r) {
        std::size_t g =
          static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(c0);
        for (int c = c0; c <= c1; ++c, ++g) {
          if (!LM[g]) continue;
          const int hx = X_h_of_cell[g];
          if (hx != -1 && suitable_h_flag[hx]) {
            total += 1.0;
          }
        }
      }
      result[sp] = total;
      continue;
    }

    // Build transit/habitat masks within ROI
    std::size_t roiH = 0;
    for (int r = r0; r <= r1; ++r) {
      const std::size_t rowOff =
        static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W);

      if (row_last_land[r] < c0 || row_first_land[r] > c1) {
        // whole row is sea in this ROI
        std::fill_n(&habitat[rowOff], W, 0u);
        std::fill_n(&transit[rowOff], W, 0u);
        // no need to touch steps: they will never be used in this row
        continue;
      }

      std::fill_n(&habitat[rowOff], W, 0u);
      std::fill_n(&transit[rowOff], W, 0u);
      // note: we will initialise steps[] only where transit==1

      const int cc0 = std::max(c0, row_first_land[r]);
      const int cc1 = std::min(c1, row_last_land[r]);

      std::size_t g =
        static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(cc0);

      std::size_t idx = rowOff + static_cast<std::size_t>(cc0 - c0);

      for (int c = cc0; c <= cc1; ++c, ++idx, ++g) {

        uint8_t hmask = 0u;
        uint8_t tmask = 0u;

        const int hx = X_h_of_cell[g];
        if (hx != -1 && suitable_h_flag[hx]) {
          hmask = 1u;
          tmask = 1u; // X counts as transit too
        }
        if (!tmask) {
          const int he = E_h_of_cell[g];
          if (he != -1 && suitable_h_flag[he]) {
            tmask = 1u; // suitable E as transit only
          }
        }

        habitat[idx] = hmask;
        transit[idx] = tmask;

        if (tmask) {
          steps[idx] = STEP_MAX;  // initialise steps *only* for transit cells
        }

        roiH += hmask;
      }
    }

    if (roiH == 0) {
      result[sp] = 0.0;
      continue;
    }

    // Seed frontier
    frontier.clear();
    for (std::size_t k : S.seed_idxs) {
      const int r = cell_r[k];
      const int c = cell_c[k];
      if (r < r0 || r > r1 || c < c0 || c > c1) continue;

      const std::size_t idx =
        static_cast<std::size_t>(r - r0) * static_cast<std::size_t>(W)
        + static_cast<std::size_t>(c - c0);

      const std::size_t g =
      static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x)
        + static_cast<std::size_t>(c);

      if (!LM[g]) continue;
      if (!transit[idx]) continue;
      if (steps[idx] == 0) continue;

      steps[idx] = 0;
      frontier.push_back(idx);
    }

    if (frontier.empty()) {
      result[sp] = 0.0;
      continue;
    }

    // BFS with on-the-fly counting and early exit
    const auto& offsets = offsets_cache[disp];
    double g_sp = 0.0;
    std::size_t remaining = roiH;

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

          if (habitat[j]) {
            g_sp += 1.0;
            if (--remaining == 0) {
              // all reachable habitat counted; can stop BFS
              frontier.clear();
              break;
            }
          }

          next.push_back(j);
        }

        if (frontier.empty()) break; // early break from outer loop too
      }

      frontier.swap(next);
    }
    result[sp] = g_sp;
  }
  } // end omp parallel

  return result;
}

