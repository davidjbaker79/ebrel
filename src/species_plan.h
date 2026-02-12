//-------------------- Pre-compute species plan of modelling --------------------

#ifndef SPECIES_PLAN_H
#define SPECIES_PLAN_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <cmath>

// ---- Seeds stored as global cell indices ----
using SeedT = uint32_t;

// ---- Plan for each species, stored as structure-of-arrays ----
struct SpeciesPlan {
  int n_s = 0;
  int n_h = 0;

  // per species ROI + settings
  std::vector<uint8_t> active;     // 0/1
  std::vector<uint8_t> fast_path;  // 0/1 (disp_raw >= universal_disp_thres)
  std::vector<int16_t> disp;       // min(disp_raw, universal_disp_thres)
  std::vector<int32_t> r0, r1, c0, c1;
  std::vector<int32_t> W, H;

  // Target as number of cells to gain
  std::vector<int32_t> O_n;

  // suitability flags: [sp*n_h + h] -> 0/1
  std::vector<uint8_t> suitable;   // size n_s*n_h

  // flat seed pool: seed_start[sp], seed_len[sp] index into seed_pool
  std::vector<uint32_t> seed_start;
  std::vector<uint32_t> seed_len;
  std::vector<SeedT>    seed_pool;

};

#endif  // SPECIES_PLAN_H

SpeciesPlan build_species_plan(
    const std::vector<double>& SD,
    const std::vector<double>& SxH,
    const std::vector<int>&    D,
    const std::vector<double>& O,
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
);
