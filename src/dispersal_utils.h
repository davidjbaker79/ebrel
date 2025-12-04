//------------------------------ Dispersal-aware G function --------------------

#ifndef DISPERSAL_UTILS_H
#define DISPERSAL_UTILS_H

#include <vector>
#include <cstdint>

// Structure for storing fragment of ROI
struct Fragment {
  int r0, r1, c0, c1;                  // ROI bounds for this fragment
  std::vector<std::size_t> seed_idxs;  // seeds (global cell indices)
};

/* Precalculate information on dispersal for each species */
struct SpeciesDispData {
  bool active;                // false if no seeds or ROI empty
  int sp;                     // species index
  int disp_raw;               // D[sp]
  int disp;                   // min(disp_raw, universal_disp_thres)
  std::vector<uint8_t> suitable_h_flag; // size n_h: 0/1
  // Seeds (global cell indices) satisfying SD==1 && E suitable
  std::vector<std::size_t> seed_idxs;
  // Collect the ROI fragments
  std::vector<Fragment>    fragments;   // per-fragment ROIs
  // ROI bounds for this species (clipped to land)
  int r0, r1, c0, c1;         // inclusive bounds
  // Fast path flag: use "global dispersal" branch only
  bool use_fast_path;
};

// Pre-compute species dispersal information
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
    const std::vector<int>& cell_c,
    int cluster_gap_cells
);

/* Function to calculate G accounting for dispersal through the landscape of habitat patches
 * - uses BFS approach for speed
 */
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
    const std::vector<SpeciesDispData>& species_info);

#endif  // DISPERSAL_UTILS_H
