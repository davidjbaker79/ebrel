//------------------------------ Dispersal-aware G function --------------------

#ifndef DISPERSAL_UTILS_H
#define DISPERSAL_UTILS_H

#include <vector>
#include <cstdint>

/* Function to calculate G accounting for dispersal through the landscape of habitat patches
 * - uses BFS approach for speed
 */
std::vector<double> compute_G(const std::vector<double>& X,
                 const std::vector<double>& E,
                 const std::vector<double>& SD,
                 const std::vector<double>& SxH,
                 const std::vector<int>& D,
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
                 const std::vector<int>& col_last_land);

#endif  // DISPERSAL_UTILS_H
