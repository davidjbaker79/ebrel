//------------------------------ Dispersal-aware G function --------------------

#ifndef DISPERSAL_UTILS_H
#define DISPERSAL_UTILS_H

#include <vector>
#include <cstdint>

// Function to calculate G accounting for dispersal through the landscape of habitat patches
std::vector<double> compute_G(const std::vector<double>& X,
                 const std::vector<double>& E,
                 const std::vector<double>& SD,
                 const std::vector<double>& SxH,
                 const std::vector<int>& D,
                 int n_h,
                 int dim_x,
                 int dim_y,
                 int max_disp_thres,
                 int disp_boundary);

#endif  // DISPERSAL_UTILS_H
