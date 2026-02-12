#ifndef EBREL_BUILDER_H
#define EBREL_BUILDER_H

#include <vector>
#include <cstdint>

#include "run_ebrel.h"

// Fully prepared inputs + derived quantities for an EBREL run.
struct BuiltEBREL {
  RunEBRELInput in;
  RunEBRELOptions opt;
  double sigma = 0.0;          // derived as 1/mean(D>0)
  std::vector<double> W;       // optional: n_h * n_cells, habitat-major
};

// Build a complete RunEBRELInput (including caches + species plan),
// derive sigma, and optionally precompute proposal weights W.
BuiltEBREL build_ebrel(
    std::vector<int8_t>  E,    // size: n_cells
    std::vector<double>  C,    // size: n_h * n_cells
    std::vector<double>  SD,   // size: n_s * n_cells
    std::vector<int>     D,    // size: n_s
    std::vector<double>  SxH,  // size: n_h * n_s (idx = h*n_s + s)
    std::vector<double>  O,    // size: n_s
    std::vector<uint8_t> LM,   // size: n_cells
    int dim_x,
    int dim_y,
    int n_h,
    int n_s,
    double sentinel,
    double sigma_in,
    int universal_disp_thres = 20,
    int max_disp_steps = 10,
    int roi_cap = 100,
    bool precompute_W = true
);

#endif // EBREL_BUILDER_H
