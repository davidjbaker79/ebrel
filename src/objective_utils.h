#ifndef OBJECTIVE_UTILS_H
#define OBJECTIVE_UTILS_H

#include <vector>
#include <cmath>
#include "dispersal_utils.h"  // compute_G

// Return structure for H(x)
struct HResult {
  double H;
  double Fx;
  double gx;
  double F1;
  double F2;
  std::vector<double> g;
};

// Declare F1 objective function
double compute_F1(const std::vector<double>& X,
                  const std::vector<double>& C,
                  int n_h,
                  int dim_x,
                  int dim_y);

// Declare F2 objective function
double compute_F2(const std::vector<double>& X,
                  const std::vector<double>& E,
                  int n_h,
                  int dim_x,
                  int dim_y);

// Declare internal geometry helper functions
double euclidean_distance(const std::vector<double>& a,
                          const std::vector<double>& b);

double pairwise_distances(const std::vector<std::vector<double>>& A,
                          const std::vector<std::vector<double>>& B);

// Declare function to compute F(x)
double compute_F(const std::vector<double>& X,
                 const std::vector<double>& C,
                 const std::vector<double>& E,
                 double alpha,
                 double beta,
                 int n_h,
                 int dim_x,
                 int dim_y);

// Declare function to compute H(x)
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
                  int n_h,
                  int dim_x,
                  int dim_y,
                  int max_disp_thres,
                  int disp_boundary);



#endif // OBJECTIVE_UTILS_H
