//--------------------------- Ebrel setup functions ----------------------------

// These functions are for creating an ebrel object with validation. The ebrel
// object is then suppled to the run_ebrel function.

// Memory layout note for reference:
// - Habitat-major 3D fields (E, C): idx = h * (dim_x*dim_y) + cell
// - Species-major 3D field (SD):     idx = s * (dim_x*dim_y) + cell
// - SxH (n_h * n_s):                 idx = h * n_s + s
// - O is species-level:              size == n_s

#ifndef EBREL_SETUP_H
#define EBREL_SETUP_H

#include <vector>
#include <stdexcept>
#include <string>
#include <cstddef>

// Check that each cell has at most one habitat present
// - at present this is binary 0/1 but may relax later to accommodate partial habitat cover
void validate_existing_habitat(const std::vector<double>& E,
                               int dim_x,
                               int dim_y,
                               int n_h);

// Validate input dimensions for creating the Ebrel class object
inline void validate_dimensions(const std::vector<double>& E,
                                const std::vector<double>& C,
                                const std::vector<double>& SD,
                                const std::vector<int>& D,
                                const std::vector<double>& SxH,
                                const std::vector<double>& O,
                                int dim_x,
                                int dim_y,
                                int n_h,
                                int n_s) {
  if (dim_x <= 0 || dim_y <= 0 || n_h <= 0 || n_s <= 0) {
    throw std::runtime_error("dim_x, dim_y, n_h, and n_s must all be positive.");
  }

  const std::size_t cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  const std::size_t Hsz   = cells * static_cast<std::size_t>(n_h);
  const std::size_t Ssz   = cells * static_cast<std::size_t>(n_s);
  const std::size_t SxHsz = static_cast<std::size_t>(n_h) * static_cast<std::size_t>(n_s);

  if (E.size()  != Hsz)   throw std::runtime_error("E size mismatch: got " + std::to_string(E.size())  + ", expected " + std::to_string(Hsz));
  if (C.size()  != Hsz)   throw std::runtime_error("C size mismatch: got " + std::to_string(C.size())  + ", expected " + std::to_string(Hsz));
  if (SD.size() != Ssz)   throw std::runtime_error("SD size mismatch: got " + std::to_string(SD.size()) + ", expected " + std::to_string(Ssz));
  if (D.size()  != static_cast<std::size_t>(n_s))
    throw std::runtime_error("D size mismatch: got " + std::to_string(D.size())  + ", expected " + std::to_string(n_s));
  if (SxH.size()!= SxHsz) throw std::runtime_error("SxH size mismatch: got " + std::to_string(SxH.size()) + ", expected " + std::to_string(SxHsz));
  if (O.size()  != static_cast<std::size_t>(n_s))
    throw std::runtime_error("O size mismatch: got " + std::to_string(O.size()) + ", expected " + std::to_string(n_s) + " (species-level).");
}

// Create unavailable mask from C
[[nodiscard]] std::vector<double> get_unavailable_mask(const std::vector<double>& C,
                                                       int n_h,
                                                       int dim_x,
                                                       int dim_y,
                                                       double sentinel = 1e10);

// Build Ebrel class object outputs
inline void create_ebrel_class_object(const std::vector<double>& E,
                                      const std::vector<double>& C,
                                      const std::vector<double>& SD,
                                      const std::vector<int>& D,
                                      const std::vector<double>& SxH,
                                      const std::vector<double>& O,
                                      int dim_x,
                                      int dim_y,
                                      int n_h,
                                      int n_s,
                                      double sentinel,
                                      std::vector<double>& U_out) {
  validate_dimensions(E, C, SD, D, SxH, O, dim_x, dim_y, n_h, n_s);
  validate_existing_habitat(E, dim_x, dim_y, n_h);
  U_out = get_unavailable_mask(C, n_h, dim_x, dim_y, sentinel);
}

#endif  // EBREL_SETUP_H
