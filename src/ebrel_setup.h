//--------------------------- Ebrel setup functions ----------------------------

// These functions are for creating an Ebrel object with validation. The Ebrel
// object is then supplied to the run_ebrel function.

// Memory layout note for reference:
// - Habitat-major 3D fields (E, C, P): idx = h * (dim_x*dim_y) + cell
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
                                const std::vector<uint8_t>& LM,
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

  if (E.size()  != Hsz)
    throw std::runtime_error("E size mismatch: got " + std::to_string(E.size())  +
                             ", expected " + std::to_string(Hsz));
  if (C.size()  != Hsz)
    throw std::runtime_error("C size mismatch: got " + std::to_string(C.size())  +
                             ", expected " + std::to_string(Hsz));
  if (SD.size() != Ssz)
    throw std::runtime_error("SD size mismatch: got " + std::to_string(SD.size()) +
                             ", expected " + std::to_string(Ssz));
  if (D.size()  != static_cast<std::size_t>(n_s))
    throw std::runtime_error("D size mismatch: got " + std::to_string(D.size())  +
                             ", expected " + std::to_string(n_s));
  if (SxH.size()!= SxHsz)
    throw std::runtime_error("SxH size mismatch: got " + std::to_string(SxH.size()) +
                             ", expected " + std::to_string(SxHsz));
  if (O.size()  != static_cast<std::size_t>(n_s))
    throw std::runtime_error("O size mismatch: got " + std::to_string(O.size()) +
                             ", expected " + std::to_string(n_s) + " (species-level).");
  if (LM.size() != cells)
    throw std::runtime_error("LM size mismatch: got " + std::to_string(LM.size()) +
                             ", expected " + std::to_string(cells) + " (one per cell).");
}

// Create unavailable mask from C
[[nodiscard]] std::vector<double> get_unavailable_mask(const std::vector<double>& C,
                                                       int n_h,
                                                       int dim_x,
                                                       int dim_y,
                                                       double sentinel = 1e10);

// Compute land spans by row and col for skipping sea
inline void compute_land_spans(const std::vector<uint8_t>& LM,
                               int dim_x,
                               int dim_y,
                               std::vector<int>& row_first_land,
                               std::vector<int>& row_last_land,
                               std::vector<int>& col_first_land,
                               std::vector<int>& col_last_land)
{
  // One entry per row / column
  row_first_land.assign(dim_y, dim_x);  // "no land" => first = dim_x
  row_last_land.assign(dim_y, -1);      // "no land" => last = -1

  col_first_land.assign(dim_x, dim_y);  // "no land" => first = dim_y
  col_last_land.assign(dim_x, -1);      // "no land" => last = -1

  // Scan LM row-major: g = r * dim_x + c
  for (int r = 0; r < dim_y; ++r) {
    std::size_t base = static_cast<std::size_t>(r) * static_cast<std::size_t>(dim_x);
    for (int c = 0; c < dim_x; ++c) {
      std::size_t g = base + static_cast<std::size_t>(c);
      if (!LM[g]) continue;  // sea

      if (c < row_first_land[r]) row_first_land[r] = c;
      if (c > row_last_land[r])  row_last_land[r]  = c;

      if (r < col_first_land[c]) col_first_land[c] = r;
      if (r > col_last_land[c])  col_last_land[c]  = r;
    }
  }
}

// Build Ebrel class object outputs
inline void create_ebrel_class_object(const std::vector<double>& E,
                                      const std::vector<double>& C,
                                      const std::vector<double>& SD,
                                      const std::vector<int>&    D,
                                      const std::vector<double>& SxH,
                                      const std::vector<double>& O,
                                      const std::vector<uint8_t>& LM,
                                      int dim_x,
                                      int dim_y,
                                      int n_h,
                                      int n_s,
                                      double sentinel,
                                      double sigma,
                                      std::vector<double>& U_out,
                                      std::vector<int>& row_first_land,
                                      std::vector<int>& row_last_land,
                                      std::vector<int>& col_first_land,
                                      std::vector<int>& col_last_land)
{
  // Check sizes, including LM
  validate_dimensions(E, C, SD, D, SxH, O, LM,
                      dim_x, dim_y, n_h, n_s);

  validate_existing_habitat(E, dim_x, dim_y, n_h);

  // Derive land spans from the LM mask
  compute_land_spans(LM, dim_x, dim_y,
                     row_first_land, row_last_land,
                     col_first_land, col_last_land);

  // Existing unavailable-mask logic
  U_out = get_unavailable_mask(C, n_h, dim_x, dim_y, sentinel);
}

#endif  // EBREL_SETUP_H
