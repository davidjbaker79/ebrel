//------------------------------ Ebrel builder ---------------------------------

#include "ebrel_builder.h"

#include "ebrel_setup.h"
#include "species_plan.h"
#include "optimisation_utils.h"  // compute_distance_weights (or wherever it lives)

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>   // std::move

// ------------------------- Local helper functions ----------------------------

namespace {

inline void ensure(bool ok, const char* msg) {
  if (!ok) throw std::invalid_argument(msg);
}

} // namespace

//---------------------------- Main function -----------------------------------

BuiltEBREL build_ebrel(
    std::vector<int8_t>  E,
    std::vector<double>  C,
    std::vector<double>  SD,
    std::vector<int>     D,
    std::vector<double>  SxH,
    std::vector<double>  O,
    std::vector<uint8_t> LM,
    int dim_x,
    int dim_y,
    int n_h,
    int n_s,
    double sentinel,
    double sigma_in,
    int universal_disp_thres,
    int max_disp_steps,
    int roi_cap,
    bool precompute_W
) {
  BuiltEBREL b;

  // ---- basic sanity ----
  ensure(dim_x > 0 && dim_y > 0, "dim_x and dim_y must be positive");
  ensure(n_h   > 0 && n_s   > 0, "n_h and n_s must be positive");

  const std::size_t cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  const std::size_t Hsz   = cells * static_cast<std::size_t>(n_h);
  const std::size_t Ssz   = cells * static_cast<std::size_t>(n_s);
  const std::size_t SxHsz = static_cast<std::size_t>(n_h) * static_cast<std::size_t>(n_s);

  ensure(E.size()  == cells, "E size mismatch (expected dim_x*dim_y)");
  ensure(LM.size() == cells, "LM size mismatch (expected dim_x*dim_y)");
  ensure(C.size()  == Hsz,   "C size mismatch (expected dim_x*dim_y*n_h)");
  ensure(SD.size() == Ssz,   "SD size mismatch (expected dim_x*dim_y*n_s)");
  ensure(SxH.size()== SxHsz, "SxH size mismatch (expected n_h*n_s)");
  ensure(O.size()  == static_cast<std::size_t>(n_s), "O size mismatch (expected n_s)");
  ensure(D.size()  == static_cast<std::size_t>(n_s), "D size mismatch (expected n_s)");

  // range checks on E
  for (const int8_t v : E) {
    ensure(v >= -1 && v < n_h, "E contains value outside [-1, n_h-1]");
  }

  // ---- fill core fields ----
  b.in.dim_x = dim_x;
  b.in.dim_y = dim_y;
  b.in.n_h   = n_h;
  b.in.n_s   = n_s;

  b.in.universal_disp_thres = universal_disp_thres;
  b.in.max_disp_steps       = max_disp_steps;
  b.in.roi_cap              = roi_cap;

  b.in.E  = std::move(E);
  b.in.C  = std::move(C);
  b.in.SD = std::move(SD);
  b.in.D  = std::move(D);
  b.in.SxH= std::move(SxH);
  b.in.O  = std::move(O);
  b.in.LM = std::move(LM);

  // ---- setup derived state (U, spans, coords, tiles, rowruns, sigma) ----
  create_ebrel_class_object(
    b.in.E, b.in.C, b.in.SD, b.in.D, b.in.SxH, b.in.O, b.in.LM,
    b.in.universal_disp_thres,
    dim_x, dim_y, n_h, n_s,
    sentinel,
    sigma_in,
    b.sigma,
    b.in.U,
    b.in.row_first_land, b.in.row_last_land,
    b.in.col_first_land, b.in.col_last_land,
    b.in.cell_r, b.in.cell_c,
    b.in.Etiles_per_h,
    b.in.rowruns_cache
  );

  // ---- species plan ----
  b.in.species_plan = build_species_plan(
    b.in.SD, b.in.SxH, b.in.D, b.in.O,
    n_h, n_s, dim_x, dim_y,
    b.in.universal_disp_thres,
    b.in.max_disp_steps,
    b.in.roi_cap,
    b.in.row_first_land, b.in.row_last_land,
    b.in.col_first_land, b.in.col_last_land,
    b.in.E,
    b.in.cell_r, b.in.cell_c
  );

  // ---- proposal weights W (distance decay) ----
  if (precompute_W) {
    b.W = compute_distance_weights(b.in.E, b.in.U, n_h, dim_x, dim_y, b.sigma);
  } else {
    b.W.clear();
  }

  return b;
}
