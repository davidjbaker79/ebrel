#' Create a C++-friendly object for EBREL optimization
#'
#' Prepares raster and matrix data as flattened vectors suitable for passing to C++.
#' This function extracts and flattens relevant input rasters and matrices into a list
#' of vectors used by EBREL's C++ backend.
#'
#' @details
#' Let \code{n_cells <- dim_x * dim_y}, \code{n_h <- nlyr(E_rast)}, and \code{n_s <- nlyr(SD_rast)}.
#' Flattening is **layer-major** (column-major on the \code{values(...)} matrix): all cells of layer 1,
#' then all cells of layer 2, etc. This yields vectors of lengths \code{n_cells * n_h} for habitat/cost
#' rasters and \code{n_cells * n_s} for species rasters.
#'
#' @param E_rast A `SpatRaster` with \strong{n_h} layers: habitat suitability/categorical habitat per habitat option.
#' @param C_rast A `SpatRaster` with \strong{n_h} layers: conversion costs/constraints per habitat option.
#'   Values \eqn{\ge} 1e11 indicate unavailable cells (kept as-is; interpreted downstream).
#' @param SD_rast A `SpatRaster` with \strong{n_s} layers: species/feature baseline distributions or constraints.
#' @param D_vec Numeric vector of length \strong{n_s}: maximum dispersal (in grid cells) for each species/feature.
#' @param SxH_mat Numeric/integer/logical matrix of shape \strong{n_s × n_h}: species–habitat associations (ideally 0/1).
#' @param O_vec Numeric vector of length \strong{n_s}: species/feature targets (e.g., proportional increases).
#'
#' @return A list with flattened numeric or integer vectors:
#' \describe{
#'   \item{E}{length \code{n_cells * n_h}}
#'   \item{C}{length \code{n_cells * n_h}}
#'   \item{SD}{length \code{n_cells * n_s}}
#'   \item{D}{length \code{n_s}}
#'   \item{SxH}{length \code{n_s * n_h}}
#'   \item{O}{length \code{n_s}}
#' }
#'
#' @examples
#' \dontrun{
#' obj <- prepare_ebrel_r_to_cpp(E_rast, C_rast, SD_rast, D_vec, SxH_mat, O_vec)
#' }
#' @export
prepare_ebrel_r_to_cpp <- function(E_rast, C_rast, SD_rast, D_vec, SxH_mat, O_vec) {
  # ---- Type checks ----
  if (!inherits(E_rast, "SpatRaster")) stop("E_rast must be a terra::SpatRaster.")
  if (!inherits(C_rast, "SpatRaster")) stop("C_rast must be a terra::SpatRaster.")
  if (!inherits(SD_rast, "SpatRaster")) stop("SD_rast must be a terra::SpatRaster.")

  # ---- Layer counts ----
  n_cells <- terra::ncell(E_rast)
  n_h <- terra::nlyr(E_rast)
  n_hC <- terra::nlyr(C_rast)
  n_s <- terra::nlyr(SD_rast)

  if (n_h < 1L) stop("E_rast must have at least one layer (n_h >= 1).")
  if (n_s < 1L) stop("SD_rast must have at least one layer (n_s >= 1).")
  if (n_hC != n_h) {
    stop(sprintf("C_rast must have the same number of layers as E_rast (n_h). Found: C_rast=%d, E_rast=%d.", n_hC, n_h))
  }

  # ---- Vector/matrix argument checks ----
  if (!is.numeric(D_vec) || is.matrix(D_vec) || is.list(D_vec)) {
    stop("D_vec must be a numeric vector (length n_s).")
  }
  if (length(D_vec) != n_s) {
    stop(sprintf("D_vec length (%d) must equal the number of species layers in SD_rast (n_s = %d).", length(D_vec), n_s))
  }
  if (!is.matrix(SxH_mat)) stop("SxH_mat must be a matrix with shape n_s x n_h.")
  if (nrow(SxH_mat) != n_s || ncol(SxH_mat) != n_h) {
    stop(sprintf("SxH_mat must have dimensions n_s x n_h (%d x %d); found %d x %d.",
                 n_s, n_h, nrow(SxH_mat), ncol(SxH_mat)))
  }
  if (!(is.numeric(SxH_mat) || is.integer(SxH_mat) || is.logical(SxH_mat))) {
    stop("SxH_mat must be numeric, integer, or logical.")
  }
  # Warn (not stop) if not strictly binary
  unique_sxh <- unique(as.vector(SxH_mat))
  if (!all(unique_sxh %in% c(0, 1, NA))) {
    warning("SxH_mat contains values other than 0/1/NA; values will be passed through as-is.")
  }

  if (!is.numeric(O_vec) || is.matrix(O_vec) || is.list(O_vec)) {
    stop("O_vec must be a numeric vector (length n_s).")
  }
  if (length(O_vec) != n_s) {
    stop(sprintf("O_vec length (%d) must equal the number of species layers in SD_rast (n_s = %d).", length(O_vec), n_s))
  }

  # ---- Flatten rasters (layer-major: all cells of layer 1, then layer 2, ...) ----
  # Using mat=TRUE ensures an ncell x nlayer matrix; as.vector() (column-major) gives layer-major flattening.
  E_mat  <- terra::values(E_rast, mat = TRUE)
  C_mat  <- terra::values(C_rast, mat = TRUE)
  SD_mat <- terra::values(SD_rast, mat = TRUE)

  if (!is.matrix(E_mat) || nrow(E_mat) != n_cells || ncol(E_mat) != n_h) {
    stop("Internal error extracting E_rast values: unexpected matrix shape.")
  }
  if (!is.matrix(C_mat) || nrow(C_mat) != n_cells || ncol(C_mat) != n_h) {
    stop("Internal error extracting C_rast values: unexpected matrix shape.")
  }
  if (!is.matrix(SD_mat) || nrow(SD_mat) != n_cells || ncol(SD_mat) != n_s) {
    stop("Internal error extracting SD_rast values: unexpected matrix shape.")
  }

  E  <- as.vector(E_mat)
  C  <- as.vector(C_mat)
  SD <- as.vector(SD_mat)

  # ---- Final sanity checks on lengths ----
  if (length(E)  != n_cells * n_h) stop("E length mismatch after flattening.")
  if (length(C)  != n_cells * n_h) stop("C length mismatch after flattening.")
  if (length(SD) != n_cells * n_s) stop("SD length mismatch after flattening.")

  # Optional heads-up about unavailability sentinel in costs
  if (any(is.finite(C) & C >= 1e11)) {
    warning("C contains values >= 1e11; these will be treated as unavailable by downstream code.")
  }

  # ---- Return list formatted for C++ use ----
  list(
    E   = E,
    C   = C,
    SD  = SD,
    D   = as.numeric(D_vec),
    SxH = as.vector(as.numeric(SxH_mat)),
    O   = as.numeric(O_vec)
  )
}
