#' Run EBREL optimisation from terra rasters
#'
#' High-level user-facing wrapper:
#' - extracts raster values in R (terra-friendly)
#' - passes flat vectors to C++
#' - C++ builds caches (U, spans, caches, sigma, W) and runs optimisation
#' - returns results as rasters + diagnostics
#'
#' @param E_rast SpatRaster (1 layer) with habitat id per cell: -1..(n_h-1)
#'        (If you still have multi-layer E, you’ll need a conversion step before this.)
#' @param C_rast SpatRaster with n_h layers (cost per habitat option)
#' @param SD_rast SpatRaster with n_s layers
#' @param D_vec numeric/integer vector length n_s
#' @param SxH_mat matrix n_s x n_h
#' @param O_vec numeric vector length n_s
#' @param LM_rast SpatRaster (1 layer), land mask 0/1
#' @param opt list of SA/run options (passed through to C++)
#' @param sentinel numeric, cost sentinel encoding unavailability (default 1e10)
#' @param return_maps logical, if TRUE return rasters for best solution etc.
#' @export
run_ebrel_R <- function(E_rast,
                        C_rast,
                        SD_rast,
                        D_vec,
                        SxH_mat,
                        O_vec,
                        LM_rast,
                        opt = list(),
                        sentinel = 1e10) {

  # ---- Check alignment of rasters ----
  stopifnot(terra::compareGeom(C_rast, LM_rast, stopOnError = FALSE))
  stopifnot(terra::compareGeom(C_rast, E_rast,  stopOnError = FALSE))
  stopifnot(terra::compareGeom(C_rast, SD_rast, stopOnError = FALSE))

  # ---- Checks ----
  stopifnot(inherits(E_rast, "SpatRaster"),
            inherits(C_rast, "SpatRaster"),
            inherits(SD_rast, "SpatRaster"),
            inherits(LM_rast, "SpatRaster"))

  dim_x <- terra::ncol(C_rast)
  dim_y <- terra::nrow(C_rast)
  n_cells <- terra::ncell(C_rast)
  n_h <- terra::nlyr(C_rast)
  n_s <- terra::nlyr(SD_rast)

  # ---- Checks ----
  if (terra::nlyr(E_rast) != 1)
    stop("E_rast must be a single-layer raster with values -1..(n_h-1)")
  if (length(D_vec) != n_s)
    stop("D_vec length must equal n_s")
  if (!is.matrix(SxH_mat) || nrow(SxH_mat) != n_s || ncol(SxH_mat) != n_h)
    stop("SxH_mat must be n_s x n_h")
  if (length(O_vec) != n_s)
    stop("O_vec length must equal n_s")

  # ---- Extract & flatten ---
  E  <- as.integer(terra::values(E_rast, mat = FALSE))
  C  <- as.vector(terra::values(C_rast, mat = TRUE))
  SD <- as.vector(terra::values(SD_rast, mat = TRUE))
  LM <- as.integer(terra::values(LM_rast, mat = FALSE))

  if (any(is.na(E))) stop("E contains NA")
  if (any(E < -1L | E >= n_h))
    stop("E values must be in -1..(n_h-1)")
  if (!all(LM %in% c(0L, 1L)))
    stop("LM must be 0/1")

  # ---- Call C++ one-shot runner ----
  res <- run_ebrel_cpp(
    E_int  = E,
    C  = C,
    SD = SD,
    D  = as.integer(D_vec),
    SxH = as.numeric(t(SxH_mat)),
    O   = as.numeric(O_vec),
    LM_int  = LM,
    dim_x = dim_x,
    dim_y = dim_y,
    n_h   = n_h,
    n_s   = n_s,
    sentinel = sentinel,
    opt = opt
  )

  # --- reconstruct raster ---
  X_best_rast <- terra::rast(E_rast)
  terra::values(X_best_rast) <- res$X_best

  res$X_best_rast <- X_best_rast
  res
}
