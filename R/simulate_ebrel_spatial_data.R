#' Simulate spatial habitat and species distribution data for EBREL
#'
#' Generates virtual habitat configurations, land–sea masks, conversion costs,
#' and species distributions for testing spatial conservation and species
#' distribution models. The simulation supports scalable grid dimensions,
#' flexible numbers of habitats and species, grid-scaled species dispersal
#' distances drawn from a Gamma distribution, and a tunable proportion of
#' rare (low-prevalence) species.
#'
#' Habitat is generated as autocorrelated raster layers, optionally overlaid
#' with a simple land–sea mask. Species are assigned to one or more habitat
#' types, seeded into suitable cells, and their ranges are grown using a
#' breadth-first search constrained by habitat and species-specific dispersal
#' distances. Dispersal distances are drawn from a Gamma distribution whose
#' mean and variance are scaled to the grid dimension \code{dim_x} and
#' controlled by a discrete long-tail class (see \code{disp_longtail} and
#' \code{.map_longtail_to_gamma()}).
#'
#' @param dim_x Integer. Number of rows in the spatial grid (e.g. x/easting).
#' @param dim_y Integer. Number of columns in the spatial grid (e.g. y/northing).
#' @param n_h Integer. Number of distinct habitat types (the last type is treated
#'   as non-selectable/urban in the cost construction).
#' @param n_s Integer. Number of species.
#' @param disp_max Integer. Maximum allowable dispersal distance (in cells) for
#'   species. Gamma-drawn distances are truncated at this value.
#' @param disp_longtail Integer (0–3). Dispersal long-tail class passed to
#'   \code{.map_longtail_to_gamma()}. Controls the typical mean and spread of
#'   species' dispersal distances relative to grid size:
#'   \itemize{
#'     \item \code{0}: very short-tailed dispersal
#'     \item \code{1}: short-tailed dispersal
#'     \item \code{2}: intermediate dispersal
#'     \item \code{3}: long-tailed (more long-distance dispersers)
#'   }
#'   Larger values yield species with, on average, longer and more variable
#'   dispersal distances. (See \code{.map_longtail_to_gamma()} for exact mapping.)
#' @param rarity_bias Numeric (≥0). Controls how many species are rare
#'   (low-prevalence). A value of \code{0} produces more equal prevalence across
#'   species; higher values (e.g. 1–3) generate more strongly right-skewed
#'   prevalence distributions, with many rare species. Implemented via a
#'   Beta-distributed rarity weight.
#' @param fixed_O Optional numeric (>0). If supplied, applies the same occupancy
#'   target (proportional increase in occupancy) to all species. If \code{NULL},
#'   species-specific occupancy targets are drawn at random.
#' @param maskLandArea Logical. If \code{TRUE}, habitat and species layers are
#'   masked to the land area defined by the internally generated land mask.
#'   If \code{FALSE}, the land mask is still generated but returned as a separate
#'   layer and not used to mask habitat or species.
#' @param prop_sea Numeric in \eqn{(0, 1)}. Proportion of the grid height that
#'   can be occupied by sea when generating the land–sea mask. Passed to the
#'   internal coastline/land-mask generator.
#' @param convert_to_cpp_format Logical. If \code{TRUE} (default), returns data
#'   converted into the vector/matrix formats required by the EBREL C++ code via
#'   \code{prepare_ebrel_r_to_cpp()}. If \code{FALSE}, returns higher-level R
#'   objects (mostly \code{SpatRaster}s and matrices) suitable for inspection
#'   and further manipulation in R.
#' @param seed Optional integer. If provided, used to seed the random number
#'   generator to make results reproducible.
#'
#' @import terra
#' @import gstat
#'
#' @return
#' If \code{convert_to_cpp_format = FALSE}, returns a list:
#' \describe{
#'   \item{\code{dim_x}}{Number of rows.}
#'   \item{\code{dim_y}}{Number of columns.}
#'   \item{\code{n_h}}{Number of habitat types.}
#'   \item{\code{n_s}}{Number of species.}
#'   \item{\code{E}}{\code{SpatRaster} with \code{n_h} layers of binary habitat
#'     (0/1) after handling NAs/unavailable cells.}
#'   \item{\code{C}}{\code{SpatRaster} with \code{n_h} layers of habitat conversion
#'     costs (large values for unavailable/urban).}
#'   \item{\code{SD}}{\code{SpatRaster} with \code{n_s} layers of binary species
#'     presence/absence.}
#'   \item{\code{SxH}}{Matrix \code{[n_s, n_h]} indicating which species use which
#'     habitats (1 = used, 0 = not used).}
#'   \item{\code{D}}{Numeric vector of length \code{n_s} giving species dispersal
#'     distances (in cells) after truncation by \code{disp_max}.}
#'   \item{\code{O}}{Numeric vector of length \code{n_s} of species occupancy
#'     targets.}
#'   \item{\code{sigma}}{Scalar penalty weight based on the mean dispersal distance
#'     (\code{1 / mean(D)}).}
#'   \item{\code{LM}}{\code{SpatRaster} land mask (1 = land, 0 = non-land/sea).}
#' }
#'
#' If \code{convert_to_cpp_format = TRUE}, returns the list produced by
#' \code{prepare_ebrel_r_to_cpp()}, containing the same information restructured
#' into simple vectors/matrices suitable for passing directly to the EBREL C++
#' optimisation routines (see \code{prepare_ebrel_r_to_cpp()} for details).
#'
#' @examples
#' \dontrun{
#' # Basic example with 3 habitats and 15 species on a 50x50 grid
#' sim <- simulate_ebrel_spatial_data(
#'   dim_x = 50, dim_y = 50,
#'   n_h = 3, n_s = 15,
#'   disp_max = 20,
#'   rarity_bias = 2,
#'   convert_to_cpp_format = FALSE,
#'   seed = 123
#' )
#'
#' # Dispersal distances
#' hist(sim$D, main = "Species dispersal distances")
#'
#' # Plot first habitat layer
#' plot(sim$E[[1]], main = "Habitat 1")
#'
#' # Average prevalence per species
#' avg_prev <- mean(global(sim$SD, fun = "sum", na.rm = TRUE)) / (50 * 50)
#' avg_prev
#' }
#'
#' @export
simulate_ebrel_spatial_data <- function(
    dim_x,
    dim_y,
    n_h, # Number of selectable opportunity habitat types is n_h - 1
    n_s,
    disp_max,
    disp_longtail = 1,
    rarity_bias = 1,
    fixed_O = NULL,
    maskLandArea = FALSE,
    prop_sea = 0.2,
    convert_to_cpp_format = TRUE,
    seed = NULL) {
  # --- Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # --- Number of cells
  n_cells <- dim_x * dim_y

  # --- Max dispersal
  dim_max <- max(dim_x, dim_y)

  # --- Dispersal params
  gamma_params <- .map_longtail_to_gamma(disp_longtail, dim_max)
  long_disp_shape <- gamma_params$shape
  long_disp_scale <- gamma_params$scale

  # --- Generate an autocorrelated response variable
  hab_rast <- .generate_habitat_rast(
    dim_x = dim_x,
    dim_y = dim_y,
    n_h = n_h,
    unavail_hab_prop = 0.2,
    prop_sea = prop_sea
  )

  # --- Mask land or reset landMask
  landMask <- hab_rast$LM
  if (maskLandArea) {
    hab_rast[[2]] <- mask(hab_rast[[2]], landMask)
    hab_rast[[3]] <- mask(hab_rast[[3]], landMask)
  } else {
    # Blank out landMask if no masking occurred
    landMask <- setValues(landMask, NA)
  }

  # --- Generate rarity weights using beta distribution
  rarity_weights <- sort(rbeta(n_s, shape1 = 1, shape2 = 3 * rarity_bias + 0.1), decreasing = TRUE)

  # --- Simulated species probability distributions (pd)
  pd <- array(0, dim = c(dim_x, dim_y, n_s))

  # ---  Habitat association matrix
  SxH <- matrix(0, nrow = n_s, ncol = n_h)
  rownames(SxH) <- paste0("sp_", 1:n_s)
  colnames(SxH) <- paste0("hab_", 1:n_h)

  # --- Species dispersal distances ---
  D <- floor(rgamma(n_s, shape = long_disp_shape, scale = long_disp_scale))
  D <- pmax(1, pmin(D, disp_max))

  # --- Binary SD array
  SD_rast <- vector("list", n_s)

  # - Try to create species in landscape
  for (s in 1:n_s) {
    npres <- 0
    tries <- 0L

    # cap habitats drawn and align probs
    k_max <- min(4L, n_h)
    hab_probs <- c(0.5, 0.25, 0.125, 0.125)[seq_len(k_max)]

    repeat {
      tries <- tries + 1
      if (tries > 200L) {
        stop(
          "simulate_ebrel_spatial_data(): failed to place species ", s,
          " after 200 attempts (no eligible habitat cells)."
        )
      }

      # draw how many habitats and which ones
      k <- sample.int(k_max, 1, prob = hab_probs)
      hs <- sample.int((n_h), k, replace = FALSE)

      # cells with any of the chosen habitats
      hs_max <- terra::app(hab_rast$E[[hs]], "max")
      # sample seed cells ONLY from >0
      ssize <- sample.int(5L, 1)
      cells <- .safe_sample_cells(hs_max, ssize)
      if (length(cells) == 0L) next # try new hs

      # seed raster
      occ_seed_rast <- terra::setValues(hs_max, 0)
      occ_seed_rast[cells] <- 1

      # grow with BFS, constrained by hs_max
      SD_s <- .bfs_fill(hs_max, occ_seed_rast, width = D[s])

      # any presences?
      n_pres_cells <- as.numeric(terra::global(SD_s, "sum", na.rm = TRUE)[1, 1])
      if (!is.na(n_pres_cells) && n_pres_cells > 0) {
        SxH[s, hs] <- 1
        SD_rast[[s]] <- SD_s
        npres <- 1
        break
      }
      # else retry with different hs
    }
  }
  SD_rast <- do.call(c, SD_rast)

  # --- Targets
  if (is.null(fixed_O)) {
    O <- round(runif(n_s, min = 0.1, max = 0.3), 1)
  } else {
    O <- rep(fixed_O, n_s)
  }

  # --- Sigma for weights
  sigma <- 1 / mean(D)

  # --- Costs

  # - Set base costs for land management type per pixel
  hab_base_cost <- .approx_doubling(n_h, include_zero = FALSE, round_digits = 0)

  # - Set landscape ramp to scale base costs by geographic area
  x_ramp_vec <- seq(0.1, 1, length.out = dim_y)
  geo_C_ramp <- setValues(hab_rast$HO[[1]], rep(x_ramp_vec, dim_x))

  # - Assign costs
  C_rast <- lapply(seq_len(n_h-1), function(i) {
    C_i <- hab_rast$HO[[i]] # Get opportunity
    BC_i <- hab_base_cost[[i]] # Get base cost
    C_i <- ifel(C_i == 1, BC_i, 0)
    C_i <- C_i * geo_C_ramp
    C_i[C_i == 0] <- 1e10
    C_i
  })
  C_rast <- do.call(c, C_rast)
  # Add sentinel cost for urban (i.e. unavailable everywhere)
  C_urb <- setValues(C_rast[[1]], 1e10)
  C_rast <- c(C_rast, C_urb)
  # Name
  names(C_rast) <- paste0("C", 1:(n_h))

  # --- Sort out NAs
  E_rast <- ifel(is.na(hab_rast$E), 0, hab_rast$E)
  SD_rast <- ifel(is.na(SD_rast), 0, SD_rast)
  C_rast <- ifel(is.na(C_rast), 1e10, C_rast)
  LM_rast <- ifel(is.na(landMask), 0, landMask)

  # --- Simulated data in R formats
  sim_r <- list(
    dim_x = dim_x,
    dim_y = dim_y,
    n_h = n_h,
    n_s = n_s,
    E = E_rast,
    C = C_rast,
    SD = SD_rast,
    SxH = SxH,
    D = D,
    O = O,
    sigma = sigma,
    LM = LM_rast
  )

  # --- Return
  if (convert_to_cpp_format) {
    # Convert R data formats to cpp vector formats required for ebrel c++
    sim_cpp <-
      prepare_ebrel_r_to_cpp(
        E_rast = E_rast,
        C_rast = C_rast,
        SD_rast = SD_rast,
        D_vec = D,
        SxH_mat = SxH,
        O_vec = O,
        sigma = sigma,
        LM_rast = LM_rast
      )
    return(sim_cpp)
  } else {
    # Return in R formats, but will require running 'prepare_ebrel_r_to_cpp()' later
    sim_r <- list(
      dim_x = dim_x,
      dim_y = dim_y,
      n_h = n_h,
      n_s = n_s,
      E = E_rast,
      C = C_rast,
      SD = SD_rast,
      SxH = SxH,
      D = D,
      O = O,
      sigma = sigma,
      LM = LM_rast
    )
    return(sim_r)
  }
}

#' Generate autocorrelated habitat rasters with availability mask
#'
#' Creates a set of \code{SpatRaster} layers representing (i) continuous habitat
#' suitability surfaces (\code{HP}), (ii) a binary one-hot existing habitat stack
#' (\code{E}), and (iii) a binary habitat opportunity stack (\code{HO}) after
#' masking out an "urban/unavailable" fraction of the landscape.
#'
#' Internally, unconditional Gaussian fields are simulated with \pkg{gstat},
#' rescaled to \eqn{[0,1]}, and the lowest \code{unavail_hab_prop} quantile is
#' treated as unavailable. \code{E} is built by one-hot encoding the
#' layer-wise \code{which.max(HP)} and multiplying by a probabilistic
#' conversion mask; the function retries that conversion until every habitat
#' layer in \code{E} has at least one positive cell.
#'
#' @param dim_x Integer. Number of rows in the grid. Default \code{50}.
#' @param dim_y Integer. Number of columns in the grid. Default \code{50}.
#' @param n_h Integer. Number of habitat layers to simulate. Default \code{4}.
#' @param unavail_hab_prop Numeric in \eqn{[0,1)}. Proportion of cells treated as
#'   unavailable (e.g., urban) and masked out of \code{HP}/\code{E}/\code{HO}.
#'   Default \code{0.25}.
#'
#' @return A list with three \code{SpatRaster} stacks (each with \code{n_h} layers):
#' \describe{
#'   \item{\code{HP}}{Continuous habitat potentials, names \code{hp1..hpN}.}
#'   \item{\code{E}}{Binary existing habitat (one-hot), names \code{E1..EN}.}
#'   \item{\code{HO}}{Binary habitat opportunities after masking, names \code{HO1..HON}.}
#' }
#'
#' @details
#' This function is stochastic. For reproducibility, set R’s seed (e.g., \code{set.seed(1)})
#' before calling. The internal retry loop ensures \emph{every} \code{E} layer has at least
#' one positive cell.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- generate_habitats_rast(dim_x = 30, dim_y = 40, n_h = 3, unavail_hab_prop = 0.2)
#' sim$HP; sim$E; sim$HO
#' terra::global(sim$E, "sum", na.rm = TRUE)  # cells per habitat
#' }
#'
#' @importFrom gstat vgm gstat
#' @importFrom terra rast values app setValues which.max ifel global minmax ncell names
#' @noRd
.generate_habitat_rast <- function(dim_x = 50, dim_y = 50, n_h = 4,
                                  unavail_hab_prop = 0.1, prop_sea = 0.2) {

  # --- Create autocorrelated surface
  xy <- expand.grid(x = 1:dim_x, y = 1:dim_y)
  varioMod <- vgm(psill = 1, range = 40, model = 'Exp', nugget = 0.0)
  zDummy <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE,
                  beta = 0, model = varioMod, nmax = 50)

  # -- Simulate surface for "urban"
  urban_sim <- predict(zDummy, newdata = xy, nsim = 1, debug.level = 0)
  urban <- .rescale01(rast(urban_sim, type = "xyz"))
  urban_q <- quantile(values(urban), prob = unavail_hab_prop)
  urban_binary <- ifel(urban < urban_q, 1, NA)

  # -- Simulate surfaces for "habitat"
  habitat_sim <- predict(zDummy, newdata = xy, nsim = (n_h-1), debug.level = 0)
  habitat <- rast(habitat_sim, type = "xyz")
  for (x in seq_len(n_h-1)) {
    habitat[[x]] <- .rescale01(habitat[[x]])
  }
  habitat <- ifel(urban_binary == 1, NA, habitat)

  # -- Establish habitat opportunities (HO)
  HO <- habitat
  for (h in 1:(n_h-1)) {
    qX <- quantile(values(HO[[h]]), prob = pmin(rnorm(1, 0.90,0.01), 0.99), na.rm = TRUE
    ) # areas mod to small
    HO[[h]] <- ifel(HO[[h]] < qX, 0, 1) # also mask urban
  }

  # - Add zeroed H0 for urban (i.e. no opportunities)
  HO <- c(HO, setValues(HO[[1]],0))

  # -- Establish existing habitat
  # - One-hot encode E (i.e. only one habitat per tile)
  habitat_max <- which.max(habitat)
  E_hab <- rast(lapply(seq_len(n_h-1), function(k) 1 * (habitat_max == k)))
  E_hab <- mask(E_hab, urban_binary, inverse = TRUE)
  E_hab[is.na(E_hab)] <- 0

  # - Convert HO to existing habitat by converting land
  layers_zero <- 0
  while (layers_zero == 0) {

    # Simulate converted land (e.g. agriculture)
    convert_sim <- predict(zDummy, newdata = xy, nsim = 1, debug.level = 0)
    convert <- (.rescale01(rast(convert_sim, type = "xyz"))^3) * 0.5

    # Convert to binary (presence/absence) and mask out of E
    convert_binary <- .prob_to_bin(convert)
    E <- E_hab * convert_binary

    # Check that all habitats exist
    # per-layer max (vector, one per layer)
    mx <- terra::global(E, "max", na.rm = TRUE)[, 1]
    zero_max_idx <- which(is.na(mx) | mx == 0)
    if (length(zero_max_idx) == 0)
      layers_zero <- 1
  }

  # - Add back in urban (no opportunity but can be existing)
  E <- c(E, urban_binary)
  E[is.na(E)] <- 0

  # - Set names
  names(habitat) <- paste0("hp", 1:(n_h-1))
  names(E) <- paste0("E", 1:n_h)
  names(HO) <- paste0("HO", 1:n_h)

  # - Make sea mask
  LM <- setValues(habitat[[1]], 0)
  LM <- .make_land_mask(LM, prop_sea = prop_sea, smooth_window = 5)

  list(
    HP = habitat,
    E = E,
    HO = HO,
    LM = LM
  )
}

#' Generate a simple land–sea mask with a smoothed coastline
#'
#' Creates a single-layer \code{SpatRaster} mask distinguishing land from sea,
#' based on a template raster. The coastline is generated by drawing a random
#' shoreline height (in rows from the bottom of the raster) per column,
#' smoothing it with a moving-average filter, and then filling cells below the
#' coastline as sea.
#'
#' The mask is aligned to the geometry (extent, resolution, CRS) of the input
#' raster \code{r}. Cell values are set to \code{1} on land and \code{NA} over
#' sea.
#'
#' @param r A \code{terra::SpatRaster} object providing the target geometry
#'   (extent, resolution, projection, and dimensions) for the mask.
#' @param prop_sea Numeric in \eqn{(0, 1]}. Approximate maximum proportion of
#'   the raster height (from the bottom) that can be occupied by sea. This
#'   determines the upper bound on the randomly drawn coastline heights.
#' @param smooth_window Integer \eqn{\ge 1}. Width of the moving-average filter
#'   used to smooth the randomly drawn coastline heights across columns. Larger
#'   values produce smoother coastlines.
#' @param seed Optional integer. If provided, passed to \code{set.seed()} to
#'   make the coastline generation reproducible.
#'
#' @return A single-layer \code{terra::SpatRaster} with the same geometry as
#'   \code{r}, named \code{"land"}, where:
#'   \itemize{
#'     \item \code{1} indicates land.
#'     \item \code{NA} indicates sea.
#'   }
#'
#' @details
#' Let \eqn{n_r} and \eqn{n_c} be the number of rows and columns in \code{r}.
#' The function:
#' \enumerate{
#'   \item Computes \code{max_sea_rows = round(nr * prop_sea)}, the maximum
#'     number of rows from the bottom that may be sea.
#'   \item Draws a random coastline height (in rows from the bottom) for each
#'     column, uniformly between \code{1} and \code{max_sea_rows}.
#'   \item Smooths these heights across columns using a moving-average filter
#'     of width \code{smooth_window}.
#'   \item For each column, marks all cells below the smoothed coastline as
#'     sea (\code{NA}) and the remaining cells as land (\code{1}).
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Example template raster (e.g. 100 x 100 grid)
#' r <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#'
#' # Generate a land–sea mask with ~30% sea at the bottom and a smooth coastline
#' land_mask <- make_land_mask(r, prop_sea = 0.3, smooth_window = 7, seed = 42)
#'
#' plot(land_mask)
#' }
#'
#' @keywords internal
#' @noRd
.make_land_mask <- function(r, prop_sea = 0.3, smooth_window = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  nr <- nrow(r)
  nc <- ncol(r)

  # Target maximum sea height (in rows from the bottom)
  max_sea_rows <- max(1, round(nr * prop_sea))

  # Random coastline height per column (from 1 to max_sea_rows)
  coast_raw <- sample(1:max_sea_rows, nc, replace = TRUE)

  # Smooth coastline
  k <- rep(1 / smooth_window, smooth_window)
  coast_smooth <- stats::filter(coast_raw, k, sides = 2)

  # replace NAs at edges with original values
  coast_smooth[is.na(coast_smooth)] <- coast_raw[is.na(coast_smooth)]
  coast <- pmax(1, pmin(max_sea_rows, round(coast_smooth)))

  # Build a matrix where 1 = land, NA = sea
  mat <- matrix(1, nrow = nr, ncol = nc)
  for (j in seq_len(nc)) {
    sea_rows <- (nr - coast[j] + 1):nr   # rows that should be NA for column j
    mat[sea_rows, j] <- NA               # row first, then column
  }

  # Create SpatRaster with same geometry as input
  land_r <- r
  values(land_r) <- as.vector(t(mat))
  names(land_r) <- "land"

  return(land_r)
}


#' Buffer binary cells on a raster
#'
#' Internal function that applies a square buffer of radius \code{width} and
#' returns a spatraster where any neighbor within the window is converted to 1.
#'
#' @param r A \code{SpatRaster} (single layer) with binary values 0/1.
#' @param width Integer \eqn{\ge 0}. Neighborhood radius; window size is
#'   \code{(2*width + 1) x (2*width + 1)}. Default \code{2}.
#'
#' @return A \code{SpatRaster} (single layer) of 0/1 after function
#'
#' @keywords internal
#' @noRd
#' @importFrom terra focal ifel
.buffer_cells <- function(r, width = 2) {
  w <- matrix(1, nrow = 2*width + 1, ncol = 2*width + 1)
  out <- focal(r, w = w, fun = max, na.policy = "omit", pad = TRUE, padValue = 0)
  out <- ifel(out > 0, 1, 0)
}

#' BFS-style region growth constrained by suitability and hop width
#'
#' Internal function, that, tarting from seed cells (\code{occ_seed_rast > 0}),
#' repeatedly grows the visited set using \code{buffer_cells()} with radius
#' \code{width}, but only admits cells where \code{hs_max > 0}. Stops when no
#' new cells are added.
#'
#' @param hs_max \code{SpatRaster} (single layer) with 0/1 suitability mask.
#' @param occ_seed_rast \code{SpatRaster} (single layer) with 0/1 initial seeds.
#' @param width Integer \eqn{\ge 0}. Maximum hop size per iteration. Default \code{2}.
#'
#' @return A \code{SpatRaster} (single layer) of 0/1 indicating the filled region.
#'
#' @keywords internal
#' @noRd
#' @importFrom terra ifel global
.bfs_fill <- function(hs_max, occ_seed_rast, width = 2) {

  # ensure binary 0/1
  hs <- ifel(hs_max > 0, 1, 0)
  visited <- ifel(occ_seed_rast > 0, 1, 0)

  repeat {
    expanded <- .buffer_cells(visited, width = width)

    # candidate cells reachable in one hop AND suitable
    candidates <- ifel((expanded == 1) & (hs == 1), 1, 0)
    new <- ifel((candidates == 1) & (visited == 0), 1, 0)

    n_new <- as.numeric(global(new, "sum", na.rm = TRUE)[1, 1])
    if (is.na(n_new) || n_new == 0) break

    # add new frontier to visited
    visited <- ifel(new == 1, 1, visited)
  }

  # final 0/1 raster of filled area (restricted to suitable habitat)
  ifel(visited == 1, 1, 0)
  visited
}

#' Rescale to [0,1]
#'
#' Linearly rescales a \code{SpatRaster} to the unit interval using its layer-wise
#' minimum and maximum; safeguards against zero range. Internal function
#'
#' @param x A \code{SpatRaster} (single layer).
#' @return A \code{SpatRaster} on \eqn{[0,1]}.
#'
#' @keywords internal
#' @noRd
#' @importFrom terra minmax
.rescale01 <- function(x) {
  rng <- terra::minmax(x); (x - rng[1]) / max(1e-9, (rng[2] - rng[1]))
}

#' Bernoulli thresholding of a probability raster
#'
#' Converts a probability raster \code{p} into a binary raster by drawing
#' \code{U(0,1)} per cell and returning \code{1*(U < p)}. Internal function.
#'
#' @param prob_rast A \code{SpatRaster} (single layer) with values in \eqn{[0,1]}.
#' @param seed Optional integer seed; if provided, sets R’s RNG state before sampling.
#'
#' @return A \code{SpatRaster} (single layer) of 0/1 with name \code{"binary"}.
#'
#' @keywords internal
#' @noRd
#' @importFrom terra rast values ncell
.prob_to_bin <- function(prob_rast, seed = NULL) {

  # draw uniform(0,1) for each cell
  rnd <- rast(prob_rast)
  values(rnd) <- runif(ncell(rnd))

  out <- (rnd < prob_rast) * 1
  names(out) <- "binary"
  out
}

#' Generate an approximately doubling cost scale
#'
#' Creates a short geometric sequence that ends at \code{max_val} and (optionally)
#' starts with 0. Useful for setting tiered costs where each successive tier is
#' about twice the previous one.
#'
#' @param n Integer (>= 1). Number of values to return.
#' @param max_val Numeric (> 0). The maximum (final) positive value in the sequence.
#' @param include_zero Logical. If \code{TRUE}, the first element is 0 and the
#'   remaining \code{n-1} values form a halving series up to \code{max_val}.
#'   If \code{FALSE}, all \code{n} values are positive.
#' @param round_digits Integer or \code{NULL}. If not \code{NULL}, round the
#'   resulting values to this many decimal places.
#'
#' @return A numeric vector of length \code{n}.
#' \itemize{
#'   \item If \code{include_zero = TRUE} and \code{n > 1}:
#'         \code{c(0, max_val/2^{m-1}, ..., max_val/2, max_val)} with \code{m = n-1}.
#'   \item If \code{include_zero = TRUE} and \code{n = 1}: \code{0}.
#'   \item If \code{include_zero = FALSE}:
#'         \code{max_val/2^{n-1}}, \code{...}, \code{max_val/2}, \code{max_val}.
#' }
#'
#' @examples
#' approx_doubling(4, max_val = 800, include_zero = TRUE)
#' # 0 200 400 800
#'
#' approx_doubling(4, max_val = 800, include_zero = FALSE)
#' # 100 200 400 800
#'
#' approx_doubling(5, max_val = 1000, include_zero = TRUE, round_digits = 0)
#'
#' @keywords internal
#' @noRd
.approx_doubling <- function(n, max_val = 1000, include_zero = TRUE, round_digits = NULL) {
  stopifnot(n >= 1)

  if (include_zero) {
    if (n == 1) return(0)
    m <- n - 1                 # number of positive values
    vals <- max_val / 2^((m - 1):0)  # geometric, ends at max_val
    out <- c(0, vals)
  } else {
    m <- n
    vals <- max_val / 2^((m - 1):0)
    out <- vals
  }

  if (!is.null(round_digits)) out <- round(out, round_digits)
  out
}

#' Map a long-tail dispersal class to Gamma parameters (grid-scaled)
#'
#' Converts a discrete long-tail class into parameters of a
#' \eqn{\mathrm{Gamma}(\text{shape}, \text{scale})} distribution used to draw
#' dispersal distances. The mapping is defined in terms of the maximum possible
#' distance on the grid (\code{max_dist}), allowing dispersal distances to scale
#' automatically with grid size (e.g., 25×25 vs 200×200).
#'
#' Each \code{longtail} class corresponds to a pair of fractions
#' (\code{mean_frac}, \code{sd_frac}) that determine the mean and standard
#' deviation of the Gamma distribution as:
#' \deqn{\mu = \text{mean\_frac} \times \text{max\_dist}}
#' \deqn{\sigma = \text{sd\_frac} \times \text{max\_dist}}
#'
#' These are then converted to Gamma parameters via:
#' \deqn{\text{shape} = (\mu^2) / \sigma^2}
#' \deqn{\text{scale} = \sigma^2 / \mu}
#'
#' Larger \code{longtail} values correspond to increasingly long-tailed
#' dispersal distributions (larger mean fraction, larger SD fraction).
#'
#' @param longtail Integer (0–3). Dispersal class where:
#'   \itemize{
#'     \item \code{0}: very short-tailed (\code{mean_frac = 0.05}, \code{sd_frac = 0.10})
#'     \item \code{1}: short-tailed    (\code{mean_frac = 0.15}, \code{sd_frac = 0.15})
#'     \item \code{2}: intermediate    (\code{mean_frac = 0.25}, \code{sd_frac = 0.25})
#'     \item \code{3}: long-tailed     (\code{mean_frac = 0.35}, \code{sd_frac = 0.35})
#'   }
#'
#' @param max_dist Numeric. The maximum possible cell-to-cell distance on the grid
#'   (typically the grid diagonal, e.g., \code{sqrt(nx^2 + ny^2)}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{shape}: Gamma shape parameter.
#'   \item \code{scale}: Gamma scale parameter.
#' }
#'
#' @examples
#' # Suppose a 100x100 grid:
#' max_dist <- sqrt(100^2 + 100^2)
#'
#' .map_longtail_to_gamma(0, max_dist)
#' .map_longtail_to_gamma(1, max_dist)
#' .map_longtail_to_gamma(3, max_dist)
#'
#' @keywords internal
#' @noRd
.map_longtail_to_gamma <- function(longtail, max_dist) {
  # choose relative behaviour by 'longtail'
  base <- switch(as.character(longtail),
                 "0" = list(mean_frac = 0.05, sd_frac = 0.1),
                 "1" = list(mean_frac = 0.15, sd_frac = 0.15),
                 "2" = list(mean_frac = 0.25, sd_frac = 0.25),
                 "3" = list(mean_frac = 0.35, sd_frac = 0.35),
                 stop("Unknown longtail value")
  )

  mean_D <- base$mean_frac * max_dist
  sd_D   <- base$sd_frac   * max_dist

  shape <- (mean_D^2) / (sd_D^2)
  scale <- (sd_D^2) / mean_D

  list(shape = shape, scale = scale)
}


#' Sample eligible cell indices from a mask raster
#'
#' Returns uniformly sampled cell indices from the set of cells where the mask
#' has finite values greater than 0. If there are fewer eligible cells than
#' requested, the sample size is reduced accordingly. If no cells are eligible,
#' an empty integer vector is returned.
#'
#' @param mask_rast A `terra::SpatRaster` mask; cells with `values(mask_rast) > 0`
#'   (and finite) are considered eligible.
#' @param size Integer (>= 1). Desired number of cells to sample.
#'
#' @return An integer vector of 1-based cell indices (length `<= size`).
#'
#' @examples
#' \dontrun{
#' r <- terra::rast(nrows = 5, ncols = 5)
#' terra::values(r) <- c(0,1,2,NA,0)[(1:25 - 1) %% 5 + 1]
#' safe_sample_cells(r, size = 3)  # indices where values > 0
#' }
#'
#' @importFrom terra values
#' @keywords internal
#' @noRd
.safe_sample_cells <- function(mask_rast, size) {
  v <- terra::values(mask_rast)
  idx <- which(is.finite(v) & v > 0)
  if (length(idx) == 0L) return(integer(0))
  size <- min(size, length(idx))
  sample(idx, size)
}
