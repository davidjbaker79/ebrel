#' Simulate Spatial Habitat and Species Distribution Data
#'
#' Generates virtual habitat configurations, conversion costs, and species distributions
#' for testing spatial conservation and species distribution models.
#' The simulation supports scalable grid dimensions, flexible numbers of habitats and species,
#' and realistic species dispersal distances drawn from a gamma distribution. It also supports
#' generating a high proportion of rare (low-prevalence) species.
#'
#' @param dim1 Integer. Number of rows in the spatial grid.
#' @param dim2 Integer. Number of columns in the spatial grid.
#' @param n_h Integer. Number of distinct habitat types.
#' @param n_s Integer. Number of species.
#' @param disp_max Integer. maximum allowable dispersal distances for species.
#' @param disp_longtail Numeric (0–1). Controls the proportion of long-distance dispersers. Low values generate mostly short dispersers; high values generate more long-tailed (dispersive) species. Default is 0.5.
#' @param rarity_bias Numeric (≥0). Controls how many species are rare (low-prevalence). Set to 0 for equal prevalence across species; higher values (e.g. 1–3) produce more rare species. Default is 1.
#' @param seed Integer. Random seed for reproducibility. Default is 1.
#'
#' @import terra
#' @import gstat
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{dim}}{Vector of grid dimensions \code{(dim1, dim2)}.}
#'   \item{\code{n_h}}{Number of habitats.}
#'   \item{\code{n_s}}{Number of species.}
#'   \item{\code{E}}{3D array \code{[dim1, dim2, n_h]} of binary habitat layers.}
#'   \item{\code{C}}{3D array \code{[dim1, dim2, n_h]} of habitat conversion costs.}
#'   \item{\code{U}}{3D array \code{[dim1, dim2, n_h]} indicating unavailable cells (1 = unavailable).}
#'   \item{\code{Y}}{3D array \code{[dim1, dim2, n_h]} of habitat configurations after accounting for availability.}
#'   \item{\code{SD}}{3D array \code{[dim1, dim2, n_s]} of binary species presence/absence values.}
#'   \item{\code{pd}}{3D array \code{[dim1, dim2, n_s]} of species occurrence probabilities.}
#'   \item{\code{S}}{Matrix \code{[n_h, n_s]} indicating which species use which habitats.}
#'   \item{\code{Disp}}{Numeric vector of species dispersal distances.}
#'   \item{\code{O}}{Numeric vector of species occupancy targets.}
#'   \item{\code{Alpha}}{Vector of penalty multipliers.}
#'   \item{\code{sigma}}{Scalar penalty weight based on mean dispersal distance.}
#'   \item{\code{seed}}{The seed used for random number generation.}
#' }
#'
#' @examples
#' sim <- simulate_ebrel_spatial_data(dim1 = 50, dim2 = 50, n_h = 3, n_s = 15, rarity_bias = 2)
#' hist(sim$Disp)  # Dispersal distances
#' image(sim$E[,,1])  # Plot first habitat layer
#' mean(colSums(sim$SD, dims = 2)) / (50 * 50)  # Average prevalence per species
#'
#' @export
simulate_ebrel_spatial_data <- function(
    dim1,
    dim2,
    n_h,
    n_s,
    disp_max,
    disp_longtail = 0.5,
    rarity_bias = 1,
    fixed_O = NULL) {

  # --- Number of cells
  n_cells <- dim1 * dim2

  gamma_params <- .map_longtail_to_gamma(disp_longtail)
  long_disp_shape <- gamma_params$shape
  long_disp_scale <- gamma_params$scale

  # --- Generate an autocorrelated response variable
  hab_rast <- .generate_habitat_rast(
    dim1 = dim1,
    dim2 = dim2,
    n_h = n_h,
    unavail_hab_prop = 0.25
  )

  # --- Generate rarity weights using beta distribution
  rarity_weights <- sort(rbeta(n_s, shape1 = 1, shape2 = 3 * rarity_bias + 0.1), decreasing = TRUE)

  # --- Simulated species probability distributions (pd)
  pd <- array(0, dim = c(dim1, dim2, n_s))

  # ---  Habitat association matrix
  SxH <- matrix(0, nrow = n_s, ncol = n_h)
  rownames(SxH) <- paste0("sp_", 1:n_s)
  colnames(SxH) <- paste0("hab_", 1:n_h)

  # --- Species dispersal distances ---
  D <- floor(rgamma(n_s, shape = long_disp_shape, scale = long_disp_scale))
  D <- pmax(1, pmin(D, disp_max))

  # --- Binary SD array
  SD_rast <- vector("list", n_s)

  for (s in 1:n_s) {
    npres <- 0
    tries <- 0L

    # cap habitats drawn and align probs
    k_max <- min(4L, n_h)
    hab_probs <- c(0.5, 0.25, 0.125, 0.125)[seq_len(k_max)]

    repeat {
      tries <- tries + 1
      if (tries > 200L) stop("simulate_ebrel_spatial_data(): failed to place species ", s,
                             " after 200 attempts (no eligible habitat cells).")

      # draw how many habitats and which ones
      k  <- sample.int(k_max, 1, prob = hab_probs)
      hs <- sample.int(n_h, k, replace = FALSE)

      # cells with any of the chosen habitats
      hs_max <- terra::app(hab_rast$E[[hs]], "max")
      # sample seed cells ONLY from >0
      ssize <- sample.int(5L, 1)
      cells <- .safe_sample_cells(hs_max, ssize)
      if (length(cells) == 0L) next  # try new hs

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

  # Set base costs for land management type per pixel
  hab_base_cost <- .approx_doubling(n_h, include_zero = FALSE, round_digits = 0)

  # Set landscape ramp to scale base costs by geographic area
  x_ramp_vec <- seq(0.1, 1, length.out = dim2)
  geo_C_ramp <- setValues(hab_rast$HO[[1]], rep(x_ramp_vec, dim1))

  # Assign costs
  C_rast <- lapply(seq_len(n_h), function(i) {
    C_i <- hab_rast$HO[[i]] # Get opportunity
    BC_i <- hab_base_cost[[i]] # Get base cost
    C_i <- ifel(C_i == 1, BC_i, 0)
    C_i <- C_i * geo_C_ramp
    C_i[C_i == 0] <- 1e10
    C_i
  })
  C_rast <- do.call(c, C_rast)
  names(C_rast) <- paste("C", 1:n_h)

  # Sort out NAs
  E_rast <- ifel(is.na(hab_rast$E), 0, hab_rast$E)
  SD_rast <- ifel(is.na(SD_rast), 0, SD_rast)
  C_rast <- ifel(is.na(C_rast), 1e10, C_rast)

  # Return
  list(
    dim = c(dim1, dim2),
    n_h = n_h,
    n_s = n_s,
    E = E_rast,
    C = C_rast,
    SD = SD_rast,
    SxH = SxH,
    D = D,
    O = O,
    sigma = sigma
  )
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
#' @param dim1 Integer. Number of rows in the grid. Default \code{50}.
#' @param dim2 Integer. Number of columns in the grid. Default \code{50}.
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
#' sim <- generate_habitats_rast(dim1 = 30, dim2 = 40, n_h = 3, unavail_hab_prop = 0.2)
#' sim$HP; sim$E; sim$HO
#' terra::global(sim$E, "sum", na.rm = TRUE)  # cells per habitat
#' }
#'
#' @importFrom gstat vgm gstat
#' @importFrom terra rast values app setValues which.max ifel global minmax ncell names
#' @export
.generate_habitat_rast <- function(dim1 = 50, dim2 = 50, n_h = 4,
                                   unavail_hab_prop = 0.25) {

  # --- Create autocorrelated surface
  xy <- expand.grid(x = 1:dim1, y = 1:dim2)
  varioMod <- vgm(psill = 1, range = 40, model = 'Exp', nugget = 0.0)
  zDummy <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE,
                beta = 0, model = varioMod, nmax = 50)

  # -- Simulate surface for "urban"
  urban_sim <- predict(zDummy, newdata = xy, nsim = 1, debug.level = 0)
  urban <- .rescale01(rast(urban_sim, type = "xyz"))
  urban_q <- quantile(values(urban), prob = unavail_hab_prop)
  urban_binary <- ifel(urban < urban_q, 1, 0)

  # -- Simulate surfaces for "habitat"
  habitat_sim <- predict(zDummy, newdata = xy, nsim = n_h, debug.level = 0)
  habitat <- rast(habitat_sim, type = "xyz")
  for (x in seq_len(n_h)) {
    habitat[[x]] <- .rescale01(habitat[[x]])
  }
  habitat <- ifel(urban_binary == 1, NA, habitat)

  # -- Establish habitat opportunities (HO)
  HO <- habitat
  for (h in 1:(n_h)) {
    qX <- quantile(values(HO[[h]]), prob = pmin(rnorm(1, 0.90,0.01), 0.99), na.rm = TRUE
                   ) # areas mod to small
    HO[[h]] <- ifel(HO[[h]] < qX, 0, 1) # also mask urban
  }

  # -- Establish existing habitat
  # - One-hot encode E (i.e. only one habitat per tile)
  habitat_max <- which.max(habitat)
  E_hab <- rast(lapply(seq_len(n_h), function(k) 1 * (habitat_max == k)))

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

  # - Set names
  names(habitat) <- paste0("hp", 1:n_h)
  names(E) <- paste0("E", 1:n_h)
  names(HO) <- paste0("HO", 1:n_h)

  list(
    HP = habitat,
    E = E,
    HO = HO
  )
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

#' Map a long-tail control to Gamma parameters
#'
#' Converts a unit-interval "long-tail" control into parameters of a
#' \eqn{\mathrm{Gamma}(\text{shape}, \text{scale})} distribution used to draw
#' dispersal distances. Larger \code{longtail} yields smaller shape (heavier
#' tail) and larger scale.
#'
#' Mapping:
#' \deqn{\text{shape} = \text{base\_shape} \cdot (1 - \text{longtail}) + 0.5}
#' \deqn{\text{scale} = \text{base\_scale} \cdot (1 + 4 \cdot \text{longtail})}
#'
#' @param longtail Numeric in \[0, 1]; 0 = short-tailed, 1 = long-tailed.
#' @param base_shape Positive numeric; baseline shape when \code{longtail = 0}.
#' @param base_scale Positive numeric; baseline scale when \code{longtail = 0}.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{shape}: gamma shape parameter.
#'   \item \code{scale}: gamma scale parameter.
#' }
#'
#' @examples
#' map_longtail_to_gamma(0)   # near base_shape + 0.5, base_scale
#' map_longtail_to_gamma(0.5) # intermediate shape/scale
#' map_longtail_to_gamma(1)   # shape ~ 0.5, scale ~ 5 * base_scale
#'
#' @keywords internal
#' @noRd
.map_longtail_to_gamma <- function(longtail, base_shape = 5, base_scale = 1) {
  inv <- 1 - longtail
  shape <- base_shape * inv + 0.5
  scale <- base_scale * (1 + 4 * longtail)
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
