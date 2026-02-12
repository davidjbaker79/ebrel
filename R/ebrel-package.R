#' @keywords internal
#' @useDynLib ebrel, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

#' Generate naive starting point
#'
#' @export
generate_X0_A_R <- function(U, n_h, dim_x, dim_y, base_prob = 0.85, seed) {
  .Call("_ebrel_generate_X0_A_R", PACKAGE = "ebrel", U, n_h, dim_x, dim_y, base_prob, seed)
}

