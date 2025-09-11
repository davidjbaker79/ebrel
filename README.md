# **_ebrel_: spatial prioritisation for nature recovery**

R package for running the Ebrel spatial prioritisation model with multiple options for habitat creation. Given multiple candidate habitat types and costs, ebrel identifies a configuration that maximises biodiversity objectives subject to constraints, while accommodating species with different dispersal capacities and connectivity needs.
________________________________________
Overview
ebrel provides fast, reproducible tools to run spatial prioritisation using the Ebrel model, as described in:

> Unnithan Kumar, S., Baker, D. J., Maclean, I. M. D., & Gaston, K. J. (2025). Spatial prioritisation for nature recovery with multiple options for habitat creation. Journal of Applied Ecology, 00, 1–13. https://doi.org/10.1111/1365-2664.70144

The package implements the same objective functions as Unnithan Kumar et al. (2025), while optimising parts of the implementation for computational efficiency. In particular, ebrel:
* Rewrites key routines in C++ for speed (via Rcpp).
* Uses a Breadth-First Search (BFS)-based accessibility computation instead of the original hierarchical clustering approach. This reduces typical complexity from ~O(n²) to ~O(V+E) per species on the underlying adjacency graph, and scales better on large grids.
* Provides convenience functions to simulate data, prepare inputs, and run the optimiser end-to-end.
________________________________________
# Installation
### Install from GitHub

install.packages("remotes")
remotes::install_github("davidjbaker79/ebrel", build_vignettes = TRUE)

System requirements
* R ≥ 4.2 recommended
* *C++17-capable toolchain (used indirectly by R)
* Windows: Rtools (matching your R version)

________________________________________
# Quick start

See vignette for example.

________________________________________
# How it works (high level)
* Objective functions: Identical to Unnithan Kumar et al. (2025), combining habitat suitability, species distributions, conversion costs, and connectivity terms into a single score HH.
* Connectivity & accessibility: Computed via BFS on a rook/queen adjacency graph to assess reachability of proposed habitat from current species’ distributions within species-specific dispersal thresholds.
* Dispersal tiers: Supports species with differing dispersal capacities (e.g., highly limited, moderate/stepping-stone, and good dispersers) through configurable connectivity logic.
* Search algorithm: A simulated-annealing-like metaheuristic proposes local changes (cell ↦ habitat type) using distance-weighted probabilities, with a Metropolis acceptance rule and early-stopping heuristics.
________________________________________
# Citing
If you use _ebrel_, please cite the paper and, if appropriate, the software:

### Paper
> Unnithan Kumar, S., Baker, D. J., Maclean, I. M. D., & Gaston, K. J. (2025). Spatial prioritisation for nature recovery with multiple options for habitat creation. Journal of Applied Ecology. https://doi.org/10.1111/1365-2664.70144

### Software
> Baker, D. J., et al. (2025). ebrel: Spatial prioritisation for nature recovery with multiple habitat options. R package version X.Y.Z. GitHub: https://github.com/davidjbaker79/ebrel
> 
@article{unnithankumar2025ebrel,
  title   = {Spatial prioritisation for nature recovery with multiple options for habitat creation},
  author  = {Unnithan Kumar, S. and Baker, D. J. and Maclean, I. M. D. and Gaston, K. J.},
  journal = {Journal of Applied Ecology},
  year    = {2025},
  doi     = {10.1111/1365-2664.70144}
}

@manual{ebrelRpkg,
  title  = {ebrel: Spatial prioritisation for nature recovery},
  author = {Baker, D. J. and collaborators},
  year   = {2025},
  note   = {R package version X.Y.Z},
  url    = {https://github.com/davidjbaker79/ebrel}
}
________________________________________
# Contributing
Bug reports and pull requests are welcome via GitHub Issues/PRs. For larger features, please open an issue first to discuss scope/design.
________________________________________
# License
This project’s license is provided in the LICENSE file of the repository.
________________________________________
# Contact
Questions or suggestions? Please open an issue on GitHub or contact @davidjbaker79.
