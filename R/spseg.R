#' Spatial and Network Segregation Measures
#'
#' Calculates spatial segregation measures from point coordinates or [sf::sf]
#' geometries, or a supplied network, using local population environments.
#'
#' `spseg()` is the primary interface for the package. Local environments
#' and segregation indices are produced by a unified C++ engine. The default
#' calculation uses normalized biweight local environments with `power = 3`
#' and radius neighborhoods.
#'
#' By default, `spseg()` returns indices over default distance bands. These
#' default bands are distance quantiles from a reproducible population-weighted
#' sample of support points; see [spseg_default_bands()] for the quantile defaults.
#' If `output = "indices"`, coordinates, input data, and local environments
#' are not stored in the result.
#'
#' Set `bands = 0` to compute indices without spatial weighting. In that
#' case, local environments contain the observed composition at coincident
#' support points, or just the focal unit when coordinates are distinct.
#'
#' When `neighbors = "knn"`, `bands` is interpreted as population-count
#' thresholds before distance weighting. Each variable
#' selects the requested eligible population, truncating the focal unit when
#' needed and including equidistant boundary units proportionately. Neighboring
#' units are added by centroid distance. An unattainable threshold raises an error.
#'
#' Each named variable has a result in `result$analyses`. Its multigroup
#' measures are stored in `indices$multigroup` and requested pairwise
#' measures in `indices$pairwise`.
#' Pairwise scalar measures are returned as one symmetric matrix per
#' bandwidth with one row and column for each population group. Units with
#' no population in either group for a given pair do not contribute to that
#' pair's accumulators.
#'
#' The `surface` argument controls population-surface construction before
#' local environments are calculated. `surface = "grid"` redistributes
#' areal counts over a square grid. `surface = "pycno"` applies internal
#' groupwise pycnophylactic interpolation: each group is smoothed separately
#' by averaging each grid cell with neighboring cells, then rescaled within
#' each source polygon after every iteration so known group totals are
#' preserved.
#'
#' When computing the spatial information theory index, zero proportions are
#' evaluated using the limiting value $x \\log(x) = 0$ as $x$ approaches
#' zero. This avoids artificial perturbation of the data and prevents
#' `0 * log(0)` from producing `NaN`.
#'
#' @param x A numeric matrix or data frame with two coordinate columns, or
#'   an object of class [sf::sf]. For polygon [sf::sf] inputs, centroids
#'   are used as the focal coordinates. With `variables`, a plain data frame
#'   may contain both attributes and the columns selected by `coords`.
#'   With `network`, supply an observation table; coordinates are optional.
#' @param data Optional attribute data with one row per coordinate or feature.
#'   This may contain category memberships, observed
#'   numerical values, and counts selected by the specifications. If omitted,
#'   attributes come from `x` (excluding sf geometry).
#' @param measures A character vector naming measures to compute: "all",
#'   "exposure", "information", "diversity", "dissimilarity", or "square_root".
#'   `"all"` includes all five families. Aliases `"relative_diversity"` and
#'   `"variance_ratio"` select diversity.
#' @param bands Neighborhood thresholds. With `neighbors = "radius"`, these
#'   are distance bandwidths in the units of `x`. With `neighbors = "knn"`,
#'   these are population-count thresholds. If omitted for radius
#'   neighborhoods, [spseg_default_bands()] is used.
#' @param weighting Spatial weighting scheme used within each neighborhood.
#'   One of "biweight", "unweighted", "inverse", or "exponential".
#' @param power Power parameter for the biweight and inverse-distance
#'   weighting schemes.
#' @param normalize Logical. If `TRUE`, distances are divided by the
#'   bandwidth before weights are calculated. The bandwidth still defines
#'   the search radius when `neighbors = "radius"`.
#' @param neighbors Neighborhood strategy. "radius" uses distance
#'   bandwidths. "knn" uses count-based adaptive neighborhoods.
#' @param search Neighbor-search backend. "boost" uses the exact
#'   Boost.Geometry R-tree or Boost.Graph traversal. "brute" computes Euclidean
#'   distances directly and is mainly useful for validation.
#' @param surface Population surface used before local environments are
#'   calculated. "raw" uses the input support directly. "grid"
#'   redistributes areal counts over a square grid. "pycno" applies
#'   groupwise pycnophylactic smoothing on that grid while preserving source
#'   polygon group totals.
#' @param surface_audit Logical. If `TRUE`, retain auditable surface
#'   intermediates in `result$surface_info$audit`, including source counts and
#'   geometry, source-zone assignments for retained grid cells, fallback
#'   contributions, the initial gridded counts, and pycnophylactically smoothed
#'   counts when applicable. The default retains only lightweight grid and
#'   fallback metadata.
#' @param scope Whether to calculate "multigroup" measures, "pairwise"
#'   two-group measures, or "both" for categorical outcomes. Ordered outcomes
#'   always return their cumulative-threshold summaries; exposure is separate.
#' @param output Return mode. "indices" stores only segregation indices.
#'   "full" stores both indices and local environments. "localenv"
#'   stores local environments without calculating indices.
#'   Categorical environments contain category shares; ordered variables default
#'   to compact local quantiles of global ranks. Customize the retained summary
#'   with `localenv` in the variable constructors, independently of the indices.
#' @param verbose Logical. If `TRUE`, report progress messages.
#' @param variables Required uniquely named list of [spseg_variables]
#'   specifications, including for a single segregation dimension. A count
#'   specification can select several columns. All variables share spatial
#'   discovery; kNN uses each variable's eligible population and adaptive
#'   distance weights. Each specification can override `measures`.
#'   Grid and pycnophylactic surfaces require count specifications for all
#'   variables and polygonal sf input; raw support also accepts observed values.
#'   See `vignette("spseg-multiple-variables")` for formulas and sources.
#' @param coords With plain data frames,
#'   two coordinate column names or positions. For sf input, omit this argument
#'   to use the active geometry. Distances are Euclidean in coordinate units;
#'   use an appropriate projected CRS for metric bandwidths.
#' @param network Optional [spseg_network()] specification. If omitted,
#'   neighborhoods use Euclidean distances. Supplying a network uses shortest
#'   paths exclusively, regardless of any geometry. Network calculations require
#'   explicit `bands`, `search = "boost"`, and `surface = "raw"`. Radius bands
#'   are in hops or cost units; kNN bands remain variable-specific populations.
#'   Unreachable population does not count toward kNN; an unattainable target
#'   raises an error identifying the variable and focal observation. The global
#'   reference distribution still includes all supplied observations.
#' @param node_id Single attribute column name or position mapping each
#'   observation to a network node. Required with `network`; otherwise omit.
#' @param population Default population-weight column or numeric vector for
#'   observed-value specifications.
#'   Descriptor `weights` override it; count specifications use their row totals.
#'   kNN defaults to one per observed-value row and selects each variable's
#'   own eligible population after omissions. For radius neighborhoods,
#'   consistent count/weight totals can supply the default; specify this
#'   reference population explicitly when those totals differ. It also controls
#'   population-based default distance bands for radius neighborhoods.
#'   Network observed-value specifications default to one per row unless weights
#'   are supplied; count specifications always use their own row totals.
#' @param ... Additional population-surface options, such as `cellsize`,
#'   `celldim`, `nrow`, `ncol`, `max_iter`, `converge`,
#'   `surface_geometry`. Grid and pycnophylactic surfaces use
#'   [terra::rasterize()] to assign source polygons to generated grid-cell
#'   centroids. Non-empty polygons that receive no centroid-assigned cell
#'   contribute their counts to the cell containing an interior point of the
#'   polygon. Set `surface_geometry = "points"` to store centroid point
#'   geometries, or `surface_geometry = "polygons"`/`TRUE` to construct and
#'   store square grid-cell geometries for mapping. By default, generated
#'   surfaces store only centroid coordinates and population counts.
#'
#' @return A [spseg_result] with named variable results in `analyses`, plus
#'   shared `bands`, `neighbors`, `diagnostics`, and `surface_info`. Input data,
#'   coordinates, geometry, and environments are retained for `output = "full"`
#'   or `"localenv"`. Surface audits are retained when requested.
#'
#' @note
#' The exposure/isolation index, P, is presented in matrix form. The spatial
#' exposure of group 'm' to group 'n' is located in row 'm' and column 'n'.
#' The matrix is rarely symmetric in practice.
#'
#' @references
#' Reardon, S. F. and O'Sullivan, D. (2004) Measures of spatial segregation.
#'   *Sociological Methodology*, **34**, 121-162.
#'
#' Reardon, S. F., Farrell, C. R., Matthews, S. A., O'Sullivan, D., Bischoff,
#'   K., and Firebaugh, G. (2009) Race and space in the 1990s: Changes in the
#'   geographic scale of racial residential segregation, 1990-2000.
#'   *Social Science Research*, **38**, 55-70.
#'
#' @seealso [spseg_default_bands()], [as.data.frame.spseg_result()]
#'
#' @examples
#' data(spseg_data)
#' xy <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
#' variables <- list(ethnicity = spseg_categorical(counts = c("A1", "A2")))
#' result <- spseg(xy, variables, data = spseg_data, bands = c(0, 2), output = "full")
#' as.data.frame(result)
#' head(as.data.frame(result, what = "env"))
#' result$analyses$ethnicity$env[[1]]
#'
#' geom <- sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10),
#'   c(0,10), c(0,0)))) |> sf::st_make_grid(cellsize = 1)
#' blocks <- sf::st_sf(spseg_data, geometry = geom)
#' gridded <- spseg(blocks, variables, surface = "grid", nrow = 10, ncol = 10,
#'                  bands = 2, output = "localenv")
#' sf::st_as_sf(gridded)
#'
#' @export
spseg <- function(
  x,
  variables,
  coords = NULL,
  data = NULL,
  population = NULL,
  measures = "all",
  bands = NULL,
  weighting = c("biweight", "unweighted", "inverse", "exponential"),
  power = 3,
  normalize = TRUE,
  neighbors = c("radius", "knn"),
  search = c("boost", "brute"),
  surface = c("raw", "grid", "pycno"),
  surface_audit = FALSE,
  scope = c("both", "multigroup", "pairwise"),
  output = c("indices", "full", "localenv"),
  verbose = FALSE,
  network = NULL,
  node_id = NULL,
  ...
) {
  if (missing(variables))
    stop("Supply a named 'variables' list of spseg_variable specifications", call. = FALSE)
  .spseg_compute(x, data, variables, coords, population, measures, bands,
                 match.arg(weighting), power, normalize, neighbors, search,
                 surface, surface_audit, scope, output, verbose, match.call(), list(...),
                 network, node_id)
}
