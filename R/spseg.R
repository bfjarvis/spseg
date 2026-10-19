#' Spatial Segregation Measures
#'
#' Calculates spatial segregation measures from point coordinates or [sf::sf]
#' geometries using local population environments.
#'
#' `spseg()` is the primary interface for the package. Local environments
#' and segregation indices are produced by a unified C++ engine. The default
#' calculation uses normalized biweight local environments with `power = 3`
#' and radius neighborhoods.
#'
#' By default, `spseg()` returns indices over default distance bands. These
#' default bands are the first through fifth deciles of distances among a
#' reproducible population-weighted sample of points; see [default_bands()].
#' If `output = "indices"`, coordinates, input data, and local environments
#' are not stored in the result.
#'
#' Set `bands = 0` to compute indices without spatial weighting. In that
#' case, each local environment is a copy of the observed population values
#' for the focal unit.
#'
#' When `neighbors = "knn"`, `bands` is interpreted as population-count
#' thresholds. The focal unit is included fully. If the focal unit's
#' population exceeds a threshold, its observed composition is used for that
#' threshold. Otherwise, neighboring units are added in order of centroid
#' distance until the threshold is reached, with the final boundary unit
#' included fractionally when needed.
#'
#' Multigroup measures are stored in `result$indices$multigroup`. Pairwise
#' measures, when requested, are stored in `result$indices$pairwise`.
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
#' @param data Optional attribute data with one row per coordinate or feature.
#'   Legacy calls require numeric population counts in at least two columns.
#'   With `variables`, this may instead contain category memberships, observed
#'   numerical values, and counts selected by the specifications. If omitted,
#'   attributes come from `x` (excluding sf geometry).
#' @param measures A character vector naming measures to compute: "all",
#'   "exposure", "information", "diversity", "dissimilarity", or "square_root".
#'   See `variables` for aliases and the legacy/new meaning of `"all"`.
#' @param bands Neighborhood thresholds. With `neighbors = "radius"`, these
#'   are distance bandwidths in the units of `x`. With `neighbors = "knn"`,
#'   these are population-count thresholds. If omitted for radius
#'   neighborhoods, [default_bands()] is used.
#' @param weighting Spatial weighting scheme used within each neighborhood.
#'   One of "biweight", "unweighted", "inverse", or "exponential".
#' @param power Power parameter for the biweight and inverse-distance
#'   weighting schemes.
#' @param normalize Logical. If `TRUE`, distances are divided by the
#'   bandwidth before weights are calculated. The bandwidth still defines
#'   the search radius when `neighbors = "radius"`.
#' @param neighbors Neighborhood strategy. "radius" uses distance
#'   bandwidths. "knn" uses count-based adaptive neighborhoods.
#' @param search Neighbor-search backend. "kdtree" uses the exact
#'   nanoflann kd-tree backend. "brute" computes candidate distances
#'   directly and is mainly useful for validation.
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
#' @param verbose Logical. If `TRUE`, report progress messages.
#' @param variables Optional uniquely named list of [seg_variables]
#'   specifications. All variables share neighborhood searches and spatial
#'   weights. The new interface supports `surface = "raw"`. The measure
#'   `"square_root"` and aliases `"relative_diversity"`/`"variance_ratio"` are
#'   available. With `variables`, `"all"` includes the root family; legacy
#'   `"all"` retains its existing four-family selection. Each specification can
#'   override `measures`. See `vignette("seg-multiple-variables")` for formulas,
#'   estimator conventions, sources, and properties.
#' @param coords For the multiple-variable interface with plain data frames,
#'   two coordinate column names or positions. For sf input, omit this argument
#'   to use the active geometry. Distances are Euclidean in coordinate units;
#'   use an appropriate projected CRS for metric bandwidths.
#' @param population Common neighborhood population column or numeric vector
#'   for the multiple-variable interface. Defaults to one per individual row,
#'   or consistent count/weight totals. Supply it explicitly when variable
#'   populations differ. Each variable retains its own eligible population for
#'   index averaging and its reference distribution. kNN thresholds use this
#'   common population, not a sum across variable specifications.
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
#' @return An S3 object of class `seg_result`. It is a list with fields
#'   including `bands`, `indices`, `env`, `coords`, `data`,
#'   `geometry`, `neighbors`, `surface_info`, and `crs`. Some fields are
#'   `NULL` depending on `output`. Calls with `variables` return a
#'   [seg_multi_result] inheriting from `seg_result`, with named results under
#'   `analyses` and shared geometry/neighborhood metadata.
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
#' @seealso [default_bands()], [as.data.frame.seg_result()]
#'
#' @examples
#' data(segdata)
#'
#' coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
#' test_data <- segdata[, 1:2]
#'
#' # Default behavior: indices across default bands.
#' idx <- spseg(coords, data = test_data)
#' idx
#' as.data.frame(idx)
#'
#' # Store local environments for selected distance bands.
#' env <- spseg(coords, data = test_data, bands = c(0, 2), output = "localenv")
#' head(as.data.frame(env, what = "env", shape = "long"))
#'
#' # Store indices and local environments together.
#' full <- spseg(coords, data = test_data, bands = c(1, 2), output = "full")
#' as.data.frame(full, what = "exposure")
#'
#' # sf inputs preserve geometry in local-environment outputs.
#' geom <- sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10),
#'   c(0,10), c(0,0)))) |>
#'   sf::st_make_grid(cellsize = 1)
#' grd_sf <- sf::st_sf(test_data, geometry = geom)
#'
#' sf_env <- spseg(grd_sf, bands = 2, output = "localenv")
#' sf::st_as_sf(sf_env)
#'
#' @export
spseg <- function(
  x,
  data = NULL,
  measures = "all",
  bands = NULL,
  weighting = c("biweight", "unweighted", "inverse", "exponential"),
  power = 3,
  normalize = TRUE,
  neighbors = c("radius", "knn"),
  search = c("kdtree", "brute"),
  surface = c("raw", "grid", "pycno"),
  surface_audit = FALSE,
  scope = c("both", "multigroup", "pairwise"),
  output = c("indices", "full", "localenv"),
  verbose = FALSE,
  variables = NULL,
  coords = NULL,
  population = NULL,
  ...
) {
  call <- match.call()
  if (!is.null(variables)) {
    return(.spseg_multiple(
      x, data, variables, coords, population, measures, bands,
      if (missing(weighting)) "biweight" else weighting,
      power, normalize, neighbors, search, surface, surface_audit,
      scope, output, verbose, call, list(...)
    ))
  }
  if (!is.null(coords) || !is.null(population)) {
    stop("'coords' and 'population' require 'variables'", call. = FALSE)
  }
  output <- match.arg(output)
  scope <- match.arg(scope)
  neighbors <- match.arg(neighbors)
  search <- match.arg(search)
  surface <- match.arg(surface)
  surface_audit <- isTRUE(surface_audit)
  dots <- list(...)
  weighting <- if (missing(weighting)) "biweight" else match.arg(weighting)
  normalize <- isTRUE(normalize)

  if (neighbors == "knn" && is.null(bands)) {
    stop(
      "'bands' must be supplied for count-based kNN neighborhoods",
      call. = FALSE
    )
  }

  if (!is.null(bands)) {
    bands <- as.numeric(bands)
  }
  if (!is.null(bands) && any(!is.finite(bands))) {
    stop("'bands' must contain finite numeric values", call. = FALSE)
  }
  if (!is.null(bands) && neighbors == "radius" && any(bands < 0)) {
    stop("'bands' must be greater than or equal to 0", call. = FALSE)
  }
  if (!is.null(bands) && neighbors == "knn" && any(bands <= 0)) {
    stop("'bands' must be greater than 0 for kNN neighborhoods", call. = FALSE)
  }

  checked <- if (verbose) {
    if (is.null(data)) chksegdata(x) else chksegdata(x, data)
  } else {
    if (is.null(data)) {
      suppressMessages(chksegdata(x))
    } else {
      suppressMessages(chksegdata(x, data))
    }
  }
  pop_surface <- spseg_surface(
    x = x,
    coords = checked$coords,
    data = checked$data,
    surface = surface,
    args = dots,
    audit = surface_audit,
    verbose = verbose
  )
  pop_surface$coords <- as.matrix(pop_surface$coords)
  pop_surface$data <- as.matrix(pop_surface$data)

  if (is.null(bands)) {
    bands <- default_bands(pop_surface$coords, pop_surface$data)
  }

  measures_to_compute <- if (identical(output, "localenv")) {
    character()
  } else {
    spseg_measures(measures)
  }
  engine <- if ("square_root" %in% measures_to_compute) {
    .seg_legacy_engine(pop_surface$coords, pop_surface$data, bands, power,
                       weighting, normalize, measures_to_compute, scope,
                       neighbors, search, output)
  } else seg_engine_coords(
    coords = pop_surface$coords,
    data = pop_surface$data,
    bands = bands,
    power = power,
    weighting = weighting,
    normalize = normalize,
    measures = measures_to_compute,
    scope = scope,
    neighbors = neighbors,
    search = search,
    keep_env = !identical(output, "indices"),
    keep_indices = !identical(output, "localenv")
  )
  indices <- if (identical(output, "localenv")) {
    list()
  } else {
    spseg_indices_from_engine(engine$indices, bands)
  }
  result <- spseg_result(
    coords = pop_surface$coords,
    data = pop_surface$data,
    env = engine$env,
    bands = bands,
    indices = indices,
    measures = if (!identical(output, "localenv")) {
      spseg_measures(measures)
    } else {
      character()
    },
    scope = scope,
    weighting = weighting,
    power = power,
    normalize = normalize,
    crs = st_crs(checked$proj4string),
    output = output,
    call = call,
    geometry = if (identical(surface, "raw")) {
      checked$geometry
    } else {
      pop_surface$geometry
    },
    neighbors = neighbors,
    search = search,
    surface = surface,
    surface_info = pop_surface$surface_info
  )
  result
}
