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
#' Each named variable has a result in `result$analyses`. Its overall
#' measures are stored in `indices$overall` and requested pairwise
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
#' @param scope Whether to calculate "overall" measures, "pairwise"
#'   two-group measures, or "both" for categorical/ordinal categories and resolved
#'   continuous global groups. Ordinal pairwise indices compare individual
#'   categories using the same pair-conditioned formulas as categorical input.
#'   Overall indices use the full categorical distribution, ordinal cumulative
#'   thresholds, or continuous rank distribution. Continuous pairwise indices
#'   use binary group shares within the original full-variable neighborhoods,
#'   independently of summary mode. Exposure is selected through `measures`
#'   independently of scope. This global default can be overridden per variable.
#' @param output Return mode. "indices" stores only segregation indices.
#'   "full" stores both indices and local environments. "localenv"
#'   stores local environments without calculating indices.
#'   Categorical and ordinal environments default to category shares; continuous
#'   variables default to shares of resolved global quartile groups. Use
#'   `summary = list(type = "quantiles")` for synchronized local/grouped rank
#'   summaries, with `mean_rank = TRUE` in that list for a separate mean.
#'   Customize continuous summaries with `summary`, or ordinal summaries with
#'   `localenv`, independently of full-resolution scalar indices.
#'   This global default can be overridden per variable. Shared input geometry,
#'   coordinates, data, and population are retained once if any variable stores
#'   environments.
#' @param keep_scores Logical. Retain per-location H/R/S local variation and
#'   D discrepancies, for the requested measures and scopes. Extract with
#'   `as.data.frame(result, what = "scores")`.
#' @param keep_contributions Logical. Retain additive per-location index
#'   contributions, extracted with `what = "contributions"`. They sum to the
#'   corresponding regional index, using original focal-population weights and
#'   pair-specific weights for pairwise results. Both retention options default
#'   to `FALSE` and work independently of `output`; only requested arrays are
#'   allocated. Exposure and rank profiles have separate output.
#'   Retention memory grows with locations, bands, requested scalar measures,
#'   and group pairs; each enabled option stores its own arrays.
#'   Both retention options can be overridden per variable, including by `FALSE`.
#' @param verbose Logical. If `TRUE`, report progress messages.
#' @param variables Column names for a single categorical variable, or a uniquely
#'   named list of [spseg_variables] specifications. A single column name is
#'   interpreted as observed category memberships; two or more column
#'   names are interpreted as category counts (aggregated data). A message reports
#'   this assumption. The result is named after the single column, or `categorical`
#'   for multiple count columns. Use an explicit specification for other variable
#'   types, a custom result name, or population weights for long-form counts.
#'   Single-column count input is not supported. All variables share spatial
#'   discovery; kNN uses each variable's eligible population and adaptive
#'   distance weights. Each specification can override `measures`.
#'   Grid and pycnophylactic surfaces require count specifications for all
#'   variables and polygonal sf input; raw support also accepts observed values.
#'   See `vignette("spseg-measure-definitions")` for formulas and sources,
#'   and `vignette("spseg-multiple-variables")` for workflows.
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
#' In shares mode the exposure/isolation measure, P, is presented in matrix form. The spatial
#' exposure of group 'm' to group 'n' is located in row 'm' and column 'n'.
#' The matrix is rarely symmetric in practice.
#' In continuous quantiles mode, `measures = "exposure"` requests focal-group
#' averages of local rank statistics, stored as `rank_summary` rather than P.
#' These descriptive profiles are neither P matrices nor segregation indices.
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
#' variables <- list(ethnicity = spseg_categorical(vars = c("A1", "A2")))
#' result <- spseg(xy, variables, data = spseg_data, bands = c(0, 2), output = "full")
#' as.data.frame(result)
#' head(as.data.frame(result, what = "env"))
#' result$analyses$ethnicity$env[[1]]
#'
#' # Shorthand for one categorical variable stored as category counts:
#' spseg(xy, c("A1", "A2"), data = spseg_data, bands = c(0, 2))
#' # For microdata, use variables = "race" to select observed memberships.
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
  scope = c("both", "overall", "pairwise"),
  output = c("indices", "full", "localenv"),
  keep_scores = FALSE,
  keep_contributions = FALSE,
  verbose = FALSE,
  network = NULL,
  node_id = NULL,
  ...
) {
  if (missing(variables)) {
    stop("Supply column names or a named 'variables' list of spseg_variable specifications",
         call. = FALSE)
  }
  if (is.character(variables) && is.null(dim(variables))) {
    if (!length(variables) || anyNA(variables) || any(!nzchar(variables)) ||
        anyDuplicated(variables)) {
      stop("'variables' column names must be non-empty, non-missing, and unique",
           call. = FALSE)
    }
    columns <- unname(variables)
    if (length(columns) == 1L) {
      variables <- setNames(list(spseg_categorical(vars = columns)), columns)
      message("Assuming categorical observed memberships: '", columns,
              "' contains observed category memberships; result named '",
              columns, "'.")
    } else {
      variables <- list(categorical = spseg_categorical(vars = columns))
      message("Assuming aggregated categorical data: selected columns contain ",
              "category counts for one variable; result named 'categorical'.")
    }
  }
  call <- match.call()
  dots <- list(...)
  if (!is.list(variables) || !length(variables) || is.null(names(variables)) ||
      anyNA(names(variables)) || any(!nzchar(names(variables))) ||
      anyDuplicated(names(variables)) ||
      !all(vapply(variables, inherits, logical(1), "spseg_variable")))
    stop("'variables' must be a uniquely named list of spseg_variable specifications",
         call. = FALSE)
  neighbors <- match.arg(neighbors, c("radius", "knn"))
  search <- match.arg(search, c("boost", "brute"))
  weighting <- match.arg(weighting, c("biweight", "unweighted", "inverse", "exponential"))
  surface <- match.arg(surface, c("raw", "grid", "pycno"))
  relational <- !is.null(network)
  if (relational && surface != "raw")
    stop("Network input requires surface = 'raw'", call. = FALSE)
  scope <- match.arg(scope, c("both", "overall", "pairwise"))
  output <- match.arg(output, c("indices", "full", "localenv"))
  retention <- list(keep_scores = keep_scores,
                    keep_contributions = keep_contributions)
  for (option in names(retention)) {
    value <- retention[[option]]
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop("'", option, "' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(surface_audit) || length(surface_audit) != 1L || is.na(surface_audit))
    stop("'surface_audit' must be TRUE or FALSE", call. = FALSE)
  if (surface == "raw" && length(dots))
    stop("Surface options require surface = 'grid' or 'pycno'", call. = FALSE)
  if (surface != "raw" && any(vapply(variables, function(z) is.null(z$counts), logical(1))))
    stop("Grid and pycnophylactic surfaces require category counts in multiple 'vars' columns for every variable",
         call. = FALSE)
  if (length(power) != 1L || !is.numeric(power) || !is.finite(power) || power <= 0)
    stop("'power' must be a positive finite number", call. = FALSE)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
    stop("'normalize' must be TRUE or FALSE", call. = FALSE)
  prepared_network <- NULL
  if (relational) {
    if (search != "boost") stop("Network input requires search = 'boost'", call. = FALSE)
    if (is.null(bands)) stop("Supply explicit 'bands' for network neighborhoods", call. = FALSE)
    prepared_network <- .spseg_network_input(x, data, coords, network, node_id)
    checked <- prepared_network$checked
  } else {
    if (!is.null(node_id)) stop("'node_id' requires a network", call. = FALSE)
    checked <- .spseg_geometry(x, data, coords)
  }
  count_data <- lapply(variables, function(spec) {
    if (!is.null(spec$counts)) .spseg_counts(spec, checked$data) else NULL
  })
  supplied_population <- !is.null(population)
  population <- if (neighbors == "knn" || relational) {
    .spseg_mass(population, checked$data, nrow(checked$data))
  } else .spseg_population(population, variables, checked$data, count_data)
  if (neighbors == "radius" && (!is.finite(sum(population)) || sum(population) <= 0))
    stop("The common population total must be positive and finite", call. = FALSE)
  compiled <- Map(function(spec, name, counts) {
    .spseg_compile(spec, checked$data, population, measures, variable = name,
                   counts = counts)
  }, variables, names(variables), count_data)
  surface_info <- list(type = "raw", fallback = list(count = 0L), audit = NULL)
  if (surface != "raw") {
    # Gridding transports all count tables in one surface-construction pass.
    # The first eligible population is a sufficient reference for kNN because
    # its thresholds use each descriptor's mass, not this reference vector.
    if (neighbors == "knn" && !supplied_population) population <- compiled[[1L]]$mass
    prepared <- .spseg_surface_variables(x, checked, variables, compiled, population,
                                       measures, surface, dots, surface_audit, verbose)
    checked <- prepared$checked
    population <- prepared$population
    compiled <- prepared$compiled
    surface_info <- prepared$surface_info
  }
  if (is.null(bands)) {
    if (neighbors == "knn") stop("Supply population-count 'bands' for kNN", call. = FALSE)
    bands <- spseg_default_bands(checked$coords, population = population)
  }
  if (!is.numeric(bands) || !length(bands) || any(!is.finite(bands)) ||
      any(bands < 0) || (neighbors == "knn" && any(bands <= 0)))
    stop("Supply finite nonnegative radius bands or positive population-count kNN bands",
         call. = FALSE)
  if (neighbors == "knn") {
    for (name in names(compiled)) {
      available <- sum(compiled[[name]]$mass)
      # Allow accumulated rounding of fractional weights in R's sum; the
      # native selector separately verifies the population actually selected.
      tolerance <- (length(compiled[[name]]$mass) + 16) *
        .Machine$double.eps * max(bands)
      if (max(bands) > available + tolerance)
        stop("Variable '", name, "' has eligible population ", format(available),
             ", below the requested kNN threshold ", format(max(bands)), call. = FALSE)
    }
  }
  if (verbose) message("Calculating ", length(compiled), " variables with shared neighborhood searches")
  compiled <- Map(function(config, spec) {
    config$settings <- list(scope = spec$scope %||% scope,
      output = spec$output %||% output,
      keep_scores = spec$keep_scores %||% keep_scores,
      keep_contributions = spec$keep_contributions %||% keep_contributions)
    config
  }, compiled, variables)
  store_input <- any(vapply(compiled, function(z)
    z$settings$output != "indices", logical(1)))
  store_local_keys <- store_input || any(vapply(compiled, function(z)
    z$settings$keep_scores || z$settings$keep_contributions, logical(1)))
  engine <- spseg_engine_cpp(checked$coords %||% matrix(numeric(), 0L, 2L),
                         compiled, as.numeric(bands),
                         power, .spseg_weighting_id(weighting), normalize,
                         match(neighbors, c("radius", "knn")) - 1L,
                         match(search, c("boost", "brute")) - 1L,
                         if (relational) prepared_network$graph else list())
  analyses <- Map(.spseg_prepare_result, engine$analyses, compiled,
                  MoreArgs = list(bands = bands))
  result <- list(analyses = analyses, bands = as.numeric(bands),
                 coords = if (store_input) checked$coords else NULL,
                 data = if (store_input) checked$data else NULL,
                 geometry = if (store_input) checked$geometry else NULL,
                 population = if (store_input) population else NULL,
                 crs = checked$crs, weighting = weighting, power = power,
                 normalize = normalize, surface = surface, output = output,
                 keep_scores = keep_scores, keep_contributions = keep_contributions,
                 surface_info = surface_info,
                 network = if (relational) prepared_network$summary else NULL,
                 node_id = if (relational && store_local_keys) prepared_network$ids else NULL,
                 neighbors = list(type = neighbors, values = bands,
                                  space = if (relational) "network" else "euclidean",
                                  units = if (neighbors == "knn") "population" else
                                    if (relational) prepared_network$summary$units else "distance",
                                  distance_units = if (relational) prepared_network$summary$units else "coordinate",
                                  population = if (neighbors == "knn") "variable" else NULL,
                                  engine = search, scope = scope),
                 diagnostics = engine$diagnostics, call = call)
  class(result) <- "spseg_result"
  result
}
