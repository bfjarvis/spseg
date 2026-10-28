.spseg_geometry <- function(x, data, coords) {
  if (inherits(x, "sf")) {
    if (!is.null(coords)) stop("For sf input, use the active geometry instead of 'coords'",
                              call. = FALSE)
    geometry <- sf::st_geometry(x)
    xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_make_valid(geometry))))
    source <- if (is.null(data)) sf::st_drop_geometry(x) else as.data.frame(data)
    crs <- sf::st_crs(x)
  } else {
    geometry <- NULL
    crs <- sf::st_crs(NA_character_)
    source <- if (is.null(data)) as.data.frame(x) else as.data.frame(data)
    if (is.null(coords)) {
      xy <- as.matrix(x)
      if (ncol(xy) != 2L)
        stop("Supply two coordinate columns, or select them using 'coords'",
             call. = FALSE)
    } else xy <- as.matrix(.spseg_select(as.data.frame(x), coords))
  }
  if (!is.numeric(xy) || ncol(xy) != 2L || nrow(xy) < 1L ||
      nrow(xy) != nrow(source) || any(!is.finite(xy)))
    stop("Supply finite two-dimensional coordinates with one location per data row",
         call. = FALSE)
  storage.mode(xy) <- "double"
  list(coords = xy, data = source, geometry = geometry, crs = crs)
}

.spseg_population <- function(population, variables, data) {
  n <- nrow(data)
  if (!is.null(population)) return(.spseg_mass(population, data, n))
  candidates <- lapply(variables, function(spec) {
    if (!is.null(spec$counts)) rowSums(.spseg_counts(spec, data)) else
      if (!is.null(spec$weights)) .spseg_mass(spec$weights, data, n) else NULL
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) return(rep(1, n))
  first <- candidates[[1L]]
  if (!all(vapply(candidates, function(z) isTRUE(all.equal(z, first,
                    check.attributes = FALSE)), logical(1))))
    stop("Variable populations differ; supply the common neighborhood 'population' explicitly",
         call. = FALSE)
  first
}

.spseg_restore_analysis <- function(result, config, bands, keep_env) {
  for (group in c("multigroup", "pairwise")) {
    result$indices[[group]] <- lapply(result$indices[[group]], function(z) {
      names(z) <- as.character(bands)
      z
    })
  }
  for (key in names(result$indices$pairwise)) {
    result$indices$pairwise[[key]] <- lapply(result$indices$pairwise[[key]], function(z) {
      dimnames(z) <- list(config$pair_labels, config$pair_labels)
      z
    })
  }
  if (length(result$rank_summary)) {
    names(result$rank_summary) <- as.character(bands)
    result$rank_summary <- lapply(result$rank_summary, function(z) {
      dimnames(z) <- list(config$exposure_labels, config$metadata$localenv$labels)
      z
    })
  }
  if (keep_env) {
    names(result$env) <- as.character(bands)
    result$env <- lapply(result$env, function(z) {
      colnames(z) <- config$metadata$localenv$labels
      z
    })
    if (length(result$interpolation)) {
      names(result$interpolation) <- as.character(bands)
      result$interpolation <- lapply(result$interpolation, function(z) {
        colnames(z) <- c("bottom_category_share", "top_category_share",
                         "interpolated_share")
        z
      })
    }
  }
  metadata <- config$metadata
  # Millions of exact values are calculation inputs, not lightweight metadata.
  if (config$type == 3L && !keep_env) {
    metadata$distinct_values <- sum(config$totals > 0)
    metadata$labels <- metadata$scores <- metadata$shares <- NULL
  }
  result$definition <- metadata
  if (!is.null(result$knn)) result$knn <- data.frame(band = bands, result$knn)
  names(result$status) <- as.character(bands)
  result
}

.spseg_compute <- function(x, data, variables, coords, population, measures,
                            bands, weighting, power, normalize, neighbors,
                            search, surface, surface_audit, scope, output,
                            verbose, call, dots, network = NULL, node_id = NULL) {
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
  scope <- match.arg(scope, c("both", "multigroup", "pairwise"))
  output <- match.arg(output, c("indices", "full", "localenv"))
  if (!is.logical(surface_audit) || length(surface_audit) != 1L || is.na(surface_audit))
    stop("'surface_audit' must be TRUE or FALSE", call. = FALSE)
  if (surface == "raw" && length(dots))
    stop("Surface options require surface = 'grid' or 'pycno'", call. = FALSE)
  if (surface != "raw" && any(vapply(variables, function(z) is.null(z$counts), logical(1))))
    stop("Grid and pycnophylactic surfaces require 'counts' for every variable",
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
  supplied_population <- !is.null(population)
  population <- if (neighbors == "knn" || relational) {
    .spseg_mass(population, checked$data, nrow(checked$data))
  } else .spseg_population(population, variables, checked$data)
  if (neighbors == "radius" && (!is.finite(sum(population)) || sum(population) <= 0))
    stop("The common population total must be positive and finite", call. = FALSE)
  compiled <- Map(function(spec, name) {
    .spseg_compile(spec, checked$data, population, measures, variable = name)
  }, variables, names(variables))
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
  engine <- spseg_engine_cpp(checked$coords %||% matrix(numeric(), 0L, 2L),
                         compiled, as.numeric(bands),
                         power, .spseg_weighting_id(weighting), normalize,
                         match(neighbors, c("radius", "knn")) - 1L,
                         match(search, c("boost", "brute")) - 1L,
                         if (relational) prepared_network$graph else list(),
                         as.integer(.spseg_scope_flags(scope)), output != "indices",
                         output != "localenv")
  analyses <- Map(.spseg_restore_analysis, engine$analyses, compiled,
                  MoreArgs = list(bands = bands, keep_env = output != "indices"))
  result <- list(analyses = analyses, bands = as.numeric(bands),
                 coords = if (output != "indices") checked$coords else NULL,
                 data = if (output != "indices") checked$data else NULL,
                 geometry = if (output != "indices") checked$geometry else NULL,
                 population = if (output != "indices") population else NULL,
                 crs = checked$crs, weighting = weighting, power = power,
                 normalize = normalize, surface = surface, output = output,
                 surface_info = surface_info,
                 network = if (relational) prepared_network$summary else NULL,
                 node_id = if (relational && output != "indices") prepared_network$ids else NULL,
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

.spseg_surface_variables <- function(x, checked, variables, compiled, population,
                                    measures, surface, dots, audit, verbose) {
  sizes <- vapply(compiled, function(z) ncol(z$counts), integer(1))
  ends <- cumsum(sizes)
  columns <- Map(seq.int, ends - sizes + 1L, ends)
  combined <- do.call(cbind, lapply(compiled, function(z) z$counts))
  colnames(combined) <- make.unique(unlist(lapply(compiled, function(z) z$metadata$labels)))
  reference <- which(vapply(compiled, function(z) identical(z$mass, population), logical(1)))
  if (!length(reference)) {
    reference_column <- ncol(combined) + 1L
    combined <- cbind(combined, population)
    colnames(combined)[reference_column] <- make.unique(c(colnames(combined)[-reference_column],
                                                         ".population"))[reference_column]
  }
  prepared <- .spseg_surface(x, checked$coords, combined, surface, dots, audit, verbose)
  population <- if (length(reference)) {
    rowSums(prepared$data[, columns[[reference[1L]]], drop = FALSE])
  } else prepared$data[, reference_column]
  recompiled <- Map(function(spec, config, selected) {
    counts <- prepared$data[, selected, drop = FALSE]
    colnames(counts) <- config$metadata$labels
    spec$counts <- seq_len(ncol(counts))
    if (config$type == 3L) spec$scores <- config$scores
    .spseg_compile(spec, as.data.frame(counts), population, measures,
                   resolution = config$summary_config)
  }, variables, compiled, columns)
  checked$coords <- as.matrix(prepared$coords)
  checked$data <- as.data.frame(prepared$data)
  checked$geometry <- prepared$geometry
  prepared$surface_info$variables <- setNames(columns, names(variables))
  list(checked = checked, population = population, compiled = recompiled,
       surface_info = prepared$surface_info)
}
