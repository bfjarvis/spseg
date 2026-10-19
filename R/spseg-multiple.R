.seg_geometry <- function(x, data, coords) {
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
      if (is.null(data) || ncol(xy) != 2L)
        stop("Data-frame input requires named 'coords', or separate two-column coordinates and 'data'",
             call. = FALSE)
    } else xy <- as.matrix(.seg_select(as.data.frame(x), coords))
  }
  if (!is.numeric(xy) || ncol(xy) != 2L || nrow(xy) < 1L ||
      nrow(xy) != nrow(source) || any(!is.finite(xy)))
    stop("Supply finite two-dimensional coordinates with one location per data row",
         call. = FALSE)
  storage.mode(xy) <- "double"
  list(coords = xy, data = source, geometry = geometry, crs = crs)
}

.seg_population <- function(population, variables, data) {
  n <- nrow(data)
  if (!is.null(population)) return(.seg_mass(population, data, n))
  candidates <- lapply(variables, function(spec) {
    if (!is.null(spec$counts)) rowSums(.seg_counts(spec, data)) else
      if (!is.null(spec$weights)) .seg_mass(spec$weights, data, n) else NULL
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

.seg_restore_analysis <- function(result, config, bands, keep_env) {
  for (group in c("multigroup", "pairwise")) {
    result$indices[[group]] <- lapply(result$indices[[group]], function(z) {
      names(z) <- as.character(bands)
      z
    })
  }
  for (key in names(result$indices$pairwise)) {
    result$indices$pairwise[[key]] <- lapply(result$indices$pairwise[[key]], function(z) {
      dimnames(z) <- list(config$metadata$labels, config$metadata$labels)
      z
    })
  }
  if (keep_env) {
    names(result$env) <- as.character(bands)
    if (config$type != 3L) result$env <- lapply(result$env, function(z) {
      colnames(z) <- config$metadata$labels
      z
    })
  }
  metadata <- config$metadata
  # Millions of exact values are calculation inputs, not lightweight metadata.
  if (config$type == 3L && !keep_env) {
    metadata$distinct_values <- sum(config$totals > 0)
    metadata$labels <- metadata$scores <- metadata$shares <- NULL
  }
  result$definition <- metadata
  names(result$status) <- as.character(bands)
  result
}

.spseg_multiple <- function(x, data, variables, coords, population, measures,
                            bands, weighting, power, normalize, neighbors,
                            search, surface, surface_audit, scope, output,
                            verbose, call, dots) {
  if (!is.list(variables) || !length(variables) || is.null(names(variables)) ||
      anyNA(names(variables)) || any(!nzchar(names(variables))) ||
      anyDuplicated(names(variables)) ||
      !all(vapply(variables, inherits, logical(1), "seg_variable")))
    stop("'variables' must be a uniquely named list of seg_variable specifications",
         call. = FALSE)
  neighbors <- match.arg(neighbors, c("radius", "knn"))
  search <- match.arg(search, c("kdtree", "brute"))
  weighting <- match.arg(weighting, c("biweight", "unweighted", "inverse", "exponential"))
  surface <- match.arg(surface, c("raw", "grid", "pycno"))
  scope <- match.arg(scope, c("both", "multigroup", "pairwise"))
  output <- match.arg(output, c("indices", "full", "localenv"))
  if (surface != "raw" || isTRUE(surface_audit) || length(dots))
    stop("The multiple-variable interface currently supports surface = 'raw' without surface options; legacy surface calls remain available",
         call. = FALSE)
  if (length(power) != 1L || !is.numeric(power) || !is.finite(power) || power <= 0)
    stop("'power' must be a positive finite number", call. = FALSE)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
    stop("'normalize' must be TRUE or FALSE", call. = FALSE)
  checked <- .seg_geometry(x, data, coords)
  population <- .seg_population(population, variables, checked$data)
  if (!is.finite(sum(population)) || sum(population) <= 0)
    stop("The common population total must be positive and finite", call. = FALSE)
  compiled <- lapply(variables, .seg_compile, data = checked$data,
                     population = population, measures = measures)
  if (is.null(bands)) {
    if (neighbors == "knn") stop("Supply population-count 'bands' for kNN", call. = FALSE)
    bands <- default_bands(checked$coords, cbind(population, 0))
  }
  if (!is.numeric(bands) || !length(bands) || any(!is.finite(bands)) ||
      any(bands < 0) || (neighbors == "knn" && any(bands <= 0)))
    stop("Supply finite nonnegative radius bands or positive population-count kNN bands",
         call. = FALSE)
  if (verbose) message("Calculating ", length(compiled), " variables with shared neighborhood searches")
  engine <- seg_multi_cpp(checked$coords, population, compiled, as.numeric(bands),
                         power, spseg_weighting_id(weighting), normalize,
                         match(neighbors, c("radius", "knn")) - 1L,
                         match(search, c("kdtree", "brute")) - 1L,
                         as.integer(spseg_scope_flags(scope)), output != "indices",
                         output != "localenv")
  analyses <- Map(.seg_restore_analysis, engine$analyses, compiled,
                  MoreArgs = list(bands = bands, keep_env = output != "indices"))
  result <- list(analyses = analyses, bands = as.numeric(bands),
                 coords = if (output != "indices") checked$coords else NULL,
                 data = if (output != "indices") checked$data else NULL,
                 geometry = if (output != "indices") checked$geometry else NULL,
                 population = if (output != "indices") population else NULL,
                 crs = checked$crs, weighting = weighting, power = power,
                 normalize = normalize, surface = surface, output = output,
                 neighbors = list(type = neighbors, values = bands,
                                  units = if (neighbors == "knn") "population" else "distance",
                                  engine = search, scope = scope),
                 diagnostics = engine$diagnostics, call = call)
  class(result) <- c("seg_multi_result", "seg_result")
  result
}

.seg_legacy_engine <- function(coords, data, bands, power, weighting, normalize,
                               measures, scope, neighbors, search, output) {
  config <- .seg_compile(seg_categorical(counts = seq_len(ncol(data))),
                         as.data.frame(data), rowSums(data), spseg_measures(measures))
  config$legacy_env <- TRUE
  engine <- seg_multi_cpp(coords, rowSums(data), list(config), bands, power,
                         spseg_weighting_id(weighting), normalize,
                         match(neighbors, c("radius", "knn")) - 1L,
                         match(search, c("kdtree", "brute")) - 1L,
                         as.integer(spseg_scope_flags(scope)), output != "indices", TRUE)
  restored <- .seg_restore_analysis(engine$analyses[[1L]], config, bands, output != "indices")
  if (output != "indices") restored$env <- lapply(restored$env, function(z) {
    rownames(z) <- rownames(data)
    z
  })
  indices <- restored$indices
  for (key in c("d", "r", "h", "p")) if (is.null(indices$multigroup[[key]]))
    indices$multigroup[[key]] <- if (key == "p") list() else numeric()
  list(indices = indices, env = restored$env)
}
