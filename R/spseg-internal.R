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

.spseg_population <- function(population, variables, data, count_data) {
  n <- nrow(data)
  if (!is.null(population)) return(.spseg_mass(population, data, n))
  candidates <- Map(function(spec, counts) {
    if (!is.null(counts)) rowSums(counts) else
      if (!is.null(spec$weights)) .spseg_mass(spec$weights, data, n) else NULL
  }, variables, count_data)
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) return(rep(1, n))
  first <- candidates[[1L]]
  if (!all(vapply(candidates, function(z) isTRUE(all.equal(z, first,
                    check.attributes = FALSE)), logical(1))))
    stop("Variable populations differ; supply the common neighborhood 'population' explicitly",
         call. = FALSE)
  first
}

.spseg_prepare_result <- function(result, config, bands, keep_env) {
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
