`%||%` <- function(x, y) if (is.null(x)) y else x

#' Spatial segregation results
#'
#' [spseg()] returns a `seg_result`. Each element of `analyses` contains `indices`, optional
#' `env`, `definition`, `normalizers`, and a status for each band. Stored
#' categorical/ordinal/grouped environments are category-share matrices;
#' continuous environments are sparse lists of exact values, global midranks,
#' and local population shares. Inputs and environments are omitted in the
#' default `output = "indices"` mode.
#' Each definition records its `eligible_population`. For kNN, each analysis
#' also contains a `knn` table: requested population, minimum/maximum selected
#' population before distance weighting, and minimum/maximum neighborhood
#' radius, by band across focal rows with positive eligible population.
#'
#' @param x A `seg_result`.
#' @param variables Optional names of analyses to extract.
#' @param what Extract scalar `"indices"`, categorical `"pairwise"` indices,
#'   `"exposure"` matrices, or stored `"env"` distributions.
#' @param bands Optional subset of neighborhood bands.
#' @param row.names,optional Standard data-frame generic arguments.
#' @param ... Reserved for generic methods. Extraction rejects unused arguments.
#' @return Printing returns `x` invisibly. Data-frame extraction returns a long
#'   table with the analysis name, type, method, and band. Environment tables
#'   identify the input row as `location`. Spatial extraction returns this
#'   environment table with the corresponding input geometry attached.
#' @name seg_result
NULL

#' @rdname seg_result
#' @export
print.seg_result <- function(x, ...) {
  cat("\nSpatial segregation across", length(x$analyses), "variables\n")
  cat("Neighborhoods:", x$neighbors$type, "| bands:", paste(x$bands, collapse = ", "), "\n")
  for (name in names(x$analyses)) {
    a <- x$analyses[[name]]
    cat(" ", name, ": ", a$definition$type, " (", a$definition$method, ")\n", sep = "")
    bad <- unique(a$status[a$status != "ok"])
    if (length(bad)) cat("    Status:", paste(bad, collapse = ", "), "\n")
  }
  invisible(x)
}

.seg_analysis_subset <- function(x, variables, bands) {
  if (is.null(variables)) variables <- names(x$analyses)
  if (!is.character(variables) || !length(variables) || anyNA(variables) ||
      !all(variables %in% names(x$analyses)))
    stop("'variables' must name stored analyses", call. = FALSE)
  b <- if (is.null(bands)) seq_along(x$bands) else match(bands, x$bands)
  if (!length(b) || anyNA(b)) stop("'bands' must match stored bands", call. = FALSE)
  list(variables = variables, bands = b)
}

#' @rdname seg_result
#' @export
as.data.frame.seg_result <- function(x, row.names = NULL, optional = FALSE,
    ..., what = c("indices", "pairwise", "exposure", "env"), variables = NULL,
    bands = NULL) {
  .seg_unused_arguments(list(...))
  what <- match.arg(what)
  selected <- .seg_analysis_subset(x, variables, bands)
  rows <- list()
  keys <- c(d = "dissimilarity", h = "information", r = "diversity", s = "square_root")
  for (name in selected$variables) {
    analysis <- x$analyses[[name]]
    definition <- analysis$definition
    if (what != "env" && x$output == "localenv")
      stop("Indices were not calculated; use output = 'indices' or 'full'", call. = FALSE)
    if (what == "env" && is.null(analysis$env))
      stop("Local environments were not stored; use output = 'full' or 'localenv'", call. = FALSE)
    for (b in selected$bands) {
      base <- list(variable = name, type = definition$type, method = definition$method,
                   band = x$bands[b])
      if (what == "indices") {
        for (key in intersect(names(keys), names(analysis$indices$multigroup))) {
          rows[[length(rows) + 1L]] <- data.frame(base, measure = unname(keys[key]),
            value = analysis$indices$multigroup[[key]][b], status = analysis$status[b],
            row.names = NULL)
        }
      } else if (what == "pairwise") {
        for (key in intersect(names(keys), names(analysis$indices$pairwise))) {
          mat <- analysis$indices$pairwise[[key]][[b]]
          index <- which(upper.tri(mat), arr.ind = TRUE)
          if (nrow(index)) rows[[length(rows) + 1L]] <- data.frame(base,
            measure = unname(keys[key]), group_a = rownames(mat)[index[, 1]],
            group_b = colnames(mat)[index[, 2]], value = mat[index], row.names = NULL)
        }
      } else if (what == "exposure") {
        mat <- analysis$indices$multigroup$p[[b]]
        if (!is.null(mat) && length(mat)) rows[[length(rows) + 1L]] <- data.frame(base,
          expand.grid(from_group = rownames(mat), to_group = colnames(mat),
                      stringsAsFactors = FALSE), exposure = as.vector(mat))
      } else if (definition$type == "continuous") {
        env <- analysis$env[[b]]
        for (i in seq_along(env)) if (length(env[[i]]$share)) {
          rows[[length(rows) + 1L]] <- data.frame(base, location = i,
            category = NA_character_, value = env[[i]]$value,
            rank = env[[i]]$rank, share = env[[i]]$share)
        }
      } else {
        env <- analysis$env[[b]]
        index <- which(!is.na(env), arr.ind = TRUE)
        if (nrow(index)) rows[[length(rows) + 1L]] <- data.frame(base,
          location = index[, 1], category = colnames(env)[index[, 2]],
          value = NA_real_, rank = NA_real_, share = env[index])
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  rownames(out) <- row.names
  out
}

#' @rdname seg_result
#' @export
st_as_sf.seg_result <- function(x, ..., variables = NULL, bands = NULL) {
  .seg_unused_arguments(list(...))
  if (is.null(x$coords)) stop("Coordinates were not stored; use output = 'full' or 'localenv'",
                              call. = FALSE)
  out <- as.data.frame(x, what = "env", variables = variables, bands = bands)
  if (!nrow(out)) stop("No stored environment rows to extract", call. = FALSE)
  if (!is.null(x$geometry)) sf::st_sf(out, geometry = x$geometry[out$location]) else {
    out$x <- x$coords[out$location, 1]
    out$y <- x$coords[out$location, 2]
    sf::st_as_sf(out, coords = c("x", "y"), crs = x$crs)
  }
}

#' @rdname seg_result
#' @export
as.list.seg_result <- function(x, ...) unclass(x)

.seg_unused_arguments <- function(dots) {
  if (length(dots)) stop("Unused extraction arguments: ",
                         paste(names(dots), collapse = ", "), call. = FALSE)
}
