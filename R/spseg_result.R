`%||%` <- function(x, y) if (is.null(x)) y else x

#' Spatial segregation results
#'
#' [spseg()] returns a `spseg_result`. Each element of `analyses` contains `indices`, optional
#' `env`, `definition`, `normalizers`, and a status for each band. Stored
#' environments are compact matrices: categorical shares, or local quantiles
#' of global ranks for ordered variables by default. The ordered constructors'
#' `localenv` option selects quantiles or category shares; settings are recorded
#' in `definition$localenv`. The regional reference distribution is retained
#' once when environments are stored. Inputs and environments are omitted in the
#' default `output = "indices"` mode.
#' Each definition records its `eligible_population`. For kNN, each analysis
#' also contains a `knn` table: requested population, minimum/maximum selected
#' population before distance weighting, and minimum/maximum neighborhood
#' radius, by band across focal rows with positive eligible population.
#'
#' @param x A `spseg_result`.
#' @param variables Optional names of analyses to extract.
#' @param what Extract scalar `"indices"`, categorical `"pairwise"` indices,
#'   `"exposure"` matrices, or stored `"env"` distributions.
#' @param bands Optional subset of neighborhood bands.
#' @param row.names,optional Standard data-frame generic arguments.
#' @param ... Reserved for generic methods. Extraction rejects unused arguments.
#' @return Printing returns `x` invisibly. Data-frame extraction returns a long
#'   table with the analysis name, type, continuous input `representation`
#'   (`NA` for other types), method, and band. Environment tables
#'   identify the input row as `location`. Spatial extraction returns this
#'   environment table with the corresponding input geometry attached. Quantile
#'   rows have `summary = "quantiles"`, `prob`, and `rank`; share rows have
#'   `summary = "shares"`, `category`, and `share`. Unused columns are `NA`.
#'   Empty environments retain rows with missing summaries.
#'   Network environments also identify `node`. Network-only results support
#'   ordinary extraction; spatial extraction requires stored geometry or coordinates.
#' @name spseg_result
NULL

#' @rdname spseg_result
#' @method print spseg_result
#' @export
print.spseg_result <- function(x, ...) {
  cat("\nSpatial segregation across", length(x$analyses), "variables\n")
  cat("Neighborhoods:", x$neighbors$type, "| bands:", paste(x$bands, collapse = ", "), "\n")
  cat("Space:", x$neighbors$space, "| search:", x$neighbors$engine, "\n")
  for (name in names(x$analyses)) {
    a <- x$analyses[[name]]
    cat(" ", name, ": ", a$definition$type, " (", a$definition$method, ")\n", sep = "")
    bad <- unique(a$status[a$status != "ok"])
    if (length(bad)) cat("    Status:", paste(bad, collapse = ", "), "\n")
  }
  invisible(x)
}

.spseg_analysis_subset <- function(x, variables, bands) {
  if (is.null(variables)) variables <- names(x$analyses)
  if (!is.character(variables) || !length(variables) || anyNA(variables) ||
      !all(variables %in% names(x$analyses)))
    stop("'variables' must name stored analyses", call. = FALSE)
  b <- if (is.null(bands)) seq_along(x$bands) else match(bands, x$bands)
  if (!length(b) || anyNA(b)) stop("'bands' must match stored bands", call. = FALSE)
  list(variables = variables, bands = b)
}

#' @rdname spseg_result
#' @method as.data.frame spseg_result
#' @export
as.data.frame.spseg_result <- function(x, row.names = NULL, optional = FALSE,
    ..., what = c("indices", "pairwise", "exposure", "env"), variables = NULL,
    bands = NULL) {
  .spseg_unused_arguments(list(...))
  what <- match.arg(what)
  selected <- .spseg_analysis_subset(x, variables, bands)
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
      base <- list(variable = name, type = definition$type,
                   representation = definition$representation %||% NA_character_,
                   method = definition$method,
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
      } else {
        env <- analysis$env[[b]]
        location <- rep(seq_len(nrow(env)), times = ncol(env))
        column <- rep(seq_len(ncol(env)), each = nrow(env))
        quantiles <- definition$localenv$type == "quantiles"
        if (length(location)) rows[[length(rows) + 1L]] <- data.frame(base,
          location = location, summary = definition$localenv$type,
          prob = if (quantiles) definition$localenv$probs[column] else NA_real_,
          rank = if (quantiles) as.vector(env) else NA_real_,
          category = if (quantiles) NA_character_ else colnames(env)[column],
          share = if (quantiles) NA_real_ else as.vector(env))
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (what == "env" && nrow(out) && !is.null(x$node_id))
    out$node <- x$node_id[out$location]
  rownames(out) <- row.names
  out
}

#' @rdname spseg_result
#' @method st_as_sf spseg_result
#' @export
st_as_sf.spseg_result <- function(x, ..., variables = NULL, bands = NULL) {
  .spseg_unused_arguments(list(...))
  if (is.null(x$coords) && is.null(x$geometry))
    stop("No stored geometry or coordinates; network-only results cannot be mapped",
         call. = FALSE)
  out <- as.data.frame(x, what = "env", variables = variables, bands = bands)
  if (!nrow(out)) stop("No stored environment rows to extract", call. = FALSE)
  if (!is.null(x$geometry)) sf::st_sf(out, geometry = x$geometry[out$location]) else {
    out$x <- x$coords[out$location, 1]
    out$y <- x$coords[out$location, 2]
    sf::st_as_sf(out, coords = c("x", "y"), crs = x$crs)
  }
}

#' @rdname spseg_result
#' @method as.list spseg_result
#' @export
as.list.spseg_result <- function(x, ...) unclass(x)

.spseg_unused_arguments <- function(dots) {
  if (length(dots)) stop("Unused extraction arguments: ",
                         paste(names(dots), collapse = ", "), call. = FALSE)
}
