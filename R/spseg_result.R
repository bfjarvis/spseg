`%||%` <- function(x, y) if (is.null(x)) y else x

#' Spatial segregation results
#'
#' [spseg()] returns a `spseg_result`. Each element of `analyses` contains `indices`, optional
#' `env`, `definition`, `normalizers`, and a status for each band. Stored
#' environments are compact matrices: categorical and ordinal shares, or local
#' shares of global groups for continuous variables by default. Continuous
#' `summary` selects synchronized shares or rank summaries; ordinal `localenv`
#' selects quantiles or category shares. Settings are recorded
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
#'   `"exposure"` group summaries, or stored `"env"` distributions. Continuous
#'   pairwise results compare resolved global groups using binary categorical
#'   formulas within the original full-variable neighborhoods.
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
#'   Grouped extraction is explicitly typed: shares have `statistic = "share"`
#'   and retain the `exposure` column. Quantile profiles have `summary = "quantiles"`,
#'   `statistic = "quantile_rank"` or `"mean_rank"`, `rank`, and no target group.
#'   Both have `value`; mean-rank rows have missing `prob`. Rank profiles are stored
#'   in `analyses[[name]]$rank_summary`, separately from the share-valued P matrix.
#'   Local extraction also identifies `statistic`. For interval input, retained
#'   environments include `bottom_category_share`, `top_category_share`, and
#'   `interpolated_share`, stored by band in `analyses[[name]]$interpolation`.
#'   These are shares of the weighted local population in the first category,
#'   last category, and any category apportioned across requested groups. They
#'   are missing for other input types or empty local environments, and repeated
#'   across summary rows for the same location/band. Group allocation fractions
#'   are recorded in `definition$summary$groups$allocation`. A single positive
#'   interval can provide interpolated rank summaries under this stated model.
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
  extract <- switch(what, indices = .spseg_indices_rows,
    pairwise = .spseg_pairwise_rows, exposure = .spseg_exposure_rows,
    env = .spseg_env_rows)
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
      rows <- c(rows, extract(analysis, base, b, keys))
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!length(rows) && what == "env") {
    out <- data.frame(variable = character(), type = character(),
      representation = character(), method = character(), band = numeric(),
      location = integer(), summary = character(), statistic = character(),
      prob = numeric(), rank = numeric(), category = character(), share = numeric(),
      bottom_category_share = numeric(), top_category_share = numeric(),
      interpolated_share = numeric())
  }
  if (!length(rows) && what == "exposure") {
    out <- data.frame(variable = character(), type = character(),
      representation = character(), method = character(), band = numeric(),
      from_group = character(), to_group = character(), summary = character(),
      statistic = character(), prob = numeric(), exposure = numeric(),
      rank = numeric(), value = numeric())
  }
  if (what == "env" && nrow(out) && !is.null(x$node_id))
    out$node <- x$node_id[out$location]
  rownames(out) <- row.names
  out
}

.spseg_indices_rows <- function(analysis, base, b, keys) {
  rows <- list()
  for (key in intersect(names(keys), names(analysis$indices$multigroup))) {
    rows[[length(rows) + 1L]] <- data.frame(base, measure = unname(keys[key]),
      value = analysis$indices$multigroup[[key]][b], status = analysis$status[b],
      row.names = NULL)
  }
  rows
}

.spseg_pairwise_rows <- function(analysis, base, b, keys) {
  rows <- list()
  for (key in intersect(names(keys), names(analysis$indices$pairwise))) {
    mat <- analysis$indices$pairwise[[key]][[b]]
    index <- which(upper.tri(mat), arr.ind = TRUE)
    if (nrow(index)) rows[[length(rows) + 1L]] <- data.frame(base,
      measure = unname(keys[key]), group_a = rownames(mat)[index[, 1]],
      group_b = colnames(mat)[index[, 2]], value = mat[index], row.names = NULL)
  }
  rows
}

.spseg_exposure_rows <- function(analysis, base, b, keys) {
  rows <- list()
  definition <- analysis$definition
  if (length(analysis$rank_summary)) {
    mat <- analysis$rank_summary[[b]]
    stats <- definition$localenv$statistics
    if (length(mat)) {
      column <- rep(seq_len(ncol(mat)), each = nrow(mat))
      rows[[length(rows) + 1L]] <- data.frame(base,
        from_group = rep(rownames(mat), times = ncol(mat)),
        to_group = NA_character_, summary = "quantiles",
        statistic = stats$statistic[column], prob = stats$prob[column],
        exposure = NA_real_, rank = as.vector(mat), value = as.vector(mat))
    }
  } else {
    mat <- analysis$indices$multigroup$p[[b]]
    if (!is.null(mat) && length(mat)) {
      rows[[length(rows) + 1L]] <- data.frame(base,
        expand.grid(from_group = rownames(mat), to_group = colnames(mat),
                    stringsAsFactors = FALSE), summary = "shares",
        statistic = "share", prob = NA_real_, exposure = as.vector(mat),
        rank = NA_real_, value = as.vector(mat))
    }
  }
  rows
}

.spseg_env_rows <- function(analysis, base, b, keys) {
  rows <- list()
  definition <- analysis$definition
  env <- analysis$env[[b]]
  location <- rep(seq_len(nrow(env)), times = ncol(env))
  column <- rep(seq_len(ncol(env)), each = nrow(env))
  quantiles <- definition$localenv$type == "quantiles"
  stats <- definition$localenv$statistics
  if (length(location)) {
    block <- data.frame(base,
    location = location, summary = definition$localenv$type,
    statistic = if (!is.null(stats)) stats$statistic[column] else
      if (quantiles) "quantile_rank" else "share",
    prob = if (!is.null(stats)) stats$prob[column] else
      if (quantiles) definition$localenv$probs[column] else NA_real_,
    rank = if (quantiles) as.vector(env) else NA_real_,
    category = if (quantiles) NA_character_ else colnames(env)[column],
    share = if (quantiles) NA_real_ else as.vector(env))
    for (key in c("bottom_category_share", "top_category_share",
                   "interpolated_share")) {
      block[[key]] <- NA_real_
    }
    if (length(analysis$interpolation)) {
      diagnostic <- as.data.frame(analysis$interpolation[[b]][location, , drop = FALSE])
      block[names(diagnostic)] <- diagnostic
    }
    rows[[length(rows) + 1L]] <- block
  }
  rows
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
