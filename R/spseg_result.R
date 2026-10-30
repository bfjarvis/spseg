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
#' once when environments are stored. Each analysis records effective `scope`,
#' `output`, `keep_scores`, and `keep_contributions` in `settings`. Global result
#' settings record defaults; variable-specific settings take precedence. Shared
#' inputs are retained once when any analysis stores environments.
#' Each definition records its `eligible_population`. For kNN, each analysis
#' also contains a `knn` table: requested population, minimum/maximum selected
#' population before distance weighting, and minimum/maximum neighborhood
#' radius, by band across focal rows with positive eligible population.
#'
#' @param x A `spseg_result`.
#' @param variables Optional names of analyses to extract.
#' @param what Table to extract: `"indices"` for overall scalar indices,
#'   `"pairwise"` for two-group indices, `"exposure"` for focal-group summaries,
#'   `"env"` for stored local summaries, `"diagnostics"` for per-band status and
#'   kNN diagnostics, or `"interpolation"` for stored neighborhood interpolation
#'   diagnostics. `"scores"` extracts retained local variation/deviation scores;
#'   `"contributions"` extracts their additive index contributions. Both identify
#'   `scope`, `location`, `measure`, and pair members where applicable, and require
#'   their corresponding `keep_scores` or `keep_contributions` option in [spseg()].
#'   Categorical and ordinal pairwise indices compare original categories;
#'   continuous pairwise indices compare resolved global groups.
#'   With `variables = NULL`, mixed-output results include analyses with the
#'   requested stored results. Explicitly selecting an analysis without those
#'   results raises an error; so does requesting results unavailable everywhere.
#' @param bands Optional subset of neighborhood bands.
#' @param row.names,optional Standard data-frame generic arguments.
#' @param ... Reserved for generic methods. Extraction rejects unused arguments.
#' @return Printing returns `x` invisibly. Extraction returns a long data frame
#'   with identifiers and one numeric `value` column for results:
#'   \itemize{
#'     \item Indices: `variable`, `band`, `measure`, `value`.
#'     \item Pairwise indices additionally identify `group_a` and `group_b`.
#'     \item Share exposure: `variable`, `band`, `from_group`, `to_group`, `value`.
#'     \item Quantile exposure: `variable`, `band`, `from_group`, `statistic`,
#'       `prob`, `value`.
#'     \item Local shares: `variable`, `band`, `location`, `category`, `value`.
#'     \item Local ranks: `variable`, `band`, `location`, `statistic`, `prob`, `value`.
#'     \item Scores/contributions: `variable`, `band`, `scope`, `location`,
#'       `measure`, `value`, plus pair members when pairwise results are requested.
#'   }
#'   Columns are chosen from the selected analyses' summary specifications, not
#'   observed missingness. Mixed share/rank tables contain the union of their
#'   required columns; `statistic = "share"` identifies share rows in those tables.
#'   Quantile rows use `statistic = "quantile_rank"`; optional mean ranks use
#'   `"mean_rank"` and have no quantile probability. Undefined results retain
#'   their rows and missing `value`; empty tables retain their typed schema.
#'
#'   With `what = "exposure"`, shares are exposure/isolation probabilities, but
#'   quantiles are focal-population-weighted averages of local global-rank
#'   statistics, not probabilities of encountering a target group. Mean ranks
#'   summarize the full local rank distribution. These rank profiles are stored
#'   separately in `analyses[[name]]$rank_summary`, not in the P matrix. The
#'   summary specification remains in `analyses[[name]]$definition`.
#'
#'   Type, representation, method, normalizers, and calculation settings remain
#'   in the result object rather than being repeated in result tables.
#'   Retained scores and contributions are stored separately in each analysis,
#'   grouped by scope, measure, and band. `st_as_sf()` supports `what = "env"`,
#'   `"scores"`, and `"contributions"` when `output = "full"` or `"localenv"`
#'   has retained geometry or coordinates.
#'   `what = "diagnostics"` returns one row per variable/band with `status` and
#'   any kNN population/radius diagnostics. `what = "interpolation"` returns
#'   one row per stored interval neighborhood/band with `location`,
#'   `bottom_category_share`, `top_category_share`, and `interpolated_share`.
#'   These are weighted local population shares in the first category, last
#'   category, and categories split across requested groups. They are not
#'   uncertainty estimates; allocation fractions remain in
#'   `definition$summary$groups$allocation`. Non-interval analyses contribute
#'   no interpolation rows. Interpolation diagnostics require retained environments.
#'
#'   Network local and interpolation tables also identify `node`. Spatial
#'   extraction returns the compact local table with corresponding geometry.
#'   Network-only results require stored geometry or coordinates for mapping.
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
    cat("    Scope:", a$settings$scope, "| output:", a$settings$output, "\n")
    if (a$definition$type == "continuous") {
      interpretation <- if (a$definition$localenv$type == "quantiles") {
        "focal-group averages of local global-rank statistics, not exposure probabilities"
      } else "exposure/isolation probabilities"
      cat("    Summary:", a$definition$localenv$type, "-", interpretation, "\n")
    }
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
    ..., what = c("indices", "pairwise", "exposure", "env", "diagnostics",
                 "interpolation", "scores", "contributions"), variables = NULL, bands = NULL) {
  .spseg_unused_arguments(list(...))
  what <- match.arg(what)
  extract <- switch(what, indices = .spseg_indices_rows,
    pairwise = .spseg_pairwise_rows, exposure = .spseg_exposure_rows,
    env = .spseg_env_rows, diagnostics = .spseg_diagnostic_rows,
    interpolation = .spseg_interpolation_rows,
    scores = function(...) .spseg_score_rows(..., what = "scores"),
    contributions = function(...) .spseg_score_rows(..., what = "contributions"))
  selected <- .spseg_analysis_subset(x, variables, bands)
  if (is.null(variables)) {
    available <- vapply(x$analyses[selected$variables], function(a) switch(what,
      indices = , pairwise = , exposure = a$settings$output != "localenv",
      env = !is.null(a$env),
      scores = !is.null(a$scores), contributions = !is.null(a$contributions),
      interpolation = !identical(a$definition$representation, "intervals") ||
        !is.null(a$env), TRUE), logical(1))
    if (any(available)) selected$variables <- selected$variables[available]
  }
  analyses <- x$analyses[selected$variables]
  schema <- .spseg_result_schema(analyses, what)
  rows <- list()
  keys <- c(d = "dissimilarity", h = "information", r = "diversity", s = "square_root")
  for (name in selected$variables) {
    analysis <- x$analyses[[name]]
    if (what %in% c("scores", "contributions") && is.null(analysis[[what]])) {
      stop("Local ", what, " were not retained; set keep_", what, " = TRUE",
           call. = FALSE)
    }
    if (what %in% c("indices", "pairwise", "exposure") &&
        analysis$settings$output == "localenv") {
      stop("Indices were not calculated; use output = 'indices' or 'full'", call. = FALSE)
    }
    if (what == "env" && is.null(analysis$env)) {
      stop("Local environments were not stored; use output = 'full' or 'localenv'",
           call. = FALSE)
    }
    if (what == "interpolation" &&
        identical(analysis$definition$representation, "intervals") &&
        is.null(analysis$env)) {
      stop("Interpolation diagnostics were not stored; use output = 'full' or 'localenv'",
           call. = FALSE)
    }
    for (b in selected$bands) {
      base <- list(variable = name, band = x$bands[b])
      rows <- c(rows, extract(analysis, base, b, keys))
    }
  }
  # Fill only structurally inapplicable fields when combining different summaries.
  rows <- lapply(rows, function(block) {
    for (column in setdiff(names(schema), names(block))) {
      block[[column]] <- schema[[column]][rep(NA_integer_, nrow(block))]
    }
    block[names(schema)]
  })
  out <- if (length(rows)) do.call(rbind, rows) else schema
  if (what %in% c("env", "interpolation", "scores", "contributions") && !is.null(x$node_id)) {
    out$node <- x$node_id[out$location]
  }
  rownames(out) <- row.names
  out
}

.spseg_result_schema <- function(analyses, what) {
  base <- list(variable = character(), band = numeric())
  schema <- switch(what,
    indices = c(base, list(measure = character(), value = numeric())),
    pairwise = c(base, list(measure = character(), group_a = character(),
      group_b = character(), value = numeric())),
    exposure = {
      ranks <- vapply(analyses, function(a) {
        a$definition$type == "continuous" && a$definition$localenv$type == "quantiles"
      }, logical(1))
      c(base, list(from_group = character()),
        if (any(!ranks)) list(to_group = character()),
        if (any(ranks)) list(statistic = character(), prob = numeric()),
        list(value = numeric()))
    },
    env = {
      ranks <- vapply(analyses, function(a) a$definition$localenv$type == "quantiles",
                      logical(1))
      c(base, list(location = integer()),
        if (any(!ranks)) list(category = character()),
        if (any(ranks)) list(statistic = character(), prob = numeric()),
        list(value = numeric()))
    },
    diagnostics = {
      knn <- unique(unlist(lapply(analyses, function(a) setdiff(names(a$knn), "band"))))
      c(base, list(status = character()), setNames(rep(list(numeric()), length(knn)), knn))
    },
    interpolation = c(base, list(location = integer(),
      bottom_category_share = numeric(), top_category_share = numeric(),
      interpolated_share = numeric())),
    scores = , contributions = {
      has_pairs <- any(vapply(analyses, function(a)
        length(a[[what]]$pairwise) > 0, logical(1)))
      c(base, list(scope = character(), location = integer(), measure = character()),
        if (has_pairs) list(group_a = character(), group_b = character()),
        list(value = numeric()))
    }
  )
  as.data.frame(schema, stringsAsFactors = FALSE)
}

.spseg_indices_rows <- function(analysis, base, b, keys) {
  lapply(intersect(names(keys), names(analysis$indices$overall)), function(key) {
    data.frame(base, measure = unname(keys[key]),
      value = unname(analysis$indices$overall[[key]][b]), row.names = NULL)
  })
}

.spseg_score_rows <- function(analysis, base, b, keys, what) {
  rows <- list()
  for (scope in c("overall", "pairwise")) {
    store <- analysis[[what]][[scope]]
    for (key in intersect(names(keys), names(store))) {
      values <- store[[key]][[b]]
      if (scope == "overall") {
        rows[[length(rows) + 1L]] <- data.frame(base, scope = scope,
          location = seq_along(values), measure = unname(keys[key]), value = values)
      } else if (ncol(values)) {
        pair <- rep(seq_len(ncol(values)), each = nrow(values))
        rows[[length(rows) + 1L]] <- data.frame(base, scope = scope,
          location = rep(seq_len(nrow(values)), times = ncol(values)),
          measure = unname(keys[key]), group_a = analysis$score_pairs$group_a[pair],
          group_b = analysis$score_pairs$group_b[pair], value = as.vector(values))
      }
    }
  }
  rows
}

.spseg_pairwise_rows <- function(analysis, base, b, keys) {
  rows <- list()
  for (key in intersect(names(keys), names(analysis$indices$pairwise))) {
    mat <- analysis$indices$pairwise[[key]][[b]]
    index <- which(upper.tri(mat), arr.ind = TRUE)
    if (nrow(index)) {
      rows[[length(rows) + 1L]] <- data.frame(base,
        measure = unname(keys[key]), group_a = rownames(mat)[index[, 1]],
        group_b = colnames(mat)[index[, 2]], value = mat[index], row.names = NULL)
    }
  }
  rows
}

.spseg_exposure_rows <- function(analysis, base, b, keys) {
  if (length(analysis$rank_summary)) {
    mat <- analysis$rank_summary[[b]]
    if (!length(mat)) return(list())
    stats <- analysis$definition$localenv$statistics
    column <- rep(seq_len(ncol(mat)), each = nrow(mat))
    return(list(data.frame(base,
      from_group = rep(rownames(mat), times = ncol(mat)),
      statistic = stats$statistic[column], prob = stats$prob[column],
      value = as.vector(mat))))
  }
  mat <- analysis$indices$overall$p[[b]]
  if (is.null(mat) || !length(mat)) return(list())
  list(data.frame(base,
    expand.grid(from_group = rownames(mat), to_group = colnames(mat),
                stringsAsFactors = FALSE), statistic = "share", value = as.vector(mat)))
}

.spseg_env_rows <- function(analysis, base, b, keys) {
  env <- analysis$env[[b]]
  if (!length(env)) return(list())
  definition <- analysis$definition$localenv
  location <- rep(seq_len(nrow(env)), times = ncol(env))
  column <- rep(seq_len(ncol(env)), each = nrow(env))
  if (definition$type == "quantiles") {
    stats <- definition$statistics
    return(list(data.frame(base, location = location,
      statistic = if (is.null(stats)) "quantile_rank" else stats$statistic[column],
      prob = if (is.null(stats)) definition$probs[column] else stats$prob[column],
      value = as.vector(env))))
  }
  list(data.frame(base, location = location, category = colnames(env)[column],
    statistic = "share", value = as.vector(env)))
}

.spseg_diagnostic_rows <- function(analysis, base, b, keys) {
  block <- data.frame(base, status = unname(analysis$status[b]))
  if (!is.null(analysis$knn)) {
    block <- cbind(block, analysis$knn[b, setdiff(names(analysis$knn), "band"), drop = FALSE])
  }
  list(block)
}

.spseg_interpolation_rows <- function(analysis, base, b, keys) {
  if (!length(analysis$interpolation)) return(list())
  mat <- analysis$interpolation[[b]]
  if (!nrow(mat)) return(list())
  list(data.frame(base, location = seq_len(nrow(mat)), mat))
}

#' @rdname spseg_result
#' @method st_as_sf spseg_result
#' @export
st_as_sf.spseg_result <- function(x, ..., variables = NULL, bands = NULL,
                                what = c("env", "scores", "contributions")) {
  .spseg_unused_arguments(list(...))
  if (is.null(x$coords) && is.null(x$geometry))
    stop("No stored geometry or coordinates; network-only results cannot be mapped",
         call. = FALSE)
  what <- match.arg(what)
  out <- as.data.frame(x, what = what, variables = variables, bands = bands)
  if (!nrow(out)) stop("No stored local result rows to extract", call. = FALSE)
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
