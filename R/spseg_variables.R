#' Describe variables for a shared spatial segregation calculation
#'
#' These constructors describe data columns; they do not search neighborhoods.
#' Supply exactly one of `values` (a single column) and `counts` (columns in
#' category/value order). Column selectors are names or integer positions.
#' Continuous variables accept exact observations or interval membership/counts
#' through `representation`. Interval CDFs are interpolated directly; individual
#' values are not imputed.
#'
#' @param values Column containing category membership or observed values.
#' @param counts Columns containing nonnegative population counts.
#' @param levels Category order for observed ordinal or interval membership.
#'   Ordered factors supply their own order. For categorical values, factor
#'   levels or sorted distinct values are used when this is omitted.
#' @param measures Optional measure selection overriding [spseg()].
#' @param weights Population-weight column or numeric vector for `values`.
#'   Defaults to the `population` supplied to [spseg()], or its inferred default.
#' @param na.action Missing-value policy: `"fail"` or `"omit"`. Omission applies
#'   only to this variable, including its kNN stopping population. Other
#'   variables and radius bandwidths are unaffected.
#' @param representation Continuous input representation: `"exact"` (default)
#'   for observed numeric values or exact-value counts, or `"intervals"` for
#'   interval membership or ordered interval counts. Membership requires an
#'   ordered factor or explicit `levels`; count columns supply interval order.
#' @param interpolation For `representation = "intervals"`, CDF interpolation.
#'   `NULL` selects `"linear"` in overall percentile space, not in the original
#'   measurement units. Exact input does not accept this option.
#' @param scores Exact numerical values corresponding to continuous count
#'   columns with `representation = "exact"`. These are observed values, not
#'   interval midpoints. Interval input does not accept scores.
#' @param summary Continuous summary specification, defaulting to
#'   `list(type = "shares", probs = c(.25, .5, .75))`. Shares use the same
#'   global groups for retained local shares and group-averaged exposure.
#'   `type = "quantiles"` returns local global-rank quantiles and their
#'   focal-group averages instead; `mean_rank = TRUE` additionally reports the
#'   mean of the full local rank distribution. `probs` is one strictly increasing
#'   vector in [0, 1]; optional 0/1 delimiters define tails but are not evaluated
#'   as local extrema. Exact shares may use finite absolute `cuts` instead of
#'   `probs`, with implicit infinite tails. Cuts require exact values and shares
#'   mode. Tied values remain together, so achieved group sizes may differ.
#'   Interval input retains the requested grid. Categories crossing requested
#'   boundaries are proportionally apportioned in global-percentile space, using
#'   the same fractions for original focal populations and local target shares.
#'   This assumes the same within-category rank distribution at every location;
#'   it does not require monetary tail bounds. Allocation is resolved once before
#'   surface construction and recorded in `definition$summary$groups$allocation`.
#'   One warning identifies when interpolation is required. Retained neighborhoods
#'   report `bottom_category_share`, `top_category_share`, and `interpolated_share`
#'   rather than printing numerical details in the warning. The last is the local
#'   population share in input categories split across requested groups, not an
#'   uncertainty estimate. Choose probabilities on actual cumulative interval
#'   boundaries to avoid apportionment. Summary choices do not change continuous
#'   scalar indices.
#' @param exposure Deprecated for continuous variables: use `summary` instead.
#'   Explicit legacy calls retain the previous behavior with a warning and cannot
#'   be combined with `summary`. For exact continuous variables, a list with either
#'   `probs` (global weighted quantile boundaries, including 0 and 1) or `breaks`
#'   (value boundaries covering the observations). Defaults to deciles.
#'   Equal values stay together; empty groups are removed. These cuts do not
#'   affect the rank-order indices. They also supply the default groups for
#'   `localenv = "shares"` when no local cuts are specified.
#'   Interval input uses the supplied intervals for exposure and does not accept
#'   this option; custom local presentation cuts belong in `localenv`.
#' @param localenv Deprecated for continuous variables in favor of `summary`;
#'   explicit legacy calls retain independent local/exposure summaries. For
#'   ordinal variables, the current retained summary defaults to `"shares"`.
#'   Only the deprecated explicit continuous path defaults to `"quantiles"`
#'   (including when only legacy `exposure` is supplied). Select either summary,
#'   or a list with `type` and optional `probs`/`breaks`. Legacy continuous and
#'   ordinal quantiles default to `seq(.1, .9, .1)` and report local
#'   weighted quantiles of global population ranks. Ordinal and exact continuous
#'   observations use weighted midranks and the inverse weighted empirical CDF;
#'   grouped intervals use the inverse linearly interpolated CDF. Probabilities
#'   must increase strictly within [0, 1]; endpoints select the local support
#'   limits. Shares default to the supplied ordinal categories or intervals,
#'   or the continuous variable's exposure groups (deciles unless overridden).
#'   For shares, `probs` must increase from 0 to 1 and define global quantile
#'   groups; grouped intervals are apportioned using their interpolated CDF.
#'   Only exact continuous variables accept original-value `breaks`.
#'   These settings affect retained environments only, not indices or exposure.
#'   Environments are stored by `output = "full"` or `"localenv"` in [spseg()].
#' @return A `spseg_variable` specification for `spseg(variables = ...)`.
#' @seealso [spseg()]
#' @examples
#' spseg_categorical(values = "ethnicity")
#' spseg_ordinal(counts = c("primary", "secondary", "tertiary"))
#' spseg_continuous(counts = c("low", "middle", "high"), representation = "intervals")
#' spseg_continuous(values = "income", summary = list(probs = c(.2, .8)))
#' spseg_continuous(values = "income",
#'                  summary = list(type = "quantiles", mean_rank = TRUE))
#' spseg_ordinal(values = "education", localenv = "quantiles")
#' @name spseg_variables
NULL

.spseg_variable <- function(type, values, counts, levels, measures, weights,
                          na.action, interpolation = NULL, scores = NULL,
                          exposure = NULL, localenv = "shares", representation = NULL) {
  if (is.null(values) == is.null(counts)) {
    stop("Supply exactly one of 'values' and 'counts'", call. = FALSE)
  }
  if (!is.null(values) && length(values) != 1L) {
    stop("'values' must select one column", call. = FALSE)
  }
  if (!is.null(counts) && length(counts) == 0L) {
    stop("'counts' must select at least one column", call. = FALSE)
  }
  if (!is.null(counts) && !is.null(weights)) {
    stop("'weights' applies to 'values'; count columns already contain mass",
         call. = FALSE)
  }
  structure(list(type = type, values = values, counts = counts,
                 levels = levels, measures = measures, weights = weights,
                 na.action = match.arg(na.action, c("fail", "omit")),
                 interpolation = interpolation, scores = scores,
                 exposure = exposure, representation = representation,
                 localenv = .spseg_localenv_options(localenv)),
            class = "spseg_variable")
}

#' @rdname spseg_variables
#' @export
spseg_categorical <- function(values = NULL, counts = NULL, levels = NULL,
                            measures = NULL, weights = NULL,
                            na.action = c("fail", "omit")) {
  .spseg_variable("categorical", values, counts, levels, measures, weights,
                match.arg(na.action))
}

#' @rdname spseg_variables
#' @export
spseg_ordinal <- function(values = NULL, counts = NULL, levels = NULL,
                       measures = NULL, weights = NULL,
                       na.action = c("fail", "omit"), localenv = "shares") {
  .spseg_variable("ordinal", values, counts, levels, measures, weights,
                match.arg(na.action), localenv = localenv)
}

#' @rdname spseg_variables
#' @export
spseg_continuous <- function(values = NULL, counts = NULL,
                          representation = c("exact", "intervals"),
                          levels = NULL, scores = NULL,
                          measures = NULL, weights = NULL,
                          na.action = c("fail", "omit"), exposure = NULL,
                          localenv = "quantiles", interpolation = NULL,
                          summary = NULL) {
  legacy <- !missing(exposure) || !missing(localenv)
  if (legacy && !missing(summary)) {
    stop("Do not combine 'summary' with deprecated 'localenv' or 'exposure'",
         call. = FALSE)
  }
  representation <- match.arg(representation)
  if (representation == "intervals") {
    if (!is.null(scores))
      stop("'scores' applies only to representation = 'exact'", call. = FALSE)
    if (!is.null(exposure))
      stop("Interval exposure uses the supplied intervals; customize 'localenv' for presentation",
           call. = FALSE)
    interpolation <- interpolation %||% "linear"
    if (!identical(interpolation, "linear"))
      stop("Only 'linear' interpolation is supported", call. = FALSE)
  } else {
    if (!is.null(levels))
      stop("Continuous 'levels' applies only to representation = 'intervals'", call. = FALSE)
    if (!is.null(interpolation))
      stop("'interpolation' applies only to representation = 'intervals'", call. = FALSE)
  }
  if (!is.null(values) && !is.null(scores)) {
    stop("'scores' applies only to continuous count columns", call. = FALSE)
  }
  out <- .spseg_variable("continuous", values, counts, levels, measures, weights,
                match.arg(na.action), scores = scores, exposure = exposure,
                localenv = localenv, representation = representation,
                interpolation = interpolation)
  out$summary <- if (legacy) NULL else .spseg_summary_options(summary)
  if (legacy) {
    warning("Continuous 'localenv' and 'exposure' are deprecated; use 'summary'. ",
            "Explicit legacy settings retain their previous behavior.", call. = FALSE)
  }
  out
}

.spseg_localenv_options <- function(options) {
  if (is.character(options) && length(options) == 1L && !is.na(options))
    options <- list(type = options)
  if (!is.list(options) || is.null(names(options)) || anyDuplicated(names(options)) ||
      any(!names(options) %in% c("type", "probs", "breaks")) ||
      !is.character(options$type) || length(options$type) != 1L ||
      is.na(options$type) || !options$type %in% c("quantiles", "shares"))
    stop("'localenv' must select 'quantiles' or 'shares', optionally with probs or breaks",
         call. = FALSE)
  if (!is.null(options$probs) && !is.null(options$breaks))
    stop("Supply only one of localenv 'probs' and 'breaks'", call. = FALSE)
  if (options$type == "quantiles") {
    if (!is.null(options$breaks))
      stop("Local quantiles accept 'probs', not 'breaks'", call. = FALSE)
    options$probs <- options$probs %||% seq(.1, .9, .1)
  }
  p <- options$probs
  if (!is.null(p) && (!is.numeric(p) || !length(p) || any(!is.finite(p)) ||
      any(p < 0 | p > 1) || is.unsorted(p, strictly = TRUE)))
    stop("Local-environment probabilities must increase strictly within [0, 1]",
         call. = FALSE)
  if (options$type == "shares" && !is.null(p) &&
      (length(p) < 2L || p[1L] != 0 || p[length(p)] != 1))
    stop("Local share probabilities must increase from 0 to 1", call. = FALSE)
  b <- options$breaks
  if (!is.null(b) && (!is.numeric(b) || length(b) < 2L || anyNA(b) ||
      is.unsorted(b, strictly = TRUE)))
    stop("Local share breaks must increase strictly", call. = FALSE)
  options
}

.spseg_localenv <- function(spec, labels, totals, scores) {
  options <- spec$localenv
  intervals <- identical(spec$representation, "intervals")
  exact <- spec$type == "continuous" && !intervals
  if (!is.null(options$breaks) && !exact)
    stop("Local share 'breaks' require exact continuous values; use 'probs' for intervals",
         call. = FALSE)
  if (options$type == "quantiles") {
    return(list(type = 1L, probs = as.numeric(options$probs), breaks = numeric(),
      map = integer(), definition = list(type = "quantiles", probs = options$probs,
        labels = paste0("q", format(options$probs, trim = TRUE)), scale = "global_rank",
        method = if (intervals) "inverse_linear_cdf" else "inverse_midrank_cdf")))
  }
  if (intervals && !is.null(options$probs)) {
    p <- options$probs
    return(list(type = 2L, probs = numeric(), breaks = as.numeric(p), map = integer(),
      definition = list(type = "shares", probs = p,
        labels = paste0("(", p[-length(p)], ",", p[-1L], "]"),
        scale = "global_rank", method = "linear_cdf")))
  }
  cuts <- options[intersect(names(options), c("probs", "breaks"))]
  if (exact || length(cuts)) {
    if (!length(cuts)) cuts <- spec$exposure
    groups <- .spseg_exposure(if (exact) scores else seq_along(labels),
                              totals, cuts)
    definition <- c(list(type = "shares"), groups$definition)
    map <- groups$map
  } else {
    definition <- list(type = "shares", labels = labels, method = "supplied_categories")
    map <- seq_along(labels) - 1L
  }
  list(type = 0L, probs = numeric(), breaks = numeric(), map = as.integer(map),
       definition = definition)
}

.spseg_select <- function(data, columns) {
  if (is.character(columns)) {
    if (anyNA(columns) || !all(columns %in% names(data)))
      stop("Unknown column: ", paste(setdiff(columns, names(data)), collapse = ", "),
           call. = FALSE)
    index <- match(columns, names(data))
  } else if (is.numeric(columns) && !anyNA(columns) &&
             all(is.finite(columns)) && all(columns >= 1 & columns <= ncol(data)) &&
             all(columns == as.integer(columns))) {
    index <- as.integer(columns)
  } else stop("Column selectors must be names or valid integer positions", call. = FALSE)
  if (anyDuplicated(index)) stop("Columns must not be repeated", call. = FALSE)
  data[, index, drop = FALSE]
}

.spseg_mass <- function(value, data, n, default = rep(1, n)) {
  if (is.null(value)) return(default)
  if (is.character(value) && length(value) == 1L)
    value <- .spseg_select(data, value)[[1L]]
  if (!is.numeric(value) || length(value) != n ||
      any(!is.finite(value) | value < 0)) {
    stop("Population weights must be finite, nonnegative, and have one value per row",
         call. = FALSE)
  }
  as.numeric(value)
}

.spseg_counts <- function(spec, data) {
  out <- .spseg_select(data, spec$counts)
  if (!all(vapply(out, is.numeric, logical(1))))
    stop("Count columns must be numeric", call. = FALSE)
  out <- as.matrix(out)
  bad <- !stats::complete.cases(out)
  if (any(bad)) {
    if (spec$na.action == "fail") stop("Missing counts; choose na.action = 'omit' explicitly",
                                       call. = FALSE)
    out[bad, ] <- 0
  }
  if (any(!is.finite(out) | out < 0))
    stop("Counts must be finite and nonnegative", call. = FALSE)
  storage.mode(out) <- "double"
  out
}

.spseg_measures <- function(measures) {
  measures[measures %in% c("relative_diversity", "variance_ratio")] <- "diversity"
  valid <- c("exposure", "information", "diversity", "dissimilarity", "square_root")
  if (!is.character(measures) || !length(measures) || anyNA(measures) ||
      any(!measures %in% c(valid, "all")))
    stop("Unknown or empty measure selection", call. = FALSE)
  if ("all" %in% measures) valid else unique(measures)
}

.spseg_exposure <- function(scores, totals, options) {
  if (is.null(options)) options <- list(probs = seq(0, 1, .1))
  if (!is.list(options) || length(options) != 1L || is.null(names(options)) ||
      !all(names(options) %in% c("probs", "breaks")))
    stop("'exposure' must contain exactly one of 'probs' or 'breaks'", call. = FALSE)
  present <- totals > 0
  v <- scores[present]
  mass <- totals[present]
  if (!length(v)) return(list(map = rep(-1L, length(scores)), labels = character(),
                              definition = list(requested = options)))
  if (!is.null(options$probs)) {
    p <- options$probs
    if (!is.numeric(p) || length(p) < 2L || any(!is.finite(p)) ||
        is.unsorted(p, strictly = TRUE) || p[1L] != 0 || p[length(p)] != 1)
      stop("Exposure probabilities must increase strictly from 0 to 1", call. = FALSE)
    cumulative <- cumsum(mass) / sum(mass)
    inner <- v[vapply(p[-c(1L, length(p))],
                       function(q) which(cumulative >= q - 8*.Machine$double.eps)[1L], integer(1))]
    breaks <- c(-Inf, unique(inner[inner < max(v)]), Inf)
  } else {
    breaks <- options$breaks
    if (!is.numeric(breaks) || length(breaks) < 2L || anyNA(breaks) ||
        is.unsorted(breaks, strictly = TRUE) ||
        breaks[1L] > min(v) || breaks[length(breaks)] < max(v))
      stop("Exposure breaks must increase strictly and cover the observed values",
           call. = FALSE)
  }
  groups <- cut(scores, breaks, include.lowest = TRUE, right = TRUE)
  original <- as.integer(groups)
  shares <- vapply(seq_len(length(breaks) - 1L),
                   function(k) sum(totals[which(original == k)]), numeric(1))
  keep <- which(shares > 0)
  list(map = match(original, keep) - 1L, labels = levels(groups)[keep],
       definition = list(requested = options, breaks = breaks, right = TRUE,
                         labels = levels(groups)[keep], shares = shares[keep] / sum(shares)))
}

.spseg_compile <- function(spec, data, population, measures, resolution = NULL,
                           variable = NULL) {
  n <- nrow(data)
  intervals <- identical(spec$representation, "intervals")
  exact <- spec$type == "continuous" && !intervals
  selected <- .spseg_measures(spec$measures %||% measures)
  counts <- matrix(numeric(), 0L, 0L)
  code <- rep(-1L, n)
  scores <- numeric()
  if (!is.null(spec$counts)) {
    counts <- .spseg_counts(spec, data)
    labels <- colnames(counts)
    mass <- rowSums(counts)
    totals <- colSums(counts)
    if (exact) {
      scores <- spec$scores
      if (!is.numeric(scores) || length(scores) != ncol(counts) ||
          any(!is.finite(scores)))
        stop("Supply finite 'scores' matching the continuous count columns", call. = FALSE)
      if (anyDuplicated(scores)) {
        # Collapse duplicate exact values, retaining the population mass.
        keys <- sort(unique(scores))
        counts <- vapply(keys, function(z) rowSums(counts[, scores == z, drop = FALSE]),
                         numeric(n))
        dim(counts) <- c(n, length(keys))
        scores <- keys
      } else {
        order <- order(scores)
        scores <- scores[order]
        counts <- counts[, order, drop = FALSE]
      }
      totals <- colSums(counts)
      labels <- as.character(scores)
    }
  } else {
    values <- .spseg_select(data, spec$values)[[1L]]
    mass <- .spseg_mass(spec$weights, data, n, population)
    bad <- is.na(values)
    if (any(bad) && spec$na.action == "fail")
      stop("Missing observed values; choose na.action = 'omit' explicitly", call. = FALSE)
    mass[bad] <- 0
    if (exact) {
      if (!is.numeric(values) || any(!is.finite(values[!bad])))
        stop("Continuous values must be finite numeric observations", call. = FALSE)
      scores <- sort(unique(values[mass > 0]))
      labels <- as.character(scores)
      code <- match(values, scores) - 1L
    } else {
      levels <- spec$levels
      if (is.null(levels)) {
        if (spec$type != "categorical" && !is.ordered(values))
          stop("Supply ordered 'levels' or an ordered factor", call. = FALSE)
        levels <- if (is.factor(values)) levels(values) else sort(unique(values[!bad]))
      }
      if (anyNA(levels) || anyDuplicated(levels) || !length(levels))
        stop("Category levels must be nonempty, distinct, and nonmissing", call. = FALSE)
      labels <- as.character(levels)
      code <- match(values, levels) - 1L
      if (any(is.na(code) & !bad))
        stop("Some observed categories are absent from 'levels'", call. = FALSE)
    }
    code[is.na(code)] <- -1L
    totals <- numeric(length(labels))
    ok <- code >= 0L & mass > 0
    if (any(ok)) {
      sums <- rowsum(mass[ok], code[ok] + 1L, reorder = FALSE)
      totals[as.integer(rownames(sums))] <- sums[, 1L]
    }
  }
  total <- sum(mass)
  if (!is.finite(total) || total <= 0)
    stop("Each variable must have a positive finite eligible population", call. = FALSE)
  shares <- totals / total
  type <- switch(spec$type, categorical = 0L, ordinal = 1L,
                 continuous = if (intervals) 2L else 3L)
  method <- c("categorical", "ordinal_cutoffs", "linear_cdf", "empirical_midrank")[type + 1L]
  synchronized <- !is.null(spec$summary)
  resolved <- NULL
  exposure <- if (synchronized) {
    resolved <- resolution %||% .spseg_resolve_summary(spec$summary, labels, totals,
                                                      scores, intervals, variable)
    list(map = resolved$map, labels = resolved$labels,
         definition = resolved$definition$groups)
  } else if (type == 3L) {
    if ("exposure" %in% selected) .spseg_exposure(scores, totals, spec$exposure) else
      list(map = integer(), labels = character(), definition = NULL)
  } else list(map = seq_along(labels) - 1L, labels = labels, definition = NULL)
  localenv <- if (synchronized) {
    definition <- resolved$definition
    list(type = if (definition$type == "shares") 0L else 1L,
      probs = definition$effective_probs %||% numeric(), breaks = numeric(),
      map = resolved$map,
      definition = list(type = definition$type, probs = definition$effective_probs,
        labels = definition$statistics$label, statistics = definition$statistics,
        scale = if (definition$type == "quantiles") "global_rank" else "group_share"))
  } else .spseg_localenv(spec, labels, totals, scores)
  list(type = type, counts = counts, code = as.integer(code), mass = mass,
       synchronized = synchronized, summary_config = resolved,
       group_fractions = if (synchronized) resolved$fractions else matrix(numeric(), 0, 0),
       split_intervals = if (synchronized) {
         as.integer(resolved$definition$interpolation$split_intervals - 1L)
       } else integer(),
       env_mean = synchronized && spec$summary$mean_rank,
       pair_labels = if (synchronized) resolved$labels else labels,
       totals = totals, shares = shares, scores = scores,
       z = cumsum(shares) - shares / 2,
       exposure_map = as.integer(exposure$map), exposure_labels = exposure$labels,
       env_type = localenv$type, env_probs = localenv$probs, env_breaks = localenv$breaks,
       env_map = localenv$map, env_columns = length(localenv$definition$labels),
       flags = as.integer(c("exposure", "information", "diversity", "dissimilarity",
                            "square_root") %in% selected),
       metadata = list(type = spec$type, representation = spec$representation,
                       interpolation = spec$interpolation,
                       method = method, measures = selected,
                       eligible_population = total,
                       labels = labels, shares = shares, scores = scores,
                       ties = if (type == 3L) "weighted_midrank" else NULL,
                       na.action = spec$na.action,
                       summary = if (synchronized) resolved$definition else NULL,
                       exposure = exposure$definition, localenv = localenv$definition))
}
