#' Describe variables for a shared spatial segregation calculation
#'
#' These constructors describe data columns; they do not search neighborhoods.
#' Supply exactly one of `values` (a single column) and `counts` (columns in
#' category/value order). Column selectors are names or integer positions.
#'
#' @param values Column containing category membership or observed values.
#' @param counts Columns containing nonnegative population counts.
#' @param levels Category order for observed ordinal or interval membership.
#'   Ordered factors supply their own order. For categorical values, factor
#'   levels or sorted distinct values are used when this is omitted.
#' @param measures Optional measure selection overriding [spseg()].
#' @param weights Population-weight column or numeric vector for `values`.
#'   Defaults to the common `population` supplied to [spseg()], or one per row.
#' @param na.action Missing-value policy: `"fail"` or `"omit"`. Omission applies
#'   only to this variable, leaving the common neighborhood population intact.
#' @param interpolation Grouped continuous CDF interpolation, currently `"linear"`
#'   in overall percentile space, not in the original measurement units.
#' @param scores Exact numerical values corresponding to continuous count
#'   columns. These are observed values, not interval midpoints.
#' @param exposure For observed continuous variables, a list with either
#'   `probs` (global weighted quantile boundaries, including 0 and 1) or `breaks`
#'   (value boundaries covering the observations). Defaults to deciles.
#'   Equal values stay together; empty groups are removed. These cuts affect
#'   exposure only, not the rank-order indices.
#' @return A `seg_variable` specification for `spseg(variables = ...)`.
#' @seealso [spseg()]
#' @examples
#' seg_categorical(values = "ethnicity")
#' seg_ordinal(counts = c("primary", "secondary", "tertiary"))
#' seg_grouped(counts = c("low", "middle", "high"))
#' seg_continuous(values = "income", exposure = list(probs = c(0, .2, .8, 1)))
#' @name seg_variables
NULL

.seg_variable <- function(type, values, counts, levels, measures, weights,
                          na.action, interpolation = NULL, scores = NULL,
                          exposure = NULL) {
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
                 exposure = exposure), class = "seg_variable")
}

#' @rdname seg_variables
#' @export
seg_categorical <- function(values = NULL, counts = NULL, levels = NULL,
                            measures = NULL, weights = NULL,
                            na.action = c("fail", "omit")) {
  .seg_variable("categorical", values, counts, levels, measures, weights,
                match.arg(na.action))
}

#' @rdname seg_variables
#' @export
seg_ordinal <- function(values = NULL, counts = NULL, levels = NULL,
                       measures = NULL, weights = NULL,
                       na.action = c("fail", "omit")) {
  .seg_variable("ordinal", values, counts, levels, measures, weights,
                match.arg(na.action))
}

#' @rdname seg_variables
#' @export
seg_grouped <- function(values = NULL, counts = NULL, levels = NULL,
                       measures = NULL, weights = NULL,
                       na.action = c("fail", "omit"),
                       interpolation = "linear") {
  if (!identical(interpolation, "linear")) {
    stop("Only 'linear' interpolation is supported", call. = FALSE)
  }
  .seg_variable("grouped", values, counts, levels, measures, weights,
                match.arg(na.action), interpolation = interpolation)
}

#' @rdname seg_variables
#' @export
seg_continuous <- function(values = NULL, counts = NULL, scores = NULL,
                          measures = NULL, weights = NULL,
                          na.action = c("fail", "omit"), exposure = NULL) {
  if (!is.null(values) && !is.null(scores)) {
    stop("'scores' applies only to continuous count columns", call. = FALSE)
  }
  .seg_variable("continuous", values, counts, NULL, measures, weights,
                match.arg(na.action), scores = scores, exposure = exposure)
}

.seg_select <- function(data, columns) {
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

.seg_mass <- function(value, data, n, default = rep(1, n)) {
  if (is.null(value)) return(default)
  if (is.character(value) && length(value) == 1L)
    value <- .seg_select(data, value)[[1L]]
  if (!is.numeric(value) || length(value) != n ||
      any(!is.finite(value) | value < 0)) {
    stop("Population weights must be finite, nonnegative, and have one value per row",
         call. = FALSE)
  }
  as.numeric(value)
}

.seg_counts <- function(spec, data) {
  out <- .seg_select(data, spec$counts)
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

.seg_measures <- function(measures) {
  measures[measures %in% c("relative_diversity", "variance_ratio")] <- "diversity"
  valid <- c("exposure", "information", "diversity", "dissimilarity", "square_root")
  if (!is.character(measures) || !length(measures) || anyNA(measures) ||
      any(!measures %in% c(valid, "all")))
    stop("Unknown or empty measure selection", call. = FALSE)
  if ("all" %in% measures) valid else unique(measures)
}

.seg_exposure <- function(scores, totals, options) {
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

.seg_compile <- function(spec, data, population, measures) {
  n <- nrow(data)
  selected <- .seg_measures(spec$measures %||% measures)
  counts <- matrix(numeric(), 0L, 0L)
  code <- rep(-1L, n)
  scores <- numeric()
  if (!is.null(spec$counts)) {
    counts <- .seg_counts(spec, data)
    labels <- colnames(counts)
    mass <- rowSums(counts)
    totals <- colSums(counts)
    if (spec$type == "continuous") {
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
    values <- .seg_select(data, spec$values)[[1L]]
    mass <- .seg_mass(spec$weights, data, n, population)
    bad <- is.na(values)
    if (any(bad) && spec$na.action == "fail")
      stop("Missing observed values; choose na.action = 'omit' explicitly", call. = FALSE)
    mass[bad] <- 0
    if (spec$type == "continuous") {
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
  type <- match(spec$type, c("categorical", "ordinal", "grouped", "continuous")) - 1L
  method <- c("categorical", "ordinal_cutoffs", "linear_cdf", "empirical_midrank")[type + 1L]
  exposure <- if (type == 3L) {
    if ("exposure" %in% selected) .seg_exposure(scores, totals, spec$exposure) else
      list(map = integer(), labels = character(), definition = NULL)
  } else list(map = seq_along(labels) - 1L, labels = labels, definition = NULL)
  list(type = type, counts = counts, code = as.integer(code), mass = mass,
       totals = totals, shares = shares, scores = scores,
       z = cumsum(shares) - shares / 2,
       exposure_map = as.integer(exposure$map), exposure_labels = exposure$labels,
       flags = as.integer(c("exposure", "information", "diversity", "dissimilarity",
                            "square_root") %in% selected),
       metadata = list(type = spec$type, method = method, measures = selected,
                       labels = labels, shares = shares, scores = scores,
                       ties = if (type == 3L) "weighted_midrank" else NULL,
                       na.action = spec$na.action,
                       exposure = exposure$definition))
}
