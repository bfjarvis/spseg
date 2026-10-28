.spseg_summary_options <- function(summary) {
  if (is.null(summary)) {
    summary <- list()
  }
  if (is.character(summary) && length(summary) == 1L) {
    summary <- list(type = summary)
  }
  if (!is.list(summary) ||
      (length(summary) && (is.null(names(summary)) || anyNA(names(summary)) ||
       any(!nzchar(names(summary))) || anyDuplicated(names(summary)) ||
       any(!names(summary) %in% c("type", "probs", "cuts", "mean_rank"))))) {
    stop("'summary' must contain type, probs or cuts, and optional mean_rank",
         call. = FALSE)
  }
  type <- summary$type %||% "shares"
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("shares", "quantiles")) {
    stop("summary 'type' must be 'shares' or 'quantiles'", call. = FALSE)
  }
  if (!is.null(summary$probs) && !is.null(summary$cuts)) {
    stop("Supply only one of summary 'probs' and 'cuts'", call. = FALSE)
  }
  mean_rank <- summary$mean_rank %||% FALSE
  if (!is.logical(mean_rank) || length(mean_rank) != 1L || is.na(mean_rank) ||
      (mean_rank && type != "quantiles")) {
    stop("summary 'mean_rank' must be TRUE or FALSE and requires quantiles mode",
         call. = FALSE)
  }
  if (!is.null(summary$cuts)) {
    cuts <- summary$cuts
    if (type != "shares") {
      stop("summary 'cuts' requires shares mode; use 'probs' for rank quantiles",
           call. = FALSE)
    }
    if (!is.numeric(cuts) || !length(cuts) || any(!is.finite(cuts)) ||
        is.unsorted(cuts, strictly = TRUE)) {
      stop("summary 'cuts' must be finite and strictly increasing", call. = FALSE)
    }
    return(list(type = type, probs = NULL, cuts = cuts, mean_rank = mean_rank))
  }
  probs <- summary$probs %||% c(.25, .5, .75)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || is.unsorted(probs, strictly = TRUE)) {
    stop("summary 'probs' must increase strictly within [0, 1]", call. = FALSE)
  }
  list(type = type, probs = as.numeric(probs[probs > 0 & probs < 1]),
       cuts = NULL, mean_rank = mean_rank)
}

.spseg_resolve_summary <- function(options, labels, totals, scores, intervals,
                                   variable = NULL) {
  requested <- options$probs
  effective <- requested
  shares <- as.numeric(totals / sum(totals))
  cumulative <- cumsum(shares)
  reasons <- character()
  fractions <- matrix(numeric(), 0L, 0L)
  split_intervals <- integer()
  if (intervals) {
    if (!is.null(options$cuts)) {
      stop("Absolute summary 'cuts' require exact values, not interval labels",
           call. = FALSE)
    }
    group_labels <- paste0("(", utils::head(c(0, effective, 1), -1L), ",",
                           utils::tail(c(0, effective, 1), -1L), "]")
    lower <- c(0, utils::head(cumulative, -1L))
    edges <- c(0, requested, 1)
    fractions <- matrix(0, length(labels), length(group_labels),
                         dimnames = list(labels, group_labels))
    # Allocate population mass, not income values; open monetary tails need no bounds.
    for (j in seq_along(group_labels)) {
      overlap <- pmax(0, pmin(cumulative, edges[j + 1L]) - pmax(lower, edges[j]))
      fractions[shares > 0, j] <- overlap[shares > 0] / shares[shares > 0]
    }
    positive <- which(shares > 0)
    allocated <- rowSums(fractions)
    representable <- positive[allocated[positive] > 0]
    fractions[representable, ] <- fractions[representable, , drop = FALSE] /
      allocated[representable]
    # Tiny positive tails can collapse to one rank at floating-point precision.
    collapsed <- positive[allocated[positive] == 0]
    if (length(collapsed)) {
      groups <- findInterval(cumulative[collapsed], edges, all.inside = TRUE)
      fractions[cbind(collapsed, groups)] <- 1
    }
    split_intervals <- which(rowSums(fractions > 8 * .Machine$double.eps) > 1L)
    if (length(split_intervals)) {
      reasons <- "linear_interpolation"
    }
    map <- integer()
    group_definition <- list(labels = group_labels,
      shares = as.numeric(colSums(fractions * shares)),
      interval_ids = lapply(seq_along(group_labels), function(j) which(fractions[, j] > 0)),
      interval_labels = labels, zero_mass_intervals = which(totals == 0),
      source_shares = shares, allocation = fractions,
      method = "linear_global_rank")
  } else {
    exposure <- .spseg_rank_groups(scores, totals,
      if (is.null(options$cuts)) list(probs = c(0, requested, 1)) else
        list(breaks = c(-Inf, options$cuts, Inf)))
    map <- exposure$map
    group_labels <- exposure$labels
    group_definition <- exposure$definition
  }
  group_definition$requested_groups <- if (is.null(options$cuts)) {
    length(requested) + 1L
  } else length(options$cuts) + 1L
  group_definition$effective_groups <- length(group_labels)
  group_definition$boundaries <- c(0,
    utils::head(cumsum(group_definition$shares), -1L), 1)
  statistics <- if (options$type == "shares") {
    data.frame(statistic = "share", prob = NA_real_, label = group_labels)
  } else {
    data.frame(statistic = c(rep("quantile_rank", length(effective)),
                            if (options$mean_rank) "mean_rank" else character()),
               prob = c(effective, if (options$mean_rank) NA_real_),
               label = c(if (length(effective)) paste0("q", effective) else character(),
                         if (options$mean_rank) "mean_rank"))
  }
  definition <- list(type = options$type, requested_probs = requested,
    effective_probs = effective, cuts = options$cuts, mean_rank = options$mean_rank,
    groups = group_definition, reasons = reasons,
    interpolation = list(needed = length(split_intervals) > 0L,
                         split_intervals = split_intervals),
    status = if (length(split_intervals)) "interpolated" else "ok",
    statistics = statistics)
  if (length(reasons)) {
    warning("Interval summary interpolation",
      if (!is.null(variable)) paste0(" for '", variable, "'"),
      ": requested groups proportionally split income categories in global-rank space. ",
      "The same fractions apply at every location; see what = 'interpolation' ",
      "and definition$summary$groups$allocation.", call. = FALSE)
  }
  list(map = as.integer(map), fractions = fractions, labels = group_labels,
       definition = definition)
}
