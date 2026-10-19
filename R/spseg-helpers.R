#' Default Distance Bands
#'
#' Calculates default distance bands from sampled pairwise distances.
#'
#' This helper samples up to `n` points, calculates pairwise distances among
#' the sampled points, and returns the requested quantiles. By default,
#' sampled points are selected with probability proportional to their total
#' population, and only the lower half of the distance distribution is used.
#' It is used by [spseg()] when no bandwidths are supplied.
#'
#' @param x A numeric coordinate matrix, data frame of coordinates, or
#'   [sf::sf] object accepted by [spseg()].
#' @param data Optional population data. If omitted and `x` contains attached
#'   data, those data are used.
#' @param n Maximum number of points to sample.
#' @param probs Probabilities passed to [stats::quantile()]. The default
#'   returns the first through fifth deciles.
#' @param weighted Logical. If TRUE, sample points using row sums of `data`
#'   as population weights.
#'
#' @return A numeric vector of distance bands.
#'
#' @seealso [spseg()]
#'
#' @export
default_bands <- function(
  x,
  data = NULL,
  n = 500,
  probs = c(0.01,0.05,seq(0.1,0.7,0.1)),
  weighted = TRUE
) {
  tmp <- if (is.null(data)) {
    suppressMessages(chksegdata(x))
  } else {
    suppressMessages(chksegdata(x, data))
  }
  coords <- tmp$coords
  data <- tmp$data
  n_points <- nrow(coords)
  sample_n <- min(n_points, n)

  weights <- NULL
  if (weighted) {
    weights <- rowSums(data)
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) {
      weights <- NULL
    }
  }

  if (n_points > sample_n) {
    if (!is.null(weights)) sample_n <- min(sample_n, sum(weights > 0))
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    set.seed(1)
    on.exit(
      {
        if (has_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    coords <- coords[
      sample.int(n_points, sample_n, prob = weights),
      ,
      drop = FALSE
    ]
  }

  d <- as.numeric(dist(coords))
  d <- d[d > 0]
  if (length(d) == 0) {
    stop(
      "failed to calculate default bands from the input coordinates",
      call. = FALSE
    )
  }

  as.numeric(quantile(d, probs = probs, names = FALSE, type = 8))
}

spseg_weighting_id <- function(weighting) {
  match(weighting, c("unweighted", "biweight", "inverse", "exponential")) - 1L
}

spseg_scope_flags <- function(scope) {
  c(
    multigroup = scope %in% c("multigroup", "both"),
    pairwise = scope %in% c("pairwise", "both")
  )
}
