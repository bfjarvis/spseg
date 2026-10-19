#' Default Distance Bands
#'
#' Calculates distance quantiles from a reproducible sample of support points.
#' This is the default bandwidth calculation used by [spseg()]. Sampling uses
#' population weights when supplied, otherwise equal weights, and preserves
#' the caller's random-number state. Coincident-point distances are omitted.
#'
#' @param x A numeric two-column coordinate matrix, a data frame with columns
#'   selected by `coords`, or an [sf::sf] object. Polygon inputs use centroids.
#' @param coords Two coordinate column names or positions for a data frame.
#'   Omit for sf input or a two-column coordinate matrix.
#' @param population Optional population-weight column name or numeric vector.
#'   Defaults to one per support point. For aggregated variables, supply the
#'   relevant count totals. At least one weight must be positive.
#' @param n Maximum number of support points to sample; an integer of at least 2.
#' @param probs Requested distance quantiles, in [0, 1]. The defaults are 0.01,
#'   0.05, and 0.1 through 0.7 in increments of 0.1.
#' @return A numeric vector of distance bands, in coordinate units.
#' @seealso [spseg()]
#' @examples
#' points <- data.frame(x = 1:10, y = 0, residents = 1:10)
#' default_bands(points, coords = c("x", "y"), population = "residents")
#' @export
default_bands <- function(x, coords = NULL, population = NULL, n = 500,
                          probs = c(0.01, 0.05, seq(0.1, 0.7, 0.1))) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2 || n != floor(n))
    stop("'n' must be an integer of at least 2", call. = FALSE)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1))
    stop("'probs' must contain probabilities in [0, 1]", call. = FALSE)
  checked <- .seg_geometry(x, NULL, coords)
  coords <- checked$coords
  n_points <- nrow(coords)
  weights <- .seg_mass(population, checked$data, n_points)
  if (!is.finite(sum(weights)) || !any(weights > 0))
    stop("Population total must be positive and finite", call. = FALSE)
  sample_n <- min(n_points, n)

  if (n_points > sample_n) {
    sample_n <- min(sample_n, sum(weights > 0))
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
