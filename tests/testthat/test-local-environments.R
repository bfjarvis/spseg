test_that("weighted local quantiles use global midranks and preserve ties", {
  d <- data.frame(x = c(0, 1, 3, 4), y = 0,
    v = c(0, 0, 20, 100), w = c(1, 3, 2, 4))
  summary <- list(type = "quantiles", probs = c(.1, .5, .9))
  run <- function(d, spec) {
    spseg(d, list(v = spec), coords = c("x", "y"), bands = c(1, 4),
      weighting = "unweighted", output = "full")
  }
  r <- run(d, spseg_continuous(vars = "v", weights = "w", summary = summary))
  expect_equal(unname(r$analyses$v$env[[1]][1, ]), rep(.2, 3))
  expect_equal(unname(r$analyses$v$env[[1]][3, ]), c(.5, .8, .8))
  expect_equal(unname(r$analyses$v$env[[2]][1, ]), c(.2, .5, .8))
  expect_identical(r$analyses$v$definition$localenv$scale, "global_rank")
  counts <- transform(d, low = ifelse(v == 0, w, 0),
    mid = ifelse(v == 20, w, 0), high = ifelse(v == 100, w, 0))
  aggregate <- run(counts, spseg_continuous(vars = c("low", "mid", "high"),
    scores = c(0, 20, 100), summary = summary))
  expect_equal(r$analyses$v$env, aggregate$analyses$v$env, tolerance = 1e-12)
  expect_equal(r$analyses$v$indices, aggregate$analyses$v$indices, tolerance = 1e-12)
  d$v <- exp(d$v / 100)
  transformed <- run(d, spseg_continuous(vars = "v", weights = "w", summary = summary))
  expect_equal(r$analyses$v$env, transformed$analyses$v$env)
})

test_that("interval quantiles and shares use the same linear CDF model", {
  d <- data.frame(x = 1:3, y = 0, a = c(2, 0, 0), b = c(0, 1, 1), c = c(0, 1, 5))
  run <- function(d, columns, summary) {
    suppressWarnings(spseg(d, list(v = spseg_continuous(vars = columns,
      representation = "intervals", summary = summary)), coords = c("x", "y"),
      bands = 0, weighting = "unweighted", output = "full"))
  }
  q <- run(d, c("a", "b", "c"), list(type = "quantiles", probs = c(.25, .5, .75)))
  s <- run(d, c("a", "b", "c"), list(probs = c(.1, .5)))
  expect_equal(unname(q$analyses$v$env[[1]]),
    rbind(c(.05, .1, .15), c(.3, .4, .7), c(.46, .64, .82)), tolerance = 1e-12)
  expect_equal(unname(s$analyses$v$env[[1]]),
    rbind(c(.5, .5, 0), c(0, 7/12, 5/12), c(0, 11/36, 25/36)), tolerance = 1e-12)
  scalar <- c("d", "h", "r", "s")
  expect_equal(q$analyses$v$indices$overall[scalar],
    s$analyses$v$indices$overall[scalar], tolerance = 1e-12)
  d$b1 <- d$b * .3
  d$b2 <- d$b * .7
  split <- run(d, c("a", "b1", "b2", "c"),
    list(type = "quantiles", probs = c(.25, .5, .75)))
  expect_equal(q$analyses$v$env, split$analyses$v$env, tolerance = 1e-12)
})

test_that("local summaries follow kernels and proportional kNN boundaries", {
  d <- data.frame(x = c(0, 1, -1, 3), y = 0, v = 1:4, w = c(1, 2, 4, 3))
  ranks <- (cumsum(d$w) - d$w / 2) / sum(d$w)
  probabilities <- c(.25, .5, .75)
  inverse <- function(weights) {
    present <- which(weights > 0)
    cumulative <- cumsum(weights[present]) / sum(weights)
    vapply(probabilities, function(p) ranks[present[which(cumulative >= p)[1]]],
      numeric(1))
  }
  spec <- list(v = spseg_continuous(vars = "v", weights = "w",
    summary = list(type = "quantiles", probs = probabilities)))
  radius <- spseg(d, spec, coords = c("x", "y"), bands = 4,
    weighting = "exponential", output = "full")
  expect_equal(unname(radius$analyses$v$env[[1]][1, ]),
    inverse(d$w * exp(-abs(d$x) / 4)))
  nearest <- spseg(d, spec, coords = c("x", "y"), neighbors = "knn", bands = 4,
    weighting = "unweighted", output = "full")
  expect_equal(unname(nearest$analyses$v$env[[1]][1, ]), inverse(c(1, 1, 2, 0)))
  expect_equal(nearest$analyses$v$knn$min_population, 4)
})

test_that("summary choices preserve scalar indices and ordinal behavior", {
  d <- data.frame(x = 1:6, y = 0, v = c(1, 1, 2, 3, 5, 8),
    o = ordered(c("a", "a", "b", "c", "c", "c")))
  specs <- list(default = spseg_continuous(vars = "v"),
    quantiles = spseg_continuous(vars = "v", summary = list(type = "quantiles")),
    custom = spseg_continuous(vars = "v", summary = list(cuts = 2)),
    ordinal = spseg_ordinal(vars = "o"),
    ordinal_quantiles = spseg_ordinal(vars = "o", localenv = "quantiles"))
  for (mode in c("radius", "knn")) {
    r <- spseg(d, specs, coords = c("x", "y"), bands = c(1, 3),
      neighbors = mode, output = "full")
    indices <- spseg(d, specs, coords = c("x", "y"), bands = c(1, 3), neighbors = mode)
    expect_equal(as.data.frame(r), as.data.frame(indices), tolerance = 1e-12)
    for (name in c("quantiles", "custom")) {
      scalar <- c("d", "h", "r", "s")
      expect_equal(r$analyses$default$indices$overall[scalar],
        r$analyses[[name]]$indices$overall[scalar])
    }
    expect_equal(r$analyses$ordinal$indices, r$analyses$ordinal_quantiles$indices)
    expect_equal(unname(r$analyses$ordinal$env[[1]]),
      diag(3)[as.integer(d$o), , drop = FALSE])
    expect_equal(ncol(r$analyses$ordinal_quantiles$env[[1]]), 9L)
    expect_equal(ncol(r$analyses$custom$env[[1]]), 2L)
    expect_equal(ncol(r$analyses$default$env[[1]]), 4L)
  }
})

test_that("local shares do not require exposure calculation", {
  d <- data.frame(x = 1:10, y = 0, v = 1:10, o = ordered(rep(letters[1:2], each = 5)))
  r <- spseg(d, list(v = spseg_continuous(vars = "v", summary = list(probs = .5)),
    o = spseg_ordinal(vars = "o", localenv = list(type = "shares", probs = c(0, .5, 1)))),
    coords = c("x", "y"), measures = "information", bands = 20,
    weighting = "unweighted", output = "full")
  for (a in r$analyses) {
    expect_equal(unname(a$env[[1]]), matrix(.5, 10, 2))
    expect_null(a$indices$overall$p)
  }
})

test_that("empty local summaries and mixed sf extraction retain focal rows", {
  d <- data.frame(x = c(0, 10, 20), y = 0, g = c("a", "b", "b"), v = c(NA, 1, 2))
  x <- sf::st_as_sf(d, coords = c("x", "y"), crs = 3857, sf_column_name = "location")
  r <- spseg(x, list(g = spseg_categorical(vars = "g"),
    v = spseg_continuous(vars = "v", na.action = "omit",
      summary = list(type = "quantiles"))), bands = 0, output = "localenv")
  expect_true(all(is.na(r$analyses$v$env[[1]][1, ])))
  table <- as.data.frame(r, what = "env")
  expect_equal(nrow(table), 3 * (2 + 3))
  quantiles <- subset(table, statistic == "quantile_rank")
  expect_equal(sort(unique(quantiles$prob)), c(.25, .5, .75))
  expect_false(any(c("share", "rank", "summary") %in% names(quantiles)))
  expect_true(all(is.na(quantiles$category)))
  expect_s3_class(sf::st_as_sf(r, variables = "v"), "sf")
  expect_equal(nrow(sf::st_as_sf(r, variables = "v")), 9L)
})

test_that("continuous environment storage stays compact at large bandwidths", {
  d <- data.frame(x = seq(0, 1, length.out = 500), y = 0, v = seq_len(500))
  r <- spseg(d, list(v = spseg_continuous(vars = "v",
    summary = list(type = "quantiles", probs = seq(.1, .9, .1)))),
    coords = c("x", "y"), bands = c(1, 2, 3), weighting = "unweighted", output = "full")
  expect_equal(dim(r$analyses$v$env[[3]]), c(500L, 9L))
  expect_lt(as.numeric(object.size(r)), 1024^2)
  expect_equal(r$analyses$v$definition$scores, as.numeric(d$v))
})

test_that("ordinal summary options reject invalid specifications", {
  for (p in list(numeric(), c(.5, .5), c(.8, .2), c(-.1, .5), NA_real_, "median")) {
    expect_error(spseg_ordinal(vars = "o", localenv = list(type = "quantiles", probs = p)),
      "probabilities")
  }
  expect_error(spseg_ordinal(vars = "o", localenv = "unknown"), "localenv")
  expect_error(spseg_ordinal(vars = "o", localenv = list(type = "shares", breaks = 0:1)),
    "localenv")
})
