summary_reference_data <- function() {
  data.frame(
    x = c(0, 10), y = 0,
    value_10 = c(2, 0), value_20 = c(1, 1),
    value_30 = c(0, 2), value_40 = c(1, 1)
  )
}

run_reference_summary <- function(summary) {
  spseg(
    summary_reference_data(),
    list(income = spseg_continuous(
      vars = c("value_10", "value_20", "value_30", "value_40"),
      scores = c(10, 20, 30, 40), summary = summary
    )),
    coords = c("x", "y"), bands = 0, weighting = "unweighted",
    output = "full"
  )
}

test_that("continuous shares and ranks match independent two-unit references", {
  shares <- run_reference_summary(list(type = "shares", probs = 0.5))
  endpoint_delimited <- run_reference_summary(list(
    type = "shares", probs = c(0, 0.5, 1)
  ))
  expect_equal(shares$analyses$income$env[[1]],
    rbind(c(0.75, 0.25), c(0.25, 0.75)), ignore_attr = TRUE)
  expect_equal(shares$analyses$income$env[[1]],
    endpoint_delimited$analyses$income$env[[1]], ignore_attr = TRUE)

  ranks <- run_reference_summary(list(
    type = "quantiles", probs = 0.5, mean_rank = TRUE
  ))
  expect_equal(ranks$analyses$income$env[[1]],
    rbind(c(0.125, 0.375), c(0.625, 0.625)), ignore_attr = TRUE)
  expect_equal(ranks$analyses$income$rank_summary[[1]],
    rbind(c(0.25, 0.4375), c(0.5, 0.5625)), ignore_attr = TRUE)
  expect_null(ranks$analyses$income$indices$multigroup$p)

  grouped <- as.data.frame(ranks, what = "exposure")
  expect_equal(grouped$statistic,
    rep(c("quantile_rank", "mean_rank"), each = 2))
  expect_equal(grouped$value, c(0.25, 0.5, 0.4375, 0.5625))
  expect_false(any(c("to_group", "exposure", "rank") %in% names(grouped)))

  # Direct pair-conditioned categorical formula, using the two unit tables.
  expect_equal(unname(shares$analyses$income$indices$pairwise$d[[1]][1, 2]),
    0.5)
})

test_that("exact value cuts agree with hand-counted ties", {
  result <- run_reference_summary(list(type = "shares", cuts = 20))
  expect_equal(result$analyses$income$env[[1]],
    rbind(c(0.75, 0.25), c(0.25, 0.75)), ignore_attr = TRUE)
  expect_error(spseg_continuous(vars = "v",
    summary = list(type = "quantiles", cuts = 20)), "requires shares")
})

test_that("interval populations split proportionally across requested rank bins", {
  d <- data.frame(x = c(0, 10), y = 0,
    a = c(4, 0), b = c(0, 2), c = c(0, 4), d = c(0, 0))
  specs <- list(income = spseg_continuous(
    representation = "intervals", vars = c("a", "b", "c", "d"),
    summary = list(type = "shares", probs = c(0.25, 0.5, 0.75))
  ))
  warning_messages <- character()
  result <- withCallingHandlers(
    spseg(d, specs, coords = c("x", "y"), bands = c(0, 2),
      weighting = "unweighted", output = "full"),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warning_messages, 1)
  expect_match(warning_messages, "Interval summary interpolation")

  analysis <- result$analyses$income
  definition <- analysis$definition$summary
  expected_allocation <- rbind(
    c(5 / 8, 3 / 8, 0, 0),
    c(0, 1 / 2, 1 / 2, 0),
    c(0, 0, 3 / 8, 5 / 8),
    c(0, 0, 0, 0)
  )
  expect_equal(definition$requested_probs, c(0.25, 0.5, 0.75))
  expect_equal(definition$effective_probs, c(0.25, 0.5, 0.75))
  expect_equal(definition$groups$allocation, expected_allocation,
    tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(rowSums(definition$groups$allocation)[1:3], rep(1, 3),
    tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(rowSums(definition$groups$allocation)[4], 0,
    ignore_attr = TRUE)

  # Global interval mass is (.4, .2, .4, 0); each requested rank bin
  # therefore receives exactly one quarter of the total population.
  global_interval_shares <- c(0.4, 0.2, 0.4, 0)
  expect_equal(colSums(definition$groups$allocation * global_interval_shares),
    rep(0.25, 4), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(definition$groups$shares, rep(0.25, 4), tolerance = 1e-12)

  expected_local <- rbind(
    c(5 / 8, 3 / 8, 0, 0),
    c(0, 1 / 6, 5 / 12, 5 / 12)
  )
  expect_equal(analysis$env[[1]], expected_local,
    tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(analysis$indices$multigroup$p[[1]], rbind(
    c(5 / 8, 3 / 8, 0, 0),
    c(3 / 8, 7 / 24, 1 / 6, 1 / 6),
    c(0, 1 / 6, 5 / 12, 5 / 12),
    c(0, 1 / 6, 5 / 12, 5 / 12)
  ), tolerance = 1e-12, ignore_attr = TRUE)
  # Groups 1 and 2 have global masses .25 each. Their focal masses are
  # (2.5, 0) and (1.5, 1), with local pair shares .625 and 0. Thus
  # D = [.8*|.625-.5| + .2*|0-.5|] / [.5] = .4.
  expect_equal(unname(analysis$indices$pairwise$d[[1]][1, 2]), 0.4,
    tolerance = 1e-12)
})

test_that("fractional interval allocation handles few and one populated bins", {
  run_shares <- function(d) {
    spseg(d, list(income = spseg_continuous(
      representation = "intervals", vars = c("a", "b", "c", "d"),
      summary = list(type = "shares", probs = c(0.25, 0.5, 0.75))
    )), coords = c("x", "y"), bands = 0, weighting = "unweighted",
    output = "full")
  }

  # Two populated intervals separated by two zero-mass intervals.
  few_data <- data.frame(x = 0:1, y = 0,
    a = c(4, 0), b = c(0, 0), c = c(0, 6), d = c(0, 0))
  warning_messages <- character()
  few <- withCallingHandlers(run_shares(few_data), warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warning_messages, 1)
  allocation <- few$analyses$income$definition$summary$groups$allocation
  expect_equal(allocation, rbind(
    c(5 / 8, 3 / 8, 0, 0), rep(0, 4),
    c(0, 1 / 6, 5 / 12, 5 / 12), rep(0, 4)
  ), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(colSums(allocation * c(0.4, 0, 0.6, 0)),
    rep(0.25, 4), tolerance = 1e-12, ignore_attr = TRUE)

  # A single populated interval still spans the complete rank distribution.
  one_data <- data.frame(x = 0:1, y = 0,
    a = c(0, 0), b = c(0, 3), c = c(0, 0), d = c(0, 0))
  warning_messages <- character()
  one <- withCallingHandlers(run_shares(one_data), warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warning_messages, 1)
  allocation <- one$analyses$income$definition$summary$groups$allocation
  expect_equal(allocation, rbind(rep(0, 4), rep(0.25, 4),
    rep(0, 4), rep(0, 4)), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(unname(one$analyses$income$env[[1]][2, ]), rep(0.25, 4),
    tolerance = 1e-12)
  expect_true(all(is.na(one$analyses$income$env[[1]][1, ])))
  expect_equal(unname(one$analyses$income$indices$pairwise$d[[1]][1, 2]), 0,
    tolerance = 1e-12)
  expect_equal(one$analyses$income$definition$summary$effective_probs,
    c(0.25, 0.5, 0.75))
})

test_that("interval summary grids leave full continuous scalar indices unchanged", {
  d <- data.frame(x = 0:2, y = 0,
    low = c(4, 0, 2), middle = c(0, 5, 0), high = c(0, 3, 4))
  variables <- list(
    quartiles = spseg_continuous(representation = "intervals",
      vars = c("low", "middle", "high"),
      summary = list(type = "shares", probs = c(0.25, 0.5, 0.75))),
    custom = spseg_continuous(representation = "intervals",
      vars = c("low", "middle", "high"),
      summary = list(type = "shares", probs = c(0.2, 0.8)))
  )
  result <- suppressWarnings(spseg(d, variables, coords = c("x", "y"),
    bands = c(0, 2), weighting = "unweighted", output = "full"))
  scalar_names <- c("d", "h", "r", "s")
  expect_equal(result$analyses$quartiles$indices$multigroup[scalar_names],
    result$analyses$custom$indices$multigroup[scalar_names], tolerance = 1e-12)
  expect_equal(result$analyses$quartiles$normalizers,
    result$analyses$custom$normalizers, tolerance = 1e-12)
})

test_that("summary grids do not change continuous scalar indices or search counts", {
  d <- data.frame(x = 0:5, y = 0, income = c(1, 1, 2, 3, 5, 8))
  specs <- list(
    default = spseg_continuous(vars = "income"),
    median = spseg_continuous(vars = "income",
      summary = list(type = "quantiles", probs = 0.5)),
    custom = spseg_continuous(vars = "income",
      summary = list(type = "shares", probs = c(0.2, 0.8)))
  )
  result <- spseg(d, specs, coords = c("x", "y"), bands = c(0, 2),
    weighting = "unweighted", output = "full")
  scalar_names <- c("d", "h", "r", "s")
  reference <- result$analyses$default
  for (name in names(specs)[-1]) {
    expect_equal(result$analyses[[name]]$indices$multigroup[scalar_names],
      reference$indices$multigroup[scalar_names], tolerance = 1e-12)
    expect_equal(result$analyses[[name]]$normalizers, reference$normalizers,
      tolerance = 1e-12)
  }
  expect_equal(result$diagnostics$focal_searches, nrow(d))

  exposure_only <- spseg(d, list(income = specs$median), coords = c("x", "y"),
    bands = c(0, 2), weighting = "unweighted", measures = "exposure",
    output = "full")
  expect_equal(as.data.frame(exposure_only, what = "exposure"),
    as.data.frame(spseg(d, list(income = specs$median), coords = c("x", "y"),
      bands = c(0, 2), weighting = "unweighted", output = "full"),
      what = "exposure"))
})

test_that("individual and count inputs give the same local and grouped summaries", {
  counts <- summary_reference_data()
  people <- data.frame(
    x = c(rep(0, 4), rep(10, 4)), y = 0,
    income = c(10, 10, 20, 40, 20, 30, 30, 40)
  )
  individual <- spseg(people, list(income = spseg_continuous(vars = "income",
    summary = list(type = "quantiles", probs = 0.5, mean_rank = TRUE))),
    coords = c("x", "y"), bands = 0, weighting = "unweighted", output = "full")
  aggregate <- spseg(counts, list(income = spseg_continuous(
    vars = c("value_10", "value_20", "value_30", "value_40"),
    scores = c(10, 20, 30, 40),
    summary = list(type = "quantiles", probs = 0.5, mean_rank = TRUE)
  )), coords = c("x", "y"), bands = 0, weighting = "unweighted", output = "full")
  expect_equal(individual$analyses$income$env[[1]][c(1, 5), ],
    aggregate$analyses$income$env[[1]], tolerance = 1e-12)
  expect_equal(individual$analyses$income$rank_summary[[1]],
    aggregate$analyses$income$rank_summary[[1]], tolerance = 1e-12)
  expect_equal(individual$analyses$income$indices, aggregate$analyses$income$indices,
    tolerance = 1e-12)
})

test_that("radius and kNN summaries agree with the brute search reference", {
  d <- data.frame(x = c(0, 1, 3, 4), y = 0,
    income = c(0, 0, 20, 100), weight = c(1, 3, 2, 4))
  for (neighborhood in c("radius", "knn")) {
    bands <- if (neighborhood == "radius") c(1, 4) else c(2, 4)
    args <- list(x = d, variables = list(income = spseg_continuous(
      vars = "income", weights = "weight",
      summary = list(type = "quantiles", probs = c(0.25, 0.5, 0.75),
        mean_rank = TRUE))), coords = c("x", "y"), bands = bands,
      neighbors = neighborhood, weighting = "unweighted", output = "full")
    fast <- do.call(spseg, args)
    args$search <- "brute"
    reference <- do.call(spseg, args)
    expect_equal(fast$analyses, reference$analyses, tolerance = 1e-12)
    expect_equal(fast$diagnostics$focal_searches, nrow(d))
  }
})
