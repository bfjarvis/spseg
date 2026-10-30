output_example <- function(variables, ...) {
  d <- data.frame(x = 1:4, y = 0, a = c(4, 2, 1, 0), b = c(0, 2, 3, 4),
    c = c(1, 1, 1, 1), value = c(1, 2, 5, 10))
  spseg(d, variables, coords = c("x", "y"), bands = c(0, 2.5),
    population = rowSums(d[c("a", "b", "c")]), ...)
}

test_that("variable settings inherit global defaults and override explicitly", {
  constructors <- list(spseg_categorical, spseg_ordinal, spseg_continuous)
  for (constructor in constructors) {
    spec <- constructor(vars = c("a", "b", "c"))
    if (identical(constructor, spseg_continuous)) spec$scores <- c(1, 2, 3)
    inherit <- output_example(list(v = spec), scope = "overall", output = "full",
      keep_scores = TRUE, keep_contributions = TRUE)$analyses$v
    expect_equal(inherit$settings, list(scope = "overall", output = "full",
      keep_scores = TRUE, keep_contributions = TRUE))
    expect_length(inherit$indices$pairwise, 0)
    expect_length(inherit$env, 2)
    spec$scope <- "pairwise"
    spec$output <- "indices"
    spec$keep_scores <- FALSE
    spec$keep_contributions <- FALSE
    override <- output_example(list(v = spec), scope = "overall", output = "full",
      keep_scores = TRUE, keep_contributions = TRUE)$analyses$v
    expect_equal(override$settings, list(scope = "pairwise", output = "indices",
      keep_scores = FALSE, keep_contributions = FALSE))
    expect_null(override$env)
    expect_null(override$scores)
    expect_null(override$contributions)
    expect_named(override$indices$overall, "p")
    expect_named(override$indices$pairwise, c("d", "r", "h", "s"))
  }
})

test_that("mixed storage retains shared input once and extracts available results", {
  variables <- list(
    full = spseg_categorical(vars = c("a", "b", "c"), output = "full",
      scope = "both", keep_scores = TRUE),
    indices = spseg_ordinal(vars = c("a", "b", "c"), scope = "overall",
      keep_contributions = TRUE),
    env = spseg_continuous(vars = "value", output = "localenv"))
  result <- output_example(variables, scope = "pairwise", output = "indices")
  for (name in names(variables)) {
    separate <- output_example(variables[name], scope = "pairwise", output = "indices")
    expect_equal(result$analyses[[name]], separate$analyses[[name]])
  }
  expect_equal(nrow(result$coords), 4)
  expect_equal(nrow(result$data), 4)
  expect_null(result$analyses$indices$env)
  expect_length(result$analyses$env$indices$overall, 0)
  expect_setequal(unique(as.data.frame(result, what = "env")$variable), c("full", "env"))
  expect_setequal(unique(as.data.frame(result)$variable), c("full", "indices"))
  expect_equal(unique(as.data.frame(result, what = "scores")$variable), "full")
  expect_equal(unique(as.data.frame(result, what = "contributions")$variable), "indices")
  expect_error(as.data.frame(result, what = "env", variables = "indices"), "not stored")
  expect_error(as.data.frame(result, variables = "env"), "not calculated")
  expect_error(as.data.frame(result, what = "scores", variables = "indices"), "not retained")
  expect_equal(result$diagnostics$focal_searches, 4)
  expect_equal(result$diagnostics$local_environments, 4 * 2 * 3)
})

test_that("variable output overrides control geometry and surface recompilation", {
  shapes <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 3, ymax = 1), crs = sf::st_crs(3857))),
    n = c(3, 1))
  x <- sf::st_sf(a = c(4, 2, 1), b = c(1, 2, 4), geometry = shapes)
  specs <- list(stored = spseg_categorical(vars = c("a", "b"), output = "full"),
    lean = spseg_ordinal(vars = c("a", "b"), output = "indices", scope = "overall"))
  for (surface in c("raw", "grid")) {
    args <- if (surface == "grid")
      list(celldim = .5, surface_geometry = "points") else list()
    result <- do.call(spseg, c(list(x = x, variables = specs, bands = 1.5,
      surface = surface, output = "localenv"), args))
    expect_equal(result$analyses$stored$settings$output, "full")
    expect_equal(result$analyses$lean$settings$output, "indices")
    expect_length(result$analyses$lean$indices$pairwise, 0)
    expect_null(result$analyses$lean$env)
    spatial <- sf::st_as_sf(result)
    expect_equal(unique(spatial$variable), "stored")
    expect_equal(length(result$geometry), nrow(result$coords))
    specs$stored$output <- "indices"
    lean <- do.call(spseg, c(list(x = x, variables = specs, bands = 1.5,
      surface = surface, output = "full"), args))
    expect_null(lean$coords)
    expect_null(lean$geometry)
    expect_null(lean$data)
    specs$stored$output <- "full"
  }
})

test_that("constructor overrides validate choices and nullable logical flags", {
  for (constructor in list(spseg_categorical, spseg_ordinal, spseg_continuous)) {
    expect_error(constructor(vars = c("a", "b"), scope = "multigroup"), "arg")
    expect_error(constructor(vars = c("a", "b"), output = "unknown"), "arg")
    expect_error(constructor(vars = c("a", "b"), keep_scores = NA), "keep_scores")
    expect_error(constructor(vars = c("a", "b"), keep_contributions = 1), "keep_contributions")
    spec <- constructor(vars = c("a", "b"), scope = "both", output = "full",
      keep_scores = FALSE, keep_contributions = TRUE)
    expect_equal(spec$scope, "both")
    expect_false(spec$keep_scores)
    expect_true(spec$keep_contributions)
  }
})

test_that("pairwise-only ordered input preserves exposure independently", {
  variables <- list(
    ordinal = spseg_ordinal(vars = c("a", "b", "c")),
    intervals = spseg_continuous(vars = c("a", "b", "c"),
      representation = "intervals", summary = list(probs = c(.35, .8))),
    exact = spseg_continuous(vars = "value"),
    ranks = spseg_continuous(vars = "value",
      summary = list(type = "quantiles")))
  both <- output_example(variables, scope = "both")
  pairwise <- output_example(variables, scope = "pairwise",
    keep_scores = TRUE, keep_contributions = TRUE)
  for (name in names(variables)) {
    a <- pairwise$analyses[[name]]
    b <- both$analyses[[name]]
    expect_equal(a$indices$pairwise, b$indices$pairwise)
    expect_null(a$indices$overall$h)
    expect_length(a$scores$overall, 0)
    expect_length(a$contributions$overall, 0)
    expect_equal(a$indices$overall$p, b$indices$overall$p)
    expect_equal(a$rank_summary, b$rank_summary)
  }
  expect_equal(nrow(as.data.frame(pairwise)), 0)
})

test_that("network extraction retains keys for variable-specific local output", {
  d <- data.frame(node = c("c", "a", "b"), a = c(0, 2, 1), b = c(2, 0, 1))
  net <- spseg_network(data.frame(from = c("a", "b"), to = c("b", "c")))
  result <- spseg(d, list(
    scores = spseg_categorical(vars = c("a", "b"), keep_scores = TRUE),
    env = spseg_ordinal(vars = c("a", "b"), output = "localenv")),
    network = net, node_id = "node", bands = 1)
  for (what in c("scores", "env")) {
    rows <- as.data.frame(result, what = what)
    expect_equal(rows$node, d$node[rows$location])
    expect_equal(unique(rows$variable), what)
  }
  expect_equal(unique(as.data.frame(result)$variable), "scores")
})
