score_example <- function(neighbors = "radius", ...) {
  d <- data.frame(x = 1:5, y = 0, a = c(4, 0, 2, 1, 0),
    b = c(1, 3, 2, 0, 0), c = c(0, 2, 1, 4, 0),
    income = c(1, 2, 4, 8, 16))
  vars <- c("a", "b", "c")
  spseg(d, coords = c("x", "y"), variables = list(
    category = spseg_categorical(vars = vars),
    ordinal = spseg_ordinal(vars = vars),
    exact = spseg_continuous(vars = "income", weights = rowSums(d[vars])),
    intervals = spseg_continuous(vars = vars, representation = "intervals",
      summary = list(probs = c(.35, .65)))),
    population = rowSums(d[vars]), neighbors = neighbors,
    bands = if (neighbors == "radius") c(0, 2.5) else c(5, 12), ...)
}

test_that("score retention is optional and preserves streaming results", {
  for (neighbors in c("radius", "knn")) {
    ordinary <- score_example(neighbors, output = "full")
    retained <- score_example(neighbors, output = "full", keep_scores = TRUE,
                              keep_contributions = TRUE)
    expect_equal(retained$diagnostics$local_environments, 5 * 2 * 4)
    expect_equal(retained$diagnostics$focal_searches, 5)
    for (name in names(ordinary$analyses)) {
      a <- ordinary$analyses[[name]]
      b <- retained$analyses[[name]]
      expect_null(a$scores)
      expect_null(a$contributions)
      expect_equal(a$indices, b$indices)
      expect_equal(a$env, b$env)
      for (key in names(b$contributions$overall)) {
        for (band in 1:2) expect_equal(
          sum(b$contributions$overall[[key]][[band]]),
          unname(b$indices$overall[[key]][band]), tolerance = 1e-12)
      }
      for (key in names(b$contributions$pairwise)) {
        pairs <- utils::combn(nrow(b$indices$pairwise[[key]][[1]]), 2)
        for (band in 1:2) expect_equal(
          unname(colSums(b$contributions$pairwise[[key]][[band]])),
          unname(b$indices$pairwise[[key]][[band]][t(pairs)]), tolerance = 1e-12)
      }
    }
  }
})

test_that("raw categorical scores follow the vector primitives", {
  result <- score_example(output = "indices", keep_scores = TRUE)
  a <- result$analyses$category
  p <- c(4, 1, 0) / 5
  P <- c(7, 6, 7) / 20
  expected <- c(h = -sum(p[p > 0] * log(p[p > 0])),
    r = sum(p * (1 - p)), s = 0, d = sum(abs(p - P)))
  for (key in names(expected)) expect_equal(
    a$scores$overall[[key]][[1]][1], unname(expected[key]))
  expect_null(a$env)
  expect_null(a$contributions)
  expect_true(all(is.na(vapply(a$scores$overall,
    function(z) z[[1]][5], numeric(1)))))
  rows <- as.data.frame(result, what = "scores")
  expect_true(all(c("location", "scope", "measure", "value") %in% names(rows)))
  expect_true(all(c("overall", "pairwise") %in% rows$scope))
  expect_error(as.data.frame(result, what = "contributions"), "keep_contributions")
})

test_that("contributions can be retained without indices or environments", {
  result <- score_example(output = "localenv", keep_contributions = TRUE,
                          scope = "overall", measures = "information")
  reference <- score_example(scope = "overall", measures = "information")
  for (name in names(result$analyses)) {
    a <- result$analyses[[name]]
    expect_length(a$indices$overall, 0)
    expect_null(a$scores)
    expect_named(a$contributions$overall, "h")
    expect_length(a$contributions$pairwise, 0)
    expect_equal(vapply(a$contributions$overall$h, sum, numeric(1)),
      reference$analyses[[name]]$indices$overall$h, tolerance = 1e-12)
  }
  indices_only <- score_example(keep_contributions = TRUE)
  expect_null(indices_only$analyses$category$env)
  expect_null(indices_only$analyses$category$scores)
  expect_gt(nrow(as.data.frame(indices_only, what = "contributions")), 0)
  expect_error(score_example(keep_scores = NA), "keep_scores")
  expect_error(score_example(keep_contributions = 1), "keep_contributions")
})

test_that("absent groups and zero-population focal units have defined handling", {
  d <- data.frame(x = 1:4, y = 0, a = c(4, 0, 0, 0),
    b = c(0, 3, 0, 0), c = c(0, 0, 2, 0), absent = 0)
  result <- spseg(d, list(g = spseg_categorical(vars = names(d)[3:6])),
    coords = c("x", "y"), bands = 0, keep_scores = TRUE,
    keep_contributions = TRUE, output = "full")
  a <- result$analyses$g
  for (key in names(a$contributions$pairwise)) {
    values <- a$contributions$pairwise[[key]][[1]]
    expect_true(all(is.na(values[, c(3, 5, 6)])))
    expect_equal(unname(values[4, c(1, 2, 4)]), c(0, 0, 0))
    expect_equal(values[3, 1], 0)
    expect_true(is.na(a$scores$pairwise[[key]][[1]][3, 1]))
  }
  pair_only <- spseg(d, list(g = spseg_categorical(vars = c("a", "b"))),
    coords = c("x", "y"), bands = 0, scope = "pairwise", keep_scores = TRUE)
  expect_length(pair_only$analyses$g$scores$overall, 0)
  expect_named(pair_only$analyses$g$indices$overall, "p")
})

test_that("retained local results preserve geographic and network row keys", {
  x <- sf::st_as_sf(data.frame(x = 1:3, y = 0,
    a = c(2, 1, 0), b = c(0, 1, 2)), coords = c("x", "y"), crs = 3857)
  result <- spseg(x, list(g = spseg_categorical(vars = c("a", "b"))),
    bands = c(0, 2), output = "full", keep_contributions = TRUE)
  rows <- as.data.frame(result, what = "contributions", bands = 2)
  spatial <- sf::st_as_sf(result, what = "contributions", bands = 2)
  expect_equal(sf::st_geometry(spatial), sf::st_geometry(x)[rows$location])
  net <- spseg_network(data.frame(from = c("a", "b"), to = c("b", "c")))
  d <- data.frame(node = c("c", "a", "b"), a = c(0, 2, 1), b = c(2, 0, 1))
  result <- spseg(d, list(g = spseg_categorical(vars = c("a", "b"))),
    network = net, node_id = "node", bands = 1, keep_scores = TRUE)
  rows <- as.data.frame(result, what = "scores")
  expect_equal(rows$node, d$node[rows$location])
})

test_that("constant variables and exposure-only requests keep compact arrays", {
  d <- data.frame(x = 1:3, y = 0, a = c(1, 2, 3), b = 0)
  run <- function(...) spseg(d,
    list(g = spseg_categorical(vars = c("a", "b"))), coords = c("x", "y"),
    bands = 0, keep_scores = TRUE, keep_contributions = TRUE, ...)
  constant <- run()
  for (values in constant$analyses$g$contributions$overall)
    expect_true(all(is.na(values[[1]])))
  exposure <- run(measures = "exposure")
  expect_length(exposure$analyses$g$scores$overall, 0)
  expect_length(exposure$analyses$g$scores$pairwise, 0)
  expect_equal(nrow(as.data.frame(exposure, what = "scores")), 0)
})
