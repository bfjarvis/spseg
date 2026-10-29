compact_fixture <- function(output = "full", measures = "all") {
  d <- data.frame(x = c(0, 1, 3, 5), y = 0,
    a = c(8, 6, 1, 2), b = c(2, 4, 9, 8), income = c(1, 2, 4, 8))
  spseg(d, list(race = spseg_categorical(vars = c("a", "b")),
    income = spseg_continuous(vars = "income",
      summary = list(type = "quantiles", probs = .5, mean_rank = TRUE))),
    coords = c("x", "y"), population = rep(10, 4), bands = c(0, 2),
    output = output, measures = measures)
}

test_that("compact tables expose results rather than duplicated metadata", {
  r <- compact_fixture()
  expect_identical(names(as.data.frame(r)), c("variable", "band", "measure", "value"))
  expect_identical(names(as.data.frame(r, what = "pairwise")),
    c("variable", "band", "measure", "group_a", "group_b", "value"))
  expect_identical(names(as.data.frame(r, what = "exposure", variables = "race")),
    c("variable", "band", "from_group", "to_group", "value"))
  expect_identical(names(as.data.frame(r, what = "exposure", variables = "income")),
    c("variable", "band", "from_group", "statistic", "prob", "value"))
  expect_identical(names(as.data.frame(r, what = "env", variables = "race")),
    c("variable", "band", "location", "category", "value"))
  expect_identical(names(as.data.frame(r, what = "env", variables = "income")),
    c("variable", "band", "location", "statistic", "prob", "value"))
  env <- as.data.frame(r, what = "env", variables = "race", bands = 0)
  expect_equal(env$value, as.vector(r$analyses$race$env[[1]]))
  exposure <- as.data.frame(r, what = "exposure", variables = "income", bands = 2)
  expect_equal(exposure$value, as.vector(r$analyses$income$rank_summary[[2]]))
  expect_identical(r$analyses$income$definition$representation, "exact")
  expect_true(length(r$analyses$income$normalizers) > 0)
  expect_output(print(r), "not exposure probabilities")
})

test_that("mixed summaries include only their semantic union of columns", {
  r <- compact_fixture()
  env <- as.data.frame(r, what = "env")
  expect_identical(names(env),
    c("variable", "band", "location", "category", "statistic", "prob", "value"))
  expect_true(all(env$statistic[env$variable == "race"] == "share"))
  expect_true(all(is.na(env$category[env$variable == "income"])))
  expect_true(all(is.na(env$prob[env$variable == "race"])))
  exposure <- as.data.frame(r, what = "exposure")
  expect_identical(names(exposure),
    c("variable", "band", "from_group", "to_group", "statistic", "prob", "value"))
  expect_true(all(exposure$statistic[exposure$variable == "race"] == "share"))
  expect_true(all(is.na(exposure$to_group[exposure$variable == "income"])))
  expect_true(all(is.na(exposure$prob[exposure$statistic == "mean_rank"])))
})

test_that("empty and all-missing tables retain meaningful typed columns", {
  r <- compact_fixture(measures = "information")
  empty <- as.data.frame(r, what = "exposure", variables = "race")
  expect_identical(names(empty), c("variable", "band", "from_group", "to_group", "value"))
  expect_equal(nrow(empty), 0L)
  expect_type(empty$value, "double")
  empty_rank <- as.data.frame(r, what = "exposure", variables = "income")
  expect_identical(names(empty_rank),
    c("variable", "band", "from_group", "statistic", "prob", "value"))
  expect_type(empty_rank$prob, "double")
  expect_type(empty_rank$value, "double")
  constant <- spseg(data.frame(x = 0:1, y = 0, a = 1, b = 0),
    list(race = spseg_categorical(vars = c("a", "b"))),
    coords = c("x", "y"), bands = 0)
  indices <- as.data.frame(constant)
  expect_identical(names(indices), c("variable", "band", "measure", "value"))
  expect_true(all(is.na(indices$value)))
  expect_true(nrow(indices) > 0)
  expect_identical(as.data.frame(constant, what = "diagnostics")$status, "constant")
  empty_env <- spseg(data.frame(x = 0:1, y = 0, v = 1:2),
    list(v = spseg_continuous(vars = "v",
      summary = list(type = "quantiles", probs = c(0, 1)))),
    coords = c("x", "y"), bands = 0, output = "localenv")
  table <- as.data.frame(empty_env, what = "env")
  expect_equal(nrow(table), 0L)
  expect_identical(names(table),
    c("variable", "band", "location", "statistic", "prob", "value"))
  expect_type(table$location, "integer")
  expect_type(table$value, "double")
})

test_that("diagnostics are separated at the appropriate row granularity", {
  r <- compact_fixture()
  expect_identical(names(as.data.frame(r, what = "diagnostics")),
    c("variable", "band", "status"))
  expect_equal(nrow(as.data.frame(r, what = "diagnostics")), 4L)
  expect_equal(nrow(as.data.frame(r, what = "interpolation")), 0L)
  d <- data.frame(x = 0:2, y = 0, a = c(1, 0, 3), b = c(1, 2, 0))
  variables <- list(income = spseg_continuous(vars = c("a", "b"),
    representation = "intervals"))
  r <- suppressWarnings(spseg(d, variables, coords = c("x", "y"),
    bands = c(1, 3), neighbors = "knn", output = "full"))
  diagnostics <- as.data.frame(r, what = "diagnostics", bands = 3)
  expect_equal(nrow(diagnostics), 1L)
  expect_equal(diagnostics$requested, 3)
  expect_equal(diagnostics$min_population, 3)
  interpolation <- as.data.frame(r, what = "interpolation")
  expect_equal(nrow(interpolation), 6L)
  expect_identical(names(interpolation), c("variable", "band", "location",
    "bottom_category_share", "top_category_share", "interpolated_share"))
  expect_equal(interpolation$top_category_share[interpolation$band == 1],
    r$analyses$income$interpolation[[1]][, "top_category_share"], ignore_attr = TRUE)
  expect_false("interpolated_share" %in% names(as.data.frame(r, what = "env")))
  indices <- suppressWarnings(spseg(d, variables, coords = c("x", "y"), bands = 1))
  expect_error(as.data.frame(indices, what = "interpolation"), "not stored")
  expect_equal(nrow(as.data.frame(indices, what = "diagnostics")), 1L)
  expect_equal(nrow(as.data.frame(compact_fixture(output = "localenv"),
    what = "diagnostics")), 4L)
})

test_that("compact spatial and network tables preserve join identifiers", {
  d <- data.frame(x = 0:2, y = 0, group = c("a", "b", "a"), node = letters[1:3])
  points <- sf::st_as_sf(d, coords = c("x", "y"), crs = 3857)
  r <- spseg(points, list(group = spseg_categorical(vars = "group")),
    bands = c(0, 1), output = "full")
  mapped <- sf::st_as_sf(r, bands = 1)
  expect_identical(names(mapped), c("variable", "band", "location", "category", "value", "geometry"))
  expect_equal(sf::st_geometry(mapped), sf::st_geometry(points)[mapped$location])
  network <- spseg_network(data.frame(from = c("a", "b"), to = c("b", "c")))
  net <- spseg(d, list(group = spseg_categorical(vars = "group")),
    network = network, node_id = "node", bands = c(0, 1), output = "full")
  env <- as.data.frame(net, what = "env", bands = 1)
  expect_equal(env$node, d$node[env$location])
  expect_identical(tail(names(env), 1), "node")
})
