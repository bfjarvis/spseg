test_that("count preparation is reused during population inference", {
  d <- data.frame(x = 1:3, y = 0, a = c(1, NA, 3), b = c(3, 2, 1),
    c = c(2, 1, 2), e = c(2, 3, 2))
  original <- spseg:::.spseg_counts
  calls <- 0L
  testthat::local_mocked_bindings(.spseg_counts = function(spec, data) {
    calls <<- calls + 1L
    original(spec, data)
  }, .package = "spseg")
  # Omitted count rows still determine each variable's eligible population.
  specs <- list(first = spseg_categorical(counts = c("a", "b"), na.action = "omit"),
    second = spseg_categorical(counts = c("c", "e")))
  expect_error(spseg(d, specs, coords = c("x", "y"), bands = 0),
    "populations differ")
  expect_identical(calls, 2L)
  calls <- 0L
  result <- spseg(d, specs, coords = c("x", "y"), population = rep(4, 3),
    bands = 0, output = "full")
  expect_identical(calls, 2L)
  expect_equal(result$analyses$first$definition$eligible_population, 8)
  expect_true(all(is.na(result$analyses$first$env[[1]][2, ])))
})

test_that("surface setup validates the current option vocabulary", {
  geom <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 2, ymax = 1))), n = c(2, 1))
  x <- sf::st_sf(a = c(4, 1), b = c(1, 4), geometry = geom)
  run <- function(...) spseg(x, list(race = spseg_categorical(counts = c("a", "b"))),
    bands = 0, surface = "pycno", ...)
  expect_error(run(keep_geometry = TRUE), "surface options")
  expect_error(run(surface_geometry = "polygon"), "surface_geometry")
  expect_error(run(nrow = 2), "both 'nrow' and 'ncol'")
  expect_error(run(celldim = Inf), "positive finite")
  expect_error(run(celldim = NA_real_), "positive finite")
  expect_error(run(nrow = 2, ncol = 4, max_iter = 0), "positive integer")
  expect_error(run(nrow = 2, ncol = 4, converge = NA_real_), "finite numeric")
})
