polygon_point_fixture <- function() {
  concave <- sf::st_polygon(list(rbind(
    c(0, 0), c(4, 0), c(4, 1), c(1, 1), c(1, 4), c(0, 4), c(0, 0)
  )))
  square <- sf::st_polygon(list(rbind(
    c(5, 0), c(6, 0), c(6, 1), c(5, 1), c(5, 0)
  )))
  sf::st_sf(A = c(9, 1), B = c(1, 9),
            geometry = sf::st_sfc(concave, square))
}

test_that("polygon representatives default to interior points with a centroid option", {
  x <- polygon_point_fixture()
  variables <- list(group = spseg_categorical(vars = c("A", "B")))
  run <- function(points) spseg(x, variables, polygon_points = points,
    bands = 4.5, weighting = "unweighted", output = "full")
  interior <- run("point_on_surface")
  centroid <- run("centroid")
  default <- spseg(x, variables, bands = 4.5, weighting = "unweighted", output = "full")
  expect_equal(default$coords, interior$coords)
  expect_equal(default$analyses, interior$analyses)
  expect_identical(default$polygon_points, "point_on_surface")
  expect_identical(centroid$polygon_points, "centroid")
  expect_equal(interior$geometry, sf::st_geometry(x))
  points <- sf::st_as_sf(as.data.frame(interior$coords), coords = c("X", "Y"))
  expect_true(all(diag(sf::st_within(points, x, sparse = FALSE))))
  centroids <- suppressWarnings(sf::st_centroid(x))
  expect_false(sf::st_within(centroids[1, ], x[1, ], sparse = FALSE)[1, 1])
  expect_equal(centroid$coords, sf::st_coordinates(centroids))
  # The radius lies between the two representative-point distances.
  # Interior locations retain own-area composition; centroids combine both areas.
  expect_equal(unname(interior$analyses$group$env[[1]][, 1]), c(.9, .1))
  expect_equal(unname(centroid$analyses$group$env[[1]][, 1]), c(.5, .5))
})

test_that("automatic distance bands use the selected polygon representatives", {
  x <- polygon_point_fixture()
  interior <- suppressWarnings(sf::st_point_on_surface(x))
  centroid <- suppressWarnings(sf::st_centroid(x))
  interior_distance <- as.numeric(dist(sf::st_coordinates(interior)))
  centroid_distance <- as.numeric(dist(sf::st_coordinates(centroid)))
  expect_gt(interior_distance, centroid_distance)
  expect_equal(spseg_default_bands(x, probs = .5), interior_distance)
  expect_equal(spseg_default_bands(x, probs = .5, polygon_points = "centroid"),
               centroid_distance)
  variables <- list(group = spseg_categorical(vars = c("A", "B")))
  expect_equal(spseg(x, variables)$bands, rep(interior_distance, 9))
  expect_equal(spseg(x, variables, polygon_points = "centroid")$bands,
               rep(centroid_distance, 9))
})

test_that("point data and gridded surfaces do not change with polygon point selection", {
  x <- polygon_point_fixture()
  variables <- list(group = spseg_categorical(vars = c("A", "B")))
  points <- suppressWarnings(sf::st_point_on_surface(x))
  first <- spseg(points, variables, bands = 1, output = "full")
  second <- spseg(points, variables, bands = 1, output = "full",
                   polygon_points = "centroid")
  expect_equal(first$coords, second$coords)
  expect_equal(first$analyses, second$analyses)
  for (surface in c("grid", "pycno")) {
    run <- function(method) spseg(x, variables, bands = 1, surface = surface,
      polygon_points = method, nrow = 8, ncol = 12, output = "full")
    first <- run("point_on_surface")
    second <- run("centroid")
    expect_equal(first$coords, second$coords)
    expect_equal(first$data, second$data)
    expect_equal(first$analyses, second$analyses)
  }
  expect_error(spseg(x, variables, polygon_points = "unknown"), "arg")
  expect_error(spseg_default_bands(x, polygon_points = "unknown"), "arg")
})
