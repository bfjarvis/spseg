auto_representation_data <- function() {
  data.frame(x = c(0, 0, 10, 10), y = 0,
    value = c(10, 20, 20, 30),
    band = ordered(c("low", "middle", "middle", "high"),
                   levels = c("low", "middle", "high")),
    weight = c(2, 1, 3, 2))
}

test_that("numeric and ordered microdata select their continuous representations", {
  d <- auto_representation_data()
  result <- suppressWarnings(spseg(d, list(
    numeric_auto = spseg_continuous(vars = "value", weights = "weight"),
    numeric_exact = spseg_continuous(vars = "value", weights = "weight",
                                    representation = "exact"),
    ordered_auto = spseg_continuous(vars = "band", weights = "weight"),
    ordered_intervals = spseg_continuous(vars = "band", weights = "weight",
                                        representation = "intervals")
  ), coords = c("x", "y"), bands = c(0, 20), output = "full"))
  expect_equal(result$analyses$numeric_auto, result$analyses$numeric_exact)
  expect_equal(result$analyses$ordered_auto, result$analyses$ordered_intervals)
  expect_identical(result$analyses$numeric_auto$definition$representation, "exact")
  expect_identical(result$analyses$ordered_auto$definition$representation, "intervals")
  expect_identical(result$analyses$ordered_auto$definition$interpolation, "linear")
})

test_that("count columns select intervals unless exact scores are supplied", {
  d <- data.frame(x = c(0, 10), y = 0, low = c(2, 0),
                  middle = c(1, 3), high = c(0, 2))
  columns <- c("low", "middle", "high")
  result <- suppressWarnings(spseg(d, list(
    auto = spseg_continuous(vars = columns),
    intervals = spseg_continuous(vars = columns, representation = "intervals"),
    scores = spseg_continuous(vars = columns, scores = c(10, 20, 30)),
    exact = spseg_continuous(vars = columns, scores = c(10, 20, 30),
                            representation = "exact")
  ), coords = c("x", "y"), bands = c(0, 20), output = "full"))
  expect_equal(result$analyses$auto, result$analyses$intervals)
  expect_equal(result$analyses$scores, result$analyses$exact)
  # The aggregated table represents the same weighted population as the microdata.
  micro <- suppressWarnings(spseg(auto_representation_data(), list(
    income = spseg_continuous(vars = "band", weights = "weight")
  ), coords = c("x", "y"), bands = c(0, 20), output = "full"))
  expect_equal(result$analyses$auto$indices, micro$analyses$income$indices)
  expect_error(spseg(d, list(income = spseg_continuous(vars = columns,
    representation = "exact")), coords = c("x", "y"), bands = 0), "scores")
})

test_that("explicit membership order provides an interval interpretation", {
  d <- auto_representation_data()
  d$band <- factor(d$band, levels = c("high", "low", "middle"), ordered = FALSE)
  spec <- spseg_continuous(vars = "band", levels = c("low", "middle", "high"))
  result <- suppressWarnings(spseg(d, list(income = spec), coords = c("x", "y"),
                                   bands = 0))
  expect_identical(result$analyses$income$definition$representation, "intervals")
  expect_identical(result$analyses$income$definition$labels, c("low", "middle", "high"))
  expect_error(spseg(d, list(income = spseg_continuous(vars = "band")),
    coords = c("x", "y"), bands = 0), "ordered factor or explicit 'levels'")
  expect_error(spseg_continuous(vars = "band", representation = "exact",
    levels = c("low", "high")), "levels")
  expect_error(spseg(d, list(income = spseg_continuous(vars = "value",
    interpolation = "linear")), coords = c("x", "y"), bands = 0), "interpolation")
  expect_error(spseg_continuous(vars = c("low", "high"),
    representation = "intervals", scores = c(1, 2)), "scores")
})

test_that("automatic interval selection is preserved through gridding", {
  geometry <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 2, ymax = 1))), n = c(2, 1))
  d <- sf::st_sf(low = c(2, 1), high = c(1, 2), geometry = geometry)
  for (surface in c("grid", "pycno")) {
    result <- suppressWarnings(spseg(d, list(
      auto = spseg_continuous(vars = c("low", "high")),
      explicit = spseg_continuous(vars = c("low", "high"),
                                  representation = "intervals")
    ), bands = c(0, 2), surface = surface, nrow = 2, ncol = 4, output = "full"))
    expect_equal(result$analyses$auto, result$analyses$explicit)
    expect_identical(result$analyses$auto$definition$representation, "intervals")
  }
})
