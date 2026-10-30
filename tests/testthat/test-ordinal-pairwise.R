test_that("ordinal pairs and exposure reuse categorical category shares", {
  d <- data.frame(x = 1:5, y = 0, low = c(4, 0, 1, 0, 0),
    middle = c(1, 0, 2, 0, 0), high = c(0, 0, 1, 5, 0), absent = 0)
  vars <- c("low", "middle", "high", "absent")
  for (neighbors in c("radius", "knn")) {
    result <- spseg(d, coords = c("x", "y"), variables = list(
      nominal = spseg_categorical(vars = vars),
      ordered = spseg_ordinal(vars = vars),
      quantiles = spseg_ordinal(vars = vars, localenv = "quantiles")),
      neighbors = neighbors,
      bands = if (neighbors == "radius") c(0, 2.5) else c(5, 10),
      output = "full", keep_scores = TRUE, keep_contributions = TRUE)
    reference <- result$analyses$nominal
    for (name in c("ordered", "quantiles")) {
      actual <- result$analyses[[name]]
      expect_named(actual$indices$pairwise, c("d", "r", "h", "s"))
      expect_equal(actual$indices$pairwise, reference$indices$pairwise)
      expect_equal(actual$indices$overall$p, reference$indices$overall$p)
      expect_equal(actual$scores$pairwise, reference$scores$pairwise)
      expect_equal(actual$contributions$pairwise, reference$contributions$pairwise)
      # A group's absent counterpart remains undefined, including in local arrays.
      expect_true(all(is.na(actual$indices$pairwise$h[[1]][, "absent"])))
    }
    expect_equal(result$diagnostics$focal_searches, nrow(d))
    expect_equal(result$diagnostics$local_environments, nrow(d) * 2 * 3)
    pairs <- as.data.frame(result, what = "pairwise", variables = "ordered")
    expect_equal(nrow(pairs), choose(length(vars), 2) * 2 * 4)
    exposure <- as.data.frame(result, what = "exposure", variables = "ordered")
    expect_equal(nrow(exposure), length(vars)^2 * 2)
  }
})

test_that("two ordinal categories yield the corresponding two-group indices", {
  d <- data.frame(x = 1:3, y = 0, low = c(4, 2, 1), high = c(1, 2, 4))
  run <- function(scope) spseg(d, coords = c("x", "y"), variables = list(
    education = spseg_ordinal(vars = c("low", "high"))),
    bands = c(0, 2.5), scope = scope)
  both <- run("both")$analyses$education
  for (key in c("h", "r", "s", "d")) {
    expect_equal(unname(both$indices$overall[[key]]),
      vapply(both$indices$pairwise[[key]], function(z) unname(z["low", "high"]),
        numeric(1)), ignore_attr = TRUE)
  }
  overall <- run("overall")$analyses$education
  expect_length(overall$indices$pairwise, 0)
  expect_equal(overall$indices$overall, both$indices$overall)
  pairs <- run("pairwise")$analyses$education
  expect_equal(pairs$indices$pairwise, both$indices$pairwise)
  expect_named(pairs$indices$overall, "p")
})

test_that("weighted ordinal memberships retain category-pair labels", {
  d <- data.frame(x = c(0, 0, 1, 2), y = 0, n = c(3, 1, 2, 4),
    level = ordered(c("low", "high", "middle", "high"),
      levels = c("low", "middle", "high")))
  result <- spseg(d, coords = c("x", "y"), variables = list(
    ordinal = spseg_ordinal(vars = "level", weights = "n"),
    category = spseg_categorical(vars = "level", weights = "n")),
    bands = c(0, 2), keep_contributions = TRUE)
  a <- result$analyses$ordinal
  expect_equal(a$indices$pairwise, result$analyses$category$indices$pairwise)
  expect_equal(a$indices$overall$p,
    result$analyses$category$indices$overall$p)
  expect_equal(rownames(a$indices$pairwise$h[[1]]), levels(d$level))
})
