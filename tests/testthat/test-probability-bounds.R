test_that("grouped local probabilities respect exact bounds across kernels", {
  n <- 128L
  d <- data.frame(x = seq(0, 1, length.out = n), y = 0, value = seq_len(n))
  specs <- list(income = spseg_continuous(vars = "value"))
  for (neighbors in c("radius", "knn")) {
    for (weighting in c("biweight", "unweighted", "inverse", "exponential")) {
      result <- spseg(d, specs, coords = c("x", "y"), neighbors = neighbors,
        bands = if (neighbors == "radius") c(.03, .08, .15) else c(5, 13, 25),
        weighting = weighting, output = "full")
      a <- result$analyses$income
      for (env in a$env) {
        expect_true(all(env >= 0 & env <= 1))
        expect_equal(unname(rowSums(env)), rep(1, n), tolerance = 1e-12)
        expect_true(any(env == 0))
        expect_true(any(env > 1 - 1e-12))
      }
      for (exposure in a$indices$overall$p) {
        expect_true(all(exposure >= 0 & exposure <= 1))
        expect_equal(unname(rowSums(exposure)), rep(1, 4), tolerance = 1e-12)
      }
      env_rows <- as.data.frame(result, what = "env")
      expect_true(all(env_rows$value >= 0 & env_rows$value <= 1))
    }
  }
})

test_that("bounded summaries preserve interior precision and undefined values", {
  d <- data.frame(x = c(0, 10, 20), y = 0,
    low = c(1, 0, 0), high = c(2, 1, 0), absent = 0)
  result <- spseg(d, list(
    category = spseg_categorical(vars = c("low", "high", "absent")),
    ordinal = spseg_ordinal(vars = c("low", "high", "absent")),
    exact = spseg_continuous(vars = c("low", "high", "absent"), scores = c(1, 2, 3),
      summary = list(cuts = 1))), coords = c("x", "y"), bands = 0,
    output = "full")
  for (a in result$analyses) {
    expect_equal(unname(a$env[[1]][1, 1]), 1 / 3, tolerance = 0)
    expect_true(all(is.na(a$env[[1]][3, ])))
    finite <- a$env[[1]][is.finite(a$env[[1]])]
    expect_true(all(finite >= 0 & finite <= 1))
  }
  expect_true(all(is.na(result$analyses$category$indices$overall$p[[1]]["absent", ])))
})

test_that("rank and interpolated probability outputs retain their bounds", {
  n <- 128L
  d <- data.frame(x = seq(0, 1, length.out = n), y = 0,
    value = seq_len(n), low = rep(c(1, 0), each = n / 2),
    high = rep(c(0, 1), each = n / 2))
  result <- suppressWarnings(spseg(d, list(
    rank = spseg_continuous(vars = "value",
      summary = list(type = "quantiles", mean_rank = TRUE)),
    intervals = spseg_continuous(vars = c("low", "high"), representation = "intervals")),
    coords = c("x", "y"), bands = c(.03, .15), output = "full"))
  for (a in result$analyses) {
    values <- unlist(c(a$env, a$rank_summary, a$interpolation), use.names = FALSE)
    expect_true(all(values >= 0 & values <= 1))
  }
  expect_true(all(unlist(result$analyses$intervals$indices$overall$p) <= 1))
})

test_that("the engine corrects roundoff but rejects invalid computed probabilities", {
  # Inject invalid internal allocations to exercise guards unreachable with valid input.
  d <- data.frame(low = c(1, 0), high = c(0, 1))
  spec <- spseg_continuous(vars = c("low", "high"), representation = "intervals",
    summary = list(probs = .5))
  config <- spseg:::.spseg_compile(spec, d, rep(1, 2), "all")
  config$settings <- list(scope = "both", output = "localenv",
    keep_scores = FALSE, keep_contributions = FALSE)
  run <- function(allocation) {
    input <- config
    input$group_fractions <- allocation
    spseg:::spseg_engine_cpp(cbind(c(0, 1), 0), list(v = input), 0, 3, 1L,
      TRUE, 0L, 1L, list())$analyses$v$env[[1]]
  }
  tiny <- diag(2) * (1 + .Machine$double.eps)
  expect_equal(unname(run(tiny)), diag(2), tolerance = 0)
  tiny_negative <- diag(2)
  tiny_negative[1, 2] <- -.Machine$double.eps
  expect_equal(unname(run(tiny_negative)), diag(2), tolerance = 0)
  expect_error(run(diag(2) * 1.000001), "outside.*roundoff allowance")
  negative <- diag(2)
  negative[1, 2] <- -1e-6
  expect_error(run(negative), "outside.*roundoff allowance")
  nonfinite <- diag(2)
  nonfinite[1, 1] <- Inf
  expect_error(run(nonfinite), "non-finite probability")
})
