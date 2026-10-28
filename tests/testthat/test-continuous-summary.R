summary_fixture <- function() {
  data.frame(x = c(0, 10), y = 0, a = c(2, 0), b = c(1, 1),
             c = c(0, 2), d = c(1, 1))
}

summary_run <- function(d, summary = NULL, ...) {
  spseg(d, list(income = spseg_continuous(counts = c("a", "b", "c", "d"),
    scores = c(10, 20, 30, 40), summary = summary)), coords = c("x", "y"),
    bands = 0, weighting = "unweighted", output = "full", ...)
}

test_that("continuous summaries default to synchronized quartile shares", {
  d <- summary_fixture()
  result <- summary_run(d)
  endpoints <- summary_run(d, list(probs = c(0, .25, .5, .75, 1)))
  expect_identical(spseg_continuous(values = "income")$summary$type, "shares")
  expect_equal(result$analyses, endpoints$analyses)
  expect_equal(result$analyses$income$env[[1]],
               as.matrix(d[c("a", "b", "c", "d")]) / 4, ignore_attr = TRUE)
  expect_equal(result$diagnostics$focal_searches, 2L)
  expect_true(all(as.data.frame(result, what = "exposure")$statistic == "share"))
})

test_that("local and grouped rank statistics use the same focal-weighted quantities", {
  result <- summary_run(summary_fixture(),
    list(type = "quantiles", probs = .5, mean_rank = TRUE))
  analysis <- result$analyses$income
  expect_equal(analysis$env[[1]], cbind(c(1/8, 5/8), c(3/8, 5/8)),
               ignore_attr = TRUE)
  expect_equal(analysis$rank_summary[[1]], cbind(c(1/4, 1/2), c(7/16, 9/16)),
               ignore_attr = TRUE)
  expect_null(analysis$indices$multigroup$p)
  out <- as.data.frame(result, what = "exposure")
  expect_setequal(out$statistic, c("quantile_rank", "mean_rank"))
  expect_true(all(is.na(out$to_group)))
  expect_true(all(is.na(out$exposure)))
  expect_true(all(is.na(out$prob[out$statistic == "mean_rank"])))
  expect_equal(out$value, out$rank)
  env <- as.data.frame(result, what = "env")
  expect_equal(env$rank[env$statistic == "mean_rank"], c(3/8, 5/8))
})

test_that("continuous pairwise groups and scalar calculations are independent of presentation", {
  d <- summary_fixture()
  shares <- summary_run(d, list(probs = .5))
  ranks <- summary_run(d, list(type = "quantiles", probs = .5))
  expect_equal(shares$analyses$income$indices$pairwise,
               ranks$analyses$income$indices$pairwise, tolerance = 1e-12)
  for (measure in c("h", "r", "d", "s")) {
    expect_equal(shares$analyses$income$indices$multigroup[[measure]],
                 ranks$analyses$income$indices$multigroup[[measure]])
  }
  expect_equal(unname(shares$analyses$income$indices$pairwise$d[[1]][1, 2]), .5)
  expect_equal(shares$analyses$income$indices$multigroup$p[[1]],
               matrix(c(5/8, 3/8, 3/8, 5/8), 2), ignore_attr = TRUE)
  other <- summary_run(d, list(type = "quantiles", probs = c(.2, .8)))
  for (measure in c("h", "r", "d", "s")) {
    expect_equal(shares$analyses$income$indices$multigroup[[measure]],
                 other$analyses$income$indices$multigroup[[measure]])
  }
  expect_equal(shares$analyses$income$normalizers, other$analyses$income$normalizers)
})

test_that("exact cuts preserve ties and reject ambiguous specifications", {
  result <- summary_run(summary_fixture(), list(cuts = 20))
  expect_equal(result$analyses$income$env[[1]],
               rbind(c(.75, .25), c(.25, .75)), ignore_attr = TRUE)
  expect_error(spseg_continuous(values = "v", summary = list(cuts = 20, probs = .5)),
               "only one")
  expect_error(spseg_continuous(values = "v",
    summary = list(type = "quantiles", cuts = 20)), "requires shares")
  for (p in list(c(.5, .25), c(.5, .5), NA_real_, c(-.1, .5))) {
    expect_error(spseg_continuous(values = "v", summary = list(probs = p)),
                 "increase strictly")
  }
})

test_that("interval allocation retains requested groups and is auditable", {
  d <- data.frame(x = 0:1, y = 0, a = c(10, 0), b = c(0, 10),
                  c = c(30, 40), d = c(10, 0))
  vars <- list(income = spseg_continuous(counts = c("a", "b", "c", "d"),
    representation = "intervals", summary = list(type = "quantiles")))
  expect_warning(result <- spseg(d, vars, coords = c("x", "y"), bands = c(0, 2),
                                 output = "full"), "Interval summary interpolation")
  definition <- result$analyses$income$definition$summary
  expect_equal(definition$requested_probs, c(.25, .5, .75))
  expect_equal(definition$effective_probs, c(.25, .5, .75))
  expect_equal(definition$groups$effective_groups, 4L)
  expect_equal(colSums(definition$groups$allocation * c(.1, .1, .7, .1)),
               rep(.25, 4), ignore_attr = TRUE)
  expect_equal(ncol(result$analyses$income$env[[1]]), 3L)
  expect_error(spseg(d, list(income = spseg_continuous(counts = c("a", "b"),
    representation = "intervals", summary = list(cuts = 20))),
    coords = c("x", "y"), bands = 0), "require exact")
})

test_that("single interval summaries explicitly use the rank interpolation model", {
  d <- data.frame(x = 0:1, y = 0, a = c(2, 4), b = 0)
  run <- function(summary) spseg(d, list(income = spseg_continuous(
    counts = c("a", "b"), representation = "intervals", summary = summary)),
    coords = c("x", "y"), bands = 0, output = "full")
  expect_warning(ranks <- run(list(type = "quantiles")), "Interval summary interpolation")
  expect_equal(unname(ranks$analyses$income$env[[1]]),
               matrix(rep(c(.25, .5, .75), each = 2), 2))
  expect_equal(nrow(as.data.frame(ranks, what = "env")), 6L)
  expect_equal(nrow(as.data.frame(ranks, what = "exposure")), 12L)
  expect_warning(mean <- run(list(type = "quantiles", mean_rank = TRUE)),
                 "Interval summary interpolation")
  expect_equal(unname(mean$analyses$income$env[[1]][, 4]), c(.5, .5))
  expect_warning(shares <- run(list(type = "shares")), "Interval summary interpolation")
  expect_equal(unname(shares$analyses$income$indices$multigroup$p[[1]]), matrix(.25, 4, 4))
})

test_that("open upper categories apportion deciles and report neighborhood diagnostics", {
  d <- data.frame(x = c(0, 10), y = 0, bottom = c(70, 100), top = c(30, 0))
  spec <- spseg_continuous(counts = c("bottom", "top"), representation = "intervals",
    summary = list(probs = seq(.1, .9, .1)))
  expect_warning(result <- spseg(d, list(income = spec), coords = c("x", "y"),
    bands = c(0, 20), weighting = "unweighted", output = "full"),
    "Interval summary interpolation")
  analysis <- result$analyses$income
  expect_equal(unname(analysis$definition$summary$groups$allocation[2, 9:10]),
               c(1/3, 2/3))
  expect_equal(unname(analysis$env[[1]][, 10]), c(.2, 0))
  expect_equal(unname(analysis$indices$multigroup$p[[1]][, 10]),
               c(rep(7/85, 8), 12/85, 1/5))
  expect_equal(unname(analysis$interpolation[[1]]),
               cbind(c(.7, 1), c(.3, 0), c(1, 1)))
  expect_equal(unname(analysis$interpolation[[2]]),
               cbind(c(.85, .85), c(.15, .15), c(1, 1)))
  env <- as.data.frame(result, what = "env", bands = 0)
  expect_equal(unique(env$top_category_share[env$location == 1]), .3)
  expect_true(all(env$interpolated_share == 1))
  mapped <- sf::st_as_sf(result, bands = 0)
  expect_equal(mapped$interpolated_share, env$interpolated_share)
  aligned <- spseg_continuous(counts = c("bottom", "top"), representation = "intervals",
    summary = list(probs = .85))
  expect_no_warning(unsplit <- spseg(d, list(income = aligned), coords = c("x", "y"),
    bands = c(0, 20), weighting = "unweighted", output = "full"))
  expect_true(all(unsplit$analyses$income$interpolation[[1]][, 3] == 0))
  for (measure in c("h", "r", "d", "s")) {
    expect_equal(analysis$indices$multigroup[[measure]],
                 unsplit$analyses$income$indices$multigroup[[measure]])
  }
  expect_equal(analysis$normalizers, unsplit$analyses$income$normalizers)
  expect_warning(compact <- spseg(d, list(income = spec), coords = c("x", "y"),
    bands = 0, weighting = "unweighted"), "Interval summary interpolation")
  expect_length(compact$analyses$income$interpolation, 0L)
  expect_null(compact$analyses$income$env)
})

test_that("removed continuous arguments fail instead of selecting an old path", {
  expect_error(spseg_continuous(values = "v", localenv = "quantiles"),
               "unused argument")
  expect_error(spseg_continuous(values = "v", exposure = list(probs = c(0, 1))),
               "unused argument")
  expect_error(spseg_continuous(values = "v", localenv = "shares",
                                summary = list(type = "shares")), "unused argument")
  expect_no_warning(spseg_ordinal(values = "v", localenv = "quantiles"))
})
