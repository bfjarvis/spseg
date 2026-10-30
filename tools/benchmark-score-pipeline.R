# Run in separate fresh R sessions for baseline and candidate libraries:
# Rscript tools/benchmark-score-pipeline.R /path/to/library /tmp/results.rds 20
# Compare median elapsed times and the saved numerical results across runs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply a package library and an output RDS path")
.libPaths(c(normalizePath(args[1]), .libPaths()))
library(spseg)
reps <- if (length(args) > 2L) as.integer(args[3]) else 20L
stopifnot(length(reps) == 1L, !is.na(reps), reps > 0)

set.seed(20261030)
n <- 3000L
columns <- paste0("g", 1:6)
counts <- matrix(sample(0:20, n * 6, replace = TRUE), n, 6,
  dimnames = list(NULL, columns))
d <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), counts,
  income = exp(rnorm(n)), node = seq_len(n))
population <- rowSums(counts)
category <- list(group = spseg_categorical(vars = columns))
mixed <- list(
  group = spseg_categorical(vars = columns),
  education = spseg_ordinal(vars = columns),
  income = spseg_continuous(vars = "income"),
  intervals = spseg_continuous(vars = columns, representation = "intervals",
    summary = list(probs = head(cumsum(colSums(counts) / sum(counts)), -1)))
)
make_case <- function(data = d, variables = category, bands = c(30, 75, 150),
                      output = "indices", ...) {
  c(list(x = data, variables = variables, coords = c("x", "y"),
    population = rowSums(data[columns]), bands = bands, output = output), list(...))
}
cases <- list(
  categorical_zero = make_case(bands = 0),
  categorical_radius = make_case(),
  mixed_radius = make_case(variables = mixed),
  mixed_knn = make_case(variables = mixed, neighbors = "knn", bands = c(100, 500, 1500)),
  mixed_full = make_case(variables = mixed, output = "full"),
  dense_radius = make_case(data = d[1:1000, ], bands = c(300, 750, 1500)),
  network = make_case(data = d[1:1000, ], bands = c(1, 5, 15), node_id = "node",
    network = spseg_network(data.frame(from = 1:999, to = 2:1000)))
)

data(spseg_data, package = "spseg")
xy <- as.matrix(expand.grid(x = seq(.5, 9.5), y = seq(.5, 9.5)))
cases$toy <- list(x = xy, data = spseg_data,
  variables = list(group = spseg_categorical(vars = names(spseg_data)[5:6])),
  bands = c(0, 1, 2, 4), output = "full")
auck <- sf::st_read(system.file("extdata", "auck", "auck.shp", package = "spseg"),
  quiet = TRUE)
auck$Korean <- auck$KR06
auck$Non_Korean <- auck$POP06 - auck$KR06
cases$auckland <- list(x = auck,
  variables = list(group = spseg_categorical(vars = c("Korean", "Non_Korean"))),
  bands = seq(500, 10000, 500), output = "full")

results <- timings <- list()
for (name in names(cases)) {
  call <- cases[[name]]
  result <- do.call(spseg, call)
  # Geometry and input tables are unchanged; retain numerical engine outputs.
  results[[name]] <- lapply(result$analyses, function(a)
    a[c("indices", "env", "normalizers", "rank_summary", "interpolation", "status", "knn")])
  batches <- switch(name, toy = 100L, categorical_zero = 20L,
                    network = 5L, 1L)
  invisible(do.call(spseg, call))
  elapsed <- numeric(reps)
  for (r in seq_len(reps)) {
    gc()
    elapsed[r] <- system.time(for (j in seq_len(batches))
      invisible(do.call(spseg, call)))[["elapsed"]] / batches
  }
  timings[[name]] <- elapsed
  cat(name, sprintf("median %.6f seconds\n", median(elapsed)))
}
saveRDS(list(timings = timings, results = results, session = sessionInfo()), args[2])
