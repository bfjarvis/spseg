# Run with Rscript tools/benchmark-multiple.R after installing seg, or load
# the development package before sourcing. No external data are required.
library(seg)
set.seed(881)
n <- as.integer(Sys.getenv("SEG_BENCHMARK_N", "5000"))
d <- data.frame(x=runif(n), y=runif(n), ethnicity=sample(letters[1:4], n, TRUE),
                education=ordered(sample(1:4, n, TRUE)),
                age=sample(18:90, n, TRUE), income=round(exp(rnorm(n, 10, 1))))
variables <- list(ethnicity=seg_categorical(values="ethnicity"),
                   education=seg_ordinal(values="education"),
                   age=seg_continuous(values="age"),
                   income=seg_continuous(values="income"))
bands <- c(.025,.05,.1) * sqrt(5000/n)
measure <- function() {
  joint <- system.time(result <- spseg(d, coords=c("x","y"), variables=variables,
    bands=bands, scope="multigroup"))["elapsed"]
  separate <- system.time(invisible(lapply(seq_along(variables), function(i)
    spseg(d, coords=c("x","y"), variables=variables[i], bands=bands,
           scope="multigroup"))))["elapsed"]
  stopifnot(result$diagnostics$focal_searches==n)
  c(joint=unname(joint), separate=unname(separate), speedup=unname(separate/joint))
}
print(replicate(3, measure()))
