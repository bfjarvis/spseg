# T003 — Implementation assessment handoff

Status: in review
Owner: developer. Reviewers: validation/testing and management.
Starting revision / worktree: `de48d12a3e9d29ca09eb289873f1877ba5582117`,
`codex/relational-neighborhoods`, `/Users/bjarvis/Documents/Develop/R/spseg`.

## Brief

- Read current `AGENTS.md`, `project/GOALS.md`, and `project/TEAM.md` in `spseg`.
  Explicitly use the renamed checkout if the chat retains the old `seg` directory.
- Reuse your recent codebase scan; do not repeat the full scan or test suite.
- Record the implementation map, public-interface limitations, and important gaps
  against the documented geographic/network and multi-variable workflows.
- Prioritize a concise handoff of the network-test abort: tested revision, copied
  sources, build/test commands, environment, observed crash, and evidence paths.
  Distinguish observations from hypotheses; validation owns independent reproduction.
- Note relevant interoperability concerns for tidyverse, sf, sfnetworks edge-list
  preparation, and data.table without claiming support that was not checked.
- Edit only this report. No package fixes, generated-file changes, or vignette edits;
  continue leaving research manuscripts outside this assignment.

## Results

Prepared 2026-10-26 from the earlier scan and test evidence. Only this report was
edited. No package fixes, new test runs, broad rescan, or manuscript review.

### Provenance

- Assessed/tested package revision: `de48d12a3e9d29ca09eb289873f1877ba5582117`.
  At the scan, unrelated edits were `.Rbuildignore` and the getting-started
  vignette, plus untracked papers. They were excluded from the test copy.
- Report completion HEAD: `e2b8acac868cdca8d0119f79c5439557e6a41dff`, same branch.
  `git diff --name-status de48d12a3e9d29ca09eb289873f1877ba5582117 HEAD -- DESCRIPTION NAMESPACE R src tests data inst man`
  returned no changes. At resumption, unrelated dirty files were `project/PLAN.md`
  and `project/PROGRESS.md`; concurrent specialist report edits appeared during
  completion. All were preserved.
- The earlier temporary source copy remains identified by
  `/private/var/folders/t5/w6p00r9118l6rtbthpcvyf180000gp/T/RtmpK6NNyy/spseg-codebase-scan-6bc85094ec9c`.
  Before the interruption, an MD5 comparison confirmed all 43 non-binary files
  under the copied entries below matched the canonical checkout. Temporary
  artifacts can disappear; recreate from the assessed revision if necessary.

### Implementation map

| Area | Implementation and responsibilities |
| --- | --- |
| Public API and orchestration | `R/spseg.R` and `R/spseg-internal.R`: required named variable specifications, validation, shared surface preparation, native call, S3 result assembly |
| Variable definitions | `R/spseg_variables.R`: categorical/ordinal/continuous constructors, eligibility, exact-value weighted midranks, interval CDFs, exposure cuts, retained summaries |
| Measurement | `src/spseg_engine.cpp`: one search per focal location, reused across bands/variables; local distributions, index accumulators, optional compact environment storage |
| Selection and kernels | `src/spseg_neighbors.h`: radius/population-kNN contributions, fractional focal inclusion, proportional tied boundaries, distance kernels |
| Search | `src/spseg_euclidean.h`: Boost.Geometry R-tree/brute force; `src/spseg_graph.cpp`: sparse Boost.Graph BFS/Dijkstra and reusable query buffers |
| Surfaces | `R/spseg_surface.R`: terra centroid rasterization and fallback counts; `src/spseg_surface.cpp`: groupwise averaging and effective-zone total correction |
| Results | `R/spseg_result.R`: S3 print/list/data-frame/sf extraction, variable/band selection, long-form share or global-rank-quantile summaries |

Index averaging uses each variable's original eligible focal population.
Categorical formulas use local shares; ordinal formulas sum variation across
cumulative cutoffs; continuous formulas integrate CDF variation. Exact continuous
values use an empirical weighted-midrank convention; intervals use linear CDF
interpolation in global percentile space. These are distinct interpretations.

Candidate discovery is shared across variables and bands; kNN contributions are
also shared when eligible-population vectors coincide. kNN selects population
before kernel weighting, including fractional support units and tied boundary
shells. Pairwise scope does not select a new neighborhood for each pair.
Nanoflann/kd-tree is no longer used. Measurement remains single-threaded.

Pycno receives fallback contributions before smoothing. Each iteration averages
the cell and its valid eight neighbors, then corrects group totals within
effective zones. The stopping threshold is maximum cell change below
`max(initial) * 10^(-converge)`, subject to `max_iter`. Fallback polygons need not
remain separate constraints. Preserve the warnings and optional surface audits.

Roxygen comments are authoritative; regenerate help/NAMESPACE rather than editing
generated files. Regenerate native wrappers with `Rcpp::compileAttributes()`.

### Network crash: reproduction and observations

The earlier checks used separate R processes and a temporary source copy. Installed
dependency libraries were reused, not a private dependency library. Canonical
compiled artifacts and the installed spseg package were not modified. Equivalent
source-copy setup and the original test command:

```r
source_dir <- "/Users/bjarvis/Documents/Develop/R/spseg"
scan_dir <- tempfile("spseg-codebase-scan-")
dir.create(scan_dir)
entries <- c("DESCRIPTION", "NAMESPACE", "R", "src", "tests", "data", "inst", "man")
stopifnot(all(file.copy(file.path(source_dir, entries), scan_dir, recursive = TRUE)))
unlink(list.files(file.path(scan_dir, "src"), pattern = "[.](o|so|dll)$",
                  full.names = TRUE))
devtools::test(scan_dir, reporter = "summary", stop_on_failure = TRUE)
```

The temporary copy compiled; the suite reached network tests and exited 134,
without a normal test summary. Two fresh-process network-only follow-ups also
exited 134. The last explicitly attached Rcpp:

```r
library(Rcpp)
pkgload::load_all(scan_dir, recompile = FALSE, helpers = FALSE)
testthat::test_file(file.path(scan_dir, "tests/testthat/test-network.R"),
                    reporter = "progress")
```

Loading completed. The last displayed progress was 145 passing network
expectations, with no ordinary assertion failure. A separate fresh-process run
after `pkgload::load_all(scan_dir, recompile = FALSE)` passed all six non-network
test files with exit 0:

```r
filter <- paste0("^(boost-neighbors|core-regression|default-bands|",
                 "local-environments|multiple-variables|unified-interface)$")
testthat::test_local(scan_dir, filter = filter, reporter = "summary",
                     stop_on_failure = TRUE)
```

Local crash reports are under `/Users/bjarvis/Library/Logs/DiagnosticReports/`:
`R-2026-10-26-104350.ips` (full suite), `R-2026-10-26-104409.ips` (network-only),
and `R-2026-10-26-104431.ips` (explicit Rcpp attachment). The last records
`EXC_CRASH`/`SIGABRT` with `free_tiny_botch` and `malloc_zone_error` above
`exception_to_condition_template<Rcpp::exception>` (`exceptions.h:366`),
`rcpp_exception_to_r_condition` (`exceptions.h:370`), and
`_spseg_spseg_engine_cpp` (`RcppExports.cpp:33`). No standalone console/build log
was saved. Inspect raw crash reports for machine-specific details before sharing.

The unreachable-population error test at `tests/testthat/test-network.R:148` was
a candidate to isolate, not an established minimal reproduction from this scan.
Its first call requests two people at an isolated node containing one. The exact
failing expression was not identified here. The current plan separately reports
Validation's isolation to that test block; Validation owns that independent
evidence and its durable report.

The stack establishes a failure during exception conversion, not its cause.
Earlier corruption, exception/runtime incompatibility, and build configuration
remain hypotheses; do not attribute the failure to graph traversal or Rcpp alone.

### Recorded environment

| Component | Version/configuration |
| --- | --- |
| R | 4.6.0 (2026-04-24), `aarch64-apple-darwin23` |
| Rcpp / RcppArmadillo | 1.1.2 / 15.6.2.1 |
| BH / sf / terra | 1.90.0.1 / 1.1.3 / 1.9.50 |
| testthat / devtools / pkgload | 3.3.2 / 2.5.2 / 1.5.2 |
| Configured CXX17 | `/opt/homebrew/bin/g++-16`, Homebrew GCC 16.2.0 |
| Configured CXX17FLAGS | `-falign-functions=64 -Wall -g -O2 -O3 -march=native -fPIC -fopenmp` |

`library(Rcpp)` warned it was built under R 4.6.1. Compiler settings came from
`R CMD config`, not a retained build transcript; development-build overrides
were not recorded. `otool -L` on the snapshot's `src/spseg.so` confirmed R 4.6
BLAS/LAPACK/libR and Homebrew GCC `libstdc++`/`libgomp` linkage. These observations
are inputs for controlled comparison, not proof of a toolchain defect or internal
parallel execution.

### Interface and interoperability limits

- Named `variables` are mandatory, even for one outcome. Categorical scope
  controls multigroup/pairwise scalars; ordered summaries and directional exposure
  have different meanings. Stored environments are shares or rank quantiles,
  not local count matrices. Indices use full distributions, independently of
  retained presentation settings.
- Default indices-only output omits input/support/environment storage. Environment
  extraction is long-form only; historical `shape` arguments are rejected.
  Data frames identify rows by `location` and optional network `node`, without
  automatically binding coordinates or arbitrary source identifiers.
  `st_as_sf()` attaches geometry or reconstructs coordinate points.
- Network workflows require explicit bands, `search = "boost"`, and raw support.
  Users prepare topology, adjacency/snapping, costs, node IDs, and population
  allocation externally. Geometry does not define relationships. The sfnetworks
  vignette demonstrates an edge-list handoff, but its optional preparation chunk
  was not executed in this scan.
- Grid/pycno requires polygon sf input and count specifications for every variable;
  observed-value specifications cannot be gridded directly. Geometry and full
  audits are optional to control memory. Radius zero combines coincident points.
- sf and mixed-variable extraction have regression evidence. Tibble,
  grouped-tibble, and data.table workflows were not independently checked; grouping,
  row identity, and absence of caller mutation must not be advertised as validated.
  Check these plus nonstandard geometry columns when changing input/extraction.
- `DESCRIPTION` declares Rcpp/RcppArmadillo/BH in `LinkingTo`, but not Rcpp in
  `Imports`. Namespace/loading requirements merit packaging review; this is not
  evidence that the declaration caused the abort.

### Next actions and limitations

1. Validation supplies its minimal network error-path reproduction and an isolated
   source/private-library baseline; compare matched R/Rcpp/compiler configurations
   before assigning a cause. Development then implements the agreed targeted repair
   and regression test, followed by independent network-suite verification.
2. Support the analyst-owned getting-started vignette and subsequent small income
   and network examples in the order specified by the current plan. Preserve the
   user's narrative; route scientific interpretation to Science and independent
   expected results to Validation.
3. Keep optimization and new methods outside this report. No new suite, benchmark,
   vignette render, `R CMD check`, Windows check, or scientific-source audit was
   performed for T003. Earlier non-network success is a local baseline, not release
   readiness or cross-platform/scientific validation.
