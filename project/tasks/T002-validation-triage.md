# T002 — Validation inventory and network-test triage

Status: in review (2026-10-26)
Owner: validation/testing. Reviewers: developer and management.
Starting revision / worktree: `de48d12a3e9d29ca09eb289873f1877ba5582117`,
canonical checkout. The assessment snapshot also included tracked working-tree
contents for 58 files; pre-existing edits were `.Rbuildignore`, `.gitignore`, and
`vignettes/spseg-getting-started.Rmd`, with untracked `AGENTS.md` and `project/`.
The temporary package source was copied from tracked package files after those
edits; no package or test source was edited for this assessment.

## Results

### Network abort reproduced

I copied the working package sources to `/private/tmp/spseg-T002-er5dx_ot/source`,
built and installed them to the private library at
`/private/tmp/spseg-T002-er5dx_ot/library`, then ran tests in fresh `Rscript
--vanilla` processes. `R CMD INSTALL --preclean --library=/private/tmp/spseg-T002-er5dx_ot/library /private/tmp/spseg-T002-er5dx_ot/source`
completed successfully. Environment: macOS Sequoia 15.6.1, arm64, R 4.6.0,
Rcpp 1.1.2, RcppArmadillo 15.6.2.1, BH 1.90.0.1, sf 1.1.3, terra 1.9.50,
testthat 3.3.2; package compilation used configured Homebrew GCC 16 and
`-march=native -fopenmp` flags. `sfnetworks` is unavailable.

The complete `test-network.R` run exits with code 134 (SIGABRT). Running each of
its eight top-level `test_that()` expressions separately isolates the abort to
the sixth test, “unreachable eligible populations produce informative errors.”
The other seven network test blocks pass. The failure reduces to a two-node
directed graph with one reachable and one unreachable node, a categorical
variable, and a kNN target of 2. The unreachable focal node requires an R error;
that process aborts with code 134. The same fixture with an undirected, reachable
path or with a larger global kNN target survives and reports/catches a regular R
error. The source path throws `SearchComplete` to stop graph traversal and catches
that type in `GraphSearch::Impl::find`; later, the engine calls Rcpp `stop()` when
the eligible population cannot meet the target. This narrows the observed failure
to the unreachable-population error path, but does not prove whether exception
handling, Rcpp long-jump behavior, or another native interaction causes SIGABRT.
No package fix was attempted.

Reproduction commands (from a disposable copy):

```sh
R CMD INSTALL --preclean --library=/private/tmp/spseg-T002-er5dx_ot/library /private/tmp/spseg-T002-er5dx_ot/source
Rscript --vanilla /private/tmp/spseg-T002-er5dx_ot/run-tests.R file test-network.R
Rscript --vanilla /private/tmp/spseg-T002-er5dx_ot/minimal.R unreachable
```

The test harness scripts, manifest, exact status records, and logs are under
`/private/tmp/spseg-T002-er5dx_ot/`; in particular, `install.log`,
`test-status.json`, `network-batch.log`, `block-1.log` through `block-8.log`, and
`minimal-unreachable.log`, `minimal-reachable.log`, `minimal-global.log`.
Running the test file in a single process gives exit 134 after the test reporter
prints its dots; per-test isolation identifies the sixth test. The other six
existing test files (56 test blocks total) pass in fresh processes. The first
harness attempt had an output-serialization bug after passing tests; it was fixed
in the temporary harness and those files were rerun successfully. The full
network test file was rerun after harness repair and still aborted.

### Evidence inventory and limits

The current suite contains seven test files and 64 `test_that()` blocks. Useful
independent correctness evidence includes:

- `helper-index-reference.R` calculates categorical indices from stored local
  compositions using separate formulas; the multiple-variable tests include
  numerical integration references for ordinal, grouped, and empirical
  continuous measures.
- Network tests include a small Floyd-Warshall distance reference and a separate
  R implementation of radius/kNN neighborhood composition. They compare directed
  and undirected, weighted and unweighted paths, several kernels, ties, zero-cost
  edges, transit nodes, parallel edges, shared observation nodes, and edge/row
  ordering.
- Coordinate neighbor tests compare Boost search with brute-force search, and
  regression tests cover output modes, local composition invariants, surface
  count conservation, and input validation.

These assertions are evidence for the listed fixtures and implemented conventions,
not proof of the full Reardon measure definitions or general network validity.
The network path-versus-Euclidean test intentionally checks equivalent fixtures;
the independent path/composition reference is stronger evidence for traversal and
composition, but index results remain supported by the common measure reference
suite. No release-level scientific validation is inferred from passing tests.

Material unverified areas for this assessment:

- The native abort blocks normal completion of the network file, so the observed
  test run cannot establish a clean full network-suite pass.
- Only macOS arm64 with one R/compiler/dependency stack was exercised. No Linux,
  Windows, alternate compiler, or older supported R check was run; compiler flags
  include host-specific `-march=native`.
- `sfnetworks` is not installed. Its vignette handoff and edge-list extraction
  were not executed. The suite does not establish broad tibble, grouped-data, or
  data.table interoperability; no dedicated data.table tests were found.
- No clean-session vignette render or `R CMD check` was run. Those checks are
  outside this focused reproduction, and the getting-started vignette has
  pre-existing in-progress edits.
- Tests exercise many edge and missingness conventions, but do not by themselves
  establish all population eligibility, geographic units/CRS, neighborhood
  weighting, income/ordinal scientific assumptions, or proposed network
  generalizations. Cross-platform and realistic road-network evidence remains
  necessary before making broad portability or network-validity claims.

### Handoff / next action

Developer should repair the unreachable eligible-population error path and add a
regression that asserts an informative R error without terminating the process.
Then rerun the isolated minimal case and complete network file in fresh processes;
record the repaired revision and compiler/runtime. Management should retain
portability, `sfnetworks` handoff, clean-session vignette rendering, and release
`R CMD check` as follow-up validation items. This report records the baseline
failure and validation evidence; it does not accept a repair or validate release
readiness.
