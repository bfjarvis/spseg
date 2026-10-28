# T012 — Independently validate continuous summary implementation

Status: in review (revised independent references, full suite, and constrained package check pass)
Owner: validation/testing. Reviewers: science officer and management.

Starting revision: `e2b8acac868cdca8d0119f79c5439557e6a41dff` in the canonical
checkout. Pre-existing changes include `AGENTS.md`, `README.md`, `NEWS.md`,
project records/reports, and the user's in-progress
`vignettes/spseg-getting-started.Rmd`; all are preserved. The only T012 code
artifact is `tests/testthat/test-continuous-summary-reference.R`.

Read T011, the design, and existing T008/T010 reference examples. Prepare independent
tidyverse/base-R reference calculations in temporary files while developer owns
package source/tests. Use private libraries and copied snapshots, never reinstall
into the active framework library. Preserve all user edits, particularly getting started.

Verify on the candidate: shared quartiles, optional endpoints, cuts conflicts,
exact/count equivalence, focal versus neighbor weighting, mean/median/pooled
distinctions, interval exact-match/nearest/tie/duplicate/fewer-interval/zero/single
cases, one adjustment warning, auditable metadata, synchronized local/grouped
outputs, binary pairwise values, unchanged continuous scalar normalizers and
neighborhood/search behavior. Use independent expectations, not package internals.
Test full versus grouped-only outputs, missingness, radius/kNN, and network cases.

Coordinate with developer by reading its incremental T011 report. When interface
is ready, write independent regression tests only to
`tests/testthat/test-continuous-summary-reference.R`; developer owns other test
files. Before then, edit only this report and temporary scripts. Report observed
failures immediately in this report/final reply for management to relay. Do not
repair package source or change user's vignette. No commits/pushes.

Run appropriate existing regression checks on the final source snapshot. The prior
network unreachable-population SIGABRT is a known unresolved blocker: investigate
that failure if it persists, separating source defects from compiler/runtime issues,
and provide precise evidence to developer/management. Never skip it silently or
declare a passing suite after a process abort. Record exact tested revision/diff,
commands, warnings, failures, and checks not run. Mark in review with an acceptance
recommendation only when the candidate evidence supports it.

## Results

Pending.

## Incremental results (2026-10-28)

### Independent summary reference cases prepared

The approved T010 eight-person/two-unit fixture is encoded independently in
`/private/tmp/spseg-T012-reference.R`; it loads no package code. The script
asserts hand-derived expectations for local shares, local mean ranks, inverse
empirical-CDF medians, group averages weighted by original focal group counts,
and binary low/high dissimilarity. Values are:

| Quantity | Expected values |
| --- | --- |
| U/V local low-high shares | `(.75, .25)` / `(.25, .75)` |
| U/V local mean rank | `.375` / `.625` |
| U/V local median rank | `.125` / `.625` |
| Low/high focal-group exposure rows | `(.625, .375)` / `(.375, .625)` |
| Low/high focal-group mean rank | `.4375` / `.5625` |
| Low/high focal-group average local median | `.25` / `.50` |
| Binary low/high dissimilarity | `.5` |

The same script checks interval-grid expectations for exact quartile boundaries,
insufficient-interval fallback, duplicate nearest-boundary matches, lower-boundary
tie-breaking, and a single positive interval. `Rscript --vanilla
/private/tmp/spseg-T012-reference.R` passes. These are reference expectations,
not yet comparisons with the new package implementation.

### Independent resolution of the prior network abort

The T002 failure was reproducible in the local configured GCC 16 toolchain. The
isolated source copies were made from the `e2b8acac868cdca8d0119f79c5439557e6a41dff`
worktree before T011's later in-progress package edits. I rebuilt one isolated copy at `/private/tmp/spseg-T012-gcc-plain/source` with GCC
16.2.0, overriding the user's `~/.R/Makevars` to remove `-march=native` and
`-fopenmp` while retaining C++17 and ordinary `-O2 -g -fPIC` flags. The package
installed successfully, but the exact minimal unreachable-population case and
whole `test-network.R` still exited 134 (SIGABRT). Thus those extra optimization
and OpenMP flags do not explain this result.

The same source was then built in
`/private/tmp/spseg-T012-clang/source` with Apple Clang 17.0.0.13.5, C++17, and
`-O2 -g -fPIC`, with the installed Homebrew Fortran runtime supplied for linking
in the temporary Makevars. The unreachable two-node directed kNN fixture now
returns the expected informative R error and exits normally; reachable and
infeasible-global-threshold controls also survived with their expected outcomes.
The complete `test-network.R` passed. I then ran `testthat::test_dir()` on all
seven current test files in one fresh process; all existing test files passed,
including the network test. Log: `/private/tmp/spseg-T012-clang/all-tests.log`.

Both T002's GCC16 build and this T012 source snapshot have identical hashes for
all 27 tracked files under `R/`, `src/`, and `tests/`:
`67f4028ad70c8782bcd55c2c39ba5dab946b4be480919e24ae551f9b354f64c1` (computed
from sorted path plus file SHA-256 entries). The test itself was not weakened or
changed. The comparison makes the local GCC16/C++ runtime configuration the
leading explanation; it does not rule out source undefined behavior exposed by
that compiler. Rcpp 1.1.2's `Rcpp::stop()` implementation throws an
`Rcpp::exception`, and the generated Rcpp wrapper catches it to create an R
condition, so the earlier generic long-jump hypothesis is not supported by this
inspection. A debugger backtrace was not obtained: LLDB could not launch the
process in this environment.

The Clang build used a private makevars file to avoid the host's user GCC16
override and to provide the locally installed Fortran runtime:

```make
CC = /usr/bin/clang
CXX = /usr/bin/clang++
CXX11 = /usr/bin/clang++
CXX14 = /usr/bin/clang++
CXX17 = /usr/bin/clang++
CXX17STD = -std=gnu++17
CFLAGS = -O2 -fPIC
CXXFLAGS = -O2 -g -fPIC
CXX11FLAGS = -O2 -g -fPIC
CXX14FLAGS = -O2 -g -fPIC
CXX17FLAGS = -O2 -g -fPIC
FLIBS = -L/opt/homebrew/opt/gfortran/lib/gcc/current/gcc/aarch64-apple-darwin24/16 -L/opt/homebrew/opt/gfortran/lib/gcc/current -lemutls_w -lheapt_w -lgfortran -lquadmath
```

Reproduction used:

```sh
R_MAKEVARS_USER=/private/tmp/spseg-T012-clang/Makevars \
  R CMD INSTALL --preclean \
  --library=/private/tmp/spseg-T012-clang/library \
  /private/tmp/spseg-T012-clang/source
Rscript --vanilla \
  /private/tmp/spseg-T012-run-all-tests.R \
  /private/tmp/spseg-T012-clang/library \
  /private/tmp/spseg-T012-clang/source/tests/testthat
```

The GCC comparison used the same install command with
`/private/tmp/spseg-T012-gcc-plain/Makevars`, its own source copy and library.
Focused outcomes are under each temp root (`minimal-*.log`, `network.log`, and
`install.log`). Environment: macOS Sequoia 15.6.1 arm64, R 4.6.0, Rcpp 1.1.2,
RcppArmadillo 15.6.2.1, BH 1.90.0.1, sf 1.1.3, terra 1.9.50, testthat 3.3.2.
No package source or test file was edited. No vignettes, `R CMD check`, Windows,
Linux, alternate R version, or external network workflows were tested. The
successful makevars path is specific to this host and is evidence for a compiler
comparison, not a general installation recipe.

### Executable summary candidate validation

T011 recorded an executable core-source handoff on 2026-10-28. The isolated
candidate snapshot at `/private/tmp/spseg-T012-candidate/source` included
untracked R source and developer regression tests, plus current tracked files;
the manifest is `/private/tmp/spseg-T012-candidate/manifest.json` (snapshot
digest `4c12bffc1e46f1b485f65eb25e494d482335c22038c60993662470b6c8a3bb08`).
It was built and installed with the private Apple Clang Makevars described
above. The build/install succeeded. No vignette was rendered or edited.

Added `tests/testthat/test-continuous-summary-reference.R`, based on independent
hand arithmetic in `/private/tmp/spseg-T012-reference.R`. It now checks local
shares, endpoint-delimited grid equivalence, tied-value median midranks,
focal-weighted grouped rank summaries, typed rank extraction, binary pairwise
dissimilarity, exact cuts, interval boundary resolution, scalar invariance,
individual/count equivalence, grouped-only extraction, and radius/kNN search
agreement.

The first complete snapshot test exposed post-snapshot changes in the current
R source. I therefore copied the current complete non-project source tree
(including untracked files and the user vignette) to
`/private/tmp/spseg-T012-current/source`, rebuilt and installed it with the same
private Apple Clang Makevars, and reran only the independent file. Snapshot
digest: `a7e7b7eb1103f3a80d314b817a129bb5f0b4d3c5ae52141d50808bab760d9ac5`
(72 files). Build/install passed; the focused independent test file passed
(43 expectations; `/private/tmp/spseg-T012-current/reference-test.log`). The
test checks:

- interval exact matching; a requested grid larger than the positive-interval
  count; duplicate nearest-boundary matches; lower-boundary tie resolution; and
  a single positive interval, including one-warning behavior and effective-grid
  metadata;
- unchanged continuous scalar measures and normalizers across summary grids,
  with shared focal-search count;
- equality of local/grouped rank summaries and segregation indices for
  individual and count representations of the same population;
- equality of requested grouped profiles under `measures = "exposure"` and
  full-measure calculation; and
- radius and kNN results agreeing with the brute-search path across multiple
  bands.

The earlier full-suite run was on a previous candidate snapshot and failed in
`test-local-environments.R` at its 10-failure cap. It showed deprecated legacy
argument warnings and stale default-quantile/output expectations. That result is
superseded by the complete current-source test run below.

### Previous-candidate complete-source verification (2026-10-28)

After the developer confirmed the candidate and regression migrations were ready,
I copied the current complete non-project source tree, including untracked
sources/tests and the in-progress getting-started vignette, into
`/private/tmp/spseg-T012-final/source`. Snapshot digest:
`7018efc7a4f1b79e84327917943aab1e47cd9ccfa3c65dad389c08f467f15986`
(72 files). Compiler: Apple Clang 17.0.0.13.5; R 4.6.0; macOS Sequoia 15.6.1
arm64.

`R CMD build --no-build-vignettes --no-manual`, private-library install, and the
complete `testthat::test_dir()` all passed. The full run included the independent
reference file and network tests. Log:
`/private/tmp/spseg-T012-final/tests.log`.

The first `R CMD check --no-manual --no-build-vignettes` stopped before package
checks because suggested package `sfnetworks` is unavailable. With force-suggests
disabled, check ran package tests successfully but failed its vignette code
check: the user-owned `spseg-getting-started.Rmd` has a duplicate chunk label
`micro-indices`. The exact diagnostic is in
`/private/tmp/spseg-T012-final/spseg.Rcheck/spseg-getting-started.Rmd.log`.
Other vignette chunks checked successfully. No vignette was rendered, edited, or
regenerated by this task. A final
`_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --ignore-vignettes`
completed with exit 0 and `Status: 1 NOTE`; package tests passed and vignette
checks were explicitly skipped. The remaining NOTE is that
`.spseg_resolve_summary` uses unqualified `head` and `tail`, which R CMD check
reports as undefined globals and recommends importing from `utils`.
`sfnetworks` remains unavailable. Check log:
`/private/tmp/spseg-T012-final/check-ignore-vignettes.log`.

This passing result was for the prior interval-boundary matching/fallback
contract. It remains historical evidence for the other summary contracts but
does not establish acceptance of the revised interval allocation behavior.

### Revised interval interpolation contract — final verification

Management superseded nearest-boundary/fallback matching with interpolation and
proportional population allocation so requested quantiles remain available.
The previous candidate tests and verdict above are historical; the revised
expectations below were tested against the final candidate.

Prepared independent reference cases use global interval mass `(0.4, 0.2, 0.4,
0)`, with requested rank bins at `.25, .5, .75`. The hand-derived allocation is
`(5/8, 3/8, 0, 0)`, `(0, 1/2, 1/2, 0)`, `(0, 0, 3/8, 5/8)`, and all zeros for
the unpopulated interval. Every positive source interval conserves its mass;
each target rank bin receives `.25` of global population. Independently derived
local bin shares and focal-group exposure rows are asserted in
`tests/testthat/test-continuous-summary-reference.R`. Additional cases check
two positive intervals separated by zero-mass intervals, a single positive
interval spanning all four requested rank bins, one warning per adjusted
variable, retention of requested probabilities, and invariance of full
continuous scalar indices across summary grids.

The following final run records that evidence.

After the developer declared the revised candidate ready, I copied the complete
current non-project source tree, including untracked source and test files, into
`/private/tmp/spseg-T012-revised-final/source`. The user's unfinished
getting-started vignette was copied byte-for-byte and never edited, rendered, or
regenerated. The 72-file manifest is
`/private/tmp/spseg-T012-revised-final/manifest.json`, with aggregate SHA-256
`c50b5c0d666dd95a0b5b18b4f6c65b720c46007bc5e58e1d623bb610b86824dd`. A final
hash comparison found no differences between the snapshot and canonical files.
Key hashes: `R/spseg_summary.R` =
`4eace09fb81616189c297e3545d000d1aa313d60e9f4bbe1f1621f6cbdc5a0a7`,
`src/spseg_engine.cpp` =
`7679323159c1f91109dee1ad3e67081f3ae9d8948d589b3a83bebd822b7845b0`, and
`tests/testthat/test-continuous-summary-reference.R` =
`f54d3fb25e3acdca2e1449696e5953b9b789570967f8de9b05fefbf203d8a8ba`.

On macOS Sequoia 15.6.1 arm64, R 4.6.0, and Apple Clang 17.0.0.13.5:

- `R CMD build --no-build-vignettes --no-manual` and private-library install
  passed.
- The revised independent file passed all 48 expectations. It checks exact
  proportional allocation, per-positive-interval conservation, zero-mass
  intervals, preservation of requested probabilities, local/focal target
  shares, single-interval allocation, and full-resolution scalar invariance.
  Its fractional interval pairwise dissimilarity reference is `0.4`; with one
  populated interval, the corresponding binary pairwise dissimilarity is `0`.
  Log: `/private/tmp/spseg-T012-revised-final/independent-final.log`.
- The complete test directory passed, including all nine test files and network
  regressions. Log:
  `/private/tmp/spseg-T012-revised-final/full-tests-final.log`.
- `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --ignore-vignettes`
  completed with `Status: OK`. R CMD check reported `sfnetworks` unavailable
  for checking, but no NOTE/WARNING/ERROR. Vignette checks were explicitly
  skipped. Log: `/private/tmp/spseg-T012-revised-final/check-final.log`.

**Validation verdict:** recommend acceptance of the revised implementation test
evidence for T012. The constrained package check is clean, with `sfnetworks`
unavailable and vignette checks skipped as recorded limits. This is not release
certification or scientific validation of every published measure. The earlier
local GCC16 native exception abort remains a host/toolchain limitation; this
final suite was built and run under the private Apple Clang configuration. No
source or developer-owned test was changed by validation. No commit/push was
made; management owns that step.
