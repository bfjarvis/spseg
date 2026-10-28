# T011 — Implement synchronized continuous summaries

Status: in review (user-directed interpolation revision)
Owner: developer. Reviewers: science officer, validation/testing, analyst, management.

## Authorized outcome

The user approved implementation and a final management commit of all outstanding
non-ignored changes after useful tests/reviews pass. The user is actively editing
`vignettes/spseg-getting-started.Rmd`: do not edit, regenerate, or restore it.
Management handles staging/commit; workers do not commit or push.

Read `project/design/continuous-summaries.md` and T008–T010. Implement one `summary`
specification replacing separate continuous localenv/exposure configuration:

- Default summary type: shares, following the management recommendation now adopted
  for implementation. Quantiles remain explicit. Default probabilities in either
  type: interior quartiles `.25, .5, .75`. One vector governs grouping and local
  evaluations. Accept optional 0/1 delimiters, remove them for local evaluation.
- Exact `cuts` are mutually exclusive with explicit `probs`, with implicit infinite
  tails. Initially allow cuts in shares mode; reject quantiles/cuts with useful
  guidance. Keep ties intact; metadata distinguishes requested/achieved groups.
- Interval probabilities resolve once against original eligible global interval
  totals, before surface redistribution. Nearest interior cumulative boundaries,
  lower-boundary tie break, repeated-match merge; use every positive interval if
  requested group count exceeds interval count. Warn once per adjusted variable
  and retain full mapping/fallback metadata. Use effective grid for both outputs.
- Compute local summaries once; retain and focal-population-average the same
  quantities. Quantiles mode optionally includes separately named mean rank,
  computed from the full distribution, not sparse reported quantiles.
- Keep grouped extraction under the existing `what = "exposure"` request for both
  modes, explicitly typed/labeled by summary/statistic. `measures = "exposure"`
  requests that grouped descriptive output. Shares retain classic P semantics;
  quantile/mean-rank profiles are not P matrices or segregation indices. Store
  them separately internally rather than mislabel rank values as shares. Preserve
  existing share columns where possible; document intentional rank output changes.
- Add binary pairwise segregation between named resolved global groups for continuous
  input, in either summary mode. Use full-variable neighborhoods, pair-conditioned
  local shares, pair focal mass/reference normalization, existing categorical
  formulas. Preserve full-resolution continuous scalar indices and normalizers.
- Handle empty/one-positive-interval cases without invented quantiles. Preserve
  eligibility, weighting, network direction, and search sharing.

## Ownership and implementation practice

Own package R/C++ source, native bindings, generated help, NEWS, and implementation
tests. Analyst owns README and methods/network vignette documentation after your
interface handoff; coordinate through this report. Do not edit those files or the
getting-started vignette. Regenerate roxygen/Rcpp outputs as needed and review diffs.

Provide a short compatibility path for explicit legacy continuous `localenv` and
`exposure` calls, with deprecation and conflict diagnostics, so the user's unfinished
vignette can remain unchanged. Preserve explicit old behavior where feasible; do
not accept conflicting new/legacy settings silently. Ordinal/categorical contracts
should remain unchanged. Document intentional pre-1.0 default changes.

Use a temporary source copy/private library for builds. Only you edit package source
initially; validation writes its reference checks in temporary files until handoff.
Review the known network-error abort and existing environment/toolchain evidence
before a full suite. If it blocks required checks, diagnose and fix its narrow cause
with validation; do not weaken tests or fold in unrelated optimization.

Acceptance: tiny independent examples agree; default/cuts/matching/typed output and
pairwise tests pass; continuous scalar/search invariance, equivalent counts/individuals,
ties, empty cases, multiple variables/bands, and relevant network regressions pass.
Run useful existing checks; report limits. Preserve user's vignette byte-for-byte
relative to edits you observe. Save incremental interface/candidate evidence here
so independent reviewers can start before your final summary.

## Results

### Current user-directed revision: interval apportionment

The user rescinded whole-interval matching/fallback after the original candidate
below. Requested probabilities are now retained. Each interval spans its original
eligible global rank mass; overlap with each requested group determines a fraction.
Exactly the same fractions apportion original focal group counts and local weighted
neighbor shares. This assumes a common within-category global-rank distribution
across locations; it does not assume or infer monetary tail bounds. Full-resolution
continuous scalar calculations remain unchanged. The previous matching/single-bin
empty-profile results below are historical and superseded by this revision.

The user then requested numerical diagnostics in neighborhood output instead of
in warning text. New retained `analysis$interpolation[[band]]` matrices and
`as.data.frame(..., what = "env")`/sf rows report `bottom_category_share`,
`top_category_share`, and `interpolated_share`. The last is local weighted mass in
input categories split across groups, not a standard error or a fraction of raw
geographic rows. Diagnostics are retained only with full/localenv output to avoid
an n-by-band allocation in indices-only mode; unrelated analyses receive NA columns.
Empty environments receive NA. Bottom/top refer to original first/last categories,
which can coincide for a one-category input.

`definition$summary$groups$allocation` records the interval-by-group fractions;
`source_shares`, original labels and zero-mass interval IDs provide the global
audit. `effective_probs` equals the requested grid; no cuts merge or move and no
insufficient-interval fallback applies. `summary$interpolation` records needed and
split_intervals; status is interpolated/ok. One brief warning per interpolated
variable explains the common-fraction assumption and points to the diagnostics.
Probabilities on actual interval boundaries require no apportionment warning.
Absolute cuts still require identified exact values; labels alone supply no income
bounds. Quantile values continue using the existing linear rank-CDF model.

Core candidate is executable at canonical R/C++ sources and builds in private
`/private/tmp/spseg-interpolation-TwBwiP/lib` (Apple Clang recipe from T012).
Focused implementation tests and the eight applicable regression files passed
(`focused.log`, `regression.log`, install log under that root). The prior independent
reference file at that run still described superseded matching. Core source is ready for refreshed
validation/science review. No commit, user-library install, getting-started edit,
or independent-reference edit was performed for this revision.

The top-category fixture has global category masses `.85/.15`, with local top
shares `.30/0` and requested deciles. Fractions allocate `1/3` and `2/3` of the
top category to the ninth/tenth groups; local top-decile shares are `.20/0`.
Focal-weighted exposure to the top decile is `7/85` for the first eight groups,
`12/85` for the ninth, and `.20` for the tenth. The implementation checks these
hand values, diagnostics by band, sf extraction, no-warning aligned boundaries,
scalar/reference invariance, and omitted diagnostics in indices-only output.
README/methods/help text now describe apportionment, not boundary matching;
roxygen regenerated `spseg_result.Rd` and `spseg_variables.Rd` successfully.

The refreshed independent reference file was run against the private installed
candidate on 2026-10-28. All numerical references passed; the full directory
reported four failures only because named row/column sums of the allocation
matrix were compared against unnamed vectors (reference lines 91, 93, 98, 139).
The allocation matrix deliberately retains original category and requested-bin
labels for the audit trail. Validation owns those assertions; no independent
test was edited. Final full-suite acceptance awaits their attribute-aware checks.
Also qualified new resolver `head`/`tail` calls with `utils::` to address the
previous package-check NOTE without changing semantics. This final qualification
does not require native recompilation; a refreshed private-source check remains
for Validation. The protected getting-started vignette was not touched.

Starting HEAD: `e2b8acac868cdca8d0119f79c5439557e6a41dff`, 2026-10-28.
Existing project/docs/vignette edits are preserved. Builds used private copies.

### Initial interface handoff (historical)

- `spseg_continuous(..., summary = list(type = "shares", probs = c(.25,.5,.75)))`
  is the default. `type = "quantiles"` uses the same grid for grouping and local
  evaluation; optional `mean_rank = TRUE` adds the distinct mean statistic.
- Endpoint 0/1 delimiters are accepted and stripped from evaluation. Exact shares
  accept `cuts` instead of `probs`; quantiles/cuts and interval/cuts error.
- Explicit old continuous `localenv`/`exposure` arguments select a deprecated
  legacy path; combining them with explicit `summary` errors.
- `definition$summary` will record requested/effective grid and group resolution;
  grouped ranks will be stored separately as `analysis$rank_summary`, not in P.
- `as.data.frame(..., what = "exposure")` will return `summary`/`statistic`,
  `from_group`, `to_group`, `prob`, `exposure`, `rank`, and `value`: share rows
  retain `exposure`; rank rows use `rank`, with no false target-group category.
- Continuous pairwise matrices will use resolved global group labels and the
  original full-variable neighborhoods. Scalar continuous indices stay unchanged.
- Source ownership: developer edits R/C++, generated help, NEWS, and implementation
  tests only; independent validation owns `test-continuous-summary-reference.R`.
  Analyst owns README/vignette updates. Candidate readiness will be recorded here.

### Executable candidate handoff

Core source candidate is now present in canonical `R/spseg_summary.R`,
`R/spseg_variables.R`, `R/spseg-internal.R`, `R/spseg_result.R`, and
`src/spseg_engine.cpp`. It builds and installs in the private library at
`/private/tmp/spseg-t011-hl4yQj/lib`. The first focused implementation test run
exercised defaults, rank averages, cuts, matching, pairwise, legacy and one-interval
cases; two failures were a comparison wrongly including P in rank mode (test
corrected) and named interval probabilities (source corrected). Rebuild from
canonical sources for independent checks; the first private snapshot is stale.
Core API is ready for independent validation; generated docs and regression-suite
adjustments are still in progress. No native binding signature changed.

`definition$summary` contains `type`, `requested_probs`, `effective_probs`, `cuts`,
`mean_rank`, `groups`, `matching`, `reasons`, `status`, and `statistics`.
`groups` includes labels/shares, requested/effective counts and boundaries; interval
groups additionally record original interval IDs/labels, zero-mass IDs, and cut-after
IDs. `matching` records requested, matched, cut_after, and tie columns.
`statistics` describes local/grouped columns with statistic/prob/label fields.
`status = "single_interval"` identifies a zero-column quantile profile; extraction
returns a typed empty table, or mean-rank rows if explicitly requested.

Validation's Clang-versus-GCC diagnosis is accepted as a toolchain limitation for
these checks: subsequent builds use its private Apple Clang Makevars, not a native
network workaround or a global compiler change. README/vignette ownership and the
unfinished getting-started file remain untouched by this assignment.

Roxygen regenerated successfully using `roxygen2::roxygenise(package.dir = ".",
load_code = roxygen2::load_source)`, changing `man/spseg.Rd`,
`man/spseg_result.Rd`, and `man/spseg_variables.Rd` only. The initial sandbox write
was blocked by the obsolete workspace root; the authorized retry succeeded.
No duplicate generation handoff is needed. New focused tests pass under Clang;
the wider suite is being updated where assertions assumed the superseded defaults.

Current core/test candidate is ready for an independent rebuild/rerun. Explicit
legacy quantile regression fixtures now call a test helper that asserts the
deprecation warning and returns the old specification; modern defaults and group
behavior have separate implementation/reference cases. No network assertion was
removed or weakened. Interval fixture tests unrelated to matching explicitly use
their supplied global boundaries, or assert the expected adjustment warning.
The Clang snapshot and installed private candidate are under
`/private/tmp/spseg-t011-hl4yQj`; `manifest.tsv` records actual copied files,
including untracked new R/test files, excluding compiled artifacts. Full-suite
results are being written to `tests-clang-3.log` (earlier runs are retained).

### Historical final outcome and verification (superseded)

The user subsequently replaced the interval matching contract with proportional
interpolation. This section records the earlier candidate only; see the revised
handoff below for the current implementation and acceptance limits.

Implementation is ready for management review/commit. Continuous defaults now use
quartile shares, explicit quantiles optionally include full-distribution mean rank,
and grouped extraction has the typed schema above. Whole-interval matching is
resolved once from original eligible totals and reused through grid/pycno
recompilation. Binary continuous pairwise results use the resolved global groups
and original variable neighborhoods, not pair-specific searches. Full-resolution
continuous scalar accumulation is unchanged. Explicit legacy arguments retain
their old independent behavior with deprecation diagnostics and new/old conflicts
error; ordinal/categorical calculations retain their previous contracts.

Final build/test commands (temporary source includes all current R/src/tests/data/
inst/man files, including new untracked implementation and independent tests):

```sh
R_MAKEVARS_USER=/private/tmp/spseg-T012-clang/Makevars \
  R CMD INSTALL --preclean \
  --library=/private/tmp/spseg-t011-hl4yQj/lib \
  /private/tmp/spseg-t011-hl4yQj/source
```

```r
.libPaths(c("/private/tmp/spseg-t011-hl4yQj/lib", .libPaths()))
library(spseg)
testthat::test_dir("/private/tmp/spseg-t011-hl4yQj/source/tests/testthat",
                   reporter = "summary", stop_on_failure = TRUE)
```

- Build/install passed: `install-clang-2.log` under the experiment root.
- All **nine** test files passed, including the 43-expectation independent reference
  file, focused implementation tests, legacy compatibility regressions, and the
  complete network file. No test warnings, failures, errors, or skips were reported
  in `tests-clang-3.log` (exit 0).
- `manifest.tsv` records the tested source files and MD5 hashes, excluding native
  build artifacts. A final comparison to canonical files found only
  `R/spseg_variables.R` (roxygen comments) and `man/spseg_variables.Rd` different.
  Parsed R expressions with `keep.source = FALSE` are identical. All runtime
  source and test behavior therefore still match the passing snapshot.
- Roxygen regenerated the three affected help files using the command recorded
  above; no Rcpp wrapper regeneration was needed because native signatures did not
  change. NEWS now consistently describes quartile shares; continuous decile/
  quantile defaults are confined to explicitly deprecated legacy documentation.
- `git diff --check -- R src man NEWS.md tests` passed. T013 found no material
  scientific defect; its NEWS/help clarity requests have been addressed. The
  efficiency review's nonblocking suggestions were left out of scope. Independent
  final review/acceptance remains with the assigned reviewers and management.

Developer-owned changes: the new `R/spseg_summary.R`, constructor/compile and
surface-resolution plumbing, local-summary reuse and pairwise/rank accumulation,
typed result extraction, generated help, NEWS, and implementation/legacy-regression
tests. Validation owns `test-continuous-summary-reference.R`. README and vignette
changes visible in the shared tree belong to their existing owners and were not
edited in T011. In particular, the user's unfinished getting-started vignette was
not edited, regenerated, formatted, restored, or rendered by this assignment.

Limits: macOS arm64/R 4.6.0 with private Apple Clang configuration was exercised.
The previously isolated GCC16 native exception abort is not repaired or claimed
fixed; no test was weakened to bypass it. No user/global compiler configuration,
active R session, or framework library was modified. No Windows/Linux check,
`R CMD check`, broad performance benchmark, or vignette render was run by this
assignment. No commit/push was made; management owns the authorized final commit.

### Final revised interpolation handoff (2026-10-28)

**Source is complete and frozen pending independent validation. No commit/push.**
This section supersedes the historical final outcome above. Requested interval
probabilities are retained; overlap in eligible global-rank space allocates each
original category fractionally to requested groups. The same fractions apply to
original focal counts, local target shares, and group normalizers. Allocation is
resolved before surface redistribution. No monetary tail bounds are inferred.
Full-resolution scalar indices remain separate and unchanged.

Current candidate is the canonical working tree at HEAD
`e2b8acac868cdca8d0119f79c5439557e6a41dff` plus uncommitted changes. Principal
candidate SHA-256 hashes (not hashes of the superseded matching snapshot):

```text
4eace09fb81616189c297e3545d000d1aa313d60e9f4bbe1f1621f6cbdc5a0a7  R/spseg_summary.R
ed1b8152edf669d2bc19436bc759addb8f0df0ac0f22e0b53f59f054541d61f5  R/spseg_variables.R
4e0d946e06395061c1c66322850519fc477ceafe0689163e72054a5e50e817c0  R/spseg_result.R
674d3f92c22bbacf9d49e495f687c1148e89f5b29ca9da1a182de3de1cf09fb2  R/spseg-internal.R
e141b760f95f4b05c3763a30693b482b9ee78b4a632ebcfcdc1af7ed32d35e62  R/spseg.R
7679323159c1f91109dee1ad3e67081f3ae9d8948d589b3a83bebd822b7845b0  src/spseg_engine.cpp
47a53041d52f7c71aba9bf2d8cf458698db7475bcc38fa86dae890a65ad45d80  tests/testthat/test-continuous-summary.R
451365a9732367b170f06edb809e8a92811d5e8a325fa3e08169c757e2e56a72  tests/testthat/test-multiple-variables.R
690d20c3c0245e32652043e9d10731b0b90c03e10539271d1e573a55063ccd83  man/spseg_variables.Rd
f650864e6f31481fb9e7c5ec5fbb6a3199ae4cb8752f2ff7293e7f3b5c1b2de8  man/spseg_result.Rd
7e005faa2a89969a508c3bfd7fe3e73c724d16ed4348b57f5cd390f9e9cc6189  NEWS.md
```

Per-neighborhood/band diagnostics are `bottom_category_share`,
`top_category_share`, and `interpolated_share`, available with retained
environments in data-frame/sf extraction and internal interpolation matrices.
They describe weighted local population shares, not estimation uncertainty.
Indices-only output avoids retaining these matrices. Allocation fractions and
original-category global shares remain auditable in summary metadata. The warning
is brief and contains no neighborhood-specific numerical report.

Completed developer checks: private Apple Clang build/install passed under
`/private/tmp/spseg-interpolation-TwBwiP`; focused interpolation tests and eight
regression files passed, including full network regressions. The refreshed
independent reference file was also executed: numerical assertions passed, but
four named-versus-unnamed allocation row/column sum assertions failed. Validation
owns these tests and must resolve the attribute expectations and rerun all nine
files. The final `utils::head`/`utils::tail` qualification was parsed successfully
but is newer than the installed snapshot; a fresh complete-source build is needed
for final acceptance. Scoped `git diff --check` passed. Roxygen help was regenerated;
native signatures did not change.

Outstanding limits: no revised-candidate R CMD check or vignette render by this
developer, no Windows/Linux validation, and the prior GCC16 exception-abort issue
is not claimed fixed. Independent/science/docs acceptance and management's central
design update remain their owners' responsibility. README and methods vignette
interpolation prose were updated during this user-directed revision (unlike the
earlier candidate). The protected getting-started vignette remains untouched.
