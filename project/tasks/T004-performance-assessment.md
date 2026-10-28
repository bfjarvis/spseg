# T004 — Performance evidence and benchmark priorities

Status: in review
Owner: optimization. Reviewers: developer, validation/testing, and management.
Starting revision / worktree: `de48d12a3e9d29ca09eb289873f1877ba5582117`,
canonical checkout, assessed 2026-10-26. Existing tracked edits: `.Rbuildignore`,
`.gitignore`, and `vignettes/spseg-getting-started.Rmd`; untracked `AGENTS.md` and
`project/`. No existing package-source edits were present. Only this report is
owned and edited by this assignment.

Final documentation check: `e2b8acac868cdca8d0119f79c5439557e6a41dff`,
2026-10-26. No changes to `src/`, `R/`, or `tools/` since the assessed revision.
Existing edits to management records and T002/T005 were preserved. The report
remains in review; this check adds no executable performance evidence.

## Brief

- Read `AGENTS.md`, `project/GOALS.md`, and `project/TEAM.md`; reconcile prior
  optimization recommendations with the current code and benchmark script.
- Preserve the user's prior emphasis on `output = "full"`, all requested outputs,
  and reuse of fixed geography across simulation steps with changing occupancy.
- Separate measured findings from hypotheses. Specify workloads, memory limits,
  baseline checks, and numerical-equivalence criteria for proposed optimizations.
- Include geographic and network cases; assess cache validity when eligible
  populations change. Do not assume a speedup or omit requested environments.
- Deliver ranked opportunities and a reproducible benchmark plan, reusing existing
  evidence. Defer fresh timed suites until coordinated with validation; concurrent
  test/build workloads would confound timings. Flag any missing baseline evidence.
- Edit only this report; no engine changes, caching API implementation, or new
  shared benchmark artifacts during this first assessment.

## Results

Assessment completed 2026-10-26; management acceptance pending. Only this report
was edited; no implementation or candidate revision exists. Recommend prioritizing
reusable geography and cross-band aggregation, supported by small exact kernel
improvements. All proposed gains require new full-output measurements. Keep the
measurement engine serial: independent cases will supply external parallelism.

### Provenance and current-code reconciliation

- Earlier executable assessment: `2334d39622176b58cb89273a587e0f9ef9d72718`,
  2026-10-22. Current assessment: `de48d12a3e9d29ca09eb289873f1877ba5582117`.
  `git diff 2334d39 de48d12 -- src tools` is empty. Current source inspection
  confirms the previously identified loops and allocations still exist.
- Public preparation/presentation changed: interval inputs now use
  `spseg_continuous(representation = "intervals")`; `spseg_grouped()` was removed.
  Ordinal environments now default to category shares; exact and interval
  continuous environments default to nine rank quantiles. Future benchmarks must
  use current constructors and explicitly record presentation settings.
- [Existing benchmark](../../tools/benchmark-multiple.R) remains a joint-versus-
  separate variable comparison: seed 881, default 5,000 uniform points, four
  variables, three radius bands, multigroup scope, three repetitions. It does not
  request full output, compare candidate revisions, measure memory, or exercise
  repeated geography. Its focal-search assertion is not numerical validation.
- [Decision record](../DECISIONS.md) contains no accepted scientific decision
  that resolves the scientific status of the estimators. This assessment targets
  equivalence to current behavior, not independent scientific validation.

Historical evidence was re-read, not rerun. Local temporary artifacts remain
available at assessment time; they are not durable project dependencies:

- [Benchmark driver](/private/tmp/spseg-review-benchmark.R), SHA-256
  `4529e98abb31664d2afb2884e9e25ec75aafa561bb46d6a94e2bff85053495ae`.
- [Timing log](/var/folders/t5/w6p00r9118l6rtbthpcvyf180000gp/T/spseg-review-qk1hcanq/benchmark.log),
  SHA-256 `5803300df5e8fc109055831d5477c2a70b9ea3c30c10f91f31d3e5e43b76a671`.
- [Build log](/var/folders/t5/w6p00r9118l6rtbthpcvyf180000gp/T/spseg-review-qk1hcanq/install.log)
  and [non-network test log](/var/folders/t5/w6p00r9118l6rtbthpcvyf180000gp/T/spseg-review-qk1hcanq/tests-nonnetwork.log).

### Historical measurements and limits

The driver used seed 20261022, 4,000 uniform unit-square points, prepared R inputs,
one native-engine warm-up, and three elapsed-time repetitions. It used
`keep_env = FALSE`, `keep_indices = TRUE`, Boost search, and normalized weights.
Environment: macOS arm64, R 4.6.0, BH 1.90.0.1, GCC 16.2.0, C++17;
the build log includes local `-O3 -march=native -fopenmp` flags. These are not
portable package defaults or evidence of engine threading. CPU model, memory
peaks, dependency versions beyond those recorded, and background load were not
captured systematically.

| Historical workload | Median elapsed seconds | Interpretation |
| --- | ---: | --- |
| Four-category multigroup, power-three biweight, radius 0.25 | 0.195 | Single-band reference only |
| Same, nine bands from 0.05 to 0.25 | 0.456 | Extra aggregation/kernel work despite one search per focal |
| Exact continuous, same one/nine bands | 0.971 / 3.573 | Prioritizes investigation of distribution processing; not a profile |
| 128 observed categories, radius 0.025, all measures, multigroup/both | 0.012 / 0.065 | Additional outputs differ; not an equivalent-output speedup |
| Same, exposure only, multigroup/both | 0.012 / 0.059 | Driver's `all.equal()` check passed; exposes unnecessary scalar-pair machinery |
| Four-category kNN, one target 400 / nine targets 50–400 | 0.204 / 0.475 | Repeated population selection and weighting remain |

These timings are neither full-output benchmarks nor measured improvements from
optimized code. The exposure-only comparison is secondary to the user's all-
results workload; its log says “identical,” but the actual check was `all.equal()`
with default tolerance, not bitwise identity. The other comparisons change work
or outputs and cannot establish achievable speedups. No native profile was
obtained: the earlier process sampler lacked permission. No cache or approximately
2x improvement has been measured.

Historical non-network tests passed after the temporary build included package
data. The full suite aborted on network unreachable-population handling, also
reproduced then against the installed package. No cause is assigned here;
[Validation T002](T002-validation-triage.md) owns current reproduction and triage.
Historical passes do not certify the current revision.

### Ranked opportunities for full results

Ranks express workload relevance, not measured speedups or implementation approval.
All retain every requested focal environment, scalar index, exposure matrix, and
pairwise matrix.

| Rank | Proposed change and current code | Benefit and qualification |
| --- | --- | --- |
| 1 | Reusable fixed geography; `src/spseg_engine.cpp:417–441`, `src/spseg_euclidean.h:25–73`, `src/spseg_graph.cpp:156–205` | Search objects are rebuilt per call and every focal launches discovery. Reuse index/topology and geometric neighborhoods across simulation steps; recompute population-dependent kNN selection. Highest relevance to repeated geography, but memory and invalidation govern feasibility. |
| 2 | Cross-band aggregation; `src/spseg_engine.cpp:443–466`, `src/spseg_neighbors.h:59–123` | Process sorted bands and reuse shell/prefix totals. For unweighted count input, aggregation can fall from O(K sum(m_b)) to O(K m_max + B K), excluding index/output work. Weighted normalized kernels cannot reuse plain sums. |
| 3 | Continuous ordering/integration reuse; `src/spseg_engine.cpp:167–181,252–267,291–301` | Sort candidate value codes once per focal and variable across bands; preserve code grouping while updating weights. Reuse known global-rank endpoint indices instead of searching for them again. Exact continuous full output still needs distributions for quantiles. |
| 4 | Relevant categorical pairs; `src/spseg_engine.cpp:146–148,337–345` | Individual rows contribute only to K−1 pairs containing their category, rather than all K(K−1)/2. Sparse count rows can visit pairs incident to populated focal categories without double-counting. Dense count rows have less benefit; complete output remains quadratic. |
| 5 | Specialized powers and reusable buffers; `src/spseg_neighbors.h:22–29,61–88`, `src/spseg_engine.cpp:249–275,321` | Specialize powers 1/2/3, retain generic powers, dispatch outside hot loops, and reuse capacity. Low structural risk, but floating-point differences need checks. Do not remove compensated population summation or boundary safeguards. |
| 6 | Count layout and repeated origins; `src/spseg_engine.cpp:273`, `src/spseg_graph.cpp:181–205` | Benchmark a row-packed copy versus existing column-major gathers. Reuse radius environments at shared origins, preserving each row's focal contribution and output. Copies/caches cost memory in every external worker. kNN focal-first selection prevents blindly sharing environments. |
| 7 | Pycnophylactic buffers/passes; `src/spseg_surface.cpp:29–65,127–137` | Precompute zone counts, reuse two buffers, and combine correction with convergence scans. Retain averaging stencil, conservation and stopping semantics. Relevant only to this surface workflow. |

For rank 2, the default integer-power kernel permits distance-moment formulas
such as `(1 - (d/b)^3)^3 = 1 - 3d^3/b^3 + 3d^6/b^6 - d^9/b^9`. Cancellation near
the boundary makes this a higher-risk exact-arithmetic proposal, not a ready
replacement. Keep direct evaluation as reference and potentially as fallback.
Restore original band order and duplicate-band results after any internal sorting.

Additional bounded candidates: precompute fixed interval-to-presentation-bin
overlaps (`store_environment`, lines 218–227); avoid a full distance sort in a
single-radius path when only filtering and a minimum positive distance are needed.
Neither warrants priority without profiles. Prior exposure-only shortcuts remain
valid special cases but do not characterize all-measure full output.

Do not prioritize internal parallelism, focal sampling, dropped vacant focal rows,
or replacing continuous distributions with interval approximations. They either
conflict with the requested workload or change outputs/estimands. Any approximation
would need a separate opt-in specification and scientific review.

### Fixed-geography cache contract to assess before implementation

Separate three layers rather than caching a prior final result:

1. **Geometry:** stable site IDs and coordinates, metric/boundary convention, or
   graph nodes/edges/costs/direction and observation-node mapping; reusable index,
   distances, and complete equal-distance shells. Invalidate or explicitly remap
   on changes. Do not retain unsafe views into mutable caller-owned R data.
2. **Neighborhood selection and kernels:** fixed radius support/weights can be
   reused for identical support, bands, power and normalization. For kNN, geometric
   ordering survives population changes, but shell masses, stopping shell,
   fractional inclusion, radius, zero-distance weighting and summaries must be
   recomputed for each distinct current eligible-population vector. Unchanged
   occupancy alone is insufficient if variable missingness or weights change.
3. **Measurement:** recompute compositions, global shares/ranks, normalizers,
   exposure groups, local summaries and indices from the current attributes.
   A geometry cache does not authorize reusing continuous value order or ranks
   across changing attributes; rank 3 above concerns reuse within one measurement.

Use fixed support with explicit zero population for vacant sites as the reference
simulation contract, retaining all focal environments. Row deletion is a different
input contract: radius inverse weighting currently uses the nearest positive
geometric distance, even if that site's population is zero. Deleting it can change
weights at coincident observations. Compare cached and fresh calls on the *same*
support. Explicit bands remain fixed; automatically chosen default bands can
change with population weights even when coordinates do not.

For regular Euclidean grids, ordered offset shells can replace per-origin lists;
clip offsets at edges and retain exact distance/tie semantics. Grid recognition
and translated-coordinate floating-point differences must not introduce or remove
ties. Current Euclidean input is not a toroidal metric: periodic wrapping would
be a separately specified extension, not a silent optimization.

For irregular geography, use bounded radius lists or geometry-only neighbor
prefixes, including currently vacant sites. If later eligible population requires
more neighbors, extend through a complete shell using exact search; never treat
cache exhaustion as population exhaustion. Preserve the true unreachable-target
error and its location/variable context. Do not cache a focal-sufficient singleton
as a complete geometric neighborhood for other observations at that origin.

A fixed graph may reuse topology even if observation-node mappings change, but
cached observation lists and selected neighborhoods then require remapping or
invalidation. Preserve directed reachability, transit nodes, zero-cost shells,
parallel edges and isolated nodes. No all-pairs matrix should be the default.

Measure cold preparation, first use, warm use, extensions and eviction separately.
For T measurements, compare `sum(fresh_call_times)` against
`prepare_time + sum(reused_call_times)`, including validation and key checks.
A 2x gain requires enough removable spatial work: if fraction f of baseline time
is removed, the ideal limit is `1 / (1 - f)` before cache overhead. f is unknown.

### Ranked reproducible benchmark plan (proposed; not executed)

All primary workloads use `output = "full"`, `measures = "all"`, categorical
`scope = "both"`, default power-three normalized biweight, and current local
presentation defaults. Record every output dimension. Test unweighted, inverse,
exponential, and `normalize = FALSE` as targeted variants, not a giant factorial
suite. Use explicit bands. Prespecify seeds below and save generator parameters;
generate fixtures before timing and use identical fixtures for both revisions.

| Order / seed | Workload | Variants and purpose |
| --- | --- | --- |
| W1 / 4001 | Fixed 50×50 and 100×100 unit grids; four observed categories; T=1,10,50 measurements | Radius bands 2,5,10 and population-kNN targets 10,50,100. Separate label-only swaps at fixed occupancy from 10%/30% moving vacancies; permute a fixed occupied label pool each step. Add a second variable with explicitly changing omissions. Compare fresh calls with proposed warm/cold/extension paths. |
| W2 / 4002 | Uniform and clustered unit-square points; n=2,500 and 10,000; four-category, five-level ordinal, and exact-continuous variables together | One band 0.1 versus nine from 0.02 to 0.1; kNN target 100 versus nine targets 20–100. Include tied continuous values and a near-unique distribution. Isolates band/distribution reuse with all results. |
| W3 / 4003 | n=2,500; K=4,32,128 observed categories, then dense/sparse count matrices | Sparse rows have at most two populated categories; include fractional counts and explicit common radius population. Radius 0.05, kNN target 50, and one/three bands. Isolates relevant pairs, layout, and output-memory costs. |
| W4 / 4004 | 25×25 and 50×50 lattice networks; hop and nonnegative cost versions; one/three observations per occupied node | Fixed topology, label swaps and moving vacancies. Explicit transit nodes; use hop bands 2,5 and kNN 10,50. Weighted edges use specified costs 0,1,2; separate directed, duplicate-edge and disconnected correctness fixtures. Run after Validation triage. |
| W5 / 4005 | 100×100 and 200×200 pycnophylactic grids from reproducible adjacent square polygons; four/eight count columns | Full calls plus isolated smoother timing; max_iter 200 and convergence 3/6. Include holes and zero-target zones; record actual iteration counts if later instrumentation exposes them. Separate surface construction from measurement. |

W1–W4 also need small independent reference cases before large timing. Large sizes
are conditional on memory preflight; do not silently downsize failing candidates.
For full local-presentation coverage, add explicit ordinal quantiles and continuous
shares/interval inputs using the current API. No real-world dataset is required.

Protocol for the later coordinated run:

1. Validation establishes a usable baseline and management schedules uncontended
   timing. Export pinned baseline/candidate sources into isolated directories and
   private R libraries, excluding stale objects. Record revision, dirty diff,
   source manifest, OS/CPU/RAM, R/dependency/compiler versions, flags and BLAS.
2. Preserve the serial engine; cap dependency threads at one *in benchmark worker
   subprocesses*. Warm once, then take at least five paired repetitions with
   reproducibly alternated baseline/candidate order. Batch very short calls so
   timer resolution does not dominate. Report raw observations, median, spread
   and paired ratios; increase repetitions only to resolve uncertainty.
3. Time public `spseg()` calls as the primary outcome, including preparation and
   full result construction. Use native-engine timing/profiles secondarily to
   attribute costs. Compare outputs outside timing; do not time fixture generation
   or result comparison. Retain and inspect a full result from both paths.
4. Measure peak process RSS with a platform-appropriate external monitor in
   separate runs and report its units/method; R object size misses native caches.
   Record output object bytes, cache bytes, temporary peaks and teardown behavior.
   Initially release previous full results between steps; separately assess a
   declared multi-step retention policy if the simulation actually retains them.
5. After serial measurements, run independent case throughput at 1,2,4 workers,
   capped by available cores/memory. Record cases/minute, total peak memory,
   per-worker preparation and warm-cache behavior. Avoid serializing live native
   pointers between workers; include per-worker cache replication/initialization.

### Numerical and memory acceptance criteria

These are proposed engineering gates for Developer/Validation review, not a
relaxation of existing tests or a claim of scientific validity:

- Compare all returned numerical analyses, environments, knn summaries and
  normalizers; require exact row/column/band order, names, dimensions, NA masks,
  statuses and input preservation. Permit only explicitly documented operational
  changes such as cache counters or a different recorded call, not blanket
  attribute removal. Preserve caller RNG state.
- Aim for bitwise equality when only reuse changes. For reordered arithmetic,
  report maximum absolute/relative errors, not merely `all.equal()`. Proposed
  ordinary-case bound: `abs(a-b) <= 1e-12 + 1e-12 * max(abs(a), abs(b))` for finite
  entries. Preserve tighter existing requirements and existing specially justified
  endpoint cases; any relaxation requires explicit review. Quantile crossings
  need independent threshold cases because tiny CDF errors can move discrete
  quantiles. Near-constant denominators also require separate stress cases.
- Retain compensated population stopping, focal-first inclusion, proportional
  complete boundary shells, the kNN bandwidth adjustment and zero-distance rules.
  Verify selected population against target with the current 64-epsilon-scaled
  check; do not replace it with the generic output tolerance. Preserve closed
  radius boundaries, including the square-root-rounding regression.
- Extend [neighbor tests](../../tests/testthat/test-boost-neighbors.R),
  [local-output tests](../../tests/testthat/test-local-environments.R), and
  [network references](../../tests/testthat/test-network.R) with cold/warm cache,
  forced extension/eviction, changing weights/omissions, empty origins, duplicate
  coordinates, row permutations, and changed topology/cost/bands/kernel inputs.
  Compare small Euclidean cases to brute search and small graphs to the existing
  independent all-pairs reference. Validate every step, not just the initial state.
- A minimum storage estimate for environments is `8*n*B*sum(env_columns)` bytes;
  full exposure and scalar pairwise matrices add quadratic category storage.
  Peak memory is higher because native accumulators, R inputs and output matrices
  coexist. Preflight against a proposed aggregate budget of 70% physical RAM;
  divide by worker count after accounting for coordinator memory. Report exact
  chosen byte limits before running, and record skipped cases with their reason.
- Proposed cache cap per worker: the smaller of 256 MiB and 10% of its available
  worker budget, adjustable only as a reported workload parameter. Test reaching
  the cap; exact fallback/eviction must remain correct and within total budget.
  Account for vectors, indices and overhead rather than just packed distances.
  Large memory increases must show a throughput benefit at intended concurrency;
  no unbounded all-pairs cache or hidden output omission qualifies as a gain.
- Report single-call and amortized speedups and break-even step count. Accept a
  performance claim only when repeat-to-repeat uncertainty supports it and all
  numerical gates pass. Preserve no-cache behavior for workloads where caching
  loses. Management should set deployment-specific runtime targets after baseline
  measurements; there is no evidence yet to promise 2x.

### Handoff and checks not run

Fresh work consisted of source/diff inspection and verification of the historical
driver/logs. No fresh build, package test, timed suite, memory profile, vignette
render or `R CMD check` was run. No code, cache API, shared benchmark artifact,
central record or colleague-owned document was edited.

Management: retain this as a deferred performance backlog behind the current
getting-started vignette, transparent income/network examples, and network error
repair. When optimization is scheduled, review ranking/budgets and coordinate
timing with Validation after T002. Developer should then scope geometry reuse and
exact population-dependent selection, including memory and invalidation costs.
Optimization's next executable deliverable would be a pinned full-output baseline
for W1/W2 and agreed network cases, followed by separately measured candidates.
Validation and Science retain independent correctness and scientific interpretation.
