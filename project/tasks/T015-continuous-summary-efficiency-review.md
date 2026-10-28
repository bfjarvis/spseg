# T015 — Bounded efficiency review of the summary candidate

Status: in review
Owner: optimization. Reviewer: management.

The user approved implementation, useful tests, aligned reviews, and a final commit.
Review the current T011 candidate and its approved design, especially
`R/spseg_summary.R`, preparation/recompile paths, and `src/spseg_engine.cpp`.

Confirm summary work is resolved/reused appropriately, no repeated neighborhood
searches or dense all-pairs income representations were introduced, grouped-only
outputs avoid per-location retention, and pairwise cost follows resolved group
count and requested scope. Look for material accidental time/memory regressions,
not optional future optimization. Separate complexity reasoning from measurements.

Write only this report with a concise acceptance recommendation and any concrete
blocking issue. Do not edit source/docs/tests, run timed benchmarks while other
agents build/test, restart historical optimization work, or change model settings.
Do not touch the user's getting-started vignette, commit, or push.

## Results

**Current recommendation: accept the revised interpolation/allocation path from
the bounded efficiency-review perspective, subject to the other reviews.** The
latest source inspection below supersedes the earlier matching-candidate approval.
No executable tests, timings or memory measurements were run.

### Revised interpolation candidate — 2026-10-28

The user replaced whole-interval boundary matching/fallback with interpolation
and proportional allocation preserving the requested quantiles. This instruction
supersedes the earlier contract; at inspection time the design document and T011
final handoff still described whole-interval matching. Their owners must align
those records before final acceptance. This review makes no independent scientific
choice between the contracts.

Inspected the revised uncommitted source at the same HEAD `e2b8aca`; only T015 was
edited. Snapshot SHA-256 hashes:

| File | SHA-256 |
| --- | --- |
| `R/spseg_summary.R` | `ddb3a34f0829040f5a0d8ca24adc5c0d3d03087303cc35279cc87737de0c4c17` |
| `R/spseg-internal.R` | `674d3f92c22bbacf9d49e495f687c1148e89f5b29ca9da1a182de3de1cf09fb2` |
| `R/spseg_variables.R` | `ed1b8152edf669d2bc19436bc759addb8f0df0ac0f22e0b53f59f054541d61f5` |
| `src/spseg_engine.cpp` | `7679323159c1f91109dee1ad3e67081f3ae9d8948d589b3a83bebd822b7845b0` |

Findings (n rows, K original intervals, G requested summary groups, S reported
statistics and B bands):

- **One reusable allocation, not expanded people or income pairs.**
  `R/spseg_summary.R:64–96` constructs a K-by-G fraction matrix using overlap of
  global-rank intervals. Resolution costs O(KG) time/storage, including temporary
  normalization/metadata arrays. The retained allocation documents the specified
  model; it is not an n-by-K observation expansion or K-by-K income-pair matrix.
  Exact observations continue to use integer maps and sparse local distributions.
- **Allocation is resolved once from original totals.**
  `R/spseg_variables.R:417–437` uses the supplied `resolution` when present.
  `R/spseg-internal.R:225–231` passes the original `summary_config` through surface
  recompilation. No matrix construction or probability resolution appears inside
  focal/band loops. C++ receives the fraction matrix as an Rcpp matrix member,
  with no explicit per-focal matrix copy (`src/spseg_engine.cpp:108,176–180`).
- **Allocation arithmetic has an explicit G factor.**
  `add_group_mass()` scans G fractions per original interval. Interval local and
  focal group calculations therefore cost O(KG) per row/band when required,
  in addition to existing neighborhood/count processing. Default G is four.
  Very large requested G can make the dense allocation and zero-fraction scans
  costly; this is bounded by requested grouping dimensions, not hidden expansion
  to all distinct income pairs. No new optimization requirement is proposed here.
- **Search and summary reuse remain intact.**
  The single `find()` call is still before the band/model updates
  (`src/spseg_engine.cpp:515–543`). `summarize_environment()` is called once per
  row/model/band and feeds retained local output and grouped accumulation. Shares
  mode now assigns `local_group = summary` (lines 379–380), removing the earlier
  duplicate grouping pass. Quantiles still need distinct grouped-share arithmetic
  for binary pairwise results. Both focal and target group computations reuse the
  same prepared fractions; no group-specific search was added.
- **Grouped-only mode still avoids per-location retention.**
  The interpolation list and n-by-3 diagnostic matrices are guarded by retention
  flags (`src/spseg_engine.cpp:114,170–172,320–326`). Full/localenv output adds
  three diagnostic doubles per row/band (about `24*n*B` payload bytes), alongside
  requested environments. Indices-only output retains neither those diagnostics
  nor n-by-S environments; its new persistent allocation metadata is O(KG), with
  grouped/pairwise accumulators O(BGS)/O(BG²). This is intentional model metadata.
- **Pairwise dimensions remain G, and scope still suppresses pairwise work.**
  Constructor/finish loops use resolved summary groups, not original intervals
  (`src/spseg_engine.cpp:164–168` and `pair_groups` in `finish()`). The earlier
  exposure-only/default-scope extra pair scan remains a non-blocking observation;
  no scalar-family gating fix was introduced in this revision.

No material accidental expansion, repeated resolution or search duplication was
found in the inspected revision. The proportional model necessarily adds group
allocation arithmetic; this statement does not establish its runtime or numerical
correctness. Developer/Validation should rebuild and verify the revised allocation,
surface reuse and full-versus-grouped-only behavior; earlier matching-candidate
test passes cannot validate the changed semantics. Management should reconcile
the design/handoff and recheck hashes if implementation changes before commit.
No benchmarks, tests, builds, source changes or vignette edits were made by T015.

### Earlier matching-candidate assessment (historical)

**Earlier recommendation: accept from the bounded efficiency-review perspective.** No
blocking search duplication, dense income expansion, or grouped-only retention
regression was found. Correctness/scientific acceptance remains with the other
reviewers. This is source/complexity inspection, not measured performance.

### Reviewed candidate and scope

Assessed 2026-10-28 in the canonical checkout at HEAD
`e2b8acac868cdca8d0119f79c5439557e6a41dff`, including the uncommitted T011 candidate.
Read current AGENTS/GOALS, [approved design](../design/continuous-summaries.md),
and [T011 handoff](T011-continuous-summary-implementation.md). Existing source,
test, documentation, project-record and user vignette edits were preserved.
Only this report was edited; no commit or push was made.

The principal files had these SHA-256 hashes at the final source inspection
(HEAD alone does not identify this uncommitted candidate). Concurrent R changes
added variable names to resolution warnings; those paths were re-inspected and
preserve the complexity conclusions. The engine hash remained unchanged.

| File | SHA-256 |
| --- | --- |
| `R/spseg_summary.R` | `9ac937eff30a53ad94e14e22e5c37322ed8d75eac80aba06fd9c5aff39071233` |
| `R/spseg-internal.R` | `2bd0c788ef18e537f1cb6d45f25770e1bc00bfb46b32cc2c217e3a2e345dc2a9` |
| `R/spseg_variables.R` | `4037de6917c5fdf0493d9c0d1c95261d9d3436b18dff6af129f4a7165b38fc0b` |
| `src/spseg_engine.cpp` | `f184bb0ab2c11717075e84680e6f2eaa2618d6c619cce986710a2a5260bcfd24` |

### Evidence

Let n be support rows, K original income values/intervals, G resolved global
groups, S local summary columns, B bands, and m a focal neighborhood's size.

- **Resolution is outside the focal/band loops.**
  `R/spseg_variables.R:409–430` resolves once per synchronized variable.
  `R/spseg-internal.R:217–223` passes `config$summary_config` into surface
  recompilation, where `resolution %||% ...` reuses it. Interval matching therefore
  does not repeat per surface cell, band or focal row. The resolver has O(KP)
  matching and O(KG) mapping/metadata passes for P requested probabilities;
  default P is three. No observation-by-observation distance structure is created.
- **Search sharing is preserved.**
  `R/spseg-internal.R:167–174` makes one engine call. The engine still creates
  one spatial backend, groups equal eligible-population vectors, and calls
  `find()` once per focal before the band/model updates
  (`src/spseg_engine.cpp:468–517`). Summary mode, grouped profiles and new pairwise
  outputs do not launch a second search or select pair-specific neighborhoods.
- **The same summary vector supplies local and grouped outputs.**
  `src/spseg_engine.cpp:302–308` calls `summarize_environment()` once per
  focal/model/band when retention or synchronized grouped output needs it.
  Lines 364–369 average that vector; they do not reconstruct quantiles from
  retained matrices or pool neighbors into a new distribution. Mean rank uses
  the current full local distribution in the summary function. Sparse observed
  income remains a neighborhood-local vector of code/weight pairs, O(m), rather
  than an n-by-K indicator matrix (`R/spseg_variables.R:338–339,366–398`;
  `src/spseg_engine.cpp:274–289`). Existing count input remains a count matrix.
- **Grouped-only output streams through bounded accumulators.**
  With `output = "indices"`, `keep_env` is false; the environment list has length
  zero and n-by-S matrices are not allocated (`src/spseg_engine.cpp:109,155–166`).
  The per-focal S-vector is temporary. Grouped storage is O(BGS) for rank profiles
  or O(BG²) for share exposure, plus optional O(BG²) pairwise storage. Input rows,
  coordinates and geometry are omitted from ordinary indices-only results;
  exact-value metadata is compacted in `R/spseg-internal.R:72–76`. Explicitly
  requested surface audit remains a separate intentional retention option.
- **Pairwise dimensions follow resolved groups and requested scope.**
  `src/spseg_engine.cpp:111,159–163,373–382,423–440` uses G, not K, for synchronized
  continuous pair lists, accumulators and returned matrices. Multigroup scope
  suppresses them; localenv-only output allocates no pair accumulator and returns
  before pair updates. With requested pairwise scalar outputs, worst-case update
  work is O(nBG²), and scalar flags gate their arithmetic/returned families.
  User-requested very many groups can still make G² large; this is explicit group
  output size, not an accidental all-distinct-income pair matrix. Full output
  additionally retains O(nBS) local values, as requested.

### Non-blocking observations and handoff

Two bounded extra passes remain; neither justifies delaying this candidate or
reopening the earlier optimization program:

1. `pairwise` is enabled by scope without checking whether any scalar family is
   selected. An exposure-only request with default `scope = "both"` still allocates
   the five-slot pair accumulator and scans pairs, although no scalar pairwise
   matrices are returned. This pre-existing categorical behavior now also applies
   to synchronized continuous groups. It costs O(BG²) storage and O(nBG²) pair
   visits; default G is at most four, so it is not a distinct-income explosion.
   A later narrow change could gate pair machinery on any selected scalar flag;
   verify unchanged grouped output with exposure-only requests in both summary
   modes and scopes before adopting it.
2. In synchronized shares mode, `summarize_environment()` aggregates target shares,
   then lines 349–357 aggregate `local_group` again for pairwise use. Even with
   multigroup-only exposure this latter pass runs. It is an additional O(m) sparse
   or O(K) count-distribution pass, not another search or retained dense object.
   Reusing the shares vector would be a small later cleanup, subject to equality
   checks; quantile profiles still require a separate grouped-share calculation.

No build, install, executable test, benchmark, profiler, memory measurement,
vignette render or package check was run by T015. The T011 private-build results
are the developer's evidence, not independently reproduced results here. Inspection
also found existing focused assertions for search count and summary/pairwise
equivalence; their current pass/fail status belongs to Validation. No speedup or
peak-memory claim is made. Management can proceed with aligned correctness,
science and documentation reviews; recheck these paths if the hashed source files
change before the final commit.

### Closing confirmation — revised allocation, 2026-10-28

**Accept from the bounded efficiency perspective; no material blocker found.**
Rechecked the current uncommitted candidate at HEAD `e2b8aca`. The revised
allocation review above is already saved; this confirmation covers the latest
`R/spseg_summary.R` SHA-256
`4eace09fb81616189c297e3545d000d1aa313d60e9f4bbe1f1621f6cbdc5a0a7`.
Its additional `utils::head()`/`utils::tail()` qualification does not change the
allocation structure. Internal/variable/engine hashes still match the revised
review table. The design now records interpolation as the current instruction.

Allocation remains one O(KG) interval-by-requested-group matrix reused through
surface recompilation. No additional searches, observation expansion or all-income
pair matrix appears. Indices-only mode still omits per-location environments and
interpolation diagnostics; requested pairwise work remains sized by G and scope.
The non-blocking exposure-only pair-scan observation remains. This is complexity
inspection only: no tests, builds, timing or memory measurements were run. Only
T015 was edited; management can close this review alongside the revised correctness
evidence, with no speedup claim from T015.
