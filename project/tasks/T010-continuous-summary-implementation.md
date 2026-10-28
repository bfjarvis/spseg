# T010 — Implementation map for continuous summary harmonization

Status: in review
Owner: developer. Reviewers: science officer, analyst, and management.

## Brief

Assess the user's proposed synchronized local-environment/exposure framework for
continuous variables. Shares mode should use shared global groups and report their
local shares and group-averaged exposures. Quantiles mode should retain local
quantiles in global-rank units and aggregate corresponding summaries by focal global
income groups, with mean neighbor rank considered explicitly.

Map current constructor/preparation, engine accumulation, retained output, metadata,
and extraction paths. Identify existing sharing/inheritance, exact versus interval
differences, ties, variable-specific eligible populations, and output compatibility.
Separate focal group definitions, target group definitions, and local quantile
probabilities; do not silently conflate them. Review T008/T009 when available.

For proposed pairwise measures, explain what categorical machinery could be reused,
and what requires a scientific contract (global-group binary comparison versus
dyadic rank-pair decomposition; pair-restricted versus original neighborhoods).
Do not change full-resolution continuous scalar indices due to summary settings.

## Deliverable and boundaries

Write only this report: concise current-behavior map, implementation options and
migration risks, memory/runtime implications, and a focused validation checklist.
No source/API changes, generated files, or test execution unless a small isolated
probe is needed to establish current behavior. No broad scan or benchmark suite.
Use the canonical spseg directory explicitly if the chat retains the old seg cwd.
Record inspected revision and mark in review when ready.

## Results

Assessed `e2b8acac868cdca8d0119f79c5439557e6a41dff`
(`codex/relational-neighborhoods`) on 2026-10-28. Existing workflow, vignette, and
project/report edits were preserved. Only this report changed; no code, API,
generated files, or tests changed or ran. The arithmetic below is hand-derived,
not package output. Proposed summaries are definitions for review, not newly
validated segregation measures.

Recommendation: compile one summary specification in R, compute its local
statistics once per focal unit/band, and reuse them for retention and
focal-population-weighted aggregation. Preserve the full-resolution continuous
scalar calculations, eligibility, and neighborhood searches. Keep rank summaries
separate from the classic share-valued exposure matrix and from pairwise indices.

### Current behavior and affected paths

| Path | Existing behavior and implication |
| --- | --- |
| `R/spseg_variables.R:112`, constructor | Separate `exposure` and `localenv`; continuous defaults to rank quantiles. Interval input rejects `exposure` overrides. |
| `.spseg_localenv():174` | Exact shares inherit exposure cuts only when local cuts are absent. Quantiles use independent local probabilities. Interval shares can overlap custom global-rank bins; default interval shares use supplied intervals. |
| `.spseg_exposure():262` and `.spseg_compile():299` | Exact exposure uses global weighted value cuts, ties intact, empty groups removed. Interval exposure uses one supplied interval per group. Compiles separate `env_map`/`env_breaks` and `exposure_map`; metadata likewise separate. |
| `src/spseg_engine.cpp:186`, `store_environment()` | Computes local shares or inverse-CDF rank quantiles only when environments are retained. No local mean rank. |
| `Model::update():240`, exposure block at 319 | Uses original focal group counts and local neighbor group shares. Focal and target use the same integer group map. It does not average retained local quantiles. |
| `Model::finish():348` | Divides exposure accumulators by global focal group totals; returns square matrices in `indices$multigroup$p`. Pairwise scalar accumulation is enabled only for categorical models (`type == 0`). |
| `R/spseg-internal.R:45`, `R/spseg_result.R:67` | Restores environment labels; extraction assumes every rank-summary column corresponds to one quantile probability. Exposure extraction assumes group-to-group matrices. Both need typed column/axis metadata for the proposed outputs. |

Exact input uses weighted global midranks, including exact-value count tables.
Interval input uses a linear CDF in global percentile space. Each variable has
its own eligible mass after omissions; kNN neighborhoods can therefore differ
between variables. Group definitions must use that variable's eligible population,
not an unrelated common population or a second neighborhood selection.

### A coherent proposed specification

Historical proposal, superseded by the user refinement below. Retained to explain
the initial assessment; independent public probability fields are no longer
recommended. Neither proposal is an implemented API:

```r
summary = list(
  type = "shares",
  groups = list(probs = c(0, .5, 1))
)

summary = list(
  type = "quantiles",
  groups = list(probs = c(0, .5, 1)),
  local_probs = c(.25, .5, .75),
  mean_rank = TRUE
)
```

`groups` defines global focal membership and, by default in shares mode, the
same target partition. `local_probs` evaluates each local CDF; it never defines
groups. Internally retain distinct focal and target maps even when their
definitions are shared. Do not add independent public overrides in the initial
implementation unless the user needs them. T009 uses `focal_groups` in quantile
mode; either spelling works, but a consistent `groups` field avoids a mode-specific
rename. Exact values can also support original-value boundaries for grouping.

The least ambiguous first delivery is synchronized shares with unchanged
share-exposure semantics; rank mode adds explicitly labeled grouped rank summaries.
Changing the current default from quantiles to shares is a user-visible behavior
change, not an implementation necessity. Management/user must choose the default.
Also decide calculation gating: `measures = "exposure"` currently means the
share-valued P matrix. Do not silently replace it with rank means. Prefer a
separately named grouped-rank request/output; its exact switch remains unresolved.

### User refinement: one public probability grid

Amended 2026-10-28 following the user's preference for constrained, synchronized
settings. This section supersedes the independent `groups$probs`/`local_probs`
recommendation above and its naming discussion with T009. No implementation or
test execution accompanies this amendment.

```r
summary = list(type = "quantiles", probs = c(.25, .5, .75))
summary = list(type = "shares", probs = c(.25, .5, .75))
summary = list(type = "shares", cuts = c(25000, 50000, 100000))
```

- Accept finite, strictly increasing, distinct **interior** `probs` in `(0,1)`.
  Reject explicit endpoints, missing/nonfinite entries, duplicates, unsorted or
  explicitly empty vectors with actionable messages; do not silently sort them.
  When neither field is supplied, the user-selected default for both modes is
  `c(.25, .5, .75)` (quartiles), superseding the earlier candidate deciles.
  Validate supplied fields before applying defaults, so an implicit probability
  default never conflicts with a user's explicit `cuts`.
- Compile two different internal tasks from that one vector: grouping boundaries
  are `c(0, probs, 1)`, whereas local evaluation positions are exactly `probs`.
  Thus three quartile probabilities define four global focal groups and three
  local rank summaries, without implicit minimum/maximum outputs. Shares mode
  uses the same global partition for targets. Do not add independent public
  focal/target/local probability controls.
- Keep exact-value ties intact. Requested global quartiles can produce fewer or
  unequal-sized realized groups; preserve mappings and achieved group masses in
  metadata. Continue evaluating every requested local probability, rather than
  modifying the probability grid to match the realized groups. Reuse the current
  deterministic quantile convention and empty-group policy unless separately
  revised. Coupling the requests does not equate the two operations.
- Make explicit `cuts` and `probs` mutually exclusive. For exact input in shares
  mode, interpret finite, strictly increasing cuts in original income units,
  complete with `-Inf`/`Inf`, and use right-closed bins for both focal and target
  membership. Never convert dollar thresholds directly into local probabilities.
  Recommend rejecting `cuts` in quantiles mode initially: cuts alone do not state
  local evaluation probabilities. Converting cuts to achieved global percentiles
  is a possible future contract, but introduces additional tie/empty-cut semantics
  and is unnecessary for this constrained design.
- Current interval input lacks numeric income endpoints; labels do not identify
  absolute thresholds. Reject absolute `cuts` for that representation. The user's
  subsequent interval-boundary matching decision below supersedes the earlier
  proposal to reject within-interval probability requests or fractionally allocate
  focal membership. No extra public probability knob is needed.

This coupling changes presentation/grouping only. Continuous H/R/D/S, reference
normalizers, ranks, eligibility, and spatial searches remain based on original
full-resolution inputs. Keep the earlier distinction between averaged local
quantiles, pooled quantiles, and mean rank. Add validation for endpoint completion,
three-probability/four-group dimensions, tied-group collapse, conflicting fields,
and mode/representation restrictions to the acceptance checks below.

### User refinement: deterministic interval matching

The user selects probability grouping, default quartiles, in both modes. For
interval input, resolve one effective grid per variable from its original eligible
population before any surface redistribution or neighborhood searches; carry it
through preparation rather than resolving per unit or band. This amendment
supersedes the earlier interval-overlap proposal for **group membership** below.

1. Keep supplied interval order, compute global eligible interval totals, and
   exclude zero-mass intervals from the boundary candidates and interval count
   `K`. Candidate interior boundaries are cumulative shares after each positive
   interval except the last. Retain original interval IDs and zero-mass IDs.
2. Let requested group count be `length(probs) + 1`. If `K` is smaller, use all
   actual interior boundaries: each populated supplied interval becomes a group.
   Otherwise match each requested cut to the nearest candidate in probability
   distance. Break midpoint ties toward the lower cumulative boundary; use a
   fixed small numerical tie tolerance (e.g. eight machine epsilons on this
   `[0,1]` scale), not a public policy option.
3. Merge repeated selected boundaries, preserving order. Nearest-cut matching is
   monotone; uniqueness is not guaranteed. Accept and disclose fewer effective
   groups rather than implement a global unique-assignment optimizer. Example:
   actual internal boundaries `(.10,.20,.90)` match quartiles to `(.20,.20,.90)`;
   the effective grid is `(.20,.90)`, giving three groups. This simple rule
   minimizes each cut's discrepancy and never splits a supplied interval.
4. Use `c(0, effective_probs, 1)` for global groups and exactly `effective_probs`
   for local rank-quantile evaluation. Shares use the same groups as targets.
   Do not retain unmatched requested local positions or add extrema implicitly.
   Focal memberships remain whole-interval sums. Local quantile **values** may
   still interpolate within intervals under the existing rank-CDF model; that
   does not split focal membership or identify unobserved dollar values.

Issue one warning per adjusted variable, independent of `verbose`, listing
requested-to-matched probabilities/interval IDs, any merged cuts, requested versus
effective group counts, and the insufficient-interval fallback when used. Do not
warn for an exact unchanged match. Store the complete mapping regardless of
warnings: requested/effective probabilities, complete effective endpoints,
cut-after original interval IDs, group-to-interval membership, achieved group
shares, removed zero-mass IDs, tie decisions, and collision/fallback reasons.
This allows audit without repeatedly warning per band or location.

For one positive interval, effective interior probabilities are empty: one share
group with valid exposure 1 where environments exist, and no local quantile
columns. Preserve an explicitly requested mean-rank statistic independently;
otherwise return a typed empty quantile profile/status, not an invented median.
All-zero eligible mass remains the existing actionable input error, before
normalization/matching. Zero-mass intervals introduce no repeated cut candidates
or division; their original positions remain auditable. None of this modifies
full-resolution scalar calculations, ranks, eligible populations, or searches.
No code or tests were run for this bounded refinement. Future acceptance cases
must cover exact matches, lower-boundary ties, collisions, insufficient intervals,
single/zero-positive cases, and equality of grouped/local effective grids.

### Statistical contract needed by the implementation

For one variable/band let `F_i` be the normalized local distribution after
neighbor population, kernel, and fractional kNN inclusion weights. Let `n_ia`
be the original eligible population in focal group `a` at unit `i`, and
`N_a = sum_i n_ia`. These are focal weights, not neighborhood population totals.
People in an aggregate unit share that unit's bespoke local environment.

- Shares: `s_ic = Pr_F_i(target group c)`; aggregate
  `P_ac = sum_i n_ia * s_ic / N_a`. Use the very same `s_ic` for retained local
  shares and exposure. Rows sum to one when target groups partition all eligible
  mass and environments are valid; the matrix need not be symmetric.
- Quantiles: `q_i(tau)` is the local inverse CDF of global ranks at `tau`.
  Aggregate `Q_a(tau) = sum_i n_ia * q_i(tau) / N_a`. This is an average of
  local quantiles, not a quantile of pooled neighbors and not a P matrix.
- Mean rank: `mu_i = E_F_i(global rank)`; aggregate
  `M_a = sum_i n_ia * mu_i / N_a`. A median or average of a few reported
  quantiles cannot substitute for this expectation. For exact values use
  `sum_g p_ig * z_g`; for the interval linear-CDF model use interval rank
  midpoints with local interval shares. Compute it from the full distribution.

Exact ties retain global weighted midranks and cannot be split to force equal
group sizes. Record requested cuts, realized value boundaries, achieved shares,
removed empty groups, and the existing inverse-empirical-CDF rule. Global cuts
stay fixed across bands. Directed network distances need no different averaging
formula and provide no general symmetry or pooled-neighbor interpretation.

Zero eligible focal mass contributes nothing, although its local map summary may
exist. Empty focal groups return NA/status. Recommend preserving the current
conservative failure semantics for undefined positive-mass local environments
rather than silently dropping them and changing denominators; per-group versus
whole-analysis invalidation should be an explicit decision. Do not renormalize
focal weights by selected kNN population or effective kernel mass.

### Tiny checkable example

Eight equally weighted people, four income values, two distinct units. Use raw
support and radius zero, so everyone experiences their unit's composition.

| Unit | Income 10 | Income 20 | Income 30 | Income 40 | Low focal count | High focal count |
| --- | --- | --- | --- | --- | --- | --- |
| U | 2 | 1 | 0 | 1 | 3 | 1 |
| V | 0 | 1 | 2 | 1 | 1 | 3 |

Each income has two people globally, giving ranks `(1/8, 3/8, 5/8, 7/8)`.
Global low/high cuts put 10/20 in low and 30/40 in high, four people per group.

| Local statistic | U | V |
| --- | --- | --- |
| Low/high shares | `(3/4, 1/4)` | `(1/4, 3/4)` |
| Mean neighbor rank | `3/8` | `5/8` |
| Local median rank, inverse empirical CDF | `1/8` | `5/8` |

The low group weights U/V by `(3/4, 1/4)`; high weights them by `(1/4, 3/4)`.
Consequently:

| Focal group | Average low/high shares | Average mean rank | Average local median rank |
| --- | --- | --- | --- |
| Low | `(5/8, 3/8)` | `7/16 = .4375` | `1/4 = .25` |
| High | `(3/8, 5/8)` | `9/16 = .5625` | `1/2 = .5` |

An unweighted mean over units would incorrectly give every group `(.5,.5)`
share exposure. Even pooling normalized environments with the correct low-group
focal weights gives income shares `(3/8, 1/4, 1/8, 1/4)` and median rank `3/8`,
not the average local median `1/4`. This difference does not require unequal
neighborhood sizes. Group means commute with this particular normalized mixture;
quantiles do not. Pooling raw neighbor counts would define yet another weighting
when neighborhood masses differ.

### Interval-input boundary

Historical alternative below: superseded for group membership by the user's
whole-interval matching contract above; retained as the original assessment.

Supplied intervals can immediately define shared focal/target groups, and adjacent
intervals can be merged without inventing membership. Custom rank cuts inside an
interval require fractional membership, not the current single-integer exposure
map. Under the existing linear-in-global-rank model, a source interval's contribution
to a rank group is its count times overlap width divided by its global rank width.
The same mapping must be applied to original focal counts and local target shares.

This is implementable with a sparse interval-to-group overlap map, but amounts to
model-based within-interval allocation, not observed income-group membership.
Science/user must approve it before enabling those cuts for exposure. Do not
interpret interval endpoints as income dollar values or recover dollar quantiles
without additional information. Zero-global-mass intervals require no division.
Phase-one conservative behavior can restrict focal grouping to supplied/merged
intervals while retaining explicitly modeled local rank summaries.

### Pairwise is a separate decision

Binary comparisons of named global groups can reuse the categorical accumulator
formulas after mapping distributions to shares: focal pair weight
`n_ia + n_ib`, local composition `s_ia / (s_ia + s_ib)`, and global pair reference
`N_a / (N_a + N_b)`. Skip zero focal-pair mass; define handling of empty local
pairs explicitly; globally absent groups give undefined comparisons. This yields
binned categorical indices, not the original full-resolution continuous indices.
In the example the low/high binary dissimilarity is `.5`.

Recommend retaining neighborhoods selected from the variable's full eligible
population, then conditioning composition on the pair, consistent with current
categorical semantics. A pair-only population-kNN search is a different analysis
and would need an explicit request and separate eligibility/stopping rules.

Exposure entries are directional group-to-group shares, whereas pairwise H/D/R/S
are symmetric two-group segregation indices. A dyadic rank-pair decomposition or
quantile-profile contrast would require its own scientific definition and
normalizer. Sparse quantiles do not reconstruct group shares or those indices;
this does not establish that no quantile-based segregation measures exist.
Use shares for the proposed categorical pairwise extension regardless of local
presentation mode, only if the user approves that particular meaning of pairwise.

### Implementation sequence, cost, and migration

1. Normalize the single public `probs` grid once in R (or shares-only `cuts`),
   deriving separate internal grouping boundaries and local evaluation positions,
   plus the mean flag and column descriptors. Record their
   relationship and realized partitions in one summary definition; retain the
   continuous scalar definition separately.
2. Refactor local summarization out of the storage-only routine into a small
   per-unit/band calculation reused by storage and grouped accumulation. It must
   run for requested grouped summaries even in indices-only mode, but need not
   allocate per-unit result matrices. Leave scalar accumulators unchanged.
3. Share mode needs `B x A x C` grouped accumulators; rank mode needs
   `B x A x (Q + mean)` rather than a square exposure cube. Compute local means
   in the existing distribution scan and quantiles from the same sorted CDF.
   Preserve sparse exact-observation processing; do not create an `n x K` matrix
   of every distinct income or rerun spatial searches for presentation changes.
4. Use distinct result fields/extractors for classic share exposure and grouped
   rank statistics. Add `statistic`/probability descriptors for both local and
   grouped rank columns; mean rows have no probability. The current extractor's
   all-columns-are-quantiles assumption must change. Retaining environments costs
   `O(n B C)` or `O(n B (Q + 1))`; grouped-only output need not pay that cost.
5. If approved, extract the categorical pair accumulator for mapped groups rather
   than duplicate a complete model or build one binary dataset per pair. Pairwise
   storage/work scales quadratically with group count and must remain requestable.

The mathematical contract, output gating, default mode, and naming require user
approval before coding. Both old arguments supplied with conflicting settings
must never be silently collapsed: either migrate with explicit conflict errors
and a documented transition, or announce a pre-1.0 breaking change. Retained-output
defaults and `indices$multigroup$p` meaning are particularly sensitive. Regenerate
roxygen/native wrappers only when the approved implementation requires it.

### Focused acceptance checks and peer review

- Use the eight-person table above to check local/grouped shares, mean ranks,
  quantiles, unequal focal weights, and distinction from mixture quantiles.
- Require identical scalar H/R/D/S and normalizers before/after changes to summary
  mode, cuts, local probabilities, and retention. Share exposure should change
  when its target groups change; that is not a scalar-invariance failure.
- Verify weighted individuals equal equivalent exact-value counts, ties remain
  intact, omitted values affect only their variable, and changing presentation
  does not change radius/kNN eligibility, boundary inclusion, or search results.
- Check supplied/merged intervals; if approved, test within-interval overlap maps
  for both focal and target masses against explicit hand arithmetic and totals.
- Cover empty groups/environments, directed networks, multiple bands/variables,
  grouped-only versus full retention, and long/sf extraction with typed columns.
  Pairwise checks must state the selected population and condition on named pairs.

Read T009 in review after forming the independent recommendation. Agree on one
summary framework, focal-population weighting, distinct rank outputs, and scalar
invariance. Differences for management to reconcile: consistent `groups` spelling,
exact scope of interval-group interpolation, and pairwise selected-global-group
comparisons versus cumulative-threshold comparisons. T008 still had pending
results at the final read; scientific approval remains outstanding. No literature
claims or new method validation are inferred from this implementation assessment.
