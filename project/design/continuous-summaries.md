# Continuous local environments, grouped summaries, and pairwise comparisons

Status: implemented and accepted by management on 2026-10-28 after revised
interval allocation tests and scientific/workflow/efficiency review (T011–T015).
The user revised interval handling to preserve requested quantiles through
interpolation and proportional allocation, with a warning. Earlier boundary
matching/fallback instructions in task briefs and reports are historical.
Final verification and its platform/vignette limits are recorded in T012/T014.
Inputs: [science](../tasks/T008-continuous-summary-science.md),
[analyst](../tasks/T009-continuous-summary-usability.md), and
[developer](../tasks/T010-continuous-summary-implementation.md) assessments.

## Common contract

Define global income groups once and choose a local summary. Compute each local
statistic once, retain it if requested, and average that same statistic within
each focal income group. Synchronization is the default. Both averaging levels
must retain their intended weights: neighbor/kernel/inclusion weights locally,
and original eligible focal-group population for aggregation. Do not pool raw
neighbor counts across focal people or weight geographic rows equally by accident.

| Mode | Local output | Focal-group average |
| --- | --- | --- |
| Shares | Fraction of the neighborhood in each global target group | Mean of those fractions among focal-group members; classic exposure matrix |
| Quantiles | Global rank at each requested within-neighborhood quantile | Mean of each local quantile rank among focal-group members |
| Mean rank, separately named | Mean global rank of neighbors | Mean of the local mean ranks among focal-group members |

Global grouping probabilities define who belongs to a group. Local quantile
probabilities define positions in each neighborhood distribution. These remain
distinct internal operations, but the user explicitly prefers one public `probs`
vector to govern both. Do not expose independent group/local probability controls.
Mean rank is a moment of the full local rank distribution, not an average of a
few reported quantiles.

Grouped quantile summaries are descriptive rank-context profiles, not a classic
share-valued exposure matrix or a new segregation index. Their rows need not sum
to one. A mean of local quantiles is generally different from a quantile of the
focal-weighted mixture of environments; pooling raw counts can change weights too.
Use explicit output type/statistic metadata. For implementation, retain
`measures = "exposure"` and `what = "exposure"` as the grouped-experience request
and extractor in both modes. Preserve classic P storage/semantics in shares mode;
quantiles mode has separate descriptive rank-profile storage and clearly typed
rank/probability/statistic output. Document this explicit mode-dependent contract;
never label rank profiles as classic share matrices or segregation indices.

## Approved implementation target

```r
spseg_continuous(
  values = "income",
  summary = list(type = "shares", probs = c(.25, .5, .75))
)

spseg_continuous(
  values = "income",
  summary = list(
    type = "quantiles", probs = c(.25, .5, .75), mean_rank = TRUE
  )
)

spseg_continuous(
  values = "income",
  summary = list(type = "shares", cuts = c(30000, 60000))
)
```

With `probs = c(.25, .5, .75)`, the global grouping uses `c(0, .25, .5, .75, 1)`.
Shares mode has four focal/target groups. Quantiles mode has four focal groups
and evaluates three local quantiles at `.25`, `.5`, and `.75`. Implicit endpoints
complete the groups without adding local minimum/maximum columns.

The user selects probability grouping as the default in both shares and quantiles
modes, with `probs = c(.25, .5, .75)` when omitted. Absolute `cuts` remain an
explicit alternative. This resolves the earlier candidate decile default.

Accept explicit 0/1 endpoints as equivalent delimiters for convenience/migration,
stripping them before local evaluation; interior-only input is the preferred form.
Validate numeric, finite, strictly increasing probabilities in [0,1]; do not sort
or deduplicate invalid input silently. Exact ties can change achieved group sizes
without changing the requested local evaluation probabilities.

`cuts` is an alternative to `probs`, with implicit `-Inf`/`Inf` tails and a stated
right-closed boundary convention. Reject simultaneous explicit `cuts` and `probs`;
apply defaults only after this choice is resolved. Original-value cuts require
identified original values/boundaries, not merely ordered interval labels.

For the first implementation, allow absolute `cuts` in shares
mode. Science also identifies a coherent quantiles-mode extension: translate cuts
once to achieved eligible-population cumulative probabilities F(c), then use those
for local quantiles while keeping absolute focal groups. Analyst/developer flag
that this implicit translation may surprise users and needs rules for duplicate,
empty, and endpoint probabilities. This scope remains a user decision; it must not
silently introduce an independent second probability grid. Quantiles/cuts therefore
receives an actionable error in this initial implementation.

Use shares as the default for introductory exposure workflows, with quantiles
explicitly selected. Document the intentional pre-1.0 default change and provide
a deprecated path preserving explicit legacy continuous localenv/exposure calls.
Default exact-input groups are user-specified global quartiles;
interval inputs preserve that request using the interpolation model below.

## Resolve probabilities for interval input

The user now requests interpolation and proportional allocation rather than
matching requested probabilities to observed interval boundaries. Preserve the
requested grid even when it implies more groups than populated input intervals.
Warn that within-interval positions and memberships are estimated, and disclose
the assumption in documentation and result metadata.

Resolve the allocation from the original eligible global interval population,
before surface redistribution or neighborhood search, and reuse it across bands.
Fractional target membership must conserve each interval's population; fractional
focal membership must also govern the averaging denominators and pairwise focal
weights. Do not replace these weights with whole-interval assignments or equal
weights per geographic row. Quantile summaries evaluate the requested probabilities.

The implemented model is uniform global rank within each supplied interval,
with the same conditional rank distribution at every location. For interval g,
let its global rank span be [q_prev, q_g] and its positive regional share be pi_g.
Its allocation to requested group a is the length of overlap with that group's
rank span, divided by pi_g. Apply this same allocation matrix to original focal
counts, local target shares, and regional group normalizers. Positive interval
allocation rows sum to one; zero-mass intervals contribute zero. Requested group
shares equal their requested probability widths, within numerical precision.

Local quantile values invert the full piecewise-linear rank CDF; mean rank remains
a full-distribution moment. No monetary tail bounds or imputed person identities
are needed. Finer groups do not add observed within-interval sorting information.
One populated interval yields local shares equal to requested widths, quantile
ranks equal to their probabilities, mean rank .5, and zero binary pairwise indices
in nonempty environments. These are modeled outcomes, not evidence of actual
absence of sorting. No eligible mass remains an actionable input error.

Warn once per variable when requested groups split supplied intervals; aligned
boundaries require no warning. Store requested/effective probabilities and the
allocation matrix. Retained neighborhoods include bottom_category_share,
top_category_share, and interpolated_share; the last is the population share in
split source intervals, not uncertainty or bias. Absolute-cuts scope is unchanged.
Independent validation must cover focal/target conservation, requested-grid
retention, grouped/pairwise weighting, and unchanged full continuous indices.

## Pairwise comparisons

The pairwise option compares two named global groups using binary
categorical segregation formulas. It uses their local shares, renormalized within
the pair, and original focal mass from that pair. It is independent of whether
the displayed local summary is shares or quantiles. It loses within-group income
detail, depends on the chosen cuts, and must not replace the full-resolution
continuous scalar indices.

Preserve neighborhoods selected from the full eligible variable
population, consistent with current categorical pairwise behavior. Selecting k
members of the pair instead requires a different neighbor-selection specification.
State handling of globally absent groups and locally empty pairs explicitly.

Group-to-group exposure shares are directional; binary pairwise segregation
indices are symmetric. Contrasts between group quantile profiles and dyadic
rank-pair decompositions are additional estimands, not automatic equivalents of
the binary indices. Do not claim quantile-based alternatives cannot exist; defer
new index definitions and associated properties until separately specified.

## Data rules and acceptance

- Preserve weighted global midranks and keep tied exact values together. Record
  requested versus achieved group shares, which need not be equal quantiles.
- For interval data, preserve the requested grid through documented interpolation
  and proportional focal/target allocation. Independently test conservation,
  uneven interval populations, cuts crossing one interval, and zero/single-positive
  interval cases, including warnings, metadata, and grouped/pairwise weights.
- Use variable-specific eligible population and current neighborhood semantics.
  Make absent groups and undefined environments explicit; do not silently drop
  focal residents from denominators. Directed networks remain outward summaries.
- Verify the tiny independent examples in T008/T010, individual/count equivalence,
  ties, empty cases, and typed extraction. Preserve existing continuous H/R/D/S,
  normalizers, and searches across summary changes; grouped/pairwise outputs may
  change with their defining cuts.
- T011 records the implementation decisions and ownership. Require independent
  validation, science review, and analyst documentation before management commits
  outstanding changes. The user-owned getting-started vignette stays untouched.
  Existing design arithmetic is not a passing test of an implemented feature.
