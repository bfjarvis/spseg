# T009 — User interface for synchronized continuous summaries

Status: in review
Owner: analyst. Reviewers: science officer, developer, and management.

## Brief

The user finds separate continuous `localenv` and `exposure` specifications confusing.
They want a common framework with synchronization by default:
- Shares mode: both summarize shares of the same global income groups, locally
  and averaged within focal income groups.
- Quantiles mode: both summarize global income ranks at requested within-neighborhood
  quantiles, locally and averaged within focal income groups; consider mean rank too.
- Pairwise segregation is desirable but its meaning must be specified separately.

Read current AGENTS/GOALS and relevant public documentation. Review T008 as it
becomes available; do not block your initial user-facing assessment on it.
Propose one clear compact syntax plus at most one alternative. Explain the
distinct roles of global grouping cuts and probabilities evaluating local quantiles,
while honoring the requested single-vector constraint.
Prefer one inherited/default summary definition with explicit overrides only where
useful. Show typical tidyverse/native-pipe workflows and clear extraction labels.
Separate mean rank, average local quantile, and pooled quantile descriptions.

## Deliverable and boundaries

Write only this report: a minimal before/after usage example, recommended defaults,
output tables/labels, migration concerns, and questions requiring scientific judgment.
Treat syntax as a proposal, not an accepted API. Do not edit vignettes or code;
preserve the user's getting-started edits. Keep this a bounded design assessment,
not a broad rewrite. Record inspected revision and mark in review when ready.

## Results

Assessed revision: `e2b8acac868cdca8d0119f79c5439557e6a41dff` (2026-10-28).
The current public interface and interpretation were checked in
`R/spseg_variables.R`, `man/spseg_variables.Rd`, `man/spseg_result.Rd`, and
`vignettes/spseg-multiple-variables.Rmd`. T008 and T010 were also read in their
current draft state; both still say results are pending, so this report records
an independent usability proposal and leaves method approval to science and
management. No code, vignette, or other report was changed, and no calculations
were run for this design-only assignment.

**Design refinement (2026-10-28):** the preferred proposal below supersedes the
earlier suggestion to expose separate `focal_groups` and `local_probs` fields.
Use one interior `probs` vector in both modes, with implicit 0/1 endpoints for
global grouping; quantile mode also evaluates local ranks at those same
interior probabilities. Do not automatically add local endpoint summaries.
The further user refinement sets probability grouping as the default in both
modes, with default interior probabilities `c(0.25, 0.5, 0.75)`. For interval
input, report the requested-to-effective boundary mapping; exact matching,
ties, and collisions remain for science/development review.

### User-facing assessment

The current interface makes the user's confusion understandable. For exact
continuous input, `exposure = list(probs = ...)` controls target income groups
used by the exposure matrix, while `localenv = list(type = "quantiles",
probs = ...)` controls local rank quantiles. For interval input, exposure uses
the supplied intervals and rejects a separate exposure cut specification.
Local summary choices do not change exposure or the scalar indices. Both
probability arguments are documented, but their different jobs are easy to miss
when they are separate constructor arguments.

The current split can look like this:

```r
income <- spseg_continuous(
  values = "income",
  exposure = list(probs = seq(0, 1, 0.1)),
  localenv = list(type = "quantiles", probs = c(0.25, 0.5, 0.75))
)
```

The concepts to expose in one framework are still distinct:

- **Focal groups** determine which focal people are averaged together.
- **Target groups** partition neighbors into categories whose shares are
  reported. In shares mode, using the same global cuts for focal and target
  groups is the default. The preferred constrained interface uses the same
  `probs` cutpoints in quantile mode to define focal groups and local quantile
  positions; these roles remain conceptually distinct even though users cannot
  tune them independently.
- **Local quantile probabilities** select points in each focal person's
  weighted local distribution of global income ranks. In quantile mode they
  reuse the interior global `probs` values; implicit 0/1 group endpoints do
  not imply local minimum/maximum rows.
- **Mean neighbor rank** is a weighted arithmetic mean. It complements local
  quantiles, but is not the median or another quantile.

### Proposed syntax and defaults

Use one constructor field, provisionally `summary`, to select a synchronized
local and grouped-summary mode. Within it, accept either one interior
`probs` vector or, for absolute income thresholds, one `cuts` vector; require
exactly one. Both must be nonempty, numeric, finite, and strictly increasing.
Additionally, `probs` values must lie inside `(0, 1)`. Users provide only
interior probabilities; implicit endpoints 0 and 1 complete the global
grouping partition. In shares mode, the resulting groups define target shares
and default focal groups. In quantiles mode, the same probabilities define
global focal-income groups and local rank-quantile evaluation positions. They
must not silently create local minimum/maximum rows; request endpoints only if the
API later chooses to support that separately. This is a proposal for review,
not an accepted API.

Default `probs` to `c(0.25, 0.5, 0.75)` for both summary types. With exact
input, use weighted global income quantiles as the boundaries. With interval
input, use the supplied interval boundaries as the available grid. When the
requested number of groups exceeds the number of populated intervals, fall
back to those supplied intervals rather than implying finer information.
Otherwise, map each requested global boundary to a valid cumulative interval
boundary and use the effective probabilities consistently for focal groups,
target shares, and (in quantile mode) local quantile positions. The exact
matching/tie/collision rule is unresolved and must not be guessed in the API.
Report both requested and effective probabilities so users can see the
coarsening. The requested group count is `length(probs) + 1`.

```r
income_shares <- spseg_continuous(
  values = "income",
  summary = list(
    type = "shares",
    probs = c(0.25, 0.5, 0.75)
  )
)

shares_result <- spseg(
  people,
  variables = list(income = income_shares),
  coords = c("x", "y"),
  bands = c(2, 5),
  measures = c("information", "exposure"),
  output = "full"
)

as.data.frame(shares_result, what = "env", variables = "income")
as.data.frame(shares_result, what = "exposure", variables = "income")
```

For exact income, `probs` defines global weighted quantile-bin boundaries;
the same realized groups are used as focal groups by default. Equal income
values stay together, so achieved group sizes can differ from target
proportions. For absolute thresholds, offer the mutually exclusive `cuts`
alternative in shares mode:

```r
income_cut_shares <- spseg_continuous(
  values = "income",
  summary = list(type = "shares", cuts = c(30000, 60000, 100000))
)
```

Treat `cuts` as increasing interior value thresholds, with implicit open-ended
tails, and require a nonempty numeric vector of finite, nonmissing values that
are strictly increasing.
This keeps the cut interface parallel to interior `probs` and avoids making
users repeat `-Inf` and `Inf`. It is meaningful for exact observed income.
For interval input, keep supplied intervals as the target and focal groups;
reject `cuts` unless a future contract explicitly supports cuts aligned to
interval boundaries. Splitting an interval would require an assumption about
within-interval ranks/counts.

For example, if three supplied intervals have population shares `.2`, `.5`,
and `.3`, their cumulative boundaries are `.2` and `.7`. A requested
`probs = c(.25, .75)` could map to those effective boundaries under the
matching rule approved by science/development. Suggested metadata:

```r
result$analyses$income$definition$summary$requested_probs
# .25 .75
result$analyses$income$definition$summary$effective_probs
# .20 .70
result$analyses$income$definition$summary$group_source
# "interval_boundaries"
```

This is an illustration of the notice and provenance, not a decision that
nearest-boundary matching is the correct rule. The reported mapping should
also identify fallback to supplied intervals when requested group count is
greater than the populated interval count. Whether empty intervals count as
available groups, and how equal-distance/tied boundaries and collisions are
resolved, remain explicit science/development decisions. A user-facing message
could say: “Requested probabilities 0.25 and 0.75 were represented by
interval boundaries 0.20 and 0.70; local quantiles and group summaries use the
effective probabilities.”

Quantile mode uses the same `probs` vector for both global focal grouping and
local rank-quantile positions. The simple contract should reject `cuts` in this
mode: absolute income thresholds do not directly specify local quantile
probabilities. A possible but more surprising alternative is to map each cut
once to its achieved global cumulative probability `F(c)` in the eligible,
weighted exact-income population, then use those same values for local
quantile positions while retaining the original cuts for focal groups. That
would keep a single shared grid, but the two uses still have different
interpretations (a global rank threshold versus a within-neighborhood
cumulative probability). It would require explicit rules for ties (right-closed
cut groups and the cumulative mass through the tied value), and actionable
errors when `F(c)` is 0/1 or when distinct cuts yield duplicate `F(c)` values,
since these create empty global groups or invalid/redundant local quantile
positions. It is not identifiable from the current interval input: interval
counts do not locate an absolute cut inside a bin. Recommend shares-only
`cuts`, unless the user affirmatively prefers this coupling and science agrees
to those rules. The grouped result is an average of focal people’s local rank
summaries, not the classic exposure matrix `P*`:

```r
income_ranks <- spseg_continuous(
  values = "income",
  summary = list(
    type = "quantiles",
    probs = c(0.25, 0.5, 0.75),
    mean_rank = TRUE
  )
)

rank_result <- spseg(
  people,
  variables = list(income = income_ranks),
  coords = c("x", "y"),
  bands = c(2, 5),
  output = "full"
)

as.data.frame(rank_result, what = "env", variables = "income")
as.data.frame(rank_result, what = "rank_summary", variables = "income")
```

In this mode, each focal person's local quantiles are calculated using the
same weighted empirical-CDF convention as the existing local environment
summary, and each returned quantile is expressed in global rank units. The
`rank_summary` rows average corresponding per-focal quantiles within the
specified global focal-income group. They do **not** calculate a quantile
after pooling all neighbors in that group. If `mean_rank = TRUE`, report the
per-focal weighted local mean rank and its focal-group average as a separate
statistic.

Recommended default for discussion: make the synchronized **shares** mode the
default, using the existing exposure decile cuts for exact values and the
provided intervals for grouped input. This preserves the familiar meaning of
the classic exposure output and removes the split configuration as the default.
It changes the current default retained local environment from rank quantiles
to shares. Users who need the existing local rank quantiles would request
`type = "quantiles"` explicitly. An alternative is to retain rank quantiles as
the default summary mode and make their grouped rank summaries the default
partner; however, that mode must keep its rank summaries distinct in name and
meaning from Reardon's group-to-group exposure matrix. Science and management
should choose the default after reviewing T008.

### Hand-checkable distinction among rank summaries

Take two equally weighted focal people in the same focal income group. Their
weighted local neighbor ranks are A = `(.375, .875)` and B = `(.125, .125,
.125, .875)`. Use the existing inverse empirical-CDF rule for local quantiles.

- A's local median is `.375`; B's is `.125`. The average of the two local
  medians is `.25`.
- Pooling all six neighbors first gives a median of `.125` (the first rank
  reaching cumulative weight one half). It is not `.25`.
- A's mean local rank is `.625`; B's is `.3125`. Their equally weighted
  focal-person average is `.46875`, distinct from the averaged medians.
- Pooling all six neighbors gives mean rank `2.5 / 6`, about `.4167`, because
  B contributes four neighbors and A contributes two. This differs from
  giving each focal person equal weight before averaging.

This example is intentionally about estimands and weighting, not a segregation
index. It makes clear that the report must label whether focal persons or
pooled neighbor mass gets equal weight. Recommended grouped summaries use
focal population weights within each focal group, applied to each focal
person's local quantile or local mean. Within a focal person's local
distribution, retain the package's neighbor population weights, spatial
kernel weights, and KNN boundary fractions. Empty local environments and
empty focal groups should have explicit status and missing summaries, rather
than silently changing denominators.

In shares mode, the analogous simple table is local target shares A =
`(.75, .25)` and B = `(.25, .75)` for low/high target groups. With equal focal
weights their low-focal exposure row is `(.5, .5)`, the average of A and B's
local shares. Unequal neighborhood sizes do not alter this focal-weighted
average; pooling neighbor counts first would be a different calculation.

### Output and migration

Keep `as.data.frame(result, what = "indices")` unchanged: summary settings
must not alter full-resolution scalar segregation indices, their normalizers,
or population eligibility. In shares mode, changing global cuts necessarily
changes the categories and entries of the classic exposure/isolation matrix.
In quantiles mode, grouped rank summaries are a separate descriptive output;
they do not replace the classic exposure matrix or scalar indices. The API
must decide whether `measures = "exposure"` remains available there through
an explicit target-group definition, or is restricted to shares mode with an
actionable error. For retained local results, continue the existing long
`what = "env"` table. Share rows use `summary`,
`category`, and `share`; rank rows use `summary`, `prob`, and `rank`. Add an
explicit grouped extractor such as the proposed `what = "rank_summary"`, with
columns for `variable`, `focal_group`, `band`, `statistic`, `prob`, and `value`;
`prob` is missing for the mean. Keep classic `what = "exposure"` as the
unnormalized matrix with explicit focal and target labels. Do not label an
average quantile as classic exposure or use the same output field for a
category share and a rank.

The `summary` default change would affect retained `full`/`localenv` output
shape, so announce it and document migration: current
`localenv = "quantiles"` probabilities move to `summary$probs` in quantile
mode; current exposure quantile cuts move to `summary$probs` in shares mode.
When both legacy fields use different probability vectors, there is no exact
one-to-one migration: users must choose one vector or retain separate legacy
behavior during a transition. Current value-based exposure breaks map to
`summary$cuts` in shares mode. For a transition period, an explicit legacy
`localenv` plus `exposure` pair could remain accepted with a
deprecation notice, but conflicting old and new specifications must error with
clear guidance. Alternatively, avoid deprecation and keep the arguments as
aliases if the maintainer prefers a longer compatibility period. These are
release-policy choices, not scientific requirements.

### Pairwise question for science review

Do not make “pairwise” automatically mean quantile-distance summaries. Current
`scope = "pairwise"` is a categorical scalar-index option. For continuous data,
one defensible extension would define binary comparisons at named global
cutpoints, then use the existing binary pairwise calculation; this is a set of
threshold-specific comparisons, not a generic pairwise rank measure. Another
possible proposal is a dyadic rank contrast, but it needs its own distance,
weighting, normalization, and population interpretation and cannot be inferred
from the current categorical API. Science should decide whether continuous
pairwise means (a) binary low-versus-high contrasts at selected thresholds,
(b) pair-specific exposure using shared neighborhoods, or (c) a newly defined
rank-pair statistic. Until that decision and independent reference cases exist,
retain existing categorical pairwise behavior and avoid presenting quantile
profiles as a segregation index.

Scientific choices still requiring review are the default summary mode;
whether quantile mode includes local/grouped mean rank by default; the default
focal-group cut scheme for rank summaries; focal-person versus pooled-neighbor
weighting (recommended above as focal-person weighting); treatment and labels
for tied/empty quantile groups; and whether any interval-data focal grouping
beyond supplied intervals is allowed. No choice here is implemented or
scientifically validated.
