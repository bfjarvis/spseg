# T008 — Scientific contract for continuous summaries and pairwise measures

Status: in review
Owner: science officer. Reviewers: analyst, developer, and management.
Assessed revision: `e2b8acac868cdca8d0119f79c5439557e6a41dff` in the
canonical `spseg` checkout on 2026-10-28; the summary implementation files
inspected below match this revision. Concurrent work in `AGENTS.md`, project
records, `README.md`, `NEWS.md`, and the unfinished getting-started vignette
was preserved. This is a design proposal, not an accepted API.

## Brief

The user wants local environments and exposure summaries synchronized by default
through a common specification. Current continuous localenv defaults to within-
environment quantiles expressed as global ranks; exposure defaults to shares of
global income groups. Both are useful, but their independent syntax is confusing.

Assess the user's proposed two modes:
- Shares: define population groups (usually global quantile bins); local environments
  report group shares; exposure averages those same shares within focal groups.
- Quantiles: local environments report global ranks at chosen local quantiles;
  exposure averages corresponding local summaries within global focal income groups.
  The user also wants mean neighbor rank by focal income group considered.

Specify the mathematics, weighting, and interpretation. Distinguish global group
boundaries from local quantile evaluation probabilities, mean rank from a quantile,
and an average of local quantiles from a quantile of pooled neighborhood distributions.
Handle tied exact values and interval-valued input without silently claiming new
within-bin information; cover empty environments and asymmetric network weighting.

Assess pairwise segregation separately: binary comparisons of global groups,
pair-specific restriction/normalization and neighborhood populations, versus dyadic
rank-pair decompositions or quantile-distance contrasts. State which are established
measures, package extensions, or open research questions; do not assert there is
no quantile-based alternative simply because categorical formulas do not transfer.
Do not change continuous scalar indices merely by changing summary settings.

## Deliverable and boundaries

Write only this report. Use existing papers/current code and primary sources where
needed. Deliver concise definitions, recommended default semantics, a tiny numerical
example that distinguishes the key summaries, and unresolved scientific choices.
No proofs or literature review beyond what is necessary for the design, no package
changes, and no manuscript edits. Reuse T001 and prior research rather than restart.
Record source/revision provenance and mark in review when ready.

## Results

**Revised recommendation (supersedes the earlier independent-grid and
interval-default proposals).** Use one public vector of strictly increasing
interior probabilities `probs`, defaulting to quartiles
`c(.25, .5, .75)` for both exact and interval-valued input. Use the
resolved effective vector for both global focal-group boundaries and local
quantile evaluation. The implicit endpoints `0,1` complete the focal-group
partition; they do not add
local minimum/maximum columns. In `shares` mode those focal groups are also
the local target groups. In `quantiles` mode, report local global-rank
quantiles at the *same effective* `probs` and separately report local
**mean neighbor rank**. Average every retained local statistic within each
focal group.
This deliberate coupling makes the default compact and interpretable; do
not add a second public probability vector without an established use case.
Full-resolution continuous segregation indices remain unchanged. This is a
proposed package reporting convention, not a new Reardon index.

An optional vector of absolute income `cuts` can replace `probs`, not
accompany it. Keep the cuts as focal/target group boundaries and translate
them once to achieved global cumulative probabilities `u_j = F(c_j)` in
the eligible population; the outer income groups extend below the first
and above the last cut. In `quantiles` mode evaluate local rank quantiles
at those same `u_j`; in `shares` mode use the original absolute groups for
local targets. Thus `cuts` also produces one coupled grid. This is coherent
but may surprise users: local quantile columns are labeled by their derived
probabilities and source cuts, and fixed monetary cuts do not preserve the
rank-invariance of the group summaries when income distributions change.
The underlying continuous scalar indices remain rank-based. For
interval-valued input, the probability-grid resolution below applies to
`probs`; absolute `cuts` retain their separate observed-boundary rule.

### Definitions and interpretation

Fix one eligible population for the continuous variable, with mass `m_i` at
focal row `i`. All ranks and group boundaries refer to that population. Let
`B_a` be disjoint global income groups with focal mass `n_ia` at row `i` and
total `N_a = Σ_i n_ia`. These group boundaries are fixed across locations and
bands; their requested and achieved masses should be recorded. Let `L_i,b`
be the **normalized** local distribution of global ranks at focal row `i`
and band `b`, using the existing neighbor search, population weights, kernel,
and KNN boundary fractions. The focal unit remains in its own environment
under the current selection rule. Here `b` denotes a spatial band, not an
income group.

For shares, the local report is
`s_i,b,t = Pr_L_i,b(R in B_t)` and the group exposure is
`E_a,t,b = Σ_i n_ia s_i,b,t / N_a`. The diagonal `E_a,a,b` is an isolation
share and every valid row sums to one. These are the unnormalized exposure
semantics of [Reardon and O'Sullivan 2004](../references/Reardon%20and%20O%E2%80%99Sullivan%20-%202004%20-%20Measures%20of%20Spatial%20Segregation.pdf),
PDF pp. 17–18, with the package's stated local-environment conventions.
They are distinct from the normalized multigroup exposure *index* in
[Reardon and Firebaugh 2002](../references/Reardon_et_al_2002_Measures_of_Multigroup_Segregation.pdf).

For quantiles, use the shared **boundary/evaluation probabilities**
`u_1,...,u_h` in `(0,1)`: the supplied `probs` or probabilities derived
from `cuts`. Define
`Q_i,b(u) = inf{r: L_i,b([0,r]) >= u}` on the global-rank scale and
`M_i,b = ∫ r dL_i,b(r)`. The local report has `Q_i,b(u)` and `M_i,b` as
distinct statistics. The group summary is
`A_a,b(u) = Σ_i n_ia Q_i,b(u)/N_a` and
`A_a,b(mean) = Σ_i n_ia M_i,b/N_a`.
These are **focal-group averages of local summaries**: `A_a,b(.5)` is an
average of local median *ranks*, not a median of the group's pooled
environment. `A_a,b(mean)` is an average local mean rank, not a median.
Call this output a **rank-context summary** (or a clearly labeled
`mean_local_quantile`/`mean_local_mean_rank` table) rather than a
category-to-category exposure matrix. Each value lies on `[0,1]`; rows do
not sum to one. It describes the income-rank positions encountered by
members of focal group `a`. The averaging is linear in focal mass, not in
geographic support rows.

Reardon and O'Sullivan's exposure concept motivates averaging what each
group experiences; its published exposure output is a *share*, not these
quantile summaries. [Reardon 2011](../references/Reardon_2011_Measures_of_Income_Segregation.pdf)
establishes the use of income ranks for income segregation, not an index
based on averaging local quantiles. The quantile/mean summaries above are
descriptive package proposals. They should not inherit claims about the
bounds or exchange properties of established segregation indices.

### A four-person check

Take four equally weighted people with exact incomes 10, 20, 30, 40,
global weighted midranks `.125, .375, .625, .875`, and positions
`0, 10, 10.1, 0.1`, respectively. At an unweighted radius of `0.2`, each
person's environment is the nearby pair, including themselves. The two
focal people in the bottom-half group therefore have local rank
distributions `(.125,.875)` and `(.375,.625)`, with equal local weights.

| Bottom-half focal summary | Value | Meaning |
| --- | ---: | --- |
| Share of bottom-half neighbors | `.5` | Each local pair has one member below `.5`. |
| Average local median rank `A_bottom(.5)` | `(.125+.375)/2 = .25` | Average of the two local inverse-CDF medians. |
| Median rank after pooling their local environments | `.375` | Distinct pooled inverse-CDF operation, with masses `.25` at each of the four ranks. |
| Average local mean rank | `(.5+.5)/2 = .5` | A rank moment, distinct from either median. |

This verifies the meanings without relying on implementation output.
With unequal focal masses, replace the simple average by `n_ia/N_a`;
neighbor weights inside `L_i,b` remain a separate operation. In this
fixture the two local pairs are disjoint, so the pooled distribution
contains each of the four people once with equal weight. In general,
overlapping local environments can contribute the same neighbor more
than once to a pooled distribution; that reflects the experiences of
multiple focal residents.

### Interval-boundary checks

For four supplied intervals with eligible population shares
`(.10, .60, .20, .10)`, the actual cumulative boundaries are
`(.10, .70, .90)`. Requested quartiles `(.25, .50, .75)` independently
match `(.10, .70, .70)`; merging the last two gives effective boundaries
`(.10, .70)` and three groups with achieved shares `(.10, .60, .30)`.
Quantiles mode evaluates local ranks at `.10` and `.70`, and shares mode
uses those same three interval unions. This illustrates the resolution
lost by matching without inventing within-interval membership. Four
equal-share intervals match all requested quartiles exactly. With only
three populated intervals of shares `(.20, .50, .30)`, the fewer-interval
rule uses both actual boundaries `(.20, .70)`; neither mode claims four
groups. These checks specify the proposed contract, not tested output.

### Data and boundary rules

- **Exact values and ties.** The current code assigns one global weighted
  midrank to each exact tied value. An empirical local quantile selects the
  first distinct rank attaining `u`; it does not interpolate between exact
  values. For a `probs` boundary inside a tied atom, assign the whole atom
  to the lower group under a right-closed boundary convention. Retain the
  requested `u` for local quantile evaluation, and report achieved focal
  group masses: nominal equal-sized groups may be unequal or empty. For
  `cuts`, `F(c)` includes all values equal to `c`, so a cut never splits a
  tie. Reject `cuts` whose achieved `F(c)` values are not strictly
  increasing and interior; otherwise they would duplicate a local
  quantile column or create an empty endpoint group. Mean rank uses the
  same midranks. These are proposed boundary rules, not established code
  behavior.
- **Grouped continuous input.** The current linear CDF interpolation in
  *global percentile space* allows inverse local quantiles, local mean
  ranks, and shares of modeled rank intervals without imputing dollar
  incomes. For **focal and target groups**, first aggregate eligible mass
  over the supplied intervals across all locations. Drop zero-mass
  intervals from matching, retaining their identities in the report. If
  `m` intervals have positive mass, their `m-1` interior cumulative shares
  are the candidate boundaries. With fewer than `length(probs)+1`
  populated intervals, use every candidate boundary. Otherwise match
  each requested probability independently to its nearest candidate in
  absolute probability distance; break exact ties toward the lower
  candidate, then merge duplicate matches. This is a deterministic
  projection onto *observed* interval boundaries, so no focal membership
  is split or inferred within an interval. It can produce fewer effective
  groups even when `m` is large enough; that is preferable to forcing a
  requested cut onto a distant interval boundary. **Evaluate local
  quantiles at the effective boundaries**, not the original requests.
  Record requested/effective probabilities, matched interval identities,
  achieved group shares, and every fallback or merge; notify the user
  whenever the effective grid differs. Local quantile *values* can still
  use the existing within-interval rank-CDF interpolation, which models
  rank positions but does not assign focal people to artificial groups.
  With one populated interval, there are no interior boundaries or local
  quantile columns; retain the separately labeled mean-rank summary and
  flag the absence of resolved income groups. With no eligible mass,
  grouping is undefined. Absolute `cuts` should initially be accepted
  only at supplied interval boundaries; a cut inside an interval needs
  an income-within-bin model that rank interpolation does not identify.
  Do not infer monetary endpoints from labels, including a top-coded bin.
- **Missing and empty environments.** Only eligible mass enters the
  variable's focal group totals. An absent focal group has an undefined
  row. If any positive-mass focal row lacks a local environment, the
  existing exposure behavior marks the band undefined; mirror that
  policy for quantile summaries rather than silently averaging only
  the remaining focal residents. A later conditional-on-coverage
  option would need both its revised denominator and coverage reported.
- **Asymmetric access.** On directed networks `L_i,b` describes the
  outward reachable environment from `i`. Group averaging is still
  mathematically defined, but `E_a,t,b` need not equal a reversed
  exposure, and the quantile summary is likewise directional.
  Disconnected observations remain in the eligible regional ranking
  while local environments contain only reachable people. KNN
  infeasibility follows the current error policy. Published results
  derived for symmetric proximity should not be claimed automatically
  for this extension.

### Pairwise segregation and the summary mode

`Pairwise` has at least three distinct meanings. First, a binary
segregation index can compare two selected **global groups** `B_a,B_t`.
For each focal environment, restrict its group shares to those two and
renormalize by their combined local share; weight focal rows by their
population in either group and normalize against the combined regional
population. This is the established category-pair interpretation used
by the current *categorical* pairwise machinery. If the focal group
pair has no eligible local residents at some focal row, the convention
for that row must be explicit. The construction depends on **shares**
of the selected groups, whatever summary mode is used for presenting
local environments. It may be an optional discretized comparison for
continuous input; requesting it must not coarsen or replace the
continuous information/diversity/root/dissimilarity scalars.

Second, a `KNN-100` neighborhood selected from *all eligible income
records* is different from one selected from the union of two chosen
groups. The former is the current categorical pairwise convention;
the latter needs pair-specific neighbor selection to deserve the
`KNN-100` label for that pair. Both are coherent estimands if named
and documented. Third, a decomposition into narrow rank-slice pairs,
or a contrast between groups' averaged local quantiles, addresses
different scientific questions. A signed difference such as
`A_a,b(.5) − A_t,b(.5)` is a useful **contrast** but is not a
binary segregation index, is not normalized to `[0,1]`, and has no
automatic analogue of the categorical pairwise information or root
index. Thus a quantiles presentation mode can coexist with optional
share-based binary pairwise indices; it does not rule out future
quantile-based research measures.

**Scientific choices for management's proposal.** (1) Default both modes
to quartile probabilities; for interval input, project that grid onto
observed cumulative interval boundaries once for the eligible population
and use the effective grid everywhere. (2) Keep `probs` and absolute
`cuts` mutually exclusive, with implicit outer groups and mean neighbor
rank as a distinct quantiles-mode column. (3) Preserve ordinary scalar
indices and their reported normalizers exactly when summary settings
change. (4) If adding pairwise binary indices, choose and label whether
neighbor selection is for all eligible people or the selected group pair,
especially under KNN. (5) Decide whether an optional future `coverage`
policy is needed; the safe initial rule is undefined output for a band
with a missing positive-mass environment.

**Evidence and limits.** The design is based on the source definitions
above, the [T001 inventory](T001-method-provenance.md), and read-only
inspection of `R/spseg_variables.R`, `src/spseg_engine.cpp`,
`R/spseg_result.R`, and `vignettes/spseg-multiple-variables.Rmd` at the
assessed revision. The current code independently prepares `localenv`
and `exposure`, defaults exact exposure to weighted deciles, defaults
continuous localenv to quantiles, uses supplied intervals for grouped
exposure, and implements pairwise scalar indices only for categorical
variables. No code, test, or vignette was changed, and no package
execution or broad literature review was run. The analyst's and
developer's design reports were still pending at this assessment;
review them before accepting a unified public syntax.
