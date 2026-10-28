# T013 — Review implemented continuous summary semantics

Status: in review
Owner: science officer. Reviewers: validation/testing and management.

Read the approved implementation brief T011 and prior T008. Review the candidate
as developer records its interface/code handoff. Confirm local/grouped consistency,
focal/neighbor weighting, whole-interval matching, effective-grid reuse, unchanged
continuous scalars, and share-based binary pairwise restriction/normalization.
Rank profiles must be labeled descriptive summaries, not classic share exposure
matrices or new segregation indices. T011's practical output choice is to keep
`measures = "exposure"` and `what = "exposure"` with explicitly typed mode outputs.
Shares are the default; quantiles are explicit; quartiles are default in both.

Check the independent tiny examples and scientific documentation claims. Identify
material scientific defects versus optional future research. Do not broaden scope,
edit code/vignettes, restart proof work, or change other reports. Own only this
report; document a bounded acceptance checklist now and final candidate review
when available. Preserve private chat identities. No commits/pushes.

## Results

### Scientific acceptance checklist (before candidate review)

1. **Common grid and eligibility.** Shares is the default, quantiles is
   explicit, and omitted `probs` means interior quartiles in either mode.
   Optional 0/1 delimiters only complete groups. Compute global ranks,
   groups, and interval candidate boundaries from this variable's original
   eligible population, weighted by its population counts, once across all
   locations and bands. Exact ties stay whole; show requested and achieved
   group mass. Exact cuts work in shares mode only; quantiles/cuts must give
   a useful error rather than silently changing the estimand.
2. **Intervals (revised user decision).** Preserve requested probabilities.
   Allocate each input interval proportionally across overlapping global
   rank groups, resolving the fractions once before spatial redistribution.
   Apply the identical allocation to focal group mass, local target shares,
   and regional group normalizers. Local quantiles use the requested grid.
   Warn once when bins are split and expose the allocation and its modeling
   assumption. A single positive interval can yield modeled quantile groups;
   their identical spatial distributions reflect the assumption. No eligible
   mass remains an input error. This replaces the earlier snapping/fallback
   acceptance criterion; see the revised contract below.
3. **Local and focal weights.** Each local summary uses the existing
   neighbor/kernel/population/KNN boundary weights and a local denominator.
   Each grouped summary averages the corresponding stored local value with
   original eligible focal-group mass, not equal geographic-row weights or
   pooled neighbor counts. Undefined positive-focal environments follow an
   explicit band policy. Directed networks retain outward interpretation.
4. **Output meaning.** Share rows are group-to-group P exposures and sum to
   one where defined. Quantiles are local inverse-CDF **global ranks** at
   the effective grid; their focal-group averages are means of local
   quantiles, not quantiles of a pooled environment. Optional mean neighbor
   rank is separately named and computed from the full local rank
   distribution. Quantile/mean-rank profiles are typed descriptive output
   under the existing exposure request/extractor, never labeled P matrices
   or segregation indices.
5. **Index separation.** Switching summary mode or grid leaves the
   full-resolution continuous H/R/D/S values and normalizers unchanged.
   Optional binary pairwise indices use the selected resolved group pair,
   pair-conditioned local shares and pair focal/reference mass, with the
   already selected **full-variable** neighborhoods. A KNN pairwise label
   must not claim pair-specific neighbor selection. Empty local pairs and
   globally absent groups need explicit outcomes; pairwise outputs do not
   replace full-resolution indices.

Review the candidate, independent tiny cases, and documentation against
these points after the developer's interface/code handoff. This checklist
records a proposed acceptance standard, not an implementation pass.

### Earlier candidate scientific review (interval findings superseded)

Reviewed the executable T011 candidate in the five named R/C++ source files
and its generated help, README, methods vignette, NEWS, and focused tests on
2026-10-28. The starting tracked revision is
`e2b8acac868cdca8d0119f79c5439557e6a41dff`; current source is still
uncommitted. The private Apple Clang install used for the probes is the T012
candidate snapshot (`4c12bffc1e46f1b485f65eb25e494d482335c22038c60993662470b6c8a3bb08`),
which predates at least one subsequent small source fix. Source conclusions
refer to the canonical files inspected, not a claim that the private install
is the final package build.

**No material scientific defect found in the reviewed paths.** R resolves
the interval grid from original eligible interval totals and passes that
resolution through surface recompilation. It retains whole-interval group
maps and uses the effective probabilities for local quantiles. The C++
engine normalizes local kernel/inclusion/population contributions, stores
the selected local summary, and averages that same value with focal group
mass. Mean rank uses the full local rank distribution. Exact tied values
use one weighted midrank. Pairwise H/R/D/S use pair-conditioned local
shares, pair focal mass and pair regional normalizers, within the selected
full-variable neighborhood. The scalar continuous accumulators and
normalizers operate on the original fine-grained local distributions and
do not read the summary grid. Rank profiles have separate internal
storage and typed extraction; they are not put in P.

Independent checks support those conclusions. The T008 four-person fixture
returned bottom-group average local median rank `.25` and mean neighbor
rank `.5`, distinct from the hand-derived pooled median `.375`. A separate
three-record weighted-values fixture (masses `2,1,1`) returned grouped
mean ranks `.375,.625`, grouped median ranks `.25,.5625`, and binary
dissimilarity `.5`; an equivalent exact count-table run agreed. A skewed
four-interval fixture with shares `(.10,.60,.20,.10)` warned once and
resolved quartiles to `(.10,.70)`, with achieved group shares
`(.10,.60,.30)` and local quantile columns at those effective positions.
T012's independently derived eight-person reference test passed 12
expectations against its isolated candidate snapshot. These probes are
small reference checks, not a general property proof.

**Actionable documentation issue.** The current `NEWS.md` has a new bullet
correctly saying continuous summaries default to quartile shares, but a
later bullet in the same 0.9.0.9000 section still says retained continuous
environments default to quantiles. Clarify that later statement as history
or update it to the new default. The `localenv` help text should also say
unambiguously that its quantile default applies only on the deprecated
explicit legacy path. The README and current methods vignette otherwise
distinguish share P exposures from descriptive rank profiles and state the
whole-interval and full-variable pairwise conventions accurately.

**Remaining acceptance gate.** T012 reports that the initial full candidate
test directory failed in `test-local-environments.R` because older tests
assumed the former default and absent continuous pairwise outputs. The
developer is updating that suite. Rebuild and rerun after its handoff;
recheck source changes against this review before final acceptance. No
source, generated documentation, or vignette was edited in this review.

### Revised interval contract: proportional allocation

The user's subsequent decision **supersedes interval snapping, merged
boundaries, and fewer-interval fallback**, including the earlier candidate
approval above. Requested quartile probabilities remain quartiles even
with fewer than four observed intervals. This amendment assesses the
revised mathematics and the current source plumbing; it does not extend
the earlier executable-candidate acceptance to the new implementation.

Let input interval `g` have regional share `pi_g` and global-rank endpoints
`q_(g-1), q_g`. For requested group `a` with endpoints `u_(a-1), u_a`, define
`alpha_ga = length([q_(g-1),q_g] intersect [u_(a-1),u_a]) / pi_g`
when `pi_g > 0`. Zero-mass intervals contribute nothing. For row counts
`c_ig` and normalized local original-interval shares `p_ig`, use
`n_ia = sum_g c_ig alpha_ga`, `s_ia = sum_g p_ig alpha_ga`, and
`P_a = sum_g pi_g alpha_ga`. The allocation rows sum to one, and regional
group shares equal the requested probability widths up to numerical error.
The group exposure is `sum_i n_ia s_ib / sum_i n_ia`; grouped quantiles
and means use those same fractional focal weights. Pairwise calculations
use `n_ia+n_ib`, local share `s_ia/(s_ia+s_ib)`, regional prevalence
`P_a/(P_a+P_b)`, and total pair mass. Neighborhood selection continues to
use the full variable population; allocating groups does not change kNN.

This is scientifically coherent under the existing linear local CDF model:
conditional on an input interval, the distribution of global ranks is
uniform over its regional rank interval and identical at every location.
The same `alpha` must therefore apply to all rows and to both focal and
target mass. Quantile *values* still invert the full piecewise-linear
local rank CDF; mean rank remains a moment of that full distribution.
No monetary values, finite tail endpoints, or identities of completed
individuals are inferred. An interval-membership record is treated as
fractional population mass under this model. Nonlinear pairwise indices
are plug-in indices of that modeled population, not expected indices over
randomly imputed individual assignments or estimates of their uncertainty.

For the previous skewed shares `(.10,.60,.20,.10)`, the four quartile
allocation rows are `(1,0,0,0)`, `(1/4,5/12,1/3,0)`,
`(0,0,1/4,3/4)`, and `(0,0,0,1)`. Regional group shares are now exactly
`(.25,.25,.25,.25)`. A fresh R process sourcing only the current resolver
confirmed these hand-derived fractions and emitted one interpolation
warning. In the current C++ source, `add_group_mass()` applies the same
matrix to original focal counts, normalized local shares, and regional
normalizers; local share summaries are reused for averaging and pairwise
calculations. Full-resolution scalar paths remain separate.

**Interpretation requiring explicit documentation.** Finer requested
groups increase output resolution under the model without adding observed
within-interval sorting information. Two groups wholly inside one input
interval have identical inferred spatial distributions. With one populated
interval, local shares equal requested widths, rank quantiles equal their
probabilities, mean rank is `.5`, and binary pairwise indices are zero in
nonempty environments. These outcomes express the imposed absence of
within-interval spatial sorting, not evidence that actual people are
unsegregated. The warning and allocation metadata must say so. The proposed
`interpolated_share` diagnostic is the local share in intervals split by
requested boundaries; it measures reliance on allocation, not uncertainty
or bias. The current roxygen source states this assumption and distinction.

No material mathematical ambiguity remains once that model is stated.
Management must update the design/briefs and final acceptance evidence so
that the superseded matching contract is not presented as current.
Independent validation must now check fractional focal and target weights,
count/membership equivalence, aligned boundaries, fewer/single intervals,
pairwise normalization, and scalar invariance on the rebuilt candidate.
No source or vignette edits were made for this amendment.

### Closing verdict on revised interpolation

**Scientific acceptance: no blocking defect found in the current source or
documentation.** The allocation matrix is applied consistently to focal
counts, local target shares, regional group mass, and pair-conditioned
index weights/normalizers. Group averages retain focal population weights;
quantiles and mean rank retain their distinct meanings. Continuous scalar
calculations and full-variable neighborhood selection remain independent
of reporting groups. Diagnostics use the normalized weighted local
distribution: first/last input-category shares and the share in categories
split by the requested grid. Generated help and the methods vignette
describe these as model-reliance diagnostics rather than error estimates,
and explicitly state the common within-category rank assumption.

Reviewed canonical runtime source identity: HEAD
`e2b8acac868cdca8d0119f79c5439557e6a41dff` plus uncommitted changes on
2026-10-28. The SHA-256 of the concatenated `path SHA256(file)` lines for
`R/spseg_summary.R`, `R/spseg_variables.R`, `R/spseg-internal.R`,
`R/spseg_result.R`, and `src/spseg_engine.cpp`, in that order with newline
terminators, is
`67bf009b1e3507c702e4293b83b6b09467c2b8fcef32de9fb5c060c3dfbab754`.

**Remaining verification gate:** T012 currently records revised interval
reference tests as prepared but not yet run against a rebuilt executable
candidate. Earlier matching/fallback test passes do not cover this revised
behavior. Final package acceptance therefore needs that revised build/test
evidence; no further scientific design decision is required by this review.
Status remains in review for management disposition. Only T013 was edited.
