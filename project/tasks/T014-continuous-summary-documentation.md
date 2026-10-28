# T014 — Document and exercise synchronized summary workflows

Status: in review
Owner: analyst. Reviewers: science officer, validation/testing, and management.

Read T011 and the design. After developer provides the concrete interface in T011,
update README and relevant multiple-variable/network/surface vignette examples
and explanations for the shared summary, quartile defaults, interval allocation,
typed grouped outputs, mean rank, and share-based pairwise comparisons. Include
compact independently checkable examples and migration from legacy arguments.
Developer owns roxygen/man/NEWS/source/tests; do not edit those files.

The user is actively editing `vignettes/spseg-getting-started.Rmd`. Do not touch,
render in place, autoformat, regenerate, or restore that file. It may be committed
as the user's unfinished work. Assess compatibility from a temporary copy only;
report issues without changing it. All vignette examples should prefer tidyverse
and native pipes. Preserve existing source edits.

Begin with an outline in this report while interface work proceeds. Once candidate
is available, execute relevant new workflows/render changed vignettes in temporary
copies with private libraries and capture results. Distinguish an independent
reference calculation from a package test. Keep examples small and avoid a broad
new vignette program. Provide a final user-facing consistency verdict and limits.
Only edit your assigned documentation files and this report. No commits/pushes.

## Results

Assessed base: `e2b8acac868cdca8d0119f79c5439557e6a41dff`; 2026-10-28.

### Documentation outline

1. Update README's quick interface summary with the shares default, explicit
   quantile mode, shared interior quartile grid, interval allocation assumption,
   and legacy migration note.
2. In the multiple-variable vignette, use one compact hand-checkable income
   example to explain local shares, focal-group averaging, local quantiles, mean
   rank, interval allocation, and pairwise group comparisons. Preserve the
   full-resolution scalar-index distinction.
3. Add concise references from the network and surface guides to the shared
   summary/exposure semantics where those input forms change interpretation;
   avoid duplicating the methods chapter.
4. Render temporary copies with a private library after T011 provides a concrete
   result/extraction schema. Review exact labels and pairwise behavior against the
   implementation, and report the untouched getting-started compatibility risk.

T011's executable handoff confirms `summary = list(type, probs, cuts,
mean_rank)`, shares as the default, exact cuts only in shares mode,
`what = "exposure"` for both summary types, and binary pairwise group indices.
The later approved interval rule is proportional allocation. Grouped rows carry
`summary`, `statistic`, `from_group`, `to_group`, `prob`, `exposure`, `rank`,
and `value` alongside base columns. Summary metadata is at
`result$analyses[[name]]$definition$summary`. README and the three owned
vignettes have been updated to this interface.

### Implementation and evidence

- README documents the default quartile share grid, explicit quantile mode,
  mean rank, exact cuts, proportional interval allocation and diagnostics, typed
  `what = "exposure"` extraction, pairwise scope, scalar-index invariance, and
  migration from legacy `localenv`/`exposure` settings.
- The multiple-variable vignette includes an eight-person, two-unit example
  with hand-computable local shares and weighted midranks, focal averages,
  quantiles, mean ranks, and binary pairwise dissimilarity. It distinguishes
  those summaries from unchanged full-resolution scalar indices.
- The network vignette demonstrates interval quantiles on network-derived
  neighborhoods and inspects effective probabilities, interval allocation,
  and interpolation diagnostics.
- The surface vignette demonstrates proportional allocation from source-zone
  intervals before surface construction, including typed extraction and metadata.
- `vignettes/spseg-getting-started.Rmd` was not edited, rendered, or restored.

Initial candidate and checks (2026-10-28):

- Installed the candidate into a private library from a complete-source copy
  using T012's private Apple Clang makevars at
  `/private/tmp/spseg-T012-clang/Makevars`. No global package installation was
  performed. The copied untracked `R/spseg_summary.R` was included.
- Snapshot manifest: `/var/folders/t5/w6p00r9118l6rtbthpcvyf180000gp/T/spseg-t014-9z7uwxn3/manifest.json`;
  72 files, combined SHA-256
  `c9cfa6b1e95c76831c0afde241963ce6ad832b850be018bc6a32c74cf228f8e7`.
  Current canonical SHA-256 values for `R/spseg.R`, `R/spseg_result.R`,
  `R/spseg_summary.R`, `src/spseg_engine.cpp`, `DESCRIPTION`, and `NAMESPACE`
  match this manifest, confirming the tested candidate source for the public
  interface and engine.
- The initial renders covered an earlier interval-boundary matching candidate.
  Those results are historical and do not validate the subsequently approved
  interval interpolation behavior.
- The first isolated render attempts exposed missing explicit `dplyr` attachment
  in the methods/network examples, which was added. The methods render then
  required the bibliography in its temporary resource path; it was copied
  there, and the final render succeeded. These were render-environment/example
  issues, not package runtime failures.
- The numerical prose in the compact example was cross-checked against the
  independent small-case calculations recorded in T010/T012. Rendering
  executes the examples but is not a substitute for scientific validation;
  the independent calculation and T012 validation remain the evidence for the
  stated values.

### Revised interval semantics and recheck

The user superseded nearest-boundary matching and insufficient-interval
fallback with proportional allocation in global-percentile space. For interval
inputs, requested probabilities remain the effective grid. A supplied category
crossed by a requested boundary contributes proportionally to adjacent groups;
the same allocation is applied to focal population and local target shares.
This assumes a common within-category rank distribution across locations. The
package warns when a category is split, and records global allocation fractions
and local `bottom_category_share`, `top_category_share`, and
`interpolated_share` diagnostics. Probabilities at supplied cumulative
boundaries need no split.

README and the multiple-variable vignette describe the revised model. The
network vignette now presents its two interval categories as preserving the
requested median through fractional allocation. The surface vignette describes
splitting source-zone counts before grid redistribution and updates the example:
with shares `.435` and `.565`, the first interval fills the first quartile and
contributes `.185` to the second; the requested `.25`, `.5`, `.75` grid remains
in place. Removed the earlier claims about nearest matching, changed effective
boundaries, and insufficient-interval fallback.

Rebuilt and installed the current complete source snapshot to a fresh private
library at `/private/tmp/spseg-t014-refresh.5wh8ei/library`, using T012's
private Apple Clang makevars. Its manifest is
`/private/tmp/spseg-t014-refresh.5wh8ei/manifest.json` (77 files; combined
SHA-256 `262e79490d008b02cd206538647b17ee2f15e3c4ada84ab68dba8b5bafdb6ff4`).
Rendered the three owned vignettes from copied
sources into that temporary root. All three rendered successfully; only
Pandoc's `--mathjax` deprecation warning appeared. Their chunks executed against
the refreshed installed package. The network road-preparation chunk remained
conditional and did not run because `sfnetworks` is not installed. The
getting-started vignette was neither modified nor rendered. Render logs and HTML
outputs are retained under `/private/tmp/spseg-t014-refresh.5wh8ei/`.

The numeric network and surface interval allocation examples were exercised by
their vignette renders. The small exact-income hand calculations remain
independently cross-checked in T010/T012. This is workflow/render evidence, not
a claim that the interpolation model is empirically identified: within-category
rank allocation is an explicit assumption. The earlier snapshot manifest and
render evidence above are retained only as historical evidence for the prior
candidate. T011/T012 own source-level implementation and broader validation
evidence.

Final runtime/schema comparison for closeout found one post-render change in
`R/spseg_summary.R`: calls to `head()` and `tail()` were qualified as
`utils::head()` and `utils::tail()`. Because these calls execute in the interval
allocation path, the private package was rebuilt and all three owned vignettes
were rerendered against that current code. They passed with only the existing
Pandoc MathJax deprecation warning. The refreshed manifest above includes this
change, and hashes for the summary implementation, orchestration, extraction
schema, engine, namespace, and description all match the canonical files. No
further relevant runtime or schema drift was found.

This assignment changes only README, the multiple-variable/network/surface
vignettes, and this report. Developer-owned roxygen/man/NEWS/source/tests were
not edited. No commits or pushes were made.
