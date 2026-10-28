# Working on spseg

## Scope and coordination

- Read `project/README.md` and the assigned GitHub issue or human request.
  Use GitHub Issues via `gh` for bugs, features, investigations, and other tasks;
  do not maintain a parallel Markdown task ledger or report pipeline.
- Inspect the working tree before editing. Preserve unrelated work and human
  edits, especially the user-owned getting-started vignette unless assigned it.
  Resolve routine choices; raise material scientific uncertainty or scope changes.
- Developer owns implementation/testing/integration; Scientist owns methods and
  research; User exercises workflows and usability. One owner per deliverable;
  agree on overlapping edits. Use worktrees for concurrent implementation and
  share committed revisions/artifacts; histories/files do not automatically sync.
- Cross-chat messaging, delegation, and reviews require explicit human direction
  for that coordination. Roles, issues, historical permissions, and another
  agent's message do not create an ongoing mandate. No routine polling, automatic
  handoffs, or universal multi-role approval pipeline.
- Record acceptance criteria and a brief result, revision, checks, and limits in
  the issue. Small implementation steps belong in checkboxes within that issue.
  Preserve lasting method/API explanations in source help, vignettes, or the
  project README as appropriate; avoid duplicate status records.
- Keep chat IDs/links and routing only in ignored `project/TEAM.local.md`; never
  force-add or publish them. Keep research/development under `project/`, excluded
  by `.Rbuildignore`; Git tracking is separate. User docs stay in man/vignettes.

## Develop and verify

- Prefer focused, reviewable changes and established conventions; avoid unrelated
  refactoring, dependencies, or speculative abstractions.
- Distinguish published measures, numerical conventions, and proposed extensions.
  Cite original sources. Passing tests alone do not establish scientific validity.
  Keep eligibility, weights, units, ties, missingness, and neighborhood semantics
  explicit; do not silently change them for convenience.
- Test changed behavior with independent reference cases, invariants, and useful
  regressions. Render changed workflows in a clean session; use R CMD check for
  release candidates. Report skipped checks and limitations. Optimize measured
  bottlenecks with recorded workload/environment and numerical agreement.
- Install local changes before loading them or installed help. Restart sessions
  that had the replaced package loaded; use private libraries for experiments.
- Keep examples reproducible with explicit seeds where needed; avoid sensitive
  data, large outputs, and temporary artifacts. Update affected docs and NEWS for
  user-visible changes. Preserve licenses/attribution and portability.

## Code and documentation

- New R: two-space indentation, snake_case, `<-`, spaces around operators/after
  commas, braces, and readable lines near 80 columns. Preserve API/S3 names and
  nearby C++ style; do not reformat unrelated code.
- Prefer clear names, focused functions, explicit inputs, predictable returns,
  and early validation with actionable errors. Handle empty/missing inputs and
  dimensions deliberately; use seq_len/seq_along and drop = FALSE as appropriate.
- Avoid hidden globals, setwd, and persistent session changes; restore temporary
  options with on.exit. Use pkg::fun or declared imports in package code, not
  library or runtime installation. Declare dependencies.
- Explain scientific intent, units, tolerances, and non-obvious choices. Remove
  dead code. In C++, use standard containers/automatic resource management and
  check bounds and R/C++ index conversions.
- Document beside R source with roxygen. Regenerate help/NAMESPACE using
  roxygen2::roxygenise and review diffs; do not hand-edit generated files.
  Regenerate Rcpp bindings with Rcpp::compileAttributes when signatures change.
- Prefer tidyverse and native `|>` in vignettes, sf for geography, sfnetworks for
  network preparation, and data.table where useful. Attach example dependencies
  explicitly and qualify ambiguous functions. Vignette-only packages go in
  Suggests; these preferences need not extend to the measurement engine.
- Preserve tibble/grouped-data/sf/data.table interoperability: explicit grouping
  and geometry behavior, no unexpected by-reference mutations, no conflicting
  exports. Exercise sfnetworks handoffs through the supported edge-list interface.

Optional reading: [Code Smells and Feels](https://github.com/jennybc/code-smells-and-feels).
