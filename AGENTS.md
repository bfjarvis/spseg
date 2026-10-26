# Working on spseg

## Start with the assignment

- Read `project/GOALS.md` and your task brief. Consult `project/PLAN.md`,
  `project/TEAM.md`, and relevant decisions when needed; do not load every report.
- Treat inherited chat history as background. Check current files and decisions;
  flag contradictions rather than silently choosing a scientific interpretation.
- Inspect the working tree before editing. Preserve unrelated and in-progress work.
- Work within the assigned scope. Resolve routine choices independently; raise
  material scientific uncertainty, incompatible requirements, or scope changes.

## Coordinate changes

- Give each deliverable one owner. Coordinate overlapping edits; use separate Git
  worktrees for concurrent implementation and identify the starting revision.
- Management owns central project records; workers own assigned task reports.
  Route findings through management rather than starting competing workstreams.
- Follow the user's authorization for cross-chat messages and delegation. A role
  description alone does not authorize launching or messaging other agents.
- Keep chat links/IDs and local routing details only in ignored
  `project/TEAM.local.md`. Never force-add it or copy those details into tracked
  reports, commits, issues, or pull requests. Refer to agents by role elsewhere.
- Keep development records, reference papers, and research drafts under `project/`.
  `.Rbuildignore` excludes this directory from package artifacts; `.gitignore`
  controls Git tracking separately. Keep user documentation in `man/` and `vignettes/`.
- Share committed revisions or explicit artifacts across worktrees; files and
  conversation history do not automatically synchronize. Recheck after integration.

## Develop and verify

- Prefer small, reviewable changes and existing package conventions. Avoid unrelated
  refactoring, new dependencies, and new abstractions without a concrete need.
- Separate published measures, numerical conventions, and proposed extensions.
  Cite original sources; do not infer scientific validity from passing code tests.
- Keep population eligibility, weights, units, ties, missingness, and neighborhood
  semantics explicit. Do not silently change them to simplify implementation.
- Test changed behavior with meaningful reference cases, invariants, and regressions.
  Expected scientific results should not merely duplicate the implementation.
- Run checks appropriate to the change. Use clean-session vignette rendering for
  workflow changes and `R CMD check` for release candidates; report checks not run.
- Optimize measured bottlenecks. Record baseline/candidate revisions, workload,
  environment, runtime and memory where relevant, and numerical agreement.
- Keep examples reproducible and practical to run. Document required external data;
  do not commit sensitive data, large outputs, or temporary build artifacts.
- Update affected documentation and `NEWS.md` for user-visible changes. Preserve
  attribution and license notices; avoid unnecessary platform-specific assumptions.

## Code hygiene and style (all roles)

- Apply these conventions to package code, tests, analysis scripts, and vignette
  chunks. Prefer clear names and explicit data flow over clever or compressed code.
- In new R code, use two-space indentation, `snake_case` names, `<-` assignment,
  spaces around operators and after commas, and braces for control-flow bodies.
  Preserve established API/S3 names and nearby C++ conventions; avoid formatting
  unrelated code. Wrap long expressions for readability, usually near 80 columns.
- Keep functions focused. Extract repeated logic when it has a coherent purpose;
  avoid speculative frameworks. Use early checks/returns to reduce deep nesting.
- Pass inputs explicitly and return predictable structures. Avoid hidden globals,
  `setwd()`, and persistent session changes; restore temporary options with
  `on.exit()`. Use project-relative paths and explicit seeds in stochastic examples.
- Validate inputs at boundaries and give actionable errors. Handle empty inputs,
  missing values, and dimensions deliberately; use `seq_len()`/`seq_along()` for
  potentially empty iteration and `drop = FALSE` when dimensions must be retained.
- Use `pkg::fun()` or declared imports in package code, not `library()` or runtime
  installation. Examples may attach packages explicitly. Declare dependencies.
- Explain scientific intent, units, numerical tolerances, and non-obvious choices
  in comments. Remove dead code; use arguments rather than commented-out branches
  to control behavior. Name meaningful constants instead of unexplained literals.
- Document package functions in roxygen comments beside their R source, including
  parameters, returns, and examples. Regenerate help and `NAMESPACE` with
  `roxygen2::roxygenise()`; review generated diffs. Do not hand-edit roxygen-generated
  files. Regenerate Rcpp bindings with `Rcpp::compileAttributes()` when needed.
- In C++, prefer standard containers and automatic resource management; check
  bounds and R/C++ index conversions. Preserve numerical behavior during refactoring.

Further reading, not mandatory context for every task:
[Jenny Bryan, Code Smells and Feels](https://github.com/jennybc/code-smells-and-feels).

## Vignettes and interoperability

- Prefer tidyverse idioms and the native `|>` pipe in vignettes. Use `sf` for
  geographic workflows and `sfnetworks` for network preparation where relevant;
  use `data.table` when it offers a clear benefit for the example.
- Show realistic preparation, analysis, and plotting workflows. Attach required
  packages explicitly and qualify ambiguous functions with `pkg::fun()`.
- Keep public functions compatible with these ecosystems. Avoid conflicting export
  names; check relevant tibble, grouped-data, `sf`, and `data.table` inputs and
  extraction workflows. Make grouping/geometry behavior explicit and avoid
  unexpectedly modifying caller-owned data by reference. Verify `sfnetworks`
  handoffs through the supported edge-list interface.
- Declare vignette-only dependencies in `Suggests`; these presentation preferences
  do not require adopting tidyverse dependencies in the measurement engine.

## Report briefly

- State the outcome, changed artifacts/revision, checks and results, limitations,
  and next action. Link evidence instead of copying logs into shared records.
- Distinguish proposed, implemented, tested, and scientifically validated work.
  Do not mark an assignment complete while its acceptance criteria remain unmet.
- Keep durable rules here, priorities in `project/PLAN.md`, and detailed evidence
  in task reports. Do not repeat the same background in each document.
