# spseg development workspace

Deliver efficient, portable, usable, open-source implementations of Reardon and
colleagues' spatial, multigroup, ordinal, and income segregation measures.
Support individual and aggregate data, geographic neighborhoods, and prepared
network edge lists, including road-network workflows. Scientific correctness and
clear interpretation take precedence over speed. Network preparation/topology
repair generally belongs in established external tools.

## Roles and work tracking

The human project owner sets priorities, resolves consequential choices, and
approves releases. Three chats support the work:

| Role | Responsibility |
| --- | --- |
| Developer | Implementation, tests, performance, API/help documentation, integration. |
| Scientist | Methods, estimands, assumptions, interpretation, and exploratory innovations. |
| User | Apply the package, exercise workflows, draft examples/vignettes, identify usability problems. |

Follow [AGENTS.md](AGENTS.md). One owner per deliverable; no standing management
or mandatory multi-role review pipeline. Cross-chat coordination happens only
with explicit human direction. Chat identities live only in ignored TEAM.local.md.

[GitHub Issues](https://github.com/bfjarvis/spseg/issues) is the sole task backlog
for bugs, feature requests, documentation, testing gaps, and scientific questions.
Use `gh issue list`, `gh issue view NUMBER`, and `gh issue create` from the repo.
Each issue needs a problem/use case and acceptance criterion; close with the
implementing revision and verification result. Use a few meaningful labels,
such as bug, enhancement, documentation, methods, and performance. Small steps
can be checkboxes within one issue. No parallel tasks, plans, or progress files.

## Current scientific conventions

Continuous summaries share one specification: quartile shares by default,
explicit quantile-rank profiles, and an optional full-distribution mean rank.
One probability vector governs groups and local quantile positions, with implicit
group endpoints. Exact-value cuts are an alternative for shares; ties stay together.
Grouped outputs average the same local statistics using focal population weights.
Rank profiles are descriptive and distinct from share-valued exposure matrices.

Interval inputs retain requested probabilities through proportional allocation in
global-rank space. Assume uniform rank within each supplied interval and the same
conditional distribution at every location. Allocate by rank-span overlap using
identical fractions for focal counts, target shares, and group normalizers,
resolved before surface redistribution. Warn when categories are split and retain
allocation metadata/diagnostics. Finer groups add modeled resolution, not observed
within-interval sorting information. This replaces earlier boundary matching.

Binary pairwise indices compare resolved groups using pair-conditioned shares
and focal mass within full-variable neighborhoods. Reporting grids do not change
full-resolution continuous scalar indices or neighborhood selection. Continuous
presentation uses only `summary`; pre-1.0 syntax compatibility is not maintained.
Proposed rank-pair innovations remain separate from established package functionality.
Refer to source help/tests for the implemented contract and original papers for
published definitions; passing tests alone do not establish scientific validity.

## Files and distribution

- README.md: project context, roles, and lasting conventions.
- AGENTS.md: concise working rules; root AGENTS.md points here for discovery.
- TEAM.local.md: Git-ignored private routing; never publish or force-add.
- references/: original sources, retained separately from package documentation.
- papers/: research drafts and supporting material; proposals are not automatically
  validated package methods.

Earlier reports and detailed verification evidence remain in Git history at
`06e3792`, including `project/tasks/T012-continuous-summary-validation.md`.
There is no working-tree archive to maintain. That implementation passed all nine
test files and 48 independent assertions, three vignette renders, and a constrained
R CMD check. Verification used macOS arm64/R 4.6.0 with private Apple Clang;
vignette checks were excluded and optional sfnetworks was unavailable. It did not
resolve the local GCC16 exception abort or establish Windows/Linux release coverage.

The human-owned getting-started vignette remains unfinished; agents must leave it
untouched unless assigned. Preserve other human edits as well.

`project/` and root AGENTS.md are excluded from R package builds. Git tracking is
separate: tracked development records appear in a source checkout. Keep user help
in man/ and vignettes/. Build/check a package artifact from a known revision for
release; a separate deployment tree is unnecessary.
