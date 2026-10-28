# Development plan

Status: continuous summary implementation and revised interval allocation closed
on 2026-10-28. Accepted changes and outstanding non-ignored work are included in
the authorized closeout commit. Organization simplification remains a separate task.
The user-owned getting-started vignette must remain untouched by agents.
Owner: management.

## Current work

| ID | Deliverable | Owner | Status | Completion criterion |
| --- | --- | --- | --- | --- |
| P00 | Concise project framework and agent rules | Management | Draft for review | User reviews goals, responsibilities, workflow, and initial scope |
| [T011](tasks/T011-continuous-summary-implementation.md) | Implement shared continuous summaries and binary group-pair indices | Developer | Accepted | Candidate and useful regression checks pass |
| [T012](tasks/T012-continuous-summary-validation.md) | Independent scientific/behavioral verification | Validation/testing | Accepted | Reference cases and final regression evidence accepted |
| [T013](tasks/T013-continuous-summary-science-review.md) | Scientific contract review | Science officer | Accepted | Candidate semantics and claims accepted |
| [T014](tasks/T014-continuous-summary-documentation.md) | User workflow/docs consistency | Analyst | Accepted | Changed examples work; user's getting-started file untouched |
| [T015](tasks/T015-continuous-summary-efficiency-review.md) | Bounded candidate efficiency review | Optimization | Accepted | No material accidental search/retention/pairwise regression |
| [T006](tasks/T006-install-decompression-warning.md) | Diagnose/remedy delayed warning after `pak::pak(".")` | Developer | Closed: user confirmed resolution | Evidence-based remedy; clean install/restart/load behavior verified |
| [T007](tasks/T007-install-warning-validation.md) | Independently verify installation-warning remedy | Validation/testing | Closed: verification accepted | Reproduction matrix, forced lazy-object checks, and candidate verdict |
| [T008](tasks/T008-continuous-summary-science.md) | Continuous summary and pairwise scientific contract | Science officer | In review | Explicit estimands, weighting, small distinguishing example, and limits |
| [T009](tasks/T009-continuous-summary-usability.md) | Synchronized specification and output design | Analyst | In review | Compact proposed syntax, defaults, and unambiguous output labels |
| [T010](tasks/T010-continuous-summary-implementation.md) | Current behavior and implementation implications | Developer | In review | Code-path map, compatibility risks, and validation checklist |

New design question: synchronize continuous localenv/exposure by default through a
shared summary framework. The user requests one public `probs` vector for global
groups and local quantile evaluations, with implicit grouping endpoints. Those
operations retain distinct internal meanings without independent public controls.
An alternative absolute `cuts` vector is proposed, mutually exclusive with `probs`.
The user selects quartiles as the default probability grid for both summary types.
For interval income, the user superseded boundary matching on 2026-10-28:
preserve requested quantiles through interpolation and proportional allocation,
with a warning. The revised implementation uses common within-interval global-rank
allocation for focal and target mass. T011–T015 now cover this candidate; older
matching/fallback instructions and results are historical. The central design
records the accepted contract.

Final evidence: all nine test files pass, including 48 independent reference
assertions and fractional pairwise checks. The private Apple Clang build and
`R CMD check --no-manual --ignore-vignettes` with force-suggests disabled pass
(`Status: OK`). All three analyst-owned vignettes render. Optional `sfnetworks`
was unavailable. The user-owned getting-started vignette remains unfinished and
unchanged by agents; its duplicate chunk label prevents the broader vignette
check. The local GCC16 exception abort remains unresolved; no source workaround
or global compiler change was introduced. See T012/T014 for evidence and limits.

Management accepts T011–T015 for this scoped implementation, not as a release
certification. Commit includes outstanding public project records and the user's
current getting-started draft, as authorized; ignored private routing is excluded.
The next organizational change is a separate task. Nonblocking exposure-only
pair-accumulator overhead remains a later optional optimization (T015).

Installation warning finding: controlled pak reinstalls over a loaded package
reproduced stale help/code database failures, including the exact `fetch(key)`
warning. New processes read the same installed files cleanly. Independent pak
installation and new-process verification passed all 52 bindings, six help topics,
bundled data, and a public calculation without warnings. README/NEWS document the
safe workflow; agent rules reinforce it. No package-logic or persistent-library
repair was indicated. On 2026-10-28 the user confirmed resolution and identified
stale help files as the cause. No further installation-warning work is pending.
The original delayed callback was not directly traced. See T006/T007 for evidence.

## Immediate delivery priorities

The user's requested sequence is getting started, then small income and network
examples whose calculations can be checked with tidyverse or by hand. Prioritize
these concrete workflows over a comprehensive audit or optimization program.

| Priority | Deliverable / owner | Acceptance criteria |
| --- | --- | --- |
| 1 | Finish getting-started vignette / analyst; developer fixes package defects | Preserve the user's current narrative and edits; coherent individual-to-aggregate-to-surface workflow; tidyverse and native pipes; clear outputs/interpretation; clean-session render and scientific review |
| 2 | Small income worked example / analyst with science officer | Tiny deterministic dataset; explicit ranks, local/global CDFs, weights and empirical normalizer; independently calculate at least one measure; compare mixing/sorting and equivalent individual/count inputs; validate numerical agreement |
| 3 | Small network worked example / analyst with science officer | Hand-checkable edge list, distance and membership tables; unweighted bandwidth and population-kNN examples; transparent local compositions and index calculation; compare package output with independent results |
| Parallel blocker | Network error-path repair / developer after validation handoff | Unreachable population produces an informative R error rather than terminating R; independent regression check and relevant network-suite run on the repaired revision |

These are the next delivery briefs, not yet dispatched implementation assignments.
First resume/consolidate the relevant baseline work. Analyst owns vignette source;
developer owns package fixes; validation owns independent checks; science reviews
definitions and interpretation. Preserve any user edits made since the first snapshot.

Candidate teaching datasets: eight people with four distinct income values repeated
twice across two neighborhoods, and a four-node path graph. Start with one fully
worked measure and unweighted neighborhoods. Add other measures, kernels, weighted
edges, and interval-income assumptions only after the basic calculation is clear.
Choose exact expected values with science/validation review, not from package output.
Keep vignette reference calculations independent of package internals; reuse the
same cases in regression tests. Full benchmarks and new-method development follow
these correctness and usability priorities.

## First milestone: verified baseline

Purpose: establish what exists, which claims have evidence, and which gaps matter
most before scheduling overlapping implementation work.

| Assignment | Owner | Status | Required result |
| --- | --- | --- | --- |
| [T001: Method/provenance inventory](tasks/T001-method-provenance.md) | Science officer | In review | Source-to-behavior map; supplements and unresolved network claims |
| [T002: Validation triage](tasks/T002-validation-triage.md) | Validation/testing | In review | Reproduce reported network abort; evidence inventory and gaps |
| [T003: Implementation handoff](tasks/T003-implementation-handoff.md) | Developer | In review | Preserve recent assessment and test-crash evidence; identify interface gaps |
| [T004: Performance assessment](tasks/T004-performance-assessment.md) | Optimization | In review | Existing evidence, full-output workloads, and benchmark plan |
| [T005: User workflow assessment](tasks/T005-user-workflows.md) | Analyst | In review | Geographic and small-network walkthroughs; reproducible friction reports |
| Consolidated backlog | Management | Awaiting reports | Prioritized issues with owners, dependencies, and acceptance criteria |

Initial assignments cover assessment and reports, not package changes. Reuse
existing findings, distinguish them from fresh checks, and preserve the user's
unfinished vignette. Read-only assessments can overlap; compile/run checks in
isolated copies with separate libraries and temporary outputs. Avoid competing
timed benchmarks. Record the revision and dirty-tree scope each report covers.

Priority evidence: validation's latest chat reports independent reproduction of
network `SIGABRT`, isolated to the unreachable-population error test; the other
seven network test blocks reportedly passed. The durable report is unfinished;
the underlying cause is not established. Complete the evidence handoff before repair.

The milestone is complete when findings are linked, conflicts are resolved or
explicitly recorded, and the next implementation priorities have acceptance criteria.
This is a baseline assessment, not a declaration of release readiness.

## Later work

Sequence scientific corrections, workflow improvements, network validation, and
optimization from the baseline findings. Set performance targets and a supported
platform matrix from evidence and user priorities; no release date is committed.
