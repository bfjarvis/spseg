# Development plan

Status: framework remains reviewable; user authorized initial specialist contact
and assignments on 2026-10-26. Owner: management.

## Current work

| ID | Deliverable | Owner | Status | Completion criterion |
| --- | --- | --- | --- | --- |
| P00 | Concise project framework and agent rules | Management | Draft for review | User reviews goals, responsibilities, workflow, and initial scope |

## First milestone: verified baseline

Purpose: establish what exists, which claims have evidence, and which gaps matter
most before scheduling overlapping implementation work.

| Assignment | Owner | Status | Required result |
| --- | --- | --- | --- |
| [T001: Method/provenance inventory](tasks/T001-method-provenance.md) | Science officer | Dispatched | Source-to-behavior map; supplements and unresolved network claims |
| [T002: Validation triage](tasks/T002-validation-triage.md) | Validation/testing | Dispatched | Reproduce reported network abort; evidence inventory and gaps |
| [T003: Implementation handoff](tasks/T003-implementation-handoff.md) | Developer | Dispatched | Preserve recent assessment and test-crash evidence; identify interface gaps |
| [T004: Performance assessment](tasks/T004-performance-assessment.md) | Optimization | Dispatched | Existing evidence, full-output workloads, and benchmark plan |
| [T005: User workflow assessment](tasks/T005-user-workflows.md) | Analyst | Dispatched | Geographic and small-network walkthroughs; reproducible friction reports |
| Consolidated backlog | Management | Awaiting reports | Prioritized issues with owners, dependencies, and acceptance criteria |

Initial assignments cover assessment and reports, not package changes. Reuse
existing findings, distinguish them from fresh checks, and preserve the user's
unfinished vignette. Read-only assessments can overlap; compile/run checks in
isolated copies with separate libraries and temporary outputs. Avoid competing
timed benchmarks. Record the revision and dirty-tree scope each report covers.

Priority evidence: development reports a network-test `SIGABRT` in a temporary
build. Validation owns reproduction; the underlying cause is not yet established.

The milestone is complete when findings are linked, conflicts are resolved or
explicitly recorded, and the next implementation priorities have acceptance criteria.
This is a baseline assessment, not a declaration of release readiness.

## Later work

Sequence scientific corrections, workflow improvements, network validation, and
optimization from the baseline findings. Set performance targets and a supported
platform matrix from evidence and user priorities; no release date is committed.
