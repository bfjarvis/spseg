# T002 — Validation inventory and network-test triage

Status: in progress (2026-10-26)
Owner: validation/testing. Reviewers: developer and management.
Starting revision / worktree: `de48d12a3e9d29ca09eb289873f1877ba5582117`,
canonical checkout; pre-existing edits in `.Rbuildignore`, `.gitignore`, and
`vignettes/spseg-getting-started.Rmd`; untracked `AGENTS.md` and `project/`.

## Brief

- Read `AGENTS.md`, `project/GOALS.md`, and the current tests. Consult the developer's
  T003 handoff when available; do not wait to inventory existing evidence.
- First priority: independently reproduce the reported network-test `SIGABRT`.
  The developer reports six non-network test files passed in a temporary build;
  neither that result nor a suspected exception-handling cause is yet verified here.
- Build/run only in an isolated copy and private temporary library. Identify the
  smallest failing test or example, commands, environment, and source snapshot.
- Inventory reference cases, regressions, and scientific/portability gaps. Separate
  successful assertions from untested claims and process crashes.
- Deliver a concise report and reproducible failure recipe, or clear evidence that
  reproduction failed. Do not repair package code during this assignment.
- Edit only this report; temporary diagnostic scripts/logs stay outside the checkout.

## Results

Pending. Record exact commands/results, evidence locations, limitations, and next steps.
