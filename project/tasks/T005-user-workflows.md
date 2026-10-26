# T005 — End-user workflow assessment

Status: in progress
Owner: analyst. Reviewers: science officer, validation/testing, and management.
Starting revision: `de48d12a3e9d29ca09eb289873f1877ba5582117`; canonical checkout assessed through an isolated temporary source snapshot on 2026-10-26. Existing modifications: `.Rbuildignore`, `.gitignore`, and `vignettes/spseg-getting-started.Rmd`; `AGENTS.md` and `project/` untracked. No prior analyst artifacts were available in this chat.

## Brief

- Read `AGENTS.md`, `project/GOALS.md`, and `project/TEAM.md`. Represent a sociologist
  using the public interface and documentation, not implementation knowledge.
- Walk through the getting-started and small network examples. Record research
  objectives, data preparation, confusing behavior, failures, and interpretation gaps.
- Prefer tidyverse workflows and native `|>`; use sf/sfnetworks where relevant.
  Identify a useful data.table case only if it adds practical value.
- The getting-started vignette is actively edited by the user. Assess an isolated
  snapshot; do not repair or overwrite it. Use a small reproducible example where
  a full run is expensive, stating exactly what was and was not executed.
- Validation is investigating a reported network-test abort. Run any network example
  in a separate R process and report failures without duplicating crash diagnosis.
- Deliver prioritized friction reports with reproduction steps and suggested
  vignette improvements. Edit only this report; keep trial scripts in temporary space.

## Results

Pending. Record source/build provenance, evidence, limitations, and next actions.
