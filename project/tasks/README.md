# Task briefs and reports

Create one `TNNN-short-name.md` file for each substantial assignment; small fixes
can stay in the plan. Management writes the brief; the assigned owner appends results.
Do not create reports for every tool call or copy entire chat histories.

Use this compact structure:

```markdown
# TNNN — Outcome
Status: proposed | active | blocked | in review | complete
Owner / reviewers:
Starting revision / worktree:

## Brief
- Objective and scope:
- Inputs: relevant files, decisions, sources
- Deliverable and acceptance criteria:
- Dependencies and edit boundaries:

## Results
- Outcome and changed artifacts / resulting revision:
- Evidence: commands, environment, observed results, report links
- Limitations, checks not run, and unresolved questions:
- Next action / handoff:
```

Identify the exact tested revision; for uncommitted work, describe the diff and
preserve the reviewed artifact. Include seeds and workload parameters when relevant.
Link large logs or local data locations instead of committing them by default.
Management records acceptance and updates the plan after required review.
