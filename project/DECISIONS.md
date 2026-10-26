# Decision record

Owner: management. Record consequential choices, not routine implementation details.
Use: ID, date, status, decision, rationale/evidence, and reconsideration trigger.
Mark replaced decisions as superseded rather than erasing their rationale.

## D001 — Project coordination records

- **Date:** 2026-10-26
- **Status:** Drafting authorized; framework awaiting review.
- **Decision:** Place concise project records in `project/` and shared agent rules
  in root `AGENTS.md`; exclude both from R package builds.
- **Rationale:** Keep goals and evidence versioned beside the code without adding
  them to the distributed package. Detailed findings belong in task reports.
- **Evidence:** User requested these drafts and the `project/` directory in the
  management chat. Chat link pending.
- **Reconsider when:** Research/data work needs a distinct dependency or publication
  lifecycle, or coordination records become difficult to maintain with the package.

No scientific decisions have yet been audited or accepted through this framework.

## D002 — Keep agent identities local

- **Date:** 2026-10-26
- **Status:** Required by the user.
- **Decision:** Store chat links/IDs and local routing only in Git-ignored
  `TEAM.local.md`. Keep `TEAM.md`, plans, and task reports free of those identities.
  Never force-add or publish the local roster.
- **Rationale/evidence:** The user supplied the roster for coordination and explicitly
  prohibited tracking or uploading the file containing agent information to GitHub.
- **Reconsider when:** The user explicitly changes that instruction.

## D003 — Separate development sources from package artifacts

- **Date:** 2026-10-26
- **Status:** Directory boundary confirmed by the user; release approach recommended.
- **Decision:** Keep research and coordination under `project/`, including the moved
  `project/references/` and `project/papers/`. Exclude `project/` and `AGENTS.md`
  from R builds. Development materials may be public; private chat identities may not.
- **Rationale:** A clean installable package does not require a second deployment
  branch or repository. Git visibility and R package contents are separate controls.
- **Release approach:** Build/check a package artifact from a known revision and
  tag that revision. A GitHub source archive still includes tracked project records.
  Worktrees isolate development work; they do not define release contents.
- **Reconsider when:** A distribution channel requires a different repository layout
  or the research workspace develops an independent lifecycle.
