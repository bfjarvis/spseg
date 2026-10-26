# Team and ownership

Status: specialist chats identified; initial coordination authorized on 2026-10-26.
Chat links and local routing are stored only in ignored `TEAM.local.md`.
This file contains shareable role descriptions, not private chat identities.

The user sets scientific direction and priorities. Management coordinates the
specialists below; specialists contribute evidence within their assigned scopes.

| Role | Primary ownership |
| --- | --- |
| Management | Scope, assignments, dependencies, central records, acceptance and integration planning |
| Science officer | Method definitions, provenance, assumptions, interpretation; methods vignette drafts |
| Developer | Implementation, package maintenance, API documentation, coordinated integration |
| Validation/testing | Independent correctness evidence, regression checks, portability, clean-session workflow verification |
| Optimization | Profiles, representative benchmarks, measured improvements with numerical agreement |
| Analyst | End-user analyses, usability findings, workflow vignette drafts |

## Working relationships

- Management assigns one owner per deliverable and identifies needed reviewers.
- The analyst owns workflow vignette drafts; science reviews interpretation,
  development supports executable examples, and validation checks reproducibility.
- Development supplies implementation tests; validation independently challenges
  behavior. Optimization coordinates changes to code owned by development.
- The user resolves consequential scientific and scope choices and release decisions.
  Routine work proceeds within the authorized assignment without repeated approval.
- Role definitions are responsibilities, not evidence of completed work. Management
  sends authorized assignments and reads results in the specialist chats; workers
  need not send cross-chat replies. See `PLAN.md` for assignment status.
