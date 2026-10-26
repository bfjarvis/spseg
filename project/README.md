# Development workspace

This directory holds project coordination, research, and supporting material.
It is excluded from R package builds by the root `.Rbuildignore`. Exclusion from
the package does not make material private or prevent it from being tracked in Git.

| Material | Location | Distribution |
| --- | --- | --- |
| Goals, roles, plans, decisions, evidence | This directory and `tasks/` | May be versioned on GitHub; excluded from the R package |
| Original reference papers | `references/` | Excluded from the R package; select files for Git deliberately |
| Research manuscripts and supporting scripts | `papers/` | May be versioned on GitHub; excluded from the R package |
| Private chat identities and routing | `TEAM.local.md` | Git-ignored; never publish or force-add |
| Package help and runnable user guides | Root `man/` and `vignettes/` | Part of the package workflow; do not move here |

Start with [goals](GOALS.md), [team responsibilities](TEAM.md), and the
[current plan](PLAN.md). See [decisions](DECISIONS.md) and [progress](PROGRESS.md)
when needed; substantial assignments have individual [task reports](tasks/README.md).

## Recommended release approach

Keep development material and package sources in one repository. Use branches and
worktrees for isolated work, rather than maintaining a second copy for deployment.
Build a source package with `R CMD build`, check the resulting tarball with
`R CMD check`, and inspect its contents before release. Tag the tested revision
and distribute the checked package artifact through the chosen release channel.

A GitHub source archive or clone contains tracked development material. An R-built
package tarball respects `.Rbuildignore`. These are different artifacts; use the
latter when promising an end-user package free of the development workspace.

Release automation and full validation are still pending. A packaging-only check
that skips vignette rendering does not establish release readiness.
