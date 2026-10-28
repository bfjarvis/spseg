# T006 — Diagnose and remedy delayed installation decompression warning

Status: in review
Owner: developer. Reviewers: validation/testing and management.
Starting revision / environment: record current source, R/pak versions, library
paths, and relevant user edits before work. Preserve all unrelated changes.

## Reported behavior

The user runs `pak::pak(".")` with spseg already loaded. Pak warns that a restart
is probably needed; the user continues. It builds and installs spseg successfully,
then after a short delay R reports:

```
Warning message:
In fetch(key) : internal error 1 in R_decompress1 with libdeflate
```

Pak reports one local spseg update, version 0.9.0.9000 to the same version,
with dependencies kept. A fresh `Rscript --vanilla` running `library(spseg)`
completed without warnings during management's initial check. This does not
establish integrity of every lazy object or explain the user's active session.

## Work plan

1. Identify the database/package involved and capture warning calls/trace where
   possible. Inspect R/pak primary sources for the relevant mechanism. Separate
   stale in-session lazy-load state from persistent database corruption, cache
   trouble, and an actual package defect.
2. Reproduce safely with temporary source copies and private libraries: fresh
   installation/loading, reinstalling an already-loaded package in the same R
   session, and loading after process restart. Force relevant lazy bindings/help
   or objects to expose deferred failures. Use the supplied pak trigger.
3. Apply the narrowest supported remedy. A development workflow/documentation fix
   is appropriate if clean installs are healthy and replacing a loaded database
   explains the warning. Do not suppress warnings, change compression blindly,
   reinstall unrelated dependencies, or bundle the separate network-crash fix.
4. Provide the exact candidate artifacts and verification recipe to management
   and validation through this report. Add a focused regression check only if
   there is a package defect with meaningful behavior to protect.

## Scope and acceptance

- Own this report and necessary installation/workflow documentation or narrowly
  justified source changes. Coordinate any shared source edits through management.
  Do not edit vignettes or other agents' reports. Do not commit or push.
- Preserve the user's active R session and default/system library. Use isolated
  reproduction and remediation first. If a persistent repair is actually needed,
  identify the affected installation and propose its exact replacement with backup
  before changing it. The framework R library is outside workspace write roots.
- Acceptance: a supported causal explanation (or clearly delimited remaining
  uncertainty), a safe actionable remedy, and clean-session verification of the
  candidate. Document whether an existing user session must be restarted.
- Keep this bounded: use the existing evidence and produce a handoff rather than
  broad dependency/toolchain maintenance. Keep current model settings.

## Results

### Outcome and candidate remedy

Completed 2026-10-26. Reinstalling unchanged spseg sources while the package was
loaded reproduced the libdeflate warning. A documentation-only variation also
reproduced the exact `In fetch(key)` call and identified the private installed
`help/spseg.rdb` in the warning frames. The same installations passed fresh-process
checks without reinstalling or repairing any file. This supports stale in-process
lazy-load/cache state rather than persistent corruption for these experiments.

Implemented only an installation-workflow remedy in `README.md` and a matching
`NEWS.md` entry: install from a fresh R process before loading spseg/help, and
restart any analysis session that had spseg loaded before using the replacement.
That restart is still required when installation ran in a separate process.
Calling `library(spseg)` again is not a reload. No compression changes, warning
suppression, dependency updates, or measurement-code changes were made.

### Provenance and isolation

- Assessed HEAD: `e2b8acac868cdca8d0119f79c5439557e6a41dff`, branch
  `codex/relational-neighborhoods`. Unrelated project/report edits were preserved;
  T003 and vignettes were not edited. No commit or push.
- Experiment root: `/private/tmp/spseg-t006-ZNHyTd`, with `source/`,
  `changed-source/`, `lib/`, `cache/`, and `driver.R`. Only the private library was
  an installation destination. Dependencies were reused read-only. The active
  user session and framework library were not changed.
- Copied `DESCRIPTION`, `NAMESPACE`, `R`, `src`, `data`, `inst`, and `man` from
  canonical sources with R's `file.copy(..., recursive = TRUE)`, removing copied
  `.o`, `.so`, and `.dll` files before building. No project records or vignettes
  entered the experiment. The variation added only a temporary help topic
  `changed-source/man/a_t006_probe.Rd`; package code/version were unchanged.
- R 4.6.0 (2026-04-24), `aarch64-apple-darwin23`; pak 0.9.5; sf 1.1.3;
  terra 1.9.50; Rcpp 1.1.2; RcppArmadillo 15.6.2.1; BH 1.90.0.1.

### Executed reproduction matrix

| Step | Observed result | Evidence under experiment root |
| --- | --- | --- |
| Fresh process installs original source with pak | One local package updated; source build and installation successful | `install.log` |
| Fresh process loads original installation | All 52 namespace bindings, six help topics, bundled data, and public calculation passed without warnings | `verify-before.log` |
| Load original package/help, then reinstall unchanged sources in that process | Pak's loaded-package warning; installation successful; subsequent code/help reads gave libdeflate warnings and lazy-load errors | `unchanged-reinstall.log` |
| Fresh process reads that exact installation, without repair | All binding/help/data/calculation checks passed without warnings | `verify-after-unchanged.log` |
| Loaded package/help, then install same-version documentation variation | Exact `fetch(key)` warning from help DB; deferred code reads also failed; four warnings/errors captured overall | `changed-reinstall.log` |
| Fresh process reads the variation, without repair | Single help fetch, all 52 bindings, seven help topics, data and calculation passed; zero warnings/errors | `verify-after-changed.log` |

The calculation uses two distinct points with pure groups `(10,0)` and `(0,10)`,
`bands = 0`, `weighting = "unweighted"`, `measures = "information"`, and
`scope = "multigroup"`; expected H is 1. The driver records warnings through
`withCallingHandlers()` without muffling them, records errors, and extracts
`datafile`/`filebase` frame values. Verification mode asserts zero warnings/errors.
Expected-failure modes catch errors to allow remaining probes, so their exit code
alone is not a success criterion; consult the per-probe results.

Help index/database MD5 values changed even on the unchanged-source rebuild.
For the documentation variation they were `c34b42eb9f2792bb40a28d6f020e2b38`
(`help/spseg.rdx`) and `3e664985d032d0ae01ebada80682c59b`
(`help/spseg.rdb`), unchanged by the clean-process verification. Thus the healthy
restart check used exactly the installed files that failed in the prior process.

### Reproduction recipe for Validation

Use an independent temporary copy/library rather than modifying this experiment
or the user's installed package. The retained local driver can be inspected at
`/private/tmp/spseg-t006-ZNHyTd/driver.R`; paths/logs are temporary, not committed
fixtures. From `source/`, the commands executed were:

```sh
Rscript --vanilla ../driver.R install
Rscript --vanilla ../driver.R verify
Rscript --vanilla ../driver.R loaded-reinstall
Rscript --vanilla ../driver.R verify
```

The second pair of modes ran from `changed-source/` with logs listed above. In
each child R process the driver prepends its private library to `.libPaths()`
and points `PKG_CACHE_DIR`, `PKG_PACKAGE_CACHE_DIR`, and `PKG_METADATA_CACHE_DIR`
into `cache/`; `PKG_SYSREQS=false`. Pak installs only the local package:

```r
pak::pak(".", lib = .libPaths(), ask = FALSE, dependencies = FALSE)
```

The essential same-process reproduction, after a successful private install:

```r
.libPaths(c(private_lib, .libPaths()))
library(spseg, lib.loc = private_lib)
tools::Rd_db("spseg", lib.loc = private_lib) # Prime this process's old help cache.
pak::pak(".", lib = .libPaths(), ask = FALSE, dependencies = FALSE)
tools:::fetchRdDB(file.path(private_lib, "spseg", "help", "spseg"), key = "spseg")
```

Run the final read again in a new `Rscript --vanilla` process before changing the
installation. Use the driver's `verify` mode for all bindings/help/data and H.
The internal `fetchRdDB` call is diagnostic, not recommended public user API.

### Mechanism and source evidence

Installed R's `tools:::fetchRdDB()` creates `fetch(key)` around
`lazyLoadDBfetch()`. In the controlled warning, `datafile` was
`/private/tmp/spseg-t006-ZNHyTd/lib/spseg/help/spseg.rdb`. Deferred namespace reads
also named `lib/spseg/R/spseg.rdb`. The user's actual delayed call was not traced;
the matching help-call mechanism is demonstrated, not proof of their specific
database path or an RStudio Help-pane callback.

R's [help DB reader](https://github.com/wch/r-source/blob/R-4-6-branch/src/library/tools/R/Rd.R)
obtains offsets from the index and fetches compressed data. The
[R 4.6 lazy-load implementation](https://github.com/wch/r-source/blob/R-4-6-branch/src/main/serialize.c)
caches database bytes by pathname; its fetch operation reads indexed bytes and
decompresses them. Replacing files while retaining old in-process state can mix
database generations. Our same-version builds changed their DB checksums, so
unchanged version numbers are not evidence of interchangeable DB offsets.

Pak's [installation API](https://pak.r-lib.org/reference/pkg_install.html) accepts
`.` and explicit library paths; its [configuration](https://pak.r-lib.org/reference/pak-config.html)
supports the private caches above. Installed pak's `pkg_install()` and
`loaded_packages()` were inspected directly: installation takes place through a
subprocess while the calling process has its own namespace state. Pak warns about
loaded packages; successful installation does not make that old state current.

### Verification, limits, and next action

- Candidate files: `README.md`, `NEWS.md`, and this report. The evidence supports
  documentation/workflow remediation, not a package logic patch or persistent
  library replacement. No user-session restart was performed on their behalf.
- Clean install and restart/load checks above verify the operational remedy;
  independent T007 verification remains owned by Validation. The active user's
  exact delayed warning and GUI callback were not inspected. If a warning survives
  a genuinely fresh process, obtain the affected DB/library path before proposing
  a backed-up targeted replacement.
- No broad suite, native network-error test, benchmark, vignette render, Windows
  test, dependency reinstall, or `R CMD check` was run. Roxygen/native interfaces
  were unchanged, so no generated-file regeneration was needed.
