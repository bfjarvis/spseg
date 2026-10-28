# T007 — Independently verify installation warning remedy

Status: in review (validation recommends acceptance; original user-session call unobserved) (2026-10-26)
Owner: validation/testing. Reviewers: developer and management.
Starting revision / environment: `e2b8acac868cdca8d0119f79c5439557e6a41dff`,
canonical checkout. At inspection, T006 was still marked ready to dispatch with
no diagnosis or candidate recorded. R 4.6.0 on macOS Sequoia 15.6.1, arm64;
pak 0.9.5. `.libPaths()` contains only
`/Library/Frameworks/R.framework/Versions/4.6/Resources/library`.

## Results

### Read-only installed-package baseline

The installed package is `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/spseg`,
version `0.9.0.9000`, built for R 4.6.0 on 2026-10-26. I ran the check in a
fresh `Rscript --vanilla` process; no installation or package files were changed.

The process forced all six exported namespace bindings, loaded `spseg_data` from
the package lazy data database (100 rows by 16 columns), then called the public
`spseg()` and `spseg_categorical()` APIs on a four-row coordinate fixture with a
zero-radius band. The call returned a `spseg_result` with finite multigroup
dissimilarity 1.0. A calling warning handler recorded zero warnings, including
no `fetch(key)` / `R_decompress1` / libdeflate warning. Exit status was 0.

Reproduce with:

```sh
Rscript --vanilla /private/tmp/spseg-T007-baseline.R
```

The script and output are `/private/tmp/spseg-T007-baseline.R` and
`/private/tmp/spseg-T007-baseline.log`. The script prints R version, library
paths, package path/version/build, forced exports, dataset dimensions, the API
result, and warning count. For reference, the inspected installed lazy database
MD5s were `spseg.rdb` d0482796f59c9bf772733897aaf68fc2 and `spseg.rdx`
4b1613c436a248f2a872377cdbc2dcd2; these identify only this observation and do
not establish integrity against a known-good build.

### Limits of the original installed-package baseline

The clean-process baseline supports that the then-installed package's tested
lazy bindings, dataset, and public calculation were healthy at test time. It did
not reproduce pak's specific workflow or establish why the user's already-loaded
session later warned. The reported trigger remains `pak::pak(".")` while spseg is
loaded, accepting pak's restart recommendation and then installing version
0.9.0.9000 over the same version with dependencies kept. The delayed warning
alone is not evidence of persistent package corruption.

At the time of the initial baseline, T006 had no diagnosis or candidate, so that
check alone could not verify pak installation, restart, and reload behavior. The
candidate sequence is assessed below. No pak installation was run in the framework
library or the user's active session. No warnings were suppressed. The separate
network test failure is out of scope.

### Installed help database baseline

At the follow-up request, I separately retrieved and parsed every installed help
topic in another fresh `Rscript --vanilla` process. `tools::Rd_db("spseg")`
listed six topics; each was then retrieved through `utils::help()` and parsed
with `utils:::.getHelpFile()`. All six completed, with zero warnings and exit
status 0. The captured topic names were `spseg`, `spseg_data`,
`spseg_default_bands`, `spseg_network`, `spseg_result`, and
`spseg_variables`. This exercises the installed help database independently of
namespace bindings and lazy package data, including the help fetch path relevant
to the reported `fetch(key)` warning.

Reproduce with:

```sh
Rscript --vanilla /private/tmp/spseg-T007-help-baseline.R
```

Script and output: `/private/tmp/spseg-T007-help-baseline.R` and
`/private/tmp/spseg-T007-help-baseline.log`. This clean read does not reproduce a
help-database failure or explain a warning after pak replaces help files while
an older session is active. Developer's controlled help reinstall fixture and
candidate diagnosis remain necessary to assess that sequence.

### Private pak install and candidate verification (2026-10-26)

After T006 recorded its candidate, I independently copied the current tracked
source tree to `/private/tmp/spseg-T007-pak-check/source`, excluding compiled
artifacts, and used `/private/tmp/spseg-T007-pak-check/library` as the only pak
installation destination. The copy contains 76 tracked files; no `.o`, `.so`,
`.dylib`, `.dll`, or `.a` files were copied. Its sorted file/hash manifest is
`/private/tmp/spseg-T007-pak-check/manifest.json` (manifest digest
`d8379c40a2498f4ffb891cc6fbf7ea8dc0af719270f26a1340529bf96639fce1`). The
source revision was `e2b8acac868cdca8d0119f79c5439557e6a41dff`, with the current
working-tree README and NEWS installation-note edits included. Package code and
installed help sources were unchanged in that candidate; version stayed
`0.9.0.9000`.

The initial pak attempt without a private cache stopped before installation with
`Cannot open lock file: Operation not permitted` in pkgcache. I inspected the
installed pkgcache configuration and reran with its supported
`R_PKG_CACHE_DIR` override pointing inside my temporary root. No user cache,
framework library, or active user session was used. The successful install
process confirmed `spseg` was absent from `loadedNamespaces()` before and after
pak. Pak installed one local package, with zero R warnings and exit status 0.

Reproduce the installation and separate-process verification with:

```sh
R_PKG_CACHE_DIR=/private/tmp/spseg-T007-pak-check/cache \
  Rscript --vanilla /private/tmp/spseg-T007-pak-install.R \
  /private/tmp/spseg-T007-pak-check/source \
  /private/tmp/spseg-T007-pak-check/library
Rscript --vanilla /private/tmp/spseg-T007-pak-verify.R \
  /private/tmp/spseg-T007-pak-check/library
```

The fresh verification process loaded the exact private pak installation, forced
all 52 namespace bindings, retrieved and parsed all six installed help topics,
loaded the 100-by-16 bundled dataset, and ran the four-row public API fixture.
The `spseg_result` dissimilarity was 1.0 as expected. Captured warning count was
zero and exit status was 0. Logs and scripts are
`/private/tmp/spseg-T007-pak-install.log`,
`/private/tmp/spseg-T007-pak-verify.log`,
`/private/tmp/spseg-T007-pak-install.R`, and
`/private/tmp/spseg-T007-pak-verify.R`.

This independently verifies the README's clean-session install and new-process
load workflow on the documentation candidate. It does not independently reproduce
the stale-session warning after reinstall: Developer owns the controlled
same-session fixture. T006 now records unchanged-source and documentation-only
reinstall reproductions in private libraries, the affected database paths, and
clean reads of those exact installations in fresh processes without repair. Its
mechanism evidence shows how stale in-process lazy-load state can refer to old
DB generations; it explicitly does not claim to have traced the user's particular
Help-pane call. My independent private pak install and fresh-process checks support
the README remedy. The README also correctly advises restarting the analysis
session that had the package loaded, rather than calling `library(spseg)` again
in that session.

### Validation recommendation

Recommend accepting the documented workflow remedy for the reproduced failure
mode: install from a fresh R process, and restart any analysis process that had
spseg loaded before using the replacement. Developer reproduced stale-session
help/code database failures and warning-free fresh-process reads; my independent
pak installation and fresh-process verification passed with zero warnings. The
exact callback that produced the user's original delayed warning remains
unobserved, so the diagnosis is supported for the controlled reproduction rather
than proven for that specific GUI event. If the warning recurs after a restart,
the README's next step—to capture the warning and `find.package("spseg")` path
in the fresh session—is appropriate. No package-code or installed-library repair
is indicated by this evidence.
