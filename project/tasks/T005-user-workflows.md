# T005 — End-user workflow assessment

Status: in review
Owner: analyst. Reviewers: science officer, validation/testing, and management.
Assessed source base: `de48d12a3e9d29ca09eb289873f1877ba5582117`, with the dirty working-tree version of `vignettes/spseg-getting-started.Rmd` (SHA256 `9c79ad1d17fdd21e2f18e871917863c47a95e5556ac275a90ed0450b8f7fb0bf`) and package sources copied into a private temporary snapshot on 2026-10-26. That snapshot was installed into a temporary R library. Management committed the framework while this assessment was underway (`e2b8acac868cdca8d0119f79c5439557e6a41dff`); the working tree now contains management edits in `project/PLAN.md` and `project/PROGRESS.md`, which were not part of the tested snapshot. The user-edited vignette was read and copied, never written by this assessment.

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

The core individual-data calculation and the small supplied-network calculations
ran from the installed snapshot. The getting-started vignette did not render to
completion: its first plot chunk refers to `segdata`, while the preceding chunk
created `curves`. This is a blocking documentation defect in the assessed snapshot.
The current material also has a detached sequence of code and prose outside Rmd
chunks, and later refers to `segcurves` and `micro` without a coherent preceding
setup. Fixing the first name alone will expose more structural problems; a clean
render is required after the owner consolidates the in-progress narrative.

### Prioritized findings

1. **P1 — Getting-started stops at the first profile plot.** To reproduce, install
   the assessed source snapshot and render the copied vignette in a fresh session.
   The assessment copied `vignettes/spseg-getting-started.Rmd` to temporary space,
   changed only `N` from 5,000 to 200 and `bands` to `c(25, 100)`, then rendered
   with `rmarkdown::render(..., envir = new.env())`. Data generation, four
   simulated patterns, and the `spseg()` calls (two bands, all measures, full
   output) completed. Rendering stopped in the first plotting chunk after
   `group-seg-graph` with `object 'segdata' not found`; that prior chunk creates
   `curves`. Thus this was not a full-size performance run, and no later
   aggregation or surface chunks executed. The 5,000-person, ten-band version was
   not run. Suggested narrative repair: choose one result object name and use it
   consistently; make every code block a real R Markdown chunk; provide one
   explicit transition from results to the selected `micro` result and local
   environments. Then render the complete vignette in a clean session.

2. **P2 — The first example asks the reader to absorb many choices before a
   checkable income interpretation.** The 200-person run returned scalar indices
   for categorical group and interval income, and completed local-environment
   computation before the later render error. That establishes that this API
   path executes, not that the interval-income result is understood or
   independently validated. The main demonstration uses simulated exact income,
   decile intervals, multiple measures, four scenarios and many bands before it
   supplies no tiny hand-calculable income case. Follow the current project
   priority with eight people, four distinct incomes repeated twice and two
   neighborhoods. Show the global ranks and each neighborhood's local
   composition, calculate one named measure independently, then compare
   individual observations with equivalent count input. State whether ranks,
   tie handling, weights, and the empirical normalizer are part of that example;
   retain interval interpolation and surface assumptions for the later sections.
   No independent income index was calculated in this bounded assessment.

3. **P2 — The edge-list example is runnable and interpretable at small scale, but
   the road-preparation bridge remains unverified here.** In a separate fresh
   `Rscript --vanilla` process, the five-node supplied graph produced indices
   and environment tables for categorical counts and interval-income counts at
   hop bands 1 and 2. A second separate process ran variable-specific KNN at
   population targets 10 and 20 and returned the per-variable selected
   population/radius summaries. Both exited successfully. The vignette clearly
   discloses that nodes with no population can be transit nodes, observations
   absent from the edge list are isolated, geometry does not define network
   distance, and network results share the global reference distribution across
   components. This is useful, but readers would benefit from a tiny table of
   edges, node counts, and one hand-checked local composition before the full
   result table. The separate road-preparation chunk was not run because
   `sfnetworks` is absent from this environment; its handoff and plotting remain
   unverified. A concrete next check is a clean run with the declared optional
   dependency installed, confirming node IDs and attributes survive
   `st_network_blend()` and the supported edge-list conversion.

4. **P3 — Weighting semantics need an immediate example beside network bands.**
   The network guide correctly explains that costs set path distance while the
   kernel controls contribution, and that default biweight gives zero weight at
   the radius boundary. Consequently a one-hop biweight band can select the
   search neighborhood but contribute zero from first-order neighbors. The
   supplied example avoids this surprise with `weighting = "unweighted"`, but a
   new user may try the package default. Keep the current explanation and make
   the choice visible next to the call (or show a paired one-band comparison).
   This is an interpretation/usability recommendation, not a claim of incorrect
   measurement.

### Run provenance and limits

The isolated run used R 4.6.0 on macOS arm64, with sf 1.1.3, terra 1.9.50,
knitr 1.51, rmarkdown 2.31, dplyr 1.2.1, and ggplot2 4.0.3. The package was
installed from a private temporary copy into a private temporary library; no
canonical build products or vignette source were changed. The temporary copy
used the dirty source as captured at assessment time. Full getting-started
rendering, the detached local-environment prose after the first failure, the
polygon aggregation/surface sections, the sfnetworks road chunk, package checks,
and independent scientific validation were not completed. `data.table` adds
no clear benefit to these small teaching examples; tidyverse verbs and base
aggregation are sufficient.

Management/development next actions: repair the object/chunk flow in the
user-owned getting-started draft and run a full clean-session render; build the
small income case with science review of its hand calculation; then add the
edge/count/local-composition table to the network case and verify the optional
sfnetworks handoff in an environment with that dependency. The network
SIGABRT/error-path report remains validation's separate assignment; this report
records only successful reachable-population demonstrations and does not
reproduce the failing path.
