# spseg 0.9.0.9000

- Added an R/hexSticker logo recipe that spells `spseg` through spatial sorting
  of continuous resources with viridis colors. Run `project/create-logo.R` to
  regenerate `man/figures/logo.png` for the GitHub README; this development-only
  script is excluded from package builds.

- Data-frame and sf extraction now returns compact identifier/result tables
  with one `value` column, omitting repeated calculation metadata and parallel
  share/rank/exposure columns. Columns follow selected summary specifications,
  not observed missingness. `what = "exposure"` continues to return either
  share exposure probabilities or grouped rank profiles, with explicit
  `statistic`/`prob` fields for ranks. Inspect calculation status and kNN
  diagnostics with `what = "diagnostics"`, and retained neighborhood interpolation
  diagnostics with `what = "interpolation"`.

- Rewrote the README as a first-use guide with GitHub installation, runnable
  examples, result extraction, count data, and network neighborhoods.

- Revised the getting-started vignette around a single microdata pattern,
  count-preserving Voronoi aggregation, non-spatial and spatial analyses, and
  portable parallel comparisons across patterns using tidyverse workflows.

- Simplified calculation orchestration, surface construction, and result
  extraction. Continuous variables now accept only `summary`; removed the
  deprecated `localenv`/`exposure` interface and its separate calculation paths.
  Ordinal `localenv` remains supported. Surface options no longer accept
  `keep_geometry` or singular geometry aliases; use `surface_geometry` with
  `"none"`, `"points"`, or `"polygons"`.
- Count columns are prepared once for population inference and compilation.
  C++ configuration codes now have descriptive names. Pycnophylactic smoothing
  reuses iteration buffers and precomputes fixed zone cell counts.

- `spseg(variables = ...)` now accepts column-name shorthand for a categorical
  variable. One column selects observed memberships; multiple columns select
  category counts. A message states the assumption and the result name.

- Continuous variables now default to synchronized quartile shares through one
  `summary` specification. Explicit quantiles use the same probability grid for
  focal groups and local rank evaluations, with optional full-distribution mean
  rank. Interval categories are proportionally apportioned across the requested
  grid using a shared linear global-rank model for focal/neighbor populations,
  replacing the earlier whole-interval matching rule. One warning flags the
  assumption; per-neighborhood output records bottom/top category shares and
  the local share requiring apportionment, with allocation fractions in metadata.
  Grouped extraction remains `what = "exposure"` with typed share/rank rows;
  rank profiles are stored separately from P. Continuous binary pairwise indices
  compare resolved global groups within original neighborhoods, while scalar
  continuous indices remain full-resolution. Continuous presentation uses only
  `summary`, replacing the former `localenv`/`exposure` settings.

- Documented a fresh-session workflow for local installation and restarting after
  replacing a loaded package, avoiding stale lazy-load and installed-help state.

- Consolidated variable specifications into `spseg_categorical()`,
  `spseg_ordinal()`, and `spseg_continuous()`. Continuous specifications now
  select `representation = "exact"` (default) or `"intervals"`, including
  interval membership and ordered interval counts. Removed `spseg_grouped()`
  without an alias. Continuous results report their representation separately
  from the variable type; calculation methods and compact environments are unchanged.

- Retained local environments default to category shares for categorical and
  ordinal variables, and synchronized global-group shares for either continuous
  representation. Ordinal constructors accept `localenv`; continuous constructors
  use `summary` for shares or compact weighted rank quantiles and optional mean
  rank. Full-resolution scalar indices use the complete available distributions;
  exact local distributions are not duplicated in the result. Extraction reports
  typed `prob`/`rank` or `category`/`share` rows.

- Added a synthetic-population overview covering individual, regular-grid,
  Voronoi, and pycnophylactic workflows. The older overview is retained and
  marked as a legacy introductory guide.

- Euclidean searches now use Boost.Geometry through BH; the vendored nanoflann
  backend has been removed. Use `search = "boost"` (the default) or `"brute"`.
- Population kNN includes equidistant boundary units proportionately after the
  focal unit, avoiding row-order dependence. Radius searches include their closed
  boundary, including cases affected by square-root rounding.
- Supplied networks now support shortest-path hop or nonnegative-cost neighborhoods
  through Boost.Graph, with optional geometry, transit nodes, and variable-specific
  population KNN. Added `spseg_network()` and a network/road-preparation vignette.

* Renamed the package to `spseg`, with `spseg_*` constructors, result classes,
  default-band helper, and native bindings. Examples and internal helpers use
  the same naming convention; no compatibility aliases are retained.
* Pre-1.0 development: `spseg()` requires named categorical, ordinal, grouped-continuous, and
  observed-continuous specifications, sharing spatial discovery.
* Multiple-variable kNN selects each variable's eligible population, including
  variable-specific boundary fractions and adaptive kernel radii. Thresholds
  exceeding a variable's total population raise an error. Results report
  selected population and radius ranges, including fractional focal inclusion.
* Added root-family indices: binary Hutchens, symmetric multigroup Atkinson,
  and ordinal/rank-order square-root measures. Added explicitly specified
  ordinal and rank-distribution dissimilarity supplements.
* Observed continuous variables use deterministic weighted midranks and
  empirical normalizers; grouped CDFs use analytical linear interpolation in
  percentile space. Share exposure is unnormalized, with default global weighted
  quartile groups, retaining the requested probabilities through interpolation
  for interval input. Exact shares can alternatively use absolute income cuts.
* Added multi-analysis results and extraction methods, methodological source
  documentation, and numerical regression tests. The new interface initially
  supports raw support and shared grid/pycnophylactic construction for count specifications.
* Removed the previous calling convention, duplicate result methods, and
  categorical calculation engine. Every call uses variable specifications,
  one calculation engine, and named analyses in a `spseg_result`.
* Removed the pre-fork concentration-profile, standalone dissimilarity, and
  spatial-proximity functions, their adjustment helpers, the old data validator,
  and unused test fixtures. Default bands use the current coordinate/population
  interface; obsolete graphics imports and the `spdep` suggestion were removed.
* Stop explicitly at fractional kNN boundaries so floating-point residuals
  cannot add another unit and enlarge the kernel radius in optimized builds.
* Excluded local reference PDFs from package builds and made execution of the
  external Providence-data example opt-in for portable package checks.
