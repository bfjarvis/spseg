# spseg 0.9.0.9000

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
  percentile space. Exposure is unnormalized, with default global weighted
  decile groups for observed continuous values and configurable cuts.
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
