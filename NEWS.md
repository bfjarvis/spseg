# seg 0.9.0.9000

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
  one calculation engine, and named analyses in a `seg_result`.
* Stop explicitly at fractional kNN boundaries so floating-point residuals
  cannot add another unit and enlarge the kernel radius in optimized builds.
* Excluded local reference PDFs from package builds and made execution of the
  external Providence-data example opt-in for portable package checks.
