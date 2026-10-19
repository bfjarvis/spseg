# seg 1.0.0.9000

* `spseg()` accepts named categorical, ordinal, grouped-continuous, and
  observed-continuous specifications, sharing spatial discovery.
* Multiple-variable kNN selects each variable's eligible population, including
  variable-specific boundary fractions and adaptive kernel radii. Thresholds
  exceeding a variable's total population raise an error. Results report
  selected population and radius ranges; legacy focal inclusion is preserved.
* Added root-family indices: binary Hutchens, symmetric multigroup Atkinson,
  and ordinal/rank-order square-root measures. Added explicitly specified
  ordinal and rank-distribution dissimilarity supplements.
* Observed continuous variables use deterministic weighted midranks and
  empirical normalizers; grouped CDFs use analytical linear interpolation in
  percentile space. Exposure is unnormalized, with default global weighted
  decile groups for observed continuous values and configurable cuts.
* Added multi-analysis results and extraction methods, methodological source
  documentation, and numerical regression tests. The new interface initially
  supports raw spatial support; legacy grid/pycnophylactic calls are preserved.
* Corrected brute-force kNN focal inclusion when coordinates coincide.
* Stop explicitly at fractional kNN boundaries so floating-point residuals
  cannot add another unit and enlarge the kernel radius in optimized builds.
* Excluded local reference PDFs from package builds and made execution of the
  external Providence-data example opt-in for portable package checks.
