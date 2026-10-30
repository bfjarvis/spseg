#ifndef SPSEG_TYPES_H
#define SPSEG_TYPES_H

namespace spseg_spatial {
// Values are the R-to-C++ configuration codes; keep their order stable.
enum VariableType { Categorical = 0, Ordinal = 1, ContinuousIntervals = 2, ContinuousExact = 3 };
enum SummaryType { Shares = 0, Quantiles = 1 };
enum Measure {
  Exposure = 0, Information = 1, Diversity = 2,
  Dissimilarity = 3, Root = 4, MeasureCount = 5
};
enum Kernel { Unweighted = 0, Biweight = 1, Inverse = 2, Exponential = 3 };
enum Neighborhood { Radius = 0, PopulationNearest = 1 };
enum SearchBackend { Indexed = 0, Brute = 1 };
}

#endif
