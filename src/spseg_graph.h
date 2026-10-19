#ifndef SPSEG_GRAPH_H
#define SPSEG_GRAPH_H

namespace spseg_spatial {
class GraphSearch {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
  explicit GraphSearch(Rcpp::List network);
  ~GraphSearch();
  std::vector<Neighbor> find(int focal, double threshold, bool knn,
                            const std::vector<const double*>& populations);
};
} // namespace spseg_spatial
#endif
