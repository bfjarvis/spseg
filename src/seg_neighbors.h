#ifndef SEG_NEIGHBORS_H
#define SEG_NEIGHBORS_H

// Shared, single-threaded spatial search and kernel helpers.
namespace seg_spatial {
struct Neighbor {
  int id;
  double d;
};

struct PointCloud {
  const arma::vec& X;
  const arma::vec& Y;

  inline size_t kdtree_get_point_count() const { return X.n_elem; }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    return dim == 0 ? X[idx] : Y[idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, PointCloud>,
  PointCloud, 2, size_t>;

inline double seg_weight(double d, double bw, double power, int weighting_id)
{
  double w = 0.0;

  switch (weighting_id) {
  case 0:
    w = 1.0;
    break;
  case 1:
    w = std::pow(1.0 - std::pow(d / bw, power), power);
    if (d > bw)
      w = 0.0;
    break;
  case 2:
    w = 1.0 / std::pow(d / bw, power);
    break;
  case 3:
    w = std::exp(-d / bw);
    break;
  default:
    w = std::exp(-d / bw);
  }

  return w < 0.0 ? 0.0 : w;
}

inline bool neighbor_distance_less(const Neighbor& a, const Neighbor& b)
{
  if (a.d == b.d)
    return a.id < b.id;
  return a.d < b.d;
}

inline std::vector<Neighbor> all_neighbors(int i, const arma::vec& X,
                                    const arma::vec& Y)
{
  const int n = X.n_elem;
  std::vector<Neighbor> neighbors;
  neighbors.reserve(n);
  for (int j = 0; j < n; ++j) {
    const double dx = X[i] - X[j];
    const double dy = Y[i] - Y[j];
    neighbors.push_back({j, std::sqrt(dx * dx + dy * dy)});
  }
  std::sort(neighbors.begin(), neighbors.end(), neighbor_distance_less);
  return neighbors;
}

inline std::vector<Neighbor> radius_neighbors(const KDTree& tree, int i,
                                       const arma::vec& X,
                                       const arma::vec& Y,
                                       double radius)
{
  if (radius < 0.0)
    return all_neighbors(i, X, Y);

  double query[2] = {X[i], Y[i]};
  std::vector<nanoflann::ResultItem<size_t, double>> results;
  const double radius2 = std::nextafter(radius * radius,
                                        std::numeric_limits<double>::infinity());
  tree.radiusSearch(query, radius2, results,
                    nanoflann::SearchParameters(0.0, true));

  std::vector<Neighbor> neighbors;
  neighbors.reserve(results.size());
  for (const auto& result : results)
    neighbors.push_back({static_cast<int>(result.first),
                         std::sqrt(result.second)});

  std::sort(neighbors.begin(), neighbors.end(), neighbor_distance_less);
  return neighbors;
}

inline std::vector<Neighbor> knn_neighbors(const KDTree& tree, int i,
                                    const arma::vec& X,
                                    const arma::vec& Y,
                                    const arma::vec& T_i,
                                    double threshold)
{
  const int n = X.n_elem;
  std::vector<Neighbor> neighbors;

  if (T_i[i] >= threshold) {
    neighbors.push_back({i, 0.0});
    return neighbors;
  }

  double query[2] = {X[i], Y[i]};
  int k = std::min(n, 32);

  while (true) {
    std::vector<size_t> ids(k);
    std::vector<double> d2(k);
    const size_t found = tree.knnSearch(query, k, ids.data(), d2.data());

    neighbors.clear();
    neighbors.push_back({i, 0.0});
    for (size_t q = 0; q < found; ++q) {
      const int id = static_cast<int>(ids[q]);
      if (id == i)
        continue;
      neighbors.push_back({id, std::sqrt(d2[q])});
    }
    std::sort(neighbors.begin() + 1, neighbors.end(), neighbor_distance_less);

    double cumulative = T_i[i];
    for (std::size_t q = 1; q < neighbors.size() && cumulative < threshold; ++q)
      cumulative += T_i[neighbors[q].id];

    if (cumulative >= threshold || k >= n)
      break;
    k = std::min(n, k * 2);
  }

  return neighbors;
}


struct Contributions {
  std::vector<int> ids;
  std::vector<double> weights;
  double scale = 1.0;
};

// Scale preserves legacy radius-environment counts; it cancels from shares.
inline Contributions radius_contributions(const std::vector<Neighbor>& neighbors,
    double band, double power, int weighting, bool normalized)
{
  Contributions out;
  double bw = normalized ? band : 1.0;
  if (bw <= 0.0) bw = 1.0;
  double tolerance = 0.0;
  if (weighting == 2) {
    for (const Neighbor& neighbor : neighbors) {
      if ((band < 0.0 || neighbor.d <= band) && neighbor.d > 0.0 &&
          (tolerance <= 0.0 || neighbor.d < tolerance)) tolerance = neighbor.d;
    }
    tolerance = tolerance > 0.0 ? tolerance / 2.0 : bw;
  }
  double total = 0.0;
  for (const Neighbor& neighbor : neighbors) {
    if (band >= 0.0 && neighbor.d > band) break;
    double d = neighbor.d;
    if (weighting == 2 && d == 0.0) d = tolerance;
    const double w = seg_weight(d, bw, power, weighting);
    out.ids.push_back(neighbor.id);
    out.weights.push_back(w);
    total += w;
  }
  if (total > 0.0) out.scale = 1.0 / total;
  return out;
}

inline Contributions knn_contributions(const std::vector<Neighbor>& neighbors,
    const arma::vec& population, double threshold, double power,
    int weighting, bool normalized)
{
  Contributions out;
  if (neighbors.empty()) return out;
  const int focal = neighbors[0].id;
  out.ids.push_back(focal);
  out.weights.push_back(1.0);
  double remaining = threshold - population[focal];
  if (remaining <= 0.0) return out;
  double farthest = 0.0, tolerance = 0.0;
  for (std::size_t q = 1; q < neighbors.size() && remaining > 0.0; ++q) {
    const Neighbor& nb = neighbors[q];
    if (population[nb.id] <= 0.0) continue;
    const double fraction = std::min(1.0, remaining / population[nb.id]);
    farthest = std::max(farthest, nb.d);
    if (nb.d > 0.0 && (tolerance <= 0.0 || nb.d < tolerance)) tolerance = nb.d;
    out.ids.push_back(nb.id);
    out.weights.push_back(fraction);
    remaining -= fraction * population[nb.id];
  }
  const double bw = normalized && farthest > 0.0 ?
    farthest * (1.0 + std::sqrt(std::numeric_limits<double>::epsilon())) : 1.0;
  tolerance = tolerance > 0.0 ? tolerance / 2.0 : bw;
  // IDs are retained in the distance order, including skipped zero-mass units.
  std::size_t q = 1;
  for (std::size_t k = 1; k < out.ids.size(); ++k) {
    while (neighbors[q].id != out.ids[k]) ++q;
    double d = neighbors[q].d;
    if (weighting == 2 && d == 0.0) d = tolerance;
    out.weights[k] *= seg_weight(d, bw, power, weighting);
  }
  return out;
}

inline std::vector<Neighbor> discover_neighbors(const KDTree* tree, int i,
    const arma::vec& x, const arma::vec& y, const arma::vec& population,
    double largest, bool knn)
{
  if (tree) return knn ? knn_neighbors(*tree, i, x, y, population, largest) :
    radius_neighbors(*tree, i, x, y, largest);
  std::vector<Neighbor> out = all_neighbors(i, x, y);
  if (knn) {
    const auto focal = std::find_if(out.begin(), out.end(),
      [i](const Neighbor& nb) { return nb.id == i; });
    std::rotate(out.begin(), focal, focal + 1);
  }
  return out;
}

} // namespace seg_spatial
#endif
