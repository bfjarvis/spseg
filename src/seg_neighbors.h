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

struct Contributions {
  std::vector<int> ids;
  std::vector<double> weights;
  bool knn = false;
  double target = 0.0;
  double population = 0.0;
  double radius = 0.0;
};

// Compensate accumulation even on platforms where long double is double.
// Fractional population weights must not create spurious kNN shortfalls.
struct PopulationSum {
  long double sum = 0.0L, correction = 0.0L;
  void add(long double value) {
    const long double adjusted = value - correction;
    const long double next = sum + adjusted;
    correction = (next - sum) - adjusted;
    sum = next;
  }
};

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
  for (const Neighbor& neighbor : neighbors) {
    if (band >= 0.0 && neighbor.d > band) break;
    double d = neighbor.d;
    if (weighting == 2 && d == 0.0) d = tolerance;
    const double w = seg_weight(d, bw, power, weighting);
    out.ids.push_back(neighbor.id);
    out.weights.push_back(w);
  }
  return out;
}

inline Contributions knn_contributions(const std::vector<Neighbor>& neighbors,
    const arma::vec& population, double threshold, double power,
    int weighting, bool normalized)
{
  Contributions out;
  out.knn = true;
  out.target = threshold;
  if (neighbors.empty()) return out;
  const int focal = neighbors[0].id;
  const double focal_fraction = population[focal] > 0.0 ?
    std::min(1.0, threshold / population[focal]) : 0.0;
  out.ids.push_back(focal);
  out.weights.push_back(focal_fraction);
  out.population = focal_fraction * population[focal];
  PopulationSum selected;
  selected.add(static_cast<long double>(focal_fraction) * population[focal]);
  long double exact_remaining = static_cast<long double>(threshold) - selected.sum;
  const double tolerance = 16 * std::numeric_limits<double>::epsilon() * threshold;
  if (exact_remaining <= tolerance) return out;
  double farthest = 0.0;
  double local_tolerance = 0.0;
  for (std::size_t q = 1; q < neighbors.size() &&
       exact_remaining > tolerance; ++q) {
    const Neighbor& nb = neighbors[q];
    if (population[nb.id] <= 0.0) continue;
    const double fraction = static_cast<double>(
      std::min(1.0L, exact_remaining / population[nb.id]));
    farthest = std::max(farthest, nb.d);
    if (nb.d > 0.0 && (local_tolerance <= 0.0 || nb.d < local_tolerance))
      local_tolerance = nb.d;
    out.ids.push_back(nb.id);
    out.weights.push_back(fraction);
    selected.add(static_cast<long double>(fraction) * population[nb.id]);
    exact_remaining = static_cast<long double>(threshold) - selected.sum;
  }
  out.population = static_cast<double>(selected.sum);
  out.radius = farthest;
  const double bw = normalized && farthest > 0.0 ?
    farthest * (1.0 + std::sqrt(std::numeric_limits<double>::epsilon())) : 1.0;
  local_tolerance = local_tolerance > 0.0 ? local_tolerance / 2.0 : bw;
  // IDs are retained in the distance order, including skipped zero-mass units.
  std::size_t q = 1;
  for (std::size_t k = 1; k < out.ids.size(); ++k) {
    while (neighbors[q].id != out.ids[k]) ++q;
    double d = neighbors[q].d;
    if (weighting == 2 && d == 0.0) d = local_tolerance;
    out.weights[k] *= seg_weight(d, bw, power, weighting);
  }
  return out;
}

// One expanding candidate query serves every distinct eligible-population
// vector. Include all ties at the stopping distance so row-ID ordering agrees
// with brute force even when a tie straddles the tree's candidate boundary.
inline std::vector<Neighbor> variable_knn_neighbors(const KDTree* tree, int i,
    const arma::vec& x, const arma::vec& y,
    const std::vector<const double*>& populations, double threshold)
{
  if (!tree) {
    auto out = all_neighbors(i, x, y);
    const auto focal = std::find_if(out.begin(), out.end(),
      [i](const Neighbor& nb) { return nb.id == i; });
    std::rotate(out.begin(), focal, focal + 1);
    return out;
  }
  bool focal_sufficient = true;
  for (const double* p : populations) if (p[i] < threshold) focal_sufficient = false;
  if (focal_sufficient) return {{i, 0.0}};
  const int n = x.n_elem;
  const double query[2] = {x[i], y[i]};
  int k = std::min(n, 32);
  while (true) {
    std::vector<size_t> ids(k);
    std::vector<double> distances(k);
    const size_t found = tree->knnSearch(query, k, ids.data(), distances.data());
    std::vector<Neighbor> out = {{i, 0.0}};
    for (size_t j = 0; j < found; ++j)
      if (int(ids[j]) != i) out.push_back({int(ids[j]), std::sqrt(distances[j])});
    std::sort(out.begin() + 1, out.end(), neighbor_distance_less);
    std::vector<PopulationSum> cumulative(populations.size());
    bool sufficient = false;
    double stopping_distance = 0.0;
    const double tolerance = 16 * std::numeric_limits<double>::epsilon() * threshold;
    for (const Neighbor& nb : out) {
      sufficient = true;
      for (std::size_t g = 0; g < populations.size(); ++g) {
        cumulative[g].add(populations[g][nb.id]);
        if (cumulative[g].sum + tolerance < threshold) sufficient = false;
      }
      if (sufficient) { stopping_distance = nb.d; break; }
    }
    if (sufficient || k >= n) {
      if (sufficient && k < n && out.back().d == stopping_distance) {
        // radiusSearch uses a strict squared-distance bound. Inflate the
        // distance before squaring too, so sqrt/square rounding cannot omit
        // the boundary candidate already found by knnSearch.
        out = radius_neighbors(*tree, i, x, y,
          std::nextafter(stopping_distance, std::numeric_limits<double>::infinity()));
        const auto focal = std::find_if(out.begin(), out.end(),
          [i](const Neighbor& nb) { return nb.id == i; });
        std::rotate(out.begin(), focal, focal + 1);
      }
      return out;
    }
    k = std::min(n, k * 2);
  }
}


} // namespace seg_spatial
#endif
