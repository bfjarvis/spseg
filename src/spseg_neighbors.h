#ifndef SPSEG_NEIGHBORS_H
#define SPSEG_NEIGHBORS_H

#include "spseg_types.h"

// Shared selection and weighting, independent of how distances were found.
namespace spseg_spatial {
struct Neighbor { int id; double d; };

inline bool neighbor_distance_less(const Neighbor& a, const Neighbor& b) {
  return a.d == b.d ? a.id < b.id : a.d < b.d;
}

inline void order_neighbors(std::vector<Neighbor>& neighbors, int focal,
                            bool focal_first) {
  std::sort(neighbors.begin(), neighbors.end(), neighbor_distance_less);
  if (focal_first) {
    auto it = std::find_if(neighbors.begin(), neighbors.end(),
      [focal](const Neighbor& nb) { return nb.id == focal; });
    if (it != neighbors.end()) std::rotate(neighbors.begin(), it, it + 1);
  }
}

inline double spseg_weight(double d, double bw, double power, int weighting) {
  const double scaled = d / bw;
  switch (weighting) {
  case Unweighted: return 1.0;
  case Biweight: return d > bw ? 0.0 :
    std::max(0.0, std::pow(1.0 - std::pow(scaled, power), power));
  case Inverse: return 1.0 / std::pow(scaled, power);
  default: return std::exp(-scaled);
  }
}

struct Contributions {
  std::vector<int> ids;
  std::vector<double> weights;
  bool knn = false;
  double target = 0.0, population = 0.0, radius = 0.0;
};

// Compensated sums also work where long double has double precision.
struct PopulationSum {
  long double sum = 0.0L, correction = 0.0L;
  void add(long double value) {
    const long double adjusted = value - correction;
    const long double next = sum + adjusted;
    correction = (next - sum) - adjusted;
    sum = next;
  }
};

inline bool populations_reached(const std::vector<PopulationSum>& sums,
                                 double target) {
  const double tolerance = 16 * std::numeric_limits<double>::epsilon() * target;
  return std::all_of(sums.begin(), sums.end(), [=](const PopulationSum& value) {
    return value.sum + tolerance >= target;
  });
}

inline Contributions radius_contributions(const std::vector<Neighbor>& neighbors,
    double band, double power, int weighting, bool normalized) {
  Contributions out;
  const double bw = normalized && band > 0 ? band : 1.0;
  double tolerance = bw;
  if (weighting == Inverse) for (const auto& nb : neighbors) {
    if (nb.d > band) break;
    if (nb.d > 0) { tolerance = nb.d / 2; break; }
  }
  for (const auto& nb : neighbors) {
    if (nb.d > band) break;
    const double d = weighting == Inverse && nb.d == 0 ? tolerance : nb.d;
    out.ids.push_back(nb.id);
    out.weights.push_back(spseg_weight(d, bw, power, weighting));
  }
  return out;
}

inline Contributions knn_contributions(const std::vector<Neighbor>& neighbors,
    const arma::vec& mass, double target, double power,
    int weighting, bool normalized) {
  Contributions out;
  out.knn = true;
  out.target = target;
  if (neighbors.empty()) return out;
  const int focal = neighbors.front().id;
  const double fraction = mass[focal] > 0 ? std::min(1.0, target / mass[focal]) : 0;
  out.ids.push_back(focal);
  out.weights.push_back(fraction);
  std::vector<double> distances{0};
  PopulationSum selected;
  selected.add(static_cast<long double>(fraction) * mass[focal]);
  const double tolerance = 16 * std::numeric_limits<double>::epsilon() * target;
  double nearest_positive = 0;
  std::size_t begin = 1;
  while (target - selected.sum > tolerance && begin < neighbors.size()) {
    std::size_t end = begin + 1;
    while (end < neighbors.size() && neighbors[end].d == neighbors[begin].d) ++end;
    PopulationSum shell;
    for (std::size_t j = begin; j < end; ++j) shell.add(mass[neighbors[j].id]);
    if (shell.sum > 0) {
      const double inclusion = static_cast<double>(
        std::min(1.0L, (target - selected.sum) / shell.sum));
      out.radius = neighbors[begin].d;
      if (out.radius > 0 && nearest_positive == 0) nearest_positive = out.radius;
      for (std::size_t j = begin; j < end; ++j) {
        const int id = neighbors[j].id;
        if (mass[id] <= 0) continue;
        out.ids.push_back(id);
        out.weights.push_back(inclusion);
        distances.push_back(neighbors[j].d);
        selected.add(static_cast<long double>(inclusion) * mass[id]);
      }
    }
    begin = end;
  }
  out.population = static_cast<double>(selected.sum);
  const double bw = normalized && out.radius > 0 ?
    out.radius * (1 + std::sqrt(std::numeric_limits<double>::epsilon())) : 1.0;
  const double zero_distance = nearest_positive > 0 ? nearest_positive / 2 : bw;
  for (std::size_t j = 1; j < out.ids.size(); ++j) {
    const double d = weighting == Inverse && distances[j] == 0 ? zero_distance : distances[j];
    out.weights[j] *= spseg_weight(d, bw, power, weighting);
  }
  return out;
}
} // namespace spseg_spatial
#endif
