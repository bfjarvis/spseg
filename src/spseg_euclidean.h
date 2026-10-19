#ifndef SPSEG_EUCLIDEAN_H
#define SPSEG_EUCLIDEAN_H

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace spseg_spatial {
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class EuclideanSearch {
  using Point = bg::model::point<double, 2, bg::cs::cartesian>;
  using Value = std::pair<Point, int>;
  using Tree = bgi::rtree<Value, bgi::rstar<32>>;
  const arma::vec& x;
  const arma::vec& y;
  std::unique_ptr<Tree> tree;

  double distance(int i, int j) const {
    const double dx = x[i] - x[j], dy = y[i] - y[j];
    return std::sqrt(dx * dx + dy * dy);
  }

public:
  EuclideanSearch(const arma::vec& x_, const arma::vec& y_, bool indexed)
    : x(x_), y(y_) {
    if (indexed) {
      std::vector<Value> values;
      values.reserve(x.n_elem);
      for (std::size_t i = 0; i < x.n_elem; ++i)
        values.emplace_back(Point(x[i], y[i]), static_cast<int>(i));
      tree.reset(new Tree(values.begin(), values.end()));
    }
  }

  std::vector<Neighbor> find(int focal, double threshold, bool knn,
                            const std::vector<const double*>& populations) const {
    std::vector<Neighbor> out;
    if (!tree) {
      out.reserve(x.n_elem);
      for (std::size_t j = 0; j < x.n_elem; ++j) {
        const double d = distance(focal, j);
        if (knn || d <= threshold) out.push_back({static_cast<int>(j), d});
      }
    } else if (!knn) {
      const double infinity = std::numeric_limits<double>::infinity();
      const bg::model::box<Point> box(
        Point(std::nextafter(x[focal] - threshold, -infinity),
              std::nextafter(y[focal] - threshold, -infinity)),
        Point(std::nextafter(x[focal] + threshold, infinity),
              std::nextafter(y[focal] + threshold, infinity)));
      for (auto it = tree->qbegin(bgi::intersects(box)); it != tree->qend(); ++it) {
        const double d = distance(focal, it->second);
        if (d <= threshold) out.push_back({it->second, d});
      }
    } else {
      std::vector<PopulationSum> sums(populations.size());
      bool focal_sufficient = true;
      for (const double* p : populations)
        if (p[focal] < threshold) focal_sufficient = false;
      if (focal_sufficient) return {{focal, 0}};
      double boundary = std::numeric_limits<double>::infinity();
      const auto query = bgi::nearest(Point(x[focal], y[focal]), x.n_elem);
      for (auto it = tree->qbegin(query); it != tree->qend(); ++it) {
        const int id = it->second;
        const double d = distance(focal, id);
        if (d > boundary) break;
        out.push_back({id, d});
        for (std::size_t g = 0; g < sums.size(); ++g) sums[g].add(populations[g][id]);
        if (populations_reached(sums, threshold)) boundary = std::min(boundary, d);
      }
    }
    order_neighbors(out, focal, knn);
    return out;
  }
};
} // namespace spseg_spatial
#endif
