#include <RcppArmadillo.h>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/queue.hpp>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>
#include "spseg_neighbors.h"
#include "spseg_graph.h"

namespace spseg_spatial {
using Graph = boost::compressed_sparse_row_graph<boost::directedS,
  boost::no_property, boost::property<boost::edge_weight_t, double>>;
using Vertex = Graph::vertex_descriptor;

// Expose the public IncidenceGraph interface without global enumeration.
// Together with compact query-local indices this avoids a graph-sized heap
// initialization inside Boost's no-init Dijkstra implementation.
struct LocalGraph { const Graph& graph; };
}
namespace boost {
template<> struct graph_traits<spseg_spatial::LocalGraph>
  : graph_traits<spseg_spatial::Graph> {
  using traversal_category = incidence_graph_tag;
};
}
namespace spseg_spatial {
inline auto out_edges(Vertex v, const LocalGraph& g) { return out_edges(v, g.graph); }
inline auto out_degree(Vertex v, const LocalGraph& g) { return out_degree(v, g.graph); }
inline auto source(Graph::edge_descriptor e, const LocalGraph& g) { return source(e, g.graph); }
inline auto target(Graph::edge_descriptor e, const LocalGraph& g) { return target(e, g.graph); }

namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();
struct SearchComplete {};

// Storage is allocated once per call; a new epoch resets only the query state.
struct Scratch {
  std::vector<std::uint64_t> stamp;
  std::vector<double> distance;
  std::vector<boost::default_color_type> color;
  std::vector<std::size_t> index;
  std::uint64_t epoch = 0;
  std::size_t count = 0;
  explicit Scratch(std::size_t n) : stamp(n, 0), distance(n), color(n), index(n) {}
  void reset() {
    if (++epoch == 0) { std::fill(stamp.begin(), stamp.end(), 0); epoch = 1; }
    count = 0;
  }
  void touch(Vertex v) {
    if (stamp[v] == epoch) return;
    stamp[v] = epoch;
    distance[v] = infinity;
    color[v] = boost::white_color;
    index[v] = count++;
  }
};

struct DistanceMap {
  using key_type = Vertex;
  using value_type = double;
  using reference = double;
  using category = boost::read_write_property_map_tag;
  Scratch& scratch;
};
inline double get(DistanceMap map, Vertex v) {
  return map.scratch.stamp[v] == map.scratch.epoch ? map.scratch.distance[v] : infinity;
}
inline void put(DistanceMap map, Vertex v, double value) {
  map.scratch.touch(v);
  map.scratch.distance[v] = value;
}

struct ColorMap {
  using key_type = Vertex;
  using value_type = boost::default_color_type;
  using reference = value_type;
  using category = boost::read_write_property_map_tag;
  Scratch& scratch;
};
inline ColorMap::value_type get(ColorMap map, Vertex v) {
  return map.scratch.stamp[v] == map.scratch.epoch ? map.scratch.color[v] : boost::white_color;
}
inline void put(ColorMap map, Vertex v, ColorMap::value_type value) {
  map.scratch.touch(v);
  map.scratch.color[v] = value;
}

struct LocalIndex {
  using key_type = Vertex;
  using value_type = std::size_t;
  using reference = value_type;
  using category = boost::readable_property_map_tag;
  Scratch& scratch;
};
inline std::size_t get(LocalIndex map, Vertex v) {
  map.scratch.touch(v);
  return map.scratch.index[v];
}

struct Query {
  Scratch& scratch;
  const std::vector<std::vector<int>>& observations;
  const std::vector<const double*>& populations;
  std::vector<PopulationSum> cumulative;
  std::vector<Neighbor> neighbors;
  double threshold, boundary;
  bool knn;
  std::size_t examined = 0;

  Query(Scratch& s, const std::vector<std::vector<int>>& rows,
        const std::vector<const double*>& mass, double target, bool adaptive)
    : scratch(s), observations(rows), populations(mass), cumulative(mass.size()),
      threshold(target), boundary(adaptive ? infinity : target), knn(adaptive) {}

  void examine(Vertex v) {
    if ((examined++ % 256) == 0) Rcpp::checkUserInterrupt();
    const double d = get(DistanceMap{scratch}, v);
    if (d > boundary) throw SearchComplete{};
    for (int row : observations[v]) {
      neighbors.push_back({row, d});
      if (knn) for (std::size_t g = 0; g < populations.size(); ++g)
        cumulative[g].add(populations[g][row]);
    }
    if (knn && populations_reached(cumulative, threshold)) boundary = std::min(boundary, d);
  }
};

struct DijkstraVisitor : boost::default_dijkstra_visitor {
  Query& query;
  explicit DijkstraVisitor(Query& q) : query(q) {}
  template<class G> void examine_vertex(Vertex v, const G&) const { query.examine(v); }
};

struct BFSVisitor : boost::default_bfs_visitor {
  Query& query;
  explicit BFSVisitor(Query& q) : query(q) {}
  template<class G> void examine_vertex(Vertex v, const G&) const { query.examine(v); }
  template<class Edge, class G> void tree_edge(Edge e, const G& g) const {
    DistanceMap map{query.scratch};
    put(map, target(e, g), get(map, source(e, g)) + 1);
  }
};
} // namespace

struct GraphSearch::Impl {
  Graph graph;
  Rcpp::IntegerVector node;
  bool weighted;
  std::vector<std::vector<int>> observations;
  Scratch scratch;

  static Graph build(Rcpp::List network) {
    Rcpp::IntegerVector from = network["from"], to = network["to"];
    Rcpp::NumericVector cost = network["cost"];
    const bool directed = Rcpp::as<bool>(network["directed"]);
    const std::size_t n = Rcpp::as<Rcpp::CharacterVector>(network["labels"]).size();
    std::vector<std::pair<Vertex, Vertex>> edges;
    std::vector<double> costs;
    edges.reserve(from.size() * (directed ? 1 : 2));
    costs.reserve(edges.capacity());
    for (R_xlen_t e = 0; e < from.size(); ++e) {
      edges.emplace_back(from[e], to[e]); costs.push_back(cost[e]);
      if (!directed && from[e] != to[e]) {
        edges.emplace_back(to[e], from[e]); costs.push_back(cost[e]);
      }
    }
    return Graph(boost::edges_are_unsorted_multi_pass, edges.begin(), edges.end(),
                 costs.begin(), n);
  }

  explicit Impl(Rcpp::List network) : graph(build(network)), node(network["node"]),
    weighted(Rcpp::as<bool>(network["weighted"])), observations(num_vertices(graph)),
    scratch(num_vertices(graph)) {
    for (R_xlen_t row = 0; row < node.size(); ++row) observations[node[row]].push_back(row);
  }

  std::vector<Neighbor> find(int focal, double threshold, bool knn,
                            const std::vector<const double*>& populations) {
    if (knn) {
      bool enough = true;
      for (const double* p : populations) if (p[focal] < threshold) enough = false;
      if (enough) return {{focal, 0}};
    }
    scratch.reset();
    DistanceMap distance{scratch};
    const Vertex origin = node[focal];
    put(distance, origin, 0);
    Query query(scratch, observations, populations, threshold, knn);
    try {
      if (weighted) {
        const LocalGraph local{graph};
        boost::dijkstra_shortest_paths_no_color_map_no_init(local, origin,
          boost::dummy_property_map{}, distance, get(boost::edge_weight, graph),
          LocalIndex{scratch}, std::less<double>{}, std::plus<double>{},
          infinity, 0.0, DijkstraVisitor(query));
      } else {
        boost::queue<Vertex> queue;
        boost::breadth_first_visit(graph, origin, queue, BFSVisitor(query), ColorMap{scratch});
      }
    } catch (const SearchComplete&) {}
    order_neighbors(query.neighbors, focal, knn);
    return std::move(query.neighbors);
  }
};

GraphSearch::GraphSearch(Rcpp::List network) : impl(new Impl(network)) {}
GraphSearch::~GraphSearch() = default;
std::vector<Neighbor> GraphSearch::find(int focal, double threshold, bool knn,
    const std::vector<const double*>& populations) {
  return impl->find(focal, threshold, knn, populations);
}
} // namespace spseg_spatial
