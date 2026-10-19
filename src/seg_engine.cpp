#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// Vendored nanoflann provenance, license, and checksum are recorded in
// inst/NOTICE.md.
#include "nanoflann.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

using namespace Rcpp;
#include "seg_neighbors.h"
using namespace seg_spatial;

namespace {

struct GroupPair {
  int a;
  int b;
  double T_ab;
  double P_a;
  double P_b;
  double E_ab;
  double I_ab;
  bool valid;
};

double entropy(const arma::vec& p_m, double log_base)
{
  double out = 0.0;
  for (arma::uword m = 0; m < p_m.n_elem; ++m) {
    if (p_m[m] > 0.0)
      out -= p_m[m] * std::log(p_m[m]) / log_base;
  }
  return out;
}

double interaction(const arma::vec& p_m)
{
  return arma::accu(p_m % (1.0 - p_m));
}

void normalize_composition(const arma::vec& L_m, arma::vec& p_m)
{
  const double L = arma::accu(L_m);
  if (L > 0.0)
    p_m = L_m / L;
  else
    p_m.zeros();
}

std::vector<GroupPair> make_pairs(const arma::vec& T_m)
{
  const int M = T_m.size();
  std::vector<GroupPair> pairs;
  pairs.reserve(M * (M - 1) / 2);
  const double log2 = std::log(2.0);

  for (int a = 0; a < M - 1; ++a) {
    for (int b = a + 1; b < M; ++b) {
      GroupPair pair;
      pair.a = a;
      pair.b = b;
      pair.T_ab = T_m[a] + T_m[b];
      pair.valid = T_m[a] > 0.0 && T_m[b] > 0.0;
      pair.P_a = pair.valid ? T_m[a] / pair.T_ab : NA_REAL;
      pair.P_b = pair.valid ? T_m[b] / pair.T_ab : NA_REAL;
      pair.E_ab = pair.valid ?
        -(pair.P_a * std::log(pair.P_a) + pair.P_b * std::log(pair.P_b)) /
          log2 : NA_REAL;
      pair.I_ab = pair.valid ?
        pair.P_a * (1.0 - pair.P_a) + pair.P_b * (1.0 - pair.P_b) :
        NA_REAL;
      pairs.push_back(pair);
    }
  }

  return pairs;
}

// Nomenclature follows Reardon and O'Sullivan (2004)
// P: exposure/isolation
void accumulate_P_multigroup(arma::cube& P_acc, int b, int i,
                          const arma::mat& N,
                          const arma::vec& T_m,
                          const arma::vec& p_m)
{
  P_acc.slice(b) += (N.row(i).t() / T_m) * p_m.t();
}

// H: information theory
void accumulate_H_multigroup(arma::vec& H_acc, int b, double T_i,
                          const arma::vec& p_m, double logM)
{
  H_acc[b] += T_i * entropy(p_m, logM);
}

// R: relative diversity
void accumulate_R_multigroup(arma::vec& R_acc, int b, double T_i,
                          const arma::vec& p_m)
{
  R_acc[b] += T_i * interaction(p_m);
}

// D: dissimilarity
void accumulate_D_multigroup(arma::vec& D_acc, int b, double T_i,
                          const arma::vec& p_m,
                          const arma::vec& P_m,
                          double denominator)
{
  D_acc[b] += (T_i / denominator) * arma::accu(arma::abs(p_m - P_m));
}

void radius_environment(arma::vec& L_m, const std::vector<Neighbor>& neighbors,
                        const arma::mat& N, double band, double power,
                        int weighting_id, bool normalized)
{
  const Contributions c = radius_contributions(neighbors, band, power,
                                              weighting_id, normalized);
  L_m.zeros();
  for (std::size_t k = 0; k < c.ids.size(); ++k)
    L_m += c.weights[k] * N.row(c.ids[k]).t();
  if (c.scale != 1.0) L_m /= (1.0 / c.scale);
}

void knn_environment(arma::vec& L_m, const std::vector<Neighbor>& neighbors,
                     const arma::mat& N, const arma::vec& T_i,
                     double threshold, double power, int weighting_id,
                     bool normalized)
{
  const Contributions c = knn_contributions(neighbors, T_i, threshold, power,
                                           weighting_id, normalized);
  L_m.zeros();
  for (std::size_t k = 0; k < c.ids.size(); ++k)
    L_m += c.weights[k] * N.row(c.ids[k]).t();
}

void accumulate_pairwise(arma::mat& D_pair_acc,
                         arma::mat& R_pair_acc,
                         arma::mat& H_pair_acc,
                         int b, int q, const GroupPair& pair,
                         const arma::mat& N,
                         const arma::vec& L_m, int i,
                         const IntegerVector& measure_flags)
{
  if (!pair.valid)
    return;

  const double T_iab = N(i, pair.a) + N(i, pair.b);
  const double L_ab = L_m[pair.a] + L_m[pair.b];
  if (T_iab <= 0.0 || L_ab <= 0.0)
    return;

  const double p_a = L_m[pair.a] / L_ab;
  const double p_b = L_m[pair.b] / L_ab;
  const double log2 = std::log(2.0);

  if (measure_flags[3] == 1) {
    const double D_den = 2.0 * pair.T_ab * pair.I_ab;
    D_pair_acc(b, q) += (T_iab / D_den) *
      (std::abs(p_a - pair.P_a) + std::abs(p_b - pair.P_b));
  }

  if (measure_flags[2] == 1) {
    const double R_i = p_a * (1.0 - p_a) + p_b * (1.0 - p_b);
    R_pair_acc(b, q) += T_iab * R_i;
  }

  if (measure_flags[1] == 1) {
    double H_i = 0.0;
    if (p_a > 0.0)
      H_i -= p_a * std::log(p_a) / log2;
    if (p_b > 0.0)
      H_i -= p_b * std::log(p_b) / log2;
    H_pair_acc(b, q) += T_iab * H_i;
  }
}

NumericVector finalize_scalar_index(const arma::vec& acc,
                                    double denominator)
{
  NumericVector out(acc.n_elem);
  for (arma::uword i = 0; i < acc.n_elem; ++i)
    out[i] = 1.0 - acc[i] / denominator;
  return out;
}

List finalize_P_multigroup(const arma::cube& P_acc, int B)
{
  List out(B);
  for (int b = 0; b < B; ++b)
    out[b] = wrap(P_acc.slice(b));
  return out;
}

List finalize_pairwise_matrix(const arma::mat& acc,
                              const std::vector<GroupPair>& pairs, int B,
                              int M, int measure_id)
{
  List out(B);
  for (int b = 0; b < B; ++b) {
    NumericMatrix mat(M, M);
    std::fill(mat.begin(), mat.end(), NA_REAL);
    for (arma::uword q = 0; q < pairs.size(); ++q) {
      const GroupPair& pair = pairs[q];
      double value = NA_REAL;
      if (pair.valid) {
        const double raw = acc(b, q);
        if (measure_id == 1)
          value = 1.0 - raw / (pair.T_ab * pair.E_ab);
        else if (measure_id == 2)
          value = 1.0 - raw / (pair.T_ab * pair.I_ab);
        else if (measure_id == 3)
          value = raw;
      }
      mat(pair.a, pair.b) = value;
      mat(pair.b, pair.a) = value;
    }
    out[b] = mat;
  }
  return out;
}

List indices_from_localenvs(const arma::mat& N,
                            const std::vector<arma::mat>& environments,
                            const IntegerVector& measure_flags,
                            const IntegerVector& scope_flags)
{
  const int n = N.n_rows;
  const int M = N.n_cols;
  const int B = environments.size();
  const bool calculate_multigroup = scope_flags[0] == 1;
  const bool calculate_pairwise = scope_flags[1] == 1;

  arma::vec L_m(M);
  arma::vec p_m(M);
  arma::vec T_i = arma::sum(N, 1);
  arma::vec T_m = arma::trans(arma::sum(N, 0));
  const double T = arma::accu(T_i);
  arma::vec P_m = T_m / T;
  arma::vec D_acc(B, arma::fill::zeros);
  arma::vec R_acc(B, arma::fill::zeros);
  arma::vec H_acc(B, arma::fill::zeros);
  arma::cube P_acc(M, M, B, arma::fill::zeros);

  const double logM = std::log(static_cast<double>(M));
  const double E = entropy(P_m, logM);
  const double I = interaction(P_m);

  const std::vector<GroupPair> pairs = make_pairs(T_m);
  const int Q = pairs.size();
  arma::mat D_pair_acc(B, Q, arma::fill::zeros);
  arma::mat R_pair_acc(B, Q, arma::fill::zeros);
  arma::mat H_pair_acc(B, Q, arma::fill::zeros);

  for (int b = 0; b < B; ++b) {
    const arma::mat& L = environments[b];
    for (int i = 0; i < n; ++i) {
      L_m = L.row(i).t();
      normalize_composition(L_m, p_m);

      if (calculate_multigroup) {
        if (measure_flags[0] == 1)
          accumulate_P_multigroup(P_acc, b, i, N, T_m, p_m);

        if (measure_flags[1] == 1)
          accumulate_H_multigroup(H_acc, b, T_i[i], p_m, logM);

        if (measure_flags[2] == 1)
          accumulate_R_multigroup(R_acc, b, T_i[i], p_m);

        if (measure_flags[3] == 1)
          accumulate_D_multigroup(D_acc, b, T_i[i], p_m, P_m, 2.0 * T * I);
      }

      if (calculate_pairwise) {
        for (int q = 0; q < Q; ++q) {
          accumulate_pairwise(D_pair_acc, R_pair_acc, H_pair_acc, b, q,
                              pairs[q], N, L_m, i, measure_flags);
        }
      }
    }
  }

  NumericVector D_out;
  NumericVector R_out;
  NumericVector H_out;
  List P_out;

  if (calculate_multigroup && measure_flags[3] == 1) {
    D_out = NumericVector(B);
    for (int b = 0; b < B; ++b)
      D_out[b] = D_acc[b];
  }

  if (calculate_multigroup && measure_flags[2] == 1)
    R_out = finalize_scalar_index(R_acc, T * I);

  if (calculate_multigroup && measure_flags[1] == 1)
    H_out = finalize_scalar_index(H_acc, T * E);

  if (calculate_multigroup && measure_flags[0] == 1)
    P_out = finalize_P_multigroup(P_acc, B);

  List multigroup = List::create(
    _["d"] = D_out,
    _["r"] = R_out,
    _["h"] = H_out,
    _["p"] = P_out
  );

  List pairwise = List::create(
    _["d"] = (calculate_pairwise && measure_flags[3] == 1) ?
      finalize_pairwise_matrix(D_pair_acc, pairs, B, M, 3) : List(),
    _["r"] = (calculate_pairwise && measure_flags[2] == 1) ?
      finalize_pairwise_matrix(R_pair_acc, pairs, B, M, 2) : List(),
    _["h"] = (calculate_pairwise && measure_flags[1] == 1) ?
      finalize_pairwise_matrix(H_pair_acc, pairs, B, M, 1) : List()
  );

  return List::create(
    _["multigroup"] = multigroup,
    _["pairwise"] = pairwise
  );
}

} // namespace

// [[Rcpp::export]]
List seg_engine_cpp(NumericVector x, NumericVector y, NumericMatrix data,
                    NumericVector bands, double power, int weighting,
                    int normalize, IntegerVector measures,
                    IntegerVector scope, int keep_env, int keep_indices,
                    int neighbors, int search)
{
  // Input coordinates and Armadillo views of them
  const arma::vec X(x.begin(), x.size(), false);
  const arma::vec Y(y.begin(), y.size(), false);
  // Input counts (n by M) and Armadillo views
  arma::mat N(data.begin(), data.nrow(), data.ncol(), false);
  NumericVector BW(bands);             // bandwidth/search-radius values
  IntegerVector measure_flags(measures); // exposure, information, diversity, dissimilarity
  IntegerVector scope_flags(scope); // multigroup and pairwise switches

  const int n = X.size();              // number of focal units
  const int M = N.n_cols;              // number of groups
  const int B = BW.size();             // number of bandwidths
  const double weight_power = power;   // kernel power parameter
  const int weighting_id = weighting;  // kernel identifier from R
  const int neighbors_id = neighbors;  // 0 = radius, 1 = count-based kNN
  const int search_id = search;        // 0 = kd-tree, 1 = brute force
  const bool normalized = normalize == 1; // scale distances by bandwidth
  const bool calculate_multigroup = scope_flags[0] == 1; // compute multigroup indices
  const bool calculate_pairwise = scope_flags[1] == 1; // compute pairwise indices
  const bool store_env = keep_env == 1; // return local environments
  const bool store_indices = keep_indices == 1; // return segregation indices

  arma::vec L_m(M);                    // local environment counts by group
  arma::vec p_m(M);                    // local environment proportions by group
  arma::vec T_i = arma::sum(N, 1);     // total observed population by unit
  arma::vec T_m = arma::trans(arma::sum(N, 0)); // total observed population by group
  const double T = arma::accu(T_i);    // total observed population
  arma::vec P_m = T_m / T;             // regional group proportions
  arma::vec D_acc(B, arma::fill::zeros); // dissimilarity accumulator by band
  arma::vec R_acc(B, arma::fill::zeros); // relative diversity accumulator by band
  arma::vec H_acc(B, arma::fill::zeros); // information theory accumulator by band
  arma::cube P_acc(M, M, B, arma::fill::zeros); // exposure/isolation accumulator by band

  double max_bw = BW[0];
  for (int b = 1; b < B; ++b)
    max_bw = std::max(max_bw, BW[b]);

  PointCloud cloud{X, Y};
  std::unique_ptr<KDTree> tree;
  if (search_id == 0)
    tree.reset(new KDTree(
      2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)
    ));

  const double logM = std::log(static_cast<double>(M));
  const double E = entropy(P_m, logM);
  const double I = interaction(P_m);

  const std::vector<GroupPair> pairs = make_pairs(T_m);
  const int Q = pairs.size();
  arma::mat D_pair_acc(B, Q, arma::fill::zeros);
  arma::mat R_pair_acc(B, Q, arma::fill::zeros);
  arma::mat H_pair_acc(B, Q, arma::fill::zeros);

  List env_list;
  if (store_env) {
    env_list = List(B);
    for (int b = 0; b < B; ++b)
      env_list[b] = NumericMatrix(n, M);
  }

  for (int i = 0; i < n; ++i) {
    const std::vector<Neighbor> neighbors_i = discover_neighbors(
      tree.get(), i, X, Y, T_i, max_bw, neighbors_id == 1);

    for (int b = 0; b < B; ++b) {
      if (neighbors_id == 1) {
        knn_environment(L_m, neighbors_i, N, T_i, BW[b], weight_power,
                        weighting_id, normalized);
      } else {
        radius_environment(L_m, neighbors_i, N, BW[b], weight_power,
                           weighting_id, normalized);
      }

      if (store_env) {
        NumericMatrix env_matrix = env_list[b];
        for (int m = 0; m < M; ++m)
          env_matrix(i, m) = L_m[m];
      }

      if (!store_indices)
        continue;

      normalize_composition(L_m, p_m);

      if (calculate_multigroup) {
        if (measure_flags[0] == 1)
          accumulate_P_multigroup(P_acc, b, i, N, T_m, p_m);

        if (measure_flags[1] == 1)
          accumulate_H_multigroup(H_acc, b, T_i[i], p_m, logM);

        if (measure_flags[2] == 1)
          accumulate_R_multigroup(R_acc, b, T_i[i], p_m);

        if (measure_flags[3] == 1)
          accumulate_D_multigroup(D_acc, b, T_i[i], p_m, P_m, 2.0 * T * I);
      }

      if (calculate_pairwise) {
        for (int q = 0; q < Q; ++q) {
          accumulate_pairwise(D_pair_acc, R_pair_acc, H_pair_acc, b, q,
                              pairs[q], N, L_m, i, measure_flags);
        }
      }
    }
  }

  NumericVector D_out;
  NumericVector R_out;
  NumericVector H_out;
  List P_out;

  if (store_indices && calculate_multigroup && measure_flags[3] == 1) {
    D_out = NumericVector(B);
    for (int b = 0; b < B; ++b)
      D_out[b] = D_acc[b];
  }

  if (store_indices && calculate_multigroup && measure_flags[2] == 1)
    R_out = finalize_scalar_index(R_acc, T * I);

  if (store_indices && calculate_multigroup && measure_flags[1] == 1)
    H_out = finalize_scalar_index(H_acc, T * E);

  if (store_indices && calculate_multigroup && measure_flags[0] == 1)
    P_out = finalize_P_multigroup(P_acc, B);

  List multigroup = List::create(
    _["d"] = D_out,
    _["r"] = R_out,
    _["h"] = H_out,
    _["p"] = P_out
  );

  List pairwise = List::create(
    _["d"] = (store_indices && calculate_pairwise && measure_flags[3] == 1) ?
      finalize_pairwise_matrix(D_pair_acc, pairs, B, M, 3) : List(),
    _["r"] = (store_indices && calculate_pairwise && measure_flags[2] == 1) ?
      finalize_pairwise_matrix(R_pair_acc, pairs, B, M, 2) : List(),
    _["h"] = (store_indices && calculate_pairwise && measure_flags[1] == 1) ?
      finalize_pairwise_matrix(H_pair_acc, pairs, B, M, 1) : List()
  );

  List indices = List::create(
    _["multigroup"] = multigroup,
    _["pairwise"] = pairwise
  );

  SEXP env_out = store_env ? static_cast<SEXP>(env_list) : R_NilValue;

  return List::create(
    _["env"] = env_out,
    _["indices"] = indices
  );
}

// [[Rcpp::export]]
List seg_indices_env_cpp(NumericMatrix data, NumericMatrix env,
                         IntegerVector measures, IntegerVector scope)
{
  arma::mat N(data.begin(), data.nrow(), data.ncol(), false);
  arma::mat L(env.begin(), env.nrow(), env.ncol(), false);
  std::vector<arma::mat> environments;
  environments.push_back(L);

  return indices_from_localenvs(N, environments, measures, scope);
}
