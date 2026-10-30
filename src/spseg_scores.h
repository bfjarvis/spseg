#ifndef SPSEG_SCORES_H
#define SPSEG_SCORES_H

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include "spseg_types.h"

namespace spseg_scores {
using namespace spseg_spatial;
using Scores = std::array<double, MeasureCount>;
using Selection = std::array<bool, MeasureCount>;
constexpr double pi = 3.1415926535897932384626433832795;

inline double entropy_term(double p) {
  return p > 0 ? -p * std::log(p) : 0.0;
}

// E, I, G, Delta evaluated together. Logarithms serve both E and G.
// K is the number of positive regional shares, prepared once per analysis.
template<class Local, class Regional>
inline Scores local_scores(const Local& p, const Regional& P,
                           const Selection& selected, int K) {
  Scores a{};
  double log_product = 0;
  double product = 1;
  bool zero_product = false;
  for (std::size_t g = 0; g < p.size(); ++g) {
    const double pg = p[g];
    const bool in_root = selected[Root] && P[g] > 0;
    const double log_p = pg > 0 && (selected[Information] || (in_root && K != 2)) ?
      std::log(pg) : 0.0;
    if (selected[Information]) a[Information] -= pg * log_p;
    if (selected[Diversity]) a[Diversity] += pg * (1 - pg);
    if (selected[Dissimilarity]) a[Dissimilarity] += std::abs(pg - P[g]);
    if (in_root) {
      zero_product = zero_product || pg <= 0;
      if(K == 2) product *= pg;
      else log_product += log_p;
    }
  }
  if (selected[Root] && K > 0 && !zero_product)
    // The two-factor geometric mean has an exact square-root evaluation.
    a[Root] = K == 2 ? std::sqrt(product) : std::exp(log_product / K);
  return a;
}

template<class Vector>
inline double E(const Vector& p) {
  double out = 0;
  for (double pg : p) out += entropy_term(pg);
  return out;
}

template<class Vector>
inline double I(const Vector& p) {
  double out = 0;
  for (double pg : p) out += pg * (1 - pg);
  return out;
}

template<class Local, class Regional>
inline double G(const Local& p, const Regional& P, int K) {
  Selection selected{};
  selected[Root] = true;
  return local_scores(p, P, selected, K)[Root];
}

template<class Local, class Regional>
inline double Delta(const Local& p, const Regional& P) {
  double out = 0;
  for (std::size_t g = 0; g < p.size(); ++g) out += std::abs(p[g] - P[g]);
  return out;
}

template<class Vector>
inline Scores index_normalizers(const Vector& P, int K) {
  Scores out{};
  out[Information] = E(P);
  out[Diversity] = I(P);
  out[Root] = K > 1 ? G(P, P, K) : 0;
  out[Dissimilarity] = 2 * out[Diversity];
  return out;
}

// W = diag(1, 1/2, 2, 1/2) preserves the documented ordered-score units.
inline Scores cdf_scores(double C, double Q, const Selection& selected) {
  C = std::max(0.0, std::min(1.0, C));
  const std::array<double, 2> p{C, 1 - C}, P{Q, 1 - Q};
  Scores a = local_scores(p, P, selected, 2);
  a[Diversity] *= 0.5;
  a[Root] = C > 0 && C < 1 ? 2 * a[Root] : 0;
  a[Dissimilarity] *= 0.5;
  return a;
}

inline double cdf_score(double C, int family) {
  Selection selected{};
  selected[family] = true;
  return cdf_scores(C, 0.5, selected)[family];
}

// Integral from 0 to c of E((u, 1-u)); symmetry stabilizes the endpoints.
inline double entropy_antiderivative(double p) {
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return 0.5;
  if (p > 0.5) return 0.5 - entropy_antiderivative(1.0-p);
  if (p < 0.1) {
    double out = p*p * (0.75 - 0.5*std::log(p));
    double term = p*p*p;
    for (int n=2; n<30; ++n) {
      out -= term / (double(n)*(n-1)*(n+1));
      term *= p;
    }
    return out;
  }
  return (p - p*p*std::log(p) + (1-p)*(1-p)*std::log1p(-p))/2;
}

// Mean ordered score along a linear CDF segment, before rank-width scaling.
inline double mean_linear_score(double a, double b, int family) {
  a = std::max(0.0, std::min(1.0, a));
  b = std::max(a, std::min(1.0, b));
  if (a == b) return cdf_score(a, family);
  if (family == Diversity) return (a+b)/2 - (a*a+a*b+b*b)/3;
  if (a+b > 1.0) { const double old=a; a=1-b; b=1-old; }
  const double width = b-a;
  if (width < 1e-6 * std::max(a, 1e-300)) {
    const double mid=(a+b)/2, d=width/2*std::sqrt(3.0/5.0);
    return (5*cdf_score(mid-d,family)+8*cdf_score(mid,family)+
            5*cdf_score(mid+d,family))/18;
  }
  if (family == Information)
    return (entropy_antiderivative(b)-entropy_antiderivative(a))/width;
  const double ta=std::atan2(std::sqrt(a),std::sqrt(1-a));
  const double tb=std::atan2(std::sqrt(b),std::sqrt(1-b));
  const double dt=std::atan2(std::sqrt(b)*std::sqrt(1-a)-
                           std::sqrt(a)*std::sqrt(1-b),
                           std::sqrt(a*b)+std::sqrt((1-a)*(1-b)));
  const double x=2*dt;
  const double one_minus_sinc = std::abs(x)<0.01 ?
    x*x/6 - std::pow(x,4)/120 + std::pow(x,6)/5040 : 1-std::sin(x)/x;
  const double angle=ta+tb;
  return dt/2 * (2*std::pow(std::sin(angle),2) +
                std::cos(2*angle)*one_minus_sinc)/width;
}

inline double mean_linear_absolute(double a, double b) {
  if (a*b >= 0) return (std::abs(a)+std::abs(b))/2;
  return (a*a+b*b)/(2*(std::abs(a)+std::abs(b)));
}

inline Scores integrate_linear_scores(double C_left, double C_right,
    double Q_left, double Q_right, const Selection& selected) {
  Scores a{};
  const double width = Q_right - Q_left;
  if (width <= 0) return a;
  for (int f : {Information, Diversity, Root}) if (selected[f])
    a[f] = width * mean_linear_score(C_left, C_right, f);
  if (selected[Dissimilarity]) a[Dissimilarity] = width *
    mean_linear_absolute(C_left-Q_left, C_right-Q_right);
  return a;
}

inline void add_scores(Scores& out, const Scores& a) {
  for (int f : {Information, Diversity, Root, Dissimilarity}) out[f] += a[f];
}
} // namespace spseg_scores
#endif
