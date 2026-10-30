#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH)]]
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <numeric>
#include "spseg_neighbors.h"
#include "spseg_euclidean.h"
#include "spseg_graph.h"
#include "spseg_scores.h"

using namespace Rcpp;
using namespace spseg_spatial;
using namespace spseg_scores;

namespace {
struct VariableAnalysis {
  // N: focal counts; T_i: focal masses; P: regional category shares.
  // Q: regional CDF; p: current local shares; A: sum_i T_i * a_i by band.
  int type, M, n_reporting_groups, n_bands;
  NumericMatrix N_input;
  NumericMatrix group_fractions;
  arma::mat N;
  IntegerVector code, exposure_map, flags, env_map;
  IntegerVector split_intervals;
  NumericVector T_i, P, z, env_probs;
  int env_type, env_columns;
  bool continuous, env_mean;
  CharacterVector labels;
  std::vector<double> Q, F_prefix;
  double T;
  Scores normalizers{};
  Selection selected{};
  int K=0;
  arma::mat A;
  arma::cube exposure_acc;
  arma::cube rank_acc;
  std::vector<double> P_reporting;
  std::vector<std::pair<int,int>> pairs;
  arma::cube A_pair;
  List environments;
  List interpolation;
  bool keep_env, indices, multigroup, pairwise;
  std::vector<bool> empty;
  NumericVector knn_target, knn_min, knn_max, radius_min, radius_max;
  bool has_knn = false;
  int local_environment_count = 0;
  bool keep_scores, keep_contributions;
  List scores_overall, scores_pairwise, contributions_overall, contributions_pairwise;
  std::vector<Scores> pair_normalizers;

  VariableAnalysis(List config, int n, int b, bool store, bool calculate, bool multi,
                   bool pair, bool store_scores, bool store_contributions)
    : type(as<int>(config["type"])), n_bands(b),
      N_input(as<NumericMatrix>(config["counts"])),
      group_fractions(config["group_fractions"]),
      N(N_input.begin(), N_input.nrow(), N_input.ncol(), false),
      code(config["code"]), exposure_map(config["exposure_map"]),
      flags(config["flags"]), T_i(config["mass"]), P(config["shares"]),
      z(config["z"]), continuous(type==ContinuousIntervals || type==ContinuousExact),
      labels(config["exposure_labels"]),
      A(b,MeasureCount,arma::fill::zeros), environments(store ? b : 0),
      interpolation(store && group_fractions.nrow()>0 ? b : 0),
      keep_env(store), indices(calculate), multigroup(multi || type!=Categorical),
      pairwise(pair && (type==Categorical || continuous)), empty(b,false),
      knn_target(b), knn_min(b,R_PosInf), knn_max(b,0.0),
      radius_min(b,R_PosInf), radius_max(b,0.0),
      keep_scores(store_scores), keep_contributions(store_contributions)
  {
    env_map = as<IntegerVector>(config["env_map"]);
    env_probs = as<NumericVector>(config["env_probs"]);
    env_type = as<int>(config["env_type"]);
    env_columns = as<int>(config["env_columns"]);
    env_mean = as<bool>(config["env_mean"]);
    split_intervals = as<IntegerVector>(config["split_intervals"]);
    M=P.size(); n_reporting_groups=labels.size();
    for (int f = 0; f < MeasureCount; ++f) selected[f] = flags[f] != 0;
    P_reporting.assign(n_reporting_groups, 0.0);
    for (int g=0; g<M; ++g) {
      add_group_mass(P_reporting,g,P[g]);
    }
    T=std::accumulate(T_i.begin(),T_i.end(),0.0);
    double cumulative=0;
    for(int g=0;g<M;++g) {
      cumulative+=P[g]; Q.push_back(std::min(1.0,cumulative));
      if(P[g]>0) ++K;
    }
    if(M) Q.back()=1.0;
    F_prefix.assign(M,0.0);
    for(int g=1;g<M;++g) F_prefix[g]=F_prefix[g-1]+(z[g]-z[g-1])*Q[g-1];
    initialize_normalizers();
    if(indices && flags[Exposure]) {
      if(continuous && env_type==Quantiles) rank_acc.zeros(n_reporting_groups,env_columns,b);
      else exposure_acc.zeros(n_reporting_groups,n_reporting_groups,b);
    }
    if((indices || keep_scores || keep_contributions) && pairwise) {
      const int pair_groups=continuous ? n_reporting_groups : M;
      for(int a=0;a<pair_groups-1;++a) for(int c=a+1;c<pair_groups;++c) {
        pairs.push_back({a,c});
        const double pa = continuous ? P_reporting[a] : P[a];
        const double pc = continuous ? P_reporting[c] : P[c];
        const double sum = pa + pc;
        const std::array<double, 2> P_ab{sum > 0 ? pa/sum : 0,
                                        sum > 0 ? pc/sum : 0};
        pair_normalizers.push_back(index_normalizers(P_ab, pa > 0 && pc > 0 ? 2 : 1));
      }
      if(indices) A_pair.zeros(pairs.size(),MeasureCount,b);
    }
    if(keep_env) for(int band=0;band<b;++band) {
      environments[band]=NumericMatrix(n,env_columns);
      if(group_fractions.nrow()>0) interpolation[band]=NumericMatrix(n,3);
    }
    if(keep_scores) {
      if(multigroup) scores_overall = allocate_local_output(n, false);
      if(pairwise) scores_pairwise = allocate_local_output(n, true);
    }
    if(keep_contributions) {
      if(multigroup) contributions_overall = allocate_local_output(n, false);
      if(pairwise) contributions_pairwise = allocate_local_output(n, true);
    }
  }

  void initialize_normalizers() {
    if(type == Categorical) normalizers = index_normalizers(P, K);
    else if(type == ContinuousIntervals) {
      normalizers[Information] = 0.5;
      normalizers[Diversity] = 1.0/6;
      normalizers[Root] = pi/4;
    } else for(int g = 0; g < M-1; ++g) {
      const double width = type == Ordinal ? 1.0 : z[g+1]-z[g];
      for(int f : {Information, Diversity, Root})
        normalizers[f] += width * cdf_score(Q[g], f);
    }
    normalizers[Dissimilarity] = 2 * normalizers[Diversity];
  }

  List allocate_local_output(int n, bool pair_output) const {
    List out;
    const char* keys[MeasureCount] = {"p", "h", "r", "d", "s"};
    for(int f : {Information, Diversity, Root, Dissimilarity}) if(selected[f]) {
      List by_band(n_bands);
      for(int b = 0; b < n_bands; ++b) {
        if(pair_output) {
          NumericMatrix values(n, pairs.size());
          std::fill(values.begin(), values.end(), NA_REAL);
          by_band[b] = values;
        } else by_band[b] = NumericVector(n, NA_REAL);
      }
      out[keys[f]] = by_band;
    }
    return out;
  }

  void retain_scores(int i, int b, const Scores& a_i, double focal_population,
                     double population, const Scores& refs, bool constant,
                     List& score_output, List& contribution_output, int pair = -1) {
    const char* keys[MeasureCount] = {"p", "h", "r", "d", "s"};
    for(int f : {Information, Diversity, Root, Dissimilarity}) if(selected[f]) {
      if(keep_scores) {
        List by_band = score_output[keys[f]];
        if(pair < 0) { NumericVector values = by_band[b]; values[i] = a_i[f]; }
        else { NumericMatrix values = by_band[b]; values(i,pair) = a_i[f]; }
      }
      if(keep_contributions) {
        double value = NA_REAL;
        if(!constant && refs[f] > 0) {
          const double normalized = a_i[f] / refs[f];
          value = (focal_population / population) *
            (f == Dissimilarity ? normalized : 1-normalized);
        }
        List by_band = contribution_output[keys[f]];
        if(pair < 0) { NumericVector values = by_band[b]; values[i] = value; }
        else { NumericMatrix values = by_band[b]; values(i,pair) = value; }
      }
    }
  }

  void add_group_mass(std::vector<double>& out, int g, double value) const {
    if(group_fractions.nrow()>0) {
      for(int a=0;a<n_reporting_groups;++a) out[a]+=value*group_fractions(g,a);
    } else if(g<exposure_map.size() && exposure_map[g]>=0) {
      out[exposure_map[g]]+=value;
    }
  }

  double focal(int i,int g) const {
    return N.n_rows ? N(i,g) : (code[i]==g ? T_i[i] : 0);
  }

  double regional_cdf_antiderivative(double x) const {
    if(M==0 || x<=z[0]) return 0;
    auto it=std::upper_bound(z.begin(),z.end(),x);
    int g=int(it-z.begin())-1;
    return F_prefix[g]+(x-z[g])*Q[g];
  }

  double integrate_step_deviation(double a,double b,double p) const {
    if(a>=b) return 0;
    const auto it=std::upper_bound(Q.begin(),Q.end(),p);
    const double cross=it==Q.end() ? 1.0 : z[int(it-Q.begin())];
    const double split=std::max(a,std::min(b,cross));
    const double left=regional_cdf_antiderivative(split)-regional_cdf_antiderivative(a);
    const double right=regional_cdf_antiderivative(b)-regional_cdf_antiderivative(split);
    return std::max(0.0,p*(split-a)-left+right-p*(b-split));
  }

  // Retain only the requested presentation. Index calculations below continue
  // to use the complete p/sparse distribution for this focal row and band.
  std::vector<double> summarize_environment(const std::vector<double>& p,
                         const std::vector<std::pair<int,double>>& sparse,
                         double local) {
    std::vector<double> summary(env_columns,0.0);
    if(local <= 0) {
      std::fill(summary.begin(),summary.end(),NA_REAL);
      return summary;
    }
    if(env_type==Quantiles) {
      int target=0;
      double cumulative=0, last=NA_REAL;
      const double tolerance=8*std::numeric_limits<double>::epsilon();
      auto add = [&](int g, double share) {
        if(share <= 0) return;
        const double previous=cumulative;
        cumulative+=share;
        last=type==ContinuousIntervals ? Q[g] : z[g];
        while(target<env_probs.size() && env_probs[target]<1 &&
              env_probs[target]<=cumulative+tolerance) {
          if(type==ContinuousIntervals) {
            const double lower=g==0 ? 0 : Q[g-1];
            const double fraction=std::max(0.0,std::min(1.0,
              (env_probs[target]-previous)/share));
            summary[target]=lower+fraction*(Q[g]-lower);
          } else summary[target]=z[g];
          ++target;
        }
      };
      if(type==ContinuousExact) for(const auto& item:sparse) add(item.first,item.second);
      else for(int g=0;g<M;++g) add(g,p[g]);
      while(target<env_probs.size()) summary[target++]=last;
      if(env_mean) {
        double mean=0;
        if(type==ContinuousExact) for(const auto& item:sparse) mean+=item.second*z[item.first];
        else for(int g=0;g<M;++g) mean+=p[g]*z[g];
        summary.back()=mean;
      }
    } else {
      if(group_fractions.nrow()>0) {
        for(int g=0;g<M;++g) add_group_mass(summary,g,p[g]);
      } else if(type==ContinuousExact) {
        for(const auto& item:sparse) {
          const int g=env_map[item.first];
          if(g>=0) summary[g]+=item.second;
        }
      } else for(int g=0;g<M;++g)
        if(env_map[g]>=0) summary[env_map[g]]+=p[g];
    }
    return summary;
  }

  double build_local_distribution(const NeighborhoodWeights& c,
      std::vector<double>& p, std::vector<std::pair<int,double>>& sparse) const {
    double local=0;
    if(type==ContinuousExact && N.n_rows==0) {
      sparse.reserve(c.ids.size());
      for(std::size_t a=0;a<c.ids.size();++a) {
        int id=c.ids[a], g=code[id];
        double w=c.weights[a]*T_i[id];
        if(g>=0 && w>0) {sparse.push_back({g,w});local+=w;}
      }
      std::sort(sparse.begin(),sparse.end(),[](const std::pair<int,double>& a,
                 const std::pair<int,double>& b){return a.first<b.first;});
      std::size_t end=0;
      for(const auto& item:sparse) {
        if(end>0 && sparse[end-1].first==item.first) sparse[end-1].second+=item.second;
        else sparse[end++]=item;
      }
      sparse.resize(end);
      if(local>0) for(auto& item:sparse) item.second/=local;
    } else {
      p.assign(M,0);
      for(std::size_t a=0;a<c.ids.size();++a) {
        int id=c.ids[a]; double w=c.weights[a];
        if(N.n_rows) {
          for(int g=0;g<M;++g) p[g]+=w*N(id,g);
        } else if(code[id]>=0) p[code[id]]+=w*T_i[id];
      }
      local=std::accumulate(p.begin(),p.end(),0.0);
      if(local>0) for(double& v:p) v/=local;
      if(type==ContinuousExact) for(int g=0;g<M;++g) if(p[g]>0) sparse.push_back({g,p[g]});
    }
    return local;
  }

  Scores categorical_scores(const std::vector<double>& p) const {
    return local_scores(p, P, selected, K);
  }

  Scores ordinal_scores(const std::vector<double>& p) const {
    Scores a_i{};
    double C = 0;
    for(int g = 0; g < M-1; ++g) {
      C = std::min(1.0, C+p[g]);
      add_scores(a_i, cdf_scores(C, Q[g], selected));
    }
    return a_i;
  }

  Scores integrate_step_scores(double C, double z_left, double z_right) const {
    Scores a = cdf_scores(C, 0.5, selected);
    const double width = z_right-z_left;
    for(int f : {Information, Diversity, Root}) a[f] *= width;
    if(selected[Dissimilarity])
      a[Dissimilarity] = integrate_step_deviation(z_left, z_right, C);
    return a;
  }

  Scores continuous_scores(const std::vector<double>& p,
                          const std::vector<std::pair<int,double>>& p_sparse) const {
    Scores a_i{};
    double C = 0, previous = 0;
    if(type == ContinuousExact) {
      for(const auto& item : p_sparse) {
        const double next = z[item.first];
        add_scores(a_i, integrate_step_scores(C, previous, next));
        C = std::min(1.0, C+item.second);
        previous = next;
      }
      if(selected[Dissimilarity])
        a_i[Dissimilarity] += integrate_step_deviation(previous, 1.0, 1.0);
    } else {
      double Q_left = 0;
      for(int g = 0; g < M; ++g) {
        C = std::min(1.0, C+p[g]);
        add_scores(a_i, integrate_linear_scores(previous, C, Q_left, Q[g], selected));
        previous = C;
        Q_left = Q[g];
      }
    }
    return a_i;
  }

  void accumulate_scores(const Scores& a_i, double focal_population, int b,
                         int pair = -1) {
    if(pair < 0) {
      for(int f : {Information, Diversity, Root, Dissimilarity}) if(selected[f])
        A(b,f) += focal_population * a_i[f];
    } else {
      for(int f : {Information, Diversity, Root, Dissimilarity}) if(selected[f])
        A_pair(pair,f,b) += focal_population * a_i[f];
    }
  }

  void process_neighborhood(int i,int band,const NeighborhoodWeights& c) {
    if(c.knn && T_i[i]>0) {
      has_knn=true;
      knn_target[band]=c.target;
      knn_min[band]=std::min(knn_min[band],c.population);
      knn_max[band]=std::max(knn_max[band],c.population);
      radius_min[band]=std::min(radius_min[band],c.radius);
      radius_max[band]=std::max(radius_max[band],c.radius);
    }
    std::vector<double> p;
    std::vector<std::pair<int,double>> sparse;
    const double local = build_local_distribution(c, p, sparse);
    ++local_environment_count;
    std::vector<double> summary;
    if(keep_env && group_fractions.nrow()>0) {
      NumericMatrix diagnostic=interpolation[band];
      double apportioned=0;
      for(int g:split_intervals) apportioned+=p[g];
      diagnostic(i,0)=local>0 ? p.front() : NA_REAL;
      diagnostic(i,1)=local>0 ? p.back() : NA_REAL;
      diagnostic(i,2)=local>0 ? apportioned : NA_REAL;
    }
    const bool need_scores = indices || keep_scores || keep_contributions;
    if(keep_env || (need_scores && continuous &&
                    (flags[Exposure] || (pairwise && env_type==Shares)))) {
      summary=summarize_environment(p,sparse,local);
      if(keep_env) {
        NumericMatrix env=environments[band];
        for(int j=0;j<env_columns;++j) env(i,j)=summary[j];
      }
    }
    if(!need_scores || (T_i[i]<=0 && !keep_scores && !keep_contributions)) return;
    if(local<=0) {
      if(T_i[i]>0) empty[band]=true;
      else if(keep_contributions) {
        // A location with zero focal population has zero additive weight.
        const char* keys[MeasureCount] = {"p", "h", "r", "d", "s"};
        for(int f : {Information, Diversity, Root, Dissimilarity}) if(selected[f]) {
          if(multigroup && normalizers[f]>0 &&
              (type==ContinuousIntervals || K>1)) {
            List by_band = contributions_overall[keys[f]];
            NumericVector values = by_band[band];
            values[i] = 0;
          }
          if(pairwise) {
            List by_band = contributions_pairwise[keys[f]];
            NumericMatrix values = by_band[band];
            for(std::size_t pair=0; pair<pairs.size(); ++pair)
              if(pair_normalizers[pair][f]>0) values(i,pair) = 0;
          }
        }
      }
      return;
    }
    if(multigroup) {
      const Scores a_i = type == Categorical ? categorical_scores(p) :
        type == Ordinal ? ordinal_scores(p) : continuous_scores(p, sparse);
      if(indices && T_i[i]>0) accumulate_scores(a_i, T_i[i], band);
      if(keep_scores || keep_contributions)
        retain_scores(i, band, a_i, T_i[i], T, normalizers,
          type != ContinuousIntervals && K<2, scores_overall, contributions_overall);
    }
    std::vector<double> local_group, focal_group;
    if(flags[Exposure] || (pairwise && continuous)) {
      local_group.assign(n_reporting_groups,0); focal_group.assign(n_reporting_groups,0);
      // Observed continuous rows stay sparse; count input already has K columns.
      if(continuous && env_type==Shares) {
        local_group=summary;
      } else if(type==ContinuousExact) {
        for(const auto& item:sparse) {
          add_group_mass(local_group,item.first,item.second);
        }
      } else for(int g=0;g<M;++g) add_group_mass(local_group,g,p[g]);
      if(N.n_rows) {
        for(int g=0;g<M;++g) add_group_mass(focal_group,g,N(i,g));
      } else if(code[i]>=0) add_group_mass(focal_group,code[i],T_i[i]);
      // Regional group totals are prepared once; accumulate focal counts here.
      if(indices && flags[Exposure]) for(int a=0;a<n_reporting_groups;++a) if(focal_group[a]>0) {
        if(continuous && env_type==Quantiles) {
          for(int j=0;j<env_columns;++j) rank_acc(a,j,band)+=focal_group[a]*summary[j];
        } else {
          for(int b=0;b<n_reporting_groups;++b) exposure_acc(a,b,band)+=
            focal_group[a]*(continuous ? summary[b] : local_group[b]);
        }
      }
    }
    if(pairwise) for(std::size_t index=0;index<pairs.size();++index) {
      int a=pairs[index].first,b=pairs[index].second;
      double weight=continuous ? focal_group[a]+focal_group[b] : focal(i,a)+focal(i,b);
      double local_pair=continuous ? local_group[a]+local_group[b] : p[a]+p[b];
      const double pa = continuous ? P_reporting[a] : P[a];
      const double pb = continuous ? P_reporting[b] : P[b];
      if(pa<=0 || pb<=0) continue;
      if(weight<=0 && !keep_scores && !keep_contributions) continue;
      Scores a_i{};
      if(local_pair>0) {
        const double u_ab = (continuous ? local_group[a] : p[a])/local_pair;
        const std::array<double, 2> p_ab{u_ab, 1-u_ab};
        const std::array<double, 2> P_ab{pa/(pa+pb), pb/(pa+pb)};
        a_i = local_scores(p_ab, P_ab, selected, 2);
        if(indices && weight>0) accumulate_scores(a_i, weight, band, index);
      } else if(weight>0) continue;
      if(keep_scores || keep_contributions) {
        retain_scores(i, band, a_i, weight, T*(pa+pb), pair_normalizers[index],
          false, scores_pairwise, contributions_pairwise, index);
        if(local_pair<=0 && keep_scores) {
          const char* keys[MeasureCount] = {"p", "h", "r", "d", "s"};
          for(int f : {Information, Diversity, Root, Dissimilarity}) if(selected[f]) {
            List by_band = scores_pairwise[keys[f]];
            NumericMatrix values = by_band[band];
            values(i,index) = NA_REAL;
          }
        }
      }
    }
  }

  List finalize_results() {
    List multigroup_out, pairwise_out, rank_summary;
    const char* names[MeasureCount]={"p","h","r","d","s"};
    CharacterVector status(n_bands);
    const bool constant=type==ContinuousIntervals ? false : K<2;
    for(int b=0;b<n_bands;++b) status[b]=empty[b] ? "empty_environment" :
      (constant ? "constant" : "ok");
    if(indices && multigroup) for(int f:{Dissimilarity,Diversity,Information,Root}) if(flags[f]) {
      NumericVector out(n_bands,NA_REAL);
      for(int b=0;b<n_bands;++b) if(!empty[b] && !constant && normalizers[f]>0) {
        const double value=A(b,f)/(T*normalizers[f]);
        out[b]=f==Dissimilarity ? value : 1-value;
      }
      multigroup_out[names[f]]=out;
    }
    if(indices && flags[Exposure] && continuous && env_type==Quantiles) {
      rank_summary=List(n_bands);
      for(int b=0;b<n_bands;++b) {
        NumericMatrix out(n_reporting_groups,env_columns);
        for(int a=0;a<n_reporting_groups;++a) for(int j=0;j<env_columns;++j)
          out(a,j)=!empty[b] && P_reporting[a]>0 ?
            rank_acc(a,j,b)/(T*P_reporting[a]) : NA_REAL;
        rank_summary[b]=out;
      }
    } else if(indices && flags[Exposure]) {
      List matrices(n_bands);
      std::vector<double> totals(n_reporting_groups,0);
      for(int a=0;a<n_reporting_groups;++a) totals[a]=P_reporting[a]*T;
      for(int b=0;b<n_bands;++b) {
        NumericMatrix out(n_reporting_groups,n_reporting_groups);
        for(int a=0;a<n_reporting_groups;++a) for(int c=0;c<n_reporting_groups;++c)
          out(a,c)=!empty[b] && totals[a]>0 ? exposure_acc(a,c,b)/totals[a] : NA_REAL;
        out.attr("dimnames")=List::create(labels,labels); matrices[b]=out;
      }
      multigroup_out["p"]=matrices;
    }
    if(indices && pairwise) for(int f:{Dissimilarity,Diversity,Information,Root}) if(flags[f]) {
      List matrices(n_bands);
      for(int b=0;b<n_bands;++b) {
        const int pair_groups=continuous ? n_reporting_groups : M;
        NumericMatrix out(pair_groups,pair_groups); std::fill(out.begin(),out.end(),NA_REAL);
        for(std::size_t index=0;index<pairs.size();++index) {
          int a=pairs[index].first,c=pairs[index].second;
          const double pa=continuous ? P_reporting[a] : P[a];
          const double pc=continuous ? P_reporting[c] : P[c];
          if(pa<=0 || pc<=0 || empty[b]) continue;
          const double t=T*(pa+pc);
          const double denominator=pair_normalizers[index][f];
          const double v=A_pair(index,f,b)/(t*denominator);
          out(a,c)=out(c,a)=f==Dissimilarity ? v : 1-v;
        }
        matrices[b]=out;
      }
      pairwise_out[names[f]]=matrices;
    }
    List refs=List::create(_["information"]=normalizers[Information],_["diversity"]=normalizers[Diversity],
                          _["dissimilarity"]=normalizers[Dissimilarity],_["square_root"]=normalizers[Root]);
    List knn_summary;
    if(has_knn) knn_summary=List::create(_["requested"]=knn_target,
      _["min_population"]=knn_min,_["max_population"]=knn_max,
      _["min_radius"]=radius_min,_["max_radius"]=radius_max);
    List out = List::create(_["indices"]=List::create(_["multigroup"]=multigroup_out,
                        _["pairwise"]=pairwise_out),
                        _["env"]=keep_env ? static_cast<SEXP>(environments) : R_NilValue,
                        _["rank_summary"]=rank_summary,
                        _["interpolation"]=interpolation,
                        _["status"]=status,_["normalizers"]=refs,
                        _["knn"]=has_knn ? static_cast<SEXP>(knn_summary) : R_NilValue);
    if(keep_scores) out["scores"] = List::create(
      _["multigroup"]=scores_overall, _["pairwise"]=scores_pairwise);
    if(keep_contributions) {
      for(List store : {contributions_overall, contributions_pairwise}) {
        for(int f=0; f<store.size(); ++f) {
          List by_band = store[f];
          for(int b=0; b<n_bands; ++b) if(empty[b]) {
            NumericVector values = by_band[b];
            std::fill(values.begin(), values.end(), NA_REAL);
          }
        }
      }
      out["contributions"] = List::create(
        _["multigroup"]=contributions_overall, _["pairwise"]=contributions_pairwise);
    }
    return out;
  }
};
} // namespace

// [[Rcpp::export(rng = false)]]
List spseg_engine_cpp(NumericMatrix coords, List variables,
    NumericVector bands, double power, int weighting, bool normalize,
    int neighbors, int search, List network, IntegerVector scope, bool keep_env,
    bool keep_indices, bool keep_scores = false, bool keep_contributions = false)
{
  const int n=as<NumericVector>(as<List>(variables[0])["mass"]).size(), B=bands.size();
  arma::vec X(coords.begin(),coords.nrow(),false), Y(coords.begin()+coords.nrow(),coords.nrow(),false);
  const double largest=*std::max_element(bands.begin(),bands.end());
  const bool relational=network.size()>0;
  std::unique_ptr<EuclideanSearch> euclidean;
  std::unique_ptr<GraphSearch> graph;
  if(relational) graph.reset(new GraphSearch(network));
  else euclidean.reset(new EuclideanSearch(X,Y,search==Indexed));
  const CharacterVector names=variables.names();
  std::vector<std::unique_ptr<VariableAnalysis>> models;
  for(int v=0;v<variables.size();++v)
    models.emplace_back(new VariableAnalysis(as<List>(variables[v]),n,B,keep_env,
      keep_indices,scope[0],scope[1],keep_scores,keep_contributions));
  const bool knn = neighbors==PopulationNearest;
  // Identical eligible populations share both the stopping rule and kernels.
  // Views refer to analysis-owned R vectors; no n-by-variable copy is needed.
  std::vector<const double*> eligible;
  std::vector<std::vector<int>> members;
  if(knn) for(int v=0;v<int(models.size());++v) {
    const double* p=models[v]->T_i.begin();
    std::size_t group=0;
    while(group<eligible.size() && !std::equal(p,p+n,eligible[group])) ++group;
    if(group==eligible.size()) { eligible.push_back(p); members.emplace_back(); }
    members[group].push_back(v);
  }
  int focal_searches=0;
  for(int i=0;i<n;++i) {
    if(i%256==0) checkUserInterrupt();
    const auto candidates=relational ? graph->find(i,largest,knn,eligible) :
                                      euclidean->find(i,largest,knn,eligible);
    ++focal_searches;
    for(int b=0;b<B;++b) {
      if(knn) {
        for(std::size_t g=0;g<eligible.size();++g) {
          arma::vec T_i(const_cast<double*>(eligible[g]),n,false);
          const NeighborhoodWeights c=population_knn_weights(candidates,T_i,bands[b],power,
                                                  weighting,normalize);
          if(std::abs(c.population-bands[b]) >
             64*std::numeric_limits<double>::epsilon()*bands[b]) {
            const std::string name=as<std::string>(names[members[g][0]]);
            std::string location="focal row "+std::to_string(i+1);
            if(relational) {
              const IntegerVector node=network["node"];
              const CharacterVector labels=network["labels"];
              location+=" (node '"+as<std::string>(labels[node[i]])+"')";
            }
            stop("Variable '%s' cannot reach population %g at %s; reachable eligible population is %g",
                 name.c_str(),bands[b],location.c_str(),c.population);
          }
          for(int v:members[g]) models[v]->process_neighborhood(i,b,c);
        }
        continue;
      }
      const NeighborhoodWeights c=radius_weights(candidates,bands[b],power,weighting,normalize);
      for(auto& model:models) model->process_neighborhood(i,b,c);
    }
  }
  List out(variables.size());
  int local_environments = 0;
  for(int v=0;v<variables.size();++v) {
    out[v]=models[v]->finalize_results();
    local_environments += models[v]->local_environment_count;
  }
  out.attr("names")=variables.attr("names");
  return List::create(_["analyses"]=out,
    _["diagnostics"]=List::create(_["focal_searches"]=focal_searches,
                                   _["variables"]=variables.size(),
                                   _["local_environments"]=local_environments));
}
