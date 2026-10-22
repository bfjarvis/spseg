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

using namespace Rcpp;
using namespace spseg_spatial;

namespace {
constexpr double pi = 3.1415926535897932384626433832795;

double variation(double p, int family) {
  p = std::max(0.0, std::min(1.0, p));
  if (p == 0.0 || p == 1.0) return 0.0;
  if (family == 1) return -p * std::log(p) - (1.0-p) * std::log1p(-p);
  if (family == 2) return p * (1.0-p);
  return 2.0 * std::sqrt(p * (1.0-p));
}

double entropy_primitive(double p) {
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return 0.5;
  if (p > 0.5) return 0.5 - entropy_primitive(1.0-p);
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

// Mean of f along a linear CDF segment. Symmetry and narrow-interval limits
// avoid cancellation in entropy/arcsine divided differences.
double segment_mean(double a, double b, int family) {
  a = std::max(0.0, std::min(1.0, a));
  b = std::max(a, std::min(1.0, b));
  if (a == b) return variation(a, family);
  if (family == 2) return (a+b)/2 - (a*a+a*b+b*b)/3;
  if (a+b > 1.0) { const double old=a; a=1-b; b=1-old; }
  const double width = b-a;
  if (width < 1e-6 * std::max(a, 1e-300)) {
    // Three-point Gauss rule in the smooth, narrow-interval limit.
    const double mid=(a+b)/2, d=width/2*std::sqrt(3.0/5.0);
    return (5*variation(mid-d,family)+8*variation(mid,family)+
            5*variation(mid+d,family))/18;
  }
  if (family == 1) return (entropy_primitive(b)-entropy_primitive(a))/width;
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

double abs_segment(double a, double b) {
  if (a*b >= 0) return (std::abs(a)+std::abs(b))/2;
  return (a*a+b*b)/(2*(std::abs(a)+std::abs(b)));
}

struct Model {
  int type, k, groups, bands;
  NumericMatrix input;
  arma::mat counts;
  IntegerVector code, exposure_map, flags, env_map;
  NumericVector mass, shares, scores, z, env_probs, env_breaks;
  int env_type, env_columns;
  CharacterVector labels;
  std::vector<double> q, prefix;
  double total, reference[5] = {0,0,0,0,0};
  int active=0;
  std::vector<int> active_groups;
  arma::mat accum;
  arma::cube exposure_acc;
  std::vector<std::pair<int,int>> pairs;
  arma::cube pair_acc;
  List environments;
  bool keep_env, indices, multigroup, pairwise;
  std::vector<bool> empty;
  NumericVector knn_target, knn_min, knn_max, radius_min, radius_max;
  bool has_knn = false;

  Model(List config, int n, int b, bool store, bool calculate, bool multi, bool pair)
    : type(as<int>(config["type"])), bands(b),
      input(as<NumericMatrix>(config["counts"])),
      counts(input.begin(), input.nrow(), input.ncol(), false),
      code(config["code"]), exposure_map(config["exposure_map"]),
      flags(config["flags"]), mass(config["mass"]), shares(config["shares"]),
      scores(config["scores"]), z(config["z"]), labels(config["exposure_labels"]),
      accum(b,5,arma::fill::zeros), environments(store ? b : 0),
      keep_env(store), indices(calculate), multigroup(multi || type!=0),
      pairwise(pair && type==0), empty(b,false),
      knn_target(b), knn_min(b,R_PosInf), knn_max(b,0.0),
      radius_min(b,R_PosInf), radius_max(b,0.0)
  {
    env_map = as<IntegerVector>(config["env_map"]);
    env_probs = as<NumericVector>(config["env_probs"]);
    env_breaks = as<NumericVector>(config["env_breaks"]);
    env_type = as<int>(config["env_type"]);
    env_columns = as<int>(config["env_columns"]);
    k=shares.size(); groups=labels.size();
    total=std::accumulate(mass.begin(),mass.end(),0.0);
    double cumulative=0;
    for(int g=0;g<k;++g) {
      cumulative+=shares[g]; q.push_back(std::min(1.0,cumulative));
      if(shares[g]>0) {++active; active_groups.push_back(g);}
    }
    if(k) q.back()=1.0;
    prefix.assign(k,0.0);
    for(int g=1;g<k;++g) prefix[g]=prefix[g-1]+(z[g]-z[g-1])*q[g-1];
    if(type==0) {
      double logroot=0;
      for(int g:active_groups) {
        reference[1]-=shares[g]*std::log(shares[g]);
        reference[2]+=shares[g]*(1-shares[g]);
        logroot+=std::log(shares[g]);
      }
      reference[4]=active>1 ? std::exp(logroot/active) : 0;
    } else if(type==2) {
      reference[1]=0.5; reference[2]=1.0/6; reference[4]=pi/4;
      // A single populated interval still represents unobserved variation.
    } else {
      for(int g=0;g<k-1;++g) {
        const double width=type==1 ? 1.0 : z[g+1]-z[g];
        for(int f:{1,2,4}) reference[f]+=width*variation(q[g],f);
      }
    }
    reference[3]=2*reference[2];
    if(indices && flags[0]) exposure_acc.zeros(groups,groups,b);
    if(indices && pairwise) {
      for(int a=0;a<k-1;++a) for(int c=a+1;c<k;++c) pairs.push_back({a,c});
      pair_acc.zeros(pairs.size(),5,b);
    }
    if(keep_env) for(int band=0;band<b;++band)
      environments[band]=NumericMatrix(n,env_columns);
  }

  double focal(int i,int g) const {
    return counts.n_rows ? counts(i,g) : (code[i]==g ? mass[i] : 0);
  }

  double root(const std::vector<double>& p) const {
    double out=0;
    for(int g:active_groups) {
      if(p[g]<=0) return 0;
      out+=std::log(p[g]);
    }
    return active ? std::exp(out/active) : 0;
  }

  double integral_global(double x) const {
    if(k==0 || x<=z[0]) return 0;
    auto it=std::upper_bound(z.begin(),z.end(),x);
    int g=int(it-z.begin())-1;
    return prefix[g]+(x-z[g])*q[g];
  }

  double integral_abs(double a,double b,double p) const {
    if(a>=b) return 0;
    const auto it=std::upper_bound(q.begin(),q.end(),p);
    const double cross=it==q.end() ? 1.0 : z[int(it-q.begin())];
    const double split=std::max(a,std::min(b,cross));
    const double left=integral_global(split)-integral_global(a);
    const double right=integral_global(b)-integral_global(split);
    return std::max(0.0,p*(split-a)-left+right-p*(b-split));
  }

  // Retain only the requested presentation. Index calculations below continue
  // to use the complete p/sparse distribution for this focal row and band.
  void store_environment(int i, int band, const std::vector<double>& p,
                         const std::vector<std::pair<int,double>>& sparse,
                         double local) {
    NumericMatrix env = environments[band];
    if(local <= 0) {
      for(int j=0;j<env_columns;++j) env(i,j)=NA_REAL;
      return;
    }
    if(env_type == 1) {
      int target=0;
      double cumulative=0, last=NA_REAL;
      const double tolerance=8*std::numeric_limits<double>::epsilon();
      auto add = [&](int g, double share) {
        if(share <= 0) return;
        const double previous=cumulative;
        cumulative+=share;
        last=type==2 ? q[g] : z[g];
        while(target<env_columns && env_probs[target]<1 &&
              env_probs[target]<=cumulative+tolerance) {
          if(type==2) {
            const double lower=g==0 ? 0 : q[g-1];
            const double fraction=std::max(0.0,std::min(1.0,
              (env_probs[target]-previous)/share));
            env(i,target)=lower+fraction*(q[g]-lower);
          } else env(i,target)=z[g];
          ++target;
        }
      };
      if(type==3) for(const auto& item:sparse) add(item.first,item.second);
      else for(int g=0;g<k;++g) add(g,p[g]);
      while(target<env_columns) env(i,target++)=last;
    } else if(env_type == 2) {
      // Apportion grouped intervals to requested global-rank bins using the
      // same linear interpolant as the grouped segregation indices.
      for(int g=0;g<k;++g) {
        const double lower=g==0 ? 0 : q[g-1], width=q[g]-lower;
        if(width<=0 || p[g]<=0) continue;
        for(int j=0;j<env_columns;++j) {
          const double overlap=std::min(q[g],env_breaks[j+1])-
                               std::max(lower,env_breaks[j]);
          if(overlap>0) env(i,j)+=p[g]*overlap/width;
        }
      }
    } else {
      if(type==3) {
        for(const auto& item:sparse) {
          const int g=env_map[item.first];
          if(g>=0) env(i,g)+=item.second;
        }
      } else for(int g=0;g<k;++g)
        if(env_map[g]>=0) env(i,env_map[g])+=p[g];
    }
  }

  void update(int i,int band,const Contributions& c) {
    if(c.knn && mass[i]>0) {
      has_knn=true;
      knn_target[band]=c.target;
      knn_min[band]=std::min(knn_min[band],c.population);
      knn_max[band]=std::max(knn_max[band],c.population);
      radius_min[band]=std::min(radius_min[band],c.radius);
      radius_max[band]=std::max(radius_max[band],c.radius);
    }
    std::vector<double> p;
    std::vector<std::pair<int,double>> sparse;
    double local=0;
    if(type==3 && counts.n_rows==0) {
      sparse.reserve(c.ids.size());
      for(std::size_t a=0;a<c.ids.size();++a) {
        int id=c.ids[a], g=code[id];
        double w=c.weights[a]*mass[id];
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
      p.assign(k,0);
      for(std::size_t a=0;a<c.ids.size();++a) {
        int id=c.ids[a]; double w=c.weights[a];
        if(counts.n_rows) {
          for(int g=0;g<k;++g) p[g]+=w*counts(id,g);
        } else if(code[id]>=0) p[code[id]]+=w*mass[id];
      }
      local=std::accumulate(p.begin(),p.end(),0.0);
      if(local>0) for(double& v:p) v/=local;
      if(type==3) for(int g=0;g<k;++g) if(p[g]>0) sparse.push_back({g,p[g]});
    }
    if(keep_env) store_environment(i,band,p,sparse,local);
    if(!indices || mass[i]<=0) return;
    if(local<=0) {empty[band]=true;return;}
    if(multigroup) {
      if(type==0) {
        for(int g=0;g<k;++g) {
          if(flags[1] && p[g]>0) accum(band,1)-=mass[i]*p[g]*std::log(p[g]);
          if(flags[2]) accum(band,2)+=mass[i]*p[g]*(1-p[g]);
          if(flags[3]) accum(band,3)+=mass[i]*std::abs(p[g]-shares[g]);
        }
        if(flags[4]) accum(band,4)+=mass[i]*root(p);
      } else if(type==3) {
        double cumulative=0,previous=0;
        for(std::size_t a=0;a<sparse.size();++a) {
          double next=z[sparse[a].first];
          if(flags[3]) accum(band,3)+=mass[i]*integral_abs(previous,next,cumulative);
          for(int f:{1,2,4}) if(flags[f])
            accum(band,f)+=mass[i]*(next-previous)*variation(cumulative,f);
          cumulative=std::min(1.0,cumulative+sparse[a].second);
          previous=next;
        }
        if(flags[3]) accum(band,3)+=mass[i]*integral_abs(previous,1.0,1.0);
      } else {
        double cumulative=0,previous=0,global_previous=0;
        for(int g=0;g<k;++g) {
          cumulative=std::min(1.0,cumulative+p[g]);
          if(type==1 && g<k-1) {
            for(int f:{1,2,4}) if(flags[f]) accum(band,f)+=mass[i]*variation(cumulative,f);
            if(flags[3]) accum(band,3)+=mass[i]*std::abs(cumulative-q[g]);
          } else if(type==2) {
            const double width=q[g]-global_previous;
            for(int f:{1,2,4}) if(flags[f] && width>0)
              accum(band,f)+=mass[i]*width*segment_mean(previous,cumulative,f);
            if(flags[3]) accum(band,3)+=mass[i]*width*
              abs_segment(previous-global_previous,cumulative-q[g]);
            previous=cumulative; global_previous=q[g];
          }
        }
      }
    }
    if(flags[0]) {
      std::vector<double> local_group(groups,0), focal_group(groups,0);
      // Observed continuous rows stay sparse; count input already has K columns.
      if(type==3) {
        for(const auto& item:sparse) {
          int g=exposure_map[item.first]; if(g>=0) local_group[g]+=item.second;
        }
      } else for(int g=0;g<k;++g) local_group[exposure_map[g]]+=p[g];
      if(counts.n_rows) {
        for(int g=0;g<k;++g) if(exposure_map[g]>=0)
          focal_group[exposure_map[g]]+=counts(i,g);
      } else if(code[i]>=0 && exposure_map[code[i]]>=0)
        focal_group[exposure_map[code[i]]]=mass[i];
      // Group totals are precomputed in finish(); accumulate population counts here.
      for(int a=0;a<groups;++a) if(focal_group[a]>0)
        for(int b=0;b<groups;++b) exposure_acc(a,b,band)+=focal_group[a]*local_group[b];
    }
    if(pairwise) for(std::size_t index=0;index<pairs.size();++index) {
      int a=pairs[index].first,b=pairs[index].second;
      double weight=focal(i,a)+focal(i,b), local_pair=p[a]+p[b];
      if(weight<=0 || local_pair<=0) continue;
      double v=p[a]/local_pair;
      for(int f:{1,2,4}) if(flags[f]) pair_acc(index,f,band)+=weight*variation(v,f);
      double q0=shares[a]/(shares[a]+shares[b]);
      if(flags[3]) pair_acc(index,3,band)+=weight*std::abs(v-q0);
    }
  }

  List finish() {
    List multigroup, pairwise_out;
    const char* names[5]={"p","h","r","d","s"};
    CharacterVector status(bands);
    const bool constant=type==2 ? false : active<2;
    for(int b=0;b<bands;++b) status[b]=empty[b] ? "empty_environment" :
      (constant ? "constant" : "ok");
    if(indices && multigroup) for(int f:{3,2,1,4}) if(flags[f]) {
      NumericVector out(bands,NA_REAL);
      for(int b=0;b<bands;++b) if(!empty[b] && !constant && reference[f]>0) {
        const double value=accum(b,f)/(total*reference[f]);
        out[b]=f==3 ? value : 1-value;
      }
      multigroup[names[f]]=out;
    }
    if(indices && flags[0]) {
      List matrices(bands);
      std::vector<double> totals(groups,0);
      for(int g=0;g<k;++g) if(exposure_map[g]>=0)
        totals[exposure_map[g]]+=shares[g]*total;
      for(int b=0;b<bands;++b) {
        NumericMatrix out(groups,groups);
        for(int a=0;a<groups;++a) for(int c=0;c<groups;++c)
          out(a,c)=!empty[b] && totals[a]>0 ? exposure_acc(a,c,b)/totals[a] : NA_REAL;
        out.attr("dimnames")=List::create(labels,labels); matrices[b]=out;
      }
      multigroup["p"]=matrices;
    }
    if(indices && pairwise) for(int f:{3,2,1,4}) if(flags[f]) {
      List matrices(bands);
      for(int b=0;b<bands;++b) {
        NumericMatrix out(k,k); std::fill(out.begin(),out.end(),NA_REAL);
        for(std::size_t index=0;index<pairs.size();++index) {
          int a=pairs[index].first,c=pairs[index].second;
          if(shares[a]<=0 || shares[c]<=0 || empty[b]) continue;
          const double t=total*(shares[a]+shares[c]);
          const double q0=shares[a]/(shares[a]+shares[c]);
          const double denominator=f==3 ? 2*variation(q0,2) : variation(q0,f);
          const double v=pair_acc(index,f,b)/(t*denominator);
          out(a,c)=out(c,a)=f==3 ? v : 1-v;
        }
        matrices[b]=out;
      }
      pairwise_out[names[f]]=matrices;
    }
    List refs=List::create(_["information"]=reference[1],_["diversity"]=reference[2],
                          _["dissimilarity"]=reference[3],_["square_root"]=reference[4]);
    List knn_summary;
    if(has_knn) knn_summary=List::create(_["requested"]=knn_target,
      _["min_population"]=knn_min,_["max_population"]=knn_max,
      _["min_radius"]=radius_min,_["max_radius"]=radius_max);
    return List::create(_["indices"]=List::create(_["multigroup"]=multigroup,
                        _["pairwise"]=pairwise_out),
                        _["env"]=keep_env ? static_cast<SEXP>(environments) : R_NilValue,
                        _["status"]=status,_["normalizers"]=refs,
                        _["knn"]=has_knn ? static_cast<SEXP>(knn_summary) : R_NilValue);
  }
};
} // namespace

// [[Rcpp::export(rng = false)]]
List spseg_engine_cpp(NumericMatrix coords, List variables,
    NumericVector bands, double power, int weighting, bool normalize,
    int neighbors, int search, List network, IntegerVector scope, bool keep_env, bool keep_indices)
{
  const int n=as<NumericVector>(as<List>(variables[0])["mass"]).size(), B=bands.size();
  arma::vec X(coords.begin(),coords.nrow(),false), Y(coords.begin()+coords.nrow(),coords.nrow(),false);
  const double largest=*std::max_element(bands.begin(),bands.end());
  const bool relational=network.size()>0;
  std::unique_ptr<EuclideanSearch> euclidean;
  std::unique_ptr<GraphSearch> graph;
  if(relational) graph.reset(new GraphSearch(network));
  else euclidean.reset(new EuclideanSearch(X,Y,search==0));
  const CharacterVector names=variables.names();
  std::vector<std::unique_ptr<Model>> models;
  for(int v=0;v<variables.size();++v)
    models.emplace_back(new Model(as<List>(variables[v]),n,B,keep_env,keep_indices,scope[0],scope[1]));
  const bool knn = neighbors==1;
  // Identical eligible populations share both the stopping rule and kernels.
  // Views refer to Model-owned R vectors; no n-by-variable copy is needed.
  std::vector<const double*> eligible;
  std::vector<std::vector<int>> members;
  if(knn) for(int v=0;v<int(models.size());++v) {
    const double* p=models[v]->mass.begin();
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
          arma::vec mass(const_cast<double*>(eligible[g]),n,false);
          const Contributions c=knn_contributions(candidates,mass,bands[b],power,
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
          for(int v:members[g]) models[v]->update(i,b,c);
        }
        continue;
      }
      const Contributions c=radius_contributions(candidates,bands[b],power,weighting,normalize);
      for(auto& model:models) model->update(i,b,c);
    }
  }
  List out(variables.size());
  for(int v=0;v<variables.size();++v) out[v]=models[v]->finish();
  out.attr("names")=variables.attr("names");
  return List::create(_["analyses"]=out,
    _["diagnostics"]=List::create(_["focal_searches"]=focal_searches,
                                   _["variables"]=variables.size()));
}
