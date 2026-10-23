/*!
 * \brief Spatial-order diagnostics. Alternative assemblies are TEST ONLY.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../SU2_CFD/include/gradients/computeGradientsLeastSquares.hpp"

using namespace PrimitiveEuler2D;
using namespace PrimitiveEulerTest;

namespace {
constexpr unsigned groups = 6;
const std::array<const char*, groups> names{{"all", "boundary", "ring1", "deep3", "ring2", "corners"}};
struct Grid {
  unsigned n;
  std::vector<Direction> xy;
  std::vector<Triangle> triangles;
  Grid(unsigned count, bool distorted) : n(count) {
    const double h = 1.0/n;
    for (unsigned j=0; j<=n; ++j) for (unsigned i=0; i<=n; ++i) {
      Direction x{{i*h,j*h}};
      if (distorted && i>0 && i<n && j>0 && j<n) {
        x[0] += 0.18*h*std::sin(2.1*i+1.3*j);
        x[1] += 0.16*h*std::cos(1.7*i-0.9*j);
      }
      xy.push_back(x);
    }
    for (unsigned j=0; j<n; ++j) for (unsigned i=0; i<n; ++i) {
      const unsigned long a=j*(n+1)+i, b=a+1, d=a+n+1, c=d+1;
      triangles.push_back({{a,b,c}}); triangles.push_back({{a,c,d}});
    }
  }
  unsigned layer(unsigned long point) const {
    const unsigned i=point%(n+1), j=point/(n+1);
    return std::min({i,j,n-i,n-j});
  }
  bool in(unsigned long point, unsigned group) const {
    const unsigned l=layer(point), i=point%(n+1), j=point/(n+1);
    if (group==0) return true;
    if (group==1) return l==0;
    if (group==2) return l==1;
    if (group==3) return l>=3;
    if (group==4) return l==2;
    return (i==0 || i==n) && (j==0 || j==n);
  }
};

State W(const Direction& x) {
  return {{2+0.15*std::sin(x[0])*std::cos(x[1]),
           0.5+0.1*std::cos(1.3*x[0])*std::sin(0.7*x[1]),
           0.3+0.08*std::sin(0.8*x[0])*std::cos(1.1*x[1]),
           3+0.2*std::cos(0.9*x[0])*std::cos(1.2*x[1])}};
}
Gradient G(const Direction& x) {
  return {{{{0.15*std::cos(x[0])*std::cos(x[1]), -0.15*std::sin(x[0])*std::sin(x[1])}},
           {{-0.13*std::sin(1.3*x[0])*std::sin(0.7*x[1]), 0.07*std::cos(1.3*x[0])*std::cos(0.7*x[1])}},
           {{0.064*std::cos(0.8*x[0])*std::cos(1.1*x[1]), -0.088*std::sin(0.8*x[0])*std::sin(1.1*x[1])}},
           {{-0.18*std::sin(0.9*x[0])*std::cos(1.2*x[1]), -0.24*std::cos(0.9*x[0])*std::sin(1.2*x[1])}}}};
}
// TEST-ONLY analytic Hessians isolate physical-trace truncation; no fitted data.
using Hessian = std::array<std::array<double,3>,4>;
// Hessian rows: xx, xy, yy for each primitive component.
Hessian H(const Direction& x) {
  return {{{{-0.15*std::sin(x[0])*std::cos(x[1]), -0.15*std::cos(x[0])*std::sin(x[1]), -0.15*std::sin(x[0])*std::cos(x[1])}},
           {{-0.169*std::cos(1.3*x[0])*std::sin(0.7*x[1]), -0.091*std::sin(1.3*x[0])*std::cos(0.7*x[1]), -0.049*std::cos(1.3*x[0])*std::sin(0.7*x[1])}},
           {{-0.0512*std::sin(0.8*x[0])*std::cos(1.1*x[1]), -0.0704*std::cos(0.8*x[0])*std::sin(1.1*x[1]), -0.0968*std::sin(0.8*x[0])*std::cos(1.1*x[1])}},
           {{-0.162*std::cos(0.9*x[0])*std::cos(1.2*x[1]), 0.216*std::sin(0.9*x[0])*std::sin(1.2*x[1]), -0.288*std::cos(0.9*x[0])*std::cos(1.2*x[1])}}}};
}
State Product(const State& w, const Gradient& g) {
  return {{w[1]*g[0][0]+w[2]*g[0][1]+w[0]*(g[1][0]+g[2][1]),
           w[1]*g[1][0]+w[2]*g[1][1]+g[3][0]/w[0],
           w[1]*g[2][0]+w[2]*g[2][1]+g[3][1]/w[0],
           w[1]*g[3][0]+w[2]*g[3][1]+1.4*w[3]*(g[1][0]+g[2][1])}};
}
const CPrimitiveEulerSolver::BoundaryState exactBoundary =
    [](const Direction& x,double,const State&,const Direction&,unsigned short) { return W(x); };
const CPrimitiveEulerSolver::Source exactSource = [](const Direction& x,double) { return Product(W(x),G(x)); };
std::string Options(const std::string& method, const std::string& volume) {
  return "TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\nMUSCL_FLOW=YES\nSLOPE_LIMITER_FLOW=NONE\n"
         "NUM_METHOD_GRAD_RECON="+method+"\nPRIMITIVE_VOLUME_QUADRATURE="+volume+"\n";
}

struct Norms {
  std::array<double,groups> squared{}, weight{};
  void add(const Grid& grid, unsigned long i, double volume, double squaredError) {
    for (unsigned g=0; g<groups; ++g) if (grid.in(i,g)) {
      squared[g] += volume*squaredError; weight[g] += volume;
    }
  }
  double rms(unsigned g) const { return std::sqrt(squared[g]/weight[g]); }
};
double Squared(const State& state) { double sum=0; for (auto value:state) sum+=value*value; return sum; }
void Print(const char* kind, bool distorted, const std::string& method, const std::string& variant,
           unsigned n, const Norms& norms, const std::array<double,groups>& previous) {
  for (unsigned g=0; g<groups; ++g) {
    const double error=norms.rms(g);
    REQUIRE(std::isfinite(error));
    const double p=n==8 || previous[g]==0 || error==0 ? 0 : std::log2(previous[g]/error);
    std::cout << std::setprecision(12) << "ORDER_AUDIT," << kind << ',' << (distorted?"distorted":"regular") << ','
              << method << ',' << variant << ',' << n << ',' << names[g] << ',' << norms.weight[g] << ','
              << error << ',' << std::sqrt(norms.squared[g]/norms.weight[0]) << ',' << p << '\n';
  }
}

// Independent TEST-ONLY reassembly. Supplying analytic derivatives removes one
// error mechanism at a time; no option or alternative is added to the solver.
CPrimitiveEulerSolver::SpatialResult Assemble(const Grid& grid, const CAggregatedGeometry& metrics,
    const std::vector<Gradient>& gradient, const std::string& volume, bool exactTrace=false,
    const std::vector<Hessian>* hessian=nullptr) {
  CPrimitiveEulerSolver::SpatialResult result;
  for (auto* a : {&result.faces,&result.volume,&result.boundary,&result.source,&result.rhs}) a->resize(grid.xy.size());
  const auto reconstruct=[&](unsigned long i,const Direction& q,bool physical) {
    Direction d{{q[0]-grid.xy[i][0],q[1]-grid.xy[i][1]}};
    auto value=Reconstruct(W(grid.xy[i]),gradient[i],d);
    if (physical && exactTrace) return W(q);
    if (physical && hessian) for (unsigned k=0;k<4;++k) {
      const auto& h=(*hessian)[i][k];
      value[k]+=0.5*h[0]*d[0]*d[0]+h[1]*d[0]*d[1]+0.5*h[2]*d[1]*d[1];
    }
    return value;
  };
  for (unsigned long e=0;e<metrics.GetnEdge();++e) {
    const auto& f=metrics.GetFace(e);
    const auto waves=ComputeFluctuations(reconstruct(f.first,f.point,false),reconstruct(f.second,f.point,false),f.normal,1.4);
    for (unsigned k=0;k<4;++k) { result.faces[f.first][k]+=waves.minus[k]; result.faces[f.second][k]+=waves.plus[k]; }
  }
  for (const auto& f:metrics.GetBoundaryFaces()) {
    const auto wave=ComputeFluctuations(reconstruct(f.node,f.point,true),W(f.point),f.normal,1.4).minus;
    for (unsigned k=0;k<4;++k) result.boundary[f.node][k]+=wave[k];
  }
  for (unsigned long i=0;i<grid.xy.size();++i) {
    const auto& cell=metrics.GetCell(i);
    auto mean=W(grid.xy[i]);
    if (volume=="MEAN_COEFF") mean=Reconstruct(mean,gradient[i],cell.offset);
    auto product=Product(mean,gradient[i]);
    result.source[i]=Product(W(grid.xy[i]),G(grid.xy[i]));
    for (unsigned k=0;k<4;++k) {
      result.volume[i][k]=cell.volume*product[k];
      result.rhs[i][k]=-(result.faces[i][k]+result.boundary[i][k])/cell.volume-product[k]+result.source[i][k];
    }
  }
  return result;
}
}

TEST_CASE("Regional MMS separates boundary rings deep interior and spatial contributions", "[primitive-solver][primitive-order-audit]") {
  std::cout << "ORDER_AUDIT,kind,family,gradient,variant,N,group,volume,L2_group,L2_global_contribution,order\n";
  for (bool distorted:{false,true}) for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"}) {
    for (const std::string mode:{"NODE_POINT","MEAN_COEFF"}) {
      std::map<std::string,std::array<double,groups>> previous;
      for (unsigned n:{8u,16u,32u,64u,128u}) {
        Grid grid(n,distorted);
        MeshFixture fixture(grid.xy,grid.triangles,"PRIMITIVE_EULER",Options(method,mode));
        CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
        std::vector<Gradient> raw(grid.xy.size()), exact(grid.xy.size()), boundaryExact(grid.xy.size());
        std::vector<Hessian> hessian(grid.xy.size());
        for (unsigned long i=0;i<grid.xy.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,W(grid.xy[i]));
        const auto actual=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactBoundary,exactSource);
        const auto& metrics=solver.GetGeometryMetrics();
        for (unsigned long i=0;i<grid.xy.size();++i) {
          REQUIRE(solver.GetPositivityScale(i)==1);
          raw[i]=solver.GetLimitedGradient(i); exact[i]=G(grid.xy[i]); hessian[i]=H(grid.xy[i]);
          boundaryExact[i]=grid.layer(i)==0?exact[i]:raw[i];
        }
        const auto copy=Assemble(grid,metrics,raw,mode);
        double difference=0;
        for (unsigned long i=0;i<grid.xy.size();++i) for (unsigned k=0;k<4;++k)
          difference=std::max(difference,std::abs(copy.rhs[i][k]-actual.rhs[i][k]));
        REQUIRE(difference<2e-12); // Validate the diagnostic against the production assembly first.
        std::map<std::string,CPrimitiveEulerSolver::SpatialResult> experiments;
        experiments.emplace("native",actual);
        experiments.emplace("exact_boundary_G",Assemble(grid,metrics,boundaryExact,mode));
        experiments.emplace("exact_all_G",Assemble(grid,metrics,exact,mode));
        experiments.emplace("exact_G_trace",Assemble(grid,metrics,exact,mode,true));
        experiments.emplace("exact_G_quadratic_trace",Assemble(grid,metrics,exact,mode,false,&hessian));
        for (const auto& experiment:experiments) {
          Norms norms;
          for (unsigned long i=0;i<grid.xy.size();++i) norms.add(grid,i,metrics.GetCell(i).volume,Squared(experiment.second.rhs[i]));
          Print("residual",distorted,method,mode+":"+experiment.first,n,norms,previous[experiment.first]);
          for (unsigned g=0;g<groups;++g) previous[experiment.first][g]=norms.rms(g);
        }
        for (const std::string part:{"faces","boundary","volume_defect","gradient_error"}) {
          Norms norms;
          for (unsigned long i=0;i<grid.xy.size();++i) {
            const double v=metrics.GetCell(i).volume;
            State error{};
            double squared=0;
            if (part=="gradient_error") {
              for (unsigned k=0;k<4;++k) for (unsigned d=0;d<2;++d) squared+=std::pow(raw[i][k][d]-exact[i][k][d],2);
            } else {
              for (unsigned k=0;k<4;++k)
                error[k]=part=="faces"?-actual.faces[i][k]/v:
                    (part=="boundary"?-actual.boundary[i][k]/v:actual.source[i][k]-actual.volume[i][k]/v);
              squared=Squared(error);
            }
            norms.add(grid,i,v,squared);
          }
          Print("component",distorted,method,mode+":"+part,n,norms,previous[part]);
          for (unsigned g=0;g<groups;++g) previous[part][g]=norms.rms(g);
        }
      }
    }
  }
}

TEST_CASE("Native LS and WLS quadratic gradient accuracy is measured without residual assembly", "[primitive-solver][primitive-order-audit]") {
  const Gradient linear{{{{0.2,-0.1}},{{0.12,0.07}},{{-0.08,0.09}},{{0.3,0.11}}}};
  const Hessian quadratic{{{{0.3,0.17,-0.12}},{{-0.14,0.09,0.23}},{{0.18,-0.13,0.11}},{{-0.21,0.15,0.26}}}};
  for (bool distorted:{false,true}) for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"}) {
    std::array<double,groups> previous{};
    for (unsigned n:{8u,16u,32u,64u,128u}) {
      Grid grid(n,distorted);
      MeshFixture fixture(grid.xy,grid.triangles,"PRIMITIVE_EULER",Options(method,"NODE_POINT"));
      CPrimitiveEulerVariable field({{2,0.5,0.3,3}},grid.xy.size(),fixture.config.get());
      for (unsigned long i=0;i<grid.xy.size();++i) {
        const auto& x=grid.xy[i]; State w{{2,0.5,0.3,3}};
        for (unsigned k=0;k<4;++k) w[k]+=linear[k][0]*x[0]+linear[k][1]*x[1]+
            0.5*quadratic[k][0]*x[0]*x[0]+quadratic[k][1]*x[0]*x[1]+0.5*quadratic[k][2]*x[1]*x[1];
        field.SetPrimitiveState(i,w);
      }
      computeGradientsLeastSquares(nullptr,MPI_QUANTITIES::SOLUTION_GRAD_REC,PERIODIC_NONE,
          *fixture.geometry,*fixture.config,method=="WEIGHTED_LEAST_SQUARES",field.GetSolution(),0,4,1,
          field.GetGradient_Reconstruction(),field.GetRmatrix());
      Norms norms;
      for (unsigned long i=0;i<grid.xy.size();++i) {
        double squared=0;
        for (unsigned k=0;k<4;++k) {
          const Direction exact{{linear[k][0]+quadratic[k][0]*grid.xy[i][0]+quadratic[k][1]*grid.xy[i][1],
                                 linear[k][1]+quadratic[k][1]*grid.xy[i][0]+quadratic[k][2]*grid.xy[i][1]}};
          for (unsigned d=0;d<2;++d) squared+=std::pow(SU2_TYPE::GetValue(field.GetGradient_Reconstruction()(i,k,d))-exact[d],2);
        }
        norms.add(grid,i,SU2_TYPE::GetValue(fixture.geometry->nodes->GetVolume(i)),squared);
      }
      if (!distorted) REQUIRE(norms.rms(3)<1e-10); // Symmetric deep stencils reproduce quadratic gradients.
      Print("quadratic_gradient",distorted,method,"native",n,norms,previous);
      for (unsigned g=0;g<groups;++g) previous[g]=norms.rms(g);
    }
  }
}

namespace {
struct Segment { Direction q, s, a, b; };
std::vector<std::vector<Segment>> GeometricSegments(const Grid& grid,const CGeometry& geometry,const CAggregatedGeometry& metrics) {
  std::vector<std::vector<Segment>> segments(metrics.GetnEdge());
  for (const auto& triangle:grid.triangles) {
    Direction center{};
    for (auto node:triangle) for (unsigned d=0;d<2;++d) center[d]+=grid.xy[node][d]/3;
    for (unsigned k=0;k<3;++k) {
      const auto i=triangle[k], j=triangle[(k+1)%3];
      const long e=geometry.FindEdge(i,j,false);
      REQUIRE(e>=0); REQUIRE(static_cast<unsigned long>(e)<metrics.GetnEdge());
      const auto& f=metrics.GetFace(e);
      Direction s{{center[1]-f.point[1],-(center[0]-f.point[0])}};
      if (s[0]*(grid.xy[f.second][0]-grid.xy[f.first][0])+s[1]*(grid.xy[f.second][1]-grid.xy[f.first][1])<0)
        for (auto& value:s) value=-value;
      segments[e].push_back({{{0.5*(center[0]+f.point[0]),0.5*(center[1]+f.point[1])}},s,f.point,center});
    }
  }
  return segments;
}
}

TEST_CASE("Aggregated face moments are compared with exact geometric segment integrals", "[primitive-solver][primitive-order-audit]") {
  std::cout << "FACE_MOMENT,family,N,incidence,mean_linear_error,max_linear_error,mean_quadratic_error,mean_linear_over_h2,mean_quadratic_over_h3\n";
  for (bool distorted:{false,true}) for (unsigned n:{8u,16u,32u,64u,128u}) {
    Grid grid(n,distorted);
    MeshFixture fixture(grid.xy,grid.triangles,"PRIMITIVE_EULER",Options("WEIGHTED_LEAST_SQUARES","NODE_POINT"));
    CAggregatedGeometry metrics(*fixture.geometry,*fixture.config);
    const auto segments=GeometricSegments(grid,*fixture.geometry,metrics);
    std::array<double,2> linear{},maximum{},quadratic{},count{};
    for (unsigned long e=0;e<metrics.GetnEdge();++e) {
      const auto& f=metrics.GetFace(e);
      const unsigned group=segments[e].size()==1?0:1;
      REQUIRE(segments[e].size()>=1); REQUIRE(segments[e].size()<=2);
      Direction sum{};
      // First and second normal-weighted moments about the aggregated point.
      // Two Gauss points give exact integrals for these polynomial integrands.
      double m1[2][2]{},m2[2][3]{};
      for (const auto& segment:segments[e]) {
        for (unsigned d=0;d<2;++d) sum[d]+=segment.s[d];
        for (double xi:{-1/std::sqrt(3.0),1/std::sqrt(3.0)}) {
          const Direction offset{{segment.q[0]+0.5*xi*(segment.b[0]-segment.a[0])-f.point[0],
                                  segment.q[1]+0.5*xi*(segment.b[1]-segment.a[1])-f.point[1]}};
          for (unsigned d=0;d<2;++d) {
            for (unsigned a=0;a<2;++a) m1[d][a]+=0.5*segment.s[d]*offset[a];
            m2[d][0]+=0.5*segment.s[d]*offset[0]*offset[0];
            m2[d][1]+=0.5*segment.s[d]*offset[0]*offset[1];
            m2[d][2]+=0.5*segment.s[d]*offset[1]*offset[1];
          }
        }
      }
      for (unsigned d=0;d<2;++d) REQUIRE(std::abs(sum[d]-f.normal[d])<2e-14);
      double l=0,q=0;
      for (unsigned d=0;d<2;++d) {
        for (unsigned a=0;a<2;++a) l+=m1[d][a]*m1[d][a];
        for (unsigned a=0;a<3;++a) q+=m2[d][a]*m2[d][a];
      }
      l=std::sqrt(l); q=std::sqrt(q);
      linear[group]+=l; maximum[group]=std::max(maximum[group],l); quadratic[group]+=q; ++count[group];
    }
    for (unsigned g=0;g<2;++g) {
      REQUIRE(count[g]>0);
      const double l=linear[g]/count[g], q=quadratic[g]/count[g], h=1.0/n;
      std::cout << std::setprecision(12) << "FACE_MOMENT," << (distorted?"distorted":"regular") << ',' << n << ',' << g+1 << ','
                << l << ',' << maximum[g] << ',' << q << ',' << l/(h*h) << ',' << q/(h*h*h) << '\n';
      if (!distorted && g==1) CHECK(maximum[g]<2e-14);
      else CHECK(l/(h*h)>1e-4); // Counterexample to affine normal-moment exactness.
    }
  }
}

