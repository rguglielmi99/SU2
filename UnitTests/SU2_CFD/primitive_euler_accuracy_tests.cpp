/*!
 * \brief Spatial-order diagnostics. Alternative assemblies are TEST ONLY.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../SU2_CFD/include/gradients/computeGradientsLeastSquares.hpp"
#include <Eigen/QR>
#include <set>

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

namespace {
// Diagnostic only: a scaled two-ring quadratic fit supplies GRADIENTS ONLY.
// No fitted Hessian is stored, exposed to, or used by the production solver.
std::vector<Gradient> DiagnosticFitGradients(const Grid& grid, const CGeometry& geometry, bool weighted) {
  std::vector<Gradient> result(grid.xy.size());
  for (unsigned long i=0;i<grid.xy.size();++i) {
    std::set<unsigned long> stencil;
    for (auto j:geometry.nodes->GetPoints(i)) {
      stencil.insert(j);
      for (auto k:geometry.nodes->GetPoints(j)) stencil.insert(k);
    }
    stencil.erase(i);
    REQUIRE(stencil.size()>=5);
    double length=0;
    for (auto j:stencil) length=std::max(length,std::hypot(grid.xy[j][0]-grid.xy[i][0],grid.xy[j][1]-grid.xy[i][1]));
    Eigen::MatrixXd design(stencil.size(),5), data(stencil.size(),4);
    const auto center=W(grid.xy[i]);
    unsigned row=0;
    for (auto j:stencil) {
      const double x=(grid.xy[j][0]-grid.xy[i][0])/length, y=(grid.xy[j][1]-grid.xy[i][1])/length;
      const double weight=weighted?1/std::hypot(x,y):1;
      design.row(row)<<weight*x,weight*y,weight*0.5*x*x,weight*x*y,weight*0.5*y*y;
      const auto sample=W(grid.xy[j]);
      for (unsigned k=0;k<4;++k) data(row,k)=weight*(sample[k]-center[k]);
      ++row;
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
    REQUIRE(qr.rank()==5);
    const Eigen::MatrixXd coefficients=qr.solve(data);
    REQUIRE(coefficients.allFinite());
    for (unsigned k=0;k<4;++k) result[i][k]={{coefficients(0,k)/length,coefficients(1,k)/length}};
  }
  return result;
}

CPrimitiveEulerSolver::SpatialResult DiagnosticSegmentAssembly(const Grid& grid, const CAggregatedGeometry& metrics,
    const std::vector<std::vector<Segment>>& segments, const std::vector<Gradient>& gradient) {
  // Begin with the identical NODE_POINT volume, physical boundary and source.
  auto result=Assemble(grid,metrics,gradient,"NODE_POINT");
  std::fill(result.faces.begin(),result.faces.end(),State{});
  for (unsigned long e=0;e<metrics.GetnEdge();++e) {
    const auto& face=metrics.GetFace(e);
    for (const auto& segment:segments[e]) {
      // True median-dual segment, two-point Gauss quadrature, normal/weight once.
      for (double xi:{-1/std::sqrt(3.0),1/std::sqrt(3.0)}) {
        const Direction q{{segment.q[0]+0.5*xi*(segment.b[0]-segment.a[0]),
                           segment.q[1]+0.5*xi*(segment.b[1]-segment.a[1])}};
        const Direction normal{{0.5*segment.s[0],0.5*segment.s[1]}};
        const auto trace=[&](unsigned long i) {
          return Reconstruct(W(grid.xy[i]),gradient[i],{{q[0]-grid.xy[i][0],q[1]-grid.xy[i][1]}});
        };
        const auto waves=ComputeFluctuations(trace(face.first),trace(face.second),normal,1.4);
        for (unsigned k=0;k<4;++k) {
          result.faces[face.first][k]+=waves.minus[k]; result.faces[face.second][k]+=waves.plus[k];
        }
      }
    }
  }
  for (unsigned long i=0;i<grid.xy.size();++i) for (unsigned k=0;k<4;++k)
    result.rhs[i][k]=result.source[i][k]-(result.faces[i][k]+result.boundary[i][k]+result.volume[i][k])/metrics.GetCell(i).volume;
  return result;
}
}

TEST_CASE("Distorted MMS separates gradient accuracy from geometric segment integration", "[primitive-solver][primitive-irregular-audit]") {
  std::cout<<"IRREGULAR_MMS,family,gradient,experiment,N,group,volume,L2,order\n";
  std::cout<<"IRREGULAR_GRAD,family,gradient,experiment,N,group,volume,L2,order\n";
  for (const std::string family:{"regular","distorted","smooth_distorted"})
    for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"}) {
      std::map<std::string,std::array<double,groups>> previous;
      for (unsigned n:{8u,16u,32u,64u,128u}) {
        Grid grid(n,family=="distorted");
        if (family=="smooth_distorted") for (auto& x:grid.xy) {
          // A fixed smooth coordinate map, distinct from grid-scale jitter.
          const auto original=x;
          x[0]+=0.035*std::sin(2*3.141592653589793*original[0])*std::sin(2*3.141592653589793*original[1]);
          x[1]+=0.03*std::sin(2*3.141592653589793*original[0])*std::sin(3.141592653589793*original[1]);
        }
        MeshFixture fixture(grid.xy,grid.triangles,"PRIMITIVE_EULER",Options(method,"NODE_POINT"));
        CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
        for (unsigned long i=0;i<grid.xy.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,W(grid.xy[i]));
        const auto actual=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactBoundary,exactSource);
        const auto& metrics=solver.GetGeometryMetrics();
        const auto segments=GeometricSegments(grid,*fixture.geometry,metrics);
        for (unsigned long e=0;e<metrics.GetnEdge();++e) {
          Direction sum{};
          for (const auto& segment:segments[e]) for (unsigned d=0;d<2;++d) sum[d]+=segment.s[d];
          for (unsigned d=0;d<2;++d) REQUIRE(std::abs(sum[d]-metrics.GetFace(e).normal[d])<2e-14);
        }
        std::vector<Gradient> native(grid.xy.size()), exact(grid.xy.size());
        double minAreaRatio=1e10;
        for (const auto& t:grid.triangles) {
          const auto& a=grid.xy[t[0]];const auto& b=grid.xy[t[1]];const auto& c=grid.xy[t[2]];
          minAreaRatio=std::min(minAreaRatio,((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))*n*n);
        }
        REQUIRE(minAreaRatio>0.2); // Families stay shape regular; no collapsed triangles.
        for (unsigned long i=0;i<grid.xy.size();++i) {
          REQUIRE(solver.GetPositivityScale(i)==1);
          native[i]=solver.GetLimitedGradient(i);exact[i]=G(grid.xy[i]);
        }
        const auto fit=DiagnosticFitGradients(grid,*fixture.geometry,method=="WEIGHTED_LEAST_SQUARES");
        for (const std::string type:{"native_gradient","fit_gradient"}) {
          Norms norms; const auto& measured=type=="native_gradient"?native:fit;
          for (unsigned long i=0;i<grid.xy.size();++i) {
            double squared=0;
            for (unsigned k=0;k<4;++k) for (unsigned d=0;d<2;++d) squared+=std::pow(measured[i][k][d]-exact[i][k][d],2);
            norms.add(grid,i,metrics.GetCell(i).volume,squared);
          }
          for (unsigned g=0;g<groups;++g) {
            const double error=norms.rms(g),before=previous[type][g];
            REQUIRE(std::isfinite(error));
            std::cout<<std::setprecision(12)<<"IRREGULAR_GRAD,"<<family<<','<<method<<','<<type<<','<<n<<','
                     <<names[g]<<','<<norms.weight[g]<<','<<error<<','<<(before>0?std::log2(before/error):0)<<'\n';
            previous[type][g]=error;
          }
        }
        std::map<std::string,CPrimitiveEulerSolver::SpatialResult> experiments;
        experiments.emplace("native_aggregated",actual);
        experiments.emplace("exact_aggregated",Assemble(grid,metrics,exact,"NODE_POINT"));
        experiments.emplace("fit_gradient_aggregated",Assemble(grid,metrics,fit,"NODE_POINT"));
        experiments.emplace("exact_segmented",DiagnosticSegmentAssembly(grid,metrics,segments,exact));
        experiments.emplace("native_segmented",DiagnosticSegmentAssembly(grid,metrics,segments,native));
        const auto reproduced=Assemble(grid,metrics,native,"NODE_POINT");
        double discrepancy=0;
        for (unsigned long i=0;i<grid.xy.size();++i) for (unsigned k=0;k<4;++k)
          discrepancy=std::max(discrepancy,std::abs(reproduced.rhs[i][k]-actual.rhs[i][k]));
        REQUIRE(discrepancy<2e-12);
        for (const auto& experiment:experiments) {
          Norms norms;
          for (unsigned long i=0;i<grid.xy.size();++i) norms.add(grid,i,metrics.GetCell(i).volume,Squared(experiment.second.rhs[i]));
          for (unsigned g=0;g<groups;++g) {
            const double error=norms.rms(g), before=previous[experiment.first][g];
            REQUIRE(std::isfinite(error));
            std::cout<<std::setprecision(12)<<"IRREGULAR_MMS,"<<family<<','<<method<<','<<experiment.first<<','<<n<<','
                     <<names[g]<<','<<norms.weight[g]<<','<<error<<','<<(before>0?std::log2(before/error):0)<<'\n';
            previous[experiment.first][g]=error;
          }
        }
      }
    }
}

TEST_CASE("First order contact separates centered dissipative geometric and boundary moments", "[primitive-solver][primitive-irregular-audit]") {
  const auto contact=[](const Direction& x){return State{{2+0.1*x[0]+0.2*x[1],0.5,0.3,3}};};
  const State source{{0.11,0,0,0}};
  std::cout<<"CONTACT_SPLIT,patch,h,V,centered,dissipative,total,source,rhs,integrated_rhs,segmented,exact_G_rhs,Mxx,Mxy,Myx,Myy\n";
  for (const std::string patch:{"symmetric_fan","asymmetric_fan","asymmetric_six"}) {
    std::vector<Direction> original;
    std::vector<Triangle> cells;
    unsigned long center;
    if (patch=="asymmetric_six") {
      Grid grid(4,true); original=grid.xy;cells=grid.triangles;center=12;
    } else {
      original={{{0,0}},{{2,0}},{{2,1}},{{0,1}},patch=="symmetric_fan"?Direction{{1,0.5}}:Direction{{0.8,0.3}}};
      cells={{{0,1,4}},{{1,2,4}},{{2,3,4}},{{3,0,4}}}; center=4;
    }
    double first=0;
    for (double h:{1.0,0.5,0.25,0.125}) {
      Grid grid(2,false);grid.xy=original;grid.triangles=cells;
      for (auto& x:grid.xy) for (unsigned d=0;d<2;++d) x[d]=h*(x[d]-original[center][d]);
      auto options=Options("WEIGHTED_LEAST_SQUARES","NODE_POINT");
      options.replace(options.find("MUSCL_FLOW=YES"),14,"MUSCL_FLOW=NO");
      MeshFixture fixture(grid.xy,cells,"PRIMITIVE_EULER",options);
      CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
      for (unsigned long i=0;i<grid.xy.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,contact(grid.xy[i]));
      const auto result=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
          [&](const Direction& x,double,const State&,const Direction&,unsigned short){return contact(x);},
          [&](const Direction&,double){return source;});
      const auto& metrics=solver.GetGeometryMetrics();
      std::vector<State> centered(grid.xy.size()), dissipative(grid.xy.size());
      const auto scatter=[&](unsigned long i,unsigned long j,const Direction& s) {
        const auto left=contact(grid.xy[i]),right=contact(grid.xy[j]); State jump{};
        for (unsigned k=0;k<4;++k) jump[k]=right[k]-left[k];
        const auto path=Apply(PathMatrix(left,right,s,1.4),jump), visc=Apply(AbsolutePathMatrix(left,right,s,1.4),jump);
        for (unsigned k=0;k<4;++k) {
          centered[i][k]+=0.5*path[k];centered[j][k]+=0.5*path[k];
          dissipative[i][k]-=0.5*visc[k];dissipative[j][k]+=0.5*visc[k];
        }
      };
      double moment[2][2]{};
      for (unsigned long e=0;e<metrics.GetnEdge();++e) {
        const auto& face=metrics.GetFace(e);scatter(face.first,face.second,face.normal);
        if (face.first==center || face.second==center) {
          const auto other=face.first==center?face.second:face.first;
          const double sign=face.first==center?1:-1;
          for (unsigned d=0;d<2;++d) for (unsigned a=0;a<2;++a)
            moment[d][a]+=0.5*sign*face.normal[d]*(grid.xy[other][a]-grid.xy[center][a])/metrics.GetCell(center).volume;
        }
      }
      for (const auto& face:metrics.GetBoundaryFaces()) {
        const auto inside=contact(grid.xy[face.node]),outside=contact(face.point);State jump{};
        for (unsigned k=0;k<4;++k) jump[k]=outside[k]-inside[k];
        const auto path=Apply(PathMatrix(inside,outside,face.normal,1.4),jump),visc=Apply(AbsolutePathMatrix(inside,outside,face.normal,1.4),jump);
        for (unsigned k=0;k<4;++k) {centered[face.node][k]+=0.5*path[k];dissipative[face.node][k]-=0.5*visc[k];}
      }
      for (unsigned long i=0;i<grid.xy.size();++i) for (unsigned k=0;k<4;++k) {
        CHECK(result.volume[i][k]==0);
        CHECK(result.rhs[i][k]==Approx(source[k]-(centered[i][k]+dissipative[i][k])/metrics.GetCell(i).volume).margin(2e-12));
      }
      const auto segments=GeometricSegments(grid,*fixture.geometry,metrics);
      double segmented=0;
      for (unsigned long e=0;e<metrics.GetnEdge();++e) {
        const auto& face=metrics.GetFace(e);
        if (face.first!=center && face.second!=center) continue;
        const auto other=face.first==center?face.second:face.first;const double sign=face.first==center?1:-1;
        for (const auto& segment:segments[e]) {
          const double speed=sign*(0.5*segment.s[0]+0.3*segment.s[1]);
          segmented+=std::min(speed,0.0)*(contact(grid.xy[other])[0]-contact(grid.xy[center])[0]);
        }
      }
      const double v=metrics.GetCell(center).volume,c=centered[center][0]/v,d=dissipative[center][0]/v;
      CHECK(result.boundary[center]==State{});
      CHECK(c==Approx(0.5*(0.1*moment[0][0]+0.2*moment[0][1])+0.3*(0.1*moment[1][0]+0.2*moment[1][1])).margin(1e-13));
      for (unsigned a=0;a<2;++a) for (unsigned b=0;b<2;++b) CHECK(moment[a][b]==Approx(a==b?1.0:0.0).margin(1e-12));
      if (h==1) first=result.rhs[center][0];else CHECK(result.rhs[center][0]==Approx(first).margin(2e-12));
      // Exact contact gradient with affine MUSCL gives matching traces at ANY Q,
      // and the NODE_POINT volume is the exact 0.11; include full source/boundary.
      Gradient exact{};exact[0]={{0.1,0.2}};
      double exactRhs=0;
      for (unsigned long e=0;e<metrics.GetnEdge();++e) {
        const auto& f=metrics.GetFace(e);
        const auto trace=[&](unsigned long i){return Reconstruct(contact(grid.xy[i]),exact,{{f.point[0]-grid.xy[i][0],f.point[1]-grid.xy[i][1]}});};
        const auto waves=ComputeFluctuations(trace(f.first),trace(f.second),f.normal,1.4);
        if (f.first==center) exactRhs-=waves.minus[0]/v;
        if (f.second==center) exactRhs-=waves.plus[0]/v;
      }
      exactRhs+=source[0]-Product(contact(grid.xy[center]),exact)[0];
      CHECK(std::abs(exactRhs)<2e-12);
      std::cout<<std::setprecision(12)<<"CONTACT_SPLIT,"<<patch<<','<<h<<','<<v<<','<<c<<','<<d<<','<<c+d<<",0.11,"
               <<result.rhs[center][0]<<','<<v*result.rhs[center][0]<<','<<segmented/v<<','<<exactRhs;
      for (const auto& row:moment) for (auto value:row) std::cout<<','<<value;
      std::cout<<'\n';
    }
  }
}

TEST_CASE("First order affine contact measures strong and weak consistency on fixed domains", "[primitive-solver][primitive-irregular-audit]") {
  std::cout<<"CONTACT_WEAK,family,N,L2_deep,weighted_signed_defect,weak_zero_trace,integrated_L1\n";
  const auto contact=[](const Direction& x){return State{{2+0.1*x[0]+0.2*x[1],0.5,0.3,3}};};
  for (const std::string family:{"regular","distorted","smooth_distorted"}) {
    for (unsigned n:{8u,16u,32u,64u,128u}) {
      Grid grid(n,family=="distorted");
      if (family=="smooth_distorted") for (auto& x:grid.xy) {
        const auto old=x;
        x[0]+=0.035*std::sin(2*3.141592653589793*old[0])*std::sin(2*3.141592653589793*old[1]);
        x[1]+=0.03*std::sin(2*3.141592653589793*old[0])*std::sin(3.141592653589793*old[1]);
      }
      auto options=Options("WEIGHTED_LEAST_SQUARES","NODE_POINT");
      options.replace(options.find("MUSCL_FLOW=YES"),14,"MUSCL_FLOW=NO");
      MeshFixture fixture(grid.xy,grid.triangles,"PRIMITIVE_EULER",options);
      CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
      for (unsigned long i=0;i<grid.xy.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,contact(grid.xy[i]));
      const auto result=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
          [&](const Direction& x,double,const State&,const Direction&,unsigned short){return contact(x);},
          [](const Direction&,double){return State{{0.11,0,0,0}};});
      Norms norms;double signedDefect=0,weak=0,l1=0;
      for (unsigned long i=0;i<grid.xy.size();++i) {
        const double v=solver.GetGeometryMetrics().GetCell(i).volume,e=result.rhs[i][0];
        const auto& x=grid.xy[i];
        const double phi=256*std::pow(x[0]*(1-x[0])*x[1]*(1-x[1]),2);
        norms.add(grid,i,v,e*e);signedDefect+=v*e;weak+=v*e*phi;l1+=v*std::abs(e);
        for (unsigned k=1;k<4;++k) CHECK(std::abs(result.rhs[i][k])<1e-12);
      }
      CHECK(std::isfinite(weak));
      std::cout<<std::setprecision(12)<<"CONTACT_WEAK,"<<family<<','<<n<<','<<norms.rms(3)<<','<<signedDefect<<','<<weak<<','<<l1<<'\n';
    }
  }
}

namespace {
State ExtremumW(const Direction& x) {
  State w{{2,0.5,0.3,3}};
  const double r=(x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5);
  for (unsigned k=0;k<4;++k) w[k]+=(k+1)*0.1*r;
  return w;
}
Gradient ExtremumG(const Direction& x) {
  Gradient g{};
  for (unsigned k=0;k<4;++k) for (unsigned d=0;d<2;++d) g[k][d]=2*(k+1)*0.1*(x[d]-0.5);
  return g;
}
std::string LimiterOptions(const std::string& method,const std::string& limiter) {
  auto options=Options(method,"NODE_POINT");
  options.replace(options.find("SLOPE_LIMITER_FLOW=NONE"),std::string("SLOPE_LIMITER_FLOW=NONE").size(),"SLOPE_LIMITER_FLOW="+limiter);
  return options;
}
void LimiterActivity(const Grid& grid,const CPrimitiveEulerSolver& solver,const std::string& family,
    const std::string& method,const std::string& limiter,const std::string& field) {
  for (unsigned group=0;group<groups;++group) for (unsigned component=0;component<5;++component) {
    std::vector<double> values;
    double active=0,lt99=0,lt90=0,lt50=0,above=0,deviation=0,thetaActive=0,minTheta=1;
    for (unsigned long i=0;i<grid.xy.size();++i) if (grid.in(i,group)) {
      double phi=component<4?SU2_TYPE::GetValue(solver.GetPrimitiveNodes().GetLimiter(i,component)):
                            std::numeric_limits<double>::infinity();
      if (component==4) for (unsigned k=0;k<4;++k) phi=std::min(phi,SU2_TYPE::GetValue(solver.GetPrimitiveNodes().GetLimiter(i,k)));
      REQUIRE(std::isfinite(phi));REQUIRE(phi>=0);
      values.push_back(phi);active+=phi<1-1e-8;lt99+=phi<0.99;lt90+=phi<0.9;lt50+=phi<0.5;
      above+=phi>1+1e-8;deviation+=std::abs(1-phi);
      const double theta=solver.GetPositivityScale(i);
      thetaActive+=theta<1-1e-12;minTheta=std::min(minTheta,theta);
    }
    std::sort(values.begin(),values.end());REQUIRE_FALSE(values.empty());
    const double count=values.size();
    const auto quantile=[&](double p){return values[static_cast<size_t>(p*(values.size()-1))];};
    std::cout<<std::setprecision(12)<<"LIMITER_ACTIVITY,"<<family<<','<<method<<','<<limiter<<','<<field<<','<<grid.n<<','<<names[group]<<','<<component<<','
             <<values.size()<<','<<active/count<<','<<values.front()<<','<<quantile(0.1)<<','<<quantile(0.5)<<','<<quantile(0.9)<<','
             <<lt99/count<<','<<lt90/count<<','<<lt50/count<<','<<above/count<<','<<deviation/count<<','<<thetaActive/count<<','<<minTheta<<'\n';
  }
}
}

TEST_CASE("Production affine MUSCL limiter MMS and activity are measured without imposed order", "[primitive-solver][primitive-limiter-audit]") {
  std::cout<<"LIMITER_MMS,family,gradient,limiter,field,N,group,L2,order,gradient_L2\n";
  std::cout<<"LIMITER_EXTREMA,family,gradient,limiter,N,nodes,local_L2,order,min_phi,center_phi_rho,center_phi_u,center_phi_v,center_phi_p\n";
  std::cout<<"LIMITER_ACTIVITY,family,gradient,limiter,field,N,group,component,nodes,f_active,min_phi,p10,p50,p90,f_lt099,f_lt09,f_lt05,f_gt1,mean_abs_1_phi,f_theta,min_theta\n";
  for (const std::string family:{"regular","smooth_distorted"})
    for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"})
      for (const std::string limiter:{"NONE","VENKATAKRISHNAN","BARTH_JESPERSEN"})
        for (const std::string field:{"mms","extremum"}) {
          std::array<double,groups> previous{};
          double previousCore=0;
          for (unsigned n:{8u,16u,32u,64u,128u}) {
            Grid grid(n,false);
            if (family=="smooth_distorted") for (auto& x:grid.xy) {
              const auto old=x;
              x[0]+=0.035*std::sin(2*3.141592653589793*old[0])*std::sin(2*3.141592653589793*old[1]);
              x[1]+=0.03*std::sin(2*3.141592653589793*old[0])*std::sin(3.141592653589793*old[1]);
            }
            const auto state=[&](const Direction& x){return field=="mms"?W(x):ExtremumW(x);};
            const auto gradient=[&](const Direction& x){return field=="mms"?G(x):ExtremumG(x);};
            MeshFixture fixture(grid.xy,grid.triangles,"PRIMITIVE_EULER",LimiterOptions(method,limiter));
            CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
            for (unsigned long i=0;i<grid.xy.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,state(grid.xy[i]));
            const auto result=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
                [&](const Direction& x,double,const State&,const Direction&,unsigned short){return state(x);},
                [&](const Direction& x,double){return Product(state(x),gradient(x));});
            Norms norms,gradientNorms;
            for (unsigned long i=0;i<grid.xy.size();++i) {
              const double v=solver.GetGeometryMetrics().GetCell(i).volume;
              norms.add(grid,i,v,Squared(result.rhs[i]));
              double squared=0;const auto expected=gradient(grid.xy[i]);
              for (unsigned k=0;k<4;++k) for (unsigned d=0;d<2;++d)
                squared+=std::pow(solver.GetLimitedGradient(i)[k][d]-expected[k][d],2);
              gradientNorms.add(grid,i,v,squared);
              REQUIRE(solver.GetPositivityScale(i)==1);
            }
            for (unsigned group=0;group<groups;++group) {
              const double error=norms.rms(group);REQUIRE(std::isfinite(error));
              std::cout<<std::setprecision(12)<<"LIMITER_MMS,"<<family<<','<<method<<','<<limiter<<','<<field<<','<<n<<','<<names[group]<<','
                       <<error<<','<<(n==8?0:std::log2(previous[group]/error))<<','<<gradientNorms.rms(group)<<'\n';
              previous[group]=error;
            }
            if (field=="extremum") {
              double squared=0,weight=0,minPhi=1;unsigned count=0;
              for (unsigned j=n/2-2;j<=n/2+2;++j) for (unsigned i=n/2-2;i<=n/2+2;++i) {
                const unsigned long point=j*(n+1)+i;const double v=solver.GetGeometryMetrics().GetCell(point).volume;
                squared+=v*Squared(result.rhs[point]);weight+=v;++count;
                for (unsigned k=0;k<4;++k) minPhi=std::min(minPhi,SU2_TYPE::GetValue(solver.GetPrimitiveNodes().GetLimiter(point,k)));
              }
              const double error=std::sqrt(squared/weight);REQUIRE(std::isfinite(error));
              std::cout<<std::setprecision(12)<<"LIMITER_EXTREMA,"<<family<<','<<method<<','<<limiter<<','<<n<<','<<count<<','
                       <<error<<','<<(n==8?0:std::log2(previousCore/error))<<','<<minPhi;
              const unsigned long center=(n/2)*(n+1)+n/2;
              for (unsigned k=0;k<4;++k) std::cout<<','<<SU2_TYPE::GetValue(solver.GetPrimitiveNodes().GetLimiter(center,k));
              std::cout<<'\n';previousCore=error;
            }
            LimiterActivity(grid,solver,family,method,limiter,field);
            const auto valid=[](const State& w){REQUIRE(w[0]>0);REQUIRE(w[3]>0);};
            for (unsigned long e=0;e<solver.GetGeometryMetrics().GetnEdge();++e) {
              const auto& f=solver.GetGeometryMetrics().GetFace(e);
              valid(solver.ReconstructAt(f.first,f.point));valid(solver.ReconstructAt(f.second,f.point));
            }
            for (const auto& f:solver.GetGeometryMetrics().GetBoundaryFaces()) valid(solver.ReconstructAt(f.node,f.point));
          }
        }
}

TEST_CASE("All production limiters preserve constants and admissible affine data", "[primitive-solver][primitive-limiter-audit]") {
  const Gradient g{{{{0.2,-0.1}},{{0.12,0.07}},{{-0.08,0.09}},{{0.3,0.11}}}};
  const auto affine=[&](const Direction& x){return Reconstruct({{2,0.5,0.3,3}},g,x);};
  for (bool distorted:{false,true}) for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"})
    for (const std::string limiter:{"NONE","VENKATAKRISHNAN","BARTH_JESPERSEN"}) {
      Grid grid(16,distorted);
      MeshFixture fixture(grid.xy,grid.triangles,"PRIMITIVE_EULER",LimiterOptions(method,limiter));
      CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
      const auto constant=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
          [](const Direction&,double,const State& w,const Direction&,unsigned short){return w;});
      for (const auto& rhs:constant.rhs) for (auto value:rhs) CHECK(std::abs(value)<1e-12);
      for (unsigned long i=0;i<grid.xy.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,affine(grid.xy[i]));
      const auto result=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
          [&](const Direction& x,double,const State&,const Direction&,unsigned short){return affine(x);},
          [&](const Direction& x,double){return Product(affine(x),g);});
      for (unsigned long i=0;i<grid.xy.size();++i) {
        CHECK(solver.GetPositivityScale(i)==1);
        for (unsigned k=0;k<4;++k) {
          CHECK(SU2_TYPE::GetValue(solver.GetPrimitiveNodes().GetLimiter(i,k))==Approx(1).margin(2e-11));
          CHECK(std::abs(result.rhs[i][k])<2e-11);
          for (unsigned d=0;d<2;++d) CHECK(solver.GetLimitedGradient(i)[k][d]==Approx(g[k][d]).margin(2e-12));
        }
      }
    }
}

TEST_CASE("Active production limiters and positivity survive renumbering and boundary reversal", "[primitive-solver][primitive-limiter-audit]") {
  std::cout<<"LIMITER_STRESS,gradient,limiter,f_theta,min_theta,min_phi\n";
  Grid grid(4,true);std::vector<Direction> permuted(grid.xy.size());
  const auto permutation=[](unsigned long i){return (7*i+3)%25;};
  for (unsigned long i=0;i<grid.xy.size();++i) permuted[permutation(i)]=grid.xy[i];
  auto cells=grid.triangles;
  for (auto& t:cells) for (auto& i:t) i=permutation(i);
  std::reverse(cells.begin(),cells.end());
  for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"})
    for (const std::string limiter:{"NONE","VENKATAKRISHNAN","BARTH_JESPERSEN"}) {
      MeshFixture fa(grid.xy,grid.triangles,"PRIMITIVE_EULER",LimiterOptions(method,limiter));
      MeshFixture fb(permuted,cells,"PRIMITIVE_EULER",LimiterOptions(method,limiter));
      for (unsigned m=0;m<fb.geometry->GetnMarker();++m)
        for (unsigned long e=0;e<fb.geometry->GetnElem_Bound(m);++e) fb.geometry->bound[m][e]->Change_Orientation();
      CPrimitiveEulerSolver a(fa.geometry.get(),fa.config.get(),0),b(fb.geometry.get(),fb.config.get(),0);
      std::vector<State> initial(grid.xy.size());
      for (unsigned long i=0;i<grid.xy.size();++i) {
        initial[i]={{i==7?100.0:1.0,0.5,0.3,i==17?100.0:1.0}};
        a.GetPrimitiveNodes().SetPrimitiveState(i,initial[i]);b.GetPrimitiveNodes().SetPrimitiveState(permutation(i),initial[i]);
      }
      const auto extrapolate=[](const Direction&,double,const State& w,const Direction&,unsigned short){return w;};
      const auto ra=a.ComputeSpatialResidual(fa.geometry.get(),fa.config.get(),0,extrapolate);
      const auto rb=b.ComputeSpatialResidual(fb.geometry.get(),fb.config.get(),0,extrapolate);
      double active=0,minTheta=1,minPhi=1;
      for (unsigned long i=0;i<grid.xy.size();++i) {
        CHECK(a.GetPrimitiveNodes().GetPrimitiveState(i)==initial[i]);
        CHECK(a.GetPositivityScale(i)==Approx(b.GetPositivityScale(permutation(i))).margin(2e-12));
        active+=a.GetPositivityScale(i)<1-1e-12;minTheta=std::min(minTheta,a.GetPositivityScale(i));
        for (unsigned k=0;k<4;++k) {
          const double phi=SU2_TYPE::GetValue(a.GetPrimitiveNodes().GetLimiter(i,k));minPhi=std::min(minPhi,phi);
          CHECK(phi==Approx(SU2_TYPE::GetValue(b.GetPrimitiveNodes().GetLimiter(permutation(i),k))).epsilon(2e-10).margin(2e-12));
          CHECK(ra.rhs[i][k]==Approx(rb.rhs[permutation(i)][k]).epsilon(2e-10).margin(2e-10));
          for (unsigned d=0;d<2;++d)
            CHECK(a.GetLimitedGradient(i)[k][d]==Approx(b.GetLimitedGradient(permutation(i))[k][d]).epsilon(2e-10).margin(2e-11));
        }
      }
      const auto positive=[](const State& w){CHECK(w[0]>0);CHECK(w[3]>0);};
      for (unsigned long e=0;e<a.GetGeometryMetrics().GetnEdge();++e) {
        const auto& f=a.GetGeometryMetrics().GetFace(e);
        for (auto i:{f.first,f.second}) {positive(a.ReconstructAt(i,f.point));positive(b.ReconstructAt(permutation(i),f.point));}
      }
      for (const auto& f:a.GetGeometryMetrics().GetBoundaryFaces()) {
        positive(a.ReconstructAt(f.node,f.point));positive(b.ReconstructAt(permutation(f.node),f.point));
      }
      if (limiter=="NONE") CHECK(active>0);else CHECK(minPhi<0.99);
      std::cout<<std::setprecision(12)<<"LIMITER_STRESS,"<<method<<','<<limiter<<','<<active/grid.xy.size()<<','<<minTheta<<','<<minPhi<<'\n';
    }
}
