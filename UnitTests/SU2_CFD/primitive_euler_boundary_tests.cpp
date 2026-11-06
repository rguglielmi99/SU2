/*!
 * \brief Physical boundary invariants and constrained-slip operator tests.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../SU2_CFD/include/numerics/flow/primitive_euler_boundary.hpp"
#include "../../SU2_CFD/include/gradients/computeGradientsLeastSquares.hpp"
#include "../../SU2_CFD/include/gradients/computeGradientsGreenGauss.hpp"
#include <set>

using namespace PrimitiveEuler2D;
using namespace PrimitiveEulerTest;
namespace {
constexpr double gasGamma=1.4;
auto Precise(double value) { return Approx(value).epsilon(64*std::numeric_limits<double>::epsilon()).margin(2e-13); }
void Same(const State& a,const State& b,double tolerance=2e-12) {
  for (unsigned k=0;k<4;++k) CHECK(a[k]==Precise(b[k]).margin(tolerance));
}
State Oriented(double un,double ut,const Direction& n) { return {{1,un*n[0]-ut*n[1],un*n[1]+ut*n[0],1}}; }
double Plus(const State& w,const Direction& n) { return Boundary::NormalVelocity(w,n)+2*Boundary::Sound(w,gasGamma)/.4; }
double Minus(const State& w,const Direction& n) { return Boundary::NormalVelocity(w,n)-2*Boundary::Sound(w,gasGamma)/.4; }
}

namespace {
using GlobalFive=std::array<double,5>; // M, Px, Py, E, H; test-only diagnostics.
GlobalFive GlobalDensity(const State& w) {
  const double speed2=w[1]*w[1]+w[2]*w[2];
  return {{w[0],w[0]*w[1],w[0]*w[2],w[3]/.4+.5*w[0]*speed2,
           w[0]*(std::log(w[3])-1.4*std::log(w[0]))}};
}
GlobalFive GlobalChain(const State& w,const State& r) {
  return {{r[0],w[1]*r[0]+w[0]*r[1],w[2]*r[0]+w[0]*r[2],
           r[3]/.4+.5*(w[1]*w[1]+w[2]*w[2])*r[0]+w[0]*(w[1]*r[1]+w[2]*r[2]),
           (std::log(w[3])-1.4*std::log(w[0])-1.4)*r[0]+w[0]*r[3]/w[3]}};
}
struct BalanceBox {
  unsigned n;
  std::vector<Direction> xy;
  std::vector<Triangle> cells;
  std::vector<BoundaryMarker> markers{{"south",{}},{"east",{}},{"north",{}},{"west",{}}};
  explicit BalanceBox(unsigned count):n(count) {
    const auto id=[&](unsigned i,unsigned j){return static_cast<unsigned long>(j*(n+1)+i);};
    for (unsigned j=0;j<=n;++j) for (unsigned i=0;i<=n;++i) xy.push_back({{double(i)/n,double(j)/n}});
    for (unsigned j=0;j<n;++j) for (unsigned i=0;i<n;++i) {
      const auto a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);
      cells.push_back({{a,b,c}}); cells.push_back({{a,c,d}});
    }
    for (unsigned k=0;k<n;++k) {
      markers[0].edges.push_back({{id(k,0),id(k+1,0)}});
      markers[1].edges.push_back({{id(n,k),id(n,k+1)}});
      markers[2].edges.push_back({{id(k+1,n),id(k,n)}});
      markers[3].edges.push_back({{id(0,k+1),id(0,k)}});
    }
  }
};
const std::string balanceOptions="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
    "MUSCL_FLOW=YES\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\nSLOPE_LIMITER_FLOW=NONE\n"
    "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\n";
State BoxInitial(const Direction& x) {
  const double pi=std::acos(-1.),cx=std::cos(pi*x[0]),cy=std::cos(pi*x[1]);
  return {{1+.03*cx*cy,.1*std::sin(pi*x[0])*cy,-.1*cx*std::sin(pi*x[1]),1+.04*cx+.03*cy}};
}
struct GlobalBalance {
  GlobalFive value{},rate{};
  Direction force{},closure{};
  double minRho=1e300,minP=1e300,slip=0,maxRhs=0;
};
GlobalBalance MeasureGlobal(CPrimitiveEulerSolver& s,const CPrimitiveEulerSolver::SpatialResult& result) {
  GlobalBalance b;
  for (unsigned long i=0;i<result.rhs.size();++i) {
    const auto w=s.GetPrimitiveNodes().GetPrimitiveState(i);
    REQUIRE(w[0]>0); REQUIRE(w[3]>0);
    const auto value=GlobalDensity(w),rate=GlobalChain(w,result.rhs[i]);
    const double volume=s.GetGeometryMetrics().GetCell(i).volume;
    for (unsigned k=0;k<5;++k) {
      REQUIRE(std::isfinite(value[k])); REQUIRE(std::isfinite(rate[k]));
      b.value[k]+=volume*value[k]; b.rate[k]+=volume*rate[k];
    }
    for (auto r:result.rhs[i]) b.maxRhs=std::max(b.maxRhs,std::abs(r));
    b.minRho=std::min(b.minRho,w[0]); b.minP=std::min(b.minP,w[3]);
  }
  for (const auto& face:s.GetGeometryMetrics().GetBoundaryFaces()) {
    const auto trace=s.ReconstructAt(face.node,face.point);
    const auto w=s.GetPrimitiveNodes().GetPrimitiveState(face.node);
    const double length=std::hypot(face.normal[0],face.normal[1]);
    b.slip=std::max(b.slip,std::abs(w[1]*face.normal[0]+w[2]*face.normal[1])/length);
    for (unsigned d=0;d<2;++d) {
      b.force[d]+=trace[3]*face.normal[d]; b.closure[d]+=face.normal[d];
    }
  }
  CHECK(b.slip<2e-13);
  for (auto value:b.closure) CHECK(value==Approx(0).margin(2e-14));
  return b;
}
void PrintGlobal(const char* kind,unsigned n,unsigned step,double time,double dt,const GlobalBalance& b,
                 const GlobalFive& initial,const GlobalFive& integrated={},const Direction& impulse={}) {
  std::cout << std::setprecision(17) << "GLOBAL_BALANCE," << kind << ',' << n << ',' << step << ',' << time << ',' << dt;
  for (auto value:b.value) std::cout << ',' << value;
  for (auto value:b.rate) std::cout << ',' << value;
  std::cout << ',' << b.force[0] << ',' << b.force[1] << ',' << (b.value[0]-initial[0])/initial[0]
            << ',' << (b.value[3]-initial[3])/initial[3] << ',' << b.value[4]-initial[4]
            << ',' << b.value[1]-initial[1]+impulse[0] << ',' << b.value[2]-initial[2]+impulse[1]
            << ',' << impulse[0] << ',' << impulse[1];
  for (auto value:integrated) std::cout << ',' << value;
  for (unsigned k=0;k<5;++k) std::cout << ',' << b.value[k]-initial[k]-integrated[k];
  std::cout << ',' << b.minRho << ',' << b.minP << ',' << b.slip << ',' << b.maxRhs << '\n';
}
}

TEST_CASE("Global balance chain diagnostics agree with independent directional derivatives", "[primitive-global-balance]") {
  const State w{{1.13,.27,-.18,1.21}},direction{{.07,-.11,.13,-.09}};
  const auto rate=GlobalChain(w,direction);
  for (double epsilon:{1e-4,1e-5}) {
    auto plus=w,minus=w;
    for (unsigned k=0;k<4;++k) { plus[k]+=epsilon*direction[k]; minus[k]-=epsilon*direction[k]; }
    const auto a=GlobalDensity(plus),b=GlobalDensity(minus);
    for (unsigned k=0;k<5;++k) CHECK(rate[k]==Approx((a[k]-b[k])/(2*epsilon)).epsilon(2e-8).margin(5e-11));
  }
}

TEST_CASE("Closed square measures equilibrium and nonuniform semidiscrete physical balances", "[primitive-global-balance]") {
  std::cout << "GLOBAL_BALANCE,kind,N,step,time,dt,M,Px,Py,E,H,dM,dPx,dPy,dE,dH,Fx,Fy,M_relative,E_relative,H_delta,Px_impulse_defect,Py_impulse_defect,impulse_x,impulse_y,int_dM,int_dPx,int_dPy,int_dE,int_dH,time_M,time_Px,time_Py,time_E,time_H,min_rho,min_p,slip,max_rhs\n";
  for (unsigned n:{8u,16u,32u}) {
    BalanceBox grid(n);
    MeshFixture f(grid.xy,grid.cells,"PRIMITIVE_EULER",balanceOptions,"MARKER_EULER",grid.markers);
    CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0);
    if (n==8) {
      for (unsigned long i=0;i<grid.xy.size();++i) s.GetPrimitiveNodes().SetPrimitiveState(i,{{1,0,0,1}});
      const auto r=s.ComputeSpatialResidual(f.geometry.get(),f.config.get());
      const auto quiet=MeasureGlobal(s,r);
      CHECK(quiet.maxRhs<2e-13); CHECK(quiet.value[0]==Precise(1)); CHECK(quiet.value[3]==Precise(2.5));
      for (auto rate:quiet.rate) CHECK(rate==Approx(0).margin(2e-13));
      PrintGlobal("quiet",n,0,0,0,quiet,quiet.value);
      for (unsigned step=0;step<5;++step) s.AdvanceExplicitEuler(f.geometry.get(),f.config.get(),step*.001,.001);
      const auto rFinal=s.ComputeSpatialResidual(f.geometry.get(),f.config.get(),.005);
      const auto final=MeasureGlobal(s,rFinal);
      for (unsigned k=0;k<5;++k) CHECK(final.value[k]==Approx(quiet.value[k]).margin(2e-13));
      PrintGlobal("quiet_final",n,5,.005,.001,final,quiet.value);
    }
    // Direct API diagnostic at t=0, independent of the quiet evolution clock.
    CPrimitiveEulerSolver nonuniform(f.geometry.get(),f.config.get(),0);
    for (unsigned long i=0;i<grid.xy.size();++i) nonuniform.GetPrimitiveNodes().SetPrimitiveState(i,BoxInitial(grid.xy[i]));
    const auto r=nonuniform.ComputeSpatialResidual(f.geometry.get(),f.config.get());
    const auto b=MeasureGlobal(nonuniform,r);
    // Analytic pressure at the same quadrature points distinguishes force
    // quadrature from reconstruction and from the semidiscrete momentum defect.
    Direction analyticForce{};
    for (const auto& face:nonuniform.GetGeometryMetrics().GetBoundaryFaces())
      for (unsigned d=0;d<2;++d) analyticForce[d]+=BoxInitial(face.point)[3]*face.normal[d];
    CHECK(analyticForce[0]==Precise(-.08)); CHECK(analyticForce[1]==Precise(-.06));
    PrintGlobal("initial",n,0,0,0,b,b.value);
  }
}

TEST_CASE("Closed square separates semidiscrete balance from nonlinear Forward Euler drift", "[primitive-global-balance]") {
  for (unsigned n:{8u,16u,32u}) for (unsigned divisor:{1u,2u}) {
    if (divisor==2 && n!=16) continue;
    BalanceBox grid(n);
    MeshFixture f(grid.xy,grid.cells,"PRIMITIVE_EULER",balanceOptions,"MARKER_EULER",grid.markers);
    CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0);
    for (unsigned long i=0;i<grid.xy.size();++i) s.GetPrimitiveNodes().SetPrimitiveState(i,BoxInitial(grid.xy[i]));
    const unsigned steps=divisor*static_cast<unsigned>(std::ceil(.05/s.EstimateTimeStep(f.geometry.get(),f.config.get(),0,.1)));
    const double dt=.05/steps;
    GlobalFive initial{},integrated{}; Direction impulse{};
    for (unsigned step=0;step<=steps;++step) {
      const auto r=s.ComputeSpatialResidual(f.geometry.get(),f.config.get(),step*dt);
      const auto b=MeasureGlobal(s,r);
      if (step==0) initial=b.value;
      if (step==0 || step==steps/2 || step==steps)
        PrintGlobal(divisor==1?"evolve":"half_dt",n,step,step*dt,dt,b,initial,integrated,impulse);
      if (step==steps) break;
      // Independent exact FE polynomial remainders for momentum and energy.
      GlobalFive predicted{};
      for (unsigned long i=0;i<grid.xy.size();++i) {
        const auto w=s.GetPrimitiveNodes().GetPrimitiveState(i); const auto& d=r.rhs[i];
        const double volume=s.GetGeometryMetrics().GetCell(i).volume;
        predicted[1]+=volume*dt*dt*d[0]*d[1]; predicted[2]+=volume*dt*dt*d[0]*d[2];
        predicted[3]+=volume*(dt*dt*(.5*w[0]*(d[1]*d[1]+d[2]*d[2])+d[0]*(w[1]*d[1]+w[2]*d[2]))+
                              .5*dt*dt*dt*d[0]*(d[1]*d[1]+d[2]*d[2]));
      }
      s.AdvanceExplicitEuler(f.geometry.get(),f.config.get(),step*dt,dt);
      GlobalFive after{};
      for (unsigned long i=0;i<grid.xy.size();++i) {
        const auto q=GlobalDensity(s.GetPrimitiveNodes().GetPrimitiveState(i));
        for (unsigned k=0;k<5;++k) after[k]+=s.GetGeometryMetrics().GetCell(i).volume*q[k];
      }
      for (unsigned k=0;k<4;++k) CHECK(after[k]-b.value[k]-dt*b.rate[k]==Approx(predicted[k]).margin(2e-12));
      for (unsigned k=0;k<5;++k) integrated[k]+=dt*b.rate[k];
      for (unsigned d=0;d<2;++d) impulse[d]+=dt*b.force[d];
    }
  }
}

namespace {
// Phase 11: uniformly refine a hexagonal fan; only new boundary vertices
// are projected onto the circle. No manufactured forcing or boundary callback.
struct RotatingDisk {
  std::vector<Direction> xy{{{0,0}}};
  std::vector<Triangle> cells;
  std::set<std::pair<unsigned long,unsigned long>> boundary;
  explicit RotatingDisk(unsigned levels) {
    const double pi=std::acos(-1.);
    for (unsigned k=0;k<6;++k) xy.push_back({{std::cos(k*pi/3),std::sin(k*pi/3)}});
    for (unsigned k=0;k<6;++k) {
      const unsigned long a=k+1,b=(k+1)%6+1;
      cells.push_back({{0,a,b}}); boundary.insert(std::minmax(a,b));
    }
    for (unsigned level=0;level<levels;++level) {
      std::map<std::pair<unsigned long,unsigned long>,unsigned long> midpoints;
      std::set<std::pair<unsigned long,unsigned long>> nextBoundary;
      const auto midpoint=[&](unsigned long a,unsigned long b) {
        const auto key=std::make_pair(std::min(a,b),std::max(a,b));
        const auto old=midpoints.find(key);
        if (old!=midpoints.end()) return old->second;
        Direction q{{.5*(xy[a][0]+xy[b][0]),.5*(xy[a][1]+xy[b][1])}};
        const auto node=static_cast<unsigned long>(xy.size());
        if (boundary.count(key)) {
          const double r=std::hypot(q[0],q[1]); for (auto& v:q) v/=r;
          nextBoundary.insert(std::minmax(a,node)); nextBoundary.insert(std::minmax(node,b));
        }
        xy.push_back(q); midpoints.emplace(key,node); return node;
      };
      std::vector<Triangle> next;
      for (const auto& t:cells) {
        const auto a=midpoint(t[0],t[1]),b=midpoint(t[1],t[2]),c=midpoint(t[2],t[0]);
        next.push_back({{t[0],a,c}}); next.push_back({{a,t[1],b}});
        next.push_back({{c,b,t[2]}}); next.push_back({{a,b,c}});
      }
      cells=std::move(next); boundary=std::move(nextBoundary);
    }
    for (const auto& t:cells) {
      const auto& a=xy[t[0]]; const auto& b=xy[t[1]]; const auto& c=xy[t[2]];
      REQUIRE((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])>0);
    }
  }
  double h() const {
    double value=0;
    for (const auto& t:cells) for (unsigned k=0;k<3;++k)
      value=std::max(value,std::hypot(xy[t[k]][0]-xy[t[(k+1)%3]][0],xy[t[k]][1]-xy[t[(k+1)%3]][1]));
    return value;
  }
};
State RotatingExact(const Direction& x) { return {{1,-.2*x[1],.2*x[0],1+.02*(x[0]*x[0]+x[1]*x[1])}}; }
const std::string diskOptions="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
    "MUSCL_FLOW=YES\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\nSLOPE_LIMITER_FLOW=NONE\n"
    "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\n";
void InitializeDisk(CPrimitiveEulerSolver& s,const RotatingDisk& grid) {
  for (unsigned long i=0;i<grid.xy.size();++i) {
    const auto& x=grid.xy[i];
    const Gradient g{{{{0,0}},{{0,-.2}},{{.2,0}},{{.04*x[0],.04*x[1]}}}};
    for (auto value:VolumeProduct(RotatingExact(x),g,1.4)) CHECK(value==Approx(0).margin(2e-17));
    s.GetPrimitiveNodes().SetPrimitiveState(i,RotatingExact(x));
  }
}
struct DiskMeasure {
  std::array<State,3> rhs{},error{}; // all, wall, interior: volume-normalized RMS.
  double slip=0,trace=0,minRho=1e300,minP=1e300,minTangent=1e300;
  double wallVelocityGradientError=0,wallPressureNormalError=0;
  double total() const { double r=0; for (auto v:rhs[0]) r+=v*v; return std::sqrt(r); }
};
DiskMeasure MeasureDisk(CPrimitiveEulerSolver& s,const RotatingDisk& grid,
                        const CPrimitiveEulerSolver::SpatialResult& result) {
  DiskMeasure m;
  std::vector<Direction> normal(grid.xy.size());
  for (const auto& f:s.GetGeometryMetrics().GetBoundaryFaces()) {
    for (unsigned d=0;d<2;++d) normal[f.node][d]+=f.normal[d];
    const auto trace=s.ReconstructAt(f.node,f.point);
    m.trace=std::max(m.trace,std::abs(trace[1]*f.normal[0]+trace[2]*f.normal[1])/std::hypot(f.normal[0],f.normal[1]));
  }
  std::array<double,3> weight{};
  for (unsigned long i=0;i<grid.xy.size();++i) {
    const auto w=s.GetPrimitiveNodes().GetPrimitiveState(i),exact=RotatingExact(grid.xy[i]);
    const double area=s.GetGeometryMetrics().GetCell(i).volume, length=std::hypot(normal[i][0],normal[i][1]);
    const unsigned region=length>0?1:2;
    for (unsigned group:{0u,region}) {
      weight[group]+=area;
      for (unsigned k=0;k<4;++k) {
        REQUIRE(std::isfinite(result.rhs[i][k]));
        m.rhs[group][k]+=area*result.rhs[i][k]*result.rhs[i][k];
        m.error[group][k]+=area*(w[k]-exact[k])*(w[k]-exact[k]);
      }
    }
    REQUIRE(w[0]>0); REQUIRE(w[3]>0);
    m.minRho=std::min(m.minRho,w[0]); m.minP=std::min(m.minP,w[3]);
    if (length>0) {
      m.slip=std::max(m.slip,std::abs(w[1]*normal[i][0]+w[2]*normal[i][1])/length);
      m.minTangent=std::min(m.minTangent,std::abs(-w[1]*normal[i][1]+w[2]*normal[i][0])/length);
      const auto& gradient=s.GetLimitedGradient(i);
      const double gradientError=gradient[1][0]*gradient[1][0]+std::pow(gradient[1][1]+.2,2)+
          std::pow(gradient[2][0]-.2,2)+gradient[2][1]*gradient[2][1];
      m.wallVelocityGradientError=std::max(m.wallVelocityGradientError,std::sqrt(gradientError));
      const double normalError=((gradient[3][0]-.04*grid.xy[i][0])*normal[i][0]+
                                (gradient[3][1]-.04*grid.xy[i][1])*normal[i][1])/length;
      m.wallPressureNormalError=std::max(m.wallPressureNormalError,std::abs(normalError));
    }
  }
  CHECK(m.slip<2e-13);
  for (unsigned group=0;group<3;++group) for (unsigned k=0;k<4;++k) {
    m.rhs[group][k]=std::sqrt(m.rhs[group][k]/weight[group]);
    m.error[group][k]=std::sqrt(m.error[group][k]/weight[group]);
  }
  return m;
}
void PrintDisk(const char* kind,const RotatingDisk& grid,unsigned step,double time,double dt,const DiskMeasure& m) {
  const std::array<const char*,3> regions{{"all","wall","interior"}};
  for (unsigned group=0;group<3;++group) {
    std::cout << std::setprecision(16) << "CURVED_WALL," << kind << ',' << grid.boundary.size() << ','
              << grid.xy.size() << ',' << grid.h() << ',' << step << ',' << time << ',' << dt << ',' << regions[group];
    for (auto value:m.rhs[group]) std::cout << ',' << value;
    for (auto value:m.error[group]) std::cout << ',' << value;
    std::cout << ',' << m.total() << ',' << m.slip << ',' << m.trace << ',' << m.minRho << ',' << m.minP << ',' << m.minTangent
              << ',' << m.wallVelocityGradientError << ',' << m.wallPressureNormalError << '\n';
  }
}
}

TEST_CASE("Rotating disk measures curved-wall spatial consistency without an imposed order", "[primitive-curved-wall]") {
  std::cout << "CURVED_WALL,kind,segments,nodes,h,step,time,dt,region,rhs_rho,rhs_u,rhs_v,rhs_p,error_rho,error_u,error_v,error_p,rhs_total,slip,trace_normal,min_rho,min_p,min_tangent,wall_grad_velocity_error,wall_grad_p_normal_error\n";
  for (unsigned level:{1u,2u,3u}) {
    RotatingDisk grid(level);
    MeshFixture f(grid.xy,grid.cells,"PRIMITIVE_EULER",diskOptions,"MARKER_EULER");
    CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0); InitializeDisk(s,grid);
    const auto result=s.ComputeSpatialResidual(f.geometry.get(),f.config.get());
    const auto m=MeasureDisk(s,grid,result);
    CHECK(m.minTangent==Approx(.2).margin(2e-14));
    PrintDisk("spatial",grid,0,0,0,m);
  }
}

TEST_CASE("Primitive LS and WLS recover wall rotation without flat symmetry corrections", "[primitive-wall-gradients]") {
  for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"}) {
    double previous=1e300;
    for (unsigned level:{1u,2u,3u}) {
      RotatingDisk grid(level);
      auto options=diskOptions;
      options.replace(options.find("WEIGHTED_LEAST_SQUARES"),std::string("WEIGHTED_LEAST_SQUARES").size(),method);
      MeshFixture f(grid.xy,grid.cells,"PRIMITIVE_EULER",options,"MARKER_EULER");
      CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0); InitializeDisk(s,grid);
      const auto result=s.ComputeSpatialResidual(f.geometry.get(),f.config.get());
      const auto m=MeasureDisk(s,grid,result);
      CHECK(m.wallVelocityGradientError<2e-12);
      CHECK(m.wallPressureNormalError<previous);
      previous=m.wallPressureNormalError;
      for (unsigned long v=0;v<f.geometry->GetnVertex(0);++v) {
        const auto i=f.geometry->vertex[0][v]->GetNode();
        const auto& g=s.GetLimitedGradient(i);
        const double normalPressure=g[3][0]*grid.xy[i][0]+g[3][1]*grid.xy[i][1];
        CHECK(normalPressure>.005); // Nonzero, not assumed exact for quadratic p.
      }
      std::cout << std::setprecision(16) << "WALL_GRADIENT," << method << ',' << grid.boundary.size()
                << ',' << m.wallVelocityGradientError << ',' << m.wallPressureNormalError << '\n';
    }
  }
}

TEST_CASE("Euler shear survives while true symmetry and generic default retain correction", "[primitive-wall-gradients]") {
  const std::vector<Direction> xy{{{0,0}},{{2,0}},{{2,1}},{{0,1}},{{.8,.3}}};
  const std::vector<Triangle> cells{{{0,1,4}},{{1,2,4}},{{2,3,4}},{{3,0,4}}};
  const std::vector<BoundaryMarker> markers{{"wall",{{{0,1}}}},{"open",{{{1,2}},{{2,3}},{{3,0}}}}};
  for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES","GREEN_GAUSS"})
   for (const std::string kind:{"MARKER_EULER","MARKER_SYM"}) {
    INFO(method << " " << kind);
    const std::string options="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
        "MUSCL_FLOW=YES\nSLOPE_LIMITER_FLOW=NONE\nNUM_METHOD_GRAD_RECON="+method+"\n"+kind+"=(wall)\nMARKER_CUSTOM=(open)\n";
    MeshFixture f(xy,cells,"PRIMITIVE_EULER",options,"",markers);
    CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0);
    for (unsigned i=0;i<xy.size();++i) s.GetPrimitiveNodes().SetPrimitiveState(i,{{1+.1*xy[i][1],.4+.3*xy[i][1],0,1+.2*xy[i][1]}});
    const CPrimitiveEulerSolver::BoundaryState extrapolate=[](const Direction&,double,const State& w,const Direction&,unsigned short){return w;};
    s.ComputeSpatialResidual(f.geometry.get(),f.config.get(),0,extrapolate);
    auto& field=s.GetPrimitiveNodes(); auto& gradient=field.GetGradient_Reconstruction();
    const auto compute=[&](bool legacyDefault) {
      if (method=="GREEN_GAUSS") {
        if (legacyDefault) computeGradientsGreenGauss(nullptr,MPI_QUANTITIES::SOLUTION_GRAD_REC,PERIODIC_NONE,
            *f.geometry,*f.config,field.GetSolution(),0,4,1,gradient);
        else computeGradientsGreenGauss(nullptr,MPI_QUANTITIES::SOLUTION_GRAD_REC,PERIODIC_NONE,
            *f.geometry,*f.config,field.GetSolution(),0,4,1,gradient,false);
      } else {
        if (legacyDefault) computeGradientsLeastSquares(nullptr,MPI_QUANTITIES::SOLUTION_GRAD_REC,PERIODIC_NONE,
            *f.geometry,*f.config,method=="WEIGHTED_LEAST_SQUARES",field.GetSolution(),0,4,1,gradient,field.GetRmatrix());
        else computeGradientsLeastSquares(nullptr,MPI_QUANTITIES::SOLUTION_GRAD_REC,PERIODIC_NONE,
            *f.geometry,*f.config,method=="WEIGHTED_LEAST_SQUARES",field.GetSolution(),0,4,1,gradient,field.GetRmatrix(),false);
      }
    };
    compute(false);
    for (unsigned i:{0u,1u}) {
      for (unsigned k=0;k<4;++k) for (unsigned d=0;d<2;++d)
        CHECK(s.GetLimitedGradient(i)[k][d]==Precise(SU2_TYPE::GetValue(gradient(i,k,d))));
      if (kind=="MARKER_SYM") CHECK(s.GetLimitedGradient(i)[1][1]==Precise(0));
      else {
        CHECK(s.GetLimitedGradient(i)[1][1]>.05);
        if (method!="GREEN_GAUSS") CHECK(s.GetLimitedGradient(i)[1][1]==Precise(.3));
      }
    }
    // Omit the optional argument, exactly as unchanged native call sites do.
    compute(true);
    for (unsigned i:{0u,1u}) {
      CHECK(gradient(i,1,1)==Precise(0));
      CHECK(gradient(i,3,1)==Precise(0));
    }
   }
}

TEST_CASE("Rotating disk separates short evolution from Forward Euler time sensitivity", "[primitive-curved-wall]") {
  for (unsigned level:{1u,2u}) {
   std::vector<State> coarse;
   for (unsigned divisor:{1u,2u}) {
    RotatingDisk grid(level);
    MeshFixture f(grid.xy,grid.cells,"PRIMITIVE_EULER",diskOptions,"MARKER_EULER");
    CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0); InitializeDisk(s,grid);
    const unsigned steps=divisor*static_cast<unsigned>(std::ceil(.1/s.EstimateTimeStep(f.geometry.get(),f.config.get(),0,.05)));
    const double dt=.1/steps;
    for (unsigned step=0;step<=steps;++step) {
      if (step==0 || step==steps/2 || step==steps) {
        const auto result=s.ComputeSpatialResidual(f.geometry.get(),f.config.get(),step*dt);
        PrintDisk(divisor==1?"evolve":"half_dt",grid,step,step*dt,dt,MeasureDisk(s,grid,result));
      }
      if (step<steps) s.AdvanceExplicitEuler(f.geometry.get(),f.config.get(),step*dt,dt);
    }
    if (divisor==1) {
      for (unsigned long i=0;i<grid.xy.size();++i) coarse.push_back(s.GetPrimitiveNodes().GetPrimitiveState(i));
    } else {
      const auto result=s.ComputeSpatialResidual(f.geometry.get(),f.config.get(),.1);
      auto difference=MeasureDisk(s,grid,result);
      difference.error={};
      std::array<double,3> weight{};
      std::set<unsigned long> wall;
      for (const auto& edge:grid.boundary) { wall.insert(edge.first); wall.insert(edge.second); }
      for (unsigned long i=0;i<grid.xy.size();++i) {
        const auto w=s.GetPrimitiveNodes().GetPrimitiveState(i);
        const double volume=s.GetGeometryMetrics().GetCell(i).volume;
        for (unsigned group:{0u,wall.count(i)?1u:2u}) {
          weight[group]+=volume;
          for (unsigned k=0;k<4;++k) difference.error[group][k]+=volume*(w[k]-coarse[i][k])*(w[k]-coarse[i][k]);
        }
      }
      for (unsigned group=0;group<3;++group) for (unsigned k=0;k<4;++k)
        difference.error[group][k]=std::sqrt(difference.error[group][k]/weight[group]);
      PrintDisk("dt_difference",grid,steps,.1,dt,difference);
    }
   }
  }
}

TEST_CASE("Rotating disk monitors a bounded global explicit pseudo-time experiment", "[primitive-curved-wall]") {
  RotatingDisk grid(1);
  MeshFixture f(grid.xy,grid.cells,"PRIMITIVE_EULER",diskOptions,"MARKER_EULER");
  CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0); InitializeDisk(s,grid);
  const double dt=s.EstimateTimeStep(f.geometry.get(),f.config.get(),0,.05);
  for (unsigned step=0;step<=250;++step) {
    if (step==0 || step==10 || step==50 || step==100 || step==150 || step==200 || step==250) {
      const auto result=s.ComputeSpatialResidual(f.geometry.get(),f.config.get(),step*dt);
      PrintDisk("relax",grid,step,step*dt,dt,MeasureDisk(s,grid,result));
    }
    if (step<250) s.AdvanceExplicitEuler(f.geometry.get(),f.config.get(),step*dt,dt);
  }
}

TEST_CASE("Slip reflection preserves thermodynamics and tangential velocity under rotation", "[primitive-boundary]") {
  for (double angle:{0.,.37,1.8,3.14}) {
    Direction n{{std::cos(angle),std::sin(angle)}}, negative{{-n[0],-n[1]}};
    for (double un:{-.2,0.,.2}) {
      const auto w=Oriented(un,.7,n), reflected=Boundary::Reflect(w,n,gasGamma);
      CHECK(reflected[0]==w[0]); CHECK(reflected[3]==w[3]);
      CHECK(Boundary::NormalVelocity(reflected,n)==Precise(-un).margin(2e-15));
      CHECK(-reflected[1]*n[1]+reflected[2]*n[0]==Precise(.7));
      Same(Boundary::Reflect(reflected,n,gasGamma),w);
      Same(reflected,Boundary::Reflect(w,negative,gasGamma));
      if (un==0) Same(reflected,w);
    }
  }
  CHECK_THROWS_AS(Boundary::Reflect({{1,0,0,1}},{{2,0}},gasGamma),std::invalid_argument);
}

TEST_CASE("Farfield imposes only incoming characteristic data in each normal regime", "[primitive-boundary]") {
  const double c=std::sqrt(gasGamma);
  for (double angle:{0.,.37,1.8}) {
    Direction n{{std::cos(angle),std::sin(angle)}};
    for (double u:{-2*c,-.4*c,0.,.4*c,2*c}) {
      const auto w=Oriented(u,.2,n);
      Same(Boundary::Farfield(w,w,n,gasGamma),w);
      State exterior=w; exterior[0]*=1.0001; exterior[3]*=1.0002;
      const auto b=Boundary::Farfield(w,exterior,n,gasGamma);
      if (u!=0) CHECK(Boundary::Incoming(b,n,gasGamma)==Boundary::Incoming(w,n,gasGamma));
      if (u<=-c) Same(b,exterior);
      else if (u>=c) Same(b,w);
      else {
        const double cs=Boundary::Sound(w,gasGamma), deltaU=Boundary::NormalVelocity(b,n)-u;
        CHECK(b[3]-w[3]+w[0]*cs*deltaU==Precise(0).margin(2e-13));
        CHECK(b[3]-w[3]-w[0]*cs*deltaU==Precise(exterior[3]-w[3]-w[0]*cs*(Boundary::NormalVelocity(exterior,n)-u)).margin(2e-13));
        const auto& advected=u<0?exterior:w;
        CHECK(b[0]-w[0]-(b[3]-w[3])/(cs*cs)==Precise(advected[0]-w[0]-(advected[3]-w[3])/(cs*cs)).margin(2e-13));
        CHECK(-b[1]*n[1]+b[2]*n[0]==Precise(-advected[1]*n[1]+advected[2]*n[0]));
      }
    }
    // Linear acoustic eigenvectors: an outgoing wave should not create an
    // incoming characteristic to first order; this is not a zero-reflection theorem.
    for (double amplitude:{1e-3,5e-4,2.5e-4}) {
      const auto base=Oriented(.3,.2,n);
      auto w=base;
      w[0]+=amplitude; w[3]+=gasGamma*amplitude;
      w[1]+=c*amplitude*n[0]; w[2]+=c*amplitude*n[1];
      const auto b=Boundary::Farfield(w,base,n,gasGamma);
      const double cs=Boundary::Sound(w,gasGamma), du=Boundary::NormalVelocity(b,n)-Boundary::NormalVelocity(w,n);
      CHECK(std::abs(b[3]-w[3]+w[0]*cs*du)<2e-13);
      CHECK(std::abs(b[3]-w[3])<10*amplitude*amplitude);
    }
  }
  const State tangent{{1,.3,0,1}}, entropy{{1.2,.3,0,1}};
  Same(Boundary::Farfield(entropy,tangent,{{0,1}},gasGamma),entropy);
  MeshFixture fixture({{{0,0}},{{1,0}},{{0,1}}},{{{0,1,2}}},"PRIMITIVE_EULER",
      "TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\nPRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC\n","MARKER_FAR");
  CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get()),std::invalid_argument);
  CHECK_NOTHROW(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
      [](const Direction&,double,const State& w,const Direction&,unsigned short){return w;}));
}

TEST_CASE("Total and mass inlet states retain the outgoing invariant and physical data", "[primitive-boundary]") {
  for (double angle:{0.,.37,1.8}) {
    Direction n{{std::cos(angle),std::sin(angle)}};
    const auto w=Oriented(-.4,.15,n);
    const double speed=std::hypot(w[1],w[2]);
    const Direction dir{{w[1]/speed,w[2]/speed}};
    const double temperature=w[3]/w[0], tt=temperature+.4*speed*speed/(2*gasGamma);
    const double pt=w[3]*std::pow(tt/temperature,gasGamma/.4);
    Same(Boundary::TotalInlet(w,tt,pt,dir,n,gasGamma,1),w);
    Same(Boundary::MassInlet(w,w[0],speed,dir,n,gasGamma),w);
    auto perturbed=w; perturbed[3]*=1.001;
    const auto total=Boundary::TotalInlet(perturbed,tt,pt,dir,n,gasGamma,1);
    CHECK(Plus(total,n)==Precise(Plus(perturbed,n)).margin(2e-13));
    const double t=total[3]/total[0], v2=total[1]*total[1]+total[2]*total[2];
    CHECK(t+.4*v2/(2*gasGamma)==Precise(tt));
    CHECK(total[3]*std::pow(tt/t,gasGamma/.4)==Precise(pt));
    const auto mass=Boundary::MassInlet(perturbed,1,speed,dir,n,gasGamma);
    CHECK(mass[0]==1); CHECK(std::hypot(mass[1],mass[2])==Precise(speed));
    CHECK(Plus(mass,n)==Precise(Plus(perturbed,n)));
    CHECK_THROWS_AS(Boundary::TotalInlet(w,-1,pt,dir,n,gasGamma,1),std::invalid_argument);
  }
}

TEST_CASE("Characteristic reservoir separates entropy vorticity acoustic and sonic limits", "[primitive-boundary]") {
  for (double angle:{0.,.37,1.8}) {
    const Direction n{{std::cos(angle),std::sin(angle)}};
    for (double un:{-.3,.3}) {
      const auto base=Oriented(un,.2,n);
      for (unsigned mode=0;mode<2;++mode) {
        auto w=base;
        if (mode==0) w[0]+=.01;
        else { w[1]-=.01*n[1];w[2]+=.01*n[0]; }
        Same(Boundary::Farfield(w,base,n,gasGamma),un<0?base:w);
      }
    }
    for (double factor:{-1.000001,-1.,-.999999,.999999,1.,1.000001}) {
      const auto w=Oriented(factor*std::sqrt(gasGamma),.2,n);
      Same(Boundary::Farfield(w,w,n,gasGamma),w);
      const unsigned expected=factor < -1?4:(factor < 0?3:(factor<1?1:0));
      CHECK(Boundary::Incoming(w,n,gasGamma)==expected);
    }
    for (double a:{1e-3,5e-4,2.5e-4}) {
      const auto base=Oriented(.3,.2,n);auto incoming=base;
      incoming[0]+=a;incoming[3]+=gasGamma*a;
      incoming[1]-=std::sqrt(gasGamma)*a*n[0];incoming[2]-=std::sqrt(gasGamma)*a*n[1];
      const auto b=Boundary::Farfield(incoming,base,n,gasGamma);
      for (unsigned k=0;k<4;++k) CHECK(std::abs(b[k]-base[k])<5*a*a);
    }
  }
}

TEST_CASE("Pressure outlet extrapolates entropy tangent and outgoing invariant with explicit backflow rejection", "[primitive-boundary]") {
  for (double angle:{0.,.37,1.8}) {
    Direction n{{std::cos(angle),std::sin(angle)}};
    const auto w=Oriented(.4,.2,n), b=Boundary::Outlet(w,1.001,n,gasGamma);
    Same(Boundary::Outlet(w,1,n,gasGamma),w);
    CHECK(b[3]==1.001); CHECK(Plus(b,n)==Precise(Plus(w,n)));
    CHECK(b[3]/std::pow(b[0],gasGamma)==Precise(w[3]/std::pow(w[0],gasGamma)));
    CHECK(-b[1]*n[1]+b[2]*n[0]==Precise(.2));
    Same(Boundary::Outlet(Oriented(2,.2,n),2,n,gasGamma),Oriented(2,.2,n));
    CHECK_THROWS_AS(Boundary::Outlet(Oriented(-.1,.2,n),1,n,gasGamma),std::invalid_argument);
    CHECK_THROWS_AS(Boundary::SupersonicOutlet(w,n,gasGamma),std::invalid_argument);
    Same(Boundary::SupersonicInlet(Oriented(-2,.2,n),Oriented(-2.1,.3,n),n,gasGamma),Oriented(-2.1,.3,n));
  }
}

TEST_CASE("Nodal slip reaction preserves impermeability without changing interior contributions", "[primitive-boundary]") {
  const std::vector<Triangle> triangles{{{0,1,4}},{{1,2,4}},{{2,3,4}},{{3,0,4}}};
  for (bool jitter:{false,true}) for (double angle:{0.,.37}) for (const std::string marker:{"MARKER_EULER","MARKER_SYM"}) {
    std::vector<Direction> xy{{{0,0}},{{1,0}},{{1,1}},{{0,1}},{{jitter?.42:.5,jitter?.57:.5}}};
    for (auto& x:xy) x={{std::cos(angle)*x[0]-std::sin(angle)*x[1],std::sin(angle)*x[0]+std::cos(angle)*x[1]}};
    const std::string options="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
        "NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\nSLOPE_LIMITER_FLOW=NONE\n";
    const std::vector<BoundaryMarker> markers{{"south",{{{0,1}}}},{"east",{{{1,2}}}},
        {"north",{{{2,3}}}},{"west",{{{3,0}}}}};
    MeshFixture f(xy,triangles,"PRIMITIVE_EULER",options,marker,markers);
    CPrimitiveEulerSolver solver(f.geometry.get(),f.config.get(),0);
    for (unsigned i=0;i<5;++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,{{1,0,0,1+.1*xy[i][0]}});
    // This nonconstant pressure demonstrates why reflection alone is insufficient.
    const auto state=solver.GetPrimitiveNodes().GetPrimitiveState(4);
    const auto result=solver.ComputeSpatialResidual(f.geometry.get(),f.config.get());
    for (unsigned i=0;i<4;++i) {
      CHECK(std::abs(result.rhs[i][1])<1e-13); CHECK(std::abs(result.rhs[i][2])<1e-13);
      for (unsigned k=0;k<4;++k) CHECK(result.rhs[i][k]==Precise(-(result.faces[i][k]+result.volume[i][k]+result.boundary[i][k])/solver.GetGeometryMetrics().GetCell(i).volume).margin(1e-13));
    }
    CHECK(std::abs(result.rhs[4][1])>0.01); // No interior constraint.
    solver.AdvanceExplicitEuler(f.geometry.get(),f.config.get(),0,1e-3);
    for (unsigned i=0;i<4;++i) {
      const auto w=solver.GetPrimitiveNodes().GetPrimitiveState(i);
      CHECK(std::abs(w[1])<1e-15); CHECK(std::abs(w[2])<1e-15);
    }
    CHECK(solver.GetPrimitiveNodes().GetPrimitiveState(4)!=state);
  }
}

TEST_CASE("Rotated channel slip preserves tangential acceleration and mixed-marker corners", "[primitive-boundary]") {
  for (double angle:{0.,.37}) for (bool reverse:{false,true}) {
    const Direction tangent{{std::cos(angle),std::sin(angle)}},normal{{-tangent[1],tangent[0]}};
    std::vector<Direction> xy{{{0,0}},{{2,0}},{{2,1}},{{0,1}},{{.8,.3}}};
    for (auto& x:xy) x={{tangent[0]*x[0]+normal[0]*x[1],tangent[1]*x[0]+normal[1]*x[1]}};
    const std::vector<Triangle> triangles{{{0,1,4}},{{1,2,4}},{{2,3,4}},{{3,0,4}}};
    const std::vector<BoundaryMarker> markers{{"walls",{{{0,1}},{{2,3}}}},{"open",{{{1,2}},{{3,0}}}}};
    MeshFixture f(xy,triangles,"PRIMITIVE_EULER","TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
        "NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\nSLOPE_LIMITER_FLOW=NONE\n"
        "MARKER_EULER=(walls)\nMARKER_CUSTOM=(open)\nPRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\n","",markers,reverse);
    CPrimitiveEulerSolver solver(f.geometry.get(),f.config.get(),0);
    for (unsigned i=0;i<5;++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,
        {{1,.4*tangent[0],.4*tangent[1],1+.1*(xy[i][0]*tangent[0]+xy[i][1]*tangent[1])}});
    const CPrimitiveEulerSolver::BoundaryState extrapolate=[](const Direction&,double,const State& w,const Direction&,unsigned short){return w;};
    const auto r=solver.AdvanceExplicitEuler(f.geometry.get(),f.config.get(),0,.001,extrapolate);
    for (unsigned i=0;i<5;++i) {
      CHECK(r.rhs[i][1]*tangent[0]+r.rhs[i][2]*tangent[1]==Precise(-.1).margin(2e-13));
      CHECK(r.rhs[i][1]*normal[0]+r.rhs[i][2]*normal[1]==Precise(0).margin(2e-13));
      const auto w=solver.GetPrimitiveNodes().GetPrimitiveState(i);
      CHECK(w[1]*tangent[0]+w[2]*tangent[1]==Precise(.3999).margin(2e-13));
      CHECK(w[1]*normal[0]+w[2]*normal[1]==Precise(0).margin(2e-13));
    }
  }
}

TEST_CASE("Nonplanar symmetry remains explicitly unsupported", "[primitive-boundary]") {
  const std::vector<Direction> xy{{{0,0}},{{1,0}},{{0,1}}};
  MeshFixture f(xy,{{{0,1,2}}},"PRIMITIVE_EULER","TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n","MARKER_SYM");
  CPrimitiveEulerSolver solver(f.geometry.get(),f.config.get(),0);
  for (unsigned i=0;i<3;++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,{{1,0,0,1}});
  try { solver.ComputeSpatialResidual(f.geometry.get(),f.config.get());FAIL("Curved marker accepted"); }
  catch (const std::invalid_argument& e) { CHECK(std::string(e.what()).find("symmetry marker")!=std::string::npos); }
}

TEST_CASE("Nonmanifold Euler wall is rejected by geometry before initialization", "[primitive-boundary][primitive-polygon]") {
  // Two fluid regions touch only at vertex zero: four physical half-faces,
  // not a single smooth polygonal wall neighbourhood.
  const std::vector<Direction> xy{{{0,0}},{{1,0}},{{0,1}},{{-1,0}},{{0,-1}}};
  MeshFixture f(xy,{{{0,1,2}},{{0,3,4}}},"PRIMITIVE_EULER",
      "TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\nPRIMITIVE_INITIALIZATION=FREESTREAM\n","MARKER_EULER");
  // The existing geometry audit already rejects this topology: do not bypass
  // that earlier guard merely to exercise the later wall-incidence guard.
  try { CPrimitiveEulerSolver solver(f.geometry.get(),f.config.get(),0); FAIL("Nonmanifold wall accepted"); }
  catch (const std::invalid_argument& e) {
    CHECK(std::string(e.what()).find("Primitive boundary normal inconsistent with SU2")!=std::string::npos);
  }
}

TEST_CASE("Polygonal Euler wall uses one resultant normal and explicit freestream projection", "[primitive-boundary][primitive-polygon]") {
  const std::vector<Direction> xy{{{0,0}},{{1,.2}},{{2,.3}},{{2,1.3}},{{0,1}},{{.9,.7}}};
  const std::vector<Triangle> cells{{{0,1,5}},{{1,2,5}},{{2,3,5}},{{3,4,5}},{{4,0,5}}};
  const std::vector<BoundaryMarker> markers{{"wall",{{{0,1}},{{1,2}}}},
                                           {"open",{{{2,3}},{{3,4}},{{4,0}}}}};
  for (bool reverse:{false,true}) {
    MeshFixture f(xy,cells,"PRIMITIVE_EULER","TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
        "PRIMITIVE_INITIALIZATION=FREESTREAM\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\n"
        "SLOPE_LIMITER_FLOW=NONE\nMARKER_EULER=(wall)\nMARKER_CUSTOM=(open)\n","",markers,reverse);
    CPrimitiveEulerSolver solver(f.geometry.get(),f.config.get(),0);
    const State free{{1,.4,.2,1}};
    for (unsigned i=0;i<xy.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,free);
    const CPrimitiveEulerSolver::BoundaryState extrapolate=[](const Direction&,double,const State& w,const Direction&,unsigned short){return w;};
    CHECK_THROWS_AS(solver.ComputeSpatialResidual(f.geometry.get(),f.config.get(),0,extrapolate),std::invalid_argument);
    std::vector<Direction> normal(xy.size());
    for (const auto& face:solver.GetGeometryMetrics().GetBoundaryFaces()) if (face.marker==0)
      for (unsigned d=0;d<2;++d) normal[face.node][d]+=face.normal[d];
    CHECK(normal[1][0]==Precise(.15)); CHECK(normal[1][1]==Precise(-1));
    for (unsigned long v=0;v<f.geometry->GetnVertex(0);++v) {
      const auto* vertex=f.geometry->vertex[0][v];
      for (unsigned d=0;d<2;++d) CHECK(normal[vertex->GetNode()][d]==Precise(-vertex->GetNormal(d)));
    }
    for (auto& n:normal) { const double length=std::hypot(n[0],n[1]); if (length>0) for (auto& a:n) a/=length; }
    REQUIRE(solver.ProjectFreestreamWallVelocity()==3);
    for (unsigned i=0;i<xy.size();++i) {
      const auto w=solver.GetPrimitiveNodes().GetPrimitiveState(i);
      CHECK(w[0]==free[0]); CHECK(w[3]==free[3]);
      if (i>2) Same(w,free);
      else {
        const double un=free[1]*normal[i][0]+free[2]*normal[i][1];
        CHECK(w[1]==Precise(free[1]-un*normal[i][0]));
        CHECK(w[2]==Precise(free[2]-un*normal[i][1]));
        CHECK(std::hypot(w[1],w[2])>.3);
      }
    }
    for (unsigned step=0;step<5;++step) {
      const auto result=solver.AdvanceExplicitEuler(f.geometry.get(),f.config.get(),step*1e-4,1e-4,extrapolate);
      for (unsigned i=0;i<xy.size();++i) {
        const auto w=solver.GetPrimitiveNodes().GetPrimitiveState(i);
        CHECK(w[0]>0); CHECK(w[3]>0);
        for (auto value:result.rhs[i]) CHECK(std::isfinite(value));
        if (i<3) CHECK(w[1]*normal[i][0]+w[2]*normal[i][1]==Precise(0));
      }
    }
    CHECK_THROWS_AS(solver.ProjectFreestreamWallVelocity(),std::invalid_argument);
  }
}
