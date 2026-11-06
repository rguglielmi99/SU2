/*!
 * \brief Primitive Forward Euler: physical stepping, CFL scale and rejection guarantees.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include <limits>

using namespace PrimitiveEuler2D;
using namespace PrimitiveEulerTest;
namespace {
const std::vector<Direction> points{{{0,0}},{{2,0}},{{2,1}},{{0,1}},{{0.8,0.3}}};
const std::vector<Triangle> elements{{{0,1,4}},{{1,2,4}},{{2,3,4}},{{3,0,4}}};
const std::string options="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
    "NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\nPRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\n"
    "SLOPE_LIMITER_FLOW=VENKATAKRISHNAN\nCFL_NUMBER=0.1\n";
const CPrimitiveEulerSolver::BoundaryState extrapolate=
    [](const Direction&,double,const State& w,const Direction&,unsigned short){return w;};
void Check(const State& a,const State& b,double tolerance=2e-12) {
  for (unsigned k=0;k<4;++k) CHECK(a[k]==Approx(b[k]).margin(tolerance));
}
}

TEST_CASE("Local pseudo-time is synchronous stationary and consistent with equal-step Forward Euler", "[primitive-pseudo-time]") {
  MeshFixture fixture(points,elements,"PRIMITIVE_EULER",options);
  CPrimitiveEulerSolver local(fixture.geometry.get(),fixture.config.get(),0),
                        physical(fixture.geometry.get(),fixture.config.get(),0);
  for (unsigned long i=0;i<points.size();++i) local.GetPrimitiveNodes().SetPrimitiveState(i,{{1,0,0,1}});
  local.AdvanceLocalPseudoTime(fixture.geometry.get(),fixture.config.get(),.1,extrapolate);
  CHECK(local.GetPseudoDiagnostics().ratio>1.1);
  for (unsigned long i=0;i<points.size();++i) Check(local.GetPrimitiveNodes().GetPrimitiveState(i),{{1,0,0,1}},0);
  const auto source=[](const Direction& x,double t) { CHECK(t==0); return State{{.1+x[0],-.3,.2,.1+x[1]}}; };
  std::vector<State> before(points.size());
  for (unsigned long i=0;i<points.size();++i) {
    before[i]={{2+.1*points[i][0],.5+.03*points[i][1],.3,3+.2*points[i][1]}};
    local.GetPrimitiveNodes().SetPrimitiveState(i,before[i]); physical.GetPrimitiveNodes().SetPrimitiveState(i,before[i]);
  }
  const auto rhs=local.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,extrapolate,source);
  const double dt=1e-5;
  physical.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),0,dt,extrapolate,source);
  local.AdvanceLocalPseudoTime(fixture.geometry.get(),fixture.config.get(),.1,extrapolate,source,0,std::vector<double>(points.size(),dt));
  CHECK(local.GetPhysicalTime()==0); CHECK(physical.GetPhysicalTime()==dt);
  for (unsigned long i=0;i<points.size();++i) {
    Check(local.GetPrimitiveNodes().GetPrimitiveState(i),physical.GetPrimitiveNodes().GetPrimitiveState(i),0);
    auto expected=before[i]; for (unsigned k=0;k<4;++k) expected[k]+=dt*rhs.rhs[i][k];
    Check(local.GetPrimitiveNodes().GetPrimitiveState(i),expected,0);
    local.GetPrimitiveNodes().SetPrimitiveState(i,before[i]);
  }
  local.AdvanceLocalPseudoTime(fixture.geometry.get(),fixture.config.get(),.01,extrapolate,source);
  const auto& info=local.GetPseudoDiagnostics();
  for (unsigned long i=0;i<points.size();++i) {
    CHECK(info.steps[i]==Approx(.01*local.GetGeometryMetrics().GetCell(i).volume/rhs.spectralRadius[i]).epsilon(1e-14));
    auto expected=before[i]; for (unsigned k=0;k<4;++k) expected[k]+=info.steps[i]*rhs.rhs[i][k];
    Check(local.GetPrimitiveNodes().GetPrimitiveState(i),expected,0);
  }
  const auto fresh=local.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,extrapolate,source);
  for (unsigned long i=0;i<points.size();++i) for (unsigned k=0;k<4;++k)
    CHECK(local.LinSysRes(i,k)==-local.GetGeometryMetrics().GetCell(i).volume*fresh.rhs[i][k]);
  const auto frozen=local.GetPrimitiveNodes().GetSolution();
  const auto tiny=local.AdvanceLocalPseudoTime(fixture.geometry.get(),fixture.config.get(),1e-300,extrapolate,source);
  double residualNorm=0;
  for (unsigned long i=0;i<points.size();++i) for (unsigned k=0;k<4;++k) {
    CHECK(local.GetPrimitiveNodes().GetSolution(i,k)==frozen(i,k)); // Update rounds away.
    residualNorm+=tiny.rhs[i][k]*tiny.rhs[i][k];
  }
  CHECK(residualNorm>.1); // A tiny update must never masquerade as a small spatial residual.
}

TEST_CASE("Local pseudo-time backtracks globally and restores all work arrays on exhausted or exceptional attempts", "[primitive-pseudo-time]") {
  MeshFixture f(points,elements,"PRIMITIVE_EULER",options);
  CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0);
  for (unsigned long i=0;i<points.size();++i) s.GetPrimitiveNodes().SetPrimitiveState(i,{{1,0,0,1}});
  const auto sink=[](const Direction&,double) { return State{{-2,0,0,-2}}; };
  std::vector<double> steps{1,2,1,2,1};
  s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,extrapolate,sink,8,steps);
  const auto& info=s.GetPseudoDiagnostics();
  REQUIRE(info.reductions==3); CHECK(info.effectiveCfl==Approx(.0125)); CHECK(info.rejectedAttempts==3);
  for (unsigned long i=0;i<points.size();++i) {
    CHECK(info.steps[i]==steps[i]/8);
    Check(s.GetPrimitiveNodes().GetPrimitiveState(i),{{1-steps[i]/4,0,0,1-steps[i]/4}},0);
  }
  auto& nodes=s.GetPrimitiveNodes();
  const auto solution=nodes.GetSolution();
  const auto raw=nodes.GetGradient_Reconstruction(),rmatrix=nodes.GetRmatrix();
  const auto limiter=nodes.GetLimiter();
  std::vector<State> residual(points.size()),old(points.size());
  std::vector<Gradient> gradient(points.size()); std::vector<double> theta(points.size());
  for (unsigned long i=0;i<points.size();++i) {
    gradient[i]=s.GetLimitedGradient(i); theta[i]=s.GetPositivityScale(i);
    for (unsigned k=0;k<4;++k) { residual[i][k]=SU2_TYPE::GetValue(s.LinSysRes(i,k)); old[i][k]=SU2_TYPE::GetValue(nodes.GetSolution_Old(i,k)); }
  }
  CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,extrapolate,sink,1,steps),std::domain_error);
  const auto throwing=[](const Direction&,double,const State&,const Direction&,unsigned short)->State {throw std::runtime_error("callback failure");};
  CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,throwing),std::runtime_error);
  CHECK(s.GetPseudoDiagnostics().rejectedIterations==2);
  CHECK(s.GetPseudoDiagnostics().rejectedAttempts==5); CHECK(s.GetPseudoDiagnostics().accepted==1);
  CHECK(s.GetPhysicalTime()==0);
  for (unsigned long i=0;i<points.size();++i) {
    CHECK(s.GetLimitedGradient(i)==gradient[i]); CHECK(s.GetPositivityScale(i)==theta[i]);
    for (unsigned k=0;k<4;++k) {
      CHECK(nodes.GetSolution(i,k)==solution(i,k)); CHECK(nodes.GetSolution_Old(i,k)==old[i][k]);
      CHECK(s.LinSysRes(i,k)==residual[i][k]); CHECK(nodes.GetLimiter(i,k)==limiter(i,k));
      for (unsigned d=0;d<2;++d) CHECK(nodes.GetGradient_Reconstruction()(i,k,d)==raw(i,k,d));
    }
    for (unsigned a=0;a<2;++a) for (unsigned b=0;b<2;++b) CHECK(nodes.GetRmatrix(i,a,b)==rmatrix(i,a,b));
  }
  CHECK_NOTHROW(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.01,extrapolate));
}

TEST_CASE("Local pseudo-time rejects invalid CFL states geometry and step requests", "[primitive-pseudo-time]") {
  MeshFixture f(points,elements,"PRIMITIVE_EULER",options);
  CPrimitiveEulerSolver s(f.geometry.get(),f.config.get(),0);
  for (double bad:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
    CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),bad,extrapolate),std::invalid_argument);
  CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,extrapolate,{},0,{1}),std::invalid_argument);
  CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,extrapolate,{},0,std::vector<double>(points.size(),0)),std::invalid_argument);
  auto& cell=const_cast<CAggregatedGeometry::Cell&>(s.GetGeometryMetrics().GetCell(0));
  const double volume=cell.volume;
  for (double invalid:{-1.,0.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    cell.volume=invalid;
    CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,extrapolate),std::invalid_argument);
  }
  cell.volume=volume;
  s.GetPrimitiveNodes().GetSolution(0)[3]=std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,extrapolate),std::invalid_argument);
  s.GetPrimitiveNodes().SetPrimitiveState(0,{{1,.5,.3,1}});
  const auto bad=[](const Direction&,double,const State&,const Direction&,unsigned short) {return State{{1,0,0,std::numeric_limits<double>::infinity()}};};
  CHECK_THROWS_AS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,bad),std::invalid_argument);
  for (unsigned long i=0;i<points.size();++i) s.GetPrimitiveNodes().SetPrimitiveState(i,{{1e-308,0,0,1e308}});
  CHECK_THROWS(s.AdvanceLocalPseudoTime(f.geometry.get(),f.config.get(),.1,extrapolate)); // Spectral overflow, despite finite positive inputs.
}

TEST_CASE("Forward Euler applies one rhs with correct sign volume coordinates and stage time", "[primitive-solver][primitive-time]") {
  MeshFixture fixture(points,elements,"PRIMITIVE_EULER",options);
  CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
  std::vector<State> before;
  for (unsigned long i=0;i<points.size();++i) {
    State w{{2+0.1*points[i][0],0.5+0.03*points[i][1],0.3,3+0.2*points[i][1]}};
    solver.GetPrimitiveNodes().SetPrimitiveState(i,w);before.push_back(w);
  }
  const double time=0.7,dt=1e-5;
  unsigned sourceCalls=0,boundaryCalls=0;
  const auto boundary=[&](const Direction&,double t,const State& w,const Direction&,unsigned short) {
    ++boundaryCalls;CHECK(t==time);return w;
  };
  const auto source=[&](const Direction& x,double t) {
    ++sourceCalls;CHECK(t==time);return State{{0.1+x[0]*t,-0.2,0.3,0.4+x[1]}};
  };
  const auto reference=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),time,boundary,source);
  sourceCalls=boundaryCalls=0;
  // A stale previously assembled LinSysRes must never drive a time step.
  for (unsigned long i=0;i<points.size();++i) for (unsigned k=0;k<4;++k) solver.LinSysRes(i,k)=1e4;
  const auto stage=solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),time,dt,boundary,source);
  CHECK(sourceCalls==points.size());CHECK(boundaryCalls==solver.GetGeometryMetrics().GetBoundaryFaces().size());
  CHECK(solver.GetPhysicalTime()==time+dt);
  for (unsigned long i=0;i<points.size();++i) {
    Check(stage.rhs[i],reference.rhs[i]);State expected=before[i];
    for (unsigned k=0;k<4;++k) {
      expected[k]+=dt*reference.rhs[i][k];
      CHECK(solver.GetPrimitiveNodes().GetSolution_Old(i,k)==before[i][k]);
      CHECK(solver.LinSysRes(i,k)==-solver.GetGeometryMetrics().GetCell(i).volume*stage.rhs[i][k]);
    }
    Check(solver.GetPrimitiveNodes().GetPrimitiveState(i),expected);
  }
  // A second stage recomputes gradients/limiter/residual at the updated W.
  const auto fresh=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),time+dt,extrapolate);
  const auto next=solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),time+dt,dt,extrapolate);
  for (unsigned long i=0;i<points.size();++i) Check(next.rhs[i],fresh.rhs[i]);
}

TEST_CASE("CFL estimate uses integrated face speeds global minimum geometry scaling and exterior data", "[primitive-solver][primitive-time]") {
  double initialDt=0;
  for (double scale:{1.0,0.5,0.25}) {
    auto xy=points;for (auto& x:xy) for (auto& value:x) value*=scale;
    MeshFixture fixture(xy,elements,"PRIMITIVE_EULER",options);
    CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
    const auto state=solver.GetPrimitiveNodes().GetPrimitiveState(0);
    const auto result=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,extrapolate);
    const auto& metrics=solver.GetGeometryMetrics();std::vector<double> expected(xy.size());
    const auto speed=[&](const Direction& n){return std::abs(state[1]*n[0]+state[2]*n[1])+std::sqrt(1.4*state[3]/state[0])*std::hypot(n[0],n[1]);};
    for (unsigned long e=0;e<metrics.GetnEdge();++e) {
      const auto& f=metrics.GetFace(e);expected[f.first]+=speed(f.normal);expected[f.second]+=speed(f.normal);
    }
    for (const auto& f:metrics.GetBoundaryFaces()) expected[f.node]+=speed(f.normal);
    double dt=1e10;
    for (unsigned long i=0;i<xy.size();++i) {
      CHECK(result.spectralRadius[i]==Approx(expected[i]).epsilon(1e-13));
      dt=std::min(dt,0.1*metrics.GetCell(i).volume/expected[i]);
      for (unsigned k=0;k<4;++k) solver.LinSysRes(i,k)=42+i+k;
    }
    CHECK(solver.EstimateTimeStep(fixture.geometry.get(),fixture.config.get(),0,0.1,extrapolate)==Approx(dt).epsilon(1e-13));
    CHECK(solver.EstimateTimeStep(fixture.geometry.get(),fixture.config.get(),0,0.2,extrapolate)==Approx(2*dt).epsilon(1e-13));
    if (scale==1) initialDt=dt;else CHECK(dt==Approx(initialDt*scale).epsilon(1e-13));
    const auto fast=[](const Direction&,double,const State& w,const Direction&,unsigned short){auto out=w;out[1]=100;return out;};
    CHECK(solver.EstimateTimeStep(fixture.geometry.get(),fixture.config.get(),0,0.1,fast)<dt);
    for (unsigned long i=0;i<xy.size();++i) {
      Check(solver.GetPrimitiveNodes().GetPrimitiveState(i),state);
      for (unsigned k=0;k<4;++k) CHECK(solver.LinSysRes(i,k)==42+i+k);
    }
  }
}

TEST_CASE("Rejected time steps preserve all nodal states residual old state and clock and release guard", "[primitive-solver][primitive-time]") {
  MeshFixture fixture(points,elements,"PRIMITIVE_EULER",options);
  CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
  const auto initial=solver.GetPrimitiveNodes().GetPrimitiveState(0);
  for (unsigned long i=0;i<points.size();++i) for (unsigned k=0;k<4;++k) solver.LinSysRes(i,k)=i+k+7;
  const auto unchanged=[&]() {
    CHECK(solver.GetPhysicalTime()==0);
    for (unsigned long i=0;i<points.size();++i) {
      Check(solver.GetPrimitiveNodes().GetPrimitiveState(i),initial);
      for (unsigned k=0;k<4;++k) {
        CHECK(solver.LinSysRes(i,k)==i+k+7);
        CHECK(solver.GetPrimitiveNodes().GetSolution_Old(i,k)==initial[k]);
      }
    }
  };
  for (double dt:{0.0,-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),0,dt,extrapolate),std::invalid_argument);unchanged();
  }
  CHECK_THROWS_AS(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),1e300,1e-10,extrapolate),std::invalid_argument);unchanged();
  for (double cfl:{0.0,-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS(solver.EstimateTimeStep(fixture.geometry.get(),fixture.config.get(),0,cfl,extrapolate),std::invalid_argument);unchanged();
  }
  CHECK_THROWS_AS(solver.AdvanceExplicitEuler(nullptr,fixture.config.get(),0,0.01,extrapolate),std::invalid_argument);unchanged();
  for (unsigned component:{0u,3u}) {
    CHECK_THROWS_AS(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),0,1,extrapolate,
        [&](const Direction& x,double){State s{};if(x[0]>1) s[component]=-100;return s;}),std::domain_error);unchanged();
  }
  CHECK_THROWS_AS(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),0,0.01,
      [](const Direction&,double,const State&,const Direction&,unsigned short)->State{throw std::runtime_error("boundary");}),std::runtime_error);unchanged();
  CHECK_THROWS_AS(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),0,0.01,extrapolate,
      [](const Direction&,double){return State{{0,0,std::numeric_limits<double>::quiet_NaN(),0}};}),std::invalid_argument);unchanged();
  unsigned calls=0;
  CHECK_NOTHROW(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),0,0.01,extrapolate,
      [&](const Direction&,double) {
        ++calls;
        CHECK_THROWS_AS(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),0,0.01,extrapolate),std::logic_error);
        CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,extrapolate),std::logic_error);
        CHECK_THROWS_AS(solver.EstimateTimeStep(fixture.geometry.get(),fixture.config.get(),0,0.1,extrapolate),std::logic_error);
        return State{{0.2,0.1,-0.05,0.3}};
      }));
  CHECK(calls==points.size());CHECK(solver.GetPhysicalTime()==0.01);
  const double acceptedTime=solver.GetPhysicalTime();
  std::vector<State> acceptedState(points.size()),acceptedOld(points.size()),acceptedResidual(points.size());
  double maxSavedResidual=0,maxOldDifference=0;
  for (unsigned long i=0;i<points.size();++i) {
    acceptedState[i]=solver.GetPrimitiveNodes().GetPrimitiveState(i);
    for (unsigned k=0;k<4;++k) {
      acceptedOld[i][k]=SU2_TYPE::GetValue(solver.GetPrimitiveNodes().GetSolution_Old(i,k));
      acceptedResidual[i][k]=SU2_TYPE::GetValue(solver.LinSysRes(i,k));
      maxSavedResidual=std::max(maxSavedResidual,std::abs(acceptedResidual[i][k]));
      maxOldDifference=std::max(maxOldDifference,std::abs(acceptedState[i][k]-acceptedOld[i][k]));
    }
  }
  REQUIRE(maxSavedResidual>0);
  REQUIRE(maxOldDifference>0);
  std::cout<<std::setprecision(17)<<"FORWARD_EULER_ROLLBACK,max_saved_residual="<<maxSavedResidual
           <<",max_old_difference="<<maxOldDifference<<'\n';
  const auto stiff=[](const Direction& x,double){return State{{x[0]>1?-100.0:0.0,0,0,0}};};
  CHECK_THROWS_AS(solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),acceptedTime,1,extrapolate,stiff),std::domain_error);
  CHECK(solver.GetPhysicalTime()==acceptedTime);
  for (unsigned long i=0;i<points.size();++i) {
    CHECK(solver.GetPrimitiveNodes().GetPrimitiveState(i)==acceptedState[i]);
    for (unsigned k=0;k<4;++k) {
      CHECK(solver.GetPrimitiveNodes().GetSolution_Old(i,k)==acceptedOld[i][k]);
      CHECK(solver.LinSysRes(i,k)==acceptedResidual[i][k]);
    }
  }
  const double retryDt=0.001;
  CPrimitiveEulerSolver::SpatialResult retry;
  REQUIRE_NOTHROW(retry=solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),acceptedTime,retryDt,extrapolate,stiff));
  CHECK(solver.GetPhysicalTime()==acceptedTime+retryDt);
  for (unsigned long i=0;i<points.size();++i) {
    // The accepted state is spatially uniform and the boundary extrapolates it:
    // the independently known rhs is precisely the prescribed stiff source.
    Check(retry.rhs[i],stiff(points[i],acceptedTime));
    State expected=acceptedState[i];
    for (unsigned k=0;k<4;++k) {
      expected[k]+=retryDt*retry.rhs[i][k];
      CHECK(solver.GetPrimitiveNodes().GetSolution_Old(i,k)==acceptedState[i][k]);
      CHECK(solver.LinSysRes(i,k)==-solver.GetGeometryMetrics().GetCell(i).volume*retry.rhs[i][k]);
    }
    Check(solver.GetPrimitiveNodes().GetPrimitiveState(i),expected);
  }
}

TEST_CASE("Native explicit hooks use a global step recompute residual and preserve constants", "[primitive-solver][primitive-time]") {
  MeshFixture fixture(points,elements,"PRIMITIVE_EULER",options,"MARKER_FAR");
  CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);CSolver* base=&solver;
  CHECK_THROWS_AS(base->ExplicitEuler_Iteration(fixture.geometry.get(),nullptr,fixture.config.get()),std::invalid_argument);
  base->SetTime_Step(fixture.geometry.get(),nullptr,fixture.config.get(),0,0);
  const double dt=solver.GetPrimitiveNodes().GetDelta_Time(0);
  CHECK(dt>0);CHECK(solver.GetMin_Delta_Time()==dt);CHECK(solver.GetMax_Delta_Time()==dt);
  for (unsigned long i=0;i<points.size();++i) {
    CHECK(solver.GetPrimitiveNodes().GetDelta_Time(i)==dt);
    for (unsigned k=0;k<4;++k) solver.LinSysRes(i,k)=1e6;
  }
  for (unsigned step=0;step<20;++step) base->ExplicitEuler_Iteration(fixture.geometry.get(),nullptr,fixture.config.get());
  CHECK(solver.GetPhysicalTime()==Approx(20*dt).epsilon(1e-13));
  for (unsigned long i=0;i<points.size();++i) Check(solver.GetPrimitiveNodes().GetPrimitiveState(i),{{1,0.5,0.3,1}});
  solver.GetPrimitiveNodes().SetDelta_Time(1,0.5*dt);
  CHECK_THROWS_AS(base->ExplicitEuler_Iteration(fixture.geometry.get(),nullptr,fixture.config.get()),std::invalid_argument);
  CHECK_THROWS_AS(base->SetTime_Step(fixture.geometry.get(),nullptr,fixture.config.get(),1,0),std::invalid_argument);
  CHECK_THROWS_AS(base->SetTime_Step(nullptr,nullptr,nullptr,0,0),std::invalid_argument);
}

TEST_CASE("Forward Euler measures first temporal order with exact smooth time forcing", "[primitive-solver][primitive-time]") {
  const State initial{{1,0.5,0.3,1}};const double end=0.2;
  std::cout<<"TEMPORAL_MMS,gradient,limiter,steps,dt,L2,order\n";
  for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"})
    for (const std::string limiter:{"NONE","VENKATAKRISHNAN","BARTH_JESPERSEN"}) {
      double previous=0;
      for (unsigned steps:{20u,40u,80u,160u}) {
        const std::string setup="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
            "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\nNUM_METHOD_GRAD_RECON="+method+"\nSLOPE_LIMITER_FLOW="+limiter+"\n";
        MeshFixture fixture(points,elements,"PRIMITIVE_EULER",setup);
        CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
        const double dt=end/steps;
        for (unsigned step=0;step<steps;++step) {
          const double time=step*dt;
          solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),time,dt,extrapolate,
              [&](const Direction&,double t){State s=initial;for(auto& x:s)x*=std::exp(t);return s;});
        }
        double squared=0,volume=0;
        for (unsigned long i=0;i<points.size();++i) {
          const auto w=solver.GetPrimitiveNodes().GetPrimitiveState(i);const double v=solver.GetGeometryMetrics().GetCell(i).volume;
          for (unsigned k=0;k<4;++k) squared+=v*std::pow(w[k]-initial[k]*std::exp(end),2);
          volume+=v;
        }
        const double error=std::sqrt(squared/volume),order=previous>0?std::log2(previous/error):0;
        if (previous>0) {CHECK(order>0.99);CHECK(order<1.01);}
        std::cout<<std::setprecision(12)<<"TEMPORAL_MMS,"<<method<<','<<limiter<<','<<steps<<','<<dt<<','<<error<<','<<order<<'\n';
        previous=error;
      }
    }
}

TEST_CASE("Time dependent affine MMS evolves the complete primitive spatial operator", "[primitive-solver][primitive-time]") {
  const Gradient g{{{{0.2,-0.1}},{{0.12,0.07}},{{-0.08,0.09}},{{0.3,0.11}}}};
  const auto exact=[&](const Direction& x,double t) {auto w=Reconstruct({{2,0.5,0.3,3}},g,x);for(auto& v:w)v*=std::exp(t);return w;};
  std::cout<<"COUPLED_TEMPORAL_MMS,gradient,steps,dt,L2,order\n";
  for (const std::string method:{"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"}) {
    double previous=0;
    for (unsigned steps:{32u,64u,128u,256u}) {
      const std::string setup="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
          "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\nSLOPE_LIMITER_FLOW=NONE\nNUM_METHOD_GRAD_RECON="+method+"\n";
      MeshFixture fixture(points,elements,"PRIMITIVE_EULER",setup);
      CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
      for (unsigned long i=0;i<points.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,exact(points[i],0));
      const double end=0.05,dt=end/steps;
      const auto boundary=[&](const Direction& x,double t,const State&,const Direction&,unsigned short){return exact(x,t);};
      const auto source=[&](const Direction& x,double t) {
        const auto w=exact(x,t);auto gradient=g;for(auto& row:gradient)for(auto& v:row)v*=std::exp(t);
        auto s=VolumeProduct(w,gradient,1.4);for(unsigned k=0;k<4;++k)s[k]+=w[k];return s;
      };
      if (steps==32) {
        for (double time:{0.0,0.037}) {
          double scale=1,maxDefect=0,squared=0,volume=0;
          unsigned long worstNode=0;unsigned worstComponent=0;
          for (unsigned long i=0;i<points.size();++i) {
            const auto w=exact(points[i],time);
            solver.GetPrimitiveNodes().SetPrimitiveState(i,w);
            for (auto value:w) scale=std::max(scale,std::abs(value));
          }
          const auto spatial=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),time,boundary,source);
          // A fixed, well-conditioned five-node patch with O(1) states and
          // derivatives: allow 64 state-scaled roundoff units, not a truncation
          // error allowance. Use direct absolute comparisons without Approx's
          // default relative tolerance.
          const double tolerance=64*std::numeric_limits<double>::epsilon()*scale;
          for (unsigned long i=0;i<points.size();++i) {
            const auto derivative=exact(points[i],time); // d/dt exp(t)*W_base = W.
            const double v=solver.GetGeometryMetrics().GetCell(i).volume;
            volume+=v;
            for (unsigned k=0;k<4;++k) {
              const double defect=std::abs(spatial.rhs[i][k]-derivative[k]);
              if (defect>maxDefect) {maxDefect=defect;worstNode=i;worstComponent=k;}
              squared+=v*defect*defect;
              INFO(method<<" time="<<time<<" node="<<i<<" component="<<k);
              CHECK(defect<=tolerance);
            }
          }
          std::cout<<std::setprecision(17)<<"AFFINE_SEMIDISCRETE,"<<method<<",time="<<time
                   <<",max_defect="<<maxDefect<<",L2="<<std::sqrt(squared/volume)
                   <<",tolerance="<<tolerance<<",node="<<worstNode<<",component="<<worstComponent<<'\n';
        }
        // Restore the original initial data; the temporal experiment below,
        // its manufactured data, step counts and order thresholds are unchanged.
        for (unsigned long i=0;i<points.size();++i)
          solver.GetPrimitiveNodes().SetPrimitiveState(i,exact(points[i],0));
      }
      for (unsigned step=0;step<steps;++step)
        solver.AdvanceExplicitEuler(fixture.geometry.get(),fixture.config.get(),step*dt,dt,boundary,source);
      double squared=0,volume=0;
      for (unsigned long i=0;i<points.size();++i) {
        const auto actual=solver.GetPrimitiveNodes().GetPrimitiveState(i),expected=exact(points[i],end);
        const double v=solver.GetGeometryMetrics().GetCell(i).volume;
        for(unsigned k=0;k<4;++k)squared+=v*std::pow(actual[k]-expected[k],2);
        volume+=v;
      }
      const double error=std::sqrt(squared/volume),order=previous>0?std::log2(previous/error):0;
      if(previous>0){CHECK(order>0.98);CHECK(order<1.02);}
      std::cout<<std::setprecision(12)<<"COUPLED_TEMPORAL_MMS,"<<method<<','<<steps<<','<<dt<<','<<error<<','<<order<<'\n';
      previous=error;
    }
  }
}
