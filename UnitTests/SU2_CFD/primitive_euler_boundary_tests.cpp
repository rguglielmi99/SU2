/*!
 * \brief Physical boundary invariants and constrained-slip operator tests.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../SU2_CFD/include/numerics/flow/primitive_euler_boundary.hpp"

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

TEST_CASE("Curved slip markers fail explicitly instead of imposing artificial no-slip", "[primitive-boundary]") {
  const std::vector<Direction> xy{{{0,0}},{{1,0}},{{0,1}}};
  MeshFixture f(xy,{{{0,1,2}}},"PRIMITIVE_EULER","TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n","MARKER_EULER");
  CPrimitiveEulerSolver solver(f.geometry.get(),f.config.get(),0);
  for (unsigned i=0;i<3;++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,{{1,0,0,1}});
  try { solver.ComputeSpatialResidual(f.geometry.get(),f.config.get());FAIL("Curved marker accepted"); }
  catch (const std::invalid_argument& e) { CHECK(std::string(e.what()).find("curved-wall")!=std::string::npos); }
}
