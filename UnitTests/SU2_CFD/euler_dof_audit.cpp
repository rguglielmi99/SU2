/*!
 * \brief Isolated phase-6A.1 experiment on the production CEulerSolver.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CEulerSolver.hpp"
#include <iostream>
#include <stdexcept>

using namespace PrimitiveEulerTest;
using Vec = std::array<double, 4>;
constexpr double gammaGas = 1.4, ux = 0.5, uy = 0.3, finalTime = 0.05;
constexpr double pi = 3.14159265358979323846;
bool affineContact=false;  // Independent control, selected only by the standalone harness.
void Require(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
Vec Primitive(const Direction& x, double t) {
  if (affineContact) return {{1+0.2*(x[0]-ux*t)+0.15*(x[1]-uy*t),ux,uy,1}};
  return {{1 + 0.2*std::sin(2*pi*(x[0]-ux*t))*std::cos(2*pi*(x[1]-uy*t)), ux, uy, 1}};
}
Vec Conservative(const Vec& w) {
  return {{w[0],w[0]*w[1],w[0]*w[2],w[3]/(gammaGas-1)+0.5*w[0]*(w[1]*w[1]+w[2]*w[2])}};
}
struct Contact : CVerificationSolution {
  explicit Contact(CConfig* config) : CVerificationSolution(2,4,0,config) {}
  void GetSolution(const su2double* x, su2double t, su2double* u) const override {
    const auto value=Conservative(Primitive({{SU2_TYPE::GetValue(x[0]),SU2_TYPE::GetValue(x[1])}},SU2_TYPE::GetValue(t)));
    for (unsigned k=0;k<4;++k) u[k]=value[k];
  }
  void GetBCState(const su2double* x, su2double t, su2double* u) const override { GetSolution(x,t,u); }
};
struct Grid {
  unsigned n;
  std::vector<Direction> xy;
  std::vector<Triangle> triangles;
  Grid(unsigned count, const std::string& family) : n(count) {
    for (unsigned j=0;j<=n;++j) for (unsigned i=0;i<=n;++i) {
      Direction x{{double(i)/n,double(j)/n}};
      if (family=="smooth") {
        const double s=std::sin(pi*x[0])*std::sin(pi*x[1]);
        x[0]+=0.08*s; x[1]+=0.06*s;
      } else if (family=="jitter" && i && j && i<n && j<n) {
        x[0]+=0.18/n*std::sin(2.1*i+1.3*j); x[1]+=0.16/n*std::cos(1.7*i-0.9*j);
      }
      xy.push_back(x);
    }
    for (unsigned j=0;j<n;++j) for (unsigned i=0;i<n;++i) {
      const unsigned long a=j*(n+1)+i,b=a+1,d=a+n+1,c=d+1;
      triangles.push_back({{a,b,c}}); triangles.push_back({{a,c,d}});
    }
  }
  unsigned layer(unsigned long i) const {
    const unsigned x=i%(n+1),y=i/(n+1); return std::min({x,y,n-x,n-y});
  }
  bool in(unsigned long i,unsigned g) const {
    const unsigned l=layer(i);
    if (g==0) return true;
    if (g<4) return l==g-1;
    if (g==4) return l>=3;
    if (g==6) return l>=1;
    const auto& x=xy[i]; return x[0]>=0.25 && x[0]<=0.75 && x[1]>=0.25 && x[1]<=0.75;
  }
};
using Subtriangle=std::array<Direction,3>;
struct Dual {
  std::vector<Subtriangle> pieces;
  double volume=0;
  Direction offset{};
  std::array<double,3> moment{};
};
double Area(const Subtriangle& p) {
  return 0.5*std::abs((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]));
}
std::vector<std::pair<double,double>> Gauss(unsigned n) {
  std::vector<std::pair<double,double>> rule;
  for (unsigned i=0;i<n;++i) {
    double z=std::cos(pi*(i+0.75)/(n+0.5)), derivative=0;
    for (unsigned it=0;it<30;++it) {
      double a=1,b=0;
      for (unsigned k=1;k<=n;++k) { const double c=b; b=a; a=((2*k-1)*z*b-(k-1)*c)/k; }
      derivative=n*(z*a-b)/(z*z-1);
      const double next=z-a/derivative;
      if (std::abs(next-z)<2e-16) { z=next; break; } z=next;
    }
    rule.push_back({0.5*(1+z),1/((1-z*z)*derivative*derivative)});
  }
  return rule;
}
template<class Function> Vec Integrate(const Dual& cell,unsigned order,const Function& f) {
  const auto rule=Gauss(order); Vec result{};
  for (const auto& p:cell.pieces) for (const auto& r:rule) for (const auto& s:rule) {
    Direction q{};
    for (unsigned d=0;d<2;++d) q[d]=p[0][d]+r.first*(p[1][d]-p[0][d])+(1-r.first)*s.first*(p[2][d]-p[0][d]);
    const auto value=f(q);
    const double weight=2*Area(p)*(1-r.first)*r.second*s.second;
    for (unsigned k=0;k<4;++k) result[k]+=weight*value[k];
  }
  for (auto& v:result) v/=cell.volume;
  return result;
}
std::vector<Dual> BuildDual(const Grid& grid,const CGeometry& geometry,double& validation) {
  std::vector<Dual> cells(grid.xy.size());
  for (const auto& t:grid.triangles) {
    Direction c{}; for (auto i:t) for (unsigned d=0;d<2;++d) c[d]+=grid.xy[i][d]/3;
    const auto& a=grid.xy[t[0]]; const auto& b=grid.xy[t[1]]; const auto& d=grid.xy[t[2]];
    Require((b[0]-a[0])*(d[1]-a[1])-(b[1]-a[1])*(d[0]-a[0])>0.2/(grid.n*grid.n),
            "inverted or degenerate triangle");
    for (unsigned k=0;k<3;++k) {
      const auto i=t[k]; Direction a{},b{};
      for (unsigned d=0;d<2;++d) { a[d]=(grid.xy[i][d]+grid.xy[t[(k+1)%3]][d])/2; b[d]=(grid.xy[i][d]+grid.xy[t[(k+2)%3]][d])/2; }
      cells[i].pieces.push_back({{grid.xy[i],a,c}}); cells[i].pieces.push_back({{grid.xy[i],c,b}});
    }
  }
  for (unsigned long i=0;i<cells.size();++i) {
    auto& cell=cells[i];
    // Exact polynomial moments from barycentric integration; independent of Duffy quadrature.
    for (const auto& p:cell.pieces) {
      const double area=Area(p); cell.volume+=area; Direction sum{};
      for (auto x:p) for (unsigned d=0;d<2;++d) sum[d]+=x[d]-grid.xy[i][d];
      for (unsigned d=0;d<2;++d) cell.offset[d]+=area*sum[d]/3;
      for (unsigned m=0;m<3;++m) {
        const unsigned a=m==2?1:0,b=m==0?0:1; double diagonal=0;
        for (auto x:p) diagonal+=(x[a]-grid.xy[i][a])*(x[b]-grid.xy[i][b]);
        cell.moment[m]+=area*(sum[a]*sum[b]+diagonal)/12;
      }
    }
    const double native=SU2_TYPE::GetValue(geometry.nodes->GetVolume(i));
    Require(std::abs(cell.volume/native-1)<2e-12,"dual volume differs from SU2");
    for (auto& x:cell.offset) x/=cell.volume;
    for (auto& x:cell.moment) x/=cell.volume;
    for (unsigned order:{3u,5u,8u}) {
      const auto poly=Integrate(cell,order,[&](const Direction& x) {
        const double a=x[0]-grid.xy[i][0],b=x[1]-grid.xy[i][1]; return Vec{{1,a,b,a*a}};
      });
      const Vec exact{{1,cell.offset[0],cell.offset[1],cell.moment[0]}};
      for (unsigned k=0;k<4;++k) validation=std::max(validation,std::abs(poly[k]-exact[k]));
      const auto mixed=Integrate(cell,order,[&](const Direction& x) {
        const double a=x[0]-grid.xy[i][0],b=x[1]-grid.xy[i][1]; return Vec{{a*b,b*b,0,0}};
      });
      validation=std::max({validation,std::abs(mixed[0]-cell.moment[1]),std::abs(mixed[1]-cell.moment[2])});
    }
  }
  Require(validation<2e-13,"polynomial quadrature validation failed");
  return cells;
}
const std::array<const char*,7> groupNames{{"all","boundary","ring1","ring2","deep3","fixed","interior"}};
struct Stats {
  double sum=0,weighted=0,volume=0,maximum=0; unsigned count=0;
  void add(double e,double v) { sum+=e*e; weighted+=v*e*e; volume+=v; maximum=std::max(maximum,std::abs(e)); ++count; }
  double rms() const { return std::sqrt(sum/count); }
  double l2() const { return std::sqrt(weighted/volume); }
};
int Run(int argc,char** argv) {
  Require(argc>=6,"usage: audit family N point|mean dt_divisor output.csv [geometry|affine]");
  Require(argc==6 || (argc==7 && (std::string(argv[6])=="geometry" || std::string(argv[6])=="affine")),
          "unknown audit mode");
  const std::string family=argv[1],initial=argv[3]; const unsigned n=std::stoul(argv[2]),divisor=std::stoul(argv[4]);
  affineContact=argc>6 && std::string(argv[6])=="affine";
  Require(SU2_MPI::GetSize()==1,"this audit requires one MPI rank");
  Require(n>=8 && divisor>0 && (initial=="point" || initial=="mean"),"invalid arguments");
  Require(family=="regular" || family=="smooth" || family=="jitter","invalid family");
  Grid grid(n,family);
  const std::string options="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
      "MUSCL_FLOW=YES\nMUSCL_KAPPA_FLOW=0\nSLOPE_LIMITER_FLOW=NONE\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\n"
      "FLUID_MODEL=IDEAL_GAS\nGAMMA_VALUE=1.4\nREF_DIMENSIONALIZATION=DIMENSIONAL\n"
      "INIT_OPTION=TD_CONDITIONS\nFREESTREAM_OPTION=DENSITY_FS\nFREESTREAM_DENSITY=1\nFREESTREAM_PRESSURE=1\n"
      "MACH_NUMBER=0.4\nTIME_DOMAIN=YES\nTIME_MARCHING=TIME_STEPPING\nTIME_STEP=0.001\nUNST_CFL_NUMBER=0\n";
  MeshFixture fixture(grid.xy,grid.triangles,"EULER",options);
  auto* geometry=fixture.geometry.get(); auto* config=fixture.config.get();
  config->SetGlobalParam(MAIN_SOLVER::EULER,RUNTIME_FLOW_SYS);
  Require(config->GetKind_TimeIntScheme()==EULER_EXPLICIT,"active integration scheme mismatch");
  geometry->PreprocessP2PComms(geometry,config);
  geometry->PreprocessPeriodicComms(geometry,config);
  double validation=0; const auto cells=BuildDual(grid,*geometry,validation);
  std::vector<Vec> mean0(cells.size()),meanT(cells.size()),meanW(cells.size()); double quadrature=0;
  for (unsigned long i=0;i<cells.size();++i) {
    mean0[i]=Integrate(cells[i],8,[](const Direction& x){return Conservative(Primitive(x,0));});
    meanT[i]=Integrate(cells[i],8,[](const Direction& x){return Conservative(Primitive(x,finalTime));});
    meanW[i]=Integrate(cells[i],8,[](const Direction& x){return Primitive(x,finalTime);});
    for (double t:{0.0,finalTime}) {
      const auto a=Integrate(cells[i],5,[&](const Direction& x){return Conservative(Primitive(x,t));});
      const auto& b=t==0?mean0[i]:meanT[i];
      for (unsigned k=0;k<4;++k) quadrature=std::max(quadrature,std::abs(a[k]-b[k]));
    }
  }
  Require(quadrature<1e-10,"smooth quadrature unresolved");
  std::ofstream out(argv[5]); out << std::setprecision(17);
  out << "kind,family,N,initial,divisor,group,component,rms,l2,linf,extra\n";
  const auto emit=[&](const char* kind,unsigned group,unsigned k,const Stats& s,double extra=0) {
    Require(s.count>0,"empty diagnostic group");
    out << kind << ',' << family << ',' << n << ',' << initial << ',' << divisor << ',' << groupNames[group] << ',' << k
        << ',' << s.rms() << ',' << s.l2() << ',' << s.maximum << ',' << extra << '\n';
  };
  for (unsigned g=0;g<7;++g) {
    Stats offset,scaled,moment,nonlinear; std::array<Stats,3> components; unsigned asymmetric=0;
    for (unsigned long i=0;i<cells.size();++i) if (grid.in(i,g)) {
      const auto& c=cells[i]; const double d=std::hypot(c.offset[0],c.offset[1]);
      offset.add(d,c.volume); scaled.add(d*n,c.volume); asymmetric+=d*n>0.01;
      moment.add(std::sqrt(c.moment[0]*c.moment[0]+2*c.moment[1]*c.moment[1]+c.moment[2]*c.moment[2]),c.volume);
      for (unsigned m=0;m<3;++m) components[m].add(c.moment[m],c.volume);
      const auto converted=Conservative(meanW[i]); double e=0;
      for (unsigned k=0;k<4;++k) e+=std::pow(meanT[i][k]-converted[k],2);
      nonlinear.add(std::sqrt(e),c.volume);
    }
    emit("offset",g,0,offset); emit("offset_over_h",g,0,scaled,double(asymmetric)/scaled.count);
    emit("moment2",g,0,moment); emit("meanU_minus_UmeanW",g,0,nonlinear);
    for (unsigned m=0;m<3;++m) emit("moment_component",g,m,components[m]);
  }
  out << "validation," << family << ',' << n << ',' << initial << ',' << divisor << ",all,0," << validation << ',' << quadrature << ",0,0\n";
  if (argc>6 && std::string(argv[6])=="geometry") return 0;
  std::vector<Vec> result(cells.size()); double evolved=0,minrho=1e30,minp=1e30;
  const unsigned steps=static_cast<unsigned>(std::ceil(finalTime*n*n/0.25))*divisor;
  const double dt=finalTime/steps;
  {
    QuietGeometryOutput quiet;
    CEulerSolver solver(geometry,config,0);
    Require(solver.VerificationSolution==nullptr,"unexpected verification solution");
    solver.VerificationSolution=new Contact(config);
    CSolver* solvers[MAX_SOLS]{}; solvers[FLOW_SOL]=&solver;
    auto* nodes=solver.GetNodes();
    Require(std::abs(SU2_TYPE::GetValue(config->GetGamma())-gammaGas)<1e-14,"gamma mismatch");
    Require(config->GetTime_Ref()==1 && config->GetDensity_Ref()==1 && config->GetPressure_Ref()==1,"normalization mismatch");
    for (unsigned long i=0;i<cells.size();++i) {
      const auto u=initial=="point"?Conservative(Primitive(grid.xy[i],0)):mean0[i];
      for (unsigned k=0;k<4;++k) nodes->SetSolution(i,k,u[k]);
      nodes->SetDelta_Time(i,dt);
    }
    for (unsigned step=0;step<steps;++step) {
      config->SetPhysicalTime(step*dt);
      nodes->Set_OldSolution();
      solver.BC_Custom(geometry,solvers,nullptr,nullptr,config,0);
      solver.Preprocessing(geometry,solvers,config,0,0,RUNTIME_FLOW_SYS,false);
      Require(config->GetNonphysical_Points()==0,"nonphysical nodal fallback");
      solver.Upwind_Residual(geometry,solvers,nullptr,config,0);
      Require(config->GetNonphysical_Reconstr()==0,"nonphysical reconstruction fallback");
      solver.BC_Custom(geometry,solvers,nullptr,nullptr,config,0);
      solver.ExplicitEuler_Iteration(geometry,solvers,config);
      // Strong time-dependent Dirichlet constraint at the new time, outside the interior update.
      config->SetPhysicalTime((step+1)*dt);
      solver.BC_Custom(geometry,solvers,nullptr,nullptr,config,0);
    }
    Require(std::abs(SU2_TYPE::GetValue(config->GetPhysicalTime())-finalTime)<1e-15,"final time mismatch");
    for (unsigned long i=0;i<cells.size();++i) {
      for (unsigned k=0;k<4;++k) result[i][k]=SU2_TYPE::GetValue(nodes->GetSolution(i,k));
      const auto& u=result[i]; const double p=(gammaGas-1)*(u[3]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0]);
      minrho=std::min(minrho,u[0]); minp=std::min(minp,p);
      evolved=std::max(evolved,std::abs(u[0]-(initial=="point"?Conservative(Primitive(grid.xy[i],0)):mean0[i])[0]));
      for (auto v:u) Require(std::isfinite(v),"nonfinite solution");
      Require(u[0]>0 && p>0,"nonpositive solution");
    }
  }
  Require(evolved>1e-4,"solution did not evolve");
  for (unsigned g=0;g<7;++g) for (unsigned k=0;k<4;++k) {
    Stats point,mean,difference;
    for (unsigned long i=0;i<cells.size();++i) if (grid.in(i,g)) {
      const double p=Conservative(Primitive(grid.xy[i],finalTime))[k],m=meanT[i][k],v=cells[i].volume;
      point.add(result[i][k]-p,v); mean.add(result[i][k]-m,v); difference.add(m-p,v);
    }
    emit("Epoint",g,k,point); emit("Emean",g,k,mean); emit("D",g,k,difference);
  }
  out << "evolution," << family << ',' << n << ',' << initial << ',' << divisor << ",all,0," << dt << ',' << evolved << ',' << minrho << ',' << minp << '\n';
  // Field snapshots make temporal and initialization comparisons independent of error norms.
  std::ofstream field(std::string(argv[5])+".field"); field << std::setprecision(17);
  for (unsigned long i=0;i<cells.size();++i) {
    field << i << ' ' << cells[i].volume << ' ' << grid.xy[i][0] << ' ' << grid.xy[i][1];
    for (auto v:result[i]) field << ' ' << v;
    field << '\n';
  }
  return 0;
}
int main(int argc,char** argv) {
  SU2_MPI::Init(&argc,&argv); int status=0;
  try { status=Run(argc,argv); }
  catch (const std::exception& e) { std::cerr << "AUDIT FAILED: " << e.what() << '\n'; status=1; }
  SU2_MPI::Finalize(); return status;
}
