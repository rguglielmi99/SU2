/*!
 * \brief Source-free contact transport with real native Euler BCs and primitive BCs.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CEulerSolver.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../SU2_CFD/include/numerics/flow/convection/roe.hpp"
#include "../../SU2_CFD/include/iteration/primitive_euler_boundary_config.hpp"
#include "../../Common/include/toolboxes/MMS/CVerificationSolution.hpp"
#include <chrono>

using namespace PrimitiveEuler2D;
using namespace PrimitiveEulerTest;
namespace {
double rotation=0;
std::string experiment;
double amplitude=1e-5;
double Pulse(double x) { return std::abs(x)<1.5?std::pow(1-x*x/2.25,4):0; }
State Exact(const Direction& x,double t,double speed) {
  if (experiment.find("vorticity_")==0) return {{1,speed,.01*Pulse(x[0]-3.5-speed*t),1}};
  if (experiment.find("acoustic_")==0) {
    const double c=std::sqrt(1.4), u=speed;
    const bool inlet=experiment=="acoustic_total" || experiment=="acoustic_mass";
    const double direction=inlet?-1.:1., incidentSpeed=u+direction*c;
    const double center=inlet?-2.5:2.5;
    const double incident=amplitude*Pulse(x[0]-center-incidentSpeed*t);
    const double boundary=inlet?-5.:5.;
    const double reflectedSpeed=u-direction*c;
    const double delay=(x[0]-boundary)/reflectedSpeed;
    double reflection=0;
    if (experiment=="acoustic_total") reflection=-(1-u/c)/(1+u/c);
    if (experiment=="acoustic_mass" || experiment=="acoustic_wall") reflection=1;
    if (experiment=="acoustic_outlet") reflection=-1;
    const double reflected=delay>=0 && delay<=t ? reflection*amplitude*Pulse(boundary-center-incidentSpeed*(t-delay)):0;
    double contact=0;
    if (experiment=="acoustic_mass") {
      const double convectiveDelay=(x[0]+5)/u;
      if (convectiveDelay>=0 && convectiveDelay<=t)
        contact=-2*amplitude*Pulse(-5-center-incidentSpeed*(t-convectiveDelay))/(c*c);
    }
    return {{1+(incident+reflected)/(c*c)+contact,u+direction*(incident-reflected)/c,0,1+incident+reflected}};
  }
  const double X=x[0]*std::cos(rotation)+x[1]*std::sin(rotation)-3.5-speed*t;
  const double Y=-x[0]*std::sin(rotation)+x[1]*std::cos(rotation);
  const double r2=(X*X+Y*Y)/4.;
  const double bump=r2<1?std::pow(1-r2,4):0;
  return {{1+.02*bump,speed*std::cos(rotation),speed*std::sin(rotation),1}};
}
State Conservative(const State& w) { return {{w[0],w[0]*w[1],w[0]*w[2],w[3]/.4+.5*w[0]*(w[1]*w[1]+w[2]*w[2])}}; }
State Primitive(CVariable* nodes,unsigned long i) {
  const double rho=SU2_TYPE::GetValue(nodes->GetSolution(i,0));
  const double u=SU2_TYPE::GetValue(nodes->GetSolution(i,1))/rho,v=SU2_TYPE::GetValue(nodes->GetSolution(i,2))/rho;
  return {{rho,u,v,.4*(SU2_TYPE::GetValue(nodes->GetSolution(i,3))-.5*rho*(u*u+v*v))}};
}
struct AnalyticSideData final : CVerificationSolution {
  double speed;
  AnalyticSideData(CConfig* c,double u) : CVerificationSolution(2,4,0,c),speed(u) {}
  bool IsManufacturedSolution() const override { return false; }
  void GetSolution(const su2double* x,su2double t,su2double* u) const override {
    const auto state=Conservative(Exact({{SU2_TYPE::GetValue(x[0]),SU2_TYPE::GetValue(x[1])}},SU2_TYPE::GetValue(t),speed));
    for (unsigned k=0;k<4;++k) u[k]=state[k];
  }
  void GetBCState(const su2double* x,su2double t,su2double* u) const override { GetSolution(x,t,u); }
};
}
int main(int argc,char** argv) {
  SU2_MPI::Init(&argc,&argv);
  int status=0;
  try {
    if (argc!=6 && argc!=7) throw std::invalid_argument("family N mode divisor output.csv [amplitude] required");
    const bool jitter=std::string(argv[1])=="jitter";
    const unsigned n=std::stoul(argv[2]), divisor=std::stoul(argv[4]);
    const std::string mode=argv[3];
    experiment=mode;
    if (argc==7) amplitude=std::stod(argv[6]);
    rotation=mode=="wall_oblique"?.37:0;
    const bool acoustic=mode.find("acoustic_")==0;
    const bool vorticity=mode.find("vorticity_")==0;
    const double speed=mode=="supersonic" || mode=="vorticity_supersonic"?2:(mode=="acoustic_wall"?0:.4);
    const double T=acoustic?6.:.5;
    std::vector<Direction> xy;std::vector<Triangle> triangles;
    std::vector<BoundaryMarker> markers{{"lateral",{}},{"inlet",{}},{"outlet",{}}};
    const auto id=[=](unsigned i,unsigned j){return static_cast<unsigned long>(j*(n+1)+i);};
    for (unsigned j=0;j<=n;++j) for (unsigned i=0;i<=n;++i) {
      Direction x{{-5.+10.*i/n,-5.+10.*j/n}};
      if (jitter && i>0 && i<n && j>0 && j<n) {
        x[0]+=.18*10/n*std::sin(2.1*i+1.3*j);x[1]+=.16*10/n*std::cos(1.7*i-.9*j);
      }
      xy.push_back(x);
    }
    for (auto& x:xy) x={{std::cos(rotation)*x[0]-std::sin(rotation)*x[1],std::sin(rotation)*x[0]+std::cos(rotation)*x[1]}};
    for (unsigned j=0;j<n;++j) for (unsigned i=0;i<n;++i) {
      triangles.push_back({{id(i,j),id(i+1,j),id(i+1,j+1)}});
      triangles.push_back({{id(i,j),id(i+1,j+1),id(i,j+1)}});
    }
    for (unsigned i=0;i<n;++i) {markers[0].edges.push_back({{id(i,0),id(i+1,0)}});markers[0].edges.push_back({{id(i+1,n),id(i,n)}});}
    for (unsigned j=0;j<n;++j) {markers[1].edges.push_back({{id(0,j+1),id(0,j)}});markers[2].edges.push_back({{id(n,j),id(n,j+1)}});}
    std::ostringstream cfg;cfg<<std::setprecision(17)
      <<"TIME_DOMAIN=YES\nTIME_MARCHING=TIME_STEPPING\nTIME_DISCRE_FLOW=EULER_EXPLICIT\nTIME_STEP=0.001\n"
        "CONV_NUM_METHOD_FLOW=ROE\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\n"
        "SLOPE_LIMITER_FLOW=NONE\nMUSCL_FLOW=YES\nMUSCL_KAPPA_FLOW=0\n"
        "FLUID_MODEL=IDEAL_GAS\nGAMMA_VALUE=1.4\nGAS_CONSTANT=1\nREF_DIMENSIONALIZATION=DIMENSIONAL\n"
        "INIT_OPTION=TD_CONDITIONS\nFREESTREAM_OPTION=DENSITY_FS\nFREESTREAM_DENSITY=1\nFREESTREAM_PRESSURE=1\nMACH_NUMBER="<<speed/std::sqrt(1.4)<<"\n"
        "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\nPRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC\nAOA="<<rotation*180/std::acos(-1.)<<"\n";
    if (mode=="farfield") cfg<<"MARKER_FAR=(lateral,inlet,outlet)\n";
    else {
      cfg<<(vorticity?"MARKER_CUSTOM=(lateral)\n":(mode=="symmetry"?"MARKER_SYM=(lateral)\n":
          (mode=="acoustic_wall"?"MARKER_EULER=(lateral,outlet)\n":"MARKER_EULER=(lateral)\n")));
      if (mode=="wall" || mode=="wall_oblique" || mode=="symmetry" || mode=="acoustic_far" || mode=="vorticity_far") cfg<<"MARKER_FAR=(inlet,outlet)\n";
      else if (mode=="acoustic_wall") cfg<<"MARKER_FAR=(inlet)\n";
      else if (mode=="supersonic" || mode=="vorticity_supersonic") cfg<<"MARKER_SUPERSONIC_INLET=(inlet,1,1,2,0,0)\nMARKER_SUPERSONIC_OUTLET=(outlet)\n";
      else {
        if (mode=="mass" || mode=="acoustic_mass") cfg<<"INLET_TYPE=MASS_FLOW\nMARKER_INLET=(inlet,1,"<<speed<<",1,0,0)\n";
        else {
          const double tt=1+.4*speed*speed/(2*1.4);
          cfg<<"INLET_TYPE=TOTAL_CONDITIONS\nMARKER_INLET=(inlet,"<<tt<<','<<std::pow(tt,3.5)<<",1,0,0)\n";
        }
        cfg<<"MARKER_OUTLET=(outlet,1)\n";
      }
    }
    MeshFixture pf(xy,triangles,"PRIMITIVE_EULER",cfg.str(),"",markers,false,false);
    MeshFixture cf(xy,triangles,"EULER",cfg.str(),"",markers,false,false);
    pf.config->SetDensity_Ref(1);pf.config->SetPressure_Ref(1);pf.config->SetVelocity_Ref(1);pf.config->SetTime_Ref(1);
    pf.config->SetGlobalParam(MAIN_SOLVER::PRIMITIVE_EULER,RUNTIME_FLOW_SYS);
    cf.config->SetGlobalParam(MAIN_SOLVER::EULER,RUNTIME_FLOW_SYS);
    cf.geometry->PreprocessP2PComms(cf.geometry.get(),cf.config.get());
    cf.geometry->PreprocessPeriodicComms(cf.geometry.get(),cf.config.get());
    cf.geometry->ComputeModifiedSymmetryNormals(cf.config.get());
    std::unique_ptr<CEulerSolver> conservative;
    conservative.reset(new CEulerSolver(cf.geometry.get(),cf.config.get(),0));
    if (vorticity) conservative->VerificationSolution=new AnalyticSideData(cf.config.get(),speed);
    for (unsigned short m=0;m<cf.geometry->GetnMarker();++m)
      conservative->SetUniformInlet(cf.config.get(),m);
    CPrimitiveEulerSolver primitive(pf.geometry.get(),pf.config.get(),0);
    auto* nodes=conservative->GetNodes();
    CSolver* solvers[MAX_SOLS]{};solvers[FLOW_SOL]=conservative.get();
    CUpwRoe_Flow numericalBoundary(2,4,cf.config.get(),false);
    const unsigned steps=(acoustic?8:2)*n*n*divisor;const double dt=T/steps;
    cf.config->SetDelta_UnstTimeND(dt);cf.config->SetInnerIter(0);
    for (unsigned long i=0;i<xy.size();++i) {
      const auto w=Exact(xy[i],0,speed),u=Conservative(w);
      primitive.GetPrimitiveNodes().SetPrimitiveState(i,w);
      nodes->SetSolution(i,u.data());nodes->SetDelta_Time(i,dt);
    }
    const auto boundary=[&](const Direction& q,double t,const State& w,const Direction& normal,unsigned short marker) {
      if (vorticity && pf.config->GetMarker_All_KindBC(marker)==CUSTOM_BOUNDARY) return Exact(q,t,speed);
      return PrimitiveEulerBoundaryConfig::Evaluate(*pf.config,marker,q,t,w,normal);
    };
    const auto nativeBoundary=[&]() {
      for (unsigned short m=0;m<cf.geometry->GetnMarker();++m) {
        const auto kind=cf.config->GetMarker_All_KindBC(m);
        switch (kind) {
          case FAR_FIELD: conservative->BC_Far_Field(cf.geometry.get(),solvers,&numericalBoundary,nullptr,cf.config.get(),m);break;
          case INLET_FLOW: conservative->BC_Inlet(cf.geometry.get(),solvers,&numericalBoundary,nullptr,cf.config.get(),m);break;
          case OUTLET_FLOW: conservative->BC_Outlet(cf.geometry.get(),solvers,&numericalBoundary,nullptr,cf.config.get(),m);break;
          case SUPERSONIC_INLET: conservative->BC_Supersonic_Inlet(cf.geometry.get(),solvers,&numericalBoundary,nullptr,cf.config.get(),m);break;
          case SUPERSONIC_OUTLET: conservative->BC_Supersonic_Outlet(cf.geometry.get(),solvers,&numericalBoundary,nullptr,cf.config.get(),m);break;
          case CUSTOM_BOUNDARY: conservative->BC_Custom(cf.geometry.get(),solvers,nullptr,nullptr,cf.config.get(),m);break;
          default: break;
        }
      }
      for (unsigned short m=0;m<cf.geometry->GetnMarker();++m) {
        const auto kind=cf.config->GetMarker_All_KindBC(m);
        if (kind==EULER_WALL || kind==SYMMETRY_PLANE)
          conservative->BC_Sym_Plane(cf.geometry.get(),solvers,&numericalBoundary,nullptr,cf.config.get(),m);
      }
    };
    const auto start=std::chrono::steady_clock::now();
    std::ofstream output(argv[5]); output<<std::setprecision(17)
      <<"family,N,mode,divisor,time,kind,region,component,l1,l2,linf,min_rho,min_p,max_wall_un,seconds,amplitude,dt,steps\n";
    const auto measure=[&](double t,const CPrimitiveEulerSolver::SpatialResult* residual) {
      const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      for (const std::string kind:{"P","C","PC","tau_P"}) {
        if (kind=="tau_P" && !residual) continue;
        for (unsigned region=0;region<5;++region) for (unsigned k=0;k<4;++k) {
          double weight=0,l1=0,l2=0,linf=0,minrho=1e30,minp=1e30,maxnormal=0;
          for (unsigned long i=0;i<xy.size();++i) {
            const unsigned ix=i%(n+1),iy=i/(n+1),layer=std::min({ix,iy,n-ix,n-iy});
            if (region && (region==4 ? layer<3 : layer!=region-1)) continue;
            const auto p=primitive.GetPrimitiveNodes().GetPrimitiveState(i),c=Primitive(nodes,i),exact=Exact(xy[i],t,speed);
            const double v=primitive.GetGeometryMetrics().GetCell(i).volume;
            double error=kind=="PC"?p[k]-c[k]:(kind=="P"?p[k]-exact[k]:c[k]-exact[k]);
            if (kind=="tau_P") {
              // Density contact derivative is analytic; other primitives are constant.
              const double X=xy[i][0]*std::cos(rotation)+xy[i][1]*std::sin(rotation)-3.5-speed*t;
              const double Y=-xy[i][0]*std::sin(rotation)+xy[i][1]*std::cos(rotation),r2=(X*X+Y*Y)/4;
              const double wt=k==0 && r2<1 ? .04*speed*X*std::pow(1-r2,3):0;
              double referenceDerivative=wt;
              if (vorticity) {
                const double s=xy[i][0]-3.5-speed*t;
                referenceDerivative=k==2 && std::abs(s)<1.5 ? .01*speed*8*s/2.25*std::pow(1-s*s/2.25,3):0;
              }
              error=residual->rhs[i][k]-referenceDerivative;
            }
            weight+=v;l1+=v*std::abs(error);l2+=v*error*error;linf=std::max(linf,std::abs(error));
            const auto& state=kind=="C"?c:p;minrho=std::min(minrho,state[0]);minp=std::min(minp,state[3]);
            if (iy==0 || iy==n) maxnormal=std::max(maxnormal,std::abs(-state[1]*std::sin(rotation)+state[2]*std::cos(rotation)));
          }
          output<<(jitter?"jitter":"regular")<<','<<n<<','<<mode<<','<<divisor<<','<<t<<','<<kind<<','
            <<region<<','<<k<<','<<l1/weight<<','<<std::sqrt(l2/weight)<<','<<linf<<','<<minrho<<','<<minp<<','<<maxnormal<<','<<elapsed<<','<<(acoustic?amplitude:(vorticity?.01:.02))<<','<<dt<<','<<steps<<'\n';
        }
      }
    };
    conservative->Preprocessing(cf.geometry.get(),solvers,cf.config.get(),0,0,RUNTIME_FLOW_SYS,false);
    const auto residual=primitive.ComputeSpatialResidual(pf.geometry.get(),pf.config.get(),0,boundary);
    measure(0,acoustic?nullptr:&residual);
    for (unsigned step=0;step<steps;++step) {
      cf.config->SetPhysicalTime(step*dt);
      if (vorticity) conservative->BC_Custom(cf.geometry.get(),solvers,nullptr,nullptr,cf.config.get(),0);
      nodes->Set_OldSolution();
      conservative->Preprocessing(cf.geometry.get(),solvers,cf.config.get(),0,0,RUNTIME_FLOW_SYS,false);
      if (cf.config->GetNonphysical_Points()!=0) throw std::runtime_error("Native nonphysical nodal fallback");
      conservative->Upwind_Residual(cf.geometry.get(),solvers,nullptr,cf.config.get(),0);
      if (cf.config->GetNonphysical_Reconstr()!=0) throw std::runtime_error("Native nonphysical reconstruction fallback");
      nativeBoundary();
      conservative->ExplicitEuler_Iteration(cf.geometry.get(),solvers,cf.config.get());
      for (unsigned long i=0;i<xy.size();++i) Validate(Primitive(nodes,i),1.4);
      if (vorticity) {
        cf.config->SetPhysicalTime((step+1)*dt);
        conservative->BC_Custom(cf.geometry.get(),solvers,nullptr,nullptr,cf.config.get(),0);
      }
      primitive.AdvanceExplicitEuler(pf.geometry.get(),pf.config.get(),step*dt,dt,boundary);
    }
    measure(T,nullptr);
    if (acoustic) {
      std::ofstream waves(std::string(argv[5])+".waves.csv");waves<<std::setprecision(17)<<"x,y,volume,rho_P,u_P,p_P,rho_C,u_C,p_C,rho_exact,u_exact,p_exact\n";
      for (unsigned long i=0;i<xy.size();++i) {
        const auto p=primitive.GetPrimitiveNodes().GetPrimitiveState(i),c=Primitive(nodes,i),exact=Exact(xy[i],T,speed);
        waves<<xy[i][0]<<','<<xy[i][1]<<','<<primitive.GetGeometryMetrics().GetCell(i).volume<<','<<p[0]<<','<<p[1]<<','<<p[3]<<','<<c[0]<<','<<c[1]<<','<<c[3]<<','<<exact[0]<<','<<exact[1]<<','<<exact[3]<<'\n';
      }
    }
  } catch (const std::exception& e) { std::cerr<<e.what()<<'\n';status=1; }
  SU2_MPI::Finalize();return status;
}
