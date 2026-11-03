/*!
 * \brief Phase-6B.2 evolved comparison of two production Euler solvers.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "euler_mms_test_support.hpp"
#include "../../SU2_CFD/include/solvers/CEulerSolver.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include <limits>

namespace M=EulerCommonMMS;
using namespace EulerMMSAudit;

struct Manufactured : CVerificationSolution {
  mutable unsigned long sourceCalls=0,boundaryCalls=0;
  mutable double sourceTimeError=0,boundaryTimeError=0;
  double expectedTime=0;
  explicit Manufactured(CConfig* c) : CVerificationSolution(2,4,0,c) {}
  bool IsManufacturedSolution() const override { return true; }
  void GetSolution(const su2double* x,su2double t,su2double* u) const override {
    const auto value=M::Conservative(M::Evaluate(SU2_TYPE::GetValue(x[0]),SU2_TYPE::GetValue(x[1]),SU2_TYPE::GetValue(t)).w);
    for (unsigned k=0;k<4;++k) u[k]=value[k];
  }
  void GetBCState(const su2double* x,su2double t,su2double* u) const override {
    ++boundaryCalls; boundaryTimeError=std::max(boundaryTimeError,std::abs(SU2_TYPE::GetValue(t)-expectedTime));
    GetSolution(x,t,u);
  }
  void GetMMSSourceTerm(const su2double* x,su2double t,su2double* source) const override {
    ++sourceCalls; sourceTimeError=std::max(sourceTimeError,std::abs(SU2_TYPE::GetValue(t)-expectedTime));
    const auto value=M::ConservativeSource(M::Evaluate(SU2_TYPE::GetValue(x[0]),SU2_TYPE::GetValue(x[1]),SU2_TYPE::GetValue(t)));
    for (unsigned k=0;k<4;++k) source[k]=value[k];
  }
};
struct Norm {
  unsigned long count=0; double weight=0,square=0,weighted=0,absolute=0,maximum=0;
  void add(double e,double v) {
    Require(std::isfinite(e),"nonfinite diagnostic");
    ++count; weight+=v; square+=e*e; weighted+=v*e*e; absolute+=v*std::abs(e); maximum=std::max(maximum,std::abs(e));
  }
};

Vec PrimitiveFromU(const Vec& u) {
  Require(std::isfinite(u[0]) && u[0]>0,"nonpositive conservative density");
  Vec w{{u[0],u[1]/u[0],u[2]/u[0],(M::gamma-1)*(u[3]-(u[1]*u[1]+u[2]*u[2])/(2*u[0]))}};
  for (auto value:w) Require(std::isfinite(value),"nonfinite evolved state");
  Require(w[3]>0,"nonpositive conservative pressure"); return w;
}
int Run(int argc,char** argv) {
  Require(argc==5,"usage: euler_mms_evolution family N divisor output.csv");
  Require(SU2_MPI::GetSize()==1,"one MPI rank is required");
  const std::string family=argv[1]; const unsigned n=std::stoul(argv[2]),divisor=std::stoul(argv[3]);
  Require(n>=8 && divisor>0 && (family=="regular" || family=="smooth" || family=="jitter"),"invalid case");
  Grid grid(n,family);
  const std::string common="TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
      "MUSCL_FLOW=YES\nMUSCL_KAPPA_FLOW=0\nSLOPE_LIMITER_FLOW=NONE\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\n"
      "FLUID_MODEL=IDEAL_GAS\nGAMMA_VALUE=1.4\nREF_DIMENSIONALIZATION=DIMENSIONAL\n"
      "INIT_OPTION=TD_CONDITIONS\nFREESTREAM_OPTION=DENSITY_FS\nFREESTREAM_DENSITY=1\nFREESTREAM_PRESSURE=1\n"
      "MACH_NUMBER=0.4\nTIME_DOMAIN=YES\nTIME_MARCHING=TIME_STEPPING\nTIME_STEP=0.001\nUNST_CFL_NUMBER=0\n";
  MeshFixture native(grid.xy,grid.triangles,"EULER",common);
  MeshFixture primitive(grid.xy,grid.triangles,"PRIMITIVE_EULER",common+
      "PRIMITIVE_REPRESENTATION=NODE\nPRIMITIVE_FACE_QUADRATURE=AGGREGATED\n"
      "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\nPRIMITIVE_INITIAL_STATE=(1,0.5,0.3,1)\n");
  auto* ng=native.geometry.get(); auto* nc=native.config.get();
  auto* pg=primitive.geometry.get(); auto* pc=primitive.config.get();
  // The primitive API consumes W,x,t directly and does not run CEuler's nondimensionalization.
  // Declare its reference-unit metadata explicitly; no state/source rescaling is performed.
  pc->SetDensity_Ref(1); pc->SetPressure_Ref(1); pc->SetVelocity_Ref(1); pc->SetTime_Ref(1);
  nc->SetGlobalParam(MAIN_SOLVER::EULER,RUNTIME_FLOW_SYS);
  pc->SetGlobalParam(MAIN_SOLVER::PRIMITIVE_EULER,RUNTIME_FLOW_SYS);
  ng->PreprocessP2PComms(ng,nc); ng->PreprocessPeriodicComms(ng,nc);
  Require(ng->GetnPoint()==pg->GetnPoint() && ng->GetnPoint()==grid.xy.size(),"point count mismatch");
  double coordinateError=0,volumeError=0,quadratureValidation=0;
  for (unsigned long i=0;i<grid.xy.size();++i) {
    for (unsigned d=0;d<2;++d) {
      coordinateError=std::max(coordinateError,std::abs(SU2_TYPE::GetValue(ng->nodes->GetCoord(i,d)-pg->nodes->GetCoord(i,d))));
      Require(SU2_TYPE::GetValue(ng->nodes->GetCoord(i,d))==grid.xy[i][d],"node numbering differs from grid");
    }
    volumeError=std::max(volumeError,std::abs(SU2_TYPE::GetValue(ng->nodes->GetVolume(i)-pg->nodes->GetVolume(i))));
  }
  Require(coordinateError==0 && volumeError==0,"mesh identity check failed");
  const auto dual=BuildDual(grid,*ng,quadratureValidation);
  std::unique_ptr<CEulerSolver> conservative;
  std::unique_ptr<CPrimitiveEulerSolver> prim;
  {
    QuietGeometryOutput quiet;
    conservative.reset(new CEulerSolver(ng,nc,0));
    prim.reset(new CPrimitiveEulerSolver(pg,pc,0));
  }
  Require(conservative->VerificationSolution==nullptr,"unexpected native verification solution");
  auto* exact=new Manufactured(nc); conservative->VerificationSolution=exact;
  for (auto* c:{nc,pc}) {
    Require(std::abs(SU2_TYPE::GetValue(c->GetGamma())-M::gamma)<1e-14,"gamma mismatch");
    Require(c->GetTime_Ref()==1 && c->GetPressure_Ref()==1 && c->GetDensity_Ref()==1,"normalization mismatch");
    Require(c->GetKind_TimeIntScheme_Flow()==EULER_EXPLICIT && c->GetMUSCLRampValue()==1,"scheme/ramp mismatch");
  }
  Require(!nc->GetBody_Force() && !nc->GetRotating_Frame() && !nc->GetAxisymmetric() &&
          !nc->GetGravityForce() && !nc->GetVorticityConfinement() && !nc->GetPyCustomSource(),"unexpected physical source");
  CSolver* solvers[MAX_SOLS]{}; solvers[FLOW_SOL]=conservative.get();
  CNumerics* numerics[MAX_TERMS]{};  // No physical-source numerics are active; native MMS needs none.

  constexpr double T=0.1;
  const unsigned steps=static_cast<unsigned>(std::ceil(T*n*n/0.25))*divisor;
  const double dt=T/steps;
  // Disable native screen-error sampling at stale stage time; final norms are computed below at T.
  // With limiter NONE this iteration label does not alter reconstruction or the FE update.
  nc->SetInnerIter(0);pc->SetInnerIter(0);
  nc->SetDelta_UnstTimeND(dt); pc->SetDelta_UnstTimeND(dt);
  auto* nodes=conservative->GetNodes();
  for (unsigned long i=0;i<grid.xy.size();++i) {
    exact->GetInitialCondition(ng->nodes->GetCoord(i),nodes->GetSolution(i));
    prim->GetPrimitiveNodes().SetPrimitiveState(i,M::Evaluate(grid.xy[i][0],grid.xy[i][1],0.0).w);
    nodes->SetDelta_Time(i,dt);
  }
  std::vector<Vec> wp(grid.xy.size()),up(grid.xy.size()),wc(grid.xy.size()),uc(grid.xy.size());
  Vec q0p{},q0c{},linearP{},nonlinearP{},residualC{},injectionC{},qP{},qC{};
  double minRhoP=1e30,minPressureP=1e30,minRhoC=1e30,minPressureC=1e30,minTheta=1;
  double timeError=0,dirichletError=0,updateError=0,balanceMaxC=0,balanceMaxP=0,maxStageCFL=0;
  unsigned long contractions=0;
  const auto sum=[&](const std::vector<Vec>& values) {
    Vec total{};
    for (unsigned long i=0;i<values.size();++i)
      for (unsigned k=0;k<4;++k) total[k]+=SU2_TYPE::GetValue(ng->nodes->GetVolume(i))*values[i][k];
    return total;
  };
  for (unsigned long i=0;i<grid.xy.size();++i) {
    wp[i]=prim->GetPrimitiveNodes().GetPrimitiveState(i); up[i]=M::Conservative(wp[i]);
    for (unsigned k=0;k<4;++k) uc[i][k]=SU2_TYPE::GetValue(nodes->GetSolution(i,k));
    wc[i]=PrimitiveFromU(uc[i]);
    minRhoP=std::min(minRhoP,wp[i][0]);minPressureP=std::min(minPressureP,wp[i][3]);
    minRhoC=std::min(minRhoC,wc[i][0]);minPressureC=std::min(minPressureC,wc[i][3]);
    for (unsigned k=0;k<4;++k) Require(std::abs(uc[i][k]-up[i][k])<2e-14,"initial physical state mismatch");
  }
  q0p=sum(up);q0c=sum(uc);
  unsigned long sourceCallsP=0,boundaryCallsP=0;
  double stageTime=0;
  const auto boundary=[&](const Direction& q,double t,const Vec&,const Direction&,unsigned short marker) {
    Require(marker==0,"unknown marker"); ++boundaryCallsP;
    timeError=std::max(timeError,std::abs(t-stageTime)); return M::Evaluate(q[0],q[1],t).w;
  };
  const auto source=[&](const Direction& q,double t) {
    ++sourceCallsP; timeError=std::max(timeError,std::abs(t-stageTime));
    return M::PrimitiveSource(M::Evaluate(q[0],q[1],t));
  };
  const auto start=std::chrono::steady_clock::now();
  for (unsigned step=0;step<steps;++step) {
    stageTime=step*dt; const double nextTime=(step+1)*dt;
    nc->SetPhysicalTime(stageTime);pc->SetPhysicalTime(stageTime);exact->expectedTime=stageTime;
    nodes->Set_OldSolution();
    conservative->BC_Custom(ng,solvers,nullptr,nullptr,nc,0);
    conservative->Preprocessing(ng,solvers,nc,0,0,RUNTIME_FLOW_SYS,false);
    Require(nc->GetNonphysical_Points()==0,"native nodal fallback");
    conservative->Upwind_Residual(ng,solvers,nullptr,nc,0);
    Require(nc->GetNonphysical_Reconstr()==0,"native edge fallback");
    conservative->Source_Residual(ng,solvers,numerics,nc,0);
    conservative->BC_Custom(ng,solvers,nullptr,nullptr,nc,0);
    // Retain the actual residual used by Forward Euler, with constrained rows zero.
    Vec stepResidual{};
    std::vector<Vec> oldU(grid.xy.size());
    for (unsigned long i=0;i<grid.xy.size();++i) {
      const double volume=SU2_TYPE::GetValue(ng->nodes->GetVolume(i));
      Require(nodes->GetDelta_Time(i)==dt,"nonuniform native physical step");
      for (unsigned k=0;k<4;++k) {
        oldU[i][k]=SU2_TYPE::GetValue(nodes->GetSolution(i,k));
        stepResidual[k]+=dt*SU2_TYPE::GetValue(conservative->LinSysRes(i,k));
      }
    }
    conservative->ExplicitEuler_Iteration(ng,solvers,nc);
    Vec stepInjection{};
    exact->expectedTime=nextTime;nc->SetPhysicalTime(nextTime);
    // Audit the raw FE candidate before native preprocessing could restore an old state.
    for (unsigned long i=0;i<grid.xy.size();++i) {
      Vec raw{};const double volume=SU2_TYPE::GetValue(ng->nodes->GetVolume(i));
      for (unsigned k=0;k<4;++k) {
        raw[k]=SU2_TYPE::GetValue(nodes->GetSolution(i,k));
        const double expected=oldU[i][k]-dt*SU2_TYPE::GetValue(conservative->LinSysRes(i,k))/volume;
        updateError=std::max(updateError,std::abs(raw[k]-expected));
      }
      const auto state=PrimitiveFromU(raw);
      minRhoC=std::min(minRhoC,state[0]);minPressureC=std::min(minPressureC,state[3]);
      if (grid.layer(i)==0) {
        const auto prescribed=M::Conservative(M::Evaluate(grid.xy[i][0],grid.xy[i][1],nextTime).w);
        for (unsigned k=0;k<4;++k) stepInjection[k]+=volume*(prescribed[k]-raw[k]);
      }
    }
    conservative->BC_Custom(ng,solvers,nullptr,nullptr,nc,0);
    // Primitive FE evaluates a fresh full spatial stage at (W^n,t_n).
    const auto oldWp=wp;const auto oldUp=up;
    const auto result=prim->AdvanceExplicitEuler(pg,pc,stageTime,dt,boundary,source);
    timeError=std::max(timeError,std::abs(prim->GetPhysicalTime()-nextTime));
    pc->SetPhysicalTime(nextTime);
    for (unsigned long i=0;i<grid.xy.size();++i) {
      const double volume=SU2_TYPE::GetValue(ng->nodes->GetVolume(i));
      wp[i]=prim->GetPrimitiveNodes().GetPrimitiveState(i);up[i]=M::Conservative(wp[i]);
      const double theta=prim->GetPositivityScale(i);
      minTheta=std::min(minTheta,theta);contractions+=theta<1;
      maxStageCFL=std::max(maxStageCFL,dt*result.spectralRadius[i]/volume);
      minRhoP=std::min(minRhoP,wp[i][0]);minPressureP=std::min(minPressureP,wp[i][3]);
      const auto mapped=M::Apply(M::Jacobian(oldWp[i]),result.rhs[i]);
      for (unsigned k=0;k<4;++k) {
        linearP[k]+=volume*dt*mapped[k];
        nonlinearP[k]+=volume*(up[i][k]-oldUp[i][k]-dt*mapped[k]);
        uc[i][k]=SU2_TYPE::GetValue(nodes->GetSolution(i,k));
      }
      wc[i]=PrimitiveFromU(uc[i]);
      minRhoC=std::min(minRhoC,wc[i][0]);minPressureC=std::min(minPressureC,wc[i][3]);
      if (grid.layer(i)==0) {
        const auto prescribed=M::Conservative(M::Evaluate(grid.xy[i][0],grid.xy[i][1],nextTime).w);
        for (unsigned k=0;k<4;++k) dirichletError=std::max(dirichletError,std::abs(uc[i][k]-prescribed[k]));
      }
    }
    qP=sum(up);qC=sum(uc);
    for (unsigned k=0;k<4;++k) {
      residualC[k]+=stepResidual[k];injectionC[k]+=stepInjection[k];
      balanceMaxC=std::max(balanceMaxC,std::abs(qC[k]-q0c[k]+residualC[k]-injectionC[k]));
      balanceMaxP=std::max(balanceMaxP,std::abs(qP[k]-q0p[k]-linearP[k]-nonlinearP[k]));
    }
    if (step%std::max(1u,steps/4)==0)
      std::cout << family << " N=" << n << " divisor=" << divisor << " step=" << step << '/' << steps << std::endl;
  }
  Require(timeError<2e-14 && exact->sourceTimeError==0 && exact->boundaryTimeError==0,"stage time mismatch");
  Require(exact->sourceCalls==steps*grid.xy.size() && sourceCallsP==steps*grid.xy.size(),"source stage count mismatch");
  Require(exact->boundaryCalls==3*steps*ng->GetnVertex(0),"native boundary stage count mismatch");
  Require(boundaryCallsP==steps*prim->GetGeometryMetrics().GetBoundaryFaces().size(),"primitive boundary stage count mismatch");
  Require(dirichletError<2e-14 && updateError<2e-14,"Dirichlet/FE update mismatch");
  Require(balanceMaxC<2e-11 && balanceMaxP<2e-11,"discrete ledger mismatch");
  Require(std::abs(prim->GetPhysicalTime()-T)<2e-14 && std::abs(SU2_TYPE::GetValue(nc->GetPhysicalTime())-T)<2e-14,"wrong final physical time");
  std::ofstream out(argv[4]);out << std::setprecision(17);
  out << "kind,family,N,divisor,time,group,component,count,volume,rms,l1,l2,linf\n";
  for (const std::string kind:{"W_primitive","U_primitive","W_conservative","U_conservative","difference_W","difference_U"}) {
    for (unsigned g=0;g<7;++g) for (unsigned k=0;k<4;++k) {
      Norm norm;
      for (unsigned long i=0;i<grid.xy.size();++i) if (grid.in(i,g)) {
        const auto w=M::Evaluate(grid.xy[i][0],grid.xy[i][1],T).w,u=M::Conservative(w);
        const double error=kind=="W_primitive"?wp[i][k]-w[k]:kind=="U_primitive"?up[i][k]-u[k]:
            kind=="W_conservative"?wc[i][k]-w[k]:kind=="U_conservative"?uc[i][k]-u[k]:
            kind=="difference_W"?wp[i][k]-wc[i][k]:up[i][k]-uc[i][k];
        norm.add(error,SU2_TYPE::GetValue(ng->nodes->GetVolume(i)));
      }
      out << kind << ',' << family << ',' << n << ',' << divisor << ',' << T << ',' << groupNames[g] << ',' << k
          << ',' << norm.count << ',' << norm.weight << ',' << std::sqrt(norm.square/norm.count) << ','
          << norm.absolute/norm.weight << ',' << std::sqrt(norm.weighted/norm.weight) << ',' << norm.maximum << '\n';
    }
  }
  std::ofstream checks(std::string(argv[4])+".checks");
  checks << std::setprecision(17) << "family,N,divisor,kind,component,value\n";
  const auto check=[&](const char* kind,double value,unsigned k=0) {
    checks << family << ',' << n << ',' << divisor << ',' << kind << ',' << k << ',' << value << '\n';
  };
  check("dt",dt);check("steps",steps);check("T",prim->GetPhysicalTime());check("stage_time_error",timeError);
  check("rho_min_primitive",minRhoP);check("p_min_primitive",minPressureP);check("rho_min_conservative",minRhoC);check("p_min_conservative",minPressureC);
  check("theta_min",minTheta);check("contractions",contractions);check("native_fallbacks",0);
  check("native_update_error",updateError);check("Dirichlet_error",dirichletError);check("balance_max_C",balanceMaxC);check("balance_max_P",balanceMaxP);
  check("max_stage_CFL_primitive",maxStageCFL);check("coordinate_error",coordinateError);check("volume_error",volumeError);
  check("source_calls_primitive",sourceCallsP);check("source_calls_conservative",exact->sourceCalls);
  check("runtime_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
  Vec analytic{};
  for (unsigned long i=0;i<grid.xy.size();++i) {
    const auto u=M::Conservative(M::Evaluate(grid.xy[i][0],grid.xy[i][1],T).w);
    for (unsigned k=0;k<4;++k) analytic[k]+=SU2_TYPE::GetValue(ng->nodes->GetVolume(i))*u[k];
  }
  for (unsigned k=0;k<4;++k) {
    check("Q_primitive_minus_exact_nodal",qP[k]-analytic[k],k);check("Q_conservative_minus_exact_nodal",qC[k]-analytic[k],k);
    check("native_constraint_injection",injectionC[k],k);check("native_residual_integral",residualC[k],k);
    check("primitive_linear_increment",linearP[k],k);check("primitive_nonlinear_increment",nonlinearP[k],k);
  }
  std::ofstream field(std::string(argv[4])+".field");field << std::setprecision(17);
  for (unsigned long i=0;i<grid.xy.size();++i) {
    field << i << ' ' << SU2_TYPE::GetValue(ng->nodes->GetVolume(i)) << ' ' << grid.xy[i][0] << ' ' << grid.xy[i][1];
    for (const auto* values:{&wp,&up,&wc,&uc}) for (auto value:(*values)[i]) field << ' ' << value;
    field << '\n';
  }
  return 0;
}
int main(int argc,char** argv) {
  SU2_MPI::Init(&argc,&argv);int status=0;
  try {status=Run(argc,argv);} catch(const std::exception& e) {std::cerr << "MMS EVOLUTION FAILED: " << e.what() << '\n';status=1;}
  SU2_MPI::Finalize();return status;
}
