/*!
 * \brief Phase-6B.1 semidiscrete comparison of two production Euler solvers.
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

int Run(int argc,char** argv) {
  Require(argc==4,"usage: euler_mms_audit regular|smooth|jitter N output.csv");
  Require(SU2_MPI::GetSize()==1,"one MPI rank is required");
  const std::string family=argv[1]; const unsigned n=std::stoul(argv[2]);
  Require(n>=8 && (family=="regular" || family=="smooth" || family=="jitter"),"invalid mesh family/size");
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
  std::ofstream out(argv[3]); out << std::setprecision(17);
  out << "kind,family,N,time,group,component,count,volume,rms,l1,l2,linf\n";
  const auto diagnostic=[&](const char* label,double t,double value) {
    Require(std::isfinite(value),"nonfinite check");
    out << label << ',' << family << ',' << n << ',' << t << ",check,0,0,0," << value << ",0,0,0\n";
  };
  diagnostic("coordinate_error",0,coordinateError); diagnostic("volume_error",0,volumeError);
  diagnostic("quadrature_polynomial",0,quadratureValidation);
  diagnostic("gamma_conservative",0,SU2_TYPE::GetValue(nc->GetGamma()));
  diagnostic("gamma_primitive",0,SU2_TYPE::GetValue(pc->GetGamma()));
  diagnostic("rho_ref_conservative",0,SU2_TYPE::GetValue(nc->GetDensity_Ref()));
  diagnostic("rho_ref_primitive",0,SU2_TYPE::GetValue(pc->GetDensity_Ref()));
  diagnostic("p_ref_conservative",0,SU2_TYPE::GetValue(nc->GetPressure_Ref()));
  diagnostic("p_ref_primitive",0,SU2_TYPE::GetValue(pc->GetPressure_Ref()));
  diagnostic("time_ref_conservative",0,SU2_TYPE::GetValue(nc->GetTime_Ref()));
  diagnostic("time_ref_primitive",0,SU2_TYPE::GetValue(pc->GetTime_Ref()));
  for (double time:{0.0,0.05,0.1}) {
    nc->SetPhysicalTime(time); pc->SetPhysicalTime(time); exact->expectedTime=time;
    exact->sourceCalls=exact->boundaryCalls=0; exact->sourceTimeError=exact->boundaryTimeError=0;
    double initU=0,initW=0,metricVolume=0;
    for (unsigned long i=0;i<grid.xy.size();++i) {
      const auto s=M::Evaluate(grid.xy[i][0],grid.xy[i][1],time);
      const auto u=M::Conservative(s.w);
      prim->GetPrimitiveNodes().SetPrimitiveState(i,s.w);
      auto* node=conservative->GetNodes();
      if (time==0) exact->GetInitialCondition(ng->nodes->GetCoord(i),node->GetSolution(i));
      else exact->GetSolution(ng->nodes->GetCoord(i),time,node->GetSolution(i));
      metricVolume=std::max(metricVolume,std::abs(prim->GetGeometryMetrics().GetCell(i).volume-SU2_TYPE::GetValue(ng->nodes->GetVolume(i))));
      for (unsigned k=0;k<4;++k) {
        initU=std::max(initU,std::abs(SU2_TYPE::GetValue(node->GetSolution(i,k))-u[k]));
        initW=std::max(initW,std::abs(prim->GetPrimitiveNodes().GetPrimitiveState(i)[k]-s.w[k]));
      }
    }
    Require(initU<1e-15 && initW==0 && metricVolume==0,"initial state or metric mismatch");
    conservative->GetNodes()->Set_OldSolution();
    conservative->BC_Custom(ng,solvers,nullptr,nullptr,nc,0);
    conservative->Preprocessing(ng,solvers,nc,0,0,RUNTIME_FLOW_SYS,false);
    Require(nc->GetNonphysical_Points()==0,"undeclared nonphysical nodal fallback");
    double recoveredPrimitiveError=0;
    for (unsigned long i=0;i<grid.xy.size();++i) {
      const auto w=M::Evaluate(grid.xy[i][0],grid.xy[i][1],time).w;
      const auto* node=conservative->GetNodes();
      const Vec actual{{SU2_TYPE::GetValue(node->GetDensity(i)),SU2_TYPE::GetValue(node->GetVelocity(i,0)),
                        SU2_TYPE::GetValue(node->GetVelocity(i,1)),SU2_TYPE::GetValue(node->GetPressure(i))}};
      for (unsigned k=0;k<4;++k) recoveredPrimitiveError=std::max(recoveredPrimitiveError,std::abs(actual[k]-w[k]));
    }
    Require(recoveredPrimitiveError<3e-14,"conservative physical state differs from primitive reference");
    conservative->Upwind_Residual(ng,solvers,nullptr,nc,0);
    Require(nc->GetNonphysical_Reconstr()==0,"undeclared nonphysical edge fallback");
    std::vector<Vec> flux(grid.xy.size());
    for (unsigned long i=0;i<grid.xy.size();++i) for (unsigned k=0;k<4;++k)
      flux[i][k]=SU2_TYPE::GetValue(conservative->LinSysRes(i,k));
    conservative->Source_Residual(ng,solvers,numerics,nc,0);
    Require(exact->sourceCalls==grid.xy.size() && exact->sourceTimeError==0,"native source count/time mismatch");
    double sourceIntegrated=0,sourceRecovered=0;
    std::vector<Vec> tauU(grid.xy.size());
    for (unsigned long i=0;i<grid.xy.size();++i) {
      const auto s=M::Evaluate(grid.xy[i][0],grid.xy[i][1],time);
      const auto source=M::ConservativeSource(s),ut=M::Apply(M::Jacobian(s.w),s.dt);
      const double volume=SU2_TYPE::GetValue(ng->nodes->GetVolume(i));
      for (unsigned k=0;k<4;++k) {
        const double r=SU2_TYPE::GetValue(conservative->LinSysRes(i,k));
        sourceIntegrated=std::max(sourceIntegrated,std::abs(r-flux[i][k]+volume*source[k]));
        sourceRecovered=std::max(sourceRecovered,std::abs((flux[i][k]-r)/volume-source[k]));
        tauU[i][k]=-r/volume-ut[k]; // Only free nodes are sampled below; no boundary defect is inferred here.
      }
    }
    Require(sourceIntegrated<2e-14 && sourceRecovered<5e-11,"native source sign/volume/time mismatch");
    conservative->BC_Custom(ng,solvers,nullptr,nullptr,nc,0);
    Require(exact->boundaryCalls==2*ng->GetnVertex(0) && exact->boundaryTimeError==0,"native boundary count/time mismatch");
    unsigned long primitiveSourceCalls=0,primitiveBoundaryCalls=0;
    double callbackTimeError=0;
    const auto result=prim->ComputeSpatialResidual(pg,pc,time,
        [&](const Direction& q,double t,const Vec&,const Direction&,unsigned short marker) {
          ++primitiveBoundaryCalls; callbackTimeError=std::max(callbackTimeError,std::abs(t-time));
          Require(marker==0,"unexpected boundary marker"); return M::Evaluate(q[0],q[1],t).w;
        },[&](const Direction& q,double t) {
          ++primitiveSourceCalls; callbackTimeError=std::max(callbackTimeError,std::abs(t-time));
          return M::PrimitiveSource(M::Evaluate(q[0],q[1],t));
        });
    Require(primitiveSourceCalls==grid.xy.size() && primitiveBoundaryCalls==prim->GetGeometryMetrics().GetBoundaryFaces().size()
            && callbackTimeError==0,"primitive callback count/time mismatch");
    double primitiveSourceError=0,pressureGradient=0,quadrature=0,nonlinearGap=0,chainError=0;
    std::vector<Vec> tauW(grid.xy.size()),jTauW(grid.xy.size()),difference(grid.xy.size()),constraint(grid.xy.size()),gapU(grid.xy.size()),gapW(grid.xy.size());
    for (unsigned long i=0;i<grid.xy.size();++i) {
      const auto s=M::Evaluate(grid.xy[i][0],grid.xy[i][1],time);
      const auto source=M::PrimitiveSource(s),u=M::Conservative(s.w);
      Require(prim->GetPositivityScale(i)==1,"unexpected positivity contraction");
      pressureGradient=std::max(pressureGradient,std::hypot(s.dx[3],s.dy[3]));
      const auto meanU=Integrate(dual[i],8,[&](const Direction& q){return M::Conservative(M::Evaluate(q[0],q[1],time).w);});
      const auto meanW=Integrate(dual[i],8,[&](const Direction& q){return M::Evaluate(q[0],q[1],time).w;});
      const auto meanU5=Integrate(dual[i],5,[&](const Direction& q){return M::Conservative(M::Evaluate(q[0],q[1],time).w);});
      const auto meanW5=Integrate(dual[i],5,[&](const Direction& q){return M::Evaluate(q[0],q[1],time).w;});
      const auto converted=M::Conservative(meanW);
      for (unsigned k=0;k<4;++k) {
        tauW[i][k]=result.rhs[i][k]-s.dt[k];
        primitiveSourceError=std::max(primitiveSourceError,std::abs(result.source[i][k]-source[k]));
        quadrature=std::max({quadrature,std::abs(meanU[k]-meanU5[k]),std::abs(meanW[k]-meanW5[k])});
        nonlinearGap=std::max(nonlinearGap,std::abs(meanU[k]-converted[k]));
        gapU[i][k]=meanU[k]-u[k]; gapW[i][k]=meanW[k]-s.w[k];
        if (grid.layer(i)==0) {
          constraint[i][k]=SU2_TYPE::GetValue(conservative->GetNodes()->GetSolution(i,k))-u[k];
          Require(std::abs(constraint[i][k])<2e-14,"Dirichlet state differs beyond roundoff");
          Require(conservative->LinSysRes(i,k)==0,"Dirichlet residual not zero");
        }
      }
      jTauW[i]=M::Apply(M::Jacobian(s.w),tauW[i]);
      const auto mappedRhs=M::Apply(M::Jacobian(s.w),result.rhs[i]);
      const auto ut=M::Apply(M::Jacobian(s.w),s.dt);
      for (unsigned k=0;k<4;++k) {
        difference[i][k]=jTauW[i][k]-tauU[i][k];
        chainError=std::max(chainError,std::abs(jTauW[i][k]-(mappedRhs[k]-ut[k])));
      }
    }
    std::cout << "CHECK,time=" << time << ",primitive_source=" << primitiveSourceError
              << ",quadrature=" << quadrature << ",chain=" << chainError << '\n';
    Require(primitiveSourceError<2e-14,"primitive source differs beyond roundoff");
    Require(quadrature<1e-10,"analytic mean quadrature unresolved");
    Require(chainError<1e-14,"chain identity failed");
    Require(time==0?pressureGradient==0:pressureGradient>0.01,"pressure gradient signal missing");
    diagnostic("init_U_error",time,initU); diagnostic("init_W_error",time,initW);
    diagnostic("native_recovered_W_error",time,recoveredPrimitiveError);
    diagnostic("native_source_integrated_error",time,sourceIntegrated); diagnostic("native_source_recovered_error",time,sourceRecovered);
    diagnostic("native_source_time_error",time,exact->sourceTimeError);
    diagnostic("native_boundary_time_error",time,exact->boundaryTimeError);
    diagnostic("native_source_calls",time,exact->sourceCalls); diagnostic("native_boundary_calls",time,exact->boundaryCalls);
    diagnostic("primitive_source_calls",time,primitiveSourceCalls); diagnostic("primitive_boundary_calls",time,primitiveBoundaryCalls);
    diagnostic("primitive_source_error",time,primitiveSourceError); diagnostic("callback_time_error",time,callbackTimeError);
    diagnostic("quadrature_5_vs_8",time,quadrature); diagnostic("meanU_minus_UmeanW",time,nonlinearGap);
    diagnostic("pressure_gradient_max",time,pressureGradient); diagnostic("chain_identity_error",time,chainError);
    for (const auto& item:std::vector<std::pair<std::string,const std::vector<Vec>*>>{
        {"tau_W",&tauW},{"J_tau_W",&jTauW},{"tau_U",&tauU},{"difference_U",&difference},
        {"Dirichlet_constraint",&constraint},{"mean_point_U",&gapU},{"mean_point_W",&gapW}}) {
      const bool freeOnly=item.first=="tau_U" || item.first=="difference_U" || item.first=="J_tau_W";
      for (unsigned g=0;g<7;++g) {
        if ((freeOnly && g==1) || (item.first=="Dirichlet_constraint" && g!=1)) continue;
        for (unsigned k=0;k<4;++k) {
          Norm norm;
          for (unsigned long i=0;i<grid.xy.size();++i) if (grid.in(i,g) && (!freeOnly || grid.layer(i)>0))
            norm.add((*item.second)[i][k],SU2_TYPE::GetValue(ng->nodes->GetVolume(i)));
          Require(norm.count>0,"empty diagnostic region");
          const std::string group=g==0 && freeOnly?"all_free":groupNames[g];
          out << item.first << ',' << family << ',' << n << ',' << time << ',' << group << ',' << k << ',' << norm.count << ',' << norm.weight
              << ',' << std::sqrt(norm.square/norm.count) << ',' << norm.absolute/norm.weight << ',' << std::sqrt(norm.weighted/norm.weight)
              << ',' << norm.maximum << '\n';
        }
      }
    }
  }
  return 0;
}
int main(int argc,char** argv) {
  SU2_MPI::Init(&argc,&argv); int status=0;
  try { status=Run(argc,argv); } catch (const std::exception& e) { std::cerr << "MMS AUDIT FAILED: " << e.what() << '\n'; status=1; }
  SU2_MPI::Finalize(); return status;
}
