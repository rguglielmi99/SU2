/*!
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "../../include/iteration/CPrimitiveEulerIteration.hpp"
#include "../../include/iteration/primitive_euler_boundary_config.hpp"
#include "../../include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../include/output/COutput.hpp"
#include "../../../Common/include/verification/isentropic_vortex.hpp"
#include <set>

void CPrimitiveEulerIteration::ValidateDriver(const CConfig& c) {
  if (!c.GetTime_Domain() || c.GetTime_Marching()!=TIME_MARCHING::TIME_STEPPING)
    throw std::invalid_argument("PRIMITIVE_EULER driver requires TIME_DOMAIN=YES and TIME_MARCHING=TIME_STEPPING.");
  if (SU2_MPI::GetSize()!=1 || omp_get_max_threads()!=1)
    throw std::invalid_argument("PRIMITIVE_EULER supports only serial mode: one MPI rank and one OpenMP thread.");
  CPrimitiveEulerSolver::ValidateConfiguration(c);
  if (c.GetSystemMeasurements()!=SI)
    throw std::invalid_argument("PRIMITIVE_EULER physical boundary adapter currently requires SYSTEM_MEASUREMENTS=SI");
  if (c.GetnMarker_Periodic()!=0)
    throw std::invalid_argument("PRIMITIVE_EULER periodic geometry/residual coupling is not implemented");
  if (c.GetVerification_Solution()!=VERIFICATION_SOLUTION::NONE)
    throw std::invalid_argument("PRIMITIVE_EULER driver does not support native KIND_VERIFICATION_SOLUTION/MMS; use PRIMITIVE_INITIALIZATION for its analytic benchmark.");
  if (c.GetGrid_Movement() || c.GetDeform_Mesh() || c.GetDirectDiff()!=NO_DERIVATIVE ||
      c.GetContinuous_Adjoint() || c.GetDiscrete_Adjoint() || c.GetWind_Gust() || c.GetFixed_CL_Mode() ||
      c.GetPredictor() || c.GetRelaxation() || c.GetWnd_Cauchy_Crit() ||
      c.GetBoolTurbomachinery() || c.GetBody_Force() || c.GetRotating_Frame() || c.GetGravityForce() ||
      c.GetVorticityConfinement() || c.GetPyCustomSource() || c.GetInlet_Profile_From_File())
    throw std::invalid_argument("PRIMITIVE_EULER supports static meshes without motion, differentiation or additional sources/models.");
  if (c.GetRef_NonDim()!=DIMENSIONAL)
    throw std::invalid_argument("PRIMITIVE_EULER driver requires REF_DIMENSIONALIZATION=DIMENSIONAL (unit reference scales).");
  for (unsigned short file=0;file<c.GetnVolumeOutputFiles();++file)
    if (c.GetVolumeOutputFrequency(file)==0)
      throw std::invalid_argument("PRIMITIVE_EULER OUTPUT_WRT_FREQ entries must be positive.");
  if (!std::isfinite(SU2_TYPE::GetValue(c.GetGas_Constant())) || c.GetGas_Constant()<=0)
    throw std::invalid_argument("PRIMITIVE_EULER requires a finite positive GAS_CONSTANT for derived temperature output.");
  if (c.GetKind_ConvNumScheme_Flow()!=SPACE_UPWIND || c.GetKind_Upwind_Flow()!=UPWIND::ROE)
    throw std::invalid_argument("PRIMITIVE_EULER requires CONV_NUM_METHOD_FLOW=ROE as a selection token; it uses its original primitive path kernel.");
  for (const auto* option:{"ROE_KAPPA","ENTROPY_FIX_COEFF","ROE_LOW_DISSIPATION","LOW_MACH_CORR","LOW_MACH_PREC"})
    if (c.OptionIsSet(option))
      throw std::invalid_argument(std::string("PRIMITIVE_EULER does not use native conservative flux option ")+option);
  const auto& initial=c.GetPrimitiveInitialization();
  const auto& boundary=c.GetPrimitiveBoundaryMode();
  if ((initial!="UNIFORM" && initial!="FREESTREAM" && initial!="ISENTROPIC_VORTEX") ||
      (boundary!="UNIFORM" && boundary!="ISENTROPIC_VORTEX" && boundary!="LOCAL_CHARACTERISTIC"))
    throw std::invalid_argument("PRIMITIVE_INITIALIZATION supports UNIFORM, FREESTREAM or ISENTROPIC_VORTEX; boundary mode supports UNIFORM, ISENTROPIC_VORTEX or LOCAL_CHARACTERISTIC.");
  if ((initial=="ISENTROPIC_VORTEX") != (boundary=="ISENTROPIC_VORTEX"))
    throw std::invalid_argument("The vortex benchmark requires matching analytic initialization and exterior boundary mode.");
  if (initial=="ISENTROPIC_VORTEX" && std::abs(SU2_TYPE::GetValue(c.GetGamma())-1.4)>1e-14)
    throw std::invalid_argument("The phase-7 isentropic vortex benchmark requires gamma=1.4.");
  PrimitiveEuler2D::State w{};
  for (unsigned k=0;k<4;++k) w[k]=SU2_TYPE::GetValue(c.GetPrimitiveBoundaryState()[k]);
  PrimitiveEuler2D::Validate(w,SU2_TYPE::GetValue(c.GetGamma()));
}

void CPrimitiveEulerIteration::Solve(COutput* output, CIntegration****, CGeometry**** geometry,
    CSolver***** solver, CNumerics******, CConfig** config, CSurfaceMovement**, CVolumetricMovement***,
    CFreeFormDefBox***, unsigned short zone, unsigned short instance) {
  auto* c=config[zone];
  auto* g=geometry[zone][instance][MESH_0];
  auto** container=solver[zone][instance][MESH_0];
  auto* s=dynamic_cast<CPrimitiveEulerSolver*>(container[FLOW_SOL]);
  try {
    if (!s) throw std::logic_error("PRIMITIVE_EULER factory did not instantiate CPrimitiveEulerSolver.");
    if (!initialized) {
      std::set<std::string> actual,declared;
      for (unsigned short m=0;m<g->GetnMarker();++m) {
        const auto kind=c->GetMarker_All_KindBC(m);
        if (!PrimitiveEulerBoundaryConfig::Supported(kind))
          throw std::invalid_argument("PRIMITIVE_EULER unsupported physical boundary marker: "+c->GetMarker_All_TagBound(m));
        if (kind==FAR_FIELD) {
          if ((c->OptionIsSet("PRIMITIVE_BOUNDARY_MODE") && c->GetPrimitiveBoundaryMode()!="LOCAL_CHARACTERISTIC") || c->OptionIsSet("PRIMITIVE_BOUNDARY_STATE"))
            throw std::invalid_argument("MARKER_FAR requires native physical freestream data and explicit LOCAL_CHARACTERISTIC mode; move prescribed/analytic exterior data to MARKER_CUSTOM, not PRIMITIVE_BOUNDARY_STATE on FAR.");
          if (c->GetPrimitiveBoundaryMode()!="LOCAL_CHARACTERISTIC")
            throw std::invalid_argument("Primitive finite-amplitude farfield compatibility remains under review; explicitly select PRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC for the local linearized reservoir model.");
          PrimitiveEulerBoundaryConfig::Freestream(*c);
        }
        if (!actual.insert(c->GetMarker_All_TagBound(m)).second)
          throw std::invalid_argument("PRIMITIVE_EULER does not accept duplicate physical marker tags.");
      }
      for (unsigned short m=0;m<c->GetnMarker_CfgFile();++m)
        if (!declared.insert(c->GetMarker_CfgFile_TagBound(m)).second)
          throw std::invalid_argument("PRIMITIVE_EULER does not accept duplicate configured marker tags.");
      if (actual!=declared)
        throw std::invalid_argument("PRIMITIVE_EULER mesh/config boundary marker sets differ (missing or unknown marker).");
      c->SetDensity_Ref(1); c->SetPressure_Ref(1); c->SetVelocity_Ref(1); c->SetTime_Ref(1);
      c->SetGlobalParam(MAIN_SOLVER::PRIMITIVE_EULER,RUNTIME_FLOW_SYS);
      if (c->GetPrimitiveInitialization()=="ISENTROPIC_VORTEX")
        for (unsigned long i=0;i<g->GetnPoint();++i) {
          const auto* x=g->nodes->GetCoord(i);
          s->GetPrimitiveNodes().SetPrimitiveState(i,SU2IsentropicVortex::Evaluate(SU2_TYPE::GetValue(x[0]),SU2_TYPE::GetValue(x[1]),0.).w);
        }
      if (c->GetPrimitiveInitialization()=="FREESTREAM") {
        const auto state=PrimitiveEulerBoundaryConfig::Freestream(*c);
        for (unsigned long i=0;i<g->GetnPoint();++i) s->GetPrimitiveNodes().SetPrimitiveState(i,state);
      }
      c->SetPhysicalTime(0);
      c->SetDelta_UnstTimeND(c->GetDelta_UnstTime());
      output->SetHistoryOutput(g,container,c,0,0,0);
      output->SetResultFiles(g,c,container,0,true);
      initialized=true;
    }
    const double time=s->GetPhysicalTime();
    c->SetPhysicalTime(time);
    const CPrimitiveEulerSolver::BoundaryState boundary=[=](const CPrimitiveEulerSolver::Direction& q,double t,
        const CPrimitiveEulerSolver::State& inside,const CPrimitiveEulerSolver::Direction& n,unsigned short marker) {
      return PrimitiveEulerBoundaryConfig::Evaluate(*c,marker,q,t,inside,n);
    };
    double dt=SU2_TYPE::GetValue(c->GetDelta_UnstTime());
    if (c->GetUnst_CFL()>0) {
      dt=s->EstimateTimeStep(g,c,time,SU2_TYPE::GetValue(c->GetUnst_CFL()),boundary);
      dt=std::min(dt,SU2_TYPE::GetValue(c->GetMax_Time())-time);
    }
    // Single source of truth: no generic residual/BC/integration route runs here.
    s->AdvanceExplicitEuler(g,c,time,dt,boundary);
    c->SetDelta_UnstTimeND(dt);
    c->SetPhysicalTime(s->GetPhysicalTime());
    c->SetTimeIter(c->GetTimeIter()+1);  // Output index counts accepted steps; zero belongs to the initial state.
    c->SetInnerIter(0);
    output->SetHistoryOutput(g,container,c,c->GetTimeIter(),0,0);
  } catch (const std::exception& e) {
    SU2_MPI::Error(e.what(),CURRENT_FUNCTION);
  }
}
