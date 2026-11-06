/*!
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "../../include/output/CPrimitiveEulerOutput.hpp"
#include "../../include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../include/iteration/primitive_euler_boundary_config.hpp"
#include "../../include/iteration/primitive_euler_weak_te_config.hpp"
#include <fstream>
#include <iomanip>

namespace {
const std::array<std::string,4> key{{"RHO","U","V","P"}};
const std::array<std::string,4> label{{"Density","Velocity_x","Velocity_y","Pressure"}};
CPrimitiveEulerSolver::SpatialResult CurrentResidual(CPrimitiveEulerSolver& s,CGeometry* g,CConfig* c) {
  return s.ComputeSpatialResidual(g,c,s.GetPhysicalTime(),
      [c](const PrimitiveEuler2D::Direction& q,double t,const PrimitiveEuler2D::State& w,
          const PrimitiveEuler2D::Direction& n,unsigned short marker) {
        return PrimitiveEulerBoundaryConfig::Evaluate(*c,marker,q,t,w,n);
      });
}
}
CPrimitiveEulerOutput::CPrimitiveEulerOutput(CConfig* c,unsigned short dim) : CFVMOutput(c,dim,false) {
  if (nRequestedHistoryFields==0) {
    RequestCommonHistory(!c->GetPrimitiveLocalPseudoTime());
    requestedHistoryFields.emplace_back("PRIMITIVE_RHS");
    if (c->GetnMarker_Monitoring()) requestedHistoryFields.emplace_back("AERO_COEFF");
    if (c->GetPrimitiveWeakTEGlobalID()>=0 || c->GetPrimitiveLocalPseudoTime()) requestedHistoryFields.emplace_back("PRIMITIVE_STEADY_AUDIT");
    if (c->GetPrimitiveLocalPseudoTime()) requestedHistoryFields.emplace_back("LOCAL_PSEUDO_TIME");
    nRequestedHistoryFields=requestedHistoryFields.size();
  }
  if (nRequestedScreenFields==0) {
    requestedScreenFields={"TIME_ITER","CUR_TIME","RHS_RHO","RHS_P"};
    if (c->GetPrimitiveLocalPseudoTime()) requestedScreenFields={"TIME_ITER","RHS_CURRENT","DT_PSEUDO_MIN","EFFECTIVE_CFL"};
    nRequestedScreenFields=requestedScreenFields.size();
  }
  if (nRequestedVolumeFields==0) {
    requestedVolumeFields={"COORDINATES","SOLUTION","PRIMITIVE"};
    nRequestedVolumeFields=requestedVolumeFields.size();
  }
  if (convFields.empty()) convFields.emplace_back("RHS_RHO");
  volumeFilename=c->GetVolume_FileName();
  surfaceFilename=c->GetSurfCoeff_FileName();
  restartFilename=c->GetRestart_FileName();
  multiZoneHeaderString="Primitive Euler";
}
void CPrimitiveEulerOutput::SetHistoryOutputFields(CConfig* c) {
  if (c->GetPrimitiveLocalPseudoTime())
    for (const auto& field:{"DT_PSEUDO_MIN","DT_PSEUDO_MAX","DT_PSEUDO_RATIO","DT_PSEUDO_TE","EFFECTIVE_CFL","CFL_REDUCTIONS","PSEUDO_REJECTED","PSEUDO_REJECTED_ATTEMPTS","RHS_WALL_STRONG","RHS_FARFIELD"})
      AddHistoryOutput(field,field,ScreenOutputFormat::SCIENTIFIC,"LOCAL_PSEUDO_TIME","Local stationary iteration diagnostic; not physical time");
  if (c->GetPrimitiveWeakTEGlobalID()>=0 || c->GetPrimitiveLocalPseudoTime())
    for (const auto& field:{"RHS_CURRENT","RHS_INTERIOR","RHS_TE","MAX_MACH","MAX_TRACE_MACH"})
      AddHistoryOutput(field,field,ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_STEADY_AUDIT",
                       "Current-state weak-TE run diagnostic (linear norm or maximum)");
  if (c->GetnMarker_Monitoring()) {
    AddHistoryOutput("LIFT","CL",ScreenOutputFormat::SCIENTIFIC,"AERO_COEFF","Primitive reconstructed pressure lift",HistoryFieldType::COEFFICIENT);
    AddHistoryOutput("DRAG","CD",ScreenOutputFormat::SCIENTIFIC,"AERO_COEFF","Primitive reconstructed pressure drag",HistoryFieldType::COEFFICIENT);
  }
  for (unsigned k=0;k<4;++k)
    AddHistoryOutput("RHS_"+key[k],"L2_rhs_"+key[k],ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_RHS",
                     (c->GetPrimitiveWeakTEGlobalID()>=0 || c->GetPrimitiveLocalPseudoTime()) ? "Volume-weighted L2 of current-state primitive RHS (linear, all nodes)" :
                     "Volume-weighted L2 of accepted-stage primitive RHS (linear, all nodes)");
  AddHistoryOutput("MIN_RHO","Min_rho",ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_STATE","Minimum updated nodal density");
  AddHistoryOutput("MIN_P","Min_p",ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_STATE","Minimum updated nodal pressure");
  AddHistoryOutput("MIN_THETA","Min_theta",ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_STATE","Minimum reconstruction positivity scale");
}
void CPrimitiveEulerOutput::LoadHistoryData(CConfig* c,CGeometry* g,CSolver** solvers) {
  auto* s=static_cast<CPrimitiveEulerSolver*>(solvers[FLOW_SOL]);
  if (c->GetPrimitiveLocalPseudoTime()) {
    const auto& info=s->GetPseudoDiagnostics();
    SetHistoryOutputValue("DT_PSEUDO_MIN",info.minimum); SetHistoryOutputValue("DT_PSEUDO_MAX",info.maximum);
    SetHistoryOutputValue("DT_PSEUDO_RATIO",info.ratio); SetHistoryOutputValue("EFFECTIVE_CFL",info.effectiveCfl);
    SetHistoryOutputValue("CFL_REDUCTIONS",info.reductions); SetHistoryOutputValue("PSEUDO_REJECTED",info.rejectedIterations);
    SetHistoryOutputValue("PSEUDO_REJECTED_ATTEMPTS",info.rejectedAttempts);
    double teStep=0;
    if (!info.steps.empty() && c->GetPrimitiveWeakTEGlobalID()>=0)
      teStep=info.steps.at(PrimitiveEulerWeakTEConfig::LocalNodes(*g,*c).front());
    SetHistoryOutputValue("DT_PSEUDO_TE",teStep);
  }
  if (c->GetPrimitiveWeakTEGlobalID()>=0 || c->GetPrimitiveLocalPseudoTime()) {
    const auto current=CurrentResidual(*s,g,c);  // Re-evaluate, never infer convergence from dt*RHS.
    const auto tips=PrimitiveEulerWeakTEConfig::LocalNodes(*g,*c);
    const auto tip=tips.empty()?g->GetnPoint():tips.front();
    std::vector<bool> boundary(g->GetnPoint(),false);
    std::vector<bool> wall(g->GetnPoint(),false),farfield(g->GetnPoint(),false);
    const auto& metrics=s->GetGeometryMetrics();
    for (const auto& f:metrics.GetBoundaryFaces()) {
      boundary[f.node]=true;
      wall[f.node]=wall[f.node] || c->GetMarker_All_KindBC(f.marker)==EULER_WALL || c->GetMarker_All_KindBC(f.marker)==SYMMETRY_PLANE;
      farfield[f.node]=farfield[f.node] || c->GetMarker_All_KindBC(f.marker)==FAR_FIELD;
    }
    double wallSquare=0,wallVolume=0,farSquare=0,farVolume=0;
    double square=0,interior=0,volume=0,interiorVolume=0,tipSquare=0,maxMach=0,maxTrace=0;
    const double gamma=SU2_TYPE::GetValue(c->GetGamma());
    const auto mach=[gamma](const PrimitiveEuler2D::State& w) {
      return std::hypot(w[1],w[2])/std::sqrt(gamma*w[3]/w[0]);
    };
    for (unsigned long i=0;i<g->GetnPoint();++i) {
      const double v=metrics.GetCell(i).volume;
      double norm=0; for (auto r:current.rhs[i]) norm+=r*r;
      square+=v*norm; volume+=v;
      if (!boundary[i]) { interior+=v*norm; interiorVolume+=v; }
      if (i==tip) tipSquare=norm;
      else if (wall[i]) { wallSquare+=v*norm; wallVolume+=v; }
      else if (farfield[i]) { farSquare+=v*norm; farVolume+=v; }
      maxMach=std::max(maxMach,mach(s->GetPrimitiveNodes().GetPrimitiveState(i)));
    }
    for (unsigned long e=0;e<metrics.GetnEdge();++e) {
      const auto& f=metrics.GetFace(e);
      for (auto i:{f.first,f.second}) maxTrace=std::max(maxTrace,mach(s->ReconstructAt(i,f.point)));
    }
    for (const auto& f:metrics.GetBoundaryFaces()) maxTrace=std::max(maxTrace,mach(s->ReconstructAt(f.node,f.point)));
    SetHistoryOutputValue("RHS_CURRENT",std::sqrt(square/volume));
    SetHistoryOutputValue("RHS_INTERIOR",interiorVolume>0?std::sqrt(interior/interiorVolume):0);
    SetHistoryOutputValue("RHS_TE",std::sqrt(tipSquare));
    SetHistoryOutputValue("MAX_MACH",maxMach); SetHistoryOutputValue("MAX_TRACE_MACH",maxTrace);
    if (c->GetPrimitiveLocalPseudoTime()) {
      SetHistoryOutputValue("RHS_WALL_STRONG",wallVolume>0?std::sqrt(wallSquare/wallVolume):0);
      SetHistoryOutputValue("RHS_FARFIELD",farVolume>0?std::sqrt(farSquare/farVolume):0);
    }
  }
  const auto aero=s->ComputeAerodynamics(g,c);
  if (aero.enabled) {
    SetHistoryOutputValue("LIFT",aero.cl);
    SetHistoryOutputValue("DRAG",aero.cd);
  }
  std::array<double,4> square{};
  double total=0,minrho=std::numeric_limits<double>::max(),minp=minrho,theta=1;
  for (unsigned long i=0;i<g->GetnPoint();++i) {
    const double v=s->GetGeometryMetrics().GetCell(i).volume;
    total+=v;
    const auto w=s->GetPrimitiveNodes().GetPrimitiveState(i);
    minrho=std::min(minrho,w[0]); minp=std::min(minp,w[3]);theta=std::min(theta,s->GetPositivityScale(i));
    for (unsigned k=0;k<4;++k) {
      const double rhs=-SU2_TYPE::GetValue(s->LinSysRes(i,k))/v;
      square[k]+=v*rhs*rhs;
    }
  }
  for (unsigned k=0;k<4;++k) SetHistoryOutputValue("RHS_"+key[k],std::sqrt(square[k]/total));
  SetHistoryOutputValue("MIN_RHO",minrho); SetHistoryOutputValue("MIN_P",minp);SetHistoryOutputValue("MIN_THETA",theta);
  SetHistoryOutputValue("CUR_TIME",s->GetPhysicalTime());
  if (c->GetTimeIter()==0) SetHistoryOutputValue("TIME_STEP",0);
}
void CPrimitiveEulerOutput::SetVolumeOutputFields(CConfig* c) {
  AddCoordinates();
  for (unsigned k=0;k<4;++k) AddVolumeOutput(key[k],label[k],"SOLUTION","Nodal primitive "+key[k]);
  for (const auto& field:std::array<std::string,7>{{"Speed","Mach","Temperature","Entropy","Dual_Volume","Physical_Time","Global_ID"}})
    if (field!="Physical_Time" || !c->GetPrimitiveLocalPseudoTime()) AddVolumeOutput(field,field,"PRIMITIVE",field);
  if (c->GetPrimitiveLocalPseudoTime()) {
    AddVolumeOutput("Pseudo_Iteration","Pseudo_Iteration","PRIMITIVE","Stationary iteration index; not a physical time");
    AddVolumeOutput("Local_Pseudo_Step","Local_Pseudo_Step","PRIMITIVE","Accepted local pseudo-step; zero before first iteration");
    for (const auto& field:{"RhsRho","RhsU","RhsV","RhsP"})
      AddVolumeOutput(field,field,"PRIMITIVE","Current spatial RHS after wall reaction, not the update");
  }
}
void CPrimitiveEulerOutput::LoadVolumeData(CConfig* c,CGeometry* g,CSolver** solvers,unsigned long i) {
  auto* s=static_cast<CPrimitiveEulerSolver*>(solvers[FLOW_SOL]);
  const auto w=s->GetPrimitiveNodes().GetPrimitiveState(i);
  const double gamma=SU2_TYPE::GetValue(c->GetGamma()),speed=std::hypot(w[1],w[2]);
  LoadCoordinates(g->nodes->GetCoord(i),i);
  for (unsigned k=0;k<4;++k) SetVolumeOutputValue(key[k],i,w[k]);
  SetVolumeOutputValue("Speed",i,speed);
  SetVolumeOutputValue("Mach",i,speed/std::sqrt(gamma*w[3]/w[0]));
  SetVolumeOutputValue("Temperature",i,w[3]/(w[0]*c->GetGas_Constant()));
  SetVolumeOutputValue("Entropy",i,w[3]/std::pow(w[0],gamma));
  SetVolumeOutputValue("Dual_Volume",i,s->GetGeometryMetrics().GetCell(i).volume);
  if (c->GetPrimitiveLocalPseudoTime()) {
    SetVolumeOutputValue("Pseudo_Iteration",i,c->GetTimeIter());
    const auto& steps=s->GetPseudoDiagnostics().steps;
    SetVolumeOutputValue("Local_Pseudo_Step",i,steps.empty()?0:steps.at(i));
    const std::array<std::string,4> rhsKey{{"RhsRho","RhsU","RhsV","RhsP"}};
    for (unsigned k=0;k<4;++k)
      SetVolumeOutputValue(rhsKey[k],i,-s->LinSysRes(i,k)/s->GetGeometryMetrics().GetCell(i).volume);
  }
  else SetVolumeOutputValue("Physical_Time",i,s->GetPhysicalTime());
  SetVolumeOutputValue("Global_ID",i,g->nodes->GetGlobalIndex(i));
}
bool CPrimitiveEulerOutput::WriteVolumeOutput(CConfig* c,unsigned long,bool force,unsigned short file) {
  return force || c->GetTimeIter()%c->GetVolumeOutputFrequency(file)==0;
}
void CPrimitiveEulerOutput::WriteAdditionalFiles(CConfig* c,CGeometry* g,CSolver** solvers) {
  // Called once per requested volume format: write the half-face table only once per output step.
  if (c->GetnMarker_Monitoring() && lastPressureOutput!=c->GetTimeIter()) {
    const auto aero=static_cast<CPrimitiveEulerSolver*>(solvers[FLOW_SOL])->ComputeAerodynamics(g,c);
    const auto name=c->GetSurfCoeff_FileName()+"_primitive_pressure_"+std::to_string(c->GetTimeIter())+".csv";
    std::ofstream pressure(name);
    if (!pressure) SU2_MPI::Error("Cannot write primitive pressure surface table",CURRENT_FUNCTION);
    pressure<<std::setprecision(17)<<"marker,node,x,y,p,Cp,nx,ny\n";
    for (const auto& f:aero.surface)
      pressure<<std::quoted(c->GetMarker_All_TagBound(f.marker))<<','<<g->nodes->GetGlobalIndex(f.node)<<','
              <<f.point[0]<<','<<f.point[1]<<','<<f.pressure<<','<<f.cp<<','<<f.normal[0]<<','<<f.normal[1]<<'\n';
    pressure.close();
    if (!pressure) SU2_MPI::Error("Failed writing primitive pressure surface table",CURRENT_FUNCTION);
    if (c->GetPrimitiveWeakTEGlobalID()>=0) {
      auto* solver=static_cast<CPrimitiveEulerSolver*>(solvers[FLOW_SOL]);
      const auto current=CurrentResidual(*solver,g,c);
      const auto tip=PrimitiveEulerWeakTEConfig::LocalNodes(*g,*c).front();
      const auto w=solver->GetPrimitiveNodes().GetPrimitiveState(tip);
      std::ofstream te(c->GetSurfCoeff_FileName()+"_primitive_te_"+std::to_string(c->GetTimeIter())+".csv");
      if (!te) SU2_MPI::Error("Cannot write primitive TE diagnostic",CURRENT_FUNCTION);
      te<<std::setprecision(17)<<"iteration,"<<(c->GetPrimitiveLocalPseudoTime()?"pseudo_iteration":"time")
        <<",global_id,x,y,rho,u,v,p,rhs_rho,rhs_u,rhs_v,rhs_p,trace_un,minus_rho,minus_u,minus_v,minus_p\n";
      for (const auto& f:solver->GetGeometryMetrics().GetBoundaryFaces()) if (f.node==tip) {
        const double length=std::hypot(f.normal[0],f.normal[1]);
        const PrimitiveEuler2D::Direction n{{f.normal[0]/length,f.normal[1]/length}};
        const auto trace=solver->ReconstructAt(tip,f.point);
        const auto exterior=PrimitiveEuler2D::Boundary::Reflect(trace,n,SU2_TYPE::GetValue(c->GetGamma()));
        const auto minus=PrimitiveEuler2D::ComputeFluctuations(trace,exterior,f.normal,SU2_TYPE::GetValue(c->GetGamma())).minus;
        te<<c->GetTimeIter()<<','<<(c->GetPrimitiveLocalPseudoTime()?double(c->GetTimeIter()):solver->GetPhysicalTime())
          <<','<<g->nodes->GetGlobalIndex(tip)<<','<<f.point[0]<<','<<f.point[1];
        for (auto x:w) te<<','<<x;
        for (auto x:current.rhs[tip]) te<<','<<x;
        te<<','<<(trace[1]*n[0]+trace[2]*n[1]);
        for (auto x:minus) te<<','<<x;
        te<<'\n';
      }
      te.close();
      if (!te) SU2_MPI::Error("Failed writing primitive TE diagnostic",CURRENT_FUNCTION);
    }
    lastPressureOutput=c->GetTimeIter();
  }
  if (c->GetPrimitiveGeometryAudit().empty() || c->GetTimeIter()!=0) return;
  const auto& metrics=static_cast<CPrimitiveEulerSolver*>(solvers[FLOW_SOL])->GetGeometryMetrics();
  std::ofstream file(c->GetPrimitiveGeometryAudit());
  if (!file) SU2_MPI::Error("Cannot write PRIMITIVE_GEOMETRY_AUDIT file.",CURRENT_FUNCTION);
  file<<std::setprecision(17)<<"kind,first,second,x,y,nx,ny,volume\n";
  for (unsigned long i=0;i<g->GetnPoint();++i) {
    const auto& cell=metrics.GetCell(i);
    file<<"node,"<<g->nodes->GetGlobalIndex(i)<<",0,"<<g->nodes->GetCoord(i,0)<<','<<g->nodes->GetCoord(i,1)
        <<','<<cell.offset[0]<<','<<cell.offset[1]<<','<<cell.volume<<'\n';
  }
  for (unsigned long e=0;e<metrics.GetnEdge();++e) {
    const auto& f=metrics.GetFace(e);
    file<<"edge,"<<g->nodes->GetGlobalIndex(f.first)<<','<<g->nodes->GetGlobalIndex(f.second)<<','
        <<f.point[0]<<','<<f.point[1]<<','<<f.normal[0]<<','<<f.normal[1]<<",0\n";
  }
  for (const auto& f:metrics.GetBoundaryFaces())
    file<<"boundary,"<<g->nodes->GetGlobalIndex(f.node)<<','<<f.marker<<','<<f.point[0]<<','<<f.point[1]
        <<','<<f.normal[0]<<','<<f.normal[1]<<",0\n";
}
