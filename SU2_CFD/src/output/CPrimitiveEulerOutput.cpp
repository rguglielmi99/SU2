/*!
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "../../include/output/CPrimitiveEulerOutput.hpp"
#include "../../include/solvers/CPrimitiveEulerSolver.hpp"
#include <fstream>
#include <iomanip>

namespace {
const std::array<std::string,4> key{{"RHO","U","V","P"}};
const std::array<std::string,4> label{{"Density","Velocity_x","Velocity_y","Pressure"}};
}
CPrimitiveEulerOutput::CPrimitiveEulerOutput(CConfig* c,unsigned short dim) : CFVMOutput(c,dim,false) {
  if (nRequestedHistoryFields==0) {
    RequestCommonHistory(true);
    requestedHistoryFields.emplace_back("PRIMITIVE_RHS");
    if (c->GetnMarker_Monitoring()) requestedHistoryFields.emplace_back("AERO_COEFF");
    nRequestedHistoryFields=requestedHistoryFields.size();
  }
  if (nRequestedScreenFields==0) {
    requestedScreenFields={"TIME_ITER","CUR_TIME","RHS_RHO","RHS_P"};
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
  if (c->GetnMarker_Monitoring()) {
    AddHistoryOutput("LIFT","CL",ScreenOutputFormat::SCIENTIFIC,"AERO_COEFF","Primitive reconstructed pressure lift",HistoryFieldType::COEFFICIENT);
    AddHistoryOutput("DRAG","CD",ScreenOutputFormat::SCIENTIFIC,"AERO_COEFF","Primitive reconstructed pressure drag",HistoryFieldType::COEFFICIENT);
  }
  for (unsigned k=0;k<4;++k)
    AddHistoryOutput("RHS_"+key[k],"L2_rhs_"+key[k],ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_RHS",
                     "Volume-weighted L2 of accepted-stage primitive RHS (linear, all nodes)");
  AddHistoryOutput("MIN_RHO","Min_rho",ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_STATE","Minimum updated nodal density");
  AddHistoryOutput("MIN_P","Min_p",ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_STATE","Minimum updated nodal pressure");
  AddHistoryOutput("MIN_THETA","Min_theta",ScreenOutputFormat::SCIENTIFIC,"PRIMITIVE_STATE","Minimum reconstruction positivity scale");
}
void CPrimitiveEulerOutput::LoadHistoryData(CConfig* c,CGeometry* g,CSolver** solvers) {
  auto* s=static_cast<CPrimitiveEulerSolver*>(solvers[FLOW_SOL]);
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
void CPrimitiveEulerOutput::SetVolumeOutputFields(CConfig*) {
  AddCoordinates();
  for (unsigned k=0;k<4;++k) AddVolumeOutput(key[k],label[k],"SOLUTION","Nodal primitive "+key[k]);
  for (const auto& field:std::array<std::string,7>{{"Speed","Mach","Temperature","Entropy","Dual_Volume","Physical_Time","Global_ID"}})
    AddVolumeOutput(field,field,"PRIMITIVE",field);
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
  SetVolumeOutputValue("Physical_Time",i,s->GetPhysicalTime());
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
