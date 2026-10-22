/*!
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "../../include/solvers/CPrimitiveEulerSolver.hpp"

void CPrimitiveEulerSolver::ValidateConfiguration(const CConfig& config) {
  if (config.GetPrimitiveRepresentation() != "NODE" || config.GetPrimitiveFaceQuadrature() != "AGGREGATED")
    throw std::invalid_argument("Primitive Euler currently requires NODE and AGGREGATED");
  if (config.GetPrimitiveVolumeQuadrature() != "MEAN_COEFF" && config.GetPrimitiveVolumeQuadrature() != "NODE_POINT")
    throw std::invalid_argument("Primitive volume quadrature must be MEAN_COEFF or NODE_POINT");
  if (config.GetKind_FluidModel() != STANDARD_AIR && config.GetKind_FluidModel() != IDEAL_GAS)
    throw std::invalid_argument("Primitive Euler requires an ideal gas");
  if (config.GetKind_TimeIntScheme_Flow() != EULER_EXPLICIT || config.GetnMGLevels() != 0 ||
      config.GetNewtonKrylov() || config.GetRestart() || config.GetRestart_Flow() || config.GetMultizone_Problem() ||
      config.GetnZone() != 1 || config.GetAxisymmetric() || config.AddRadiation() ||
      config.GetKind_Turb_Model() != TURB_MODEL::NONE || config.GetKind_Species_Model() != SPECIES_MODEL::NONE)
    throw std::invalid_argument("Primitive Euler requires explicit Euler, MGLEVEL=0, no restart or coupled models");
  if (config.GetTime_Marching() != TIME_MARCHING::STEADY &&
      config.GetTime_Marching() != TIME_MARCHING::TIME_STEPPING)
    throw std::invalid_argument("Primitive Euler does not support dual-time or harmonic-balance integration");
}

CPrimitiveEulerSolver::CPrimitiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short mesh) : CSolver() {
  if (!geometry || !config || mesh != MESH_0) throw std::invalid_argument("Invalid primitive Euler solver setup");
  ValidateConfiguration(*config);
  metrics.reset(new PrimitiveEuler2D::CAggregatedGeometry(*geometry, *config));
  PrimitiveEuler2D::State initial{};
  for (unsigned k = 0; k < 4; ++k) initial[k] = SU2_TYPE::GetValue(config->GetPrimitiveInitialState()[k]);
  nodes.reset(new CPrimitiveEulerVariable(initial, geometry->GetnPoint(), config));
  nDim = 2;
  nVar = nPrimVar = nPrimVarGrad = nVarGrad = 4;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nMarker = geometry->GetnMarker();
  SolverName = "PRIMITIVE_EULER";
  Residual_RMS.resize(nVar, 0.0);
  Residual_Max.resize(nVar, 0.0);
  Point_Max.resize(nVar, 0);
  Point_Max_Coord.resize(nVar, nDim) = su2double(0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  SetBaseClassPointerToNodes();
}

void CPrimitiveEulerSolver::EvolutionUnavailable() {
  throw std::logic_error("Primitive Euler phase 3: residual, boundary conditions and time integration "
                         "are not implemented");
}
