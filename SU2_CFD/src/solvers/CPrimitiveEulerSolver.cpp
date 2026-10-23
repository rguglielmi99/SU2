/*!
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "../../include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"

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
  const auto limiter = config.GetKind_SlopeLimit_Flow();
  if (limiter != LIMITER::NONE && limiter != LIMITER::BARTH_JESPERSEN && limiter != LIMITER::VENKATAKRISHNAN)
    throw std::invalid_argument("Primitive Euler supports NONE, BARTH_JESPERSEN or VENKATAKRISHNAN limiters");
  if (config.GetMUSCL_Kappa_Flow() != 0)
    throw std::invalid_argument("Primitive Euler affine reconstruction requires MUSCL_KAPPA_FLOW=0");

}

CPrimitiveEulerSolver::CPrimitiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short mesh)
    : CSolver(), originalGeometry(geometry), originalConfig(config) {
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
  limitedGradient.resize(nPoint);
  positivityScale.resize(nPoint, 1.0);
  SetBaseClassPointerToNodes();
}

void CPrimitiveEulerSolver::EvolutionUnavailable() {
  throw std::logic_error("Primitive Euler phase 4: use ComputeSpatialResidual; driver integration and time "
                         "advancement are not implemented");
}

void CPrimitiveEulerSolver::PrepareReconstruction(CGeometry& geometry, const CConfig& config) {
  using namespace PrimitiveEuler2D;
  const double gamma = SU2_TYPE::GetValue(config.GetGamma());
  for (unsigned long i = 0; i < nPoint; ++i) Validate(nodes->GetPrimitiveState(i), gamma);

  auto& raw = nodes->GetGradient_Reconstruction();
  auto& limiter = nodes->GetLimiter();
  for (unsigned long i = 0; i < nPoint; ++i)
    for (unsigned k = 0; k < 4; ++k)
      for (unsigned d = 0; d < 2; ++d) raw(i, k, d) = 0;
  limiter = su2double(1.0);
  if (config.GetMUSCL_Flow()) {
    // Use the reconstruction method for BOTH faces and volume, even when the
    // separate generic gradient method differs. No flow-layout assumptions.
    const auto method = config.GetKind_Gradient_Method_Recon();
    if (method == GREEN_GAUSS) {
      computeGradientsGreenGauss(nullptr, MPI_QUANTITIES::SOLUTION_GRAD_REC, PERIODIC_NONE,
                                geometry, config, nodes->GetSolution(), 0, 4, 1, raw);
    } else if (method == LEAST_SQUARES || method == WEIGHTED_LEAST_SQUARES) {
      if (nodes->GetRmatrix().length() != nPoint)
        throw std::invalid_argument("Primitive reconstruction requires allocated least-squares storage");
      computeGradientsLeastSquares(nullptr, MPI_QUANTITIES::SOLUTION_GRAD_REC, PERIODIC_NONE,
          geometry, config, method == WEIGHTED_LEAST_SQUARES, nodes->GetSolution(), 0, 4, 1, raw, nodes->GetRmatrix());
    } else {
      throw std::invalid_argument("Unsupported primitive reconstruction gradient method");
    }
    // Do not use GetKind_SlopeLimit(): the generic runtime option is set by
    // integration machinery that intentionally remains disabled in phase 4.
    computeLimiters(config.GetKind_SlopeLimit_Flow(), nullptr, MPI_QUANTITIES::SOLUTION_LIMITER,
        PERIODIC_NONE, PERIODIC_NONE, geometry, config, 0, 4, 0.0, nodes->GetSolution(), raw,
        nodes->GetSolution_Min(), nodes->GetSolution_Max(), limiter);
  }
  for (unsigned long i = 0; i < nPoint; ++i) {
    positivityScale[i] = 1;
    for (unsigned k = 0; k < 4; ++k) {
      const double phi = SU2_TYPE::GetValue(limiter(i, k));
      if (!std::isfinite(phi) || phi < 0) throw std::invalid_argument("Invalid primitive slope limiter");
      for (unsigned d = 0; d < 2; ++d) {
        limitedGradient[i][k][d] = phi * SU2_TYPE::GetValue(raw(i, k, d));
        if (!std::isfinite(limitedGradient[i][k][d])) throw std::invalid_argument("Non-finite primitive gradient");
      }
    }
  }
  // Contract all four components by the same nodal theta. Include every point
  // used below: effective interior faces, physical half-faces and dual means.
  const auto constrain = [&](unsigned long i, const Direction& displacement) {
    const auto w = nodes->GetPrimitiveState(i);
    const auto candidate = Reconstruct(w, limitedGradient[i], displacement);
    for (unsigned k = 0; k < 4; ++k)
      if (!std::isfinite(candidate[k])) throw std::invalid_argument("Non-finite primitive reconstruction");
    for (unsigned k : {0u, 3u}) {
      positivityScale[i] = std::min(positivityScale[i], ReconstructionPositivityFactor(w[k], candidate[k]));
    }
  };
  const auto displacement = [&](unsigned long i, const Direction& q) {
    return Direction{{q[0] - SU2_TYPE::GetValue(geometry.nodes->GetCoord(i, 0)),
                      q[1] - SU2_TYPE::GetValue(geometry.nodes->GetCoord(i, 1))}};
  };
  for (unsigned long e = 0; e < metrics->GetnEdge(); ++e) {
    const auto& face = metrics->GetFace(e);
    constrain(face.first, displacement(face.first, face.point));
    constrain(face.second, displacement(face.second, face.point));
  }
  for (const auto& face : metrics->GetBoundaryFaces()) constrain(face.node, displacement(face.node, face.point));
  if (config.GetPrimitiveVolumeQuadrature() == "MEAN_COEFF")
    for (unsigned long i = 0; i < nPoint; ++i) constrain(i, metrics->GetCell(i).offset);
  for (unsigned long i = 0; i < nPoint; ++i) {
    for (auto& row : limitedGradient[i])
      for (auto& value : row) value *= positivityScale[i];
  }
}

CPrimitiveEulerSolver::State CPrimitiveEulerSolver::ReconstructAt(unsigned long point, const Direction& q) const {
  const auto w = nodes->GetPrimitiveState(point);
  const auto* x = originalGeometry->nodes->GetCoord(point);
  const Direction displacement{{q[0] - SU2_TYPE::GetValue(x[0]), q[1] - SU2_TYPE::GetValue(x[1])}};
  const auto result = PrimitiveEuler2D::Reconstruct(w, limitedGradient.at(point), displacement);
  PrimitiveEuler2D::Validate(result, SU2_TYPE::GetValue(originalConfig->GetGamma()));
  return result;
}

CPrimitiveEulerSolver::SpatialResult CPrimitiveEulerSolver::ComputeSpatialResidual(
    CGeometry* geometry, const CConfig* config, double time, const BoundaryState& boundaryState, const Source& source) {
  using namespace PrimitiveEuler2D;
  if (omp_in_parallel()) throw std::logic_error("Primitive spatial operator must be called outside OpenMP regions");
  if (evaluating.test_and_set()) throw std::logic_error("Concurrent primitive spatial evaluation on the same solver");
  struct Unlock {
    std::atomic_flag& flag;
    ~Unlock() { flag.clear(); }
  } unlock{evaluating};
  if (!geometry || !config || geometry != originalGeometry || config != originalConfig || !std::isfinite(time))
    throw std::invalid_argument("Primitive residual requires original static geometry/configuration and finite time");
  ValidateConfiguration(*config);
  for (unsigned short marker = 0; marker < geometry->GetnMarker(); ++marker) {
    const auto kind = config->GetMarker_All_KindBC(marker);
    if (kind != CUSTOM_BOUNDARY && kind != FAR_FIELD)
      throw std::invalid_argument("Primitive phase 4 boundaries support CUSTOM or FAR only");
    if (kind == CUSTOM_BOUNDARY && !boundaryState)
      throw std::invalid_argument("Primitive CUSTOM boundary requires a primitive exterior-state callback");
  }
  PrepareReconstruction(*geometry, *config);
  const double gamma = SU2_TYPE::GetValue(config->GetGamma());
  SpatialResult result;
  for (auto* field : {&result.faces, &result.volume, &result.boundary, &result.source, &result.rhs}) field->resize(nPoint);

  // Serial deterministic assembly: each edge evaluated ONCE, both signs +.
  // No conservative edge-flux reducer and no shared OpenMP edge scatter.
  for (unsigned long e = 0; e < metrics->GetnEdge(); ++e) {
    const auto& face = metrics->GetFace(e);
    const auto waves = ComputeFluctuations(ReconstructAt(face.first, face.point),
                                         ReconstructAt(face.second, face.point), face.normal, gamma);
    for (unsigned k = 0; k < 4; ++k) {
      result.faces[face.first][k] += waves.minus[k];
      result.faces[face.second][k] += waves.plus[k];
    }
  }
  for (unsigned long i = 0; i < nPoint; ++i) {
    auto mean = nodes->GetPrimitiveState(i);
    if (config->GetPrimitiveVolumeQuadrature() == "MEAN_COEFF")
      mean = Reconstruct(mean, limitedGradient[i], metrics->GetCell(i).offset);
    result.volume[i] = VolumeProduct(mean, limitedGradient[i], gamma);
    for (auto& value : result.volume[i]) value *= metrics->GetCell(i).volume;
  }
  State far{};
  for (unsigned k = 0; k < 4; ++k) far[k] = SU2_TYPE::GetValue(config->GetPrimitiveInitialState()[k]);
  for (const auto& face : metrics->GetBoundaryFaces()) {
    const auto interior = ReconstructAt(face.node, face.point);
    const double measure = std::hypot(face.normal[0], face.normal[1]);
    const Direction unit{{face.normal[0] / measure, face.normal[1] / measure}};
    const auto exterior = boundaryState ? boundaryState(face.point, time, interior, unit, face.marker) : far;
    const auto incoming = ComputeFluctuations(interior, exterior, face.normal, gamma).minus;
    for (unsigned k = 0; k < 4; ++k) result.boundary[face.node][k] += incoming[k];
  }
  std::vector<State> integrated(nPoint);
  for (unsigned long i = 0; i < nPoint; ++i) {
    const Direction x{{SU2_TYPE::GetValue(geometry->nodes->GetCoord(i, 0)),
                       SU2_TYPE::GetValue(geometry->nodes->GetCoord(i, 1))}};
    if (source) result.source[i] = source(x, time);
    const double volume = metrics->GetCell(i).volume;
    for (unsigned k = 0; k < 4; ++k) {
      result.rhs[i][k] = -(result.faces[i][k] + result.volume[i][k] + result.boundary[i][k]) / volume + result.source[i][k];
      integrated[i][k] = -volume * result.rhs[i][k];
      if (!std::isfinite(result.rhs[i][k]) || !std::isfinite(integrated[i][k]))
        throw std::invalid_argument("Non-finite primitive spatial residual/source");
    }
  }
  // Commit only a complete, finite residual. No stale additions on repeated calls.
  for (unsigned long i = 0; i < nPoint; ++i)
    for (unsigned k = 0; k < 4; ++k) LinSysRes(i, k) = integrated[i][k];
  return result;
}

void CPrimitiveEulerSolver::Upwind_Residual(CGeometry* geometry, CSolver**, CNumerics**, CConfig* config,
                                           unsigned short mesh) {
  if (mesh != MESH_0) throw std::invalid_argument("Primitive residual requires the finest mesh");
  ComputeSpatialResidual(geometry, config);
}
