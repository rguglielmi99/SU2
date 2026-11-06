/*!
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "../../include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../include/numerics/flow/primitive_euler_boundary.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"
#include <limits>
#include <map>

namespace {
struct PrimitiveEvaluationGuard {
  std::atomic_flag& flag;
  explicit PrimitiveEvaluationGuard(std::atomic_flag& value) : flag(value) {
    if (omp_in_parallel()) throw std::logic_error("Primitive operator must be called outside OpenMP regions");
    if (flag.test_and_set()) throw std::logic_error("Concurrent primitive evaluation on the same solver");
  }
  ~PrimitiveEvaluationGuard() { flag.clear(); }
};
double PrimitiveFaceSpeed(const PrimitiveEuler2D::State& left, const PrimitiveEuler2D::State& right,
                          const PrimitiveEuler2D::Direction& normal, double gamma) {
  const auto path = PrimitiveEuler2D::PathSpectrum(left, right, normal, gamma);
  double speed = std::abs(path.advection) + path.acoustic;
  for (const auto& w : {left, right}) {
    const double physical = std::abs(w[1]*normal[0] + w[2]*normal[1]) +
                            std::sqrt(gamma*w[3]/w[0])*std::hypot(normal[0],normal[1]);
    speed = std::max(speed, physical);
  }
  if (!std::isfinite(speed) || speed <= 0) throw std::invalid_argument("Invalid primitive characteristic speed");
  return speed;
}
}

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
  throw std::logic_error("Primitive Euler phase 5: use Forward Euler; this integration/residual route is unavailable");
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
    // Slip constrains nodal velocity, not flat-symmetry derivatives: retain
    // Euler-wall gradients while preserving corrections on true SYMMETRY_PLANE.
    const auto method = config.GetKind_Gradient_Method_Recon();
    if (method == GREEN_GAUSS) {
      computeGradientsGreenGauss(nullptr, MPI_QUANTITIES::SOLUTION_GRAD_REC, PERIODIC_NONE,
                                geometry, config, nodes->GetSolution(), 0, 4, 1, raw, false);
    } else if (method == LEAST_SQUARES || method == WEIGHTED_LEAST_SQUARES) {
      if (nodes->GetRmatrix().length() != nPoint)
        throw std::invalid_argument("Primitive reconstruction requires allocated least-squares storage");
      computeGradientsLeastSquares(nullptr, MPI_QUANTITIES::SOLUTION_GRAD_REC, PERIODIC_NONE,
          geometry, config, method == WEIGHTED_LEAST_SQUARES, nodes->GetSolution(), 0, 4, 1, raw, nodes->GetRmatrix(), false);
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

std::vector<std::vector<CPrimitiveEulerSolver::Direction>>
CPrimitiveEulerSolver::BuildSlipConstraints(const CConfig& config) const {
  // A same-marker Euler chain represents a smooth polygonal wall. Its
  // integrated half-face sum is the negative native SU2 vertex normal.
  // Distinct straight markers retain the explicitly declared phase-9 corners.
  struct Incident { Direction sum{}, first{}; double measure=0; unsigned count=0; bool bent=false; };
  std::map<std::pair<unsigned long,unsigned short>,Incident> incident;
  const double tol=256*std::numeric_limits<double>::epsilon();
  for (const auto& f:metrics->GetBoundaryFaces()) {
    const auto kind=config.GetMarker_All_KindBC(f.marker);
    if (kind!=EULER_WALL && kind!=SYMMETRY_PLANE) continue;
    auto& a=incident[{f.node,f.marker}];
    const double length=std::hypot(f.normal[0],f.normal[1]);
    const Direction unit{{f.normal[0]/length,f.normal[1]/length}};
    if (a.count && (std::abs(a.first[0]*unit[1]-a.first[1]*unit[0])>tol ||
                    a.first[0]*unit[0]+a.first[1]*unit[1]<0)) a.bent=true;
    if (!a.count) a.first=unit;
    ++a.count; a.measure+=length;
    for (unsigned d=0;d<2;++d) a.sum[d]+=f.normal[d];
  }
  std::vector<unsigned> markerCount(nPoint,0);
  for (const auto& item:incident) ++markerCount[item.first.first];
  std::vector<std::vector<Direction>> basis(nPoint);
  for (const auto& item:incident) {
    const auto node=item.first.first;
    const auto& a=item.second;
    if (a.count>2) throw std::invalid_argument("Primitive slip non-manifold wall vertex");
    const double length=std::hypot(a.sum[0],a.sum[1]);
    if (!std::isfinite(length) || length<=tol*a.measure)
      throw std::invalid_argument("Primitive slip degenerate resultant wall normal");
    if (a.bent && config.GetMarker_All_KindBC(item.first.second)==SYMMETRY_PLANE)
      throw std::invalid_argument("Primitive symmetry marker must remain locally straight");
    if (a.bent && markerCount[node]>1)
      throw std::invalid_argument("Primitive polygonal wall intersects another slip marker: ambiguous corner");
    Direction n{{a.sum[0]/length,a.sum[1]/length}};
    for (const auto& b:basis[node]) {
      const double dot=n[0]*b[0]+n[1]*b[1];
      for (unsigned d=0;d<2;++d) n[d]-=dot*b[d];
    }
    const double magnitude=std::hypot(n[0],n[1]);
    if (magnitude>tol && basis[node].size()<2) {
      for (auto& value:n) value/=magnitude;
      basis[node].push_back(n);
    }
  }
  return basis;
}

unsigned long CPrimitiveEulerSolver::ProjectFreestreamWallVelocity() {
  PrimitiveEvaluationGuard guard(evaluating);
  if (physicalTime!=0 || originalConfig->GetPrimitiveInitialization()!="FREESTREAM")
    throw std::invalid_argument("Wall initialization projection requires FREESTREAM at time zero");
  const auto basis=BuildSlipConstraints(*originalConfig);
  std::vector<bool> wall(nPoint,false);
  for (const auto& f:metrics->GetBoundaryFaces())
    if (originalConfig->GetMarker_All_KindBC(f.marker)==EULER_WALL) wall[f.node]=true;
  unsigned long count=0;
  for (unsigned long i=0;i<nPoint;++i) if (wall[i]) {
    auto w=nodes->GetPrimitiveState(i);
    for (const auto& n:basis[i]) {
      const double un=w[1]*n[0]+w[2]*n[1];
      for (unsigned d=0;d<2;++d) w[1+d]-=un*n[d];
    }
    nodes->SetPrimitiveState(i,w);
    ++count;
  }
  return count;
}

CPrimitiveEulerSolver::SpatialResult CPrimitiveEulerSolver::ComputeSpatialResidual(
    CGeometry* geometry, const CConfig* config, double time, const BoundaryState& boundaryState, const Source& source) {
  PrimitiveEvaluationGuard guard(evaluating);
  return EvaluateSpatialResidual(geometry, config, time, boundaryState, source, true);
}

CPrimitiveEulerSolver::SpatialResult CPrimitiveEulerSolver::EvaluateSpatialResidual(
    CGeometry* geometry, const CConfig* config, double time, const BoundaryState& boundaryState,
    const Source& source, bool commit) {
  using namespace PrimitiveEuler2D;
  if (!geometry || !config || geometry != originalGeometry || config != originalConfig || !std::isfinite(time))
    throw std::invalid_argument("Primitive residual requires original static geometry/configuration and finite time");
  ValidateConfiguration(*config);
  const auto wallBasis=BuildSlipConstraints(*config);
  for (unsigned short marker = 0; marker < geometry->GetnMarker(); ++marker) {
    const auto kind = config->GetMarker_All_KindBC(marker);
    if (kind != CUSTOM_BOUNDARY && kind != FAR_FIELD && kind != EULER_WALL && kind != SYMMETRY_PLANE &&
        kind != INLET_FLOW && kind != OUTLET_FLOW && kind != SUPERSONIC_INLET && kind != SUPERSONIC_OUTLET)
      throw std::invalid_argument("Primitive unsupported physical boundary kind");
    if (kind != FAR_FIELD && kind != EULER_WALL && kind != SYMMETRY_PLANE && !boundaryState)
      throw std::invalid_argument("Primitive CUSTOM boundary requires a primitive exterior-state callback");
    if (kind==FAR_FIELD && config->GetPrimitiveBoundaryMode()=="LOCAL_CHARACTERISTIC" && !boundaryState)
      throw std::invalid_argument("Primitive LOCAL_CHARACTERISTIC FAR requires a configured physical exterior-state callback; no prescribed-state fallback");
  }
  for (unsigned long i=0;i<nPoint;++i) for (const auto& n:wallBasis[i]) {
    const auto w=nodes->GetPrimitiveState(i);
    if (std::abs(w[1]*n[0]+w[2]*n[1])>256*std::numeric_limits<double>::epsilon()*std::max(1.,std::hypot(w[1],w[2])))
      throw std::invalid_argument("Primitive slip wall requires compatible initial nodal normal velocity; no hidden state projection");
  }
  PrepareReconstruction(*geometry, *config);
  const double gamma = SU2_TYPE::GetValue(config->GetGamma());
  SpatialResult result;
  for (auto* field : {&result.faces, &result.volume, &result.boundary, &result.source, &result.rhs}) field->resize(nPoint);

  result.spectralRadius.resize(nPoint, 0.0);

  // Serial deterministic assembly: each edge evaluated ONCE, both signs +.
  // No conservative edge-flux reducer and no shared OpenMP edge scatter.
  for (unsigned long e = 0; e < metrics->GetnEdge(); ++e) {
    const auto& face = metrics->GetFace(e);
    const auto left = ReconstructAt(face.first, face.point), right = ReconstructAt(face.second, face.point);
    const double speed = PrimitiveFaceSpeed(left, right, face.normal, gamma);
    result.spectralRadius[face.first] += speed;
    result.spectralRadius[face.second] += speed;
    const auto waves = ComputeFluctuations(left, right, face.normal, gamma);
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
    const auto kind=config->GetMarker_All_KindBC(face.marker);
    const auto exterior = kind==EULER_WALL || kind==SYMMETRY_PLANE ?
        PrimitiveEuler2D::Boundary::Reflect(interior,unit,gamma) :
        (boundaryState ? boundaryState(face.point, time, interior, unit, face.marker) : far);
    result.spectralRadius[face.node] += PrimitiveFaceSpeed(interior, exterior, face.normal, gamma);
    const auto incoming = ComputeFluctuations(interior, exterior, face.normal, gamma).minus;
    for (unsigned k = 0; k < 4; ++k) result.boundary[face.node][k] += incoming[k];
  }
  std::vector<State> integrated(nPoint);
  for (unsigned long i = 0; i < nPoint; ++i) {
    const Direction x{{SU2_TYPE::GetValue(geometry->nodes->GetCoord(i, 0)),
                       SU2_TYPE::GetValue(geometry->nodes->GetCoord(i, 1))}};
    if (source) result.source[i] = source(x, time);
    const double volume = metrics->GetCell(i).volume;
    if (!std::isfinite(result.spectralRadius[i]) || result.spectralRadius[i] <= 0)
      throw std::invalid_argument("Invalid primitive spectral sum");
    for (unsigned k = 0; k < 4; ++k) {
      result.rhs[i][k] = -(result.faces[i][k] + result.volume[i][k] + result.boundary[i][k]) / volume + result.source[i][k];
      integrated[i][k] = -volume * result.rhs[i][k];
      if (!std::isfinite(result.rhs[i][k]) || !std::isfinite(integrated[i][k]))
        throw std::invalid_argument("Non-finite primitive spatial residual/source");
    }
  }
  // Slip is a constraint on nodal velocity, not an artificial volume source.
  // Add the boundary reaction that makes the COMPLETE velocity RHS tangent.
  // Scalar equations, interior faces, volume, sources and FE remain untouched.
  for (unsigned long i=0;i<nPoint;++i) {
    const double volume=metrics->GetCell(i).volume;
    for (const auto& n:wallBasis[i]) {
      const double normalRhs=result.rhs[i][1]*n[0]+result.rhs[i][2]*n[1];
      for (unsigned d=0;d<2;++d) {
        const double correction=normalRhs*n[d];
        result.rhs[i][1+d]-=correction;
        result.boundary[i][1+d]+=volume*correction;
      }
    }
  }
  if (commit) CommitSpatialResidual(result);
  return result;
}

void CPrimitiveEulerSolver::CommitSpatialResidual(const SpatialResult& result) {
  for (unsigned long i = 0; i < nPoint; ++i)
    for (unsigned k = 0; k < 4; ++k) LinSysRes(i,k) = -metrics->GetCell(i).volume * result.rhs[i][k];
}

double CPrimitiveEulerSolver::EstimateTimeStep(CGeometry* geometry, const CConfig* config,
                                              double time, double cfl, const BoundaryState& boundary) {
  PrimitiveEvaluationGuard guard(evaluating);
  if (!std::isfinite(cfl) || cfl <= 0) throw std::invalid_argument("Primitive CFL must be finite and positive");
  const auto result = EvaluateSpatialResidual(geometry, config, time, boundary, {}, false);
  double dt = std::numeric_limits<double>::max();
  for (unsigned long i = 0; i < nPoint; ++i)
    dt = std::min(dt, cfl * (metrics->GetCell(i).volume / result.spectralRadius[i]));
  if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(time+dt) || time+dt <= time)
    throw std::invalid_argument("Primitive time step is not finite, positive and representable");
  return dt;
}

CPrimitiveEulerSolver::SpatialResult CPrimitiveEulerSolver::AdvanceExplicitEuler(
    CGeometry* geometry, const CConfig* config, double time, double dt,
    const BoundaryState& boundary, const Source& source) {
  PrimitiveEvaluationGuard guard(evaluating);
  if (!std::isfinite(time) || !std::isfinite(dt) || dt <= 0 ||
      !std::isfinite(time+dt) || time+dt <= time)
    throw std::invalid_argument("Primitive Forward Euler requires a finite positive representable step");
  auto result = EvaluateSpatialResidual(geometry, config, time, boundary, source, false);
  std::vector<PrimitiveEuler2D::State> next(nPoint);
  const double gamma = SU2_TYPE::GetValue(config->GetGamma());
  for (unsigned long i = 0; i < nPoint; ++i) {
    next[i] = nodes->GetPrimitiveState(i);
    for (unsigned k = 0; k < 4; ++k) next[i][k] += dt * result.rhs[i][k];
    try {
      PrimitiveEuler2D::Validate(next[i], gamma);
    } catch (const std::invalid_argument&) {
      throw std::domain_error("Primitive Forward Euler rejected non-admissible candidate at node " + std::to_string(i));
    }
  }
  // All candidate states validated before touching W, old solution, residual or time.
  nodes->Set_OldSolution();
  for (unsigned long i = 0; i < nPoint; ++i) nodes->SetPrimitiveState(i, next[i]);
  CommitSpatialResidual(result);
  physicalTime = time + dt;
  return result;
}

void CPrimitiveEulerSolver::SetTime_Step(CGeometry* geometry, CSolver**, CConfig* config,
                                        unsigned short mesh, unsigned long) {
  if (!config || mesh != MESH_0) throw std::invalid_argument("Primitive time step requires finest mesh/configuration");
  const double dt = EstimateTimeStep(geometry, config, physicalTime, SU2_TYPE::GetValue(config->GetCFL(mesh)));
  for (unsigned long i = 0; i < nPoint; ++i) nodes->SetDelta_Time(i, dt);
  Min_Delta_Time = Max_Delta_Time = dt;
}

void CPrimitiveEulerSolver::ExplicitEuler_Iteration(CGeometry* geometry, CSolver**, CConfig* config) {
  const double dt = SU2_TYPE::GetValue(nodes->GetDelta_Time(0));
  for (unsigned long i = 1; i < nPoint; ++i)
    if (SU2_TYPE::GetValue(nodes->GetDelta_Time(i)) != dt)
      throw std::invalid_argument("Primitive physical Forward Euler requires one global time step");
  AdvanceExplicitEuler(geometry, config, physicalTime, dt);
}

void CPrimitiveEulerSolver::Upwind_Residual(CGeometry* geometry, CSolver**, CNumerics**, CConfig* config,
                                           unsigned short mesh) {
  if (mesh != MESH_0) throw std::invalid_argument("Primitive residual requires the finest mesh");
  ComputeSpatialResidual(geometry, config);
}
