/*!
 * \brief Primitive Euler NODE/AGGREGATED spatial operator; Forward Euler time integration.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "CSolver.hpp"
#include "../variables/CPrimitiveEulerVariable.hpp"
#include "../geometry/primitive_euler_geometry.hpp"
#include "../numerics/flow/primitive_euler_reconstruction.hpp"
#include <functional>
#include <atomic>

class CPrimitiveEulerSolver final : public CSolver {
  std::unique_ptr<CPrimitiveEulerVariable> nodes;
  std::unique_ptr<PrimitiveEuler2D::CAggregatedGeometry> metrics;
  std::vector<PrimitiveEuler2D::Gradient> limitedGradient;
  std::vector<double> positivityScale;
  const CGeometry* originalGeometry;
  const CConfig* originalConfig;
  std::atomic_flag evaluating = ATOMIC_FLAG_INIT;
  double physicalTime = 0;
  CVariable* GetBaseClassPointerToNodes() override { return nodes.get(); }
  void PrepareReconstruction(CGeometry& geometry, const CConfig& config);
  std::vector<std::vector<PrimitiveEuler2D::Direction>> BuildSlipConstraints(const CConfig& config) const;

 public:
  CPrimitiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short mesh);
  CPrimitiveEulerVariable& GetPrimitiveNodes() { return *nodes; }
  const CPrimitiveEulerVariable& GetPrimitiveNodes() const { return *nodes; }
  const PrimitiveEuler2D::CAggregatedGeometry& GetGeometryMetrics() const { return *metrics; }
  static void ValidateConfiguration(const CConfig& config);
  static void EvolutionUnavailable();
  // Explicit FREESTREAM-only initialization policy; leaves rho, p and interior nodes unchanged.
  unsigned long ProjectFreestreamWallVelocity();

  using State = PrimitiveEuler2D::State;
  using Direction = PrimitiveEuler2D::Direction;
  // Boundary data are exterior primitive states at the actual half-edge midpoint.
  using BoundaryState = std::function<State(const Direction&, double, const State&, const Direction&, unsigned short)>;
  // NODE source values are sampled at primal nodes, not dual centroids.
  using Source = std::function<State(const Direction&, double)>;
  struct SpatialResult {
    std::vector<State> faces, volume, boundary, source, rhs;
    std::vector<double> spectralRadius;  // Sum of integrated face wave speeds, not divided by V.
  };
  // Serial entry point, outside OpenMP parallel regions. Rebuilds reconstruction
  // on every call, never changes W. LinSysRes receives faces+volume+boundary-V*S
  // only after all contributions have been validated; rhs=-LinSysRes/V.
  SpatialResult ComputeSpatialResidual(CGeometry* geometry, const CConfig* config, double time = 0,
                                      const BoundaryState& boundary = {}, const Source& source = {});
  // Global physical step estimate; includes actual reconstructed/exterior path speeds.
  // A heuristic CFL scale, not a positivity/stability theorem for the full operator.
  double EstimateTimeStep(CGeometry*, const CConfig*, double time, double cfl,
                          const BoundaryState& boundary = {});
  // Transactional Forward Euler. Invalid candidates reject the WHOLE step;
  // no component clipping, hidden dt reduction or source rescaling.
  SpatialResult AdvanceExplicitEuler(CGeometry*, const CConfig*, double time, double dt,
                                    const BoundaryState& boundary = {}, const Source& source = {});
  double GetPhysicalTime() const { return physicalTime; }
  State ReconstructAt(unsigned long point, const Direction& coordinate) const;
  const PrimitiveEuler2D::Gradient& GetLimitedGradient(unsigned long point) const { return limitedGradient.at(point); }
  double GetPositivityScale(unsigned long point) const { return positivityScale.at(point); }

  void Upwind_Residual(CGeometry*, CSolver**, CNumerics**, CConfig*, unsigned short) override;
  void Centered_Residual(CGeometry*, CSolver**, CNumerics**, CConfig*, unsigned short, unsigned short) override {
    EvolutionUnavailable();
  }
  void Source_Residual(CGeometry*, CSolver**, CNumerics**, CConfig*, unsigned short) override {
    EvolutionUnavailable();
  }
  void SetTime_Step(CGeometry*, CSolver**, CConfig*, unsigned short, unsigned long) override;
  void ExplicitEuler_Iteration(CGeometry*, CSolver**, CConfig*) override;
  void ExplicitRK_Iteration(CGeometry*, CSolver**, CConfig*, unsigned short) override { EvolutionUnavailable(); }
  void ClassicalRK4_Iteration(CGeometry*, CSolver**, CConfig*, unsigned short) override { EvolutionUnavailable(); }
  void ImplicitEuler_Iteration(CGeometry*, CSolver**, CConfig*) override { EvolutionUnavailable(); }
 private:
  SpatialResult EvaluateSpatialResidual(CGeometry*, const CConfig*, double, const BoundaryState&, const Source&, bool commit);
  void CommitSpatialResidual(const SpatialResult&);
};
