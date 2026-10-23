/*!
 * \brief Primitive Euler NODE/AGGREGATED spatial operator; time evolution is unavailable.
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
  CVariable* GetBaseClassPointerToNodes() override { return nodes.get(); }
  void PrepareReconstruction(CGeometry& geometry, const CConfig& config);

 public:
  CPrimitiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short mesh);
  CPrimitiveEulerVariable& GetPrimitiveNodes() { return *nodes; }
  const CPrimitiveEulerVariable& GetPrimitiveNodes() const { return *nodes; }
  const PrimitiveEuler2D::CAggregatedGeometry& GetGeometryMetrics() const { return *metrics; }
  static void ValidateConfiguration(const CConfig& config);
  static void EvolutionUnavailable();

  using State = PrimitiveEuler2D::State;
  using Direction = PrimitiveEuler2D::Direction;
  // Boundary data are exterior primitive states at the actual half-edge midpoint.
  using BoundaryState = std::function<State(const Direction&, double, const State&, const Direction&, unsigned short)>;
  // NODE source values are sampled at primal nodes, not dual centroids.
  using Source = std::function<State(const Direction&, double)>;
  struct SpatialResult {
    std::vector<State> faces, volume, boundary, source, rhs;
  };
  // Serial entry point, outside OpenMP parallel regions. Rebuilds reconstruction
  // on every call, never changes W. LinSysRes receives faces+volume+boundary-V*S
  // only after all contributions have been validated; rhs=-LinSysRes/V.
  SpatialResult ComputeSpatialResidual(CGeometry* geometry, const CConfig* config, double time = 0,
                                      const BoundaryState& boundary = {}, const Source& source = {});
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
  void SetTime_Step(CGeometry*, CSolver**, CConfig*, unsigned short, unsigned long) override {
    EvolutionUnavailable();
  }
  void ExplicitEuler_Iteration(CGeometry*, CSolver**, CConfig*) override { EvolutionUnavailable(); }
  void ExplicitRK_Iteration(CGeometry*, CSolver**, CConfig*, unsigned short) override { EvolutionUnavailable(); }
  void ClassicalRK4_Iteration(CGeometry*, CSolver**, CConfig*, unsigned short) override { EvolutionUnavailable(); }
  void ImplicitEuler_Iteration(CGeometry*, CSolver**, CConfig*) override { EvolutionUnavailable(); }
};
