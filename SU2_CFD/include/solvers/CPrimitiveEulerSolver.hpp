/*!
 * \brief Phase-3 primitive Euler storage/factory scaffold; evolution is deliberately unavailable.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "CSolver.hpp"
#include "../variables/CPrimitiveEulerVariable.hpp"
#include "../geometry/primitive_euler_geometry.hpp"

class CPrimitiveEulerSolver final : public CSolver {
  std::unique_ptr<CPrimitiveEulerVariable> nodes;
  std::unique_ptr<PrimitiveEuler2D::CAggregatedGeometry> metrics;
  CVariable* GetBaseClassPointerToNodes() override { return nodes.get(); }

 public:
  CPrimitiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short mesh);
  CPrimitiveEulerVariable& GetPrimitiveNodes() { return *nodes; }
  const CPrimitiveEulerVariable& GetPrimitiveNodes() const { return *nodes; }
  const PrimitiveEuler2D::CAggregatedGeometry& GetGeometryMetrics() const { return *metrics; }
  static void ValidateConfiguration(const CConfig& config);
  static void EvolutionUnavailable();

  void Upwind_Residual(CGeometry*, CSolver**, CNumerics**, CConfig*, unsigned short) override {
    EvolutionUnavailable();
  }
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
