/*!
 * \brief Dedicated nodal state for primitive Euler, with no conservative-state cache.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "CScalarVariable.hpp"
#include "../numerics/flow/primitive_euler_2d.hpp"

class CPrimitiveEulerVariable final : public CScalarVariable {
  double gamma;

 public:
  // CScalarVariable supplies generic gradient/limiter arrays for a vector field;
  // it imposes no thermodynamic or conservative interpretation of Solution.
  CPrimitiveEulerVariable(const PrimitiveEuler2D::State& initial, unsigned long points, const CConfig* config);
  PrimitiveEuler2D::State GetPrimitiveState(unsigned long point) const;
  void SetPrimitiveState(unsigned long point, const PrimitiveEuler2D::State& state);
  PrimitiveEuler2D::State GetConservativeState(unsigned long point) const;

  su2double GetDensity(unsigned long i) const override { return Solution(i, 0); }
  su2double GetVelocity(unsigned long i, unsigned long d) const override { return Solution(i, d + 1); }
  su2double GetPressure(unsigned long i) const override { return Solution(i, 3); }
  su2double GetSoundSpeed(unsigned long i) const override {
    return std::sqrt(gamma * SU2_TYPE::GetValue(Solution(i, 3) / Solution(i, 0)));
  }
  // SU2's GetEnergy convention is specific total energy, NOT energy density.
  su2double GetEnergy(unsigned long i) const override {
    return Solution(i, 3) / ((gamma - 1) * Solution(i, 0)) +
           0.5 * (Solution(i, 1) * Solution(i, 1) + Solution(i, 2) * Solution(i, 2));
  }
  su2double GetEnthalpy(unsigned long i) const override { return GetEnergy(i) + GetPressure(i) / GetDensity(i); }
  const MatrixType& GetPrimitive() const override { return Solution; }
  su2double GetPrimitive(unsigned long i, unsigned long k) const override { return Solution(i, k); }
  su2double* GetPrimitive(unsigned long i) override { return Solution[i]; }
};
