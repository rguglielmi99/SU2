/*!
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "../../include/variables/CPrimitiveEulerVariable.hpp"

CPrimitiveEulerVariable::CPrimitiveEulerVariable(const PrimitiveEuler2D::State& initial, unsigned long points,
                                               const CConfig* config)
    : CScalarVariable(points, 2, 4, config), gamma(SU2_TYPE::GetValue(config->GetGamma())) {
  PrimitiveEuler2D::Validate(initial, gamma);
  for (unsigned long i = 0; i < nPoint; ++i) SetPrimitiveState(i, initial);
  Set_OldSolution();
  if (config->GetTime_Domain()) Set_Solution_time_n();
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY) Set_Solution_time_n1();
  Limiter = su2double(1.0);
}

PrimitiveEuler2D::State CPrimitiveEulerVariable::GetPrimitiveState(unsigned long point) const {
  if (point >= nPoint) throw std::out_of_range("Primitive Euler point index");
  PrimitiveEuler2D::State w{};
  for (unsigned k = 0; k < 4; ++k) w[k] = SU2_TYPE::GetValue(Solution(point, k));
  return w;
}

void CPrimitiveEulerVariable::SetPrimitiveState(unsigned long point, const PrimitiveEuler2D::State& state) {
  if (point >= nPoint) throw std::out_of_range("Primitive Euler point index");
  PrimitiveEuler2D::Validate(state, gamma);
  for (unsigned k = 0; k < 4; ++k) Solution(point, k) = state[k];
}

PrimitiveEuler2D::State CPrimitiveEulerVariable::GetConservativeState(unsigned long point) const {
  const auto w = GetPrimitiveState(point);
  PrimitiveEuler2D::Validate(w, gamma);
  return {{w[0], w[0] * w[1], w[0] * w[2], w[3] / (gamma - 1) + 0.5 * w[0] * (w[1] * w[1] + w[2] * w[2])}};
}
