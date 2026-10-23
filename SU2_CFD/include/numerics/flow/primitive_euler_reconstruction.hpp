/*!
 * \brief Affine NODE reconstruction and primitive volume product.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "primitive_euler_2d.hpp"
#include <algorithm>

namespace PrimitiveEuler2D {
using Gradient = std::array<Direction, 4>;

/*! Same 0.99 positivity policy, with scaled division to avoid w-candidate overflow. */
inline double ReconstructionPositivityFactor(double value, double candidate) {
  if (!std::isfinite(value) || value <= 0 || !std::isfinite(candidate))
    throw std::invalid_argument("Invalid primitive positivity bound");
  const double floor = std::min(1e-12, 0.1 * value);
  if (candidate >= floor && candidate > 0) return 1;
  // No positive fractional floor is representable for the smallest subnormals.
  // A zero gradient retains the positive nodal value instead of rounding to zero.
  if (floor == 0) return 0;
  const double scale = std::max(value, std::abs(candidate));
  // Materialize divisions: native builds use -ffast-math, which could otherwise
  // factor the denominator back into the overflowing (value-candidate)/scale.
  const volatile double scaledValue = value / scale;
  const volatile double scaledCandidate = candidate / scale;
  return 0.99 * ((value - floor) / scale) / (scaledValue - scaledCandidate);
}

inline State Reconstruct(const State& state, const Gradient& gradient, const Direction& displacement) {
  ValidateDirection(displacement);
  State result = state;
  for (unsigned k = 0; k < 4; ++k) {
    ValidateDirection(gradient[k]);
    result[k] += gradient[k][0] * displacement[0] + gradient[k][1] * displacement[1];
  }
  return result;
}

/*! Coefficients use ALL components of the supplied mean (including 1/mean(rho)). */
inline State VolumeProduct(const State& mean, const Gradient& gradient, double gamma) {
  State dx{}, dy{};
  for (unsigned k = 0; k < 4; ++k) {
    ValidateDirection(gradient[k]);
    dx[k] = gradient[k][0];
    dy[k] = gradient[k][1];
  }
  auto result = Apply(PhysicalMatrix(mean, {{1, 0}}, gamma), dx);
  const auto y = Apply(PhysicalMatrix(mean, {{0, 1}}, gamma), dy);
  for (unsigned k = 0; k < 4; ++k) result[k] += y[k];
  return result;
}
}  // namespace PrimitiveEuler2D
