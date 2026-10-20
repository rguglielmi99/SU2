/*!
 * \file primitive_euler_2d.hpp
 * \brief Euler 2D primitive-variable kernel for a straight path in (rho,u,v,p).
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */

#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

/*!
 * \brief Geometry-independent kernel for the initial, double-precision primal solver.
 * All states and matrix rows use (rho,u,v,p), NOT SU2's CEulerVariable primitive layout.
 * No entropy fix or positivity fallback is applied. Invalid inputs are rejected.
 * AD, boundary treatment, reconstruction and time integration are outside this kernel.
 */
namespace PrimitiveEuler2D {
using State = std::array<double, 4>;
using Direction = std::array<double, 2>;
using Matrix = std::array<State, 4>;

inline void Validate(const State& w, double gamma) {
  for (auto value : w)
    if (!std::isfinite(value)) throw std::invalid_argument("Primitive Euler: non-finite state");
  if (!(w[0] > 0.0) || !(w[3] > 0.0) || !std::isfinite(gamma) || !(gamma > 1.0))
    throw std::invalid_argument("Primitive Euler requires rho>0, p>0 and gamma>1");
}

inline void ValidateDirection(const Direction& n) {
  for (auto value : n)
    if (!std::isfinite(value)) throw std::invalid_argument("Primitive Euler: non-finite direction");
}

/*! \brief Exact straight-path average of 1/rho, with a stable equal-density limit. */
inline double PathInverseDensity(double left, double right) {
  if (!std::isfinite(left) || !std::isfinite(right) || !(left > 0.0) || !(right > 0.0))
    throw std::invalid_argument("Primitive Euler path requires positive finite densities");
  if (left == right) return 1.0 / left;
  // Sort to make endpoint reversal identical even in floating-point arithmetic.
  // lo stands for the lower density, hi for the higher
  const double lo = left < right ? left : right;
  const double hi = left < right ? right : left;
  const double delta = hi - lo;
  // The two function are equivalent on paper to compute log(hi/lo), but numerically if delta is very small computing log(hi)-log(lo) can be inaccurate so log(1+delta/lo) is better. The threshold is arbitrary but works well in practice.
  return delta < 0.5 * lo ? std::log1p(delta / lo) / delta : (std::log(hi) - std::log(lo)) / delta;
}

inline State Apply(const Matrix& a, const State& x) {
  // It is simply the product of the matrix and the vector, so a * x
  State result{};
  for (unsigned short i = 0; i < 4; ++i)
    for (unsigned short j = 0; j < 4; ++j) result[i] += a[i][j] * x[j];
  return result;
}

/*! \brief nx*Ax(W)+ny*Ay(W); direction need not be normalized (x/y axes give Ax/Ay). */
inline Matrix PhysicalMatrix(const State& w, const Direction& n, double gamma) {
  Validate(w, gamma);
  ValidateDirection(n);
  const double un = w[1] * n[0] + w[2] * n[1];
  return {{{un, w[0] * n[0], w[0] * n[1], 0.0},
           {0.0, un, 0.0, n[0] / w[0]},
           {0.0, 0.0, un, n[1] / w[0]},
           {0.0, gamma * w[3] * n[0], gamma * w[3] * n[1], un}}};
}

/*! \brief Exact integral of An(WL+s*(WR-WL)) over 0<=s<=1. */
inline Matrix PathMatrix(const State& left, const State& right, const Direction& n, double gamma) {
  Validate(left, gamma);
  Validate(right, gamma);
  State mean{};
  for (unsigned short k = 0; k < 4; ++k) mean[k] = 0.5 * left[k] + 0.5 * right[k];
  auto a = PhysicalMatrix(mean, n, gamma);
  const double beta = PathInverseDensity(left[0], right[0]);
  a[1][3] = beta * n[0];
  a[2][3] = beta * n[1];
  return a;
}

struct Spectrum {
  double advection;
  double acoustic;  // c*|n|; c for unit n, zero for n=0.
  State Eigenvalues() const { return {{advection - acoustic, advection, advection, advection + acoustic}}; }
};

inline Spectrum PathSpectrum(const State& left, const State& right, const Direction& n, double gamma) {
  Validate(left, gamma);
  Validate(right, gamma);
  ValidateDirection(n);
  const double pressure = 0.5 * left[3] + 0.5 * right[3];
  return {(0.5 * left[1] + 0.5 * right[1]) * n[0] + (0.5 * left[2] + 0.5 * right[2]) * n[1],
          std::sqrt(gamma * pressure * PathInverseDensity(left[0], right[0])) * std::hypot(n[0], n[1])};
}

/*!
 * \brief Spectral absolute value, using characteristic projectors including the double eigenvalue.
 * K=(A-un*I)/(c*|n|), P0=I-K^2, P+=(K^2+K)/2, P-=(K^2-K)/2.
 * This avoids inverting a matrix of eigenvectors at sonic or zero-advection states.
 */
inline Matrix AbsolutePathMatrix(const State& left, const State& right, const Direction& n, double gamma) {
  const auto a = PathMatrix(left, right, n, gamma);
  const auto spectrum = PathSpectrum(left, right, n, gamma);
  if (spectrum.acoustic == 0.0) return Matrix{};
  const double zero = std::abs(spectrum.advection);
  const double plus = std::abs(spectrum.advection + spectrum.acoustic);
  const double minus = std::abs(spectrum.advection - spectrum.acoustic);
  Matrix k = a;
  for (unsigned short i = 0; i < 4; ++i) {
    k[i][i] -= spectrum.advection;
    for (auto& value : k[i]) value /= spectrum.acoustic;
  }
  Matrix result{};
  for (unsigned short i = 0; i < 4; ++i)
    for (unsigned short j = 0; j < 4; ++j) {
      double k2 = 0.0;
      for (unsigned short m = 0; m < 4; ++m) k2 += k[i][m] * k[m][j];
      result[i][j] = (i == j ? zero : 0.0) + 0.5 * (plus - minus) * k[i][j] +
                     (0.5 * plus + 0.5 * minus - zero) * k2;
    }
  return result;
}

struct Fluctuations {
  State minus{};
  State plus{};
};

/*!
 * \brief D-/D+ multiplied by face measure exactly once, from an integrated normal S.
 * The caller adds minus to node i and plus to node j, BOTH with positive signs.
 * Zero-area faces return zero; no volume term or division by cell volume is included.
 */
inline Fluctuations ComputeFluctuations(const State& left, const State& right, const Direction& s, double gamma) {
  Validate(left, gamma);
  Validate(right, gamma);
  ValidateDirection(s);
  const double area = std::hypot(s[0], s[1]);
  if (area == 0.0) return {};
  const Direction n{{s[0] / area, s[1] / area}};
  State jump{};
  for (unsigned short k = 0; k < 4; ++k) jump[k] = right[k] - left[k];
  const auto path = Apply(PathMatrix(left, right, n, gamma), jump);
  const auto dissipation = Apply(AbsolutePathMatrix(left, right, n, gamma), jump);
  Fluctuations result;
  for (unsigned short k = 0; k < 4; ++k) {
    result.minus[k] = 0.5 * area * (path[k] - dissipation[k]);
    result.plus[k] = 0.5 * area * (path[k] + dissipation[k]);
  }
  return result;
}
}  // namespace PrimitiveEuler2D
