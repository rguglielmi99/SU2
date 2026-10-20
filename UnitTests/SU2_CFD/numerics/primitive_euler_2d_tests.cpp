/*!
 * \file primitive_euler_2d_tests.cpp
 * \brief Straight-path primitive Euler tests, independent of solver/mesh setup.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "../../../SU2_CFD/include/numerics/flow/primitive_euler_2d.hpp"
#include <limits>
#include <algorithm>
#include <random>
#include <vector>

using namespace PrimitiveEuler2D;

namespace {
constexpr double gammaGas = 1.4;

void CheckState(const State& actual, const State& expected, double tolerance = 2e-12) {
  for (unsigned short k = 0; k < 4; ++k) {
    INFO("component " << k);
    CHECK(actual[k] == Approx(expected[k]).margin(tolerance).epsilon(tolerance));
  }
}

State Difference(const State& right, const State& left) {
  State result{};
  for (unsigned short k = 0; k < 4; ++k) result[k] = right[k] - left[k];
  return result;
}

// Mixed tolerances: an absolute floor handles exact zeros; the relative part
// uses the magnitude of the terms BEFORE subtraction for cancellation-prone
// identities. These tolerances are for moderate double-precision scales, not
// arbitrary overflow/underflow inputs. Quadrature has a separate tighter bound.
void CheckScaled(double actual, double expected, double scale, double relative = 3e-11) {
  REQUIRE(std::isfinite(actual));
  REQUIRE(std::isfinite(expected));
  CHECK(std::abs(actual - expected) <= 2e-12 + relative * std::max(std::abs(expected), scale));
}

using LongState = std::array<long double, 4>;
using LongMatrix = std::array<LongState, 4>;

// Independent Gauss-Legendre integration of the explicit PDE coefficients.
// Neither PhysicalMatrix nor beta/PathMatrix/PathSpectrum is used as an oracle.
// Geometric panels in density resolve the thin endpoint layer of 1/rho for
// large density ratios; density changes by at most a factor 1.5 per panel.
// Eight Gauss points then integrate that smooth panel to near roundoff.
LongMatrix IntegratePath(const State& left, const State& right, const Direction& n, double gamma) {
  const long double points[] = {0.18343464249564980494L, 0.52553240991632898582L,
                                0.79666647741362673959L, 0.96028985649753623168L};
  const long double weights[] = {0.36268378337836198297L, 0.31370664587788728734L,
                                 0.22238103445337447054L, 0.10122853629037625915L};
  const long double lo = std::min(left[0], right[0]), hi = std::max(left[0], right[0]);
  std::vector<long double> panels{0.0L};
  for (long double rho = lo * 1.5L; rho < hi; rho *= 1.5L) panels.push_back((rho - lo) / (hi - lo));
  panels.push_back(1.0L);
  LongMatrix integral{};
  for (std::size_t panel = 1; panel < panels.size(); ++panel) {
    const long double mid = (panels[panel - 1] + panels[panel]) / 2;
    const long double half = (panels[panel] - panels[panel - 1]) / 2;
    for (unsigned q = 0; q < 4; ++q)
      for (int sign : {-1, 1}) {
        long double s = mid + sign * half * points[q];
        if (left[0] > right[0]) s = 1 - s;
        LongState w{};
        // Convex form avoids cancellation near the low-density endpoint.
        for (unsigned k = 0; k < 4; ++k) w[k] = (1 - s) * left[k] + s * right[k];
        const long double un = w[1] * n[0] + w[2] * n[1];
        const LongMatrix a{{{{un, w[0] * n[0], w[0] * n[1], 0}},
                             {{0, un, 0, n[0] / w[0]}},
                             {{0, 0, un, n[1] / w[0]}},
                             {{0, gamma * w[3] * n[0], gamma * w[3] * n[1], un}}}};
        for (unsigned i = 0; i < 4; ++i)
          for (unsigned j = 0; j < 4; ++j) integral[i][j] += half * weights[q] * a[i][j];
      }
  }
  return integral;
}

State Sum(const State& a, const State& b) {
  State result{};
  for (unsigned k = 0; k < 4; ++k) result[k] = a[k] + b[k];
  return result;
}

// Checks all identities using an integrated normal S, including non-unit S.
void CheckKernelCase(const State& left, const State& right, const Direction& s, double gamma) {
  const auto integral = IntegratePath(left, right, s, gamma);
  const auto a = PathMatrix(left, right, s, gamma);
  const auto absolute = AbsolutePathMatrix(left, right, s, gamma);
  const auto jump = Difference(right, left);
  const auto d = ComputeFluctuations(left, right, s, gamma);
  const auto reverse = ComputeFluctuations(right, left, {{-s[0], -s[1]}}, gamma);
  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      const double reference = static_cast<double>(integral[i][j]);
      CheckScaled(a[i][j], reference, std::abs(reference), 5e-13);
      double a2 = 0, abs2 = 0, scale = 0;
      for (unsigned k = 0; k < 4; ++k) {
        a2 += a[i][k] * a[k][j];
        abs2 += absolute[i][k] * absolute[k][j];
        scale += std::abs(a[i][k] * a[k][j]) + std::abs(absolute[i][k] * absolute[k][j]);
      }
      CheckScaled(abs2, a2, scale);
    }
    long double reference = 0;
    double action = 0, actionScale = 0;
    for (unsigned j = 0; j < 4; ++j) {
      reference += integral[i][j] * (static_cast<long double>(right[j]) - left[j]);
      action += absolute[i][j] * jump[j];
      actionScale += std::abs(absolute[i][j] * jump[j]) + std::abs(a[i][j] * jump[j]);
    }
    const double scale = actionScale + std::abs(d.minus[i]) + std::abs(d.plus[i]);
    CheckScaled(d.minus[i] + d.plus[i], static_cast<double>(reference), scale);
    CheckScaled(d.plus[i] - d.minus[i], action, scale);
    CheckScaled(reverse.minus[i], d.plus[i], scale);
    CheckScaled(reverse.plus[i], d.minus[i], scale);
  }
  for (double alpha : {0.125, 2.75}) {
    const auto scaled = ComputeFluctuations(left, right, {{alpha * s[0], alpha * s[1]}}, gamma);
    for (unsigned i = 0; i < 4; ++i) {
      const double scale = alpha * (std::abs(d.minus[i]) + std::abs(d.plus[i]));
      CheckScaled(scaled.minus[i], alpha * d.minus[i], scale);
      CheckScaled(scaled.plus[i], alpha * d.plus[i], scale);
    }
  }
}

// Independent local physical eigensystem (beta=1/rho, c=sqrt(gamma*p/rho)).
// Tests the sign of |lambda| as well as its square, including the double mode.
void CheckLocalSpectrum(const State& w, const Direction& s, double gamma) {
  const auto a = PathMatrix(w, w, s, gamma);
  const auto physical = PhysicalMatrix(w, s, gamma);
  const auto absolute = AbsolutePathMatrix(w, w, s, gamma);
  const double length = std::hypot(s[0], s[1]);
  const double nx = s[0] / length, ny = s[1] / length;
  const double c = std::sqrt(gamma * w[3] / w[0]);
  const double un = w[1] * s[0] + w[2] * s[1];
  const State lambda{{un - c * length, un, un, un + c * length}};
  const std::array<State, 4> vectors{{{{w[0] / c, -nx, -ny, w[0] * c}},
                                      {{1, 0, 0, 0}}, {{0, -ny, nx, 0}},
                                      {{w[0] / c, nx, ny, w[0] * c}}}};
  const auto reported = PathSpectrum(w, w, s, gamma).Eigenvalues();
  for (unsigned i = 0; i < 4; ++i) {
    CheckScaled(reported[i], lambda[i], std::abs(un) + c * length);
    for (unsigned j = 0; j < 4; ++j) CheckScaled(a[i][j], physical[i][j], std::abs(physical[i][j]));
    for (unsigned field = 0; field < 4; ++field) {
      const auto value = Apply(absolute, vectors[field]);
      double scale = 0;
      for (unsigned j = 0; j < 4; ++j) scale += std::abs(absolute[i][j] * vectors[field][j]);
      CheckScaled(value[i], std::abs(lambda[field]) * vectors[field][i], scale);
    }
  }
  const auto constant = ComputeFluctuations(w, w, s, gamma);
  CheckState(constant.minus, State{});
  CheckState(constant.plus, State{});
}

using Gradient = std::array<Direction, 4>;

// Synthetic affine data only: this is NOT a MUSCL/geometry implementation.
State Reconstruct(const State& w, const Gradient& g, const Direction& origin, const Direction& q) {
  State result = w;
  for (unsigned k = 0; k < 4; ++k)
    result[k] += g[k][0] * (q[0] - origin[0]) + g[k][1] * (q[1] - origin[1]);
  return result;
}

double MaxDifference(const State& a, const State& b) {
  double result = 0;
  for (unsigned k = 0; k < 4; ++k) result = std::max(result, std::abs(a[k] - b[k]));
  return result;
}
}  // namespace

TEST_CASE("Primitive Euler matrices reproduce the four PDE components", "[primitive-euler]") {
  const State w{{1.3, 0.7, -0.2, 2.1}};
  const State gx{{0.12, -0.3, 0.4, -0.5}}, gy{{-0.13, 0.2, -0.6, 0.8}};
  auto actual = Apply(PhysicalMatrix(w, {{1.0, 0.0}}, gammaGas), gx);
  const auto y = Apply(PhysicalMatrix(w, {{0.0, 1.0}}, gammaGas), gy);
  for (unsigned short k = 0; k < 4; ++k) actual[k] += y[k];
  const State expected{{w[1] * gx[0] + w[2] * gy[0] + w[0] * (gx[1] + gy[2]),
                        w[1] * gx[1] + w[2] * gy[1] + gx[3] / w[0],
                        w[1] * gx[2] + w[2] * gy[2] + gy[3] / w[0],
                        w[1] * gx[3] + w[2] * gy[3] + gammaGas * w[3] * (gx[1] + gy[2])}};
  CheckState(actual, expected);
}

TEST_CASE("Straight-path inverse density is stable and symmetric", "[primitive-euler]") {
  CHECK(PathInverseDensity(2.0, 2.0) == 0.5);
  for (double relative : {1e-14, 1e-8, 0.1, 1.0, 1e6}) {
    const double left = 2.0, right = left * (1.0 + relative);
    // Long-double reference avoids cancellation in the oracle near equal densities.
    const long double delta = static_cast<long double>(right) - left;
    const double expected = static_cast<double>(std::log1p(delta / left) / delta);
    CHECK(PathInverseDensity(left, right) == Approx(expected).epsilon(2e-14));
    CHECK(PathInverseDensity(left, right) == PathInverseDensity(right, left));
  }
  const double next = std::nextafter(2.0, 3.0);
  CHECK(PathInverseDensity(2.0, next) == Approx(0.5).epsilon(1e-15));
}

TEST_CASE("Path matrix agrees with independent numerical integration", "[primitive-euler]") {
  const State left{{0.7, -0.4, 0.8, 1.1}}, right{{2.4, 1.2, -0.3, 2.2}};
  const Direction n{{0.6, 0.8}};
  Matrix integral{};
  constexpr unsigned steps = 4096;
  // Composite Simpson quadrature of the physical matrix along the path.
  for (unsigned q = 0; q <= steps; ++q) {
    const double s = double(q) / steps;
    State w{};
    for (unsigned k = 0; k < 4; ++k) w[k] = left[k] + s * (right[k] - left[k]);
    const auto a = PhysicalMatrix(w, n, gammaGas);
    const double weight = (q == 0 || q == steps) ? 1.0 : (q % 2 ? 4.0 : 2.0);
    for (unsigned i = 0; i < 4; ++i)
      for (unsigned j = 0; j < 4; ++j) integral[i][j] += weight * a[i][j] / (3.0 * steps);
  }
  const auto path = PathMatrix(left, right, n, gammaGas);
  for (unsigned i = 0; i < 4; ++i) CheckState(path[i], integral[i]);
}

TEST_CASE("Absolute Roe matrix acts correctly on all characteristic fields", "[primitive-euler]") {
  const Direction n{{0.6, 0.8}};
  for (double un : {-4.0, -1.0, 0.0, 1.0, 4.0}) {
    const State left{{0.8, un * n[0], un * n[1], 1.0}};
    const State right{{1.7, un * n[0], un * n[1], 1.4}};
    const double rho = 0.5 * (left[0] + right[0]);
    const double beta = PathInverseDensity(left[0], right[0]);
    const double c = PathSpectrum(left, right, n, gammaGas).acoustic;
    const auto a = PathMatrix(left, right, n, gammaGas);
    const auto absolute = AbsolutePathMatrix(left, right, n, gammaGas);
    // Two independent vectors for the repeated advective eigenvalue.
    const std::array<State, 4> eigenvectors{{
        {{rho / c, -n[0], -n[1], c / beta}},
        {{1.0, 0.0, 0.0, 0.0}},
        {{0.0, -n[1], n[0], 0.0}},
        {{rho / c, n[0], n[1], c / beta}}}};
    const auto eigenvalues = PathSpectrum(left, right, n, gammaGas).Eigenvalues();
    for (unsigned field = 0; field < 4; ++field) {
      State expected{}, expectedAbs{};
      for (unsigned k = 0; k < 4; ++k) {
        expected[k] = eigenvalues[field] * eigenvectors[field][k];
        expectedAbs[k] = std::abs(eigenvalues[field]) * eigenvectors[field][k];
      }
      CheckState(Apply(a, eigenvectors[field]), expected);
      CheckState(Apply(absolute, eigenvectors[field]), expectedAbs);
    }
  }
}

TEST_CASE("Fluctuations satisfy path consistency reversal and area scaling", "[primitive-euler]") {
  const State left{{0.7, -0.4, 0.8, 1.1}}, right{{2.4, 1.2, -0.3, 2.2}};
  for (const Direction s : {Direction{{1.0, 0.0}}, Direction{{0.0, -2.0}}, Direction{{1.8, 2.4}}}) {
    const auto d = ComputeFluctuations(left, right, s, gammaGas);
    const auto reversed = ComputeFluctuations(right, left, {{-s[0], -s[1]}}, gammaGas);
    const auto scaled = ComputeFluctuations(left, right, {{3.0 * s[0], 3.0 * s[1]}}, gammaGas);
    const auto expected = Apply(PathMatrix(left, right, s, gammaGas), Difference(right, left));
    State sum{}, scaledMinus{}, scaledPlus{};
    for (unsigned k = 0; k < 4; ++k) {
      sum[k] = d.minus[k] + d.plus[k];
      scaledMinus[k] = 3.0 * d.minus[k];
      scaledPlus[k] = 3.0 * d.plus[k];
    }
    CheckState(sum, expected);
    CheckState(reversed.minus, d.plus);
    CheckState(reversed.plus, d.minus);
    CheckState(scaled.minus, scaledMinus);
    CheckState(scaled.plus, scaledPlus);
  }
}

TEST_CASE("Constant states zero-area faces sonic and one-way propagation", "[primitive-euler]") {
  const State w{{1.0, 0.0, 0.0, 1.0}};
  const State zero{};
  CheckState(ComputeFluctuations(w, w, {{2.0, 3.0}}, gammaGas).minus, zero);
  CheckState(ComputeFluctuations(w, w, {{2.0, 3.0}}, gammaGas).plus, zero);
  const State other{{1.2, 0.3, 0.2, 1.5}};
  const auto degenerate = ComputeFluctuations(w, other, {{0.0, 0.0}}, gammaGas);
  CheckState(degenerate.minus, zero);
  CheckState(degenerate.plus, zero);
  for (double speed : {-5.0, 5.0}) {
    State left = w, right = other;
    left[1] = right[1] = speed;
    const auto d = ComputeFluctuations(left, right, {{1.0, 0.0}}, gammaGas);
    CheckState(speed > 0 ? d.minus : d.plus, zero);
  }
  // Exact sonic eigenvalue, including the double zero eigenvalue at rest.
  for (double speed : {-std::sqrt(gammaGas), 0.0, std::sqrt(gammaGas)}) {
    State sonic = w;
    sonic[1] = speed;
    const auto a = PathMatrix(sonic, sonic, {{1.0, 0.0}}, gammaGas);
    const auto absolute = AbsolutePathMatrix(sonic, sonic, {{1.0, 0.0}}, gammaGas);
    // Spectral absolute value must square to A^2.
    for (unsigned j = 0; j < 4; ++j) {
      State basis{};
      basis[j] = 1.0;
      CheckState(Apply(absolute, Apply(absolute, basis)), Apply(a, Apply(a, basis)));
    }
  }
}

TEST_CASE("Primitive Euler rejects nonphysical and nonfinite inputs", "[primitive-euler]") {
  const State good{{1.0, 0.0, 0.0, 1.0}};
  for (unsigned component : {0u, 3u}) {
    State bad = good;
    bad[component] = 0.0;
    CHECK_THROWS_AS(ComputeFluctuations(bad, good, {{1.0, 0.0}}, gammaGas), std::invalid_argument);
    bad[component] = -1.0;
    CHECK_THROWS_AS(PathMatrix(good, bad, {{1.0, 0.0}}, gammaGas), std::invalid_argument);
  }
  CHECK_THROWS_AS(PhysicalMatrix(good, {{1.0, 0.0}}, 1.0), std::invalid_argument);
  CHECK_THROWS_AS(PathInverseDensity(0.0, 1.0), std::invalid_argument);
  State bad = good;
  bad[1] = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_AS(PathSpectrum(bad, good, {{1.0, 0.0}}, gammaGas), std::invalid_argument);
  CHECK_THROWS_AS(ComputeFluctuations(good, good, {{INFINITY, 0.0}}, gammaGas), std::invalid_argument);
}

TEST_CASE("Equal-state matrices and stationary contact retain their physical limits", "[primitive-euler]") {
  const Direction n{{0.6, 0.8}};
  const State left{{0.8, -0.8, 0.6, 1.2}}, right{{1.7, -1.6, 1.2, 1.2}};
  const auto physical = PhysicalMatrix(left, n, gammaGas);
  const auto path = PathMatrix(left, left, n, gammaGas);
  for (unsigned i = 0; i < 4; ++i) CheckState(path[i], physical[i]);
  // Density and tangential velocity may jump; pressure and normal velocity agree.
  const auto contact = ComputeFluctuations(left, right, n, gammaGas);
  CheckState(contact.minus, State{});
  CheckState(contact.plus, State{});
}

TEST_CASE("Primitive fluctuations are covariant under rotation", "[primitive-euler]") {
  const State left{{0.7, -0.4, 0.8, 1.1}}, right{{2.4, 1.2, -0.3, 2.2}};
  const Direction s{{1.8, 2.4}};
  const double c = std::cos(0.73), t = std::sin(0.73);
  auto rotateState = [c, t](State w) {
    const double u = w[1], v = w[2];
    w[1] = c * u - t * v;
    w[2] = t * u + c * v;
    return w;
  };
  const Direction rotatedS{{c * s[0] - t * s[1], t * s[0] + c * s[1]}};
  const auto d = ComputeFluctuations(left, right, s, gammaGas);
  const auto rotated = ComputeFluctuations(rotateState(left), rotateState(right), rotatedS, gammaGas);
  CheckState(rotated.minus, rotateState(d.minus));
  CheckState(rotated.plus, rotateState(d.plus));
}

TEST_CASE("Independent path quadrature resolves large density ratios", "[primitive-euler]") {
  // Calibration of the numerical oracle, including both endpoint orientations.
  // The analytic expression is used ONLY here, not in the many-state checks.
  for (double density : {1e-3, 0.7, 1.0, 1.00000000000001, 1e3}) {
    for (bool reversed : {false, true}) {
      const State left{{reversed ? density : 1.0, 0.2, -0.3, 1.4}};
      const State right{{reversed ? 1.0 : density, -0.5, 0.6, 0.8}};
      const auto integral = IntegratePath(left, right, {{1.0, 0.0}}, gammaGas);
      const long double delta = static_cast<long double>(density) - 1;
      const double exact = density == 1.0 ? 1.0 : static_cast<double>(std::log1p(delta) / delta);
      CheckScaled(static_cast<double>(integral[1][3]), exact, std::abs(exact), 5e-14);
      CheckScaled(static_cast<double>(integral[0][1]), 0.5 * (left[0] + right[0]),
                  0.5 * (left[0] + right[0]), 5e-14);
    }
  }
}

TEST_CASE("Many physical paths obey consistency spectral and orientation identities", "[primitive-euler]") {
  // Raw mt19937 outputs, not uniform_real_distribution, make the corpus
  // reproducible across standard-library implementations. Seed: 0x5EED1234.
  std::mt19937 generator(0x5EED1234u);
  auto unit = [&generator]() { return double(generator()) / 4294967296.0; };
  const std::array<Direction, 8> normals{{{{1, 0}}, {{-2, 0}}, {{0, 0.25}}, {{0, -3}},
                                          {{0.6, 0.8}}, {{-1.7, 0.4}}, {{0.3, -2.1}}, {{-0.9, -1.3}}}};
  const double gammas[] = {1.2, 1.4, 5.0 / 3.0};
  // 48 pairs x 8 normals: quasi-coincident, moderate and six-decade jumps.
  for (unsigned sample = 0; sample < 48; ++sample) {
    State left{{std::pow(10.0, -2 + 4 * unit()), 4 * unit() - 2, 4 * unit() - 2,
                std::pow(10.0, -2 + 4 * unit())}};
    State right{{std::pow(10.0, -2 + 4 * unit()), 4 * unit() - 2, 4 * unit() - 2,
                 std::pow(10.0, -2 + 4 * unit())}};
    if (sample % 4 == 0)
      for (unsigned k = 0; k < 4; ++k) right[k] = left[k] * (1 + (2 * unit() - 1) * 1e-10);
    if (sample % 4 == 1) {
      left[0] = 1e-3;
      right[0] = 1e3;
      left[3] = 1e3;
      right[3] = 1e-3;
      if (sample % 8 == 1) std::swap(left, right);
    }
    if (sample % 4 == 2) left[1] = left[2] = right[1] = right[2] = 0;
    const double gamma = gammas[sample % 3];
    for (unsigned normal = 0; normal < normals.size(); ++normal) {
      INFO("sample=" << sample << " normal=" << normal << " gamma=" << gamma);
      CheckKernelCase(left, right, normals[normal], gamma);
      CheckLocalSpectrum(left, normals[normal], gamma);
    }
    // Exercise the explicit zero-area branch once per physical pair.
    CheckKernelCase(left, right, {{0, 0}}, gamma);
  }
}

TEST_CASE("Splitting remains accurate on both sides of sonic thresholds", "[primitive-euler]") {
  const std::array<Direction, 3> normals{{{{-2, 0}}, {{0, -0.5}}, {{-1.2, 1.6}}}};
  for (const auto& s : normals) {
    const double length = std::hypot(s[0], s[1]);
    const double nx = s[0] / length, ny = s[1] / length;
    for (double mach : {-3.0, -1.0 - 1e-10, -1.0, -1.0 + 1e-10, 0.0,
                         1.0 - 1e-10, 1.0, 1.0 + 1e-10, 3.0}) {
      INFO("normal=" << s[0] << "," << s[1] << " normal Mach=" << mach);
      State left{{0.6, 0, 0, 1.7}}, right{{2.3, 0, 0, 0.8}};
      const double beta = static_cast<double>(IntegratePath(left, right, {{1, 0}}, gammaGas)[1][3]);
      const double cPath = std::sqrt(gammaGas * 0.5 * (left[3] + right[3]) * beta);
      // Same normal velocity on either side, opposite tangential components.
      left[1] = mach * cPath * nx - 0.3 * ny;
      left[2] = mach * cPath * ny + 0.3 * nx;
      right[1] = mach * cPath * nx + 0.3 * ny;
      right[2] = mach * cPath * ny - 0.3 * nx;
      CheckKernelCase(left, right, s, gammaGas);
      State local{{1.2, 0, 0, 0.9}};
      const double cLocal = std::sqrt(gammaGas * local[3] / local[0]);
      local[1] = mach * cLocal * nx;
      local[2] = mach * cLocal * ny;
      CheckLocalSpectrum(local, s, gammaGas);
    }
  }
}

TEST_CASE("Two dual segments and their aggregate share path sum but not splitting", "[primitive-euler]") {
  // Broken dual face C1 -> M -> C2, S=(dy,-dx) on each segment.
  // M is the midpoint of the primal edge xi=(0,0), xj=(2,0).
  const Direction c1{{0.7, -0.6}}, mid{{1, 0}}, c2{{1.5, 0.8}};
  const Direction s1{{mid[1] - c1[1], c1[0] - mid[0]}};
  const Direction s2{{c2[1] - mid[1], mid[0] - c2[0]}};
  const Direction effective{{s1[0] + s2[0], s1[1] + s2[1]}};
  const State left{{0.9, 0.2, -0.4, 1.1}}, right{{1.6, -0.3, 0.5, 1.8}};
  const auto d1 = ComputeFluctuations(left, right, s1, gammaGas);
  const auto d2 = ComputeFluctuations(left, right, s2, gammaGas);
  const auto aggregated = ComputeFluctuations(left, right, effective, gammaGas);
  CheckKernelCase(left, right, s1, gammaGas);
  CheckKernelCase(left, right, s2, gammaGas);
  CheckKernelCase(left, right, effective, gammaGas);
  const auto minus = Sum(d1.minus, d2.minus), plus = Sum(d1.plus, d2.plus);
  CheckState(Sum(minus, plus), Sum(aggregated.minus, aggregated.plus));
  // This is a deliberately resolved difference, many orders above tolerance.
  CHECK(MaxDifference(minus, aggregated.minus) > 1e-4);
  CHECK(MaxDifference(plus, aggregated.plus) > 1e-4);
}

TEST_CASE("Synthetic varying segment states distinguish both quadratures", "[primitive-euler]") {
  const Direction xi{{0, 0}}, xj{{2, 0}}, mid{{1, 0}}, c1{{0.7, -0.6}}, c2{{1.5, 0.8}};
  const Direction q1{{0.5 * (c1[0] + mid[0]), 0.5 * (c1[1] + mid[1])}};
  const Direction q2{{0.5 * (c2[0] + mid[0]), 0.5 * (c2[1] + mid[1])}};
  // AGGREGATED uses the primal edge midpoint explicitly, not a weighted point.
  const Direction qEffective = mid;
  const Direction s1{{0.6, -0.3}}, s2{{0.8, -0.5}}, effective{{1.4, -0.8}};
  const State wi{{1.1, 0.2, -0.1, 1.4}}, wj{{1.5, -0.3, 0.4, 1.8}};
  const Gradient gi{{{{0.2, 0.1}}, {{-0.1, 0.3}}, {{0.05, -0.2}}, {{0.3, 0.2}}}};
  const Gradient gj{{{{-0.1, 0.2}}, {{0.2, -0.1}}, {{-0.15, 0.1}}, {{-0.2, 0.1}}}};
  const auto l1 = Reconstruct(wi, gi, xi, q1), r1 = Reconstruct(wj, gj, xj, q1);
  const auto l2 = Reconstruct(wi, gi, xi, q2), r2 = Reconstruct(wj, gj, xj, q2);
  const auto le = Reconstruct(wi, gi, xi, qEffective), re = Reconstruct(wj, gj, xj, qEffective);
  for (const auto& w : {l1, r1, l2, r2, le, re}) {
    REQUIRE(w[0] > 0);
    REQUIRE(w[3] > 0);
  }
  CHECK(MaxDifference(l1, l2) > 1e-3);
  CHECK(MaxDifference(r1, r2) > 1e-3);
  CheckKernelCase(l1, r1, s1, gammaGas);
  CheckKernelCase(l2, r2, s2, gammaGas);
  CheckKernelCase(le, re, effective, gammaGas);
  const auto d1 = ComputeFluctuations(l1, r1, s1, gammaGas);
  const auto d2 = ComputeFluctuations(l2, r2, s2, gammaGas);
  const auto dEffective = ComputeFluctuations(le, re, effective, gammaGas);
  const auto minus = Sum(d1.minus, d2.minus), plus = Sum(d1.plus, d2.plus);
  CHECK(MaxDifference(minus, dEffective.minus) > 1e-4);
  CHECK(MaxDifference(plus, dEffective.plus) > 1e-4);
  // With changing reconstructed states, even the total path sums differ.
  CHECK(MaxDifference(Sum(minus, plus), Sum(dEffective.minus, dEffective.plus)) > 1e-4);
}

TEST_CASE("An exact continuous affine field has zero jumps on both segments", "[primitive-euler]") {
  const Direction xi{{0, 0}}, xj{{2, 0}};
  const State origin{{1.2, 0.3, -0.2, 1.5}};
  const Gradient g{{{{0.1, 0.04}}, {{-0.2, 0.3}}, {{0.15, -0.1}}, {{0.2, -0.05}}}};
  const auto wi = Reconstruct(origin, g, {{0, 0}}, xi);
  const auto wj = Reconstruct(origin, g, {{0, 0}}, xj);
  const std::array<Direction, 3> points{{{{0.85, -0.3}}, {{1.25, 0.4}}, {{1, 0}}}};
  const std::array<Direction, 3> normals{{{{0.6, -0.3}}, {{0.8, -0.5}}, {{1.4, -0.8}}}};
  for (unsigned m = 0; m < points.size(); ++m) {
    // Direct analytical field value is independent of either reconstruction.
    State exact{};
    const double x = points[m][0], y = points[m][1];
    exact = {{1.2 + 0.1 * x + 0.04 * y, 0.3 - 0.2 * x + 0.3 * y,
              -0.2 + 0.15 * x - 0.1 * y, 1.5 + 0.2 * x - 0.05 * y}};
    REQUIRE(exact[0] > 0);
    REQUIRE(exact[3] > 0);
    const auto left = Reconstruct(wi, g, xi, points[m]), right = Reconstruct(wj, g, xj, points[m]);
    CheckState(left, exact);
    CheckState(right, exact);
    CheckState(left, right);
    const auto d = ComputeFluctuations(left, right, normals[m], gammaGas);
    CheckState(d.minus, State{});
    CheckState(d.plus, State{});
  }
}
