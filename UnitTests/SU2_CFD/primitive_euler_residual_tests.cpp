/*!
 * \brief Phase-4 primitive spatial operator tests on native SU2 meshes.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include <limits>

using namespace PrimitiveEuler2D;
using namespace PrimitiveEulerTest;

namespace {
const std::vector<Direction> xy{{{0, 0}}, {{2, 0}}, {{2, 1}}, {{0, 1}}, {{0.8, 0.3}}};
const std::vector<Triangle> triangles{{{0, 1, 4}}, {{1, 2, 4}}, {{2, 3, 4}}, {{3, 0, 4}}};
const std::string setup = "TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\n"
                          "NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\nSLOPE_LIMITER_FLOW=NONE\n";
const Gradient affineGradient{{{{0.2, -0.1}}, {{0.12, 0.07}}, {{-0.08, 0.09}}, {{0.3, 0.11}}}};

void Check(const State& actual, const State& expected, double margin = 2e-12) {
  for (unsigned k = 0; k < 4; ++k) CHECK(actual[k] == Approx(expected[k]).margin(margin));
}

State Affine(const Direction& x) {
  State w{{2, 0.5, 0.3, 3}};
  for (unsigned k = 0; k < 4; ++k) w[k] += affineGradient[k][0] * x[0] + affineGradient[k][1] * x[1];
  return w;
}

State ExpandedVolume(const State& w, const Gradient& g, double gamma) {
  return {{w[1] * g[0][0] + w[2] * g[0][1] + w[0] * (g[1][0] + g[2][1]),
           w[1] * g[1][0] + w[2] * g[1][1] + g[3][0] / w[0],
           w[1] * g[2][0] + w[2] * g[2][1] + g[3][1] / w[0],
           w[1] * g[3][0] + w[2] * g[3][1] + gamma * w[3] * (g[1][0] + g[2][1])}};
}

const CPrimitiveEulerSolver::BoundaryState extrapolate =
    [](const Direction&, double, const State& interior, const Direction&, unsigned short) { return interior; };
const CPrimitiveEulerSolver::BoundaryState exactAffine =
    [](const Direction& x, double, const State&, const Direction&, unsigned short) { return Affine(x); };

void SetAffine(CPrimitiveEulerSolver& solver) {
  for (unsigned long i = 0; i < xy.size(); ++i) solver.GetPrimitiveNodes().SetPrimitiveState(i, Affine(xy[i]));
}
}  // namespace

TEST_CASE("Boundary half-faces close the native dual cells and integrate affine data", "[primitive-solver][primitive-residual]") {
  MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup);
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
  const auto& metrics = solver.GetGeometryMetrics();
  REQUIRE(metrics.GetBoundaryFaces().size() == 8);
  std::vector<Direction> closure(xy.size());
  double length = 0, xMoment = 0, yMoment = 0;
  for (unsigned long e = 0; e < metrics.GetnEdge(); ++e) {
    const auto& f = metrics.GetFace(e);
    for (unsigned d = 0; d < 2; ++d) {
      closure[f.first][d] += f.normal[d];
      closure[f.second][d] -= f.normal[d];
    }
  }
  for (const auto& f : metrics.GetBoundaryFaces()) {
    REQUIRE(f.node < 4);
    CHECK(f.marker == 0);
    length += std::hypot(f.normal[0], f.normal[1]);
    xMoment += f.point[0] * f.normal[0];
    yMoment += f.point[1] * f.normal[1];
    for (unsigned d = 0; d < 2; ++d) closure[f.node][d] += f.normal[d];
  }
  CHECK(length == Approx(6));
  CHECK(xMoment == Approx(2));
  CHECK(yMoment == Approx(2));
  for (const auto& sum : closure)
    for (auto value : sum) CHECK(value == Approx(0).margin(1e-14));
}

TEST_CASE("Constant primitive fields have zero full spatial residual", "[primitive-solver][primitive-residual]") {
  for (const std::string method : {"GREEN_GAUSS", "LEAST_SQUARES", "WEIGHTED_LEAST_SQUARES"}) {
    INFO(method);
    MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", "TIME_DISCRE_FLOW=EULER_EXPLICIT\n"
        "CONV_NUM_METHOD_FLOW=ROE\nNUM_METHOD_GRAD_RECON=" + method + "\nSLOPE_LIMITER_FLOW=NONE\n");
    CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
    const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate);
    for (unsigned long i = 0; i < xy.size(); ++i) {
      Check(result.faces[i], {});
      Check(result.volume[i], {});
      Check(result.boundary[i], {});
      Check(result.rhs[i], {});
      for (const auto& row : solver.GetLimitedGradient(i))
        for (auto value : row) CHECK(value == Approx(0).margin(1e-13));
    }
  }
}

TEST_CASE("Affine NODE reconstruction has zero jumps and the specified volume coefficients", "[primitive-solver][primitive-residual]") {
  for (const std::string mode : {"NODE_POINT", "MEAN_COEFF"}) {
    for (const std::string method : {"LEAST_SQUARES", "WEIGHTED_LEAST_SQUARES"}) {
      INFO(mode << " " << method);
      MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", "TIME_DISCRE_FLOW=EULER_EXPLICIT\n"
          "CONV_NUM_METHOD_FLOW=ROE\nNUM_METHOD_GRAD=GREEN_GAUSS\nNUM_METHOD_GRAD_RECON=" + method +
          "\nSLOPE_LIMITER_FLOW=NONE\nFLUID_MODEL=IDEAL_GAS\nGAMMA_VALUE=1.6\nPRIMITIVE_VOLUME_QUADRATURE=" + mode + "\n");
      CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
      SetAffine(solver);
      const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, exactAffine);
      const auto& metrics = solver.GetGeometryMetrics();
      for (unsigned long e = 0; e < metrics.GetnEdge(); ++e) {
        const auto& f = metrics.GetFace(e);
        Check(solver.ReconstructAt(f.first, f.point), Affine(f.point));
        Check(solver.ReconstructAt(f.second, f.point), Affine(f.point));
      }
      for (unsigned long i = 0; i < xy.size(); ++i) {
        Check(solver.GetPrimitiveNodes().GetPrimitiveState(i), Affine(xy[i]));
        CHECK(solver.GetPositivityScale(i) == 1);
        for (unsigned k = 0; k < 4; ++k)
          for (unsigned d = 0; d < 2; ++d)
            CHECK(solver.GetLimitedGradient(i)[k][d] == Approx(affineGradient[k][d]).margin(1e-13));
        const auto& cell = metrics.GetCell(i);
        const Direction reference = mode == "NODE_POINT" ? xy[i] :
            Direction{{xy[i][0] + cell.offset[0], xy[i][1] + cell.offset[1]}};
        const auto product = ExpandedVolume(Affine(reference), affineGradient, 1.6);
        State integrated{}, rhs{};
        for (unsigned k = 0; k < 4; ++k) {
          integrated[k] = cell.volume * product[k];
          rhs[k] = -product[k];
          CHECK(solver.LinSysRes(i, k) == Approx(integrated[k]).margin(2e-12));
        }
        Check(result.faces[i], {});
        Check(result.boundary[i], {});
        Check(result.volume[i], integrated);
        Check(result.rhs[i], rhs);
      }
    }
  }
}

TEST_CASE("Primitive volume product follows all four expanded equations", "[primitive-solver][primitive-residual]") {
  const State mean{{2.7, -0.4, 1.1, 3.2}};
  const Gradient gradient{{{{0.7, -0.3}}, {{-0.8, 0.4}}, {{0.5, -0.6}}, {{1.3, -0.9}}}};
  for (double gamma : {1.25, 1.6}) Check(VolumeProduct(mean, gradient, gamma), ExpandedVolume(mean, gradient, gamma));
  CHECK_THROWS_AS(VolumeProduct({{0, 1, 1, 1}}, gradient, 1.4), std::invalid_argument);
  auto invalid = gradient;
  invalid[2][0] = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_AS(VolumeProduct(mean, invalid, 1.4), std::invalid_argument);
}

TEST_CASE("First order assembles both fluctuations with positive signs exactly once", "[primitive-solver][primitive-residual]") {
  MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup + "MUSCL_FLOW=NO\n");
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
  SetAffine(solver);
  const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate);
  std::vector<State> expected(xy.size());
  State totalPath{}, totalFaces{};
  const auto& metrics = solver.GetGeometryMetrics();
  for (unsigned long e = 0; e < metrics.GetnEdge(); ++e) {
    const auto& f = metrics.GetFace(e);
    const auto left = Affine(xy[f.first]), right = Affine(xy[f.second]);
    Check(solver.ReconstructAt(f.first, f.point), left);
    const auto fluctuation = ComputeFluctuations(left, right, f.normal, 1.4);
    State jump{};
    for (unsigned k = 0; k < 4; ++k) {
      expected[f.first][k] += fluctuation.minus[k];
      expected[f.second][k] += fluctuation.plus[k];
      jump[k] = right[k] - left[k];
    }
    const auto path = Apply(PathMatrix(left, right, f.normal, 1.4), jump);
    for (unsigned k = 0; k < 4; ++k) totalPath[k] += path[k];
  }
  for (unsigned long i = 0; i < xy.size(); ++i) {
    Check(result.faces[i], expected[i]);
    Check(result.volume[i], {});
    Check(result.boundary[i], {});
    for (unsigned k = 0; k < 4; ++k) totalFaces[k] += result.faces[i][k];
  }
  Check(totalFaces, totalPath);
  CHECK(std::abs(totalFaces[0]) > 0.01);  // A conservative +/- scatter would sum to zero.
  const auto repeat = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate);
  for (unsigned long i = 0; i < xy.size(); ++i) Check(repeat.rhs[i], result.rhs[i]);
}

TEST_CASE("Limiter and positivity use one common reconstruction for faces and volume", "[primitive-solver][primitive-residual]") {
  for (const std::string limiter : {"NONE", "BARTH_JESPERSEN", "VENKATAKRISHNAN"}) {
    INFO(limiter);
    MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", "TIME_DISCRE_FLOW=EULER_EXPLICIT\n"
        "CONV_NUM_METHOD_FLOW=ROE\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\nSLOPE_LIMITER_FLOW=" + limiter + "\n");
    CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
    for (unsigned long i = 0; i < xy.size(); ++i)
      solver.GetPrimitiveNodes().SetPrimitiveState(i, {{i == 1 ? 100.0 : 1.0, 0.5, 0.3, i == 3 ? 100.0 : 1.0}});
    const auto before = solver.GetPrimitiveNodes().GetSolution();
    const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate);
    const auto& metrics = solver.GetGeometryMetrics();
    bool contracted = false, limited = false;
    for (unsigned long i = 0; i < xy.size(); ++i) {
      const auto& raw = solver.GetPrimitiveNodes().GetGradient_Reconstruction();
      const auto& gradient = solver.GetLimitedGradient(i);
      const double theta = solver.GetPositivityScale(i);
      CHECK(theta > 0);
      CHECK(theta <= 1);
      contracted |= theta < 1;
      for (unsigned k = 0; k < 4; ++k) {
        CHECK(solver.GetPrimitiveNodes().GetSolution(i, k) == before(i, k));
        const double phi = SU2_TYPE::GetValue(solver.GetPrimitiveNodes().GetLimiter(i, k));
        CHECK(std::isfinite(phi));
        CHECK(phi >= 0);
        limited |= phi < 0.99;
        if (limiter == "NONE") CHECK(phi == 1);
        for (unsigned d = 0; d < 2; ++d)
          CHECK(gradient[k][d] == Approx(SU2_TYPE::GetValue(raw(i, k, d)) * phi * theta).margin(1e-13));
      }
      const auto& cell = metrics.GetCell(i);
      const auto mean = Reconstruct(solver.GetPrimitiveNodes().GetPrimitiveState(i), gradient, cell.offset);
      REQUIRE(mean[0] > 0);
      REQUIRE(mean[3] > 0);
      auto expected = ExpandedVolume(mean, gradient, 1.4);
      for (auto& value : expected) value *= cell.volume;
      Check(result.volume[i], expected, 1e-10);
    }
    if (limiter == "NONE") CHECK(contracted);
    else CHECK(limited);
    for (unsigned long e = 0; e < metrics.GetnEdge(); ++e) {
      const auto& f = metrics.GetFace(e);
      for (auto i : {f.first, f.second}) {
        const auto w = solver.ReconstructAt(i, f.point);
        CHECK(w[0] > 0);
        CHECK(w[3] > 0);
      }
    }
    for (const auto& f : metrics.GetBoundaryFaces()) CHECK_NOTHROW(solver.ReconstructAt(f.node, f.point));
  }
}

TEST_CASE("Boundary data use half-edge quadrature and only incoming primitive fluctuations", "[primitive-solver][primitive-residual]") {
  MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup + "MUSCL_FLOW=NO\n");
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
  const State interior{{1, 0.5, 0.3, 1}};
  unsigned calls = 0;
  const auto exterior = [](const Direction& x, double t) { return State{{1.2 + 0.1*x[0], -0.2, 0.1, 1.1 + 0.1*t}}; };
  const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0.7,
      [&](const Direction& x, double t, const State& w, const Direction& n, unsigned short marker) {
        ++calls;
        Check(w, interior);
        CHECK(t == 0.7);
        CHECK(marker == 0);
        CHECK(std::hypot(n[0], n[1]) == Approx(1));
        return exterior(x, t);
      });
  const auto& boundary = solver.GetGeometryMetrics().GetBoundaryFaces();
  CHECK(calls == boundary.size());
  std::vector<State> expected(xy.size());
  for (const auto& f : boundary) {
    const auto waves = ComputeFluctuations(interior, exterior(f.point, 0.7), f.normal, 1.4);
    for (unsigned k = 0; k < 4; ++k) expected[f.node][k] += waves.minus[k];
  }
  for (unsigned long i = 0; i < xy.size(); ++i) Check(result.boundary[i], expected[i]);
}

TEST_CASE("NODE sources have nodal coordinates time sign and one volume factor", "[primitive-solver][primitive-residual]") {
  MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup);
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
  const auto source = [](const Direction& x, double t) { return State{{x[0]*x[0] + t, x[1] - t, 2*t, -0.3}}; };
  unsigned calls = 0;
  const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0.4, extrapolate,
      [&](const Direction& x, double t) {
        ++calls;
        CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), t, extrapolate),
                        std::logic_error);  // Reentry fails before touching shared work arrays.
        return source(x, t);
      });
  CHECK(calls == xy.size());
  for (unsigned long i = 0; i < xy.size(); ++i) {
    Check(result.source[i], source(xy[i], 0.4));
    Check(result.rhs[i], result.source[i]);
    for (unsigned k = 0; k < 4; ++k)
      CHECK(solver.LinSysRes(i, k) == Approx(-solver.GetGeometryMetrics().GetCell(i).volume * result.source[i][k]));
  }
  SetAffine(solver);
  const auto balanced = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, exactAffine,
      [&](const Direction& x, double) {
        // For NODE_POINT this cancels the affine volume exactly; MEAN_COEFF
        // intentionally uses different coefficients, tested separately above.
        return ExpandedVolume(Affine(x), affineGradient, 1.4);
      });
  CHECK(std::abs(balanced.rhs[0][0]) > 1e-5);  // Source must NOT be silently evaluated at the centroid.
  MeshFixture nodal(xy, triangles, "PRIMITIVE_EULER", setup + "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\n");
  CPrimitiveEulerSolver nodalSolver(nodal.geometry.get(), nodal.config.get(), 0);
  SetAffine(nodalSolver);
  const auto nodalBalance = nodalSolver.ComputeSpatialResidual(nodal.geometry.get(), nodal.config.get(), 0, exactAffine,
      [](const Direction& x, double) { return ExpandedVolume(Affine(x), affineGradient, 1.4); });
  for (const auto& rhs : nodalBalance.rhs) Check(rhs, {});
}

TEST_CASE("Failed spatial evaluations preserve the previous residual and unlock the solver", "[primitive-solver][primitive-residual]") {
  MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup);
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
  SetAffine(solver);
  solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate);
  std::vector<State> before(xy.size());
  for (unsigned long i = 0; i < xy.size(); ++i)
    for (unsigned k = 0; k < 4; ++k) before[i][k] = SU2_TYPE::GetValue(solver.LinSysRes(i, k));
  const auto unchanged = [&]() {
    for (unsigned long i = 0; i < xy.size(); ++i)
      for (unsigned k = 0; k < 4; ++k) CHECK(solver.LinSysRes(i, k) == before[i][k]);
  };
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get()), std::invalid_argument);
  unchanged();
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(nullptr, fixture.config.get(), 0, extrapolate), std::invalid_argument);
  unchanged();
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(),
      std::numeric_limits<double>::infinity(), extrapolate), std::invalid_argument);
  unchanged();
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0,
      [](const Direction&, double, const State&, const Direction&, unsigned short) { return State{{1, 0, 0, -1}}; }),
      std::invalid_argument);
  unchanged();
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate,
      [](const Direction&, double) { return State{{0, 0, std::numeric_limits<double>::quiet_NaN(), 0}}; }),
      std::invalid_argument);
  unchanged();
  solver.GetPrimitiveNodes().GetSolution(0)[0] = -1;
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate),
                  std::invalid_argument);
  unchanged();
  solver.GetPrimitiveNodes().SetPrimitiveState(0, Affine(xy[0]));
  CHECK_NOTHROW(solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0, extrapolate));
  unchanged();
}

TEST_CASE("Upwind virtual entry computes spatial residual and requires a prepared Forward Euler step", "[primitive-solver][primitive-residual]") {
  MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup, "MARKER_FAR");
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
  CSolver* base = &solver;
  CHECK_NOTHROW(base->Upwind_Residual(fixture.geometry.get(), nullptr, nullptr, fixture.config.get(), 0));
  for (unsigned long i = 0; i < xy.size(); ++i)
    for (unsigned k = 0; k < 4; ++k) CHECK(solver.LinSysRes(i, k) == Approx(0).margin(1e-13));
  SetAffine(solver);
  const auto direct = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get());
  CHECK_NOTHROW(base->Upwind_Residual(fixture.geometry.get(), nullptr, nullptr, fixture.config.get(), 0));
  for (unsigned long i = 0; i < xy.size(); ++i)
    for (unsigned k = 0; k < 4; ++k)
      CHECK(solver.LinSysRes(i, k) == Approx(-direct.rhs[i][k] * solver.GetGeometryMetrics().GetCell(i).volume));
  CHECK_THROWS_AS(base->Upwind_Residual(fixture.geometry.get(), nullptr, nullptr, fixture.config.get(), 1), std::invalid_argument);
  CHECK_THROWS_AS(base->ExplicitEuler_Iteration(fixture.geometry.get(), nullptr, fixture.config.get()), std::invalid_argument);
}

TEST_CASE("Spatial residual is invariant under node and element renumbering", "[primitive-solver][primitive-residual]") {
  MeshFixture first(xy, triangles, "PRIMITIVE_EULER", setup + "MUSCL_FLOW=NO\n");
  CPrimitiveEulerSolver a(first.geometry.get(), first.config.get(), 0);
  SetAffine(a);
  const auto ra = a.ComputeSpatialResidual(first.geometry.get(), first.config.get(), 0, extrapolate);
  const std::array<unsigned long, 5> permutation{{3, 0, 4, 1, 2}};
  std::vector<Direction> permuted(5);
  for (unsigned i = 0; i < 5; ++i) permuted[permutation[i]] = xy[i];
  auto cells = triangles;
  for (auto& triangle : cells)
    for (auto& node : triangle) node = permutation[node];
  std::reverse(cells.begin(), cells.end());
  MeshFixture second(permuted, cells, "PRIMITIVE_EULER", setup + "MUSCL_FLOW=NO\n");
  CPrimitiveEulerSolver b(second.geometry.get(), second.config.get(), 0);
  for (unsigned i = 0; i < 5; ++i) b.GetPrimitiveNodes().SetPrimitiveState(permutation[i], Affine(xy[i]));
  const auto rb = b.ComputeSpatialResidual(second.geometry.get(), second.config.get(), 0, extrapolate);
  for (unsigned i = 0; i < 5; ++i) Check(rb.rhs[permutation[i]], ra.rhs[i]);
}

TEST_CASE("Unaudited reconstruction variants and boundary kinds are rejected explicitly", "[primitive-solver][primitive-residual]") {
  for (const std::string option : {"SLOPE_LIMITER_FLOW=VAN_ALBADA_EDGE\n", "MUSCL_KAPPA_FLOW=0.3\n"}) {
    INFO(option);
    MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", "TIME_DISCRE_FLOW=EULER_EXPLICIT\n"
        "CONV_NUM_METHOD_FLOW=ROE\n" + option);
    CHECK_THROWS_AS(CPrimitiveEulerSolver(fixture.geometry.get(), fixture.config.get(), 0), std::invalid_argument);
  }
  MeshFixture symmetry(xy, triangles, "PRIMITIVE_EULER", setup, "MARKER_SYM");
  CPrimitiveEulerSolver solver(symmetry.geometry.get(), symmetry.config.get(), 0);
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(symmetry.geometry.get(), symmetry.config.get(), 0, extrapolate),
                  std::invalid_argument);
}

TEST_CASE("Effective gamma agrees across primitive diagnostics and the full spatial operator", "[primitive-solver][primitive-hardening]") {
  for (const auto& gas : std::vector<std::pair<std::string, double>>{
           {"STANDARD_AIR\nGAMMA_VALUE=1.6", 1.4}, {"IDEAL_GAS\nGAMMA_VALUE=1.6", 1.6},
           {"IDEAL_GAS\nGAMMA_VALUE=1.25", 1.25}}) {
    INFO(gas.first);
    MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup + "FLUID_MODEL=" + gas.first + "\n");
    CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
    const double gamma = gas.second;
    CHECK(fixture.config->GetGamma() == Approx(gamma).epsilon(1e-14));
    for (unsigned i = 0; i < xy.size(); ++i) {
      auto w = Affine(xy[i]);
      w[0] += 0.1 * xy[i][0] * xy[i][1];
      w[1] += 0.02 * xy[i][1] * xy[i][1];
      w[3] += 0.15 * xy[i][0] * xy[i][0];
      solver.GetPrimitiveNodes().SetPrimitiveState(i, w);
      const auto& nodes = solver.GetPrimitiveNodes();
      CHECK(nodes.GetSoundSpeed(i) == Approx(std::sqrt(gamma * w[3] / w[0])).epsilon(1e-13));
      CHECK(nodes.GetEnergy(i) == Approx(w[3]/((gamma-1)*w[0]) + 0.5*(w[1]*w[1]+w[2]*w[2])).epsilon(1e-13));
      CHECK(nodes.GetEnthalpy(i) == Approx(nodes.GetEnergy(i)+w[3]/w[0]).epsilon(1e-13));
      CHECK(nodes.GetConservativeState(i)[3] == Approx(w[3]/(gamma-1) + 0.5*w[0]*(w[1]*w[1]+w[2]*w[2])).epsilon(1e-13));
      CHECK(PhysicalMatrix(w, {{1,0}}, gamma)[3][1] == Approx(gamma*w[3]).epsilon(1e-13));
    }
    const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(), fixture.config.get(), 0.2, exactAffine);
    std::vector<State> expectedFaces(xy.size()), expectedBoundary(xy.size());
    double gammaSignal = 0;
    const auto& metrics = solver.GetGeometryMetrics();
    for (unsigned long e = 0; e < metrics.GetnEdge(); ++e) {
      const auto& f = metrics.GetFace(e);
      const auto l = solver.ReconstructAt(f.first, f.point), r = solver.ReconstructAt(f.second, f.point);
      const auto waves = ComputeFluctuations(l, r, f.normal, gamma);
      const auto wrong = ComputeFluctuations(l, r, f.normal, gamma == 1.4 ? 1.6 : 1.4);
      gammaSignal += std::abs(waves.minus[3]-wrong.minus[3]);
      for (unsigned k = 0; k < 4; ++k) {
        expectedFaces[f.first][k] += waves.minus[k];
        expectedFaces[f.second][k] += waves.plus[k];
      }
    }
    for (const auto& f : metrics.GetBoundaryFaces()) {
      const auto waves = ComputeFluctuations(solver.ReconstructAt(f.node, f.point), Affine(f.point), f.normal, gamma);
      for (unsigned k = 0; k < 4; ++k) expectedBoundary[f.node][k] += waves.minus[k];
    }
    CHECK(gammaSignal > 1e-5); // The test distinguishes the requested and effective gamma.
    for (unsigned i = 0; i < xy.size(); ++i) {
      Check(result.faces[i], expectedFaces[i]);
      Check(result.boundary[i], expectedBoundary[i]);
      const auto& cell = metrics.GetCell(i);
      const auto& gradient = solver.GetLimitedGradient(i);
      const auto mean = Reconstruct(solver.GetPrimitiveNodes().GetPrimitiveState(i), gradient, cell.offset);
      auto expected = ExpandedVolume(mean, gradient, gamma);
      Check(VolumeProduct(mean, gradient, gamma), expected);
      for (auto& value : expected) value *= cell.volume;
      Check(result.volume[i], expected);
    }
  }
}

TEST_CASE("Multiple physical markers preserve half-face ownership normals and distinct data", "[primitive-solver][primitive-hardening]") {
  const std::vector<BoundaryMarker> markers{{"south_east", {{{0,1}}, {{1,2}}}},
                                           {"north_west", {{{2,3}}, {{3,0}}}}};
  struct Expected { unsigned long node; std::string marker; Direction q, normal; };
  const std::vector<Expected> expected{
      {0,"south_east",{{0.5,0}},{{0,-1}}}, {1,"south_east",{{1.5,0}},{{0,-1}}},
      {1,"south_east",{{2,0.25}},{{0.5,0}}}, {2,"south_east",{{2,0.75}},{{0.5,0}}},
      {2,"north_west",{{1.5,1}},{{0,1}}}, {3,"north_west",{{0.5,1}},{{0,1}}},
      {3,"north_west",{{0,0.75}},{{-0.5,0}}}, {0,"north_west",{{0,0.25}},{{-0.5,0}}}};
  for (bool reverse : {false, true}) {
    INFO(reverse);
    MeshFixture fixture(xy, triangles, "PRIMITIVE_EULER", setup + "MUSCL_FLOW=NO\n", "MARKER_CUSTOM", markers, reverse);
    if (reverse) {
      // Also reverse AFTER SU2's orientation repair: retain native metrics but
      // exercise our constructor on genuinely reversed boundary connectivity.
      for (unsigned m = 0; m < fixture.geometry->GetnMarker(); ++m)
        for (unsigned long b = 0; b < fixture.geometry->GetnElem_Bound(m); ++b)
          fixture.geometry->bound[m][b]->Change_Orientation();
    }
    CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), 0);
    const auto& metrics = solver.GetGeometryMetrics();
    REQUIRE(metrics.GetBoundaryFaces().size() == expected.size());
    std::vector<unsigned> coverage(expected.size(), 0);
    std::vector<Direction> closure(xy.size());
    std::map<std::pair<unsigned short,unsigned long>, Direction> nativeSums;
    for (const auto& f : metrics.GetBoundaryFaces()) {
      const auto tag = fixture.config->GetMarker_All_TagBound(f.marker);
      unsigned matched = 0;
      for (unsigned n = 0; n < expected.size(); ++n) {
        const auto& ref = expected[n];
        if (ref.node != f.node || ref.marker != tag || ref.q != f.point) continue;
        ++coverage[n]; ++matched;
        for (unsigned d = 0; d < 2; ++d) CHECK(f.normal[d] == Approx(ref.normal[d]).margin(1e-14));
        CHECK(std::hypot(f.normal[0],f.normal[1]) == Approx(std::hypot(ref.normal[0],ref.normal[1])));
      }
      CHECK(matched == 1);
      for (unsigned d = 0; d < 2; ++d) {
        closure[f.node][d] += f.normal[d];
        nativeSums[{f.marker,f.node}][d] += f.normal[d];
      }
    }
    for (auto covered : coverage) CHECK(covered == 1);
    for (unsigned m = 0; m < fixture.geometry->GetnMarker(); ++m)
      for (unsigned long v = 0; v < fixture.geometry->GetnVertex(m); ++v) {
        const auto* vertex = fixture.geometry->vertex[m][v];
        for (unsigned d = 0; d < 2; ++d)
          CHECK(nativeSums[std::make_pair(m,vertex->GetNode())][d] == Approx(-vertex->GetNormal(d)).margin(1e-14));
      }
    for (unsigned long e = 0; e < metrics.GetnEdge(); ++e) {
      const auto& f = metrics.GetFace(e);
      for (unsigned d = 0; d < 2; ++d) { closure[f.first][d] += f.normal[d]; closure[f.second][d] -= f.normal[d]; }
    }
    for (const auto& sum : closure)
      for (auto value : sum) CHECK(value == Approx(0).margin(1e-14));
    const auto exterior = [](const std::string& tag) { return tag == "south_east" ? State{{1.3,-0.2,0.1,1.4}} : State{{0.8,0.4,-0.3,0.7}}; };
    std::vector<State> boundaryExpected(xy.size());
    const State inside{{1,0.5,0.3,1}};
    for (const auto& ref : expected) {
      const auto wave = ComputeFluctuations(inside,exterior(ref.marker),ref.normal,1.4).minus;
      for (unsigned k = 0; k < 4; ++k) boundaryExpected[ref.node][k] += wave[k];
    }
    std::vector<unsigned> callbacks(expected.size(),0);
    const auto result = solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
        [&](const Direction& q,double,const State&,const Direction&,unsigned short marker) {
          const auto tag = fixture.config->GetMarker_All_TagBound(marker);
          unsigned matched = 0;
          for (unsigned n = 0; n < expected.size(); ++n)
            if (q == expected[n].q && tag == expected[n].marker) { ++callbacks[n]; ++matched; }
          CHECK(matched == 1);
          return exterior(tag);
        });
    for (auto calls : callbacks) CHECK(calls == 1);
    for (unsigned i = 0; i < xy.size(); ++i) Check(result.boundary[i],boundaryExpected[i]);
  }
}

TEST_CASE("MUSCL spatial fields are invariant under nontrivial renumbering", "[primitive-solver][primitive-hardening]") {
  const std::array<unsigned long,5> p{{3,0,4,1,2}};
  std::vector<Direction> permuted(5);
  for (unsigned i = 0; i < 5; ++i) permuted[p[i]] = xy[i];
  auto cells = triangles;
  for (auto& t : cells) for (auto& node : t) node = p[node];
  std::reverse(cells.begin(),cells.end());
  const auto close = [](double a,double b) { CHECK(a == Approx(b).epsilon(2e-11).margin(5e-13)); };
  for (const std::string method : {"LEAST_SQUARES","WEIGHTED_LEAST_SQUARES"}) {
    for (bool stress : {false,true}) {
      INFO(method << " stress=" << stress);
      const std::string options = "TIME_DISCRE_FLOW=EULER_EXPLICIT\nCONV_NUM_METHOD_FLOW=ROE\nMUSCL_FLOW=YES\n"
          "SLOPE_LIMITER_FLOW=NONE\nNUM_METHOD_GRAD_RECON=" + method + "\n";
      MeshFixture first(xy,triangles,"PRIMITIVE_EULER",options), second(permuted,cells,"PRIMITIVE_EULER",options);
      CPrimitiveEulerSolver a(first.geometry.get(),first.config.get(),0), b(second.geometry.get(),second.config.get(),0);
      for (unsigned i = 0; i < 5; ++i) {
        auto w = Affine(xy[i]);
        w[0] += 0.2*xy[i][0]*xy[i][1]; w[3] += 0.1*xy[i][0]*xy[i][0];
        if (stress) { w[0] = i == 1 ? 100 : 1; w[3] = i == 3 ? 100 : 1; }
        a.GetPrimitiveNodes().SetPrimitiveState(i,w); b.GetPrimitiveNodes().SetPrimitiveState(p[i],w);
      }
      const auto ra = a.ComputeSpatialResidual(first.geometry.get(),first.config.get(),0,exactAffine);
      const auto rb = b.ComputeSpatialResidual(second.geometry.get(),second.config.get(),0,exactAffine);
      bool contracted = false;
      for (unsigned i = 0; i < 5; ++i) {
        close(a.GetPositivityScale(i),b.GetPositivityScale(p[i]));
        contracted |= a.GetPositivityScale(i) < 1;
        for (unsigned k = 0; k < 4; ++k) {
          for (unsigned d = 0; d < 2; ++d) close(a.GetLimitedGradient(i)[k][d],b.GetLimitedGradient(p[i])[k][d]);
          for (auto field : {&CPrimitiveEulerSolver::SpatialResult::faces, &CPrimitiveEulerSolver::SpatialResult::volume,
                             &CPrimitiveEulerSolver::SpatialResult::boundary, &CPrimitiveEulerSolver::SpatialResult::rhs})
            close((ra.*field)[i][k],(rb.*field)[p[i]][k]);
          close(SU2_TYPE::GetValue(a.LinSysRes(i,k)),SU2_TYPE::GetValue(b.LinSysRes(p[i],k)));
        }
      }
      if (stress) CHECK(contracted);
      for (unsigned long e = 0; e < a.GetGeometryMetrics().GetnEdge(); ++e) {
        const auto& f = a.GetGeometryMetrics().GetFace(e);
        for (auto i : {f.first,f.second}) {
          const auto wa = a.ReconstructAt(i,f.point), wb = b.ReconstructAt(p[i],f.point);
          for (unsigned k = 0; k < 4; ++k) close(wa[k],wb[k]);
        }
      }
      for (const auto& f : a.GetGeometryMetrics().GetBoundaryFaces()) {
        const auto wa = a.ReconstructAt(f.node,f.point), wb = b.ReconstructAt(p[f.node],f.point);
        for (unsigned k = 0; k < 4; ++k) close(wa[k],wb[k]);
      }
    }
  }
}

TEST_CASE("Positivity contraction remains finite for overflowing differences and tiny nodal values", "[primitive-solver][primitive-hardening]") {
  const double large = 1e308;
  CHECK(std::isinf(large-(-large))); // Old denominator overflowed, spuriously producing theta=0.
  CHECK(ReconstructionPositivityFactor(large,-large) == Approx(0.495).epsilon(1e-14));
  CHECK(ReconstructionPositivityFactor(1e-300,-1e-298) == Approx(0.99*0.9/101).epsilon(1e-13));
  CHECK(ReconstructionPositivityFactor(1,1) == 1);
  CHECK_THROWS_AS(ReconstructionPositivityFactor(0,1),std::invalid_argument);
  const double denormal = std::numeric_limits<double>::denorm_min();
  if (denormal > 0) CHECK(ReconstructionPositivityFactor(denormal,0) == 0);
  for (double scale : {1e-12,1e-100,1e-280}) {
    INFO(scale);
    MeshFixture fixture(xy,triangles,"PRIMITIVE_EULER",setup);
    CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
    std::vector<State> before(5);
    for (unsigned i = 0; i < 5; ++i) {
      before[i] = {{scale*(i==1?100:1),0.5,0.3,scale*(i==3?100:1)}};
      solver.GetPrimitiveNodes().SetPrimitiveState(i,before[i]);
    }
    solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,extrapolate);
    bool contracted = false;
    const auto valid = [](const State& w) {
      CHECK(w[0] > 0); CHECK(w[3] > 0);
      for (auto value : w) CHECK(std::isfinite(value));
    };
    for (unsigned i = 0; i < 5; ++i) {
      const double theta = solver.GetPositivityScale(i);
      CHECK(std::isfinite(theta)); CHECK(theta >= 0); CHECK(theta <= 1);
      contracted |= theta < 1;
      CHECK(solver.GetPrimitiveNodes().GetPrimitiveState(i) == before[i]);
      valid(Reconstruct(before[i],solver.GetLimitedGradient(i),solver.GetGeometryMetrics().GetCell(i).offset));
    }
    CHECK(contracted);
    for (unsigned long e = 0; e < solver.GetGeometryMetrics().GetnEdge(); ++e) {
      const auto& f = solver.GetGeometryMetrics().GetFace(e);
      valid(solver.ReconstructAt(f.first,f.point)); valid(solver.ReconstructAt(f.second,f.point));
    }
    for (const auto& f : solver.GetGeometryMetrics().GetBoundaryFaces()) valid(solver.ReconstructAt(f.node,f.point));
  }
}

TEST_CASE("Callback failures of every component preserve residual ownership and release evaluating", "[primitive-solver][primitive-hardening]") {
  MeshFixture fixture(xy,triangles,"PRIMITIVE_EULER",setup);
  CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
  SetAffine(solver);
  const auto baseline = solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactAffine);
  const auto retry = [&]() {
    for (unsigned i = 0; i < 5; ++i) for (unsigned k = 0; k < 4; ++k)
      CHECK(solver.LinSysRes(i,k) == -solver.GetGeometryMetrics().GetCell(i).volume*baseline.rhs[i][k]);
    CHECK_NOTHROW(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactAffine));
  };
  for (double bad : {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),
                     -std::numeric_limits<double>::infinity()}) {
    for (unsigned component = 0; component < 4; ++component) {
      INFO(component << " " << bad);
      CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
          [&](const Direction&,double,const State& inside,const Direction&,unsigned short) {
            State w=inside; w[component]=bad; return w;
          }),std::invalid_argument);
      retry();
      CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactAffine,
          [&](const Direction&,double) { State s{}; s[component]=bad; return s; }),std::invalid_argument);
      retry();
      solver.GetPrimitiveNodes().GetSolution(0)[component]=bad;
      CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactAffine),std::invalid_argument);
      solver.GetPrimitiveNodes().SetPrimitiveState(0,Affine(xy[0]));
      retry();
    }
  }
  for (unsigned component : {0u,3u}) for (double bad : {0.0,-1.0}) {
    CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
        [&](const Direction&,double,const State& inside,const Direction&,unsigned short) {
          auto w=inside; w[component]=bad; return w;
        }),std::invalid_argument);
    retry();
    solver.GetPrimitiveNodes().GetSolution(0)[component]=bad;
    CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactAffine),std::invalid_argument);
    solver.GetPrimitiveNodes().SetPrimitiveState(0,Affine(xy[0]));
    retry();
  }
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
      [](const Direction&,double,const State&,const Direction&,unsigned short)->State { throw std::runtime_error("boundary failure"); }),std::runtime_error);
  retry();
  CHECK_THROWS_AS(solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactAffine,
      [](const Direction&,double)->State { throw std::runtime_error("source failure"); }),std::runtime_error);
  retry();
}

TEST_CASE("ReconstructAt has documented zero initial gradient and requires refresh after nodal writes", "[primitive-solver][primitive-hardening]") {
  MeshFixture fixture(xy,triangles,"PRIMITIVE_EULER",setup);
  CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
  const Direction q{{0.5,0.2}};
  Check(solver.ReconstructAt(0,q),solver.GetPrimitiveNodes().GetPrimitiveState(0));
  CHECK_THROWS_AS(solver.ReconstructAt(5,q),std::out_of_range);
  SetAffine(solver);
  solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,exactAffine);
  Check(solver.ReconstructAt(0,q),Affine(q));
  auto w=Affine(xy[0]); w[1]+=0.1;
  solver.GetPrimitiveNodes().SetPrimitiveState(0,w);
  auto stale=Affine(q); stale[1]+=0.1;
  Check(solver.ReconstructAt(0,q),stale); // Existing API combines current W with last prepared G.
  solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,extrapolate);
  CHECK(std::abs(solver.ReconstructAt(0,q)[1]-stale[1]) > 1e-4);
}

namespace {
State SmoothState(const Direction& x) {
  return {{2+0.15*std::sin(x[0])*std::cos(x[1]),
           0.5+0.1*std::cos(1.3*x[0])*std::sin(0.7*x[1]),
           0.3+0.08*std::sin(0.8*x[0])*std::cos(1.1*x[1]),
           3+0.2*std::cos(0.9*x[0])*std::cos(1.2*x[1])}};
}
State SmoothSource(const Direction& x) {
  const auto w=SmoothState(x);
  const double rx=0.15*std::cos(x[0])*std::cos(x[1]), ry=-0.15*std::sin(x[0])*std::sin(x[1]);
  const double ux=-0.13*std::sin(1.3*x[0])*std::sin(0.7*x[1]), uy=0.07*std::cos(1.3*x[0])*std::cos(0.7*x[1]);
  const double vx=0.064*std::cos(0.8*x[0])*std::cos(1.1*x[1]), vy=-0.088*std::sin(0.8*x[0])*std::sin(1.1*x[1]);
  const double px=-0.18*std::sin(0.9*x[0])*std::cos(1.2*x[1]), py=-0.24*std::cos(0.9*x[0])*std::sin(1.2*x[1]);
  return {{w[1]*rx+w[2]*ry+w[0]*(ux+vy), w[1]*ux+w[2]*uy+px/w[0],
           w[1]*vx+w[2]*vy+py/w[0], w[1]*px+w[2]*py+1.4*w[3]*(ux+vy)}};
}
}

TEST_CASE("Spatial MMS measures truncation without imposing an expected order", "[primitive-solver][primitive-spatial-mms]") {
  // Steady analytic W with S=Ax(W)Wx+Ay(W)Wy, so the exact rhs is zero.
  // No time stepping, no numerical derivatives in the reference source.
  std::cout << "SPATIAL_MMS,family,variant,N,L2_all,L2_interior,L2_boundary,p_all,p_interior,p_boundary\n";
  for (bool distorted : {false,true}) {
    for (const std::string variant : {"NO","NODE_POINT","MEAN_COEFF"}) {
      std::array<double,3> previous{};
      for (unsigned n : {8u,16u,32u,64u}) {
        std::vector<Direction> grid;
        std::vector<Triangle> cells;
        const double h=1.0/n;
        const auto index=[&](unsigned i,unsigned j) { return static_cast<unsigned long>(j*(n+1)+i); };
        for (unsigned j=0;j<=n;++j) for (unsigned i=0;i<=n;++i) {
          Direction x{{i*h,j*h}};
          if (distorted && i>0 && i<n && j>0 && j<n) {
            x[0]+=0.18*h*std::sin(2.1*i+1.3*j);
            x[1]+=0.16*h*std::cos(1.7*i-0.9*j);
          }
          grid.push_back(x);
        }
        for (unsigned j=0;j<n;++j) for (unsigned i=0;i<n;++i) {
          cells.push_back({{index(i,j),index(i+1,j),index(i+1,j+1)}});
          cells.push_back({{index(i,j),index(i+1,j+1),index(i,j+1)}});
        }
        const std::string options = setup + (variant=="NO" ? "MUSCL_FLOW=NO\n" :
            "MUSCL_FLOW=YES\nPRIMITIVE_VOLUME_QUADRATURE="+variant+"\n");
        MeshFixture fixture(grid,cells,"PRIMITIVE_EULER",options);
        CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
        for (unsigned long i=0;i<grid.size();++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,SmoothState(grid[i]));
        const auto result=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
            [](const Direction& x,double,const State&,const Direction&,unsigned short) { return SmoothState(x); },
            [](const Direction& x,double) { return SmoothSource(x); });
        std::array<double,3> squared{},weights{},error{},rate{};
        for (unsigned j=0;j<=n;++j) for (unsigned i=0;i<=n;++i) {
          const auto point=index(i,j);
          const unsigned group=(i==0||i==n||j==0||j==n)?2:1;
          const double volume=solver.GetGeometryMetrics().GetCell(point).volume;
          double norm=0;
          for (auto value:result.rhs[point]) norm+=value*value;
          for (unsigned g : {0u,group}) { squared[g]+=volume*norm; weights[g]+=volume; }
        }
        for (unsigned g=0;g<3;++g) {
          error[g]=std::sqrt(squared[g]/weights[g]);
          REQUIRE(std::isfinite(error[g])); REQUIRE(error[g]>0);
          rate[g]=n==8?0:std::log(previous[g]/error[g])/std::log(2.0);
        }
        std::cout << std::setprecision(10) << "SPATIAL_MMS," << (distorted?"distorted":"regular") << ',' << variant << ',' << n;
        for (auto value:error) std::cout << ',' << value;
        for (auto value:rate) std::cout << ',' << value;
        std::cout << '\n';
        previous=error;
      }
    }
  }
}

TEST_CASE("Boundary geometry rejects missing duplicate and out of range connectivity", "[primitive-solver][primitive-hardening]") {
  for (unsigned variant=0;variant<3;++variant) {
    INFO(variant);
    MeshFixture fixture(xy,triangles,"PRIMITIVE_EULER",setup);
    auto* line=fixture.geometry->bound[0][0];
    if (variant==0) { line->SetNode(0,0); line->SetNode(1,2); } // Nonexistent primal edge.
    if (variant==1) { line->SetNode(0,xy.size()); }
    if (variant==2) {
      const auto* other=fixture.geometry->bound[0][1];
      line->SetNode(0,other->GetNode(0)); line->SetNode(1,other->GetNode(1));
    }
    CHECK_THROWS_AS(CAggregatedGeometry(*fixture.geometry,*fixture.config),std::invalid_argument);
  }
}

TEST_CASE("First order has a reproducible affine contact consistency defect on an asymmetric stencil", "[primitive-solver][primitive-hardening]") {
  const auto contact=[](const Direction& x) { return State{{2+0.1*x[0]+0.2*x[1],0.5,0.3,3}}; };
  const State source{{0.11,0,0,0}}; // u*rho_x + v*rho_y = 0.5*0.1 + 0.3*0.2.
  double referenceDefect=0;
  for (double h : {1.0,0.5,0.25,0.125}) {
    std::vector<Direction> patch=xy;
    for (auto& x:patch) for (unsigned d=0;d<2;++d) x[d]=h*(x[d]-xy[4][d]);
    for (bool muscl : {false,true}) {
      MeshFixture fixture(patch,triangles,"PRIMITIVE_EULER",setup+(muscl?"MUSCL_FLOW=YES\n":"MUSCL_FLOW=NO\n"));
      CPrimitiveEulerSolver solver(fixture.geometry.get(),fixture.config.get(),0);
      for (unsigned i=0;i<5;++i) solver.GetPrimitiveNodes().SetPrimitiveState(i,contact(patch[i]));
      const auto result=solver.ComputeSpatialResidual(fixture.geometry.get(),fixture.config.get(),0,
          [&](const Direction& x,double,const State&,const Direction&,unsigned short) { return contact(x); },
          [&](const Direction&,double) { return source; });
      if (muscl) { Check(result.rhs[4],{},2e-13); continue; }
      double incoming=0;
      const auto& metrics=solver.GetGeometryMetrics();
      for (unsigned long e=0;e<metrics.GetnEdge();++e) {
        const auto& f=metrics.GetFace(e);
        if (f.first!=4 && f.second!=4) continue;
        const auto other=f.first==4?f.second:f.first;
        const double sign=f.first==4?1:-1;
        const double speed=sign*(0.5*f.normal[0]+0.3*f.normal[1]);
        incoming+=std::min(speed,0.0)*(contact(patch[other])[0]-contact(patch[4])[0]);
      }
      const double defect=0.11-incoming/metrics.GetCell(4).volume;
      CHECK(result.rhs[4][0] == Approx(defect).epsilon(1e-11).margin(1e-13));
      CHECK(std::abs(defect)>1e-3);
      if (h==1) referenceDefect=defect;
      else CHECK(defect == Approx(referenceDefect).epsilon(1e-11).margin(1e-13));
      Check(result.boundary[4],{}); // Strictly interior: no physical boundary contribution.
      Check(result.volume[4],{});
      std::cout << "CONTACT_DEFECT,h=" << h << ",rhs_rho=" << std::setprecision(12) << defect << '\n';
    }
  }
}
