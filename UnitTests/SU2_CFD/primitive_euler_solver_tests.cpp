/*!
 * \brief Phase-3 primitive state, configuration and factory integration tests.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "primitive_euler_mesh_fixture.hpp"
#include "../../SU2_CFD/include/solvers/CPrimitiveEulerSolver.hpp"
#include "../../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../../SU2_CFD/include/variables/CEulerVariable.hpp"
#include <functional>
#include <limits>

using namespace PrimitiveEuler2D;
using namespace PrimitiveEulerTest;

namespace {
const std::vector<Direction> points{{{0, 0}}, {{2, 0}}, {{0, 1}}};
const std::vector<Triangle> elements{{{0, 1, 2}}};
const std::string explicitEuler = "TIME_DISCRE_FLOW= EULER_EXPLICIT\n";

void CheckState(const State& actual, const State& expected) {
  for (unsigned k = 0; k < 4; ++k) CHECK(actual[k] == Approx(expected[k]).margin(1e-13));
}

struct FactoryContainer {
  CSolver** solvers;
  explicit FactoryContainer(MeshFixture& fixture)
      : solvers(CSolverFactory::CreateSolverContainer(fixture.config->GetKind_Solver(), fixture.config.get(),
                                                       fixture.geometry.get(), MESH_0)) {}
  ~FactoryContainer() {
    for (unsigned i = 0; i < MAX_SOLS; ++i) delete solvers[i];
    delete[] solvers;
  }
};
}  // namespace

TEST_CASE("Primitive solver config factory and nodal storage agree", "[primitive-solver]") {
  MeshFixture fixture(points, elements, "PRIMITIVE_EULER", explicitEuler);
  CHECK(fixture.config->GetFluidProblem());
  CHECK(fixture.config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  CHECK(fixture.config->GetPrimitiveRepresentation() == "NODE");
  CHECK(fixture.config->GetPrimitiveFaceQuadrature() == "AGGREGATED");
  CHECK(fixture.config->GetPrimitiveVolumeQuadrature() == "MEAN_COEFF");
  FactoryContainer container(fixture);
  auto* solver = dynamic_cast<CPrimitiveEulerSolver*>(container.solvers[FLOW_SOL]);
  REQUIRE(solver != nullptr);
  REQUIRE(solver->GetnVar() == 4);
  const auto metadata = CSolverFactory::GetSolverMeta(solver);
  CHECK(metadata.solverType == SUB_SOLVER_TYPE::PRIMITIVE_EULER);
  CHECK(metadata.integrationType == INTEGRATION_TYPE::NONE);
  for (unsigned i = 0; i < MAX_SOLS; ++i)
    if (i != FLOW_SOL) CHECK(container.solvers[i] == nullptr);
  CHECK(solver->GetGeometryMetrics().GetnPoint() == 3);
  auto& nodes = solver->GetPrimitiveNodes();
  CHECK(solver->GetNodes() == &nodes);
  for (unsigned long i = 0; i < 3; ++i) {
    CheckState(nodes.GetPrimitiveState(i), {{1, 0.5, 0.3, 1}});
    for (unsigned k = 0; k < 4; ++k) {
      CHECK(nodes.GetSolution(i, k) == nodes.GetPrimitive(i, k));
      CHECK(nodes.GetSolution_Old(i, k) == nodes.GetSolution(i, k));
      CHECK(nodes.GetLimiter(i, k) == 1.0);
      for (unsigned d = 0; d < 2; ++d) CHECK(nodes.GetGradient(i, k, d) == 0.0);
    }
    CHECK(nodes.GetDelta_Time(i) == 0.0);
    CHECK(solver->LinSysRes.GetBlock(i)[0] == 0.0);
  }
}

TEST_CASE("Primitive state thermodynamics and conservative conversion have no stale cache", "[primitive-solver]") {
  MeshFixture fixture(points, elements, "PRIMITIVE_EULER",
                      explicitEuler + "FLUID_MODEL=IDEAL_GAS\nGAMMA_VALUE=1.6\n"
                                      "PRIMITIVE_INITIAL_STATE=(2, 3, -4, 5)\n");
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), MESH_0);
  auto& nodes = solver.GetPrimitiveNodes();
  CheckState(nodes.GetPrimitiveState(0), {{2, 3, -4, 5}});
  const double energy = 5 / 0.6 + 25;
  CheckState(nodes.GetConservativeState(0), {{2, 6, -8, energy}});
  const su2double velocity[2] = {3, -4};
  CEulerVariable standardEuler(2, velocity, energy / 2, 1, 2, 4, fixture.config.get());
  for (unsigned k = 0; k < 4; ++k)
    CHECK(standardEuler.GetSolution(0, k) == Approx(nodes.GetConservativeState(0)[k]));
  CHECK(nodes.GetSolution(0, 1) != standardEuler.GetSolution(0, 1));
  CHECK(nodes.GetSolution(0, 3) != standardEuler.GetSolution(0, 3));
  CHECK(nodes.GetEnergy(0) == Approx(energy / 2));
  CHECK(nodes.GetEnthalpy(0) == Approx((energy + 5) / 2));
  CHECK(nodes.GetSoundSpeed(0) == Approx(2.0));
  const State next{{1.5, -0.2, 0.7, 2.1}};
  nodes.SetPrimitiveState(0, next);
  CheckState(nodes.GetPrimitiveState(0), next);
  CHECK(nodes.GetDensity(0) == Approx(1.5));
  CHECK(nodes.GetVelocity(0, 0) == Approx(-0.2));
  CHECK(nodes.GetVelocity(0, 1) == Approx(0.7));
  CHECK(nodes.GetPressure(0) == Approx(2.1));
  CHECK(nodes.GetSoundSpeed(0) == Approx(std::sqrt(1.6 * 2.1 / 1.5)));
  const auto conservative = nodes.GetConservativeState(0);
  const State back{{conservative[0], conservative[1] / conservative[0], conservative[2] / conservative[0],
                     0.6 * (conservative[3] - 0.5 * (conservative[1] * conservative[1] +
                                                    conservative[2] * conservative[2]) / conservative[0])}};
  CheckState(back, next);
  CheckState(nodes.GetPrimitiveState(1), {{2, 3, -4, 5}});
  CHECK_THROWS_AS(nodes.SetPrimitiveState(0, {{-1, 0, 0, 1}}), std::invalid_argument);
  CheckState(nodes.GetPrimitiveState(0), next);
  CHECK_THROWS_AS(nodes.GetPrimitiveState(3), std::out_of_range);
}

TEST_CASE("Primitive NODE_POINT and standard limiter options parse independently", "[primitive-solver]") {
  MeshFixture fixture(points, elements, "PRIMITIVE_EULER",
                      explicitEuler + "PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT\nSLOPE_LIMITER_FLOW=NONE\n");
  CHECK(fixture.config->GetPrimitiveVolumeQuadrature() == "NODE_POINT");
  CHECK(fixture.config->GetKind_SlopeLimit_Flow() == LIMITER::NONE);
  CHECK_NOTHROW(CPrimitiveEulerSolver(fixture.geometry.get(), fixture.config.get(), MESH_0));
}

TEST_CASE("Primitive solver rejects unsupported setup without constructing evolution", "[primitive-solver]") {
  for (const auto& option : {"PRIMITIVE_REPRESENTATION=DUAL_MEAN\n", "PRIMITIVE_FACE_QUADRATURE=SEGMENTED\n",
                             "PRIMITIVE_VOLUME_QUADRATURE=EXACT\n", "PRIMITIVE_INITIAL_STATE=(0,0,0,1)\n"}) {
    INFO(option);
    MeshFixture fixture(points, elements, "PRIMITIVE_EULER", explicitEuler + option);
    CHECK_THROWS_AS(CPrimitiveEulerSolver(fixture.geometry.get(), fixture.config.get(), MESH_0), std::invalid_argument);
    CHECK_THROWS_AS(CSolverFactory::CreateSolverContainer(MAIN_SOLVER::PRIMITIVE_EULER, fixture.config.get(),
                                                          fixture.geometry.get(), MESH_0), std::invalid_argument);
  }
  MeshFixture implicit(points, elements, "PRIMITIVE_EULER", "TIME_DISCRE_FLOW=EULER_IMPLICIT\n");
  CHECK_THROWS_AS(CPrimitiveEulerSolver(implicit.geometry.get(), implicit.config.get(), MESH_0), std::invalid_argument);
  CHECK_THROWS_AS(CPrimitiveEulerSolver(implicit.geometry.get(), implicit.config.get(), 1), std::invalid_argument);
}

TEST_CASE("Primitive scaffold cannot silently report a numerical simulation", "[primitive-solver]") {
  MeshFixture fixture(points, elements, "PRIMITIVE_EULER", explicitEuler);
  CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), MESH_0);
  CSolver* base = &solver;
  std::vector<State> before;
  for (unsigned long i = 0; i < points.size(); ++i) {
    before.push_back(solver.GetPrimitiveNodes().GetPrimitiveState(i));
    for (unsigned k = 0; k < 4; ++k) solver.LinSysRes.GetBlock(i)[k] = 10 * i + k + 1;
  }
  const std::vector<std::function<void()>> calls{
      [&]() { base->Upwind_Residual(nullptr, nullptr, nullptr, nullptr, 0); },
      [&]() { base->Centered_Residual(nullptr, nullptr, nullptr, nullptr, 0, 0); },
      [&]() { base->Source_Residual(nullptr, nullptr, nullptr, nullptr, 0); },
      [&]() { base->SetTime_Step(nullptr, nullptr, nullptr, 0, 0); },
      [&]() { base->ExplicitEuler_Iteration(nullptr, nullptr, nullptr); },
      [&]() { base->ExplicitRK_Iteration(nullptr, nullptr, nullptr, 0); },
      [&]() { base->ClassicalRK4_Iteration(nullptr, nullptr, nullptr, 0); },
      [&]() { base->ImplicitEuler_Iteration(nullptr, nullptr, nullptr); }};
  for (unsigned call = 0; call < calls.size(); ++call) {
    INFO("virtual evolution method " << call);
    try {
      calls[call]();
      FAIL("Evolution method returned without throwing");
    } catch (const std::logic_error& error) {
      CHECK(std::string(error.what()) == "Primitive Euler phase 3: residual, boundary conditions and time integration "
                                        "are not implemented");
    }
    for (unsigned long i = 0; i < points.size(); ++i) {
      CheckState(solver.GetPrimitiveNodes().GetPrimitiveState(i), before[i]);
      for (unsigned k = 0; k < 4; ++k) CHECK(solver.LinSysRes.GetBlock(i)[k] == 10 * i + k + 1);
    }
  }
}

TEST_CASE("Invalid primitive writes have the strong state guarantee", "[primitive-solver]") {
  MeshFixture fixture(points, elements, "PRIMITIVE_EULER", explicitEuler);
  CPrimitiveEulerVariable nodes({{2.3, -0.6, 0.9, 1.7}}, 3, fixture.config.get());
  std::vector<State> invalid;
  for (unsigned component : {0u, 3u})
    for (double value : {0.0, -0.2}) {
      State state{{2.3, -0.6, 0.9, 1.7}};
      state[component] = value;
      invalid.push_back(state);
    }
  for (unsigned component = 0; component < 4; ++component)
    for (double value : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()}) {
      State state{{2.3, -0.6, 0.9, 1.7}};
      state[component] = value;
      invalid.push_back(state);
    }
  for (unsigned attempt = 0; attempt < invalid.size(); ++attempt) {
    INFO("invalid write " << attempt);
    nodes.SetPrimitiveState(1, {{2.3, -0.6, 0.9, 1.7}});
    const auto before = nodes.GetPrimitiveState(1);
    CHECK_THROWS_AS(nodes.SetPrimitiveState(1, invalid[attempt]), std::invalid_argument);
    for (unsigned k = 0; k < 4; ++k) CHECK(nodes.GetSolution(1, k) == before[k]);
  }
  for (unsigned long index : {3ul, std::numeric_limits<unsigned long>::max()}) {
    CHECK_THROWS_AS(nodes.GetPrimitiveState(index), std::out_of_range);
    CHECK_THROWS_AS(nodes.SetPrimitiveState(index, {{1, 0, 0, 1}}), std::out_of_range);
    CHECK_THROWS_AS(nodes.GetConservativeState(index), std::out_of_range);
  }
}

TEST_CASE("Gamma-dependent primitive properties use effective configuration constants", "[primitive-solver]") {
  for (const std::string gas : {"FLUID_MODEL=IDEAL_GAS\nGAMMA_VALUE=1.25\n",
                                "FLUID_MODEL=IDEAL_GAS\nGAMMA_VALUE=1.6\n",
                                "FLUID_MODEL=STANDARD_AIR\nGAMMA_VALUE=1.6\n"}) {
    INFO(gas);
    MeshFixture fixture(points, elements, "PRIMITIVE_EULER", explicitEuler + gas);
    const double gamma = gas.find("STANDARD_AIR") != std::string::npos ? 1.4 :
                         (gas.find("1.25") != std::string::npos ? 1.25 : 1.6);
    CHECK(fixture.config->GetGamma() == Approx(gamma));
    CPrimitiveEulerVariable nodes({{2.3, -0.6, 0.9, 1.7}}, 1, fixture.config.get());
    for (const State state : {State{{2.3, -0.6, 0.9, 1.7}}, State{{0.8, 0.7, -1.2, 2.4}}}) {
      nodes.SetPrimitiveState(0, state);
      const double kinetic = 0.5 * (state[1] * state[1] + state[2] * state[2]);
      const double specific = state[3] / ((gamma - 1) * state[0]) + kinetic;
      CHECK(nodes.GetSoundSpeed(0) == Approx(std::sqrt(gamma * state[3] / state[0])));
      CHECK(nodes.GetEnergy(0) == Approx(specific));
      CHECK(nodes.GetEnthalpy(0) == Approx(specific + state[3] / state[0]));
      CheckState(nodes.GetConservativeState(0),
                 {{state[0], state[0] * state[1], state[0] * state[2], state[3] / (gamma - 1) + state[0] * kinetic}});
      CheckState(nodes.GetPrimitiveState(0), state);
    }
  }
}

TEST_CASE("Primitive arrays preserve node independence aliasing and temporal allocation", "[primitive-solver]") {
  const std::vector<Direction> xy{{{0, 0}}, {{2, 0}}, {{2, 1}}, {{0, 1}}, {{0.8, 0.3}}};
  const std::vector<Triangle> cells{{{0, 1, 4}}, {{1, 2, 4}}, {{2, 3, 4}}, {{3, 0, 4}}};
  for (bool unsteady : {false, true}) {
    const std::string time = unsteady ? "TIME_DOMAIN=YES\nTIME_MARCHING=TIME_STEPPING\nTIME_STEP=0.01\n" : "";
    MeshFixture fixture(xy, cells, "PRIMITIVE_EULER", explicitEuler + time +
                        "NUM_METHOD_GRAD=GREEN_GAUSS\nNUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES\n"
                        "SLOPE_LIMITER_FLOW=VENKATAKRISHNAN\n");
    CHECK(fixture.config->GetKind_Gradient_Method() == GREEN_GAUSS);
    CHECK(fixture.config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
    CHECK(fixture.config->GetKind_SlopeLimit_Flow() == LIMITER::VENKATAKRISHNAN);
    CPrimitiveEulerSolver solver(fixture.geometry.get(), fixture.config.get(), MESH_0);
    auto& nodes = solver.GetPrimitiveNodes();
    CHECK(nodes.GetSolution().rows() == 5);
    CHECK(nodes.GetSolution().cols() == 4);
    CHECK(nodes.GetGradient().length() == 5);
    CHECK(nodes.GetGradient().rows() == 4);
    CHECK(nodes.GetGradient().cols() == 2);
    CHECK(nodes.GetGradient_Reconstruction().length() == 5);
    CHECK(nodes.GetGradient_Reconstruction().rows() == 4);
    CHECK(nodes.GetGradient_Reconstruction().cols() == 2);
    CHECK(nodes.GetLimiter().rows() == 5);
    CHECK(nodes.GetLimiter().cols() == 4);
    CHECK(nodes.GetSolution_time_n().size() == (unsteady ? 20 : 0));
    CHECK(nodes.GetSolution_time_n1().size() == (unsteady ? 20 : 0));
    const State initial = nodes.GetPrimitiveState(0), next{{1.8, -0.4, 0.2, 2.7}};
    for (unsigned long i = 0; i < 5; ++i) {
      CHECK(nodes.GetPrimitive(i) == nodes.GetSolution(i));
      for (unsigned k = 0; k < 4; ++k) {
        CHECK(nodes.GetSolution_Old(i, k) == initial[k]);
        if (unsteady) {
          CHECK(nodes.GetSolution_time_n(i, k) == initial[k]);
          CHECK(nodes.GetSolution_time_n1(i, k) == initial[k]);
        }
        CHECK(nodes.GetLimiter(i, k) == 1);
        for (unsigned d = 0; d < 2; ++d) {
          CHECK(nodes.GetGradient(i, k, d) == 0);
          CHECK(nodes.GetGradient_Reconstruction()(i, k, d) == 0);
        }
      }
      CHECK(nodes.GetDelta_Time(i) == 0);
    }
    nodes.SetPrimitiveState(2, next);
    for (unsigned long i = 0; i < 5; ++i) CheckState(nodes.GetPrimitiveState(i), i == 2 ? next : initial);
    nodes.GetPrimitive(2)[1] = 0.33;
    CHECK(nodes.GetSolution(2, 1) == Approx(0.33));
    nodes.SetSolution(2, 2, -0.44);
    CHECK(nodes.GetPrimitive(2, 2) == Approx(-0.44));
    nodes.GetSolution(2)[3] = -1;  // Deliberate demonstration of the unchecked SU2 raw interface.
    CHECK(nodes.GetPrimitive(2, 3) == -1);
    CHECK_THROWS_AS(nodes.GetConservativeState(2), std::invalid_argument);
    nodes.SetPrimitiveState(2, next);
    for (unsigned long i = 0; i < 5; ++i)
      for (unsigned k = 0; k < 4; ++k) CHECK(nodes.GetSolution_Old(i, k) == initial[k]);
  }
}

TEST_CASE("Factory unwinds failures after geometry allocation and remains usable",
          "[primitive-solver][primitive-factory-unwind]") {
  MeshFixture bad(points, elements, "PRIMITIVE_EULER", explicitEuler + "PRIMITIVE_INITIAL_STATE=(1,0,0,-1)\n");
  // Metrics have already been constructed when variable validation throws.
  // Repetition is an exercise for memory tools, not a leak counter/proof.
  for (unsigned attempt = 0; attempt < 32; ++attempt)
    CHECK_THROWS_AS(CSolverFactory::CreateSolverContainer(MAIN_SOLVER::PRIMITIVE_EULER, bad.config.get(),
                                                          bad.geometry.get(), 0), std::invalid_argument);
  CHECK(bad.geometry->GetnPoint() == 3);
  MeshFixture good(points, elements, "PRIMITIVE_EULER", explicitEuler);
  FactoryContainer result(good);
  auto* solver = dynamic_cast<CPrimitiveEulerSolver*>(result.solvers[FLOW_SOL]);
  REQUIRE(solver != nullptr);
  CheckState(solver->GetPrimitiveNodes().GetPrimitiveState(0), {{1, 0.5, 0.3, 1}});
}

TEST_CASE("Additional unsupported configurations fail at solver validation", "[primitive-solver]") {
  for (const std::string options : {"RESTART_SOL=YES\n", "AXISYMMETRIC=YES\n", "KIND_TURB_MODEL=SA\n",
                                    "FLUID_MODEL=VW_GAS\n",
                                    "TIME_DOMAIN=YES\nTIME_MARCHING=DUAL_TIME_STEPPING-1ST_ORDER\nTIME_STEP=0.01\n"}) {
    INFO(options);
    std::stringstream input("SOLVER=PRIMITIVE_EULER\nMGLEVEL=0\nCONV_NUM_METHOD_FLOW=ROE\n" + explicitEuler + options);
    CConfig config(input, SU2_COMPONENT::SU2_CFD, false);
    CHECK_THROWS_AS(CPrimitiveEulerSolver::ValidateConfiguration(config), std::invalid_argument);
  }
  MeshFixture fixture(points, elements, "PRIMITIVE_EULER", explicitEuler);
  fixture.config->SetMGLevels(1);  // Public setter: no coarse geometry is constructed.
  CHECK(fixture.config->GetnMGLevels() == 1);
  CHECK_THROWS_AS(CPrimitiveEulerSolver::ValidateConfiguration(*fixture.config), std::invalid_argument);
  fixture.config->SetMGLevels(0);
  fixture.config->SetMultizone_Problem(true);  // Validate only; do not initialize coupled arrays.
  CHECK_THROWS_AS(CPrimitiveEulerSolver::ValidateConfiguration(*fixture.config), std::invalid_argument);
}
