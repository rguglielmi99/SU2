/*!
 * \file primitive_euler_geometry_tests.cpp
 * \brief Checks against SU2's actual median-dual construction.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "../../SU2_CFD/include/geometry/primitive_euler_geometry.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/geometry/primal_grid/CTriangle.hpp"
#include "../../Common/include/geometry/primal_grid/CQuadrilateral.hpp"
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <map>
#include <set>

#include "primitive_euler_mesh_fixture.hpp"

using namespace PrimitiveEuler2D;
using namespace PrimitiveEulerTest;

namespace {
// Independent shoelace area/centroid of an explicitly assembled dual polygon.
struct PolygonMetrics {
  double area;
  Direction centroid;
  Direction offset;  // Relative to vertices.front(), without subtracting large coordinates.
};
PolygonMetrics Polygon(const std::vector<Direction>& vertices) {
  REQUIRE(vertices.size() >= 3);
  long double twiceArea = 0, mx = 0, my = 0;
  const auto& origin = vertices.front();
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const auto& p = vertices[i];
    const auto& q = vertices[(i + 1) % vertices.size()];
    const long double px = static_cast<long double>(p[0]) - origin[0];
    const long double py = static_cast<long double>(p[1]) - origin[1];
    const long double qx = static_cast<long double>(q[0]) - origin[0];
    const long double qy = static_cast<long double>(q[1]) - origin[1];
    const long double cross = px * qy - qx * py;
    twiceArea += cross;
    mx += (px + qx) * cross;
    my += (py + qy) * cross;
  }
  REQUIRE(std::isfinite(twiceArea));
  REQUIRE(std::abs(twiceArea) > 0);
  const Direction offset{{static_cast<double>(mx / (3 * twiceArea)), static_cast<double>(my / (3 * twiceArea))}};
  return {static_cast<double>(0.5L * std::abs(twiceArea)), {{origin[0] + offset[0], origin[1] + offset[1]}}, offset};
}

void Check(double actual, double expected, double scale = 1.0) {
  REQUIRE(std::isfinite(actual));
  REQUIRE(std::isfinite(expected));
  REQUIRE(std::isfinite(scale));
  REQUIRE(scale > 0);
  CHECK(std::abs(actual - expected) <= 2e-12 * std::max(scale, std::abs(expected)));
}

const std::vector<Triangle> fan{{{0, 1, 4}}, {{1, 2, 4}}, {{2, 3, 4}}, {{3, 0, 4}}};

struct DualReference {
  std::vector<double> areas;
  std::vector<Direction> offsets;
};

DualReference PolygonReference(const std::vector<Direction>& xy, const std::vector<Triangle>& triangles) {
  DualReference reference{std::vector<double>(xy.size(), 0), std::vector<Direction>(xy.size())};
  // Every subpolygon is expressed relative to its own primal vertex. Thus the
  // reference never forms area*absolute-coordinate or subtracts two centroids.
  for (const auto& t : triangles) {
    for (unsigned k = 0; k < 3; ++k) {
      const auto i = t[k], j = t[(k + 1) % 3], previous = t[(k + 2) % 3];
      Direction nextMid{}, prevMid{}, c{};
      for (unsigned d = 0; d < 2; ++d) {
        nextMid[d] = 0.5 * (xy[j][d] - xy[i][d]);
        prevMid[d] = 0.5 * (xy[previous][d] - xy[i][d]);
        c[d] = ((xy[j][d] - xy[i][d]) + (xy[previous][d] - xy[i][d])) / 3;
      }
      const auto polygon = Polygon({{{0, 0}}, nextMid, c, prevMid});
      reference.areas[i] += polygon.area;
      for (unsigned d = 0; d < 2; ++d) reference.offsets[i][d] += polygon.area * polygon.offset[d];
    }
  }
  for (unsigned long i = 0; i < xy.size(); ++i) {
    REQUIRE(reference.areas[i] > 0);
    for (unsigned d = 0; d < 2; ++d) reference.offsets[i][d] /= reference.areas[i];
  }
  return reference;
}

void CheckPolygonReference(const std::vector<Direction>& xy, const std::vector<Triangle>& triangles,
                           const CAggregatedGeometry& cache, double length = 1.0) {
  const auto reference = PolygonReference(xy, triangles);
  for (unsigned long i = 0; i < xy.size(); ++i) {
    const auto& cell = cache.GetCell(i);
    Check(cell.volume, reference.areas[i], length * length);
    for (unsigned d = 0; d < 2; ++d) Check(cell.offset[d], reference.offsets[i][d], length);
    // Two independent affine fields. Check the correction itself so a large
    // nodal value cannot mask an incorrect small geometric mean correction.
    for (const auto& field : {State{{1.7, 0.3, -0.2, 0}}, State{{-0.4, -0.7, 0.6, 0}}}) {
      const double nodal = field[0] + field[1] * xy[i][0] + field[2] * xy[i][1];
      const double correction = field[1] * cell.offset[0] + field[2] * cell.offset[1];
      const double expected = field[1] * reference.offsets[i][0] + field[2] * reference.offsets[i][1];
      Check(correction, expected, length);
      Check(nodal + correction, nodal + expected);
    }
  }
}

// CCW rectangular boundary 0..5, with non-corner boundary nodes 1 and 4.
// Three strict interior nodes 6,7,8; their edge valences are 6,6,3.
// Node 8 splits only triangle (1,7,6), not the entire mesh into one fan.
const std::vector<Direction> irregularXY{{{0, 0}}, {{1.1, 0}}, {{3, 0}}, {{3, 2}}, {{1.7, 2}},
                                         {{0, 2}}, {{0.8, 0.8}}, {{2.2, 1.1}}, {{1.4, 0.55}}};
const std::vector<Triangle> irregularTriangles{{{0, 1, 6}}, {{1, 2, 7}}, {{1, 7, 8}}, {{7, 6, 8}},
                                               {{6, 1, 8}}, {{2, 3, 7}}, {{3, 4, 7}}, {{4, 5, 6}},
                                               {{4, 6, 7}}, {{5, 0, 6}}};

void CheckDomainMoments(const std::vector<Direction>& xy, const CAggregatedGeometry& cache, double length = 1.0) {
  // Reference from the outer polygon only, not from duals or triangle areas.
  const auto domain = Polygon({xy[0], xy[1], xy[2], xy[3], xy[4], xy[5]});
  long double area = 0, mx = 0, my = 0;
  for (unsigned long i = 0; i < xy.size(); ++i) {
    const auto& cell = cache.GetCell(i);
    area += cell.volume;
    mx += cell.volume * (static_cast<long double>(xy[i][0]) - xy[0][0] + cell.offset[0]);
    my += cell.volume * (static_cast<long double>(xy[i][1]) - xy[0][1] + cell.offset[1]);
  }
  Check(static_cast<double>(area), domain.area, length * length);
  // First moments relative to xy[0]; together with area these also determine
  // the absolute first moments, without cancellation under translations.
  Check(static_cast<double>(mx), domain.area * domain.offset[0], length * length * length);
  Check(static_cast<double>(my), domain.area * domain.offset[1], length * length * length);
  Check(static_cast<double>(mx + area * xy[0][0]), domain.area * domain.centroid[0]);
  Check(static_cast<double>(my + area * xy[0][1]), domain.area * domain.centroid[1]);
}

void CheckIrregularClosure(const CGeometry& geometry, const CAggregatedGeometry& cache) {
  std::vector<Direction> closure(cache.GetnPoint()), nativeBoundary(cache.GetnPoint());
  std::vector<Direction> expectedBoundary(cache.GetnPoint());
  std::vector<double> lengths(cache.GetnPoint(), 0);
  for (unsigned long e = 0; e < cache.GetnEdge(); ++e) {
    const auto& face = cache.GetFace(e);
    CHECK(geometry.edges->GetNode(e, 0) < geometry.edges->GetNode(e, 1));
    const double measure = std::hypot(face.normal[0], face.normal[1]);
    for (unsigned d = 0; d < 2; ++d) {
      closure[face.first][d] += face.normal[d];
      closure[face.second][d] -= face.normal[d];
    }
    lengths[face.first] += measure;
    lengths[face.second] += measure;
  }
  for (unsigned marker = 0; marker < geometry.GetnMarker(); ++marker)
    for (unsigned long v = 0; v < geometry.GetnVertex(marker); ++v) {
      const auto* vertex = geometry.vertex[marker][v];
      for (unsigned d = 0; d < 2; ++d) nativeBoundary[vertex->GetNode()][d] += vertex->GetNormal(d);
    }
  // For a CCW boundary edge the inward half-normal is (-dy,dx)/2.
  // This independent contour reference checks the actual sign, not just closure.
  for (unsigned i = 0; i < 6; ++i) {
    const unsigned j = (i + 1) % 6;
    const Direction inward{{-0.5 * (irregularXY[j][1] - irregularXY[i][1]),
                              0.5 * (irregularXY[j][0] - irregularXY[i][0])}};
    for (unsigned node : {i, j})
      for (unsigned d = 0; d < 2; ++d) expectedBoundary[node][d] += inward[d];
  }
  // Check classes separately; counts make omission of any category visible.
  const std::array<std::vector<unsigned>, 3> groups{{{6, 7, 8}, {1, 4}, {0, 2, 3, 5}}};
  for (unsigned group = 0; group < groups.size(); ++group)
    for (unsigned i : groups[group]) {
      INFO("0=interior, 1=straight boundary, 2=corner; group=" << group << " node=" << i);
      for (unsigned d = 0; d < 2; ++d) {
        Check(nativeBoundary[i][d], expectedBoundary[i][d], lengths[i]);
        Check(closure[i][d] - nativeBoundary[i][d], 0, lengths[i]);
        if (group == 0) Check(closure[i][d], 0, lengths[i]);
      }
    }
}
}  // namespace

TEST_CASE("Primitive aggregated geometry reuses actual SU2 triangular dual metrics", "[primitive-geometry]") {
  const std::vector<Direction> xy{{{0, 0}}, {{2, 0}}, {{0, 1}}};
  const std::vector<Triangle> triangles{{{0, 1, 2}}};
  MeshFixture fixture(xy, triangles);
  auto& geometry = *fixture.geometry;
  auto& config = fixture.config;
  const CAggregatedGeometry cache(geometry, *config);
  REQUIRE(cache.GetnPoint() == 3);
  REQUIRE(cache.GetnEdge() == 3);
  CheckPolygonReference(xy, triangles, cache);
  for (unsigned long i = 0; i < 3; ++i) Check(cache.GetCell(i).volume, 1.0 / 3);
  for (unsigned long e = 0; e < cache.GetnEdge(); ++e) {
    const auto& face = cache.GetFace(e);
    for (unsigned d = 0; d < 2; ++d) {
      CHECK(face.normal[d] == geometry.edges->GetNormal(e)[d]);
      Check(face.point[d], 0.5 * (xy[face.first][d] + xy[face.second][d]));
    }
  }
  // Exact normal of the dual face on primal edge 0->1: C-M=(-1/3,1/3).
  const auto& face = cache.GetFace(geometry.FindEdge(0, 1));
  Check(face.normal[0], 1.0 / 3);
  Check(face.normal[1], 1.0 / 3);
  CHECK(std::hypot(cache.GetCell(0).offset[0], cache.GetCell(0).offset[1]) > 0.1);
}

TEST_CASE("Regular and distorted triangular fans preserve area and first moments", "[primitive-geometry]") {
  for (const Direction center : {Direction{{1, 0.5}}, Direction{{0.8, 0.3}}}) {
    const std::vector<Direction> xy{{{0, 0}}, {{2, 0}}, {{2, 1}}, {{0, 1}}, center};
    MeshFixture fixture(xy, fan);
    auto& geometry = *fixture.geometry;
    auto& config = fixture.config;
    const CAggregatedGeometry cache(geometry, *config);
    CheckPolygonReference(xy, fan, cache);
    double area = 0;
    Direction moment{}, closure{};
    std::vector<Direction> cellClosure(xy.size());
    for (unsigned long i = 0; i < xy.size(); ++i) {
      const auto& cell = cache.GetCell(i);
      area += cell.volume;
      for (unsigned d = 0; d < 2; ++d) moment[d] += cell.volume * (xy[i][d] + cell.offset[d]);
    }
    Check(area, 2);
    Check(moment[0], 2);
    Check(moment[1], 1);
    // Interior dual is closed: use original integrated normals with endpoint signs.
    for (unsigned long e = 0; e < cache.GetnEdge(); ++e) {
      const auto& f = cache.GetFace(e);
      for (unsigned d = 0; d < 2; ++d) {
        cellClosure[f.first][d] += f.normal[d];
        cellClosure[f.second][d] -= f.normal[d];
        if (f.first == 4) closure[d] += f.normal[d];
        if (f.second == 4) closure[d] -= f.normal[d];
      }
    }
    Check(closure[0], 0);
    Check(closure[1], 0);
    // Boundary vertex normals in SU2 are INWARD. Subtract them to close each
    // dual, including corners where contributions from two boundary edges meet.
    for (unsigned marker = 0; marker < geometry.GetnMarker(); ++marker)
      for (unsigned long v = 0; v < geometry.GetnVertex(marker); ++v) {
        const auto* vertex = geometry.vertex[marker][v];
        for (unsigned d = 0; d < 2; ++d) cellClosure[vertex->GetNode()][d] -= vertex->GetNormal(d);
      }
    for (const auto& normal : cellClosure) {
      Check(normal[0], 0);
      Check(normal[1], 0);
    }
  }
}

TEST_CASE("Renumbering nodes reverses normals without changing physical geometry", "[primitive-geometry]") {
  const std::vector<Direction> xy{{{0, 0}}, {{2, 0}}, {{2, 1}}, {{0, 1}}, {{0.8, 0.3}}};
  const std::array<unsigned long, 5> permutation{{3, 1, 4, 0, 2}};
  std::vector<Direction> renumbered(5);
  for (unsigned i = 0; i < 5; ++i) renumbered[permutation[i]] = xy[i];
  auto reorderedTriangles = fan;
  for (auto& t : reorderedTriangles)
    for (auto& i : t) i = permutation[i];
  MeshFixture first(xy, fan), second(renumbered, reorderedTriangles);
  auto& original = *first.geometry;
  auto& reordered = *second.geometry;
  const CAggregatedGeometry a(original, *first.config), b(reordered, *second.config);
  for (unsigned i = 0; i < 5; ++i) {
    Check(a.GetCell(i).volume, b.GetCell(permutation[i]).volume);
    for (unsigned d = 0; d < 2; ++d) Check(a.GetCell(i).offset[d], b.GetCell(permutation[i]).offset[d]);
  }
  for (unsigned long e = 0; e < a.GetnEdge(); ++e) {
    const auto& f = a.GetFace(e);
    const auto& g = b.GetFace(reordered.FindEdge(permutation[f.first], permutation[f.second]));
    const double sign = permutation[f.first] < permutation[f.second] ? 1 : -1;
    for (unsigned d = 0; d < 2; ++d) {
      Check(f.point[d], g.point[d]);
      Check(f.normal[d], sign * g.normal[d]);
    }
    const State wi{{1, 0.2, -0.1, 1.3}}, wj{{1.2, -0.3, 0.2, 0.8}};
    const auto before = ComputeFluctuations(wi, wj, f.normal, 1.4);
    const auto after = sign > 0 ? ComputeFluctuations(wi, wj, g.normal, 1.4)
                                : ComputeFluctuations(wj, wi, g.normal, 1.4);
    for (unsigned k = 0; k < 4; ++k) {
      Check(before.minus[k], sign > 0 ? after.minus[k] : after.plus[k]);
      Check(before.plus[k], sign > 0 ? after.plus[k] : after.minus[k]);
    }
  }
}

TEST_CASE("Geometry snapshot rejects corrupted native normals and volumes", "[primitive-geometry]") {
  MeshFixture fixture({{{0, 0}}, {{2, 0}}, {{0, 1}}}, {{{0, 1, 2}}});
  auto& geometry = *fixture.geometry;
  auto& config = fixture.config;
  const auto native = geometry.nodes->GetVolume(0);
  geometry.nodes->SetVolume(0, 2 * native);
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
  geometry.nodes->SetVolume(0, native);
  const Direction saved{{geometry.edges->GetNormal(0)[0], geometry.edges->GetNormal(0)[1]}};
  geometry.edges->SetNormal(0, Direction{{-saved[0], -saved[1]}});
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
  geometry.edges->SetNormal(0, saved);
  geometry.edges->SetNormal(0, Direction{{1.25 * saved[0], 1.25 * saved[1]}});
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
  geometry.edges->SetNormal(0, saved);
  CHECK_NOTHROW(CAggregatedGeometry(geometry, *config));
}

TEST_CASE("Initial geometry scope rejects unsupported mesh modes", "[primitive-geometry]") {
  MeshFixture fixture({{{0, 0}}, {{2, 0}}, {{0, 1}}}, {{{0, 1, 2}}});
  auto& geometry = *fixture.geometry;
  auto& config = fixture.config;
  geometry.SetnDim(3);
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
  geometry.SetnDim(2);
  geometry.SetMGLevel(1);
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
  geometry.SetMGLevel(0);
  geometry.elem[0]->Change_Orientation();
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
  geometry.elem[0]->Change_Orientation();
  delete geometry.elem[0];
  geometry.elem[0] = new CQuadrilateral(0, 1, 2, 0);
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
  geometry.SetnPointDomain(geometry.GetnPoint() - 1);
  CHECK_THROWS_AS(CAggregatedGeometry(geometry, *config), std::invalid_argument);
}

TEST_CASE("Irregular triangular mesh with three interiors has consistent dual geometry", "[primitive-geometry]") {
  std::vector<std::set<unsigned long>> neighbors(irregularXY.size());
  double smallest = std::numeric_limits<double>::max(), largest = 0;
  for (const auto& t : irregularTriangles) {
    const auto& a = irregularXY[t[0]];
    const auto& b = irregularXY[t[1]];
    const auto& c = irregularXY[t[2]];
    const double area = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
    REQUIRE(area > 0);
    smallest = std::min(smallest, area);
    largest = std::max(largest, area);
    for (unsigned k = 0; k < 3; ++k) {
      neighbors[t[k]].insert(t[(k + 1) % 3]);
      neighbors[t[k]].insert(t[(k + 2) % 3]);
    }
  }
  CHECK(neighbors[6].size() == 6);
  CHECK(neighbors[7].size() == 6);
  CHECK(neighbors[8].size() == 3);
  CHECK(largest > 3 * smallest);
  for (unsigned i : {6u, 7u, 8u}) {
    CHECK(irregularXY[i][0] > 0);
    CHECK(irregularXY[i][0] < 3);
    CHECK(irregularXY[i][1] > 0);
    CHECK(irregularXY[i][1] < 2);
  }
  MeshFixture fixture(irregularXY, irregularTriangles);
  const CAggregatedGeometry cache(*fixture.geometry, *fixture.config);
  REQUIRE(cache.GetnPoint() == 9);
  REQUIRE(cache.GetnEdge() == 18);  // V-E+F=9-18+10=1, one connected triangulated disk.
  CheckPolygonReference(irregularXY, irregularTriangles, cache);
  CheckDomainMoments(irregularXY, cache);
  CheckIrregularClosure(*fixture.geometry, cache);
}

TEST_CASE("Relative shoelace reference is accurate far from the coordinate origin", "[primitive-geometry]") {
  const auto local = Polygon({{{0, 0}}, {{2, 0}}, {{2, 1}}, {{0, 1}}});
  const double x = 1e9, y = -1e9;
  const auto moved = Polygon({{{x, y}}, {{x + 2, y}}, {{x + 2, y + 1}}, {{x, y + 1}}});
  Check(moved.area, 2);
  Check(moved.offset[0], 1);
  Check(moved.offset[1], 0.5);
  Check(moved.area, local.area);
  for (unsigned d = 0; d < 2; ++d) Check(moved.offset[d], local.offset[d]);
}

TEST_CASE("Triangular dual metrics obey translation and uniform scaling", "[primitive-geometry]") {
  MeshFixture base(irregularXY, irregularTriangles);
  const CAggregatedGeometry original(*base.geometry, *base.config);
  // Small unshifted mesh and moderately translated large mesh: test distinct
  // physical scales without exceeding the precision of native SU2 metrics.
  for (const auto& transform : {State{{1e-3, 0, 0, 0}}, State{{8, 128, -256, 0}}}) {
    const double scale = transform[0];
    auto xy = irregularXY;
    for (auto& point : xy)
      for (unsigned d = 0; d < 2; ++d) point[d] = scale * point[d] + transform[d + 1];
    MeshFixture moved(xy, irregularTriangles);
    const CAggregatedGeometry metrics(*moved.geometry, *moved.config);
    CheckPolygonReference(xy, irregularTriangles, metrics, scale);
    CheckDomainMoments(xy, metrics, scale);
    for (unsigned long i = 0; i < xy.size(); ++i) {
      Check(metrics.GetCell(i).volume, scale * scale * original.GetCell(i).volume, scale * scale);
      for (unsigned d = 0; d < 2; ++d)
        Check(metrics.GetCell(i).offset[d], scale * original.GetCell(i).offset[d], scale);
    }
    for (unsigned long e = 0; e < metrics.GetnEdge(); ++e) {
      const auto& f = original.GetFace(e);
      const auto& g = metrics.GetFace(moved.geometry->FindEdge(f.first, f.second));
      for (unsigned d = 0; d < 2; ++d) {
        Check(g.normal[d], scale * f.normal[d], scale);
        // Compare relative Qeff, so translation does not weaken the tolerance.
        Check(g.point[d] - transform[d + 1], scale * f.point[d], scale);
      }
    }
  }
}

TEST_CASE("Geometry validation safely rejects invalid primal coordinates and edge incidence", "[primitive-geometry]") {
  MeshFixture fixture({{{0, 0}}, {{2, 0}}, {{0, 1}}}, {{{0, 1, 2}}});
  auto& geometry = *fixture.geometry;
  // Mutate an already initialized valid geometry: do not feed invalid data to
  // the mesh reader/preprocessing, which has a different error contract.
  SECTION("Collinear triangle") {
    geometry.nodes->SetCoord(2, 0, 1.0);
    geometry.nodes->SetCoord(2, 1, 0.0);
    CHECK_THROWS_WITH(CAggregatedGeometry(geometry, *fixture.config),
                      "Degenerate/inverted triangle: run SU2 orientation checks first");
  }
  SECTION("Nonfinite coordinate") {
    geometry.nodes->SetCoord(0, 0, std::numeric_limits<double>::quiet_NaN());
    CHECK_THROWS_WITH(CAggregatedGeometry(geometry, *fixture.config), "Non-finite mesh coordinate");
  }
  SECTION("Three triangles incident on one edge") {
    // Duplicate connectivity only after SU2 preprocessing; existing edge
    // lookups remain valid and the builder reaches its incidence guard safely.
    // Use a separate fixture with enough element slots; no allocation hacking.
    MeshFixture fanFixture(irregularXY, irregularTriangles);
    auto& fanGeometry = *fanFixture.geometry;
    for (unsigned slot : {1u, 2u}) {
      delete fanGeometry.elem[slot];
      fanGeometry.elem[slot] = new CTriangle(0, 1, 6);
    }
    CHECK_THROWS_WITH(CAggregatedGeometry(fanGeometry, *fanFixture.config), "Non-manifold triangle edge");
  }
}
