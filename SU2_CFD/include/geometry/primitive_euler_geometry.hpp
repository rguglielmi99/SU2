/*!
 * \file primitive_euler_geometry.hpp
 * \brief Static serial triangular geometry for primitive Euler, AGGREGATED/NODE.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../numerics/flow/primitive_euler_2d.hpp"
#include <vector>

namespace PrimitiveEuler2D {

/*!
 * \brief Read-only snapshot built AFTER SU2 SetEdges/SetControlVolume.
 * Uses SU2's integrated normals and volumes, without changing CGeometry.
 * Stores no segment quadrature. NODE states remain values at primal vertices.
 * Centroid offsets support a later mean of the affine nodal reconstruction;
 * they do not change the unknown into a DUAL_MEAN state.
 * Rebuild after any geometry change. Initial scope: static, serial, 2D triangles,
 * finest mesh, without periodic identifications, MPI halos or AD.
 */
class CAggregatedGeometry {
 public:
  struct Face {
    unsigned long first, second; // The local indices of the nodes involved
    Direction point;   // Qeff=(xi+xj)/2.
    Direction normal;  // Integrated Sij, oriented first -> second; NOT unit length.
  };
  struct Cell {
    double volume;     // Original SU2 dual volume.
    Direction offset;  // Geometric dual centroid minus primal node position.
  };
  struct BoundaryFace {
    unsigned long node;
    unsigned short marker;
    Direction point;   // Midpoint of the boundary half-edge owned by node.
    Direction normal;  // Outward integrated normal of that half-edge.
  };

 private:
  std::vector<Face> faces;
  std::vector<Cell> cells;
  std::vector<BoundaryFace> boundaryFaces;

  static void RequireClose(double actual, double reference, double scale, const char* message) {
    if (!std::isfinite(actual) || !std::isfinite(reference) ||
        std::abs(actual - reference) > 1e-11 * scale)
      throw std::invalid_argument(message);
  }

 public:
  CAggregatedGeometry(const CGeometry& geometry, const CConfig& config) {
    if (SU2_MPI::GetSize() != 1 || geometry.GetnPoint() != geometry.GetnPointDomain())
      throw std::invalid_argument("Primitive geometry requires serial geometry without halos");
    if (geometry.GetnDim() != 2 || geometry.GetMGLevel() != MESH_0 || config.GetDynamic_Grid() ||
        config.GetnMarker_Periodic() != 0 || config.GetDiscrete_Adjoint() || config.GetContinuous_Adjoint())
      throw std::invalid_argument("Primitive geometry requires static 2D finest-grid nonperiodic primal geometry");
    if (!geometry.nodes || !geometry.edges || !geometry.elem || geometry.GetnElem() == 0)
      throw std::invalid_argument("Primitive geometry requires initialized connectivity and dual metrics");

    const auto count = geometry.GetnPoint();
    std::vector<Direction> coordinates(count), moments(count), normalSum(geometry.GetnEdge());
    std::vector<double> volumes(count, 0.0), normalScale(geometry.GetnEdge(), 0.0);
    std::vector<unsigned> incidence(geometry.GetnEdge(), 0);
    for (unsigned long i = 0; i < count; ++i) {
      for (unsigned d = 0; d < 2; ++d) {
        coordinates[i][d] = SU2_TYPE::GetValue(geometry.nodes->GetCoord(i)[d]);
        if (!std::isfinite(coordinates[i][d])) throw std::invalid_argument("Non-finite mesh coordinate");
      }
    }
    faces.reserve(geometry.GetnEdge());
    for (unsigned long e = 0; e < geometry.GetnEdge(); ++e) {
      const auto i = geometry.edges->GetNode(e, 0), j = geometry.edges->GetNode(e, 1);
      if (i >= count || j >= count || i >= j) throw std::invalid_argument("Invalid SU2 local edge ordering");
      const auto* s = geometry.edges->GetNormal(e);
      faces.push_back({i, j, {{0.5 * coordinates[i][0] + 0.5 * coordinates[j][0],
                              0.5 * coordinates[i][1] + 0.5 * coordinates[j][1]}},
                      {{SU2_TYPE::GetValue(s[0]), SU2_TYPE::GetValue(s[1])}}});
    }
    for (unsigned long element = 0; element < geometry.GetnElem(); ++element) {
      const auto* triangle = geometry.elem[element];
      if (!triangle || triangle->GetVTK_Type() != TRIANGLE)
        throw std::invalid_argument("Primitive geometry currently supports triangles only");
      const auto a = triangle->GetNode(0), b = triangle->GetNode(1), c = triangle->GetNode(2);
      if (a >= count || b >= count || c >= count) throw std::invalid_argument("Invalid triangle node");
      const double dx1 = coordinates[b][0] - coordinates[a][0], dy1 = coordinates[b][1] - coordinates[a][1];
      const double dx2 = coordinates[c][0] - coordinates[a][0], dy2 = coordinates[c][1] - coordinates[a][1];
      const double cross = dx1 * dy2 - dy1 * dx2;
      if (!std::isfinite(cross) || cross <= 0.0)
        throw std::invalid_argument("Degenerate/inverted triangle: run SU2 orientation checks first");
      const Direction center{{(coordinates[a][0] + coordinates[b][0] + coordinates[c][0]) / 3,
                               (coordinates[a][1] + coordinates[b][1] + coordinates[c][1]) / 3}};
      for (unsigned f = 0; f < 3; ++f) {
        const auto i = triangle->GetNode(triangle->GetFaces(f, 0));
        const auto j = triangle->GetNode(triangle->GetFaces(f, 1));
        const long edge = geometry.FindEdge(i, j, false);
        if (edge < 0 || static_cast<unsigned long>(edge) >= faces.size())
          throw std::invalid_argument("Missing triangle edge in SU2 connectivity");
        const auto& midpoint = faces[edge].point;
        const double sign = i < j ? 1.0 : -1.0;
        const Direction contribution{{sign * (center[1] - midpoint[1]), -sign * (center[0] - midpoint[0])}};
        for (unsigned d = 0; d < 2; ++d) normalSum[edge][d] += contribution[d];
        normalScale[edge] += std::hypot(contribution[0], contribution[1]);
        if (++incidence[edge] > 2) throw std::invalid_argument("Non-manifold triangle edge");
        // Same subtriangles (node, edge midpoint, element center) as SU2.
        // First moments are accumulated relative to each primal node.
        for (auto node : {i, j}) {
          const double mx = midpoint[0] - coordinates[node][0], my = midpoint[1] - coordinates[node][1];
          const double cx = center[0] - coordinates[node][0], cy = center[1] - coordinates[node][1];
          const double area = 0.5 * std::abs(mx * cy - my * cx);
          volumes[node] += area;
          moments[node][0] += area * (mx + cx) / 3;
          moments[node][1] += area * (my + cy) / 3;
        }
      }
    }
    for (unsigned long e = 0; e < faces.size(); ++e) {
      if (incidence[e] == 0) throw std::invalid_argument("Edge without triangle incidence");
      for (unsigned d = 0; d < 2; ++d)
        RequireClose(faces[e].normal[d], normalSum[e][d], normalScale[e],
                     "SU2 integrated normal inconsistent with triangular dual");
      const auto& face = faces[e];
      const double area = std::hypot(face.normal[0], face.normal[1]);
      const double projection = face.normal[0] * (coordinates[face.second][0] - coordinates[face.first][0]) +
                                face.normal[1] * (coordinates[face.second][1] - coordinates[face.first][1]);
      if (!std::isfinite(area) || area <= 0 || !std::isfinite(projection) || projection <= 0)
        throw std::invalid_argument("Invalid first-to-second integrated normal");
    }
    cells.reserve(count);
    for (unsigned long i = 0; i < count; ++i) {
      const double nativeVolume = SU2_TYPE::GetValue(geometry.nodes->GetVolume(i));
      if (!std::isfinite(volumes[i]) || volumes[i] <= 0 || !std::isfinite(nativeVolume) || nativeVolume <= 0)
        throw std::invalid_argument("Invalid dual volume or isolated point");
      RequireClose(nativeVolume, volumes[i], volumes[i], "SU2 dual volume inconsistent with subtriangle areas");
      const Direction offset{{moments[i][0] / volumes[i], moments[i][1] / volumes[i]}};
      if (!std::isfinite(offset[0]) || !std::isfinite(offset[1]))
        throw std::invalid_argument("Invalid dual centroid");
      cells.push_back({nativeVolume, offset});
    }
    // Straight triangular boundary edges. Determine outward orientation from
    // the incident triangle, rather than assuming a boundary node ordering.
    std::vector<unsigned long> opposite(faces.size(), count);
    for (unsigned long element = 0; element < geometry.GetnElem(); ++element) {
      const auto* triangle = geometry.elem[element];
      for (unsigned k = 0; k < 3; ++k) {
        const auto i = triangle->GetNode(k), j = triangle->GetNode((k + 1) % 3);
        const auto e = geometry.FindEdge(i, j, false);
        if (e < 0 || static_cast<unsigned long>(e) >= faces.size())
          throw std::invalid_argument("Missing boundary triangle edge in SU2 connectivity");
        if (incidence[e] == 1) opposite[e] = triangle->GetNode((k + 2) % 3);
      }
    }
    std::vector<bool> seen(faces.size(), false);
    for (unsigned short marker = 0; marker < geometry.GetnMarker(); ++marker) {
      std::vector<Direction> normalByNode(count);
      for (unsigned long b = 0; b < geometry.GetnElem_Bound(marker); ++b) {
        const auto* edge = geometry.bound[marker][b];
        if (!edge || edge->GetVTK_Type() != LINE)
          throw std::invalid_argument("Primitive boundary requires straight line elements");
        const auto i = edge->GetNode(0), j = edge->GetNode(1);
        if (i >= count || j >= count) throw std::invalid_argument("Invalid primitive boundary node");
        const auto e = geometry.FindEdge(i, j, false);
        if (e < 0 || static_cast<unsigned long>(e) >= faces.size() || incidence[e] != 1 || seen[e])
          throw std::invalid_argument("Primitive boundary must cover each exterior edge once");
        seen[e] = true;
        Direction normal{{0.5 * (coordinates[j][1] - coordinates[i][1]),
                          -0.5 * (coordinates[j][0] - coordinates[i][0])}};
        const double measure = std::hypot(normal[0], normal[1]);
        if (!std::isfinite(measure) || measure <= 0)
          throw std::invalid_argument("Invalid primitive physical boundary normal");
        const auto& midpoint = faces[e].point;
        const auto& inside = coordinates[opposite[e]];
        if (normal[0] * (inside[0] - midpoint[0]) + normal[1] * (inside[1] - midpoint[1]) > 0)
          for (auto& value : normal) value = -value;
        for (auto node : {i, j}) {
          boundaryFaces.push_back({node, marker,
              {{0.5 * (coordinates[node][0] + midpoint[0]), 0.5 * (coordinates[node][1] + midpoint[1])}}, normal});
          for (unsigned d = 0; d < 2; ++d) normalByNode[node][d] += normal[d];
        }
      }
      // SU2 boundary vertex normals point INWARD, opposite to the ones above.
      for (unsigned long v = 0; v < geometry.GetnVertex(marker); ++v) {
        const auto* vertex = geometry.vertex[marker][v];
        if (!vertex || vertex->GetNode() >= count) throw std::invalid_argument("Invalid primitive boundary vertex");
        const auto node = vertex->GetNode();
        const double scale = std::hypot(normalByNode[node][0], normalByNode[node][1]);
        for (unsigned d = 0; d < 2; ++d)
          RequireClose(normalByNode[node][d], -SU2_TYPE::GetValue(vertex->GetNormal(d)), scale,
                       "Primitive boundary normal inconsistent with SU2");
      }
    }
    for (unsigned long e = 0; e < faces.size(); ++e)
      if (incidence[e] == 1 && !seen[e]) throw std::invalid_argument("Missing primitive exterior boundary edge");
  }

  const Face& GetFace(unsigned long edge) const { return faces.at(edge); }
  const Cell& GetCell(unsigned long point) const { return cells.at(point); }
  std::size_t GetnEdge() const { return faces.size(); }
  std::size_t GetnPoint() const { return cells.size(); }
  const std::vector<BoundaryFace>& GetBoundaryFaces() const { return boundaryFaces; }
};
}  // namespace PrimitiveEuler2D
