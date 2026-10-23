/*!
 * \brief Shared real-SU2 mesh fixture for primitive Euler tests.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../SU2_CFD/include/numerics/flow/primitive_euler_2d.hpp"
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <map>

namespace PrimitiveEulerTest {
using PrimitiveEuler2D::Direction;
using Triangle = std::array<unsigned long, 3>;
struct BoundaryMarker {
  std::string tag;
  std::vector<std::array<unsigned long, 2>> edges;
};
struct QuietGeometryOutput {
  std::ostringstream sink;
  std::streambuf* previous = std::cout.rdbuf(sink.rdbuf());
  ~QuietGeometryOutput() { std::cout.rdbuf(previous); }
};

// Use a private temporary SU2 mesh and the production mesh reader/preprocessing.
// CPhysicalGeometry is final; tests therefore exercise its public API directly.
struct MeshFixture {
  std::filesystem::path directory;
  std::unique_ptr<CConfig> config;
  std::unique_ptr<CPhysicalGeometry> geometry;
  MeshFixture(const std::vector<Direction>& xy, const std::vector<Triangle>& triangles,
              const std::string& solver = "EULER", const std::string& options = "",
              const std::string& markerOption = "MARKER_CUSTOM",
              std::vector<BoundaryMarker> markers = {}, bool reverseBoundary = false) {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    for (unsigned attempt = 0;; ++attempt) {
      directory = std::filesystem::temp_directory_path() /
                  ("su2-primitive-geometry-" + std::to_string(stamp) + "-" + std::to_string(attempt));
      if (std::filesystem::create_directory(directory)) break;
    }
    const auto filename = directory / "mesh.su2";
    std::map<std::pair<unsigned long, unsigned long>, std::vector<std::pair<unsigned long, unsigned long>>> edges;
    for (const auto& t : triangles)
      for (unsigned k = 0; k < 3; ++k) {
        const auto i = t[k], j = t[(k + 1) % 3];
        edges[std::minmax(i, j)].push_back({i, j});
      }
    std::ofstream mesh(filename);
    mesh << std::setprecision(17) << "NDIME= 2\nNELEM= " << triangles.size() << "\n";
    for (const auto& t : triangles) mesh << "5 " << t[0] << " " << t[1] << " " << t[2] << "\n";
    mesh << "NPOIN= " << xy.size() << "\n";
    for (unsigned long i = 0; i < xy.size(); ++i) mesh << xy[i][0] << " " << xy[i][1] << " " << i << "\n";
    if (markers.empty()) {
      markers.push_back({"boundary", {}});
      for (const auto& edge : edges)
        if (edge.second.size() == 1) markers[0].edges.push_back({{edge.second[0].first, edge.second[0].second}});
    }
    mesh << "NMARK= " << markers.size() << "\n";
    for (const auto& marker : markers) {
      mesh << "MARKER_TAG= " << marker.tag << "\nMARKER_ELEMS= " << marker.edges.size() << "\n";
      for (const auto& edge : marker.edges)
        mesh << "3 " << edge[reverseBoundary ? 1 : 0] << " " << edge[reverseBoundary ? 0 : 1] << "\n";
    }
    mesh.close();
    QuietGeometryOutput quiet;
    std::stringstream text;
    text << "SOLVER= " << solver << "\nMGLEVEL=0\n" << markerOption << "=(";
    for (unsigned m = 0; m < markers.size(); ++m) text << (m ? "," : "") << markers[m].tag;
    text << ")\nMESH_FILENAME= " << filename.string() << "\n" << options;
    config.reset(new CConfig(text, SU2_COMPONENT::SU2_CFD, false));
    {
      CPhysicalGeometry input(config.get(), 0, 1);
      geometry.reset(new CPhysicalGeometry(&input, config.get()));
    }
    geometry->SetSendReceive(config.get());
    geometry->SetBoundaries(config.get());
    geometry->SetPoint_Connectivity();
    geometry->SetElement_Connectivity();
    geometry->SetBoundVolume();
    geometry->Check_IntElem_Orientation(config.get());
    geometry->Check_BoundElem_Orientation(config.get());
    geometry->SetEdges();
    geometry->SetVertex(config.get());
    geometry->SetControlVolume(config.get(), ALLOCATE);
    geometry->SetBoundControlVolume(config.get(), ALLOCATE);
  }
  ~MeshFixture() {
    std::error_code ignored;
    std::filesystem::remove(directory / "mesh.su2", ignored);
    std::filesystem::remove(directory, ignored);
  }
};

}  // namespace PrimitiveEulerTest
