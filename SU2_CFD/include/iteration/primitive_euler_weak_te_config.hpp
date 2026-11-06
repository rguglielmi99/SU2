/*!
 * \brief Explicit serial global-to-local weak-TE selection; geometry is validated by the solver.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "../../../Common/include/geometry/CGeometry.hpp"

namespace PrimitiveEulerWeakTEConfig {
inline std::vector<unsigned long> LocalNodes(const CGeometry& geometry, const CConfig& config) {
  if (config.GetPrimitiveWeakTEGlobalID()<0) return {};
  if (config.GetKind_Solver()!=MAIN_SOLVER::PRIMITIVE_EULER || SU2_MPI::GetSize()!=1 || omp_get_max_threads()!=1)
    throw std::invalid_argument("Primitive weak TE configuration supports only serial PRIMITIVE_EULER");
  const auto global=static_cast<unsigned long>(config.GetPrimitiveWeakTEGlobalID());
  std::vector<unsigned long> local;
  for (unsigned long i=0;i<geometry.GetnPoint();++i)
    if (geometry.nodes->GetGlobalIndex(i)==global) {
      if (!geometry.nodes->GetDomain(i)) throw std::invalid_argument("Weak TE global ID is not a domain node");
      local.push_back(i);
    }
  if (local.size()!=1) throw std::invalid_argument("Weak TE global ID must identify exactly one mesh node");
  return local;
}
}
