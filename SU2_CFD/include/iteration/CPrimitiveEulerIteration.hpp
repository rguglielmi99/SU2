/*!
 * \brief Standard-driver adapter to the existing transactional primitive FE API.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "CIteration.hpp"

class CPrimitiveEulerIteration final : public CIteration {
  bool initialized=false;
 public:
  explicit CPrimitiveEulerIteration(const CConfig* config) : CIteration(config) {}
  static void ValidateDriver(const CConfig&);
  void Solve(COutput*, CIntegration****, CGeometry****, CSolver*****, CNumerics******, CConfig**,
             CSurfaceMovement**, CVolumetricMovement***, CFreeFormDefBox***, unsigned short, unsigned short) override;
};
