/*!
 * \brief Native FVM output with explicit primitive-field and RHS semantics.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "CFVMOutput.hpp"

class CPrimitiveEulerOutput final : public CFVMOutput {
 public:
  CPrimitiveEulerOutput(CConfig*, unsigned short);
  void SetHistoryOutputFields(CConfig*) override;
  void LoadHistoryData(CConfig*, CGeometry*, CSolver**) override;
  void SetVolumeOutputFields(CConfig*) override;
  void LoadVolumeData(CConfig*, CGeometry*, CSolver**, unsigned long) override;
  bool WriteVolumeOutput(CConfig*, unsigned long, bool, unsigned short) override;
  void WriteAdditionalFiles(CConfig*, CGeometry*, CSolver**) override;
};
