/*!
 * \brief Homogeneous translating vortex and independently differentiated profile.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "euler_common_mms.hpp"
#include "../../Common/include/verification/isentropic_vortex.hpp"

namespace IsentropicVortex {
using SU2IsentropicVortex::gamma;
using SU2IsentropicVortex::beta;
using SU2IsentropicVortex::finalTime;
using SU2IsentropicVortex::length;
using SU2IsentropicVortex::center0;
using EulerCommonMMS::Conservative;
using EulerCommonMMS::Jacobian;
using EulerCommonMMS::Apply;
template<class T> EulerCommonMMS::Sample<T> Evaluate(T x,T y,T t) {
  return SU2IsentropicVortex::Evaluate<T,EulerCommonMMS::Sample<T>>(x,y,t);
}
}  // namespace IsentropicVortex
