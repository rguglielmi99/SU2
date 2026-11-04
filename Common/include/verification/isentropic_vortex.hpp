/*!
 * \brief Shared analytic benchmark, no manufactured source or discrete operator.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include <array>
#include <cmath>

namespace SU2IsentropicVortex {
constexpr double gamma=1.4, beta=2, finalTime=0.5, length=10, center0=-0.25;
template<class T> struct Sample { std::array<T,4> w,dt,dx,dy; };
// Optional sample type keeps the existing phase-7 analytic API unchanged.
template<class T, class Result=Sample<T>> Result Evaluate(T x,T y,T t) {
  const T pi=std::acos(T(-1)), X=x-T(center0)-t,Y=y-T(center0)-t;
  const T f=std::exp((1-X*X-Y*Y)/2),c=T(beta)/(2*pi);
  const T a=T(0.4L)*T(beta*beta)/(8*T(1.4L)*pi*pi);
  const T theta=1-a*f*f,rho=std::pow(theta,T(2.5L)),p=std::pow(theta,T(3.5L));
  const T tx=2*a*X*f*f,ty=2*a*Y*f*f;
  Result s;
  s.w={{rho,1-c*Y*f,1+c*X*f,p}};
  s.dx={{T(2.5L)*rho*tx/theta,c*X*Y*f,c*(1-X*X)*f,T(3.5L)*p*tx/theta}};
  s.dy={{T(2.5L)*rho*ty/theta,-c*(1-Y*Y)*f,-c*X*Y*f,T(3.5L)*p*ty/theta}};
  for (unsigned k=0;k<4;++k) s.dt[k]=-s.dx[k]-s.dy[k];
  return s;
}
}  // namespace SU2IsentropicVortex
