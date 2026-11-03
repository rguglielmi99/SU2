/*!
 * \brief Shared analytic space-time MMS; independent of either discrete operator.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include <array>
#include <cmath>

namespace EulerCommonMMS {
constexpr double gamma = 1.4;
template<class T> using Vector = std::array<T,4>;
template<class T> using Matrix = std::array<Vector<T>,4>;
template<class T> struct Sample { Vector<T> w, dt, dx, dy; };

template<class T> Sample<T> Evaluate(T x,T y,T t) {
  const T k=2*std::acos(T(-1)), sx=std::sin(k*x),cx=std::cos(k*x),sy=std::sin(k*y),cy=std::cos(k*y);
  const T a=T(0.1L)*std::cos(t),b=T(0.1L)*std::sin(t),at=-b,bt=a;
  return {{{1+a*sx*sy,T(0.5L)+a*cx*sy,T(0.3L)+a*sx*cy,1+b*cx*cy}},
          {{at*sx*sy,at*cx*sy,at*sx*cy,bt*cx*cy}},
          {{k*a*cx*sy,-k*a*sx*sy,k*a*cx*cy,-k*b*sx*cy}},
          {{k*a*sx*cy,k*a*cx*cy,-k*a*sx*sy,-k*b*cx*sy}}};
}
template<class T> Vector<T> Conservative(const Vector<T>& w) {
  return {{w[0],w[0]*w[1],w[0]*w[2],w[3]/T(0.4L)+w[0]*(w[1]*w[1]+w[2]*w[2])/2}};
}
template<class T> Matrix<T> Jacobian(const Vector<T>& w) {
  return {{{{1,0,0,0}},{{w[1],w[0],0,0}},{{w[2],0,w[0],0}},
           {{(w[1]*w[1]+w[2]*w[2])/2,w[0]*w[1],w[0]*w[2],1/T(0.4L)}}}};
}
template<class T> Vector<T> Apply(const Matrix<T>& j,const Vector<T>& v) {
  Vector<T> result{};
  for (unsigned k=0;k<4;++k) for (unsigned m=0;m<4;++m) result[k]+=j[k][m]*v[m];
  return result;
}
template<class T> Vector<T> PrimitiveSource(const Sample<T>& s) {
  const auto& w=s.w; const auto& x=s.dx; const auto& y=s.dy;
  return {{s.dt[0]+w[1]*x[0]+w[2]*y[0]+w[0]*(x[1]+y[2]),
           s.dt[1]+w[1]*x[1]+w[2]*y[1]+x[3]/w[0],
           s.dt[2]+w[1]*x[2]+w[2]*y[2]+y[3]/w[0],
           s.dt[3]+w[1]*x[3]+w[2]*y[3]+T(1.4L)*w[3]*(x[1]+y[2])}};
}
template<class T> Vector<T> ConservativeSource(const Sample<T>& s) {
  return Apply(Jacobian(s.w),PrimitiveSource(s));
}
template<class T> Vector<T> Flux(const Vector<T>& w,unsigned direction) {
  const auto u=Conservative(w); const T speed=w[1+direction];
  return {{u[0]*speed,u[1]*speed+(direction==0?w[3]:T(0)),
           u[2]*speed+(direction==1?w[3]:T(0)),(u[3]+w[3])*speed}};
}
// Independent product differentiation of U_t + F_x + G_y: no J or source reuse.
template<class T> Vector<T> ConservativePDE(const Sample<T>& s) {
  const auto& w=s.w; const auto& x=s.dx; const auto& y=s.dy; const auto& t=s.dt;
  const T r=w[0],u=w[1],v=w[2],p=w[3],q=(u*u+v*v)/2,E=p/T(0.4L)+r*q;
  const auto energyDerivative=[&](const Vector<T>& d) { return d[3]/T(0.4L)+d[0]*q+r*(u*d[1]+v*d[2]); };
  return {{t[0]+x[0]*u+r*x[1]+y[0]*v+r*y[2],
           t[0]*u+r*t[1]+x[0]*u*u+2*r*u*x[1]+x[3]+y[0]*u*v+r*y[1]*v+r*u*y[2],
           t[0]*v+r*t[2]+x[0]*u*v+r*x[1]*v+r*u*x[2]+y[0]*v*v+2*r*v*y[2]+y[3],
           energyDerivative(t)+(energyDerivative(x)+x[3])*u+(E+p)*x[1]
                              +(energyDerivative(y)+y[3])*v+(E+p)*y[2]}};
}
}  // namespace EulerCommonMMS
