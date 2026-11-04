/*!
 * \brief Ideal-gas physical boundary states; no flux or interior discretization.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "primitive_euler_2d.hpp"
#include <limits>
#include <algorithm>

namespace PrimitiveEuler2D {
namespace Boundary {
inline void UnitNormal(const Direction& n) {
  if (!std::isfinite(n[0]) || !std::isfinite(n[1]) ||
      std::abs(std::hypot(n[0],n[1])-1)>64*std::numeric_limits<double>::epsilon())
    throw std::invalid_argument("Primitive boundary requires an outward unit normal");
}
inline double NormalVelocity(const State& w,const Direction& n) { return w[1]*n[0]+w[2]*n[1]; }
inline double Sound(const State& w,double gamma) {
  Validate(w,gamma);
  const double c=std::sqrt(gamma*w[3]/w[0]);
  if (!std::isfinite(c) || c<=0) throw std::invalid_argument("Primitive boundary sound speed is not representable");
  return c;
}
inline unsigned Incoming(const State& w,const Direction& n,double gamma) {
  UnitNormal(n);
  const double u=NormalVelocity(w,n), c=Sound(w,gamma);
  const double tol=64*std::numeric_limits<double>::epsilon()*std::max(c,std::abs(u));
  return (u-c < -tol ? 1u:0u)+(u < -tol ? 2u:0u)+(u+c < -tol ? 1u:0u);
}
inline State Reflect(const State& w,const Direction& n,double gamma) {
  UnitNormal(n); Validate(w,gamma);
  auto result=w; const double u=NormalVelocity(w,n);
  result[1]-=2*u*n[0]; result[2]-=2*u*n[1]; return result;
}
// Locally characteristic reservoir: project the reservoir jump onto incoming
// eigenvectors of A_n(W_int). This is explicitly a frozen-state nonlinear BC,
// not the finite-amplitude Riemann-invariant construction used by native FAR.
// In particular a tangential entropy variation must not launch acoustic waves.
inline State Farfield(const State& w,const State& reservoir,const Direction& n,double gamma) {
  UnitNormal(n);
  const double u=NormalVelocity(w,n), c=Sound(w,gamma);
  Validate(reservoir,gamma);
  const unsigned incoming=Incoming(w,n,gamma);
  if (incoming==0) return w;
  if (incoming==4) {
    if (Incoming(reservoir,n,gamma)!=4)
      throw std::invalid_argument("Primitive supersonic farfield reservoir has incompatible incoming regime");
    return reservoir;
  }
  const double du=Boundary::NormalVelocity(reservoir,n)-u, dp=reservoir[3]-w[3];
  const double acousticMinus=(dp-w[0]*c*du)/(2*c*c);
  State result=w;
  result[0]+=acousticMinus;
  result[1]-=acousticMinus*c*n[0]/w[0];
  result[2]-=acousticMinus*c*n[1]/w[0];
  result[3]+=acousticMinus*c*c;
  const double tol=64*std::numeric_limits<double>::epsilon()*std::max(c,std::abs(u));
  if (u < -tol) {
    result[0]+=reservoir[0]-w[0]-dp/(c*c);
    const double tangent=-(reservoir[1]-w[1])*n[1]+(reservoir[2]-w[2])*n[0];
    result[1]-=tangent*n[1]; result[2]+=tangent*n[0];
  }
  Validate(result,gamma);
  const unsigned resulting=Incoming(result,n,gamma);
  if (resulting==0 || resulting==4)
    throw std::invalid_argument("Primitive locally characteristic farfield crosses a sonic regime; finite-amplitude sonic boundary solve is not implemented");
  // Classification is frozen at W_int for this stage and refreshed next stage.
  // At tangency a small acoustic perturbation can legitimately change its sign.
  return result;
}
inline State Outlet(const State& w,double pressure,const Direction& n,double gamma) {
  UnitNormal(n); const double u=NormalVelocity(w,n), c=Sound(w,gamma);
  if (!std::isfinite(pressure) || pressure<=0)
    throw std::invalid_argument("Primitive outlet requires positive static pressure");
  if (u>=c) return w;
  if (u<=0) throw std::invalid_argument("Primitive pressure outlet backflow/tangency needs additional incoming data");
  const double rho=w[0]*std::pow(pressure/w[3],1/gamma);
  const double cb=std::sqrt(gamma*pressure/rho), ub=u+2*(c-cb)/(gamma-1);
  State result{{rho,w[1]+(ub-u)*n[0],w[2]+(ub-u)*n[1],pressure}};
  Validate(result,gamma);
  if (ub<=0 || ub>=cb)
    throw std::invalid_argument("Primitive pressure outlet prescribed pressure changes subsonic outflow regime");
  return result;
}
inline State SupersonicInlet(const State& w,const State& prescribed,const Direction& n,double gamma) {
  if (Incoming(w,n,gamma)!=4 || Incoming(prescribed,n,gamma)!=4)
    throw std::invalid_argument("Primitive supersonic inlet requires four incoming characteristics");
  return prescribed;
}
inline State SupersonicOutlet(const State& w,const Direction& n,double gamma) {
  if (Incoming(w,n,gamma)!=0)
    throw std::invalid_argument("Primitive supersonic outlet requires no incoming characteristics");
  return w;
}
inline Direction InletDirection(Direction dir,const Direction& n) {
  UnitNormal(n); const double mag=std::hypot(dir[0],dir[1]);
  if (!std::isfinite(mag) || mag<=0) throw std::invalid_argument("Primitive inlet has invalid flow direction");
  for (auto& x:dir) x/=mag;
  if (dir[0]*n[0]+dir[1]*n[1]>=0)
    throw std::invalid_argument("Primitive inlet direction must point into the domain");
  return dir;
}
inline void SubsonicInflow(const State& w,const Direction& n,double gamma) {
  if (Incoming(w,n,gamma)!=3)
    throw std::invalid_argument("Primitive MARKER_INLET requires subsonic inflow; use MARKER_SUPERSONIC_INLET for full static data");
}
inline State MassInlet(const State& w,double rho,double speed,Direction dir,const Direction& n,double gamma) {
  SubsonicInflow(w,n,gamma); dir=InletDirection(dir,n);
  if (!std::isfinite(rho) || rho<=0 || !std::isfinite(speed) || speed<=0)
    throw std::invalid_argument("Primitive MASS_FLOW inlet requires positive density and velocity magnitude");
  const double plus=NormalVelocity(w,n)+2*Sound(w,gamma)/(gamma-1);
  const double c=.5*(gamma-1)*(plus-speed*(dir[0]*n[0]+dir[1]*n[1]));
  if (c<=0) throw std::invalid_argument("Primitive MASS_FLOW inlet has incompatible outgoing invariant");
  State result{{rho,speed*dir[0],speed*dir[1],rho*c*c/gamma}};
  Validate(result,gamma); SubsonicInflow(result,n,gamma); return result;
}
inline State TotalInlet(const State& w,double tt,double pt,Direction dir,const Direction& n,double gamma,double gasR) {
  SubsonicInflow(w,n,gamma); dir=InletDirection(dir,n);
  if (!std::isfinite(tt) || tt<=0 || !std::isfinite(pt) || pt<=0 || !std::isfinite(gasR) || gasR<=0)
    throw std::invalid_argument("Primitive total inlet requires positive Tt, Pt and gas constant");
  const double alpha=dir[0]*n[0]+dir[1]*n[1], gm=gamma-1;
  const double plus=NormalVelocity(w,n)+2*Sound(w,gamma)/gm;
  const double ctotal2=gamma*gasR*tt;
  // (plus-alpha*q)^2 = 4/gm^2 (ctotal^2-gm*q^2/2).
  const double a=1+.5*gm*alpha*alpha,b=-gm*alpha*plus;
  const double d=.5*gm*plus*plus-2*ctotal2/gm;
  const double discriminant=b*b-4*a*d;
  if (!std::isfinite(discriminant) || discriminant<0)
    throw std::invalid_argument("Primitive total inlet has no real characteristic-compatible state");
  const double speed=(-b+std::sqrt(discriminant))/(2*a);
  const double c2=ctotal2-.5*gm*speed*speed;
  if (speed<=0 || c2<=0 || plus-alpha*speed<=0)
    throw std::invalid_argument("Primitive total inlet has no positive characteristic-compatible state");
  const double temperature=c2/(gamma*gasR), pressure=pt*std::pow(temperature/tt,gamma/gm);
  State result{{pressure/(gasR*temperature),speed*dir[0],speed*dir[1],pressure}};
  Validate(result,gamma); SubsonicInflow(result,n,gamma); return result;
}
}  // namespace Boundary
}  // namespace PrimitiveEuler2D
