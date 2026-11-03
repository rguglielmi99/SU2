/*!
 * \brief Independent analytic checks for the unforced isentropic vortex.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "euler_isentropic_vortex.hpp"
#include "../../SU2_CFD/include/numerics/flow/primitive_euler_2d.hpp"
#include <iostream>
#include <stdexcept>
#include <algorithm>

int main() {
  namespace V=IsentropicVortex;
  namespace M=EulerCommonMMS;
  double worstP=0,worstU=0,worstFD=0,worstFluxFD=0,worstEntropy=0,worstConversion=0;
  unsigned checks=0;
  const auto require=[&](bool ok) { ++checks; if (!ok) throw std::runtime_error("vortex analytic check failed"); };
  for (double t:{0.,.25,.5}) for (double x:{-5.,-.91,-.25,.17,.73,1.6,5.})
    for (double y:{-5.,-1.1,-.25,.13,.94,2.,5.}) {
      const auto s=V::Evaluate(x,y,t);
      const auto pde=M::PrimitiveSource(s),cpde=M::ConservativePDE(s),u=V::Conservative(s.w);
      const auto ax=PrimitiveEuler2D::Apply(PrimitiveEuler2D::PhysicalMatrix(s.w,{{1,0}},V::gamma),s.dx);
      const auto ay=PrimitiveEuler2D::Apply(PrimitiveEuler2D::PhysicalMatrix(s.w,{{0,1}},V::gamma),s.dy);
      for (unsigned k=0;k<4;++k) {
        worstP=std::max({worstP,std::abs(pde[k]),std::abs(s.dt[k]+ax[k]+ay[k])});
        worstU=std::max(worstU,std::abs(cpde[k]));
      }
      require(s.w[0]>0 && s.w[3]>0);
      require(s.w[3]/s.w[0]>0);  // Dimensionless temperature is p/rho.
      worstEntropy=std::max(worstEntropy,std::abs(s.w[3]/std::pow(s.w[0],V::gamma)-1));
      const double pressure=.4*(u[3]-(u[1]*u[1]+u[2]*u[2])/(2*u[0]));
      worstConversion=std::max(worstConversion,std::abs(pressure-s.w[3]));
      worstConversion=std::max({worstConversion,std::abs(u[1]/u[0]-s.w[1]),std::abs(u[2]/u[0]-s.w[2])});
      const double eps=1e-5;
      std::array<M::Sample<double>,6> q{{V::Evaluate(x+eps,y,t),V::Evaluate(x-eps,y,t),
          V::Evaluate(x,y+eps,t),V::Evaluate(x,y-eps,t),V::Evaluate(x,y,t+eps),V::Evaluate(x,y,t-eps)}};
      const auto fxp=M::Flux(q[0].w,0),fxm=M::Flux(q[1].w,0),fyp=M::Flux(q[2].w,1),fym=M::Flux(q[3].w,1);
      const auto utp=V::Conservative(q[4].w),utm=V::Conservative(q[5].w);
      for (unsigned k=0;k<4;++k) {
        worstFD=std::max({worstFD,std::abs((q[0].w[k]-q[1].w[k])/(2*eps)-s.dx[k]),
            std::abs((q[2].w[k]-q[3].w[k])/(2*eps)-s.dy[k]),
            std::abs((q[4].w[k]-q[5].w[k])/(2*eps)-s.dt[k])});
        worstFluxFD=std::max(worstFluxFD,std::abs((utp[k]-utm[k]+fxp[k]-fxm[k]+fyp[k]-fym[k])/(2*eps)));
      }
      const auto translated=V::Evaluate(x+t,y+t,t).w,initial=V::Evaluate(x,y,0.).w;
      for (unsigned k=0;k<4;++k) require(std::abs(translated[k]-initial[k])<1e-14);
    }
  for (double t:{0.,.25,.5}) {
    const auto center=V::Evaluate(V::center0+t,V::center0+t,t).w;
    require(center[1]==1 && center[2]==1);
    require(center[3]<V::Evaluate(V::center0+t+.1,V::center0+t,t).w[3]);
  }
  require(worstP<2e-14 && worstU<3e-14 && worstFD<2e-8 && worstFluxFD<2e-8);
  require(worstEntropy<2e-14 && worstConversion<2e-14);
  std::cout << "VORTEX_ANALYTIC,checks=" << checks << ",primitive=" << worstP << ",conservative=" << worstU
            << ",derivative_FD=" << worstFD << ",flux_FD=" << worstFluxFD << ",entropy=" << worstEntropy
            << ",conversion=" << worstConversion << '\n';
}
