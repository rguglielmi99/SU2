/*!
 * \brief Analytic MMS tests independent of residual assembly.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#include "catch.hpp"
#include "euler_common_mms.hpp"
#include "../../SU2_CFD/include/numerics/flow/primitive_euler_2d.hpp"
#include <algorithm>
#include <limits>
#include <iostream>
namespace M=EulerCommonMMS;

TEST_CASE("Common Euler MMS derivatives sources and transformation are independent", "[euler-common-mms]") {
  double worstFD=0,worstIdentity=0,worstKernel=0;
  for (unsigned i=0;i<11;++i) for (double td:{0.0,0.037,0.1,0.283}) {
    const long double x=0.073L+0.077L*i,y=0.119L+0.061L*((7*i)%11),t=td;
    const auto s=M::Evaluate(x,y,t);
    REQUIRE(s.w[0]>=0.9L); REQUIRE(s.w[3]>=0.9L);
    const auto sw=M::PrimitiveSource(s),su=M::ConservativeSource(s),pde=M::ConservativePDE(s);
    for (unsigned k=0;k<4;++k) {
      worstIdentity=std::max(worstIdentity,double(std::abs(su[k]-pde[k])));
      CHECK(std::abs(su[k]-pde[k])<128*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::abs(su[k])));
    }
    // Five-point derivatives of independent sampled W, U and physical fluxes.
    for (long double h:{0.001L,0.0005L,0.00025L}) {
      const auto derivative=[&](unsigned axis,const auto& function) {
        M::Vector<long double> result{};
        for (int j:{-2,-1,1,2}) {
          const auto shifted=M::Evaluate(x+(axis==0?j*h:0),y+(axis==1?j*h:0),t+(axis==2?j*h:0));
          const auto value=function(shifted.w);
          const int weight=j==-2?1:(j==-1?-8:(j==1?8:-1));
          for (unsigned k=0;k<4;++k) result[k]+=weight*value[k]/(12*h);
        }
        return result;
      };
      for (unsigned axis=0;axis<3;++axis) {
        const auto fd=derivative(axis,[](const auto& w){return w;});
        const auto& exact=axis==0?s.dx:(axis==1?s.dy:s.dt);
        for (unsigned k=0;k<4;++k) {
          worstFD=std::max(worstFD,double(std::abs(fd[k]-exact[k])));
          CHECK(std::abs(fd[k]-exact[k])<4e-10L);
        }
      }
      const auto ut=derivative(2,[](const auto& w){return M::Conservative(w);});
      const auto fx=derivative(0,[](const auto& w){return M::Flux(w,0);});
      const auto gy=derivative(1,[](const auto& w){return M::Flux(w,1);});
      for (unsigned k=0;k<4;++k) CHECK(std::abs(ut[k]+fx[k]+gy[k]-su[k])<4e-9L);
    }
    const auto d=M::Evaluate(double(x),double(y),double(t));
    const auto ax=PrimitiveEuler2D::Apply(PrimitiveEuler2D::PhysicalMatrix(d.w,{{1,0}},M::gamma),d.dx);
    const auto ay=PrimitiveEuler2D::Apply(PrimitiveEuler2D::PhysicalMatrix(d.w,{{0,1}},M::gamma),d.dy);
    const auto source=M::PrimitiveSource(d);
    for (unsigned k=0;k<4;++k) {
      worstKernel=std::max(worstKernel,std::abs(source[k]-d.dt[k]-ax[k]-ay[k]));
      CHECK(std::abs(source[k]-d.dt[k]-ax[k]-ay[k])<2e-15);
    }
    // Check every column of dU/dW independently, not just its action on S_W.
    const auto j=M::Jacobian(s.w); const long double h=1e-3L;
    for (unsigned column=0;column<4;++column) {
      auto left=s.w,right=s.w; left[column]-=h; right[column]+=h;
      const auto a=M::Conservative(left),b=M::Conservative(right);
      for (unsigned row=0;row<4;++row) CHECK(std::abs((b[row]-a[row])/(2*h)-j[row][column])<
          128*std::numeric_limits<long double>::epsilon()/h*std::max(1.0L,std::abs(j[row][column])));
    }
    CHECK(std::abs(sw[3]-s.dt[3])>0.001L);  // Nontrivial spatial operator at these samples.
  }
  std::cout << "MMS_ANALYTIC,worst_W_FD=" << worstFD << ",worst_source_identity=" << worstIdentity
            << ",worst_primitive_kernel=" << worstKernel << '\n';
}
