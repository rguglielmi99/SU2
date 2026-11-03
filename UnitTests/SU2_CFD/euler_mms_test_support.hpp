/*!
 * \brief Frozen copy of phase-6A.1 geometry/quadrature for the shared MMS harness.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "primitive_euler_mesh_fixture.hpp"
#include "euler_common_mms.hpp"
namespace EulerMMSAudit {
using namespace PrimitiveEulerTest;
using Vec=EulerCommonMMS::Vector<double>;
constexpr double pi=3.14159265358979323846;
inline void Require(bool condition,const char* message) { if (!condition) throw std::runtime_error(message); }
struct Grid {
  unsigned n;
  std::vector<Direction> xy;
  std::vector<Triangle> triangles;
  Grid(unsigned count, const std::string& family) : n(count) {
    for (unsigned j=0;j<=n;++j) for (unsigned i=0;i<=n;++i) {
      Direction x{{double(i)/n,double(j)/n}};
      if (family=="smooth") {
        const double s=std::sin(pi*x[0])*std::sin(pi*x[1]);
        x[0]+=0.08*s; x[1]+=0.06*s;
      } else if (family=="jitter" && i && j && i<n && j<n) {
        x[0]+=0.18/n*std::sin(2.1*i+1.3*j); x[1]+=0.16/n*std::cos(1.7*i-0.9*j);
      }
      xy.push_back(x);
    }
    for (unsigned j=0;j<n;++j) for (unsigned i=0;i<n;++i) {
      const unsigned long a=j*(n+1)+i,b=a+1,d=a+n+1,c=d+1;
      triangles.push_back({{a,b,c}}); triangles.push_back({{a,c,d}});
    }
  }
  unsigned layer(unsigned long i) const {
    const unsigned x=i%(n+1),y=i/(n+1); return std::min({x,y,n-x,n-y});
  }
  bool in(unsigned long i,unsigned g) const {
    const unsigned l=layer(i);
    if (g==0) return true;
    if (g<4) return l==g-1;
    if (g==4) return l>=3;
    if (g==6) return l>=1;
    const auto& x=xy[i]; return x[0]>=0.25 && x[0]<=0.75 && x[1]>=0.25 && x[1]<=0.75;
  }
};
using Subtriangle=std::array<Direction,3>;
struct Dual {
  std::vector<Subtriangle> pieces;
  double volume=0;
  Direction offset{};
  std::array<double,3> moment{};
};
inline double Area(const Subtriangle& p) {
  return 0.5*std::abs((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]));
}
inline std::vector<std::pair<double,double>> Gauss(unsigned n) {
  std::vector<std::pair<double,double>> rule;
  for (unsigned i=0;i<n;++i) {
    double z=std::cos(pi*(i+0.75)/(n+0.5)), derivative=0;
    for (unsigned it=0;it<30;++it) {
      double a=1,b=0;
      for (unsigned k=1;k<=n;++k) { const double c=b; b=a; a=((2*k-1)*z*b-(k-1)*c)/k; }
      derivative=n*(z*a-b)/(z*z-1);
      const double next=z-a/derivative;
      if (std::abs(next-z)<2e-16) { z=next; break; } z=next;
    }
    rule.push_back({0.5*(1+z),1/((1-z*z)*derivative*derivative)});
  }
  return rule;
}
template<class Function> Vec Integrate(const Dual& cell,unsigned order,const Function& f) {
  const auto rule=Gauss(order); Vec result{};
  for (const auto& p:cell.pieces) for (const auto& r:rule) for (const auto& s:rule) {
    Direction q{};
    for (unsigned d=0;d<2;++d) q[d]=p[0][d]+r.first*(p[1][d]-p[0][d])+(1-r.first)*s.first*(p[2][d]-p[0][d]);
    const auto value=f(q);
    const double weight=2*Area(p)*(1-r.first)*r.second*s.second;
    for (unsigned k=0;k<4;++k) result[k]+=weight*value[k];
  }
  for (auto& v:result) v/=cell.volume;
  return result;
}
inline std::vector<Dual> BuildDual(const Grid& grid,const CGeometry& geometry,double& validation) {
  std::vector<Dual> cells(grid.xy.size());
  for (const auto& t:grid.triangles) {
    Direction c{}; for (auto i:t) for (unsigned d=0;d<2;++d) c[d]+=grid.xy[i][d]/3;
    const auto& a=grid.xy[t[0]]; const auto& b=grid.xy[t[1]]; const auto& d=grid.xy[t[2]];
    Require((b[0]-a[0])*(d[1]-a[1])-(b[1]-a[1])*(d[0]-a[0])>0.2/(grid.n*grid.n),
            "inverted or degenerate triangle");
    for (unsigned k=0;k<3;++k) {
      const auto i=t[k]; Direction a{},b{};
      for (unsigned d=0;d<2;++d) { a[d]=(grid.xy[i][d]+grid.xy[t[(k+1)%3]][d])/2; b[d]=(grid.xy[i][d]+grid.xy[t[(k+2)%3]][d])/2; }
      cells[i].pieces.push_back({{grid.xy[i],a,c}}); cells[i].pieces.push_back({{grid.xy[i],c,b}});
    }
  }
  for (unsigned long i=0;i<cells.size();++i) {
    auto& cell=cells[i];
    // Exact polynomial moments from barycentric integration; independent of Duffy quadrature.
    for (const auto& p:cell.pieces) {
      const double area=Area(p); cell.volume+=area; Direction sum{};
      for (auto x:p) for (unsigned d=0;d<2;++d) sum[d]+=x[d]-grid.xy[i][d];
      for (unsigned d=0;d<2;++d) cell.offset[d]+=area*sum[d]/3;
      for (unsigned m=0;m<3;++m) {
        const unsigned a=m==2?1:0,b=m==0?0:1; double diagonal=0;
        for (auto x:p) diagonal+=(x[a]-grid.xy[i][a])*(x[b]-grid.xy[i][b]);
        cell.moment[m]+=area*(sum[a]*sum[b]+diagonal)/12;
      }
    }
    const double native=SU2_TYPE::GetValue(geometry.nodes->GetVolume(i));
    Require(std::abs(cell.volume/native-1)<2e-12,"dual volume differs from SU2");
    for (auto& x:cell.offset) x/=cell.volume;
    for (auto& x:cell.moment) x/=cell.volume;
    for (unsigned order:{3u,5u,8u}) {
      const auto poly=Integrate(cell,order,[&](const Direction& x) {
        const double a=x[0]-grid.xy[i][0],b=x[1]-grid.xy[i][1]; return Vec{{1,a,b,a*a}};
      });
      const Vec exact{{1,cell.offset[0],cell.offset[1],cell.moment[0]}};
      for (unsigned k=0;k<4;++k) validation=std::max(validation,std::abs(poly[k]-exact[k]));
      const auto mixed=Integrate(cell,order,[&](const Direction& x) {
        const double a=x[0]-grid.xy[i][0],b=x[1]-grid.xy[i][1]; return Vec{{a*b,b*b,0,0}};
      });
      validation=std::max({validation,std::abs(mixed[0]-cell.moment[1]),std::abs(mixed[1]-cell.moment[2])});
    }
  }
  Require(validation<2e-13,"polynomial quadrature validation failed");
  return cells;
}
const std::array<const char*,7> groupNames{{"all","boundary","ring1","ring2","deep3","fixed","interior"}};

}  // namespace EulerMMSAudit
