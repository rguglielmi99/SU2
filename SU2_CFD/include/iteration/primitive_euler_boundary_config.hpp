/*!
 * \brief Native SU2 marker/data adapter for primitive exterior states.
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md).
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */
#pragma once
#include "../../../Common/include/CConfig.hpp"
#include "../numerics/flow/primitive_euler_boundary.hpp"
#include "../../../Common/include/verification/isentropic_vortex.hpp"

namespace PrimitiveEulerBoundaryConfig {
using namespace PrimitiveEuler2D;
inline bool Supported(unsigned short kind) {
  return kind==CUSTOM_BOUNDARY || kind==FAR_FIELD || kind==EULER_WALL || kind==SYMMETRY_PLANE ||
      kind==INLET_FLOW || kind==OUTLET_FLOW || kind==SUPERSONIC_INLET || kind==SUPERSONIC_OUTLET;
}
inline State Freestream(const CConfig& c) {
  if (c.GetSystemMeasurements()!=SI)
    throw std::invalid_argument("Primitive physical boundary adapter requires SI units");
  const double p=SU2_TYPE::GetValue(c.GetPressure_FreeStream());
  const double rho=c.GetKind_FreeStreamOption()==FREESTREAM_OPTION::DENSITY_FS ?
      SU2_TYPE::GetValue(c.GetDensity_FreeStream()) :
      p/(SU2_TYPE::GetValue(c.GetGas_Constant())*SU2_TYPE::GetValue(c.GetTemperature_FreeStream()));
  const double gamma=SU2_TYPE::GetValue(c.GetGamma()), angle=SU2_TYPE::GetValue(c.GetAoA())*std::acos(-1.)/180;
  const double speed=SU2_TYPE::GetValue(c.GetMach())*std::sqrt(gamma*p/rho);
  State result{{rho,speed*std::cos(angle),speed*std::sin(angle),p}};
  Validate(result,gamma); return result;
}
inline State Evaluate(const CConfig& c,unsigned short marker,const Direction& q,double time,
                      const State& inside,const Direction& normal) {
  if (c.GetSystemMeasurements()!=SI)
    throw std::invalid_argument("Primitive physical boundary adapter requires SI units");
  const auto kind=c.GetMarker_All_KindBC(marker);
  const auto tag=c.GetMarker_All_TagBound(marker);
  const double gamma=SU2_TYPE::GetValue(c.GetGamma()), gasR=SU2_TYPE::GetValue(c.GetGas_Constant());
  switch (kind) {
    case CUSTOM_BOUNDARY: {
      if (c.GetPrimitiveBoundaryMode()=="ISENTROPIC_VORTEX")
        return SU2IsentropicVortex::Evaluate(q[0],q[1],time).w;
      State result{};
      for (unsigned k=0;k<4;++k) result[k]=SU2_TYPE::GetValue(c.GetPrimitiveBoundaryState()[k]);
      return result;
    }
    case FAR_FIELD: return Boundary::Farfield(inside,Freestream(c),normal,gamma);
    case EULER_WALL: case SYMMETRY_PLANE: return Boundary::Reflect(inside,normal,gamma);
    case OUTLET_FLOW: return Boundary::Outlet(inside,SU2_TYPE::GetValue(c.GetOutlet_Pressure(tag)),normal,gamma);
    case SUPERSONIC_OUTLET: return Boundary::SupersonicOutlet(inside,normal,gamma);
    case SUPERSONIC_INLET: {
      const double temperature=SU2_TYPE::GetValue(c.GetInlet_Temperature(tag));
      const double pressure=SU2_TYPE::GetValue(c.GetInlet_Pressure(tag));
      const auto* velocity=c.GetInlet_Velocity(tag);
      if (velocity[2]!=0) throw std::invalid_argument("Primitive 2D inlet requires zero z velocity");
      State prescribed{{pressure/(gasR*temperature),SU2_TYPE::GetValue(velocity[0]),SU2_TYPE::GetValue(velocity[1]),pressure}};
      return Boundary::SupersonicInlet(inside,prescribed,normal,gamma);
    }
    case INLET_FLOW: {
      const auto* d=c.GetInletFlowDir(tag);
      if (d[2]!=0) throw std::invalid_argument("Primitive 2D inlet requires zero z direction");
      Direction dir{{SU2_TYPE::GetValue(d[0]),SU2_TYPE::GetValue(d[1])}};
      if (c.GetInletUseNormal()) dir={{-normal[0],-normal[1]}};
      if (c.GetKind_Inlet()==INLET_TYPE::TOTAL_CONDITIONS)
        return Boundary::TotalInlet(inside,SU2_TYPE::GetValue(c.GetInletTtotal(tag)),
            SU2_TYPE::GetValue(c.GetInletPtotal(tag)),dir,normal,gamma,gasR);
      if (c.GetKind_Inlet()==INLET_TYPE::MASS_FLOW)
        return Boundary::MassInlet(inside,SU2_TYPE::GetValue(c.GetInletTtotal(tag)),
            SU2_TYPE::GetValue(c.GetInletPtotal(tag)),dir,normal,gamma);
      throw std::invalid_argument("Primitive unsupported INLET_TYPE");
    }
    default: throw std::invalid_argument("Primitive unsupported physical boundary marker "+tag);
  }
}
}  // namespace PrimitiveEulerBoundaryConfig
