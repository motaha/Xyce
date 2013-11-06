//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_ScalingVars.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/04/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_ScalingVars_h
#define Xyce_ScalingVars_h

// ---------- Standard Includes ----------
#ifdef Xyce_DEBUG_DEVICE
#include<iostream>
#endif

namespace Xyce {
namespace Device {

// ----------   Xyce Includes   ----------

// ---------- Forward Declarations ----------

//-----------------------------------------------------------------------------
// Class         : ScalingVars
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/04/08
//-----------------------------------------------------------------------------
class ScalingVars
{
  public:
    ScalingVars () :
      x0(1.0), a0(1.0), T0(1.0), V0(1.0),
      rV0(1.0), C0(1.0), D0(1.0), u0(1.0),
      R0(1.0), rR0(1.0), t0(1.0), E0(1.0),
      F0(1.0), J0(1.0), L0(1.0), k0(1.0),
      rt0(1.0), rk0(1.0)
    {};

  public:
    double x0;  // distance scaling (cm)
    double a0;  // area scaling (cm^2)
    double T0;  // temperature scaling (K)
    double V0;  // electrostatic potential scaling (V)
    double rV0; // reciprocal of V0
    double C0;  // concentration scaling (cm^-3);
    double D0;  // diffusion coefficient scaling (cm^2/s)
    double u0;  // mobility coefficient scaling (cm^2/V/s)
    double R0;  // recombination rate scaling (cm^-3/s)
    double rR0; // reciprocal of R0
    double t0;  // time scaling (s)
    double E0;  // electric field scaling (V/s)
    double F0;  // particle flux scaling (cm^-2/s)
    double J0;  // current density scaling (A/cm^2)
    double L0;  // Laplacian scaling constant

    double k0;
    double rt0;
    double rk0;

  protected:

  private:

};

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : ScalingVars::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/09
//-----------------------------------------------------------------------------
inline ostream & operator<<(ostream & os, const ScalingVars & scaleVars)
{
  os << "\n\n-----------------------------------------\n";
  os << "\tPDE Scaling Vars:\n";
  //os << "\t\tdefad                 = " << scaleVars. <<"\n";

  os << "	x0   = " << scaleVars.x0  << "\n";  // distance scaling (cm)
  os << "	a0   = " << scaleVars.a0  << "\n";  // area scaling (cm^2)
  os << "	T0   = " << scaleVars.T0  << "\n";  // temperature scaling (K)
  os << "	V0   = " << scaleVars.V0  << "\n";  // electrostatic potential scaling (V)
  os << "	rV0  = " << scaleVars.rV0 << "\n"; // reciprocal of V0
  os << "	C0   = " << scaleVars.C0  << "\n";  // concentration scaling (cm^-3)
  os << "	D0   = " << scaleVars.D0  << "\n";  // diffusion coefficient scaling (cm^2/s)
  os << "	u0   = " << scaleVars.u0  << "\n";  // mobility coefficient scaling (cm^2/V/s)
  os << "	R0   = " << scaleVars.R0  << "\n";  // recombination rate scaling (cm^-3/s)
  os << "	rR0  = " << scaleVars.rR0 << "\n"; // reciprocal of R0
  os << "	t0   = " << scaleVars.t0  << "\n";  // time scaling (s)
  os << "	E0   = " << scaleVars.E0  << "\n";  // electric field scaling (V/s)
  os << "	F0   = " << scaleVars.F0  << "\n";  // particle flux scaling (cm^-2/s)
  os << "	J0   = " << scaleVars.J0  << "\n";  // current density scaling (A/cm^2)
  os << "	L0   = " << scaleVars.L0  << "\n";  // Laplacian scaling constant

  os << "	k0   = " << scaleVars.k0  << "\n";
  os << "	rt0  = " << scaleVars.rt0 << "\n";
  os << "	rk0  = " << scaleVars.rk0 << "\n";

  os << "-----------------------------------------\n";
  os << endl;

  return os;
}

#endif

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ScalingVars N_DEV_ScalingVars;

#endif

