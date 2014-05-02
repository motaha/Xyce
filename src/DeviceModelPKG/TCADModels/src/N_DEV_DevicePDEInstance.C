//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DevicePDEInstance.C,v $
//
// Purpose        : This file contains functions that are common to both
//                  the 1D and 2D PDE devices.  Theoretically, they could
//                  be common to 3D devices as well, but those don't exist
//                  in this code.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/23/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.82 $
//
// Revision Date  : $Date: 2014/02/24 23:49:18 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DevicePDEInstance.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_SourceData.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

// Mathematical functions and derivatives.
// Most of these were lifted from the SG Framework (sgnam.cpp).

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::DevicePDEInstance
// Purpose       : constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/01/03
//-----------------------------------------------------------------------------
DevicePDEInstance::DevicePDEInstance(
  const InstanceBlock &               IB,
  ParametricData<void> &        parametric_data,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, parametric_data, factory_block),
    Temp(getDeviceOptions().temp.getImmutableValue<double>()),
    charge(1.602176565e-19),
    kb (1.3806488e-23), // boltzmann's constant
    Vt(kb*Temp/charge),
    Ut(Vt),
    e0(8.854187817e-14),
    eSi(11.8),
    eSiO2(3.9),
    eps(eSi*e0),
    Ni(1.25e10),
    h_planck (6.62606957e-34), // Planck's constant  (in J-s)
    e_mass (9.10938291e-31),  // e- mass in kg.
    // photogen section
    photogenOnFlag(false),
    xstart(0.0),
    ystart(0.0),
    xend (0.0),
    yend (0.0),
    intensity (0.0),
    photoA1 (0.0),
    photoTstart (0.0),
    photoTstop (1.0e+100),  // just something really big...
    photoTd(0.0),
    photoTr(0.0),
    photoTf(0.0),
    photoPw(0.0),
    photoPer(0.0),
    lastPeriodIndex(-1),
    photoType(_DC_DATA),
    Data_ptr(NULL),
    DataSaved_ptr(NULL),
    photoA1_old(0.0),
    photoA1_final(0.0),
    photoA1_orig(0.0),
    photoA1_ramp(0.0),
    photoA1_ramp_old(0.0),
    photoA1_Delta(0.0),
    photoA1_DeltaC(0.0),
    maxPhotoDelta(1.0e+21),
    photoContinuationFinished (false),
    // end of photogen section

    maxVoltDelta(3.0*Vt), // 3 * thermal voltage.
    enableContinuationCalled(false),
    continuationAlpha (1.0),
    mobModelName("carr"),
    fieldDependentMobility(false),
    fieldDependentMobilityGiven(false),
    bulkMaterial("si"),
    sensOn        (false),
    sensProcess   (false),
    meshSensMod   (false),
    dopingSensMod (false),
    photogenSensMod (false),
    variablesScaled(false),
    x0_user(0.0),
    C0_user(0.0),
    t0_user(0.0),
    outputName("")
{
  setupOutputName ();
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::aux1
// Purpose       : This is the first of two "auxilliary" functions, neccessary
//                 for the Scharfetter-Gummel approximation.
//
// Special Notes : To avoid under and over-flows this function is defined
//                 by equivalent or approximate functions depending upon the
//                 value of the argument.
//
//                                       x
//                         Aux1(x) =  -------
//                                    sinh(x)
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
double DevicePDEInstance::aux1(double x)
{
  if      (x < -bernSupport.bp0_MISC) x = -bernSupport.bp0_MISC;
  else if (x >  bernSupport.bp0_MISC) x =  bernSupport.bp0_MISC;

  if (x <= bernSupport.bp0_AUX1)
  {
    return(x / sinh(x));
  }
  else if (x <= bernSupport.bp1_AUX1)
  {
    return(1 - x*x/6.0*(1.0 - 7.0*x*x/60.0));
  }
  else
  {
    return(x / sinh(x));
  }
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::pd1aux1
// Purpose       : This function returns the total derivative of the aux1
//                 function.
//
// Special Notes :
// To avoid under and overflows this function is defined by equivalent or
// approximate functions depending upon the value of the argument.
//
//                    d           sinh(x) - x*cosh(x)
//                    --Aux1(x) = -------------------
//                    dx              (sinh(x))^2
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
double DevicePDEInstance::pd1aux1(double x)
{
  double y;

#if 0
  if      (x < -bernSupport.bp0_MISC) x = -bernSupport.bp0_MISC;
  else if (x >  bernSupport.bp0_MISC) x =  bernSupport.bp0_MISC;
#else
  if      (x < -700.0) x = -700.0;
  else if (x >  700.0) x =  700.0;
#endif

  if (x <= bernSupport.bp0_DAUX1)
  {
    y = sinh(x);
    return((y - x*cosh(x))/(y*y));
  }
  else if (x <= bernSupport.bp1_DAUX1)
  {
    return(-x/3.0*(1.0 - 7.0*x*x/30.0));
  }
  else
  {
    y = sinh(x);
    return((y - x*cosh(x))/(y*y));
  }

} // pd1aux1


//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::aux2
// Purpose       : This is the second of two "auxilliary" functions, neccesary
//                 for the Scharfetter-Gummel approximation.
//
// Special Notes : To avoid under and over-flows this function is defined
//                 by equivalent or approximate functions depending upon
//                 the value of the argument.
//
//                                           1
//                              Aux2(x) = -------
//                                        1 + e^x
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
double DevicePDEInstance::aux2(double x)
{
  if (x <= bernSupport.bp0_AUX2)
  {
    return(1.0);
  }
  else if (x <= bernSupport.bp1_AUX2)
  {
    return(1.0 / (1.0 + exp(x)));
  }
  else if (x <= bernSupport.bp2_AUX2)
  {
    return(exp(-x));
  }
  else
  {
    return(0.0);
  }
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::pd1aux2
// Purpose       : This function returns the total derivative of the aux2
//                 function.
//
// Special Notes : To avoid under and overflows this function is defined
//                 by equivalent or approximate functions depending upon
//                 the value of the argument.
//
//                         d             - e^x
//                         --Aux2(x) = -----------
//                         dx          (1 + e^x)^2
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
double DevicePDEInstance::pd1aux2(double x)
{
  double y,z;

  if (x <= bernSupport.bp0_DAUX2)
  {
    return(0.0);
  }
  else if (x <= bernSupport.bp1_DAUX2)
  {
    return(-exp(x));
  }
  else if (x <= bernSupport.bp2_DAUX2)
  {
    y = exp(x); z = y + 1.0; return(-y/(z*z));
  }
  else if (x <= bernSupport.bp3_DAUX2)
  {
    return(-exp(-x));
  }
  else
  {
    return(0.0);
  }

} // pd1aux2

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::Jn
// Purpose       : Electron current density between two points.
// Special Notes : Version for passing mobility along edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::Jn
(double n1, double n2, double E, double u, double h)
{
  double MU     = u;
  double dV     = E*h/(2.0*Ut);
  double n      = n1*aux2(dV)+n2*aux2(-dV);
  double dndx   = aux1(-dV)*(n2-n1)/h;
  double J      =  MU*((n*E)+(Ut*dndx));

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout().width(24);
    Xyce::dout().precision(16);
    Xyce::dout().setf(std::ios::scientific);

    Xyce::dout() << std::endl;
    Xyce::dout() << "  MU = "<< MU << std::endl;
    Xyce::dout() << "  n  = "<< n  << " dndx = "<< dndx <<" E = "<< E << std::endl;
    Xyce::dout() << "  dV = "<< dV << " Ut = "<< Ut <<" Jn = "<< J << std::endl;
    Xyce::dout() << "  n1 = "<< n1 << " n2 = "<< n2 <<" (n2-n1) = "<<(n2-n1)<< std::endl;
    Xyce::dout() << "  h  = "<< h  << std::endl;
    Xyce::dout() << "  aux1(-dV) = " << aux1(-dV) << std::endl;
    Xyce::dout() << "  aux2( dV) = " << aux2( dV) << std::endl;
    Xyce::dout() << "  aux2(-dV) = " << aux2( dV) << std::endl;
    Xyce::dout() << std::endl;
  }
#endif

  return J;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJndV1
//
// Purpose       : This function returns the derivative of electron current
//                 density between two points in space, with
//                 respect to the voltage at node 1.
//
// Special Notes : This version passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJndV1
(double n1, double n2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double n        = n1*aux2(dV)+n2*aux2(-dV);
  double dEdv1    =  1.0/h;
  double ddVdv1   =  1.0/(2.0*Ut);
  double dNdv1    = n1*( ddVdv1*pd1aux2( dV) ) + n2*(-ddVdv1*pd1aux2(-dV) );
  double dDNDXdv1 = (n2-n1)/h * (-ddVdv1) * pd1aux1(-dV);
  double dJdv1    = MU*(dNdv1*E + n*dEdv1 + Ut*dDNDXdv1);
  return dJdv1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJndV2
//
// Purpose       : This function returns the derivative of electron current
//                 density between two points in space, with
//                 respect to the voltage at node 2.
//
// Special Notes : This version passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJndV2
(double n1, double n2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double n        = n1*aux2(dV)+n2*aux2(-dV);
  double dEdv2    = -1.0/h;
  double ddVdv2   = -1.0/(2.0*Ut);
  double dNdv2    = n1*( ddVdv2*pd1aux2( dV) ) + n2*(-ddVdv2*pd1aux2(-dV) );
  double dDNDXdv2 = (n2-n1)/h * (-ddVdv2) * pd1aux1(-dV);
  double dJdv2    = MU*(dNdv2*E + n*dEdv2 + Ut*dDNDXdv2);
  return dJdv2;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJndn1
//
// Purpose       : This function returns the derivative of electron current
//                 density between two points in space, with
//                 respect to the electron density at node 1.
//
// Special Notes : This version passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJndn1
(double n1, double n2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double dNdn1    = aux2( dV);
  double dDNDXdn1 = -aux1(-dV)/h;
  double dJdn1    = MU*(dNdn1*E + Ut*dDNDXdn1);
  return dJdn1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJndn2
//
// Purpose       : This function returns the derivative of electron current
//                 density between two points in space, with
//                 respect to the electron density at node 1.
//
// Special Notes : This version passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJndn2
(double n1, double n2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double dNdn2    = aux2(-dV);
  double dDNDXdn2 =  aux1(-dV)/h;
  double dJdn2    =  MU*(dNdn2*E + Ut*dDNDXdn2);
  return dJdn2;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::Jp
// Purpose       : Hole current density between two points.
// Special Notes : This version passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::Jp
(double p1, double p2, double E, double u, double h)
{
  double MU     = u;
  double dV     = E*h/(2.0*Ut);
  double p      = p1*aux2(-dV)+p2*aux2(dV);
  double dpdx   = aux1(-dV)*(p2-p1)/h;
  double J      = MU*(p*E-Ut*dpdx);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    Xyce::dout().width(16);
    Xyce::dout().precision(8);
    Xyce::dout().setf(std::ios::scientific);

    Xyce::dout() << std::endl;
    Xyce::dout() << "  MU = "<< MU << std::endl;
    Xyce::dout() << "  p  = "<< p  << " dpdx = "<<dpdx<<" E = "<<E <<std::endl;
    Xyce::dout() << "  dV = "<< dV << " Ut = "<<Ut<<" Jpx = "<<J<<std::endl;
    Xyce::dout() << "  p1 = "<< p1 << " p2 = "<<p2<<" (p2-p1) = "<<(p2-p1)<<std::endl;
    Xyce::dout() << "  h  = "<< h  <<  std::endl;
    Xyce::dout() << "  aux1(-dV) = " << aux1(-dV) << std::endl;
    Xyce::dout() << "  aux2( dV) = " << aux2( dV) << std::endl;
    Xyce::dout() << "  aux2(-dV) = " << aux2( dV) << std::endl;
    Xyce::dout() << std::endl;
  }
#endif

  return J;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJpdV1
//
// Purpose       : This function returns the derivative of hole current
//                 density between two points in space, with
//                 respect to the voltage at node 1.
//
// Special Notes : This function passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJpdV1
(double p1, double p2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double p        = p1*aux2(-dV)+p2*aux2(dV);
  double dEdv1    =  1.0/h;
  double ddVdv1   =  1.0/(2.0*Ut);
  double dNdv1    = p1*(-ddVdv1*pd1aux2(-dV) ) + p2*(ddVdv1*pd1aux2(dV));
  double dDNDXdv1 = (p2-p1)/h * (-ddVdv1) * pd1aux1(-dV);
  double dJdv1    = MU*(dNdv1*E + p*dEdv1 - Ut*dDNDXdv1);
  return dJdv1;
}


//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJpdV2
//
// Purpose       : This function returns the derivative of hole current
//                 density between two points in space, with
//                 respect to the voltage at node 2.
//
// Special Notes : This function passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJpdV2
(double p1, double p2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double p        = p1*aux2(-dV)+p2*aux2(dV);
  double dEdv2    = -1.0/h;
  double ddVdv2   = -1.0/(2.0*Ut);
  double dNdv2    = p1*(-ddVdv2*pd1aux2(-dV) ) + p2*(ddVdv2*pd1aux2(dV));
  double dDNDXdv2 = (p2-p1)/h * (-ddVdv2) * pd1aux1(-dV);
  double dJdv2    = MU*(dNdv2*E + p*dEdv2 - Ut*dDNDXdv2);
  return dJdv2;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJpdn1
//
// Purpose       : This function returns the derivative of hole current
//                 density between two points in space, with
//                 respect to the hole density at node 1.
//
// Special Notes : Maybe I should have called this "dJpdp1"?
//               : This version passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJpdn1
(double p1, double p2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double dNdp1    = aux2(-dV);
  double dDNDXdp1 = -aux1(-dV)/h;
  double dJdn1    = MU*(dNdp1*E - Ut*dDNDXdp1);
  return dJdn1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJpdn2
//
// Purpose       : This function returns the derivative of hole current
//                 density between two points in space, with
//                 respect to the hole density at node 2.
//
// Special Notes : This version passes mobility along an edge.
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/07/03
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJpdn2
(double p1, double p2, double E, double u, double h)
{
  double MU       = u;
  double dV       = E*h/(2.0*Ut);
  double dNdp2    = aux2( dV);
  double dDNDXdp2 =  aux1(-dV)/h;
  double dJdn2    = MU*(dNdp2*E - Ut*dDNDXdp2);
  return dJdn2;
}

// ERK Note:  the "qdep" versions of these functions are intended to address
// charge-dependent (or current-dependent) mobilities, and also to be useful
// for defect carriers (i.e. carriers that have more than a single charge).

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::J_qdep
// Purpose       : Current density between two points.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/20/10
//-----------------------------------------------------------------------------
double DevicePDEInstance::J_qdep
(double n1, double n2, double E, double u, double h, int z)
{
  double charge_number = static_cast<double>(z);
  double MU     = u;
  double dV     = -E*h/(2.0*Ut);
  double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
  double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
  double J      =  MU*((n*E)-Ut*(dndx));

  return J;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::nMidpoint
// Purpose       : Carrier density between two points.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/30/12
//-----------------------------------------------------------------------------
pdeFadType DevicePDEInstance::nMidpoint
(pdeFadType & n1, pdeFadType & n2, pdeFadType & E, double h, int z)
{
  double charge_number = static_cast<double>(z);
  pdeFadType dV     = -E*h/(2.0*Ut);
  pdeFadType arg1    = charge_number*dV;
  pdeFadType arg2    = -charge_number*dV;
  pdeFadType A1 = aux2(arg1);
  pdeFadType A2 = aux2(arg2);
  pdeFadType n = charge_number*(n1*A1+n2*A2);
  return n;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdV1_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the voltage at node 1.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/20/10
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdV1_qdep
(double n1, double n2, double E, double u, double h, int z)
{
  double charge_number   = static_cast<double>(z);
  double MU       = u;
  double dV       = -E*h/(2.0*Ut);
  double n        = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
  double dEdv1    =  1.0/h;
  double ddVdv1   =  -1.0/(2.0*Ut);
  double dNdv1    = charge_number*(n1*charge_number*ddVdv1*pd1aux2(charge_number*dV)
				   -n2*charge_number*ddVdv1*pd1aux2(-charge_number*dV));
  double dDNDXdv1 = (n2-n1)/h * (-charge_number*ddVdv1) * pd1aux1(-charge_number*dV);
  double dJdv1    = MU*((dNdv1*E + n*dEdv1) - Ut*dDNDXdv1);

  return dJdv1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdV2_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the voltage at node 2.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/20/10
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdV2_qdep
(double n1, double n2, double E, double u, double h, int z)
{
  double charge_number   = static_cast<double>(z);
  double MU       = u;
  double dV       = -E*h/(2.0*Ut);
  double n        = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
  double dEdv2    =  -1.0/h;
  double ddVdv2   =  1.0/(2.0*Ut);
  double dNdv2    = charge_number*(n1*charge_number*ddVdv2*pd1aux2(charge_number*dV)
				   -n2*charge_number*ddVdv2*pd1aux2(-charge_number*dV));
  double dDNDXdv2 = (n2-n1)/h * (-charge_number*ddVdv2) * pd1aux1(-charge_number*dV);
  double dJdv2    = MU*((dNdv2*E + n*dEdv2) - Ut*dDNDXdv2);

  return dJdv2;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdV1_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the voltage at node 1.
// Special Notes : This version is used when mu is voltage-dependent.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/27/12
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdV1_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdv1=dJdV1_qdep (n1, n2, E, u.val(), h, z);
  double dudv1=u.dx(0);

  // add in chain rule from mobility derivative.
  if (dudv1 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
    double dJ      =  dudv1*((n*E)-Ut*(dndx));
    dJdv1 += dJ;
  }

  return dJdv1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdV2_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the voltage at node 2.
// Special Notes : This version is used when mu is voltage-dependent.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/27/12
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdV2_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdv2=dJdV2_qdep (n1, n2, E, u.val(), h, z);
  double dudv2=u.dx(1);

  // add in chain rule from mobility derivative.
  if (dudv2 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
    double dJ      =  dudv2*((n*E)-Ut*(dndx));
    dJdv2 += dJ;
  }

  return dJdv2;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdn1_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the density at node 1.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/20/10
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdn1_qdep
(double n1, double n2, double E, double u, double h, int z)
{
  double charge_number = static_cast<double>(z);
  double MU     = u;
  double dV     = -E*h/(2.0*Ut);
  double dNdn1    = charge_number*aux2(charge_number*dV);
  double dDNDXdn1 = -aux1(-charge_number*dV)/h;
  double dJdn1    = MU*(dNdn1*E - Ut*dDNDXdn1);

  return dJdn1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdn2_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the electron density at node 1.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/20/10
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdn2_qdep
(double n1, double n2, double E, double u, double h, int z)
{
  double charge_number = static_cast<double>(z);
  double MU     = u;
  double dV     = -E*h/(2.0*Ut);
  double dNdn2    = charge_number*aux2(-charge_number*dV);
  double dDNDXdn2 = aux1(-charge_number*dV)/h;
  double dJdn2    = MU*(dNdn2*E - Ut*dDNDXdn2);

  return dJdn2;
}

//
//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdn1_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the density at node 1.
//
// Special Notes : This version assumes that the mobility is dependent on n1
//                 as well.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/29/12
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdn1_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdn1=dJdn1_qdep (n1, n2, E, u.val(), h, z);
  double dudn1=0.0;
  if (z < 0) dudn1=u.dx(2);
  else       dudn1=u.dx(4);

  if (dudn1 !=0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;

    dJdn1 += dudn1*((n*E)-Ut*(dndx));
  }

  return dJdn1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdn2_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the electron density at node 1.
//
// Special Notes : This version assumes that the mobility is dependent on n2
//                 as well.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/29/12
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdn2_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdn2=dJdn2_qdep (n1, n2, E, u.val(), h, z);
  double dudn2=0.0;
  if (z < 0) dudn2=u.dx(3);
  else       dudn2=u.dx(5);

  if (dudn2 !=0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;

    dJdn2 += dudn2*((n*E)-Ut*(dndx));
  }

  return dJdn2;
}
//
//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdp1_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the density at node 1.
//
// Special Notes : This version assumes that the mobility is dependent on p1
//                 as well.  In this context "p" is the density of the other
//                 carrier.  The only dependence will be via a carrier-dependent
//                 mobility.  This is a little confusing.  If "J" in this
//                 function is a hole current, then "p" refers to
//                 electrons.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/29/12
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdp1_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdp1=0.0;
  double dudp1=0.0;
  if (z < 0) dudp1=u.dx(4);
  else       dudp1=u.dx(2);

  if (dudp1 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;

    dJdp1 += dudp1*((n*E)-Ut*(dndx));
  }

  return dJdp1;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdp2_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the electron density at node 1.
//
// Special Notes : This version assumes that the mobility is dependent on p2
//                 as well.  In this context "p" is the density of the other
//                 carrier.  The only dependence will be via a carrier-dependent
//                 mobility.  This is a little confusing.  If "J" in this
//                 function is a hole current, then "p" refers to
//                 electrons.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/29/12
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdp2_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdp2=0.0;
  double dudp2=0.0;
  if (z < 0) dudp2=u.dx(5);
  else       dudp2=u.dx(3);

  if (dudp2 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
    dJdp2 += dudp2*((n*E)-Ut*(dndx));
  }

  return dJdp2;
}




//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdbm1_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the dopants at node 1.
// Special Notes : This version is used when mu is dopant-dependent and dopants are variable
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 10/24/13
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdbm1_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdbm1 = 0.0;
  double dudbm1=u.dx(6);

  // add in chain rule from mobility derivative.
  if (dudbm1 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
    double dJ      =  dudbm1*((n*E)-Ut*(dndx));
    dJdbm1 += dJ;
  }

  return dJdbm1;
}


//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdbm2_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the dopants at node 1.
// Special Notes : This version is used when mu is dopant-dependent and dopants are variable
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 10/24/13
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdbm2_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdbm2=0.0;
  double dudbm2=u.dx(8);

  // add in chain rule from mobility derivative.
  if (dudbm2 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
    double dJ      =  dudbm2*((n*E)-Ut*(dndx));
    dJdbm2 += dJ;
  }

  return dJdbm2;
}


//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdpp1_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the dopants at node 1.
// Special Notes : This version is used when mu is dopant-dependent and dopants are variable
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 10/24/13
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdpp1_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdpp1=0.0;
  double dudpp1=u.dx(7);

  // add in chain rule from mobility derivative.
  if (dudpp1 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
    double dJ      =  dudpp1*((n*E)-Ut*(dndx));
    dJdpp1 += dJ;
  }

  return dJdpp1;
}


//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::dJdpp2_qdep
// Purpose       : This function returns the derivative of current
//                 density between two points in space, with
//                 respect to the dopants at node 1.
// Special Notes : This version is used when mu is dopant-dependent and dopants are variable
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 10/24/13
//-----------------------------------------------------------------------------
double DevicePDEInstance::dJdpp2_qdep
(double n1, double n2, double E, const pdeFadType & u, double h, int z)
{
  double dJdpp2=0.0;
  double dudpp2=u.dx(9);

  // add in chain rule from mobility derivative.
  if (dudpp2 != 0.0)
  {
    double charge_number = static_cast<double>(z);
    double dV     = -E*h/(2.0*Ut);
    double n      = charge_number*(n1*aux2(charge_number*dV)+n2*aux2(-charge_number*dV));
    double dndx   = aux1(-charge_number*dV)*(n2-n1)/h;
    double dJ      =  dudpp2*((n*E)-Ut*(dndx));
    dJdpp2 += dJ;
  }

  return dJdpp2;
}





//
//
#if 0
// ----------------------------------------------------------------------------
// Function      : DevicePDEInstance::nsdep
// Purpose       : This function returns an approximate deposition profile
//                 of a step implant driven in an inert environment.
// Special Notes :
//                        1       W/2 + x          W/2 - x
//        nsdep(x,W,Dt) = - (erf(---------) + erf(----------))
//                        2      2*sqrt(Dt)       2*sqrt(Dt)
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DevicePDEInstance::nsdep (double x, double W, double Dt)
{
  double D  = 2.0 * sqrt(Dt);
  double Wh = W / 2.0;
  return 0.5 * (erf((Wh + x)/D) + erf((Wh - x)/D));
}

// ----------------------------------------------------------------------------
// Function      : DevicePDEInstance::ngdep
//
// Purpose       : This function returns an approximate Gaussian deposition.
//
//                 Adapted from a similar function in SGF.
//
// Special Notes : This function is a little flakey, in that it assumes
//                 that we're using a cylindrical geomtry
//                 (x = radius, y = height.) It also assumes that the (0,0)
//                 origin is in the upper-left-hand corner of the mesh.
//
//                 So, y=0.0 is the upper surface of the device, and the
//                 implant is coming from the top (above y=0.0).  Hence,
//                 the (y>0) conditional.
//
//                 Also, the width parameter (W), is set up to be a
//                 diameter about x=0, which is the reason for the 0.5*W -
//                 half of this diameter will impact this radius.
//
//                 The implant is completely flat and constant in the
//                 x-direction, as long as fabs(x) is less than W/2.0.
//                 Beyond W/2.0, the gaussian profile kicks in.
//
//                 The parameters ax and ay are scaling parameters, and
//                 correspond to how much you want the doping to vary with
//                 space.  A typical value for either can be set up as:
//
//                 Ax = ln (Nhi / Nlo)/(Rx*Rx)
//
//                 where:
//
//                   Nhi = the max. level of doping
//                   Nlo = the min. level of doping
//                   Rx  = distance over which this doping should vary.
//
//                   Nhi/Nlo = 10^N, where N = # of orders of magnitude
//                   that should vary between x=0 and x=Rx.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DevicePDEInstance::ngdep
(double x, double y, double W, double ax, double ay)
{
  double xprime = fabs(x) - (0.5 * W);
  return ((xprime <= 0.0) ? 1.0 : exp(-ax*xprime*xprime))*
    ((y      >  0.0) ? 0.0 : exp(-ay*y*y));
}

// ----------------------------------------------------------------------------
// Function      : DevicePDEInstance::ngdep2
//
// Purpose       : This function returns an approximate Gaussian deposition.
//
// Special Notes : This function is a modification of the original ngdep
//                 (see above), and I designed it to address some of the
//                 peculiarities of ngdep.
//
//                 (1) I've gotten rid of the width, W.  I'm just
//                 going to assume that whoever is calling this function
//                 can set that(the constant region) up on their own.
//
//                 (2) I've removed the conditionals cause things to
//                 be set to zero, or one, or whatever, if you are on one
//                 side or another of the suface.  I'm assuming that
//                 whoever calls this function can do that themselves, if
//                 they need to.
//
//                 (3) I've removed the stuff that sets the retVal to zero
//                 for y>0.  Again, this is the user's problem.
//
//                 ax and ay mean the same as they did for the original
//                 ngdep.  (see above).
//
//                 It is possible to use this for the 1D case, pretty
//                 easily.  Set the xflag to false, hold y fixed at zero,
//                 and have x correspond to the 1D mesh locations. (or, set
//                 xflag to true, hold x fixed at zero, and let y
//                 correspond to 1D mesh locations - either way).
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/28/03
// ----------------------------------------------------------------------------
double DevicePDEInstance::ngdep2
(double x, double y, double ax, double ay)
{
  double retVal = exp(-ax*x*x)* exp(-ay*y*y);
  return retVal;
}
#endif

// ----------------------------------------------------------------------------
// Function      : DevicePDEInstance::erf
// Purpose       : This function returns the error function.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DevicePDEInstance::erf(double x)
{
  double t1 = 1.0 / (1.0 + 0.3275911 * fabs(x));
  double t2 = t1 * t1;
  double t3 = t2 * t1;
  double t4 = t3 * t1;
  double t5 = t4 * t1;
  double result = 1.0 - (0.254829592*t1 - 0.284496736*t2 + 1.421413741*t3 -
                         1.453152027*t4 + 1.061405429*t5) * exp(-x*x);
  return (x < 0.0) ? -result : result;
}

// ----------------------------------------------------------------------------
// Function      : DevicePDEInstance::pd1erf
// Purpose       : This function returns the derivative of the error
//                 function with respect to the first variable.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DevicePDEInstance::pd1erf(double x)
{
  double pi =  M_PI;
  return 2.0 / sqrt(pi) * exp(-x*x);
}

// ----------------------------------------------------------------------------
// Function      : DevicePDEInstance::setupOutputName
//
// Purpose       : This function takes the device instance name and creates
//                 an appropriate "outputName" to be used for file outputs.
//
//                 At this point PDE devices are all specified as "Y" devices,
//                 meaning that the device instance name will almost always
//                 start with "Y%PDE%".  Left unchanged, this has been
//                 resulting in (for example) tecplot files named thing like,
//                 "Y%PDE%DIODE1.dat".  Obviously, the Y&PDE& prefix is
//                 not needed, so this function removes it, if it exists.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/05
// ----------------------------------------------------------------------------
void DevicePDEInstance::setupOutputName ()
{

  std::string pdeString("Y%PDE%");
  std::string neutString("Y%NEUTRON%");
  std::string::size_type pos1 = getName().find(pdeString);
  std::string::size_type pos2 = getName().find(neutString);

  if (pos1 != std::string::npos)
  {
    std::string tmp1 = "";
    if (pos1 > 0) tmp1 = getName().substr(0,pos1);
    std::string tmp2 = getName().substr(pos1+6, getName().length()-1);
    outputName = tmp1 + tmp2;
  }
  else if (pos2 != std::string::npos)
  {
    std::string tmp1 = "";
    if (pos2 > 0) tmp1 = getName().substr(0,pos2);
    std::string tmp2 = getName().substr(pos2+10, getName().length()-1);
    outputName = tmp1 + tmp2;
  }
  else
  {
    outputName = getName();
  }

  // Tecplot doesn't like file names with the character, ":", so change all
  // colons to underscores.  I personally don't like "%" characters in
  // filenames, so remove those as well.
  for (int i=0;i<outputName.size();++i)
  {
    if (outputName[i]==':') outputName[i] = '_';
    if (outputName[i]=='%') outputName[i] = '_';
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "outputName = "<< outputName << std::endl;
  }
#endif

}

} // namespace Device
} // namespace Xyce
