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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DevicePDEInstance.h,v $
//
// Purpose        : This file contains the PDE device instance base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.56.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:31 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DevicePDEInstance_h
#define Xyce_N_DEV_DevicePDEInstance_h

// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif
#include <time.h>

#include <Sacado.hpp>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_MaterialSupport.h>
#include <N_DEV_BernouliSupport.h>
#include <N_DEV_Const.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_DopeInfo.h>
#include <N_DEV_ScalingVars.h>

// ---------- Forward Declarations ----------

typedef Sacado::Fad::SFad<double,10> pdeFadType;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : inverse_fermi_one_half_N
// Purpose       : inverse fermi-dirac integral.  Implemented as a functor.
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 07/01/11
//-----------------------------------------------------------------------------
class inverse_fermi_one_half_N
{
private:
  double d__1, d__2, d__3;
  double a1, a2, a3, a4, a5, a6, a7, a8, x10, y10, yp10, x20, y20, yp20, c1, c2;
  double pi, delx, dely;

public:
  inverse_fermi_one_half_N () // this stuff gets called 1x.
  {
    double c_b2 = 4.0/3.0;
    pi = 2.0*asin(1.0);

    a1 = sqrt(2.0) / 4.0;
    a2 = 0.1875 - sqrt(3.0) / 9.0;
    a3 = sqrt(2.0) * 5.0 / 48.0 + 0.125 - sqrt(6.0) / 9.0;
    a4 = sqrt(2.0) * 5.0 / 32.0 + 1585.0/6912.0 - sqrt(3.0) * 5.0/24.0 - sqrt(5.0) / 25.0;
    d__1 = sqrt(pi) * 3.0/4.0;
    a5 = pow(d__1, c_b2);
    a6 = 4.0/3.0;
    a7 = pi   * pi   / 6.0;
    a8 = 1.0/3.0;
    x10 = 7.5;
    d__3 = x10*x10;
    y10 = log(x10) + a1 * x10 + a2 * (x10*x10) + a3*(x10*x10*x10) + a4*(d__3 * d__3);
    yp10 = 1.0/x10 + a1 + a2 * 2.0 * x10 + a3 * 3.0 * (x10 * x10) + a4 * 4.0 * (x10 * x10 * x10);
    x20 = 8.5;
    y20 = sqrt(a5 * pow(x20, a6) - a7);
    yp20 = 0.5 / sqrt(a5 * pow(x20, a6) - a7) * a6 * a5 * pow(x20, a8);
    delx = 0.5;
    dely = y20 - y10;
    c1 = dely * 0.5 / (delx * delx) - yp10 * 0.75 / delx - yp20 * 0.25 / delx;
    c2 = dely * 0.5 / (delx * delx) - yp20 * 0.75 / delx - yp10 * 0.25 / delx;
  }

  template <typename ScalarT>
  ScalarT operator()(const ScalarT & ratio)
  {
    ScalarT ret_val = 0.0;
    ScalarT tempVal = 0.0;

    // Joyce-Dixon expressions as used in Medici
    if (ratio > 0.0 && ratio <= 7.5)
    {
      tempVal = ratio*ratio;
      ret_val = log(ratio) + a1 * ratio + a2*(ratio*ratio) + a3*(ratio*ratio*ratio) + a4*(tempVal*tempVal);
    }

    // These next two clauses from Sam Myers
    if (ratio > 7.5 && ratio <= 8.0)
    {
      ScalarT diff = ratio - 7.5;
      ret_val = y10 + yp10*diff + c1*(diff*diff);
    }
    if (ratio > 8. && ratio < 8.5)
    {
      ScalarT diff = 8.5-ratio;
      ret_val = y20 - yp20*diff - c2*(diff*diff);
    }
    if (ratio >= 8.5)
    {
      ret_val = sqrt(a5 * pow(ratio, a6) - a7);
    }
    return ret_val;
  }
};

//-----------------------------------------------------------------------------
// Class         : DevicePDEInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DevicePDEInstance : public DeviceInstance
{
public:
  DevicePDEInstance(
     const InstanceBlock &       IB,
     ParametricData<void> &      parametric_data,
     const FactoryBlock &        factory_block);

  virtual ~DevicePDEInstance  () {};

private:
  DevicePDEInstance(const DevicePDEInstance & right);
  DevicePDEInstance &operator=(const DevicePDEInstance & right);

public:
  // Fermi-Dirac integral
  double fermi_one_half_B(double arg)
  {
    double pi = 4.0*atan(1.0);

    double nu_eta = pow(arg, 4.0) + 50.0 +
      33.6*arg*(1.0 - 0.68*exp(-0.17*pow(arg+1.0,2)));

    double xi = 3.0*sqrt(pi)/(4.0*pow(nu_eta,0.375));

    return 1.0/(exp(-arg)+xi);
  }

  double getVoltDepHoleDens (double Vmin, double V, double Na)
  {
    return  Na * exp ( Xycemin(CONSTMAX_EXP_ARG, ((Vmin-V)/Ut)) );
  }

  double getVoltDepElecDens (double Vmax, double V, double Nd)
  {
    return  Nd * exp ( Xycemin(CONSTMAX_EXP_ARG, ((V-Vmax)/Ut)) );
  }

  double aux1 (double x);
  double aux2 (double x);


  double pd1aux1(double x);
  double pd1aux2(double x);

  double Jn (double n1, double n2, double E, double u, double h);

  double dJndV1 (double n1, double n2, double E, double u, double h);
  double dJndV2 (double n1, double n2, double E, double u, double h);
  double dJndn1 (double n1, double n2, double E, double u, double h);
  double dJndn2 (double n1, double n2, double E, double u, double h);

  double Jp (double p1, double p2, double E, double u, double h);

  double dJpdV1 (double p1, double p2, double E, double u, double h);
  double dJpdV2 (double p1, double p2, double E, double u, double h);
  double dJpdn1 (double p1, double p2, double E, double u, double h);
  double dJpdn2 (double p1, double p2, double E, double u, double h);

  // charge dependent current density calculations
  double J_qdep (double n1, double n2, double E, double u, double h, int z);

  pdeFadType aux1 (pdeFadType & x)
  {
    pdeFadType retVal=0.0;
    if      (x < -bernSupport.bp0_MISC) x = -bernSupport.bp0_MISC;
    else if (x >  bernSupport.bp0_MISC) x =  bernSupport.bp0_MISC;

    if (x <= bernSupport.bp0_AUX1) retVal=(x / sinh(x));
    else if (x <= bernSupport.bp1_AUX1) retVal=(1 - x*x/6.0*(1.0 - 7.0*x*x/60.0));
    else retVal=(x / sinh(x));

    return retVal;
  }

  pdeFadType aux2 (pdeFadType & x)
  {
    pdeFadType retVal=0.0;

    if (x <= bernSupport.bp0_AUX2) retVal=(1.0);
    else if (x <= bernSupport.bp1_AUX2) retVal=(1.0 / (1.0 + exp(x)));
    else if (x <= bernSupport.bp2_AUX2) retVal=(exp(-x));
    else retVal=(0.0);

    return retVal;
  }


  pdeFadType nMidpoint(pdeFadType & n1, pdeFadType & n2, pdeFadType & E, double h, int z);

  double J_qdep (double n1, double n2, double E, pdeFadType & u, double h, int z)
  { return J_qdep (n1, n2, E, u.val(), h, z); }

  double dJdV1_qdep (double n1, double n2, double E, double u, double h, int z);
  double dJdV2_qdep (double n1, double n2, double E, double u, double h, int z);
  double dJdn1_qdep (double n1, double n2, double E, double u, double h, int z);
  double dJdn2_qdep (double n1, double n2, double E, double u, double h, int z);

  double dJdV1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdV2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdn1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdn2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdp1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdp2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdbm1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdbm2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdpp1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  double dJdpp2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
  //

#if 0
  double nsdep(double x, double W, double Dt);

  double ngdep(double x, double y, double W, double ax, double ay);
  double ngdep2(double x, double y, double ax, double ay);
#endif

  double erf(double x);
  double pd1erf(double x);

  void setupOutputName ();

  const std::string timeDateStamp();
  const std::string tecplotTimeDateStamp();

public:
  // physical constants:
  double Temp;     // operating temperature           (K)
  double charge;   // electron charge                 (C)
  double kb;       // boltzmann's constant            (J/K)
  double Vt;       // thermal voltage
  double Ut;       // thermal voltage, scaled.
  double e0;       // permittivity of vacuum          (F/cm)
  double eSi;      // relative permittivity of Si
  double eSiO2;    // relative permittivity of SiO2
  double eps;      // permittivity of Si              (F/cm)
  double Ni;       // intrinsic concentration of Si   (cm^-3)
  double h_planck; // planck's constant
  double e_mass;   // electron mass

  // scaling variables:
  double x0_user;  // distance scaling, as set by the user (cm)
  double C0_user;  // concentration scaling, as set by the user (cm^-3)
  double t0_user;  // time scaling, as set by the user (sec)

  ScalingVars scalingVars;

  std::map<std::string, DopeInfo *> dopeInfoMap;

  // photogen, seu variables:
  bool photogenOnFlag;
  double xstart, ystart;    // starting location.
  double xend, yend;        // ending location.
  double intensity;
  double photoA1;
  double photoTstart;
  double photoTstop;
  double photoTd;
  double photoTr;
  double photoTf;
  double photoPw;
  double photoPer;
  int lastPeriodIndex;
  int photoType;
  std::string photoString;
  SourceData * Data_ptr;
  SourceData * DataSaved_ptr;

  double photoA1_old;      // old A1      (at begin of cont. loop)
  double photoA1_final;    // final A1    (at end of cont. loop)
  double photoA1_orig;     // original A1 (at begin of cont. step)
  double photoA1_ramp;     // ramped value of A1 (A1 at current cont. step)
  double photoA1_ramp_old; // ramped value of A1 (A1 at current cont. step)
  double photoA1_Delta;    // 2-level delta.(between ckt iterations)
  double photoA1_DeltaC;   // continuation delta.(between cont. solves)
  double maxPhotoDelta;    // maximum photogen delta.
  bool photoContinuationFinished;

  // continuation parameters:
  double maxVoltDelta;
  bool enableContinuationCalled;
  double continuationAlpha;

  bool sensOn;
  bool sensProcess;
  bool meshSensMod;
  bool dopingSensMod;
  bool photogenSensMod;

  std::string mobModelName;
  bool fieldDependentMobility;
  bool fieldDependentMobilityGiven;
  std::string bulkMaterial;
  bool variablesScaled;

  MaterialSupport matSupport;
  BernouliSupport bernSupport;

  std::string outputName;  // added to remove the Y%PDE% prefix.

  // inverse fermi integral function functor.
  inverse_fermi_one_half_N fdinvObj;

protected:

private:

  template <typename T> int sgn(T val)
  {
    return (val > T(0)) - (val < T(0));
  }

};

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::timeDateStamp_
// Purpose       : get current date and time and format for .PRINT output
// Special Notes : inline
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
inline const std::string DevicePDEInstance::timeDateStamp()
{
  const time_t now = time( NULL );
  char timeDate[ 80 ];

  // format for output
  strftime( timeDate, 80, "TIME='%I:%M:%S %p' DATE='%b %d, %Y' ",
            localtime( &now ) );

  return std::string( timeDate );
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::tecplotTimeDateStamp_
// Purpose       : Get current date and time and format for .PRINT output
// Special Notes : tecplot version of timeDateStamp_.
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 9/6/04
//-----------------------------------------------------------------------------
inline const std::string DevicePDEInstance::tecplotTimeDateStamp()
{
  const time_t now = time( NULL );
  char timeDate[ 80 ];

  // format for output
  strftime( timeDate, 80, "TIME= \" %I:%M:%S %p %b %d, %Y \" ",
            localtime( &now ) );

  return std::string( timeDate );
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DevicePDEInstance  N_DEV_DevicePDEInstance;

#endif // Xyce_N_DEV_DevicePDEInstance_h
