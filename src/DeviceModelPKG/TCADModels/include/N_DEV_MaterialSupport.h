//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
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
// Filename       : $RCSfile: N_DEV_MaterialSupport.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/19/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:31 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MaterialSupport_h
#define Xyce_N_DEV_MaterialSupport_h

#include <Sacado.hpp>

// ---------- Standard Includes ----------
#include <string>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_UTL_Param.h>
#include <N_ERH_ErrorMgr.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MobInfo
// Purpose       : Mobility function information.
// Special Notes :
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 06/24/03
//-----------------------------------------------------------------------------
template <typename ScalarT>
class MobInfo
{
public:
  MobInfo():
    mobModelName("carr"),
    materialName("si"),
    holeFlag(false),
    fieldDependent(false),
    N(1.0e15),
    Na(1.0e15),
    Nd(1.0e15),
    T( CONSTREFTEMP ), // 300.15 K
    refTemp( CONSTREFTEMP ), // 300.15 K
    p(static_cast<ScalarT>(1.45e10)),
    n(static_cast<ScalarT>(1.45e10)),
    epar(static_cast<ScalarT>(0.0)),
    eperp(static_cast<ScalarT>(1.5e4))
  {};

public:
  std::string mobModelName;
  std::string materialName;
  bool holeFlag;
  bool fieldDependent;
  ScalarT N;
  ScalarT Na;
  ScalarT Nd;
  double T;
  double refTemp;
  ScalarT p;
  ScalarT n;
  ScalarT epar;
  ScalarT eperp;
};

//-----------------------------------------------------------------------------
// Class         : MaterialSupport
// Purpose       : Class which contains materials-related data and
//                 functions.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class MaterialSupport
{
public:
  MaterialSupport ();

  MaterialSupport (const MaterialSupport & right);

  virtual ~MaterialSupport ();

  template <typename ScalarT>
  ScalarT calcMob (MobInfo<ScalarT> & min);

  template <typename ScalarT>
  ScalarT calcAnalyticMob (MobInfo<ScalarT> & min);

  template <typename ScalarT>
  ScalarT calcAroraMob (MobInfo<ScalarT> & min);

  template <typename ScalarT>
  ScalarT calcCarrierMobOld (MobInfo<ScalarT> & min);

  template <typename ScalarT>
  ScalarT calcCarrierMobNew (MobInfo<ScalarT> & min);

  template <typename ScalarT>
  ScalarT calcLombardiMob (MobInfo<ScalarT> & min);

  template <typename ScalarT>
  ScalarT calcPhilipsMob (MobInfo<ScalarT> & min);

  template <typename ScalarT>
  void applyHighFieldMobilityModel(MobInfo<ScalarT> & min, ScalarT & mobil);

  double workfunc(std::string & metal);
  double affin(const std::string & material);
  double bandgap(const std::string & material, double temp);

  double calcLt (bool holeFlag, double conc);

  double calcRsrh (const std::string & material,
                   double ni,
                   double n, double p,
                   double tn, double tp);

  double calcRaug (const std::string & material, double ni, double n, double p);

  double pdRsrhN (const std::string & material, double ni,
                  double n, double p,
                  double tn, double tp);

  double pdRsrhP (const std::string & material, double ni,
                  double n, double p,
                  double tn, double tp);

  double pdRaugN (const std::string & material, double ni, double n, double p);
  double pdRaugP (const std::string & material, double ni, double n, double p);

  double getNi (const std::string & material, double temp);
  double getNi_old (const std::string & material, double temp);

  double getRelPerm (const std::string & material);

  double getEffectiveMassN (const std::string & material);
  double getEffectiveMassP (const std::string & material);

  double getDOSEffectiveMassN (const std::string & material);
  double getDOSEffectiveMassP (const std::string & material);

protected:

private:

public:

protected:

private:
};

// Mobility model functions.
//
// Note, while many of these list me (Eric Keiter) as creator, I copied a lot
// from the Charon mobility models files, which were written mostly by
// Debbie Fixel.
//

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcMob
// Purpose       : This function returns the mobility of electrons and
//                 holes for various materials.
// Special Notes :
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 6/19/03
// ----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT MaterialSupport::calcMob(MobInfo<ScalarT> & min)
{
  ScalarT mobil = 0.0;

  // Makes the mobility string case-insensitive
  ExtendedString mobility = min.mobModelName;
  mobility.toLower();

  if(mobility=="analytic" || mobility=="caughey-thomas") // Barnes?
  {
    mobil=calcAnalyticMob (min);
  }
  else if(mobility=="arora")
  {
    mobil=calcAroraMob (min);
  }
  else if(mobility=="carr")
  {
    mobil=calcCarrierMobOld (min);
  }
  else if(mobility=="carrier") // this one matches Charon
  {
    mobil=calcCarrierMobNew (min);
  }
  else if(mobility=="surface" || mobility=="lombardi")
  {
    mobil=calcLombardiMob (min);
  }
  else if(mobility=="philips")
  {
    mobil= calcPhilipsMob (min);
  }
  else  // model not recognized:
  {
    Report::UserFatal0() << "Mobility model " << mobility << " not recognized.";
  }

  if (min.fieldDependent && fabs(min.epar)>0.0)
  {
    applyHighFieldMobilityModel(min,mobil);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (mobil != 0.0 && !(mobil > 0.0) && !(mobil < 0.0))
  {
    Report::DevelFatal0().in("MaterialSupport::calcMob") << "Mobility calc = nan.";
  }
#endif

  return mobil;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcAnalyticMob
// Purpose       : This function returns the mobility of electrons and
//                 holes for various materials.
// Special Notes :
//
// This model is from the reference by D.M. Caughey and R.E. Thomas
// "Carrier Mobilities in Silicon Empirically Related to
// Doping and Field", Proc. IEEE, Vol 55, pp. 2192-2193, 1967.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/12
// ----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT MaterialSupport::calcAnalyticMob (MobInfo<ScalarT> & min)
{
  ExtendedString mater = min.materialName;
  mater.toLower();
  ScalarT mobil=0.0;
  double mun_min, mup_min;
  double mun_max, mup_max;
  double nun, nup, xin, xip, nrefn, nrefp, alphan, alphap;

  if(mater=="si") // silicon
  {
    mun_min = 55.24;	  mup_min = 49.7;
    mun_max = 1429.23;	mup_max = 479.37;
    nrefn = 1.072e17;	  nrefp = 1.606e17;
    nun = -2.3;	        nup = -2.2;
    xin = -3.8;	        xip = -3.7;
    alphan = 0.733;	    alphap = 0.70;
  }
  else if(mater=="gaas") // gallium arsenide
  {
    mun_min = 0.0;	    mup_min = 0.0;
    mun_max = 8500.0;	  mup_max = 400.0;
    nrefn = 1.69e17;	  nrefp = 2.75e17;
    nun = -1.0;	        nup = -2.1;
    xin = 0.0;	        xip = 0.0;
    alphan = 0.436;	    alphap = 0.395;
  }
  else if(mater=="sio2") // silicon dioxide
  {
    mun_min = 1e1;    mup_min = 1e-5;
    mun_max = 2e1;    mup_max = 1e-5;
    nrefn = 1.072e17; nrefp = 1.606e17;
    nun = -2.3;       nup = -2.2;
    xin = -3.8;       xip = -3.7;
    alphan = 0.733;   alphap = 0.70;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    mun_min = 497.0;  mup_min = 0.0;
    mun_max = 2.41e4; mup_max = 480.0;
    nrefn = 1.0e17;   nrefp = 1.0e30;
    nun = 0.0;        nup = 0.0;
    xin = 0.0;        xip = 0.0;
    alphan = 1.0;     alphap = 1.0;

  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    mun_min = 4000.0;  mup_min = 0.0;
    mun_max = 2.73e4;  mup_max = 480.0;
    nrefn = 3.63e17;   nrefp = 1.0e30;
    nun = 0.0;         nup = 0.0;
    xin = 0.0;         xip = 0.0;
    alphan = 1.0;      alphap = 1.0;
  }
  else if (mater=="inp")
  {
    mun_min = 497.0;  mup_min = 0.0;
    mun_max = 2.41e4; mup_max = 480.0;
    nrefn = 1.0e17;   nrefp = 1.0e30;
    nun = 0.0;        nup = 0.0;
    xin = 0.0;        xip = 0.0;
    alphan = 1.0;     alphap = 1.0;
  }
  else if (mater=="ingap")
  {
    mun_min = 0.95;  mup_min = 0.0;
    mun_max = 200.0; mup_max = 150.0;
    nrefn = 1.0e17;  nrefp = 1.0e30;
    nun = 0.0;       nup = 0.0;
    xin = 0.0;       xip = 0.0;
    alphan = 1.0;    alphap = 1.0;
  }
  else
  {
    Report::UserFatal0() << "Analytic (Caughy-Thomas) mobility model not supported for " << mater;
  }

  // hole mobility
  if(min.holeFlag)
  {
    mobil = mup_min +
      (mup_max*pow((min.T/min.refTemp),nup) - mup_min)/
      (1.0 + pow((min.T/min.refTemp),xip)*pow((min.N/nrefp),alphap));
  }
  // electron mobility
  else
  {
    mobil = mun_min +
      (mun_max*pow((min.T/min.refTemp),nun) - mun_min)/
      (1.0 + pow((min.T/min.refTemp),xin)*pow((min.N/nrefn),alphan));
  }

  return mobil;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcAroraMob
// Purpose       : This function returns the mobility of electrons and
//                 holes for various materials.
// Special Notes :
//
// Arora, Hauser, and Roulston,
// "Electron and Hole Mobilities in Silicon as a Function of
// Concentration and Temperature," IEEE Transactions on Electron
// Devices, Vol. ED-29, pp.292-295, 1967.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/12
// ----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT MaterialSupport::calcAroraMob (MobInfo<ScalarT> & min)
{
  ExtendedString mater = min.materialName;
  mater.toLower();
  ScalarT mobil=0.0;
  double mun1_aro, mup1_aro, mun2_aro, mup2_aro, an_arora, ap_arora;
  double cn_arora, cp_arora, exn1_aro, exp1_aro, exn2_aro, exp2_aro;
  double exn3_aro, exp3_aro, exn4_aro, exp4_aro;
  double alphan, alphap;

  if(mater=="si") // silicon
  {
    mun1_aro = 88.0;     mup1_aro = 54.3;
    mun2_aro = 1252.0;   mup2_aro = 407.0;
    an_arora = 0.88;     ap_arora = 0.88;
    cn_arora = 1.26e17;  cp_arora = 2.35e17;
    exn1_aro = -0.57;    exp1_aro = -0.57;
    exn2_aro = -2.33;    exp2_aro = -2.33;
    exn3_aro = 2.4;      exp3_aro = 2.4;
    exn4_aro = -0.146;   exp4_aro = -0.146;
  }
  else if(mater=="gaas") // gallium arsenide
  {
    mun1_aro = 8.5e3; mup1_aro = 4e2;
    mun2_aro = 0.0;   mup2_aro = 0.0;
    an_arora = 0.0;   ap_arora = 0.0;
    cn_arora = 1.26e17; cp_arora = 2.35e17;
    exn1_aro = -5.7e-1; exp1_aro = 0.0;
    exn2_aro = 0.0;     exp2_aro = 0.0;
    exn3_aro = 0.0;     exp3_aro = 0.0;
    exn4_aro = 0.0;     exp4_aro = 0.0;
  }
  else if(mater=="sio2") // silicon dioxide
  {
    mun1_aro = 1e1;      mup1_aro = 1e-5;
    mun2_aro = 2e1;      mup2_aro = 2e-5;
    an_arora = 0.88;     ap_arora = 0.88;
    cn_arora = 1.26e17;  cp_arora = 2.35e17;
    exn1_aro = -0.57;    exp1_aro = -0.57;
    exn2_aro = -2.33;    exp2_aro = -2.33;
    exn3_aro = 2.4;      exp3_aro = 2.4;
    exn4_aro = -0.146;   exp4_aro = -0.146;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    mun1_aro = 2.41e4;	  mup1_aro = 480.0;
    mun2_aro = 0.0;	      mup2_aro = 0.0;
    an_arora = 1.0;	      ap_arora = 1.0;
    cn_arora = 1.0e20;	  cp_arora = 1.0e20;
    exn1_aro = 0.0;	      exp1_aro = 0.0;
    exn2_aro = 0.0;	      exp2_aro = 0.0;
    exn3_aro = 0.0;	      exp3_aro = 0.0;
    exn4_aro = 0.0;	      exp4_aro = 0.0;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    mun1_aro = 2.73e4;	  mup1_aro = 480.0;
    mun2_aro = 0.0;	      mup2_aro = 0.0;
    an_arora = 1.0;	      ap_arora = 1.0;
    cn_arora = 1.0e20;	  cp_arora = 1.0e20;
    exn1_aro = 0.0;	      exp1_aro = 0.0;
    exn2_aro = 0.0;	      exp2_aro = 0.0;
    exn3_aro = 0.0;	      exp3_aro = 0.0;
    exn4_aro = 0.0;	      exp4_aro = 0.0;
  }
  else if (mater=="inp")
  {
    mun1_aro = 2.41e4;	mup1_aro = 480.0;
    mun2_aro = 0.0;	    mup2_aro = 0.0;
    an_arora = 1.0;	    ap_arora = 1.0;
    cn_arora = 1.0e20;	cp_arora = 1.0e20;
    exn1_aro = 0.0;	    exp1_aro = 0.0;
    exn2_aro = 0.0;	    exp2_aro = 0.0;
    exn3_aro = 0.0;	    exp3_aro = 0.0;
    exn4_aro = 0.0;	    exp4_aro = 0.0;
  }
  else if (mater=="ingap")
  {
    mun1_aro = 200.0;  mup1_aro = 150.0;
    mun2_aro = 0.0;    mup2_aro = 0.0;
    an_arora = 1.0;    ap_arora = 1.0;
    cn_arora = 1.0e20; cp_arora = 1.0e20;
    exn1_aro = 0.0;    exp1_aro = 0.0;
    exn2_aro = 0.0;    exp2_aro = 0.0;
    exn3_aro = 0.0;    exp3_aro = 0.0;
    exn4_aro = 0.0;    exp4_aro = 0.0;
  }
  else
  {
    Report::UserFatal0() << "Arora mobility model not supported for " << mater;
  }

  alphan = an_arora*pow((min.T/min.refTemp),exn4_aro);
  alphap = ap_arora*pow((min.T/min.refTemp),exp4_aro);

  if(min.holeFlag)
  {
    mobil = mup1_aro*pow((min.T/min.refTemp),exp1_aro)+
      (mup2_aro*pow((min.T/min.refTemp),exp2_aro))/
      (1.0+pow((min.N/cp_arora*pow((min.T/min.refTemp),exp3_aro)),alphap));
  }
  else
  {
    mobil =
      mun1_aro*pow((min.T/min.refTemp),exn1_aro)+
      (mun2_aro*pow((min.T/min.refTemp),exn2_aro))/
      (1.0+pow((min.N/cn_arora*pow((min.T/min.refTemp),exn3_aro)),alphan));
  }

  return mobil;
}



// old version.

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcCarrierMobOld
// Purpose       : This function returns the mobility of electrons and
//                 holes for various materials.
// Special Notes :
//
// J. M. Dorkel and Ph. Leturcq, “Carrier Mobilities in Silicon Semi-
// Empirically Related to Temperature, Doping, and Injection Level,”
// Solid-State Electronics, 24, pp. 821-825, 1981.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/12
// ----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT MaterialSupport::calcCarrierMobOld (MobInfo<ScalarT> & min)
{
  ExtendedString mater = min.materialName;
  mater.toLower();
  ScalarT mobil=0.0;

  double Al, Bl, Ai, Bi;
  ScalarT  mul, mui;
  ScalarT muc, X;

  if(mater=="si") // silicon
  {
    if(min.holeFlag)
    {
      Al = 495.0;
      Bl = -2.2;
      Ai = 1.00e17;
      Bi = 6.25e14;
    }
    else
    {
      Al = 1430.0;
      Bl = -2.2;
      Ai = 4.61e17;
      Bi = 1.52e15;
    }

  }
  else if(mater=="gaas") // gallium arsenide
  {
    if(min.holeFlag)
    {
      Al = 4.0e2;
      Bl = 0.0;
      Ai = 1.00e17;
      Bi = 6.25e14;
    }
    else
    {
      Al = 8.50e3;
      Bl = 0.0;
      Ai = 4.61e17;
      Bi = 1.52e15;
    }
  }
  else
  {
    Report::UserFatal0() << "Carrier-carrier mobility model not supported for " << mater;
  }

  // lattice scattering term:
  mul = Al*pow((min.T/min.refTemp),Bl);

  // impurity scattering term:
  mui = (Ai*pow(min.T,1.5)/min.N)*(log(1.0+Bi*min.T*min.T/min.N)-
                                   Bi*min.T*min.T/(min.N+Bi*min.T*min.T));

  // carrier-carrier scattering term:
  //  first, make sure n and p are kosher.
  ScalarT N = fabs(min.n); if(N == 0.0) N = 1.0;
  ScalarT P = fabs(min.p); if(P == 0.0) P = 1.0;

  muc = (2.0e17*pow(min.T,1.5)/sqrt(P*N))*
    1.0/(log(1.0+8.28e8*min.T*min.T*pow(P*N,-1.0/3.0)));

  X = sqrt(6.0*mul*(mui+muc)/(mui*muc));
  mobil = mul*(1.025/(1.0+pow(X/1.68,1.43))-0.025);

#ifdef Xyce_DEBUG_DEVICE
  if (mobil != 0.0 && !(mobil > 0.0) && !(mobil < 0.0))
  {
    Xyce::dout() << "mobil is nan" << std::endl;
    Xyce::dout() << "mul = " << mul << std::endl;
    Xyce::dout() << "mui = " << mui << std::endl;
    Xyce::dout() << "muc = " << muc << std::endl;
    Xyce::dout() << "X   = " << X << std::endl;
    Xyce::dout() << "T   = " << min.T << std::endl;
    Xyce::dout() << "n   = " << min.n << std::endl;
    Xyce::dout() << "p   = " << min.p << std::endl;
  }
#endif
  return mobil;
}

// new version

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcCarrierMobNew
// Purpose       : This function returns the mobility of electrons and
//                 holes for various materials.
//
// Special Notes : This function is intended to replace the old one,
//                 calcCarrierMobOld.  This one matches what is in Charon and
//                 Taurus.  This one also has params for more materials.
//
// J. M. Dorkel and Ph. Leturcq, “Carrier Mobilities in Silicon Semi-
// Empirically Related to Temperature, Doping, and Injection Level,”
// Solid-State Electronics, 24, pp. 821-825, 1981.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/12
// ----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT MaterialSupport::calcCarrierMobNew (MobInfo<ScalarT> & min)
{
  ExtendedString mater = min.materialName;
  mater.toLower();
  ScalarT mobil=0.0;

  double a_ccs, b_ccs, a_lic, b_lic, c_lic;
  double ex_lic, mun0_lat, exn_lat;
  double an_iis, bn_iis, mup0_lat, exp_lat;
  double ap_iis, bp_iis;

  if(mater=="si") // silicon
  {
    a_ccs = 1.04e21;
    b_ccs = 7.45e13;
    a_lic = 1.0;
    b_lic = 2.126;
    c_lic = 0.0;
    ex_lic = 0.715;
    mun0_lat = 1430.0;
    exn_lat = 2.3;
    an_iis = 2.4e21;
    bn_iis = 1.37e20;
    mup0_lat = 495.0;
    exp_lat = 2.2;
    ap_iis = 5.2e20;
    bp_iis = 5.63e19;
  }
  else if(mater=="gaas") // gallium arsenide
  {
    a_ccs = 1.04e21;
    b_ccs = 7.45e13;
    a_lic = 1.0;
    b_lic = 0.0;
    c_lic = 0.0;
    ex_lic = 0.0;
    mun0_lat = 8.50e3;
    exn_lat = 0.0;
    an_iis = 2.4e21;
    bn_iis = 1.37e20;
    mup0_lat = 400.0;
    exp_lat = 0.0;
    ap_iis = 5.2e20;
    bp_iis = 5.63e19;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    a_ccs = 1.0e22;
    b_ccs = 1.0e22;
    a_lic = 1.0;
    b_lic = 0.0;
    c_lic = 0.0;
    ex_lic = 0.0;
    mun0_lat = 2.414e4;
    exn_lat = 0.0;
    an_iis = 1.0e22;
    bn_iis = 1.0e22;
    mup0_lat = 480.0;
    exp_lat = 0.0;
    ap_iis = 1.0e22;
    bp_iis = 1.0e22;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    a_ccs = 1.0e22;
    b_ccs = 1.0e22;
    a_lic = 1.0;
    b_lic = 0.0;
    c_lic = 0.0;
    ex_lic = 0.0;
    mun0_lat = 2.73e4;
    exn_lat = 0.0;
    an_iis = 1.0e22;
    bn_iis = 1.0e22;
    mup0_lat = 480.0;
    exp_lat = 0.0;
    ap_iis = 1.0e22;
    bp_iis = 1.0e22;
  }
  else if (mater=="inp")
  {
    a_ccs = 1.0e22;
    b_ccs = 1.0e22;
    a_lic = 1.0;
    b_lic = 0.0;
    c_lic = 0.0;
    ex_lic = 0.0;
    mun0_lat = 2.414e4;
    exn_lat = 0.0;
    an_iis = 1.0e22;
    bn_iis = 1.0e22;
    mup0_lat = 480.0;
    exp_lat = 0.0;
    ap_iis = 1.0e22;
    bp_iis = 1.0e22;
  }
  else if (mater=="ingap")
  {
    a_ccs = 1.0e22;
    b_ccs = 1.0e22;
    a_lic = 1.0;
    b_lic = 0.0;
    c_lic = 0.0;
    ex_lic = 0.0;
    mun0_lat = 200.0;
    exn_lat = 0.0;
    an_iis = 1.0e22;
    bn_iis = 1.0e22;
    mup0_lat = 150.0;
    exp_lat = 0.0;
    ap_iis = 1.0e22;
    bp_iis = 1.0e22;
  }
  else
  {
    Report::UserFatal0() << "Carrier-carrier mobility model not supported for " << mater;
  }

  ScalarT n_impurity = fabs(min.N);

  // make sure n and p are kosher.
  ScalarT N = fabs(min.n); if(N == 0.0) N = 1.0;
  ScalarT P = fabs(min.p); if(P == 0.0) P = 1.0;

  // Carrier-carrier scattering term
  ScalarT muc = a_ccs*pow((min.T/min.refTemp),1.5)/sqrt(P*N)*
    1.0/(log(1.0+b_ccs*pow((min.T/min.refTemp),2.0)*pow(P*N,-1.0/3.0)));

  // electron or hole
  if(min.holeFlag)
  {
    // hole mobility
    // lattice scattering term:
    ScalarT mul_h = mup0_lat*pow((min.T/min.refTemp), -exp_lat);

    // impurity scattering term:
    ScalarT gB_hole = bp_iis*pow((min.T/min.refTemp), 2.0)/(N+P);

    ScalarT mui_h = (ap_iis*pow((min.T/min.refTemp),1.5)/n_impurity)/
      (log(1.0 + gB_hole) - gB_hole/(1.0 + gB_hole));

    ScalarT muic_h = 1.0/(1.0/muc + 1.0/mui_h);

    if (std::fabs(b_lic) < std::numeric_limits<double>::epsilon())
    {
      mobil = mul_h*(a_lic - c_lic);
    }
    else
    {
      mobil = mul_h*(a_lic/(1.0+pow(b_lic*(mul_h/muic_h),ex_lic))-c_lic);
    }
  }
  else
  {
    // electron mobility
    // lattice scattering term:
    ScalarT mul_e = mun0_lat*pow((min.T/min.refTemp), -exn_lat);

    // impurity scattering term:
    ScalarT gB_elec = bn_iis*pow((min.T/min.refTemp), 2.0)/(N+P);

    ScalarT mui_e = (an_iis*pow((min.T/min.refTemp),1.5)/n_impurity)/
      (log(1.0 + gB_elec) - gB_elec/(1.0 + gB_elec));

    ScalarT muic_e = 1.0/(1.0/muc + 1.0/mui_e);

    if (std::fabs(b_lic) < std::numeric_limits<double>::epsilon())
    {
      mobil = mul_e*(a_lic - c_lic);
    }
    else
    {
      mobil = mul_e*(a_lic/(1.0+pow(b_lic*(mul_e/muic_e),ex_lic))-c_lic);
    }
  }

  return mobil;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcLombardiMob
// Purpose       : This function returns the mobility of electrons and
//                 holes for various materials.
// Special Notes :
//
// Lombardi Surface Mobility Model --
// combines mobility of the semiconductor-insulator
// interface with bulk mobility.
// Reference:  Lombardi, Manzini, Saporito, and Vanzi,
// "A physically based mobility model for numerical
// simulation of nonplanar devices", IEEE Trans. on
// Computer-Aided Design of Integrated Circuits and
// Systems, Nov. 1988, vol. 7, no. 11, p.1164-71./
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/12
// ----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT MaterialSupport::calcLombardiMob (MobInfo<ScalarT> & min)
{
  ExtendedString mater = min.materialName;
  mater.toLower();
  ScalarT mobil=0.0;

  ScalarT muac, musr;
  ScalarT mub, mumax;
  ScalarT mun0_lsm, mun1_lsm, mun2_lsm, crn_lsm, csn_lsm;
  double bn_lsm, cn_lsm, dn_lsm, exn1_lsm, exn2_lsm, exn3_lsm;
  double exn4_lsm, exn8_lsm;
  double mup0_lsm, mup1_lsm, mup2_lsm, crp_lsm, csp_lsm;
  double bp_lsm, cp_lsm, dp_lsm, exp1_lsm, exp2_lsm, exp3_lsm;
  double exp4_lsm, exp8_lsm, pc_lsm;

  if(mater=="si")
  {
    mun0_lsm = 52.2;   mup0_lsm = 44.9;
    mun1_lsm = 43.4;   mup1_lsm = 29.0;
    mun2_lsm = 1417.0; mup2_lsm = 470.5;
    crn_lsm = 9.68e16; crp_lsm = 2.23e17;
    csn_lsm = 3.43e20; csp_lsm = 6.1e20;
    bn_lsm = 4.75e7;   bp_lsm = 9.93e6;
    cn_lsm = 1.74e5;   cp_lsm = 8.84e5;
    dn_lsm = 5.82e14;  dp_lsm = 2.05e14;
    exn1_lsm = 0.680;  exp1_lsm = 0.719;
    exn2_lsm = 2.0;    exp2_lsm = 2.0;
    exn3_lsm = 2.5;    exp3_lsm = 2.2;
    exn4_lsm = 0.125;  exp4_lsm = 0.0317;
    exn8_lsm = 2.0;    exp8_lsm = 2.0;
    pc_lsm = 9.23e16;
  }
  else if(mater=="gaas")
  {
    mun0_lsm = 0.0;    mup0_lsm = 0.0;
    mun1_lsm = 0.0;    mup1_lsm = 0.0;
    mun2_lsm = 1e6;    mup2_lsm = 1.0;
    crn_lsm = 9.68e16; crp_lsm = 2.23e17;
    csn_lsm = 0.0;     csp_lsm = 0.0;
    bn_lsm = 1e10;     bp_lsm = 1e10;
    cn_lsm = 0.0;      cp_lsm = 0.0;
    dn_lsm = 1e6;      dp_lsm = 1e6;
    exn1_lsm = 0.0;    exp1_lsm = 0.0;
    exn2_lsm = 0.0;    exp2_lsm = 0.0;
    exn3_lsm = 0.0;    exp3_lsm = 0.0;
    exn4_lsm = 0.0;    exp4_lsm = 0.0;
    exn8_lsm = 0.0;    exp8_lsm = 0.0;
    pc_lsm = 0.0;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    mun0_lsm = 0.0;	    mup0_lsm = 0.0;
    mun1_lsm = 0.0;	    mup1_lsm = 0.0;
    mun2_lsm = 1.0e6;	  mup2_lsm = 1.0e6;
    crn_lsm = 1.0e20;	  crp_lsm = 1.0e20;
    csn_lsm = 0.0;	    csp_lsm = 0.0;
    bn_lsm = 1.0e10;	  bp_lsm = 1.0e10;
    cn_lsm = 0.0;	      cp_lsm = 0.0;
    dn_lsm = 1.0e6;	    dp_lsm = 1.0e6;
    exn1_lsm = 0.0;	    exp1_lsm = 0.0;
    exn2_lsm = 0.0;	    exp2_lsm = 0.0;
    exn3_lsm = 0.0;	    exp3_lsm = 0.0;
    exn4_lsm = 0.0;	    exp4_lsm = 0.0;
    exn8_lsm = 0.0;	    exp8_lsm = 0.0;
    pc_lsm = 0.0;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    mun0_lsm = 0.0;	    mup0_lsm = 0.0;
    mun1_lsm = 0.0;	    mup1_lsm = 0.0;
    mun2_lsm = 1.0e6;	  mup2_lsm = 1.0e6;
    crn_lsm = 1.0e20;	  crp_lsm = 1.0e20;
    csn_lsm = 0.0;	    csp_lsm = 0.0;
    bn_lsm = 1.0e10;	  bp_lsm = 1.0e10;
    cn_lsm = 0.0;	      cp_lsm = 0.0;
    dn_lsm = 1.0e6;	    dp_lsm = 1.0e6;
    exn1_lsm = 0.0;	    exp1_lsm = 0.0;
    exn2_lsm = 0.0;	    exp2_lsm = 0.0;
    exn3_lsm = 0.0;	    exp3_lsm = 0.0;
    exn4_lsm = 0.0;	    exp4_lsm = 0.0;
    exn8_lsm = 0.0;	    exp8_lsm = 0.0;
    pc_lsm = 0.0;
  }
  else if (mater=="inp")
  {
    mun0_lsm = 0.0;	   mup0_lsm = 0.0;
    mun1_lsm = 0.0;	   mup1_lsm = 0.0;
    mun2_lsm = 1.0e6;	 mup2_lsm = 1.0e6;
    crn_lsm = 1.0e20;	 crp_lsm = 1.0e20;
    csn_lsm = 0.0;	   csp_lsm = 0.0;
    bn_lsm = 1.0e10;	 bp_lsm = 1.0e10;
    cn_lsm = 0.0;	     cp_lsm = 0.0;
    dn_lsm = 1.0e6;	   dp_lsm = 1.0e6;
    exn1_lsm = 0.0;	   exp1_lsm = 0.0;
    exn2_lsm = 0.0;	   exp2_lsm = 0.0;
    exn3_lsm = 0.0;	   exp3_lsm = 0.0;
    exn4_lsm = 0.0;	   exp4_lsm = 0.0;
    exn8_lsm = 0.0;	   exp8_lsm = 0.0;
    pc_lsm = 0.0;
  }
  else if (mater=="ingap")
  {
    mun0_lsm = 0.0;   mup0_lsm = 0.0;
    mun1_lsm = 0.0;   mup1_lsm = 0.0;
    mun2_lsm = 1.0e6; mup2_lsm = 1.0e6;
    crn_lsm = 1.0e20; crp_lsm = 1.0e20;
    csn_lsm = 0.0;    csp_lsm = 0.0;
    bn_lsm = 1.0e10;  bp_lsm = 1.0e10;
    cn_lsm = 0.0;     cp_lsm = 0.0;
    dn_lsm = 1.0e6;   dp_lsm = 1.0e6;
    exn1_lsm = 0.0;   exp1_lsm = 0.0;
    exn2_lsm = 0.0;   exp2_lsm = 0.0;
    exn3_lsm = 0.0;   exp3_lsm = 0.0;
    exn4_lsm = 0.0;   exp4_lsm = 0.0;
    exn8_lsm = 0.0;   exp8_lsm = 0.0;
    pc_lsm = 0.0;
  }
  else
  {
    Report::UserFatal0() << "Lobardi surface mobility model not supported for " << mater;
  }

  if(min.holeFlag)
  {
    muac = bp_lsm/min.eperp + cp_lsm*pow(min.N,exp4_lsm)/(min.T*pow(min.eperp,1.0/3.0));
    mumax = mup2_lsm*pow((min.T/min.refTemp),-exp3_lsm);
    mub = mup0_lsm*exp(-pc_lsm/min.N) + mumax/(1.0 + pow(min.N/crp_lsm, exp1_lsm))-
      mup1_lsm/(1.0 + pow(csp_lsm/min.N, exp2_lsm));
    musr = dp_lsm/pow(min.eperp, exp8_lsm);
  }
  else
  {
    muac = bn_lsm/min.eperp + cn_lsm*pow(min.N,exn4_lsm)/(min.T*pow(min.eperp,1.0/3.0));
    mumax = mun2_lsm*pow((min.T/min.refTemp),-exn3_lsm);
    mub = mun0_lsm + (mumax - mun0_lsm)/(1.0 + pow(min.N/crn_lsm, exn1_lsm))-
      mun1_lsm/(1.0 + pow(csn_lsm/min.N, exn2_lsm));
    musr = dn_lsm/pow(min.eperp, exn8_lsm);
  }

  mobil = 1.0/(1.0/muac + 1.0/mub + 1.0/musr);

  return mobil;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcPhilipsMob
// Purpose       : This function returns the mobility of electrons and
//                 holes for various materials.
// Special Notes :
//
// Philips Unified Mobility model -- references by Klaassen
// "A Unified Mobility Model for Device Simulation - I",
// Solid State Electronics, Vol. 35, pp.953-959, 1992.
// "A Unified Mobility Model for Device Simulation - II",
// Solid State Electronics, Vol. 35, pp.961-967, 1992.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/12
// ----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT MaterialSupport::calcPhilipsMob (MobInfo<ScalarT> & min)
{
  ExtendedString mater = min.materialName;
  mater.toLower();
  ScalarT mobil=0.0;

  double mmnn_um, mmxn_um, nrfn_um, alpn_um, tetn_um, nrfd_um, crfd_um;
  double mmnp_um, mmxp_um, nrfp_um, alpp_um, tetp_um, nrfa_um, crfa_um;

  if(mater=="si")
  {
    mmnn_um = 52.2;    mmnp_um = 44.9;
    mmxn_um = 1.417e3; mmxp_um = 470.5;
    nrfn_um = 9.68e16; nrfp_um = 2.23e17;
    alpn_um = 0.68;    alpp_um = 0.719;
    tetn_um = 2.285;   tetp_um = 2.247;
    nrfd_um = 4.0e20;  nrfa_um = 7.2e20;
    crfd_um = 0.21;    crfa_um = 0.5;
  }
  else if(mater=="gaas")
  {
    mmnn_um = 0.0;    mmnp_um = 0.0;
    mmxn_um = 8.5e3;  mmxp_um = 400.0;
    nrfn_um = 1.0e30; nrfp_um = 1.0e30;
    alpn_um = 1.0;    alpp_um = 1.0;
    tetn_um = 0.0;    tetp_um = 0.0;
    nrfd_um = 1.0e30; nrfa_um = 1.0e30;
    crfd_um = 1.0e30; crfa_um = 1.0e30;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    mmnn_um = 0.0;     mmnp_um = 0.0;
    mmxn_um = 2.414e4; mmxp_um = 480.0;
    nrfn_um = 1.0e30;  nrfp_um = 1.0e30;
    alpn_um = 1.0;     alpp_um = 1.0;
    tetn_um = 0.0;     tetp_um = 0.0;
    nrfd_um = 1.0e30;  nrfa_um = 1.0e30;
    crfd_um = 1.0e30;  crfa_um = 1.0e30;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    mmnn_um = 0.0;     mmnp_um = 0.0;
    mmxn_um = 2.725e4; mmxp_um = 400.0;
    nrfn_um = 1.0e30;  nrfp_um = 1.0e30;
    alpn_um = 1.0;     alpp_um = 1.0;
    tetn_um = 0.0;     tetp_um = 0.0;
    nrfd_um = 1.0e30;  nrfa_um = 1.0e30;
    crfd_um = 1.0e30;  crfa_um = 1.0e30;
  }
  else if (mater=="inp")
  {
    mmnn_um = 0.0;     mmnp_um = 0.0;
    mmxn_um = 2.414e4; mmxp_um = 480.0;
    nrfn_um = 1.0e30;  nrfp_um = 1.0e30;
    alpn_um = 1.0;     alpp_um = 1.0;
    tetn_um = 0.0;     tetp_um = 0.0;
    nrfd_um = 1.0e30;  nrfa_um = 1.0e30;
    crfd_um = 1.0e30;  crfa_um = 1.0e30;
  }
  else if (mater=="ingap")
  {
    mmnn_um = 0.0;    mmnp_um = 0.0;
    mmxn_um = 200.0;  mmxp_um = 150.0;
    nrfn_um = 1.0e30; nrfp_um = 1.0e30;
    alpn_um = 1.0;    alpp_um = 1.0;
    tetn_um = 0.0;    tetp_um = 0.0;
    nrfd_um = 1.0e30; nrfa_um = 1.0e30;
    crfd_um = 1.0e30; crfa_um = 1.0e30;
  }
  else
  {
    Report::UserFatal0() << "Philips mobility model not supported for " << mater;
  }

  double me_mo = 1.0;
  double mh_mo = 1.258;
  ScalarT Nd_doping = min.Nd;
  ScalarT Na_doping = min.Na;
  ScalarT n_dens = min.n;
  ScalarT p_dens = min.p;

  if (Nd_doping <= 1.0 ) Nd_doping = 1.0;
  if (Na_doping <= 1.0 ) Na_doping = 1.0;
  if (n_dens <= 1.0 ) n_dens = 1.0;
  if (p_dens <= 1.0 ) p_dens = 1.0;

  ScalarT Ndeff = Nd_doping * ( 1.0 + 1.0/(crfd_um + pow(nrfd_um/Nd_doping,2.0)));
  ScalarT Naeff = Na_doping * ( 1.0 + 1.0/(crfa_um + pow(nrfa_um/Na_doping,2.0)));

  double vsat = 2.4e7/(1.0+0.8*exp(min.T/600.0));
  if(min.holeFlag)
  {
    //hole mobility
    double Tps = min.T/min.refTemp;
    ScalarT Nscp = Ndeff + Naeff + n_dens;
    ScalarT Pp = Tps*Tps * (1.0/( 2.459/(3.97e13*pow(Nscp,-2.0/3.0))
                                  + 3.828/(1.36e20/(n_dens+p_dens)*mh_mo) ));

    // Initial values Pnmin is the Pmin value for 300 K.
    // Initial values of FPnmin and dFPnmin are chosen to start the loop.
    ScalarT Ppmin = 0.2891;
    ScalarT FPpmin = 1.0;
    ScalarT dFPpmin = 1.0;
    // minimizing GPp here
    while( (FPpmin > 0.00001) || (FPpmin < -0.00001) )
    {
      FPpmin = 1.0 - (0.89233/pow(0.41372+Ppmin*pow(1.0/mh_mo*Tps,0.28227),0.19778))
        + (0.005978/pow(Ppmin*pow(mh_mo/Tps,0.72169),1.80618));

      dFPpmin = -0.89233*(-0.19778)*
        pow(0.41372+Ppmin*pow(1.0/mh_mo*Tps,0.28227),-1.19778)*
        pow(1.0/mh_mo*Tps,0.28227) +
        0.005978*(-1.80618)*
        pow(Ppmin*pow(mh_mo/Tps,0.72169),-2.80618)*
        pow(mh_mo/Tps,0.72169);

      Ppmin = Ppmin - FPpmin/dFPpmin;
    }

    ScalarT FPp = (0.7643*pow(Pp,0.6478) + 2.2999 + 6.5502*mh_mo/me_mo)
      /(pow(Pp,0.6478) + 2.3670 - 0.8552*mh_mo/me_mo);

    ScalarT GPp = 1.0 - (0.89233/pow(0.41372+Pp*pow(1.0/mh_mo*Tps,0.28227),0.19778))
      + (0.005978/pow(Pp*pow(mh_mo/Tps,0.72169),1.80618));

    if (Pp < Ppmin)
    {
      GPp = 1.0 - (0.89233/pow(0.41372+Ppmin*pow(1.0/mh_mo*Tps,0.28227),0.19778))
        + (0.005978/pow(Ppmin*pow(mh_mo/Tps,0.72169),1.80618));
    }

    ScalarT Nsceffp = Ndeff*GPp + Naeff + n_dens/FPp;
    ScalarT mulattp = mmxp_um * pow(Tps,-tetp_um);

    ScalarT mu1p = mmxp_um*mmxp_um/(mmxp_um-mmnp_um)*pow(Tps,3.0*alpp_um-1.5);
    ScalarT mu2p = mmxp_um*mmnp_um/(mmxp_um-mmnp_um)*pow(Tps,-0.5);

    ScalarT muDAn = mu1p*(Nscp/Nsceffp)*pow(nrfp_um/Nscp,alpp_um)
      + mu2p*(n_dens+p_dens)/Nsceffp;

    mobil = 1.0/(1.0/mulattp+1.0/muDAn);
  }
  else
  {
    //electron mobility
    double Tns = min.T/min.refTemp;
    ScalarT Nscn = Ndeff + Naeff + p_dens;
    ScalarT Pn = Tns*Tns * (1.0/(2.459/(3.97e13*pow(Nscn,-2.0/3.0))
                                 + 3.828/(1.36e20/(n_dens+p_dens)*me_mo)));

    // Initial values Pnmin is the Pmin value for 300 K.
    // Initial values of FPnmin and dFPnmin are chosen to start the loop.
    ScalarT Pnmin = 0.3246;
    ScalarT FPnmin = 1.0;
    ScalarT dFPnmin = 1.0;
    // minimizing GPn here
    while( (FPnmin > 0.00001) || (FPnmin < -0.00001) )
    {
      FPnmin = 1.0 - (0.89233/pow(0.41372+Pnmin*pow(1.0/me_mo*Tns,0.28227),0.19778))
        + (0.005978/pow(Pnmin*pow(me_mo/Tns,0.72169),1.80618));

      dFPnmin = -0.89233*(-0.19778)*
        pow(0.41372+Pnmin*pow(1.0/me_mo*Tns,0.28227),-1.19778)*
        pow(1.0/me_mo*Tns,0.28227) +
        0.005978*(-1.80618)*
        pow(Pnmin*pow(me_mo/Tns,0.72169),-2.80618)*
        pow(me_mo/Tns,0.72169);

      Pnmin = Pnmin - FPnmin/dFPnmin;
    }

    ScalarT FPn = (0.7643*pow(Pn,0.6478) + 2.2999 + 6.5502*me_mo/mh_mo)
      /(pow(Pn,0.6478) + 2.3670 - 0.8552*me_mo/mh_mo);

    ScalarT GPn = 1.0 - (0.89233/pow(0.41372+Pn*pow(1.0/me_mo*Tns,0.28227),0.19778))
      + (0.005978/pow(Pn*pow(me_mo/Tns,0.72169),1.80618));

    if (Pn < Pnmin)
    {
      GPn = 1.0 - (0.89233/pow(0.41372+Pnmin*pow(1.0/me_mo*Tns,0.28227),0.19778))
        + (0.005978/pow(Pnmin*pow(me_mo/Tns,0.72169),1.80618));
    }

    ScalarT Nsceffn = Ndeff + Naeff*GPn + p_dens/FPn;
    ScalarT mulattn = mmxn_um * pow(Tns,-tetn_um);
    ScalarT mu1n = mmxn_um*mmxn_um/(mmxn_um-mmnn_um)*pow(Tns,3.0*alpn_um-1.5);
    ScalarT mu2n = mmxn_um*mmnn_um/(mmxn_um-mmnn_um)*pow(Tns,-0.5);

    ScalarT muDAp = mu1n*(Nscn/Nsceffn)*pow(nrfn_um/Nscn,alpn_um)
      + mu2n*(n_dens+p_dens)/Nsceffn;

    mobil = 1.0/(1.0/mulattn+1.0/muDAp);
  }


#if 0  // from Myers 1D:
  //subroutine mobility(tk,dv,cd,ca,ce,ch,eu0,hu0,eu,hu)

  // tev = temperature in eV.
  // eu = electron mobility
  // hu = hole mobility
  // tk = temp
  // dv = potential drop
  // cd = donor density
  // ca = acceptor density
  // ce = electron density
  // ch = hole density

  tk300 = tk/300.0d0
    em = 1.0d0
    hm = 1.258d0
    cd1 = (1.0d0+1.0d0/(0.21d0+(4.0d20/(cd+1.0))**2))*(cd+1.0d0)
    ca1 = (1.0d0+1.0d0/(0.50d0+(7.2d20/(ca+1.0))**2))*(ca+1.0d0)
    cscn = cd1+ca1+ch+1.0d0
    cscp = cd1+ca1+ce+1.0d0
    pn = 2.459d0/3.97d13*cscn**0.666666
    pn = pn+3.828d0/1.36d20/em*(ce+ch+1.0d0)
    pn = tk300**2/pn
    pp = 2.459d0/3.97d13*cscp**0.666666
    pp = pp+3.828d0/1.36d20/hm*(ce+ch+1.0d0)
    pp = tk300**2/pp
    fpn = (0.7643d0*pn**0.6478+2.2999d0+6.5502d0*em/hm) /(pn**0.6478+2.3670d0-0.8552d0*em/hm)
    fpp = (0.7643d0*pp**0.6478+2.2999d0+6.5502d0*hm/em) /(pp**0.6478+2.3670d0-0.8552d0*hm/em)
    par1 = (0.41372d0+pn*(tk300/em)**0.28227)**0.19778
    par2 = (pn*(em/tk300)**0.72169)**1.80618
    gpn = 1.0d0-0.89223d0/par1+0.005978d0/par2
    par1 = (0.41372d0+pp*(tk300/hm)**0.28227)**0.19778
    par2 = (pp*(hm/tk300)**0.72169)**1.80618
    gpp = 1.0d0-0.89223d0/par1+0.005978d0/par2
    cscneff = cd1+ca1*gpn+ch/fpn+1.0d0
    cscpeff = ca1+cd1*gpp+ce/fpp+1.0d0
    u1n = 1417.0d0**2/(1417.0d0-52.2d0)*tk300**(3.0*0.68-1.5)
    u2n = 1417.0d0*52.2d0/(1417.0d0-52.2d0)/sqrt(tk300)
    u1p = 470.5d0**2/(470.5d0-44.9d0)*tk300**(3.0*0.719-1.5)
    u2p = 470.5d0*44.9d0/(470.5d0-44.9d0)/sqrt(tk300)
    undap = u1n*cscn/cscneff*(9.68d16/cscn)**0.68 +u2n*(ce+ch+1.0d0)/cscneff
    updan =u1p*cscp/cscpeff*(2.23d17/cscp)**0.719 +u2p*(ce+ch+1.0d0)/cscpeff
    ulattn = 1417.0d0/tk300**2.285
    ulattp = 470.5d0/tk300**2.247
    eu0 = 1.0d0/(1.0d0/ulattn+1.0d0/undap)
    hu0 = 1.0d0/(1.0d0/ulattp+1.0d0/updan)
    vsat = 2.4d7/(1.0d0+0.8d0*exp(tk/600.0d0))
    eu = 1.0d0+(eu0*abs(dv)/vsat)**2
    eu = eu0/sqrt(eu)
    hu = 1.0d0+hu0*abs(dv)/vsat
    hu = hu0/hu
#endif

    return mobil;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::applyHighFieldMobilityModel
// Purpose       : modifies a computed mobility to include high field effects.
//
// Special Notes : There are two possible approximations applied.
//
//                 (1) Caughy-Thomas approximation.
//                 D.M. Caughey and R.E. Thomas
//                 "Carrier Mobilities in Silicon Empirically Related to
//                 Doping and Field", Proc. IEEE, Vol 55, pp. 2192-2193, 1967.
//
//                 mu = mu0/(1.0 + mu0*E/vsat)
//
//                 (2) Barnes approximation.  This one is ONLY applied to III-V
//                 materials, and ONLY to electrons, not holes.  For electrons in
//                 GaAs and other III-V materials, the velocity-field characteristic
//                 has a peak and then falls off.  That gives rise to the negative
//                 differential resistance.  The Barnes expression encapsulates that behavior.
//
//                 Barnes, J J and Lomax, R J and Haddad, G I
//                 "Finite-element simulation of GaAs MESFET's with lateral
//                 doping profiles and submicron gates"
//                 IEEE Transactions on Electron Devices, vol. 23, number 9,
//                 pp. 1042-1048, Sept., 1976
//
//                 mu = (mu0+(vsat/E)*(E/E0)^4)/(1.0 + (E/E0)^4)
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/26/12
// ----------------------------------------------------------------------------
template <typename ScalarT>
void MaterialSupport::applyHighFieldMobilityModel
(MobInfo<ScalarT> & min, ScalarT & mobil)
{
  ScalarT mob0 = mobil;

  ExtendedString mater = min.materialName;
  mater.toLower();

  double vsatn, betan, eon, vsatp, betan_ha, betap, eop, betap_ha, vsn_x1, vsn_x2;
  double en_x1, en_x2;

  if(mater=="si")
  {
    vsatn = 1.035e7; betan = 2.0; eon = 4000.0;
    vsatp = 1.035e7; betan_ha = 2.0; betap = 1.0;
    eop = 4000.0; betap_ha = 2.0; vsn_x1 = 0.0;
    vsn_x2 = 0.0; en_x1 = 0.0; en_x2 = 0.0;

    // Caughy-Thomas for both carriers
    if(min.holeFlag)
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatp);
    }
    else
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatn);
    }
  }
  else if(mater=="ge") // germanium
  {
    vsatn = 1.035e7; betan = 2.0; eon = 4000.0;
    vsatp = 1.035e7; betan_ha = 2.0; betap = 1.0;
    eop = 4000.0; betap_ha = 2.0; vsn_x1 = 0.0;
    vsn_x2 = 0.0; en_x1 = 0.0; en_x2 = 0.0;

    // Caughy-Thomas for both carriers
    if(min.holeFlag)
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatp);
    }
    else
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatn);
    }
  }
  else if(mater=="gaas") // gallium arsenide
  {
    vsatn = 7.7e6; betan = 1.0; eon = 4000.0;
    vsatp = 7.7e6; betan_ha = 2.0; betap = 1.0;
    eop = 4000.0; betap_ha = 2.0; vsn_x1 = 0.0;
    vsn_x2 = 0.0; en_x1 = 0.0; en_x2 = 0.0;

    if(min.holeFlag) // Caughy-Thomas
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatp);
    }
    else  // Barnes for e-
    {
      mobil = (mob0 + (vsatn/min.epar)*pow((min.epar/eon),4.0))/
        (1.0 + pow((min.epar/eon),4.0));
    }
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    vsatn = 4.70e+06; betan = 1.0; eon = 8.40e+03;
    vsatp = 3.00e+6; betan_ha = 2.0; betap = 1.0;
    eop = 8.40e+03; betap_ha = 2.0; vsn_x1 = -0.7493;
    vsn_x2 = 0.0; en_x1 = 4.169; en_x2 = 0.0;

    if(min.holeFlag) // Caughy-Thomas
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatp);
    }
    else  // Barnes for e-
    {
      mobil = (mob0 + (vsatn/min.epar)*pow((min.epar/eon),4.0))/
        (1.0 + pow((min.epar/eon),4.0));
    }
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    vsatn = 8.40e+06; betan = 1.0; eon = 5.07e+03;
    vsatp = 4.80e+06; betan_ha = 2.0; betap = 1.0;
    eop = 5.07e+03; betap_ha = 2.0; vsn_x1 = -1.019;
    vsn_x2 = 0.709; en_x1 = -0.7956; en_x2 = 0.63;

    if(min.holeFlag) // Caughy-Thomas
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatp);
    }
    else  // Barnes for e-
    {
      mobil = (mob0 + (vsatn/min.epar)*pow((min.epar/eon),4.0))/
        (1.0 + pow((min.epar/eon),4.0));
    }
  }
  else if (mater=="inp")
  {
    vsatn = 1.3e+07; betan = 1.0; eon = 1.06e+04;
    vsatp = 6.6e6; betan_ha = 2.0; betap = 1.0;
    eop = 1.06e+04; betap_ha = 2.0; vsn_x1 = -0.332;
    vsn_x2 = 0.0; en_x1 = 5.883; en_x2 = 0.0;

    if(min.holeFlag) // Caughy-Thomas
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatp);
    }
    else  // Barnes for e-
    {
      mobil = (mob0 + (vsatn/min.epar)*pow((min.epar/eon),4.0))/
        (1.0 + pow((min.epar/eon),4.0));
    }
  }
  else if (mater=="ingap")
  {
    vsatn = 7.7e6; betan = 1.0; eon = 4000.0;
    vsatp = 7.7e6; betan_ha = 2.0; betap = 1.0;
    eop = 4000.0; betap_ha = 2.0; vsn_x1 = -0.0341;
    vsn_x2 = 0.0; en_x1 = 1.533; en_x2 = 0.0;

    if(min.holeFlag) // Caughy-Thomas
    {
      mobil = mob0/ (1.0 + mob0*min.epar/vsatp);
    }
    else  // Barnes for e-
    {
      mobil = (mob0 + (vsatn/min.epar)*pow((min.epar/eon),4.0))/
        (1.0 + pow((min.epar/eon),4.0));
    }
  }
  else
  {
    mobil = mob0; // just ignore in this case.  Code can still run.
  }
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MaterialSupport N_DEV_MaterialSupport;

#endif
