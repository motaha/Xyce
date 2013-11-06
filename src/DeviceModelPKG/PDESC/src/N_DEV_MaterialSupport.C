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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MaterialSupport.C,v $
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
// Revision Number: $Revision: 1.23.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_MaterialSupport.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : MaterialSupport::MaterialSupport
// Purpose       : constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/19/03
//-----------------------------------------------------------------------------
MaterialSupport::MaterialSupport ()
{

}

//-----------------------------------------------------------------------------
// Function      : MaterialSupport::~MaterialSupport
// Purpose       : constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/19/03
//-----------------------------------------------------------------------------
MaterialSupport::~MaterialSupport ()
{

}

//-----------------------------------------------------------------------------
// Function      : MaterialSupport::MaterialSupport
// Purpose       : copy constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/01/03
//-----------------------------------------------------------------------------
MaterialSupport::MaterialSupport
  (const MaterialSupport & right)
{

}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getEffectiveMassN
// Purpose       : returns effective mass for electrons.
// Special Notes : Relative to free space mass.
//
//                 These are from Appendix 3 of Streetman.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getEffectiveMassN (const string & material)
{
  ExtendedString mater = material;
  mater.toLower();

  double mass=0.0;

  if (mater == "si")
  {
    mass = pow((0.98*0.19*0.19),1.0/3.0); // longitudinal mass = 0.98
                                          // transverse mass = 0.19
  }
  else if (mater == "ge" )
  {
    mass = pow((1.64*0.082*0.082),1.0/3.0); // long. mass = 1.64
                                           // trans. mass = 0.082
  }
  else if (mater == "gaas")
  {
    mass = 0.067;
  }
  // don't have these yet, but what is really needed is the DOS mass, which
  // I do have (see the DOS functions, below)
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
  }
  else if (mater == "inp")
  {
  }
  else
  {
    string msg = material;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getEffectiveMassP
// Purpose       : returns effective mass for holes.
// Special Notes : Relative to free space mass.
//
//                 These are from Appendix 3 of Streetman.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getEffectiveMassP (const string & material)
{
  ExtendedString mater = material;
  mater.toLower();
  double mass=0.0;

  if (mater == "si")
  {
    double mlh = 0.16; // light hole mass
    double mhh = 0.49; // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh,1.5)),2.0/3.0);
  }
  else if (mater == "ge" )
  {
    double mlh = 0.04; // light hole mass
    double mhh = 0.28; // heavy hole mass
    mass = pow((pow(mlh, 1.5) + pow(mhh, 1.5)),2.0/3.0);
  }
  else if (mater == "gaas")
  {
    double mlh = 0.074; // light hole mass
    double mhh = 0.5;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  // don't have these yet, but what is really needed is the DOS mass, which
  // I do have (see the DOS functions, below)
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    double mlh = 0.08; // light hole mass
    double mhh = 0.6;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    double mlh = 0.05; // light hole mass
    double mhh = 0.54; // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else if (mater == "inp")
  {
    double mlh = 0.074; // light hole mass
    double mhh = 0.5;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else
  {
    string msg = material;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getDOSEffectiveMassN
// Purpose       : returns effective mass for electrons for density of states
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
double MaterialSupport::getDOSEffectiveMassN (const string & material)
{
  ExtendedString mater = material;
  mater.toLower();
  double mass = getEffectiveMassN (mater);
  double Mc = 1.0;

  if (mater == "si")
  {
    Mc = pow(6.0,2.0/3.0);
  }
  else if (mater == "gaas")
  {
    Mc = 1.0;
  }
  else if (mater == "ge")
  {
    Mc = 2.0;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    mass = 0.074;
    Mc = 1.0;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    mass = 0.041;
    Mc = 1.0;
  }
  else if (mater == "inp")
  {
    mass = 0.079;
    Mc = 1.0;
  }
  else
  {
    string msg = material;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  mass *= Mc;
  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getDOSEffectiveMassP
// Purpose       : returns effective mass for holes for density of states
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
double MaterialSupport::getDOSEffectiveMassP (const string & material)
{
  ExtendedString mater = material;
  mater.toLower();
  double mass = getEffectiveMassP (mater);

  // Mc is not applied to holes.
#if 0
  double Mc = 1.0;
  if (mater == "si")
  {
    Mc = pow(6.0,(2.0/3.0));
  }
  else if (mater == "gaas")
  {
    Mc = 1.0;
  }
  else if (mater == "ge")
  {
    Mc = pow(2.0,(4.0/3.0));
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    mass = 0.62;
    Mc = 1.0;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    mass = 0.55;
    Mc = 1.0;
  }
  else if (mater == "inp")
  {
    mass = 0.72;
    Mc = 1.0;
  }
  else
  {
    string msg = material;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  mass *= Mc;
#endif

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNi
// Purpose       : returns intrinsic electron concentration.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/24/12
// ----------------------------------------------------------------------------
double MaterialSupport::getNi (const string & material, double temp)
{
  ExtendedString mater = material;
  mater.toLower();
  double ni=0.0;

  double charge(1.602176565e-19);
  double h_planck(6.62606957e-34); // Planck's constant  (in J-s)
  double e_mass (9.10938291e-31); // e- mass in kg.
  double kb (1.3806488e-23); // boltzmann's constant (J/K)
  double kbq = 8.6173324e-5; // boltzmann's constant  (eV/K)
  double dnbnd0 = 2.0*M_PI*e_mass*kb*temp/(h_planck*h_planck);
  dnbnd0 = 2.0*pow(dnbnd0,1.5)/1.0e6;

  double bg = bandgap(mater,temp);

  double mnDOS = getDOSEffectiveMassN(mater);
  double mpDOS = getDOSEffectiveMassP(mater);
  double Nc = dnbnd0 * pow(mnDOS,1.50);
  double Nv = dnbnd0 * pow(mpDOS,1.50);
  ni = sqrt (Nc * Nv) * exp(-bg/(2.0 * kbq * temp));

#ifdef Xyce_DEBUG_DEVICE
  std::cout << "mnDOS = " <<mnDOS   <<std::endl;
  std::cout << "mpDOS = " <<mpDOS   <<std::endl;
  std::cout << "dnbnd0 = " <<  dnbnd0 <<std::endl;
  std::cout << "Nc = " <<  Nc <<std::endl;
  std::cout << "Nv = " <<  Nv <<std::endl;
  std::cout << "Ni = " <<  ni <<std::endl;
#endif
  return ni;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNi_old
// Purpose       : returns intrinsic electron concentration.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getNi_old (const string & material, double temp)
{
  ExtendedString mater = material;
  mater.toLower();
  double ni=0.0;

  if (mater == "si")
  {
    ni = 4.9e15
         * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
         * pow(6.0,0.5) * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
    // ni = 1.25e10;
  }
  else if (mater == "gaas")
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
    //ni = 2.0e6;
  }
  else if (mater == "ge")
  {
    ni = 4.9e15
         * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
         * 2.0 * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
   // ni = 2.5e13;
  }
  // for the next several, as they are all III-V materials, I copied the
  // gaas functions.  I *think* this is correct, as I *think* that Mc is
  // going to be 1.0 for all of these.
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else if (mater == "inp")
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else
  {
    string msg = "MaterialSupport::getNi:  ";
    msg += material;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  return ni;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getRelPerm
// Purpose       : returns relative permitivity
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getRelPerm (const string & material)
{
  ExtendedString mater = material;
  mater.toLower();

  double perm;
  if (mater == "si")
  {
    perm = 11.8;
  }
  else if (mater == "sio2")
  {
    perm = 3.9;
  }
  else if (mater == "ge" )
  {
    perm = 16.0;
  }
  else if (mater == "gaas")
  {
    perm = 13.2;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    perm = 12.5;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    perm = 14.0;
  }
  else if (mater == "inp")
  {
    perm = 12.6;
  }
  else
  {
    string msg = material;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  return perm;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcRsrh
// Purpose       : Calculates schockley-read-hall recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcRsrh
  (const string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pn = Ni*Ni;

  double A = (n*p-pn);
  double B = (tp*(n+Ni)+tn*(p+Ni));

  double arg = CONSTMAX_EXP_ARG;
  if (B >= exp(arg)) B = exp(arg);

  return (A/B);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRsrhN
// Purpose       : Calculates partial derivatives for schockley-read-hall
//                 recombination, with respect to electron density.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRsrhN
  (const string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pdRsrhN;
  double A1, B1, C1;
  double dAdn;
  double dBdn;

  double pn = Ni*Ni;

  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdn = (p);

  C1 = (tp*(n+Ni)+tn*(p+Ni));
  if (C1 >= exp(arg)) C1 = exp(arg);

  B1 = 1.0/C1;
  dBdn = -1.0/(C1*C1) * tp;

  pdRsrhN = dAdn * B1 + dBdn * A1;

  return pdRsrhN;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRsrhP
// Purpose       : Calculates partial derivatives for schockley-read-hall
//                 recombination, with respect to hole density.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRsrhP
  (const string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pdRsrhP;
  double A1, B1, C1;
  double dAdp;
  double dBdp;

  double pn = Ni*Ni;

  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdp = (n);

  C1 = (tp*(n+Ni)+tn*(p+Ni));
  if (C1 >= exp(arg)) C1 = exp(arg);

  B1 = 1.0/C1;
  dBdp = -1.0/(C1*C1) * tn;

  pdRsrhP = dAdp * B1 + dBdp * A1;

  return pdRsrhP;
}


// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcRaug
// Purpose       : Calculates Auger recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 I believe (but am not sure) that the constants Cn and Cp
//                 are material dependent.  That is part of why the
//                 material name is passed in as an argument.  The values
//                 here are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcRaug
  (const string & material, double ni, double n, double p)
{
  double Ni = ni;
  double Cn = 2.8e-31;
  double Cp = 1.2e-31;
  double pn = Ni*Ni;

  double A = (n*p-pn);
  double C = (Cn*n+Cp*p);

  double arg = CONSTMAX_EXP_ARG;
  if (C >= exp(arg)) C = exp(arg);

  return (A*C);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRaugN
// Purpose       : Calculates partial derivative w.r.t. electron density
//                 for Auger recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 I believe (but am not sure) that the constants Cn and Cp
//                 are material dependent.  That is part of why the
//                 material name is passed in as an argument.  The values
//                 here are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRaugN
  (const string & material, double ni, double n, double p)
{
  double Ni = ni;
  double pdRaugN;
  double A1, B1;
  double dAdn;
  double dBdn;

  double Cn = 2.8e-31;
  double Cp = 1.2e-31;
  double pn = Ni*Ni;
  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdn = (p);

  B1 = (Cn*n+Cp*p);
  if (B1 >= exp(arg)) B1 = exp(arg);

  dBdn = Cn;

  pdRaugN = dAdn*B1 + A1*dBdn;

  return pdRaugN;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRaugP
// Purpose       : Calculates partial derivative w.r.t. hole density
//                 for Auger recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 I believe (but am not sure) that the constants Cn and Cp
//                 are material dependent.  That is part of why the
//                 material name is passed in as an argument.  The values
//                 here are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRaugP
  (const string & material, double ni, double n, double p)
{
  double Ni = ni;
  double pdRaugP;
  double A1, B1;
  double dAdp;
  double dBdp;

  double Cn = 2.8e-31;
  double Cp = 1.2e-31;
  double pn = Ni*Ni;
  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdp = (n);

  B1 = (Cn*n+Cp*p);
  if (B1 >= exp(arg)) B1 = exp(arg);

  dBdp = Cp;

  pdRaugP = dAdp*B1 + A1*dBdp;

  return pdRaugP;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::workfunc
// Purpose       : This function returns the workfunction
//                 of various metals
//
// Special Notes :
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/15/03
//----------------------------------------------------------------------------
double MaterialSupport::workfunc(string & metal)
{
  double wkfunc=0.0;

  ExtendedString metalName = metal;
  metalName.toLower ();

  if (metalName=="al")
  {
    wkfunc = 4.10;   //aluminum
  }
  else if (metalName=="ppoly")
  {
    wkfunc = 5.25;   // p+-polysilicon
  }
  else if (metalName=="npoly")
  {
    wkfunc = 4.17;   // n+-polysilicon
  }
  else if (metalName=="mo")
  {
    wkfunc = 4.53;  // molybdenum
  }
  else if (metalName=="w")
  {
    wkfunc = 4.63;  // tungsten
  }
  else if (metalName=="modi")
  {
    wkfunc = 4.80;  // molybdenum disilicide
  }
  else if (metalName=="wdi")
  {
    wkfunc = 4.80;  // tungsten disilicide
  }
  else if (metalName=="cu")
  {
    wkfunc = 4.25;   // copper
  }
  else if (metalName=="pt")
  {
    wkfunc = 5.30;   // platinum
  }
  else if (metalName=="au")
  {
    wkfunc = 4.80;   // gold
  }
  else if (metalName=="neutral")
  {
    wkfunc = 0.0;
  }
  else
  {
    string msg = metalName;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  return wkfunc;
}
//----------------------------------------------------------------------------
// Function      : MaterialSupport::affin
// Purpose       : This function returns the electron affinity
//                 of various semiconductor materials
//
// Special Notes :
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/15/03
//---------------------------------------------------------------------------
double MaterialSupport::affin(const string & material)
{

  double afty=0.0;

  ExtendedString materialName = material;
  materialName.toLower();

  if (materialName=="si")
  {
    afty = 4.17;      // silicon
  }
  else if (materialName=="ge")
  {
    afty = 4.00;     // germanium
  }
  else if (materialName=="gaas")
  {
    afty = 4.07;    // gallium arsenide
  }
  else if (materialName=="sio2")
  {
    afty = 0.97;    // silicon dioxide
  }
  else if (materialName=="nitride")
  {
    afty = 0.97;    // silicon nitride
  }
  else if (materialName=="sapphire")
  {
    afty = 0.97;     // sapphire (also known as aluminum oxide)
  }
  else
  {
    string msg = materialName;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  return afty;
}
//----------------------------------------------------------------------------
// Function      : MaterialSupport::bandgap
// Purpose       : This function returns the electronic bandgap
//                 of various semiconductor materials.
//
// Special Notes : Reference for temperature-dependent semiconductor
//                 materials is "The Standard Thermodynamic Function
//                 of the Formation of Electrons and Holes in Ge, Si,
//                 GaAs, and GaP," by C. D. Thurmond, J. Electrochem. Soc.,
//                 vol. 122, p. 1133, 1975.
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/18/03
//---------------------------------------------------------------------------
double MaterialSupport::bandgap(const string & material,  double temp)
{
  double gap = 0.0;
  ExtendedString materialName = material;
  materialName.toLower();

  if (materialName=="si")   // silicon
  {
    gap = 1.17 - 4.73e-4*pow(temp,2.0)/(temp + 636.0);
  }
  else if (materialName=="ge")  // germanium
  {
    gap = 0.7437 - 4.774e-4*pow(temp,2.0)/(temp + 235);
  }
  else if (materialName=="gaas")  // gallium arsenide
  {
    gap = 1.519 - 5.405e-4*pow(temp,2.0)/(temp + 204);
  }
  else if (materialName=="sio2")  // silicon dioxide
  {
    gap = 9.00;
  }
  else if (materialName=="nitride")  // silicon nitride
  {
    gap = 4.7;
  }
  else if (materialName=="sapphire")  // sapphire
  {
    gap = 4.7;
  }
  else if (materialName=="inalas" || materialName=="alinas") // indium aluminum arsenide
  {
    double con =  1.46;
    double val =  0.0;
    gap = con-val;
  }
  else if (materialName=="ingaas" || materialName=="gainas") // indium galium arsenide
  {
    double con =  0.96;
    double val =  0.21;
    gap = con-val;
  }
  else if (materialName=="inp")  // indium phosphide
  {
    double con =  1.21;
    double val =  -0.14;
    gap = con-val;
  }
  else
  {
    string msg = materialName;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  return gap;
}
// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcLt
// Purpose       : This function calculates carrier lifetimes.
//
// Special Notes : holeFlag parameter indicates electrons or holes:
//                     holeFlag = true  -> holes
//                     holeFlag = false -> electrons
//
// This function assumes that conc is an absolute value.
//
// This function comes from this paper:
//
//     "Analysis of High-Efficiency Silicon Solar Cells",
//     IEEE Transactions on Electron Devices, by Harry T.
//     Weaver and R. D. Nasby, vol. ED-28, no. 5, May 1981.
//
// This is a function that probably should have some material dependence,
// but at the moment it doesn't.  I think all the values in this function
// are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/20/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcLt (bool holeFlag, double conc)
{
  double lt = 0.0;
  double LT0, Nref;

  conc = fabs(conc);

  if (holeFlag)
  {
    LT0   = 3.52e-5;
    Nref  = 7.1e15;
    lt = LT0 / (1.0 + conc / Nref);
  }
  else
  {
    LT0   = 3.95e-4;
    Nref  = 7.1e15;
    lt = LT0 / (1.0 + conc / Nref);
  }
  return lt;
}

} // namespace Device
} // namespace Xyce
