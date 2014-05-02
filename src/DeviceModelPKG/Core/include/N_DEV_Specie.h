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
// Filename       : $RCSfile: N_DEV_Specie.h,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 07/27/2006
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_DEV_Specie_H
#define N_DEV_Specie_H

#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <string>
#include <N_DEV_Const.h>

namespace Xyce {
namespace Device {

class Specie
{
public:
  Specie(std::string name, double diff_prefac, double act_energy, 
         int charge_state)
    : Name(name),
      DiffusionPrefactor(diff_prefac),
      ActivationEnergy(act_energy),
      ChargeState(charge_state)
  {
  } ;

  inline const std::string & getName() const {return (Name); };
  inline void setName(std::string &name) {Name = name;};
  inline int getChargeState() {return (ChargeState);};
  inline void setChargeState(int chargestate) {ChargeState=chargestate;};
  inline double getDiffPrefactor() {return(DiffusionPrefactor);} ; 
  inline void setDiffPrefactor(double p) {DiffusionPrefactor=p;};
  inline double getActEnergy() {return(ActivationEnergy);};    
  inline void setActEnergy(double Energy) {ActivationEnergy=Energy;};
  double getDiffusionCoefficient(double Temperature);
private:
  std::string Name;
  double DiffusionPrefactor;
  double ActivationEnergy;
  int ChargeState;
};

//-----------------------------------------------------------------------------
// Function      : Specie::getDiffusionCoefficient
// Purpose       : Accessor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
inline double Specie::getDiffusionCoefficient(double Temperature)
{
  return DiffusionPrefactor*exp(-ActivationEnergy/(CONSTboltz*Temperature/CONSTQ));
}

} // namespace Device
} // namespace Xyce

#endif
