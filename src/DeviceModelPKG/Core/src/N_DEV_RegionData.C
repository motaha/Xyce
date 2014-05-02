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
// Filename       : $RCSfile: N_DEV_RegionData.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 07/19/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>
#include <N_DEV_DeviceEntity.h>

#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef Xyce_DEBUG_DEVICE
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif
#endif

// ----------    Xyce Includes  ----------
#include <N_DEV_Const.h>
#include <N_DEV_RegionData.h>


namespace Xyce {
namespace Device {

template<>
ParametricData<RegionData>::ParametricData()
{
  // Set up map for double precision variables:
  addPar("AREA", 1.0e+15, &RegionData::area);
  addPar("XLOC", 0.0, &RegionData::xloc);

  // Set up map for non-double precision variables:
  addPar("NAME", std::string("none"), &RegionData::name);
  addPar("TYPE", std::string("none"), &RegionData::type);
  addPar("FILE", std::string(""), &RegionData::reactionFile);
}

ParametricData<RegionData> &RegionData::getParametricData() {
  static ParametricData<RegionData> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// RegionData functions:
//
//-----------------------------------------------------------------------------
// Function      : RegionData::RegionData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/24/06
//-----------------------------------------------------------------------------
RegionData::RegionData ():
  CompositeParam (getParametricData()),
  name("RXNREGION"),
  outName("RXNREGION"),
  type("JUNCTION"),
  reactionFile("reaction_spec_full"),
  area(1.0e-4),
  xloc(1.83e-4),
  doNothing(false)
{}

//-----------------------------------------------------------------------------
// Function      : RegionData::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/24/06
//-----------------------------------------------------------------------------
void RegionData::processParams()
{
  ParameterMap::const_iterator p_i = getParameterMap().find("TYPE");
  const Descriptor &p = *(*p_i).second;
  ExtendedString tmp = p.value<std::string>(*this);
  p.value<std::string>(*this) = tmp.toLower();
}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : RegionData::operator<<
// Purpose       : "<<" operator
// Special Notes : doesn't print everything.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/23/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const RegionData & rd)
{
  os << " Region Data: name = " << rd.name <<
    " x=" << rd.xloc <<
    " reaction file = " << rd.reactionFile <<
    " type = " << rd.type <<
    std::endl;

  return os;
}
#endif

} // namespace Device
} // namespace Xyce
