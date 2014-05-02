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
// Filename       : $RCSfile: N_DEV_PDE_Electrode.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/21/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14 $
//
// Revision Date  : $Date: 2014/02/24 23:49:18 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <map>
#include <algorithm>
#include <string>
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_DEV_PDE_Electrode.h>
#include <N_UTL_Misc.h>
#include <N_ERH_ErrorMgr.h>
#include <N_DEV_DeviceOptions.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<PDE_1DElectrode>::ParametricData()
{
    addPar("AREA", 0.0, &PDE_1DElectrode::area)
      .setGivenMember(&PDE_1DElectrode::areaGiven);
    addPar("LOCATION", 0.0, &PDE_1DElectrode::location);

    // Set up map for non-double precision variables:
    addPar("SIDE", "left", &PDE_1DElectrode::side)
      .setGivenMember(&PDE_1DElectrode::sideGiven);
    addPar("NAME", "anode", &PDE_1DElectrode::nodeName);
    addPar("BC", "dirichlet", &PDE_1DElectrode::bcName);
    addPar("MATERIAL", "neutral", &PDE_1DElectrode::material);
    addPar("OXIDEBNDRYFLAG", false, &PDE_1DElectrode::oxideBndryFlag);
}

ParametricData<PDE_1DElectrode> &PDE_1DElectrode::getParametricData() {
  static ParametricData<PDE_1DElectrode> parMap;

  return parMap;
}

template<>
ParametricData<PDE_2DElectrode>::ParametricData()
{
    // Set up map for double precision variables:
    addPar("START", 0.0, &PDE_2DElectrode::start);
    addPar("END", 0.0, &PDE_2DElectrode::end);
    addPar("OXTHICK", 0.0, &PDE_2DElectrode::oxthick);
    addPar("OXCHARGE", 0.0, &PDE_2DElectrode::oxcharge);

    // Set up map for non-double precision variables:
    addPar("NAME", "anode", &PDE_2DElectrode::nodeName);
    addPar("BC", "dirichlet", &PDE_2DElectrode::bcName);
    addPar("SIDE", "top", &PDE_2DElectrode::side)
      .setGivenMember(&PDE_2DElectrode::sideGiven);
    addPar("MATERIAL", "neutral", &PDE_2DElectrode::material);
    addPar("OXIDEBNDRYFLAG", false, &PDE_2DElectrode::oxideBndryFlag);
}

ParametricData<PDE_2DElectrode> &PDE_2DElectrode::getParametricData() {
  static ParametricData<PDE_2DElectrode> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DElectrode
// Purpose       : constructor
//
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/07/05
//-----------------------------------------------------------------------------
PDE_2DElectrode::PDE_2DElectrode ()
  : PDE_Electrode(getParametricData()),
    iA     (0),
    iB     (0),
    uLabel (0),
    start  (0.0),
    end    (0.0),
    startGiven  (false),
    endGiven    (false),
    side   ("top"),
    sideGiven   (false)
{}

//-----------------------------------------------------------------------------
// Function      : PDE_2DElectrode:processParams
// Purpose       : process raw parameter data
//
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/11/05
//-----------------------------------------------------------------------------
void PDE_2DElectrode::processParams ()
{
  // lowercase the NAME parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("NAME"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_2DElectrode>(*this, param);
    setValue<std::string, PDE_2DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }

  // lowercase the SIDE parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("SIDE"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_2DElectrode>(*this, param);
    setValue<std::string, PDE_2DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }

  // lowercase the MATERIAL parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("MATERIAL"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_2DElectrode>(*this, param);
    setValue<std::string, PDE_2DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }

  // lowercase the BC parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("BC"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_2DElectrode>(*this, param);
    setValue<std::string, PDE_2DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }
}

//-----------------------------------------------------------------------------
// Function      : PDE_1DElectrode
// Purpose       : constructor
//
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/07/05
//-----------------------------------------------------------------------------
PDE_1DElectrode::PDE_1DElectrode ()
  : PDE_Electrode (getParametricData()),
    area(1.0),
    location(0.0),
    side   ("left"),
    sideGiven (false)
{}

//-----------------------------------------------------------------------------
// Function      : PDE_1DElectrode:processParams
// Purpose       : process raw parameter data
//
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/11/05
//-----------------------------------------------------------------------------
void PDE_1DElectrode::processParams ()
{
  // lowercase the NAME parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("NAME"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_1DElectrode>(*this, param);
    setValue<std::string, PDE_1DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }

  // lowercase the SIDE parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("SIDE"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_1DElectrode>(*this, param);
    setValue<std::string, PDE_1DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }

  // lowercase the MATERIAL parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("MATERIAL"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_1DElectrode>(*this, param);
    setValue<std::string, PDE_1DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }

  // lowercase the BC parameter:
  {
    ParameterMap::const_iterator paramIter = getParameterMap().find(std::string("BC"));
    const Descriptor &param = *(*paramIter).second;

    ExtendedString tmp = getValue<std::string, PDE_1DElectrode>(*this, param);
    setValue<std::string, PDE_1DElectrode>(*this, param, static_cast<std::string>(tmp.toLower()));
  }
}

} // namespace Device
} // namespace Xyce
