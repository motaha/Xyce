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
// Filename       : $RCSfile: N_DEV_SpecieSource.C,v $
//
// Purpose        : This file contains the details of the dope info class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_DEV_DevicePDEInstance.h>
#include <N_UTL_Expression.h>
#include <N_DEV_DeviceSupport.h>
#include <N_DEV_SpecieSource.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<SpecieSource>::ParametricData()
{
  addPar("NAME", "none", &SpecieSource::name);
}

ParametricData<SpecieSource> &SpecieSource::getParametricData() {
  static ParametricData<SpecieSource> parMap;

  return parMap;
}

// ----------------------------------------------------------------------------
// Function      : SpecieSource::SpecieSource
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
SpecieSource::SpecieSource ()
  : CompositeParam (getParametricData()),
    name("V0")
{}

// ----------------------------------------------------------------------------
// Function      : SpecieSource::processParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date :
// ----------------------------------------------------------------------------
bool SpecieSource::processParam
(Param & ndParam, std::string & param, DevicePDEInstance & di)
{
  bool bsuccess = true;

  return bsuccess;
}

// ----------------------------------------------------------------------------
// Function      : SpecieSource::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date :
// ----------------------------------------------------------------------------
void SpecieSource::processParams()
{}

} // namespace Device
} // namespace Xyce
