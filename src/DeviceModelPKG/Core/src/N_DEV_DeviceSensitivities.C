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
// Filename       : $RCSfile: N_DEV_DeviceSensitivities.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.63 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <vector>
#include <list>

#include <N_DEV_DeviceSensitivities.h>

#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Algorithm.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>
#include <N_UTL_Algorithm.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::DeviceSensitivities
// Purpose       : constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
DeviceSensitivities::DeviceSensitivities(
  DeviceMgr &           device_manager,
  const DeviceOptions & device_options)
  : deviceManager_(device_manager),
    deviceOptions_(device_options)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::~DeviceSensitivities
// Purpose       : destructor
// Special Notes : De-allocates all the devices pointed  to by deviceArray
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
DeviceSensitivities::~DeviceSensitivities()
{}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivites::registerSensParams
//
// Purpose       : This function takes an option block, which contains a
//                 list of user-defined device parameters, and stores
//                 this list in a map.  The map is the "deviceEntityMap",
//                 which maps a parameter string to a device entity (a
//                 device entity is either an instance or a model).
//
// Special Notes : The map is not complete at the end of this function - it
//                 has not filled in the other side(the dePtr side) of the map.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
bool DeviceSensitivities::registerSensParams(const Util::OptionBlock &option_block)
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (deviceOptions_.sensDebugLevel > 0)
  {
    Xyce::dout() << "DeviceSensitivites::registerSensParams called!" <<std::endl;
  }
  int numSensParams = 0;
#endif

  for (std::list<Util::Param>::const_iterator iter = option_block.getParams().begin(); iter != option_block.getParams().end(); ++iter)
  {
    if ( std::string(iter->uTag(), 0, 5) == "PARAM") // this is a vector
    {
      const std::string &tag = iter->stringValue();

#ifdef Xyce_DEBUG_DEVICE
      Xyce::dout() << "name = " << iter->uTag() << "  tag = " << tag << std::endl;
      ++numSensParams;
#endif
      deviceManager_.addDeviceEntity(tag, 0);

    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (deviceOptions_.sensDebugLevel > 0)
  {
    Xyce::dout() << "number of sensitivity parameters = "<< numSensParams << std::endl;
  }
#endif

  return bsuccess;
}

} // namespace Device
} // namespace Xyce
