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
// Filename       : $RCSfile: N_DEV_DeviceSensitivities.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/15/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.25.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceSensitivities_h
#define Xyce_N_DEV_DeviceSensitivities_h

#include <N_DEV_fwd.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceSensitivities
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------

class DeviceSensitivities
{
public:
  DeviceSensitivities(
     DeviceMgr &               device_manager,
     const DeviceOptions &     device_options);
  ~DeviceSensitivities();

private:
  DeviceSensitivities(const DeviceSensitivities &);
  DeviceSensitivities &operator=(const DeviceSensitivities &);

public:
  bool registerSensParams(const Util::OptionBlock &option_block);

private:
  DeviceMgr &                 deviceManager_;
  const DeviceOptions &       deviceOptions_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceSensitivities N_DEV_DeviceSensitivities;

#endif

