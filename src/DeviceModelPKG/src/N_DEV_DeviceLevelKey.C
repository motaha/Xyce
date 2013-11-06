//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_DeviceLevelKey.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 04/20/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>
#include <iostream>

#include <N_DEV_DeviceLevelKey.h>

namespace Xyce {
namespace Device {

std::ostream &operator<<(std::ostream &os, const Xyce::Device::DeviceLevelKey &device_level_key) {
  os << "model " << device_level_key.first << ", level " << device_level_key.second;
}

} // namespace Device
} // namespace Xyce
