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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DeviceBld.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/27/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceBld_h
#define Xyce_N_DEV_DeviceBld_h

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_UTL_Registry.h>
#include <N_DEV_Device.h>

// ----------   Forward Declaration ----------
class N_IO_CmdParse;

namespace Xyce {
namespace Device {

/**
 * The DeviceBuilder class creates new devices.
 *
 * The DeviceBuilder stores the SolverState and DeviceOptions objects and injects them on to the device factories for
 * device construction.
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   3/27/00
 */
class DeviceBuilder
{
public:
  DeviceBuilder(SolverState & ss1, DeviceOptions & do1);

  virtual ~DeviceBuilder()
  {}

  Device *createDeviceByModelType(const int model_type);
  Device *createDeviceByNetlistDeviceType(const std::string &name, int level);

private:
  DeviceBuilder(const DeviceBuilder &right);
  DeviceBuilder();

private:
  SolverState &         solState_;                    ///< Solver state of the system
  DeviceOptions &       devOptions_;                  ///< Device options from the netlist
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceBuilder N_DEV_DeviceBuilder;

#endif
