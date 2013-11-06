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
// Filename       : $RCSfile: N_DEV_DeviceBld.C,v $
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
// Revision Number: $Revision: 1.91.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceBld.h>
#include <N_DEV_Factory.h>
#include <N_DEV_SolverState.h>
#ifdef Xyce_RAD_MODELS
#include <N_DEV_ExtendedModelTypes.h>
#else
#include <N_DEV_ModelTypes.h>
#endif

namespace Xyce {
namespace Device {

/** 
 * DeviceBuilder 
 *
 * 
 * @param ss1 
 * @param do1 
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Mon Aug 12 12:49:53 2013
 */
DeviceBuilder::DeviceBuilder(
  SolverState &   ss1,
  DeviceOptions & do1)
  : solState_(ss1),
    devOptions_(do1)
{}

/**
 * Creates a device based on the model type
 *
 * Allocates a derived device class and returns a base device class pointer which points to it.
 *
 * @param model_type
 *
 * @return pointer to new device
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   3/16/00
 */
Device * DeviceBuilder::createDeviceByModelType(const int model_type)
{
  Device * device_ptr = 0;
  if (Fred::exists(getXyceRegistry(), model_type))
    device_ptr = Fred::create(getXyceRegistry(), model_type)(solState_, devOptions_);

  return device_ptr;
}

/**
 * Creates a device based on the device model name and level
 *
 * Allocates a derived device class and returns a base device class pointer which points to it.
 *
 * @param name device model name
 * @param level device level
 *
 * @return pointer to new device
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Mon Aug 12 12:45:19 2013
 */
Device * DeviceBuilder::createDeviceByNetlistDeviceType(const std::string &name, int level)
{
  Device * device_ptr = 0;
  DeviceLevelKey key(name, level);

  if (Bob::exists(getXyceRegistry2(), key))
    device_ptr = Bob::create(getXyceRegistry2(), key)(solState_, devOptions_);

  return device_ptr;
}

} // namespace Device
} // namespace Xyce
