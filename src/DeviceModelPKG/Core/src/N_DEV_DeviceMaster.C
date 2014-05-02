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
// Filename       : $RCSfile: N_DEV_DeviceMaster.C,v $
//
// Purpose        : Provides templated versions of some boilerplate functions
//                  that are device-specific (so they can't easily be included
//                  in the base device, instance, or model classes).
//
// Special Notes  : Much of the functionality of the device classes, like
//                  N_DEV_Capacitor, is simply to manage STL containers
//                  of model and instance pointers.  That management is pretty
//                  generic, so templating that functionality makes sense.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/31/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:27 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_DeviceMaster.h>
#include <N_DEV_Message.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInstance.h>

namespace Xyce {
namespace Device {

/**
 * duplicate_entity_warning reports a duplication of entity names.
 *
 * Currently models and devices can share a name and the current implementation of entityMap_ results in lost
 * information.
 *
 * @param entity        const reference to the entity that is being added
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Feb  4 10:12:14 2014
 */
void duplicate_entity_warning(const DeviceEntity &entity)
{
//  DevelFatal(entity).in("DeviceMaster::addEntity") << "Attempt to add duplicate entity";
  UserWarning(entity) << "Duplicate entity, model and device share the name " << entity.getName();
}

/**
 * instance_must_reference_model_error reports that the type of instance requires that a model be specified
 *
 * @param device                const reference to the device
 * @param model_name            const reference to the model name
 * @param netlist_path          const reference to the netlist path
 * @param line_number           line number in the netlist path
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Feb  4 10:16:23 2014
 */
void instance_must_reference_model_error(const Device &device, const std::string &model_name, const std::string &netlist_path, int line_number)
{
  UserError(device).at(netlist_path, line_number) << model_name << " instance must reference a model";
}

/**
 * could_not_find_model reports that the model name is note defined
 *
 * @param device                const reference to the device
 * @param model_name            const reference to the model name
 * @param instance_name         const reference to the instance name
 * @param netlist_path          const reference to the netlist path
 * @param line_number           line number in the netlist path
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Feb  4 10:18:05 2014
 */
void could_not_find_model(const Device &device, const std::string &model_name, const std::string &instance_name, const std::string &netlist_path, int line_number)
{
  UserError(device).at(netlist_path, line_number) << "Could not find model " << model_name << " which is referenced by instance " << instance_name;
}

/**
 * duplicate_model_warning reports that the model name is duplicated.
 *
 *
 * @param device                const reference to the device
 * @param model_name            const reference to the model name
 * @param netlist_path          const reference to the netlist path
 * @param line_number           line number in the netlist path
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Feb  4 10:18:41 2014
 */
void duplicate_model_warning(const Device &device, const std::string &model_name, const std::string &netlist_path, int line_number)
{
  UserWarning(device).at(netlist_path, line_number) << "Attempted to add model " << model_name << " that already exists, ignoring all but the first definition";
}

} // namespace Device
} // namespace Xyce
