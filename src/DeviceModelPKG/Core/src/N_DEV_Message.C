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
// Filename      : $RCSfile: N_DEV_Message.C,v $
//
// Purpose       : Contains the class definition for the N_ERH_ErrorMgr
//                 class.  This class is the error handler for Xyce.
//
// Special Notes : 
//
// Creator       : Eric Keiter, SNL,  Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.2 $
//
// Revision Date  : $Date: 2014/02/27 00:52:16 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <N_DEV_Message.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_UTL_Demangle.h>

namespace Xyce {
namespace Device {

// For DeviceBlock: netlistFileName_, parsedLine_[getNumberOfNodes()].lineNumber_
// For OptionBlock: netlistFileName_, parsedLine[0].lineNumber_
// For ParameterBlock: netlistFileName_, parsedLine[0].lineNumber_

// Currently the model and instance headers are the same, so just use DeviceEntity.  If the output wants to be
// different, then copy the DeviceEntity code to make DeviceModel and DeviceInstance.

struct deviceEntityHeader {
  deviceEntityHeader(const DeviceEntity &device_entity)
    : deviceEntity_(device_entity)
  {}

  const DeviceEntity &        deviceEntity_;
};


std::ostream &operator<<(std::ostream &os, const deviceEntityHeader &x)
{
  os << "Device " << x.deviceEntity_.getEntityType() << " " << x.deviceEntity_.getName();

  return os;
}


ParamWarning::ParamWarning(const DeviceEntity &device_entity)
  : Report::UserWarning()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

ParamError::ParamError(const DeviceEntity &device_entity)
  : Report::UserError()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserWarning::UserWarning(const DeviceEntity &device_entity)
  : Report::UserWarning()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserWarning0::UserWarning0(const DeviceEntity &device_entity)
  : Report::UserWarning0()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserError::UserError(const DeviceEntity &device_entity)
  : Report::UserError()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserError0::UserError0(const DeviceEntity &device_entity)
  : Report::UserError0()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserFatal::UserFatal(const DeviceEntity &device_entity)
  : Report::UserFatal()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

UserFatal0::UserFatal0(const DeviceEntity &device_entity)
  : Report::UserFatal0()
{
  at(device_entity.netlistLocation());

  os() << deviceEntityHeader(device_entity) << ": ";
}

DevelFatal::DevelFatal(const DeviceEntity &device_entity, const char *function_name)
  : Report::DevelFatal()
{
  at(device_entity.netlistLocation());
  in(function_name);

  os() << deviceEntityHeader(device_entity) << ": ";
}

DevelFatal0::DevelFatal0(const DeviceEntity &device_entity, const char *function_name)
  : Report::DevelFatal0()
{
  at(device_entity.netlistLocation());
  in(function_name);

  os() << deviceEntityHeader(device_entity) << ": ";
}

struct deviceHeader {
  deviceHeader(const Device &device_)
    : device_(device_)
  {}

  const Device &        device_;
};

std::ostream &operator<<(std::ostream &os, const deviceHeader &x)
{
  os << "Device " << x.device_.getName();

  return os;
}

UserWarning::UserWarning(const Device &device)
  : Report::UserWarning()
{
  os() << deviceHeader(device) << ": ";
}

UserWarning0::UserWarning0(const Device &device)
  : Report::UserWarning0()
{
  os() << deviceHeader(device) << ": ";
}

UserError::UserError(const Device &device)
  : Report::UserError()
{
  os() << deviceHeader(device) << ": ";
}

UserError0::UserError0(const Device &device)
  : Report::UserError0()
{
  os() << deviceHeader(device) << ": ";
}

UserFatal::UserFatal(const Device &device)
  : Report::UserFatal()
{
  os() << deviceHeader(device) << ": ";
}

UserFatal0::UserFatal0(const Device &device)
  : Report::UserFatal0()
{
  os() << deviceHeader(device) << ": ";
}

DevelFatal::DevelFatal(const Device &device, const char *function_name)
  : Report::DevelFatal()
{
  in(function_name);

  os() << deviceHeader(device) << ": ";
}

DevelFatal0::DevelFatal0(const Device &device, const char *function_name)
  : Report::DevelFatal0()
{
  in(function_name);

  os() << deviceHeader(device) << ": ";
}

void device_assertion_error(const DeviceEntity &device_entity, const std::type_info &type, const char *label) 
{
  DevelFatal0(device_entity).in(demangle(type.name()).c_str()) << "Assertion " << label << " failed";
}

} // namespace Device
} // namespace Xyce
