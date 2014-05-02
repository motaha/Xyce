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
// Filename      : $RCSfile: N_DEV_Message.h,v $
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
// Revision Number: $Revision: 1.5.2.3 $
//
// Revision Date  : $Date: 2014/02/27 00:52:16 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Message_h
#define Xyce_N_DEV_Message_h

#include <N_DEV_fwd.h>

#include <N_ERH_Message.h>

namespace Xyce {
namespace Device {

struct ParamWarning : public Report::UserWarning
{
  ParamWarning(const DeviceEntity &device_entity);
  ParamWarning(const Device &device);
};

struct ParamError : public Report::UserError
{
  ParamError(const DeviceEntity &device_entity);
  ParamError(const Device &device);
};

struct UserWarning : public Report::UserWarning
{
  UserWarning(const DeviceEntity &device_entity);
  UserWarning(const Device &device);
};

struct UserWarning0 : public Report::UserWarning0
{
  UserWarning0(const DeviceEntity &device_entity);
  UserWarning0(const Device &device);
};

struct UserError : public Report::UserError
{
  UserError(const DeviceEntity &device_entity);
  UserError(const Device &device);
};

struct UserError0 : public Report::UserError0
{
  UserError0(const DeviceEntity &device_entity);
  UserError0(const Device &device);
};

struct UserFatal : public Report::UserFatal
{
  UserFatal(const DeviceEntity &device_entity);
  UserFatal(const Device &device);
};

struct UserFatal0 : public Report::UserFatal0
{
  UserFatal0(const DeviceEntity &device_entity);
  UserFatal0(const Device &device);
};

struct DevelFatal : public Report::DevelFatal
{
  DevelFatal(const DeviceEntity &device_entity, const char *function_name = 0);
  DevelFatal(const Device &device, const char *function_name = 0);
};

struct DevelFatal0 : public Report::DevelFatal0
{
  DevelFatal0(const DeviceEntity &device_entity, const char *function_name = 0);
  DevelFatal0(const Device &device, const char *function_name = 0);
};

void device_assertion_error(const DeviceEntity &device_entity, const std::type_info &type, const char *label);


template <class T>
void Assert(bool test, const T &device_instance, const char *label) 
{
  if (!test)
    device_assertion_error(device_instance, typeid(device_instance), label);
}

} // namespace Device
} // namespace Xyce

#define AssertLIDs(cmp) Assert(cmp, *this, #cmp)

#endif // Xyce_N_DEV_Message_h
