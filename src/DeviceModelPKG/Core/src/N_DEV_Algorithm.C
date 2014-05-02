//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Algorithm.C,v $
//
// Purpose        :
//
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_Algorithm.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>

namespace Xyce {
namespace Device {

DeviceEntity *findDeviceEntity(EntityTypeIdDeviceMap::const_iterator begin, EntityTypeIdDeviceMap::const_iterator end, const std::string &entity_name) {
  for (EntityTypeIdDeviceMap::const_iterator it = begin; it != end; ++it) {
    DeviceEntity *device_entity = (*it).second->findEntity(entity_name);
    if (device_entity)
      return device_entity;
  }

  return 0;
}

struct DeviceModelBackInsertOp: public DeviceModelOp
{
    DeviceModelBackInsertOp(std::back_insert_iterator<std::vector<DeviceModel *> > it)
      : it_(it)
    {}

    virtual bool operator()(DeviceModel *model) {
      (*it_)++ = model;

      return true;
    }

    std::back_insert_iterator<std::vector<DeviceModel *> > it_;
};


void getDeviceModels(const Device &device, std::back_insert_iterator<std::vector<DeviceModel *> > it) {
  DeviceModelBackInsertOp op(it);

  device.forEachModel(op);
}

struct DeviceInstanceBackInsertOp: public DeviceInstanceOp
{
    DeviceInstanceBackInsertOp(std::back_insert_iterator<std::vector<DeviceInstance *> > it)
      : it_(it)
    {}

    virtual bool operator()(DeviceInstance *instance) {
      (*it_)++ = instance;

      return true;
    }

    std::back_insert_iterator<std::vector<DeviceInstance *> > it_;
};

void getDeviceInstances(const Device &device, std::back_insert_iterator<std::vector<DeviceInstance *> > it) {
  DeviceInstanceBackInsertOp op(it);

  device.forEachInstance(op);
}

void getDeviceInstances(const DeviceModel &device_model, std::back_insert_iterator<std::vector<DeviceInstance *> > it) {
  DeviceInstanceBackInsertOp op(it);

  device_model.forEachInstance(op);
}

struct DeviceInstanceNameOp: public DeviceInstanceOp
{
    DeviceInstanceNameOp(std::back_insert_iterator<std::vector<std::string> > it)
      : it_(it)
    {}

    virtual bool operator()(DeviceInstance *instance) {
      (*it_)++ = instance->getName();

      return true;
    }

    std::back_insert_iterator<std::vector<std::string> > it_;
};

void getDeviceNames(const Device &device, std::back_insert_iterator<std::vector<std::string> > it) {
  DeviceInstanceNameOp op(it);

  device.forEachInstance(op);
}

} // namespace Device
} // namespace Xyce

