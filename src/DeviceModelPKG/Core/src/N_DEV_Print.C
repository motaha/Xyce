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
// Filename       : $RCSfile: N_DEV_Print.C,v $
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
// Revision Number: $Revision: 1.4 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <map>
#include <utility>
#include <vector>
#include <iostream>

#include <N_DEV_Print.h>
#include <N_DEV_Device.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMaster.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : printOutModels
// Purpose       : debugging tool
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 2/03/06
//-----------------------------------------------------------------------------
std::ostream &print(std::ostream &os, const Device &device)
{
  std::vector<DeviceModel *> models;
  getDeviceModels(device, std::back_inserter(models));

  os << std::endl
     << std::endl << section_divider << std::endl
     << "Number of " << device.getName() << " models: " << models.size() << std::endl;

  int i = 0;
  for (std::vector<DeviceModel *>::const_iterator it = models.begin(); it != models.end(); ++it, ++i)
  {
    os << i << ": name = " << (*it)->getName() << " type = " << (*it)->getType() << " level = " << (*it)->getLevel() << std::endl;

    (*it)->printOutInstances(os);
  }
  os << section_divider << std::endl;

  return os;
}

typedef std::map<std::string, std::pair<DeviceModel *, std::vector<DeviceInstance *> > > UglyMap;

struct UglyDeviceModelOp: public DeviceModelOp
{
  UglyDeviceModelOp(const Device &device, UglyMap &model_map)
    : device_(device),
      modelMap_(model_map),
      deviceCount_(0)
    {}

  virtual bool operator()(DeviceModel *model) {
    std::pair<UglyMap::iterator, bool> result = modelMap_.insert(UglyMap::value_type(model->getName(), std::make_pair(model, std::vector<DeviceInstance *>())));

    if (result.second) {
      getDeviceInstances(*model, std::back_inserter((*result.first).second.second));
      if ((*result.first).second.second.empty())
        modelMap_.erase(result.first);
      else {
        deviceCount_ += (*result.first).second.second.size();
        std::sort((*result.first).second.second.begin(), (*result.first).second.second.end(), DeviceEntityCmp());
      }
    }
    return true;
  }

  const Device &        device_;
  UglyMap &             modelMap_;
  int                   deviceCount_;
};



//-----------------------------------------------------------------------------
// Function      : DeviceMaster:printDotOpOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 5/4/12
//-----------------------------------------------------------------------------
std::ostream &printDotOpOutput (std::ostream &os, const Device &device)
{
  static const int stdWidth = 15;

// Output models
  UglyMap model_map;

  UglyDeviceModelOp op(device, model_map);
  device.forEachModel(op);

  os << std::endl
     << "Number of " << device.getName() << " models: " << model_map.size() << std::endl
     << std::setw(stdWidth) << "name ";

  if (!model_map.empty())
  {
    int column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
    {
      if (column_count%8 == 0)
        os << std::endl;
      os << std::setw(stdWidth) << (*it_model).first ;
    }
    os << std::endl;

    os << std::setw(stdWidth) << "type ";
    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
    {
      if (column_count%8 == 0)
        os << std::endl;
      os << std::setw(stdWidth) << (*it_model).second.first->getType() ;
    }
    os << std::endl;

    os << std::setw(stdWidth) << "level ";
    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
    {
      if (column_count%8 == 0)
        os << std::endl;
      os << std::setw(stdWidth) << (*it_model).second.first->getLevel() ;
    }
    os << std::endl;

    const ParameterMap &parMap = (*model_map.begin()).second.first->getParameterMap();

    column_count = 1;
    for (ParameterMap::const_iterator it_parameter = parMap.begin(); it_parameter != parMap.end(); ++it_parameter)
    {
      for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model, ++column_count)
      {
        if (it_model == model_map.begin())
        {
          if (column_count%8 == 0)
            os << std::endl;
          os << std::setw(stdWidth) << (*it_parameter).first;
        }

        os << std::setw(stdWidth);
        printParameter(os, *(*it_model).second.first, (*it_parameter).first, *(*it_parameter).second);
      }
      os << std::endl;
    }
    os << std::endl;

// Output instances
    os << "Number of " << device.getName() << " instances: " << op.deviceCount_ << std::endl
       << std::setw(stdWidth) << "name ";

    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model)
    {
      const std::vector<DeviceInstance *> &instance_list = (*it_model).second.second;

      for (std::vector<DeviceInstance *>::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance, ++column_count)
      {
        if (column_count%8 == 0)
          os << std::endl;
        os << std::setw(stdWidth) << (*it_instance)->getName();
      }
    }
    os << std::endl
       << std::setw(stdWidth) << "model ";


    const DeviceInstance &first_instance = *(*model_map.begin()).second.second.front();

    column_count = 1;
    for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model)
    {
      const DeviceModel &model = *(*it_model).second.first;
      const std::vector<DeviceInstance *> &instance_list = (*it_model).second.second;

      if ( !instance_list.empty () ) {
        for (std::vector<DeviceInstance *>::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance, ++column_count)
        {
          if (column_count%8 == 0)
            os << std::endl;
          os << std::setw(stdWidth) << model.getName();
        }
      }
    }

    os << std::endl;

    // output instance parameters:
    const ParameterMap &parameter_map = first_instance.getParameterMap();
    for (ParameterMap::const_iterator it_parameter = parameter_map.begin() ; it_parameter != parameter_map.end(); ++it_parameter)
    {
      os << std::setw(stdWidth) << it_parameter->first;

      column_count = 1;
      for (UglyMap::const_iterator it_model = model_map.begin(); it_model != model_map.end(); ++it_model)
      {
        const DeviceModel &model = *(*it_model).second.first;
        const std::vector<DeviceInstance *> &instance_list = (*it_model).second.second;

        if ( !instance_list.empty () ) {
          for (std::vector<DeviceInstance *>::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance, ++column_count)
          {
            if (column_count%8 == 0)
              os << std::endl;
            os << std::setw(stdWidth);
            printParameter(os, *(*it_instance), (*it_parameter).first, *(*it_parameter).second);
          }
        }
      }
      os << std::endl;
    }
    os << std::endl;
  }

  return os;
}


} // namespace Device
} // namespace Xyce

