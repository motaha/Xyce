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
// Filename       : $RCSfile: N_DEV_PlaceHolder.h,v $
//
// Purpose        : This is a "placeholder" device.  It doesn't do anything.
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
// Revision Number: $Revision: 1.12.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Placeholder_h
#define Xyce_N_DEV_Placeholder_h

// ---------- Standard Includes ----------
#include <list>

// ----------   Xyce Includes   ----------
#include <N_DEV_Device.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : PlaceHolder
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class PlaceHolder : public Device
{
  public:
    ~PlaceHolder();

    static Device *factory(SolverState & ss1, DeviceOptions & do1);

    virtual DeviceModel *addModel(const ModelBlock & MB)
    {}

    virtual DeviceInstance *addInstance(
      InstanceBlock &   instance_block,
      MatrixLoadData &  matrix_load_data,
      SolverState &     solver_state,
      ExternData &      extern_data,
      DeviceOptions &   device_options,
      std::string &     device_type_name)
    {}

    virtual std::ostream &printOutModels(std::ostream &os) const
    {}

    virtual std::ostream &printDotOpOutput(std::ostream &os) const
    {}

    virtual std::ostream &printOutInstances(std::ostream &os) const
    {}

    virtual bool processParams(std::string param = "")
    {}

    virtual bool processInstanceParams(std::string param = "")
    {}

    virtual DeviceModel *getDefaultModel() const {
      return 0;
    }

    virtual bool getDeviceNames(std::vector<std::string> &nameVector) const
    {}

    virtual void getInstances(std::vector<DeviceInstance *> &instanceVector) const
    {}

    virtual bool deleteInstance(const std::string & tmpName)
    {}

  private:
    PlaceHolder(SolverState & solver_state, DeviceOptions & device_options);
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::PlaceHolder N_DEV_PlaceHolder;

#endif

