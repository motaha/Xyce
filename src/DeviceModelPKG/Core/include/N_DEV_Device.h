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
// Filename       : $RCSfile: N_DEV_Device.h,v $
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
// Revision Number: $Revision: 1.101.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_Device_h
#define Xyce_N_DEV_Device_h

// ---------- Standard Includes ----------
#include <string>
#include <vector>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_NoCase.h>
#include <N_DEV_DeviceOptions.h>

// ---------- Forward Declarations ----------
class N_IO_CmdParse;
class N_LAS_Matrix;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : Device
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Device
{
  public:
    enum LinearDevice {NONLINEAR_DEVICE = false, LINEAR_DEVICE = true};

    typedef std::map<std::string, DeviceModel *, LessNoCase> ModelNameMap;
    typedef std::map<std::string, DeviceEntity *, LessNoCase> EntityNameMap;

    Device(SolverState & solver_state, DeviceOptions & device_options)
      : name_("DefaultName"),
        linearDevice_(Device::LINEAR_DEVICE),
        commandLine_(device_options.commandLine),
        numInstances_(0),
        solverState_(solver_state),
        deviceOptions_(device_options),
        modelNameMap_(),
        localEntityNameMap_()
    {}

    virtual ~Device()
    {}

  private:
    Device();
    Device(const Device &);
    Device &operator=(const Device &);

  public:
    const std::string & getName() const {
      return name_;
    }

    void setName(const std::string &name) {
      name_ = name;
    }

    bool isLinearDevice() const {
      return linearDevice_ == LINEAR_DEVICE;
    }

    void setLinearDeviceFlag(LinearDevice linearDevice) {
      linearDevice_ = linearDevice;
    }

    virtual DeviceModel *getDefaultModel() const = 0;

    virtual DeviceModel *addModel(const ModelBlock &model_block) = 0;

    virtual DeviceInstance *addInstance(
      InstanceBlock &   instance_block,
      MatrixLoadData &  matrix_load_data,
      SolverState &     solver_state,
      ExternData &      extern_data,
      DeviceOptions &   device_options,
      std::string &     device_type_name) = 0;

    // functions called by the loader, via the device mgr.
    virtual bool updateSources() {
      return true;
    }

    virtual bool updateState(double * solVec, double * staVec, double * stoVec) {
      return true;
    }

    virtual bool updateSecondaryState(double * staDerivVec, double * stoVec) {
      return true;
    }

    // load functions, residual:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ) {
      return true;
    }

    // load functions, Jacobian:
    virtual bool loadDAEMatrices(N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx) {
      return true;
    }

    // API functions (in this case currently used only by Xygra)
    virtual bool deleteInstance(const std::string & instance_name) = 0;
    virtual std::ostream &printOutModels(std::ostream &os) const = 0;
    virtual std::ostream &printDotOpOutput(std::ostream &os) const = 0;

    virtual bool getDeviceNames(std::vector<std::string> &nameVector) const  = 0;
    virtual void getInstances(std::vector<DeviceInstance *> &instanceVector) const = 0;

    DeviceEntity *findLocalEntity(const std::string &entity_name) {
      EntityNameMap::iterator it = localEntityNameMap_.find(entity_name);
      if (it != localEntityNameMap_.end())
        return (*it).second;

      return 0;
    }

    const DeviceEntity *findLocalEntity(const std::string &entity_name) const {
      EntityNameMap::const_iterator it = localEntityNameMap_.find(entity_name);
      if (it != localEntityNameMap_.end())
        return (*it).second;

      return 0;
    }

  protected:
    DeviceModel *findModel(const std::string &model_name) {
      ModelNameMap::iterator it = modelNameMap_.find(model_name);
      if (it != modelNameMap_.end())
        return (*it).second;

      return 0;
    }

    const DeviceModel *findModel(const std::string &model_name) const {
      ModelNameMap::const_iterator it = modelNameMap_.find(model_name);
      if (it != modelNameMap_.end())
        return (*it).second;

      return 0;
    }

    void insertModel(const std::string &model_name, DeviceModel *model) {
      modelNameMap_.insert(ModelNameMap::value_type(model_name, model));
    }

    void insertLocalEntity(const std::string &entity_name, DeviceEntity *entity) {
      localEntityNameMap_.insert(EntityNameMap::value_type(entity_name, entity));
    }

    int getNumInstances() const {
      return numInstances_;
    }

    void incrementInstanceCount() {
      ++numInstances_;
    }

    const SolverState & getSolverState() const {
      return solverState_;
    }

    const DeviceOptions &getDeviceOptions() const {
      return deviceOptions_;
    }

    SolverState & getSolverState() {
      return solverState_;
    }

    DeviceOptions &getDeviceOptions() {
      return deviceOptions_;
    }

    /**
     * @todo Move these to the DeviceTemplate class to make Device a pure interface.
     */
  private:
    std::string     name_;
    LinearDevice    linearDevice_;
    N_IO_CmdParse & commandLine_;
    int             numInstances_;
    SolverState &   solverState_;
    DeviceOptions   deviceOptions_;
    ModelNameMap    modelNameMap_;
    EntityNameMap   localEntityNameMap_; // Entities local to this processor
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Device N_DEV_Device;

#endif

