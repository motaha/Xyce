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
// Filename       : $RCSfile: N_DEV_DeviceTemplate.h,v $
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
// Revision Number: $Revision: 1.47.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Device_Template_h
#define Xyce_N_DEV_Device_Template_h

#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>

#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceOptions.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceTemplate
//
// Purpose       : This class contains a lot of boilerplate functions, such
//                 as addInstance, which are needed by every device.  However,
//                 as functions of this nature are device-specific, they
//                 didn't naturally fit as base class Device functions.
//
//                 This class is a template, which takes 2 typenames;
//                 model, instance and model-ptr STL vector.  A specific
//                 example of these would be CapacitorModel and
//                 CapacitorInstance
//
//                 Each derived device class will need to set up its own
//                 version of this template.
//
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
template<class Model, class Instance>
class DeviceTemplate : public Device
{
  protected:
    typedef std::vector<Instance *> InstanceVector;
    typedef std::map<std::string, Model *, LessNoCase> ModelMap;

  public:
    DeviceTemplate(
      const std::string & device_name,
      const std::string & class_name,
      const std::string & default_model_name,
      LinearDevice        linear_device,
      SolverState &       solver_state,
      DeviceOptions &     device_options,
      int                 num_ext_vars = 0)
      : Device(solver_state, device_options),
        className_(class_name),
        defaultModelName_(default_model_name + " (" + device_name + ")"),
        defaultModel_(new Model(ModelBlock(default_model_name + " (" + device_name + ")", ""), solver_state, device_options)),
        numExtVars_(num_ext_vars),
        modelMap_(),
        instanceVector_()
    {
      setName(device_name);
      setLinearDeviceFlag(linear_device);

      // add default model to the model list.
      insertModel(defaultModelName_, defaultModel_);
      modelMap_[defaultModelName_] = defaultModel_;
    }

    DeviceTemplate(
      const std::string & device_name,
      const std::string & class_name,
      const std::string & default_model_name,
      const std::string & model_type,
      LinearDevice        linear_device,
      SolverState &       solver_state,
      DeviceOptions &     device_options)
      : Device(solver_state, device_options),
        className_(class_name),
        defaultModelName_(default_model_name + " (" + device_name + ")"),
        defaultModel_(new Model(ModelBlock(default_model_name + " (" + device_name + ")", model_type), solver_state, device_options)),
        numExtVars_(0),
        modelMap_(),
        instanceVector_()
    {
      setName(device_name);
      setLinearDeviceFlag(linear_device);

      // add default model to the model list.
      insertModel(defaultModelName_, defaultModel_);
      modelMap_[defaultModelName_] = defaultModel_;
    }

    ~DeviceTemplate()
    {
      // loop over models:
      typename ModelMap::iterator firstM = modelMap_.begin();
      typename ModelMap::iterator lastM  = modelMap_.end();
      for (typename ModelMap::iterator iterM = firstM; iterM != lastM; ++iterM)
      {
        delete (*iterM).second;
      }
    }

  private:
    DeviceTemplate(const DeviceTemplate<Model, Instance> &right);
    DeviceTemplate(const Device &);
    DeviceTemplate &operator=(const DeviceTemplate &right);

  protected:
    const InstanceVector &getInstanceVector() const {
      return instanceVector_;
    }

    const ModelMap &getModelMap() const {
      return modelMap_;
    }

    bool isModelRequired() const {
      return Instance::getParametricData().getModelRequired();
    }

    DeviceModel *getModelPointer(const std::string & model_name) const {
      return findModel(model_name);
    }

// Device virtual overrides
    virtual DeviceModel *getDefaultModel() const {
      return defaultModel_;
    }

    virtual bool deleteInstance(const std::string & device_instance_name);

    virtual bool getBreakPoints( std::vector<N_UTL_BreakPoint> & breakPointTimes );

    virtual std::ostream &printOutModels(std::ostream &os) const ;

    virtual std::ostream &printDotOpOutput(std::ostream &os) const;

    virtual DeviceModel *addModel(const ModelBlock & MB);

    virtual DeviceInstance *addInstance(
      InstanceBlock &     instance_block,
      MatrixLoadData &    matrix_load_data,
      SolverState &       solver_state,
      ExternData &        extern_data,
      DeviceOptions &     device_options,
      std::string &       return_default_model_name);

    virtual bool updateSources();
    virtual bool updateState (double * solVec, double * staVec, double * stoVec);
    virtual bool updateSecondaryState (double * staDerivVec, double * stoVec);
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
    virtual bool getDeviceNames(std::vector<std::string> &devNames) const;
    virtual void getInstances(std::vector<DeviceInstance *> &divec) const;

  private:
    const std::string     className_;
    const std::string     defaultModelName_;
    Model * const         defaultModel_;
    int                   numExtVars_;
    ModelMap              modelMap_;
    InstanceVector        instanceVector_;
};

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::addInstance
// Purpose       : This function adds an instance to the list of instances.
//
// Special Notes : All device instances are contained in one of many lists.
//                 Each model owns exactly one list of instances.  If
//                 an instance is added which does not reference a model, or
//                 that references a model that does not exist, then a new
//                 model has to be added.
//
//                 For now, this function will only add a model for instances
//                 that do not specifically reference one.  If an instance
//                 references a non-existant model, this function will return a
//                 fatal error.  (This function should only be called from the
//                 device manager function, AddDeviceInstance, and that
//                 function will not call this one if the instance refers to a
//                 non-existant model.)
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
DeviceInstance * DeviceTemplate<Model, Instance>::addInstance(
  InstanceBlock &       instance_block,
  MatrixLoadData &      matrix_load_data,
  SolverState &         solver_state,
  ExternData &          extern_data,
  DeviceOptions &       device_options,
  std::string &         return_default_model_name)
{
  std::string model_name = instance_block.getModelName();

  return_default_model_name = defaultModelName_;

  if (model_name.empty())  // If no model name, use the default model.
  {
    if (isModelRequired())
    {
      std::ostringstream oss;
      oss << className_ + "::addInstance A " << getName() << " instance must reference a model";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
    }
    else
    {
      model_name = defaultModelName_;
    }
  }

  if (!findModel(model_name)) // Instance references a non-existant model
  {
    std::ostringstream oss;

    oss << className_ << "::addInstance could not find model " << model_name << " which is referenced by instance " << instance_block.getName() + ".";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str() );
  }

  Model &model = *modelMap_[model_name];
  Instance *instance = new Instance(instance_block, model, matrix_load_data, solver_state, extern_data, device_options);

  insertLocalEntity(instance_block.getName(), instance);

  instance->getModel().getInstanceVector().push_back(instance);

  instanceVector_.push_back(instance);

  incrementInstanceCount();

#ifdef Xyce_SIZEOF
  int size = sizeof(*instance);
  std::cout << "Size of device instance " << instance_block.name << " = " << size << std::endl;
#endif

  return instance;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::addModel
// Purpose       : This function adds a model to the list of device models.
//
// Special Notes : This is the  "public" version of the function, meaning that
//                 it returns a bool to indicate success or failure (rather
//                 than a pointer to the model).  Also, it checks to see if the
//                 model being added is already part of the list.  If the model
//                 already exists, it generates an error.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
DeviceModel *DeviceTemplate <Model, Instance>::addModel(const ModelBlock & model_block)
{
  // Check to make sure that there doesn't already exist a model of this name.
  if (findModel(model_block.name))
  {
    std::ostringstream oss;
    oss << className_ << "::addModel:\n\tAttempted to add model that already exists.  "
        << model_block.name << "\nIgnoring all but the first definition.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, oss.str());

    return 0;
  }

  Model *model = new Model(model_block, getSolverState(), getDeviceOptions());

  // Add to the model list.
  insertModel(model_block.name, model);
  modelMap_[model_block.name] = model;

  // Add model to the local processor entity list
  insertLocalEntity(model_block.name, model);

#ifdef Xyce_SIZEOF
  int size = sizeof(*model);
  std::cout << "Size of device model " << model_block.name << " = " << size << std::endl;
#endif

  return model;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate:printOutModels
// Purpose       : debugging tool
// Special Notes :
// Scope         : public
// Creator       : Eric                                  Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/03/06
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
std::ostream &DeviceTemplate<Model, Instance>::printOutModels(std::ostream &os) const
{
  const std::string dashedline("-----------------------------------------------------------------------------");

  os << std::endl
            << std::endl << dashedline << std::endl
            << "Number of "+ className_ + " models: " << modelMap_.size() << std::endl;

  int i = 0;
  for (typename ModelMap::const_iterator it_model = modelMap_.begin(); it_model != modelMap_.end(); ++it_model, ++i)
  {
    const Model &model = *(*it_model).second;
    os << i << ": name = " << model.getName() << " type = " << model.getType() << " level = " << model.getLevel() << std::endl;

    model.printOutInstances(os);
  }
  os << dashedline << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate:printDotOpOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric                                  Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/4/12
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
std::ostream &DeviceTemplate<Model, Instance>::printDotOpOutput (std::ostream &os) const
{
  static const int stdWidth = 15;

// Output models
  std::vector<const Model *> model_list;
  for (typename ModelMap::const_iterator it_model = modelMap_.begin(); it_model != modelMap_.end(); ++it_model)
  {
    const Model *model = (*it_model).second;
    if ( !model->getInstanceVector().empty ())
    {
      model_list.push_back(model);
    }
  }

  std::sort(model_list.begin(), model_list.end(), DeviceEntityCmp());

  ExtendedString strippedDefaultModelName = getName();

  os << std::endl
     << "Number of "+ strippedDefaultModelName + " models: " << model_list.size() << std::endl
     << std::setw(stdWidth) << "name";

  for (typename std::vector<const Model *>::const_iterator it_model = model_list.begin(); it_model != model_list.end(); ++it_model)
  {
    os << std::setw(stdWidth) << (*it_model)->getName() ;
  }
  os << std::endl;

  if (!model_list.empty())
  {
    os << std::setw(stdWidth) << "type";
    for (typename std::vector<const Model *>::const_iterator it_model = model_list.begin(); it_model != model_list.end(); ++it_model)
    {
      os << std::setw(stdWidth) << (*it_model)->getType() ;
    }
    os << std::endl;

    os << std::setw(stdWidth) << "level";
    for (typename std::vector<const Model *>::const_iterator it_model = model_list.begin(); it_model != model_list.end(); ++it_model)
    {
      os << std::setw(stdWidth) << (*it_model)->getLevel() ;
    }
    os << std::endl;

    const DeviceEntity::ParameterMap &parMap = *model_list.front()->getPMap();
    for (DeviceEntity::ParameterMap::const_iterator it_parameter = parMap.begin(); it_parameter != parMap.end(); ++it_parameter)
    {
      for (typename std::vector<const Model *>::const_iterator it_model = model_list.begin(); it_model != model_list.end(); ++it_model)
      {
        ExtendedString parName = (*it_parameter).first;
        parName.toLower ();
        if (it_model == model_list.begin())
        {
          os << std::setw(stdWidth) << parName;
        }

        os << std::setw(stdWidth);
        (*it_model)->printFormattedOutputParam(os, parName);
      }
      os << std::endl;
    }
    os << std::endl;

// Output instances
    int  numInstances = instanceVector_.size();
    os << "Number of " << strippedDefaultModelName << " instances: " << numInstances << std::endl;

    os << std::setw(stdWidth) << "name" ;
    for (typename std::vector<const Model *>::const_iterator it_model = model_list.begin(); it_model != model_list.end(); ++it_model)
    {
      InstanceVector instance_list = (*it_model)->getInstanceVector();
      std::sort(instance_list.begin(), instance_list.end(), DeviceEntityCmp());

      for (typename InstanceVector::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance)
      {
        os << std::setw(stdWidth) << (*it_instance)->getName();
      }
    }
    os << std::endl;

    const Instance &first_instance = *model_list.front()->getInstanceVector().front();

    // // output model, if it was required:
    // if (first_instance.getModelRequired())
    // {
      os << std::setw(stdWidth) << "model";
      for (typename std::vector<const Model *>::const_iterator it_model = model_list.begin(); it_model != model_list.end(); ++it_model)
      {
        const Model &model = *(*it_model);
        InstanceVector instance_list = (*it_model)->getInstanceVector();

        if ( !instance_list.empty () ) {
          std::sort(instance_list.begin(), instance_list.end(), DeviceEntityCmp());

          for (typename InstanceVector::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance)
          {
            ExtendedString modelName = model.getName();
            modelName.toLower();
            os << std::setw(stdWidth) << modelName;
          }
        }
      }

      os << std::endl;
//    }

    // output instance parameters:
    const DeviceEntity::ParameterMap &parameter_map = (*first_instance.getPMap());
    for (DeviceEntity::ParameterMap::const_iterator it_parameter = parameter_map.begin() ; it_parameter != parameter_map.end(); ++it_parameter)
    {
      ExtendedString parName = it_parameter->first;
      parName.toLower ();
      os << std::setw(stdWidth) << parName;

      for (typename std::vector<const Model *>::const_iterator it_model = model_list.begin(); it_model != model_list.end(); ++it_model)
      {
        const Model &model = *(*it_model);
        InstanceVector instance_list = model.getInstanceVector();

        if ( !instance_list.empty () ) {
          std::sort(instance_list.begin(), instance_list.end(), DeviceEntityCmp());

          for (typename InstanceVector::const_iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance)
          {
            os << std::setw(stdWidth);
            (*it_instance)->printFormattedOutputParam(os, parName);
          }
        }
      }
      os << std::endl;
    }
    os << std::endl;
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::deleteInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/03/06
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::deleteInstance (const std::string & instance_name)
{
  bool bsuccess = true;

  // loop over models:
  typename ModelMap::iterator firstM = modelMap_.begin();
  typename ModelMap::iterator lastM  = modelMap_.end();
  for (typename ModelMap::iterator iterM = firstM; iterM != lastM; ++iterM)
  {
    // loop over instances:
    typename InstanceVector::iterator firstI = (*iterM).second->getInstanceVector().begin();
    typename InstanceVector::iterator lastI = (*iterM).second->getInstanceVector().end();

    for (typename InstanceVector::iterator iterI=firstI; iterI!=lastI; ++iterI)
    {
      ExtendedString nameArg(instance_name);
      nameArg.toUpper ();
      ExtendedString nameIter((*iterI)->getName());
      nameIter.toUpper ();
      if (nameArg == nameIter)
      {
        bsuccess = true;
        delete (*iterI);
        (*iterM).second->getInstanceVector().erase(iterI);
        return bsuccess;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::getBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
// Special Notes : For most devices, this function won't get called.  The
//                 device manager only calls getBreakPoints for certain
//                 pre-selected devices.  Ultimately, this function probably
//                 isn't neccessary, as the device manager usually accesses
//                 device instances directly, rather than going through a
//                 device class middleman.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/05/06
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::getBreakPoints( std::vector<N_UTL_BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  const std::string dashedline("-----------------------------------------------------------------------------");
  const std::string dashedline2("---------------------");
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << std::endl << dashedline << std::endl;
    std::cout << className_ + "::getBreakPoints:\n";
  }
#endif

  // loop over models:
  typename ModelMap::iterator firstM = modelMap_.begin();
  typename ModelMap::iterator lastM  = modelMap_.end();
  for (typename ModelMap::iterator iterM = firstM; iterM != lastM; ++iterM)
  {
    // loop over instances:
    typename InstanceVector::iterator firstI = (*iterM).second->getInstanceVector().begin();
    typename InstanceVector::iterator lastI  = (*iterM).second->getInstanceVector().end();
    for (typename InstanceVector::iterator iterI = firstI; iterI != lastI; ++iterI)
    {
      bool tmpBool = (*iterI)->getInstanceBreakPoints(breakPointTimes);
      bsuccess = bsuccess && tmpBool;
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << dashedline2;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::updateSources
//
// Purpose       :
//
// Special Notes : Just like for the getBreakPoints function, for most devices,
//                 this function won't get called.  The
//                 device manager only calls updateSources for certain
//                 pre-selected devices.  Ultimately, this function probably
//                 isn't neccessary, as the device manager usually accesses
//                 device instances directly, rather than going through a
//                 device class middleman.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/05/06
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::updateSources()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  const std::string dashedline("----------------------------------------------------------------------------");
  const std::string dashedline2("---------------------");
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << std::endl;
    std::cout << dashedline << std::endl;
    std::cout << className_ + "::updateSources\n";
  }
#endif

  // first loop over the models:
  typename ModelMap::iterator firstM = modelMap_.begin();
  typename ModelMap::iterator lastM = modelMap_.end();
  for (typename ModelMap::iterator iterM=firstM; iterM!=lastM; ++iterM)
  {
    // loop over the instances for this model.
    typename InstanceVector::iterator firstI = (*iterM).second->getInstanceVector().begin();
    typename InstanceVector::iterator lastI   = (*iterM).second->getInstanceVector().end();
    for (typename InstanceVector::iterator iterI=firstI; iterI!=lastI; ++iterI)
    {
#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        std::cout << dashedline2 <<std::endl;
        std::cout << "  name = " << (*iterI)->getName() <<std::endl;
      }
#endif
      bool tmpBool = (*iterI)->updateSource();
      bsuccess = bsuccess && tmpBool;
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    std::cout << std::endl;
    std::cout << dashedline << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->updatePrimaryState();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::updateSecondaryState (double * staDerivVec, double * stoVec)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->updateSecondaryState ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::loadDAEVectors(double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->loadDAEFVector();
    bsuccess = bsuccess && tmpBool;
    tmpBool = (*it)->loadDAEQVector();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  bool bsuccess = true;
  for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
  {
    bool tmpBool = (*it)->loadDAEdFdx();
    bsuccess = bsuccess && tmpBool;
    tmpBool = (*it)->loadDAEdQdx();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}


//----------------------------------------------------------------------------
// Function       : DeviceTemplate::getInstances
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date  : 08/27/2008
//----------------------------------------------------------------------------
template<typename Model, typename Instance>
void DeviceTemplate<Model, Instance>::getInstances(std::vector<DeviceInstance *> & instanceVector) const
{
  instanceVector.clear();

  // loop over models:
  for (typename ModelMap::const_iterator iterM = modelMap_.begin(); iterM != modelMap_.end(); ++iterM)
  {
    const InstanceVector &instances = (*iterM).second->getInstanceVector();

    // loop over instances:
    for (typename InstanceVector::const_iterator iterI = instances.begin(); iterI != instances.end(); ++iterI)
    {
      instanceVector.push_back(*iterI);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceTemplate::getDeviceNames
//
// Purpose       : Returns a vector of names of all devices of this type
//
// Special Notes : Just like for the getBreakPoints and updateSources function,
//                 for most devices, this function won't get called.  It is
//                 used by special API calls to interface with external codes,
//                 which only happens for very limited device types.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/18/08
//-----------------------------------------------------------------------------
template<typename Model, typename Instance>
bool DeviceTemplate<Model, Instance>::getDeviceNames(std::vector<std::string> &nameVector) const
{
  nameVector.clear();

  for (typename ModelMap::const_iterator iterM = modelMap_.begin(); iterM != modelMap_.end(); ++iterM)
  {
    const InstanceVector &instances = (*iterM).second->getInstanceVector();

    // loop over instances:
    for (typename InstanceVector::const_iterator iterI = instances.begin(); iterI != instances.end(); ++iterI)
    {
      nameVector.push_back((*iterI)->getName());
    }
  }

  return true;
}

} // namespace Device
} // namespace Xyce

#endif
