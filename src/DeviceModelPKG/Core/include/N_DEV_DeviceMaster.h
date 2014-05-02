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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DeviceMaster.h,v $
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
// Revision Number: $Revision: 1.3.2.2 $
//
// Revision Date  : $Date: 2014/03/03 18:29:27 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Device_Template_h
#define Xyce_N_DEV_Device_Template_h

#include <map>
#include <string>
#include <vector>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>

namespace Xyce {
namespace Device {

void duplicate_entity_warning(const DeviceEntity &entity);
void instance_must_reference_model_error(const Device &device, const std::string &name, const std::string &netlist_path, int line_number);
void could_not_find_model(const Device &device, const std::string &model_name, const std::string &instance_name, const std::string &netlist_path, int line_number);
void duplicate_model_warning(const Device &device, const std::string &model_name, const std::string &netlist_path, int line_number);

//-----------------------------------------------------------------------------
// Class         : DeviceMaster
//
// Purpose       : This class contains a lot of boilerplate functions, such
//                 as addInstance, factory_block, which are needed by every device.  However,
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

/**
 * DeviceMaster instantiates a device as described by the device traits T.
 *
 * The Device Traits are described in in the device configuration.
 *
 * @see Xyce::Device::Configuration
 *
 * @param T     device traits class
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Feb  4 10:33:44 2014
 */
template<class T>
class DeviceMaster : public Device
{
public:
  typedef typename T::ModelType ModelType;            ///< Make the model begin defined available
  typedef typename T::InstanceType InstanceType;      ///< Make the instance being define available

protected:
  typedef std::vector<InstanceType *> InstanceVector;
  typedef std::map<std::string, ModelType *, LessNoCase> ModelMap ;
  typedef std::map<std::string, DeviceEntity *, LessNoCase> EntityMap;

public:
  /**
   *  Constructs a device
   *
   * When a device is constructed, a default model is created.
   *
   * @param configuration             const reference to device configuration
   * @param factory_block             const reference to the factory provided parameters
   * @param solver_state              const reference to the solver state to use for the device
   * @param device_options            const reference to the device options to use for the device
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:20:40 2014
   */
  DeviceMaster(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : Device(),
      deviceName_(T::name()),
      defaultModelName_(std::string(T::deviceTypeName()) + " (" + T::name() + ")"),
      configuration_(configuration),
      deviceOptions_(device_options),
      solverState_(solver_state),
      modelMap_(),
      instanceVector_(),
      entityMap_(),
      defaultModel_(new ModelType(configuration_, ModelBlock(defaultModelName_, ""), factory_block))
  {
    // add default model to the model list.
    modelMap_[defaultModelName_] = defaultModel_;
  }

  /**
   * Constructs a device
   *
   * When a device is constructed, a default model is created.
   *
   * @param configuration             const reference to the device configuration
   * @param model_type_name           const reference to the model type name inserted into the ModelBlock when creating the default model
   * @param factory_block             const reference to the factory provided parameters
   * @param solver_state              const reference to the solver state to use for the device
   * @param device_options            const reference to the device options to use for the device
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:22:41 2014
   */
  DeviceMaster(
     const std::string &       model_type_name,
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : Device(),
      deviceName_(T::name()),
      defaultModelName_(std::string(T::deviceTypeName()) + " (" + T::name() + ")"),
      configuration_(configuration),
      deviceOptions_(device_options),
      solverState_(solver_state),
      modelMap_(),
      instanceVector_(),
    entityMap_(),
    defaultModel_(new ModelType(configuration_, ModelBlock(defaultModelName_, model_type_name), factory_block))
  {
    // add default model to the model list.
    modelMap_[defaultModelName_] = defaultModel_;
  }

  /**
   * Destroys the device
   *
   * Delete the models created by this device.  The model will handle deleting the instances.
   *
   * @return
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:23:31 2014
   */
  virtual ~DeviceMaster()
  {
    for (typename ModelMap::iterator model_it = modelMap_.begin(); model_it != modelMap_.end(); ++model_it)
    {
      delete (*model_it).second;
    }
  }

private:
  DeviceMaster(const DeviceMaster &right);                ///< No copying
  DeviceMaster(const Device &);                             ///< Eh?
  DeviceMaster &operator=(const DeviceMaster &right);     ///< No assignments

public:
  /**
   * Returns the name of this device
   *
   * @return const reference to the device's name
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:25:37 2014
   */
  virtual const std::string &getName() const /* override */ 
  {
    return deviceName_;
  }

  /**
   * Returns the default model name to use if the instance being created does not specify one
   *
   * @return const reference to the default model name
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:26:10 2014
   */
  virtual const std::string &getDefaultModelName() const /* override */ 
  {
    return defaultModelName_;
  }
  /**
   * Returns a pointer to the model or instance entity with the specified name
   *
   * @param entity_name       const reference to the entity name
   *
   * @return pointer to the entity or 0 is not found
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:36:44 2014
   */
  virtual DeviceEntity *findEntity(const std::string &entity_name) /* override */ 
  {
    EntityMap::iterator it = entityMap_.find(entity_name);
    if (it != entityMap_.end())
      return (*it).second;

    return 0;
  }

  /**
   * Returns a pointer to the model or instance entity with the specified name
   *
   * @param entity_name       const reference to the entity name
   *
   * @return const pointer to the entity or 0 is not found
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:36:44 2014
   */
  virtual const DeviceEntity *findEntity(const std::string &entity_name) const /* override */ 
  {
    EntityMap::const_iterator it = entityMap_.find(entity_name);
    if (it != entityMap_.end())
      return (*it).second;

    return 0;
  }

  /**
   * Returns true if this device is a linear device
   *
   * @return the device is linear value provided via the device trait.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:51:52 2014
   */
  virtual bool isLinearDevice() const /* override */ 
  {
    return T::isLinearDevice();
  }

  /**
   * Returns true if this device is a PDE device
   *
   * @return the device is PDE value provided via the device trait.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:53:10 2014
   */
  virtual bool isPDEDevice() const /* override */ 
  {
    return T::isPDEDevice();
  }

  /**
   * Executes operator op, passing its DeviceModel pointer, for each device model
   *
   * Since the model information is managed in a map of ModelType pointers, it is not possible to return the map as an
   * association to DeviceModel pointers.  This function iterates over the models and passes it the DeviceModel
   * pointer for each.
   *
   * @param op        functor taking a DeviceModel pointer as an argument
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:53:45 2014
   */
  virtual void forEachModel(DeviceModelOp &op) const /* override */ 
  {
    for (typename ModelMap::const_iterator it = modelMap_.begin(); it != modelMap_.end(); ++it)
      op((*it).second);
  }

  /**
   * Executes operator op, passing its DeviceInstance pointer, for each device instance
   *
   * Since the instance information is managed in a map of InstanceType pointers, it is not possible to return the map as an
   * association to DeviceInstance pointers.  This function iterates over the instances and passes it the DeviceInstance
   * pointer for each.
   *
   * @param op        functor taking a DeviceInstance pointer as an argument
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:53:45 2014
   */
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */ 
  {
    for (typename InstanceVector::const_iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
      op(*it);
  }

  // virtual bool deleteInstance(DeviceInstance *instance) /* override */;

  virtual DeviceModel *addModel(const ModelBlock & MB, const FactoryBlock &factory_block) /* override */;
  virtual DeviceInstance *addInstance(const InstanceBlock &instance_block, const FactoryBlock &factory_block); /* override */

  virtual bool updateSources() /* override */;
  virtual bool updateState (double * solVec, double * staVec, double * stoVec)/* override */;
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) /* override */;
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ) /* override */;
  virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx) /* override */;

protected:
  /**
   * Returns the solver state given during device construction
   *
   * @return const reference to the solver state
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:27:14 2014
   */
  const SolverState &getSolverState() const 
  {
    return solverState_;
  }

  /**
   * Returns the device options given during device construction
   *
   * @return const reference to the device options
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:27:48 2014
   */
  const DeviceOptions &getDeviceOptions() const 
  {
    return deviceOptions_;
  }

  /**
   * Returns an iterator to the beginning of the vector of all instances created for this device
   *
   * While a device instance is created, the device model owns the pointer to the device.  The instanceVector_ gets a
   * copy so that all instances of this device may be iterated over without needing to go through the model.
   *
   * @return iterator to the beginning of the device instance vector
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:28:29 2014
   */
  typename InstanceVector::const_iterator getInstanceBegin() const 
  {
    return instanceVector_.begin();
  }

  /**
   * Returns an iterator to the ending of the vector of all instances created for this device
   *
   * While a device instance is created, the device model owns the pointer to the device.  The instanceVector_ gets a
   * copy so that all instances of this device may be iterated over without needing to go through the model.
   *
   * @return iterator to the ending of the device instance vector
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:28:29 2014
   */
  typename InstanceVector::const_iterator getInstanceEnd() const 
  {
    return instanceVector_.end();
  }

  /**
   * Returns true if the model name must be specified for each instance
   *
   * @return the device model required value provided via the device trait.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:32:10 2014
   */
  bool isModelRequired() const 
  {
    return T::modelRequired();
  }

private:
  /**
   * Adds an entity to the mapping of model and instance name to its entity
   *
   * Note that currently if the name is duplicated, the previous value is overwritten.
   *
   * @param name      const reference to the model or instance name
   * @param entity    pointer to the model or instance
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Tue Feb  4 10:35:00 2014
   */
  void addEntity(const std::string &name, DeviceEntity *entity) 
  {
    std::pair<EntityMap::iterator, bool> result = entityMap_.insert(EntityMap::value_type(name, entity));
    if (!result.second)
      duplicate_entity_warning(*entity);
  }

  virtual bool getBreakPoints(std::vector<Util::BreakPoint> & breakPointTimes);

private:
  const std::string           deviceName_;
  const std::string           defaultModelName_;
  const Configuration &       configuration_;
  const SolverState &         solverState_;
  const DeviceOptions &       deviceOptions_;
  ModelMap                    modelMap_;
  InstanceVector              instanceVector_;
  EntityMap                   entityMap_;
  ModelType * const           defaultModel_;
};

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::addInstance
// Purpose       : This function adds an instance to the list of instances.
//
// Special Notes : All device instances are contained in one of many lists.
//                 Each model owns exactly one list of instances.  If
//                 an instance is added which does not reference a model, or
//                 that references a model that does not exist, factory_block, then a new
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
template<class T>
DeviceInstance *
DeviceMaster<T>::addInstance(
   const InstanceBlock &         instance_block,
   const FactoryBlock &          factory_block)
{
  std::string model_name = instance_block.getModelName();

  if (model_name.empty())  // If no model name, use the default model.
  {
    if (isModelRequired())
    {
      instance_must_reference_model_error(*this, model_name, instance_block.netlistFileName_, instance_block.lineNumber_);
    }
    else
    {
      model_name = defaultModelName_;
    }
  }

  typename ModelMap::iterator it = modelMap_.find(model_name);
  if (it == modelMap_.end()) 
  {
    could_not_find_model(*this, model_name, instance_block.getName(), instance_block.netlistFileName_, instance_block.lineNumber_);
    return 0;
  }

  ModelType &model = *(*it).second;

  InstanceType *instance = new InstanceType(configuration_, instance_block, model, factory_block);
  instance->setDefaultParamName(T::instanceDefaultParameter());

  addEntity(instance_block.getName(), instance);

  model.addInstance(instance);

  instanceVector_.push_back(instance);

  return instance;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::addModel
// Purpose       : This function adds a model to the list of device models.
//
// Special Notes : This is the  "public" version of the function, meaning that
//                 it returns a bool to indicate success or failure (rather
//                 than a pointer to the model).  Also, factory_block, it checks to see if the
//                 model being added is already part of the list.  If the model
//                 already exists, it generates an error.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
template<class T>
DeviceModel *DeviceMaster<T>::addModel(
   const ModelBlock &    model_block,
   const FactoryBlock &  factory_block)
{
  // Check to make sure that there doesn't already exist a model of this name.
  if (modelMap_.find(model_block.name) != modelMap_.end()) 
  {
    duplicate_model_warning(*this, model_block.name, model_block.netlistFileName_, model_block.lineNumber_);

    return 0;
  }

  ModelType *model = new ModelType(configuration_, model_block, factory_block);

  // Add to the model list.
  modelMap_[model_block.name] = model;

  // Add model to the local processor entity list
  addEntity(model_block.name, model);

  return model;
}

// //-----------------------------------------------------------------------------
// // Function      : DeviceMaster::deleteInstance
// // Purpose       :
// // Special Notes :
// // Scope         : public
// // Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// // Creation Date : 2/03/06
// //-----------------------------------------------------------------------------
// template<class T>
// bool DeviceMaster<T>::deleteInstance(const std::string & instance_name)
// {
//   for (typename ModelMap::iterator it_model = modelMap_.begin(); it_model != modelMap_.end(); ++it_model)
//   {
//     InstanceVector &instance_list = (*it_model).second->getInstanceVector();

//     for (typename InstanceVector::iterator it_instance = instance_list.begin(); it_instance != instance_list.end(); ++it_instance)
//     {
//       if (equal_nocase(instance_name, (*it_instance)->getName()))
//       {
//         delete *it_instance;
//         instance_list.erase(it_instance);

//         return true;
//       }
//     }
//   }

//   return false;
// }

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::getBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
// Special Notes : For most devices, this function won't get called.  The
//                 device manager only calls getBreakPoints for certain
//                 pre-selected devices.  Ultimately, factory_block, this function probably
//                 isn't neccessary, as the device manager usually accesses
//                 device instances directly, rather than going through a
//                 device class middleman.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/05/06
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::getBreakPoints( std::vector<Util::BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  for (typename InstanceVector::iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
    bsuccess = (*it)->getInstanceBreakPoints(breakPointTimes) && bsuccess;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::updateSources
//
// Purpose       :
//
// Special Notes : Just like for the getBreakPoints function, for most devices, factory_block,
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
template<class T>
bool DeviceMaster<T>::updateSources()
{
  bool bsuccess = true;

  for (typename InstanceVector::iterator it = instanceVector_.begin(); it != instanceVector_.end(); ++it)
    bsuccess = (*it)->updateSource() && bsuccess;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMaster::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::updateState (double * solVec, double * staVec, double * stoVec)
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
// Function      : DeviceMaster::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::updateSecondaryState (double * staDerivVec, double * stoVec)
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
// Function      : DeviceMaster::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::loadDAEVectors(double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ)
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
// Function      : DeviceMaster::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, factory_block, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
template<class T>
bool DeviceMaster<T>::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
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

} // namespace Device
} // namespace Xyce

#endif
