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
// Revision Number: $Revision: 1.116.2.2 $
//
// Revision Date  : $Date: 2014/03/03 18:29:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_Device_h
#define Xyce_N_DEV_Device_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>

class N_LAS_Matrix;

namespace Xyce {
namespace Device {

struct DeviceModelOp: public std::unary_function<DeviceModel *, bool>
{
  virtual ~DeviceModelOp()
  {}

  virtual bool operator()(DeviceModel *model) = 0;
};

struct DeviceInstanceOp: public std::unary_function<DeviceInstance *, bool>
{

  virtual ~DeviceInstanceOp()
  {}

  virtual bool operator()(DeviceInstance *instance) = 0;
};

///
///  The Device class is an interface for device implementations.
/// 
///  In general, DeviceMaster is the only class that actually inherits
///  from Device.  Most devices either use the DeviceMaster class, but
///  some derive from DeviceMaster in what is known as a Master class.
/// 
///  The interfaces are unfortunately used for but derived
///  implementation as well as usage.  This should be changed at some
///  point, but this is the currently basic design of most interface
///  classes.
/// 
///  @author Eric Keiter, SNL, Parallel Computational Sciences
///  @date   3/16/00
class Device
{
public:
  Device()
  {}

  virtual ~Device()
  {}

private:
  Device(const Device &);                     ///< No copying
  Device &operator=(const Device &);          ///< No assignment

public:
  ///
  ///  Returns true if the device is linear
  /// 
  ///  @return  true if the device is linear.
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:36:27 2014
  virtual bool isLinearDevice() const = 0;

  ///
  ///  Returns true is the device is a PDE device
  /// 
  ///  @return true is the device is a PDE device.
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:37:02 2014
  virtual bool isPDEDevice() const = 0;

  ///
  ///  Returns the name given to the device
  /// 
  ///  @return const reference to the device name.
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:37:34 2014
  virtual const std::string &getName() const = 0;

  ///
  ///  Returns the name of the default model that would to used for this device
  /// 
  ///  @return const reference to the name of the default model
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:38:02 2014
  virtual const std::string &getDefaultModelName() const = 0;

  ///
  ///  Returns the device entity with the specified name
  /// 
  ///  @param entity_name       const reference to the name of the entity
  /// 
  ///  @return pointer to the device entity
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:39:23 2014
  virtual DeviceEntity *findEntity(const std::string &entity_name) = 0;

  ///
  ///  Returns the device entity with the specified name
  /// 
  ///  @param entity_name       const reference name of the entity
  /// 
  ///  @return const pointer to the device entity
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:39:23 2014
  virtual const DeviceEntity *findEntity(const std::string &entity_name) const = 0;

  ///
  ///  Creates a device model and adds it to the device's list of models
  /// 
  /// 
  ///  @param model_block       const reference to the model block describing the model to create
  ///  @param factory_block     const reference to the factory data needed to create the model
  /// 
  ///  @return pointer to the newly created device model
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:40:30 2014
  virtual DeviceModel *addModel(const ModelBlock &model_block, const FactoryBlock &factory_block) = 0;

  ///
  ///  Creates a device instance and adds to the device model's instance list
  /// 
  /// 
  ///  @param instance_block    const reference to the model block describing the instance to create
  ///  @param factory_block     const reference to the factory data needed to create the instance
  /// 
  ///  @return pointer to the newly creates device instance
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:43:08 2014
  virtual DeviceInstance *addInstance(const InstanceBlock &instance_block, const FactoryBlock &factory_block) = 0;

  ///
  ///  Updates the devices source information
  /// 
  ///  This function is called by the analysis subsystem when it is time to update the device source information.
  /// 
  ///  @return true if the update was successful
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:44:34 2014
  virtual bool updateSources() 
  {
    return true;
  }

  ///
  ///  Updates the devices state information
  /// 
  ///  This function is called by the analysis subsystem when it is time to update the device state information.
  /// 
  ///  @return true if the update was successful
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:44:34 2014
  virtual bool updateState(double * solVec, double * staVec, double * stoVec) 
  {
    return true;
  }

  ///
  ///  Updates the devices secondary state information
  /// 
  ///  This function is called by the analysis subsystem when it is time to update the device secondary state information.
  /// 
  ///  @return true if the update was successful
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:44:34 2014
  virtual bool updateSecondaryState(double * staDerivVec, double * stoVec) 
  {
    return true;
  }

  ///
  ///  Populates the device's ExternData object with these pointers
  /// 
  ///  THIS FUNCTION MUST BE CALLED PRIOR TO CALLING loadDAEMatrices.
  /// 
  ///  @param solVec            pointer to the analysis solution vector for this device
  ///  @param fVec              pointer to the analysis f vector for this device
  ///  @param qVec              pointer to the analysis q vector for this device
  ///  @param storeLeadF        pointer to the analysis lead store f vector for this device
  ///  @param storeLeadQ        poitner to the analysis load store q vector for this device
  /// 
  ///  @return true if the update was successful
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 16:47:20 2014
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ) 
  {
    return true;
  }

  ///
  ///  Populates the device's Jacobian object with these pointers
  /// 
  ///  THIS FUNCTION MUST BE CALLED AFTER CALLING loadDAEVectors.
  /// 
  ///  @param dFdx      pointer to the analysis dFdx matrix
  ///  @param dQdx      pointer to the analysis dQdx matrix
  /// 
  ///  @return
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Wed Jan 29 17:06:52 2014
  virtual bool loadDAEMatrices(N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx) 
  {
    return true;
  }

  // ///
  // ///  Delete the specified instance from the device
  // /// 
  // ///  TODO: Check this implementation in DeviceMaster...does it really work?
  // /// 
  // ///  @param instance_name     const reference to the instance name of the device instance to be deleted
  // /// 
  // ///  @return true if the delete was successful
  // /// 
  // ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  // ///  @date   Wed Jan 29 17:08:58 2014
  // virtual bool deleteInstance(DeviceInstance *instance) = 0;

  ///
  ///  Executes op on each DeviceModel pointer of the device
  /// 
  ///  To use this function, create a class which derives from DeviceModelOp.  Then implement the
  ///  operator()(DeviceModel/// ).  Sample classes are defined in N_DEV_Algorithm.C.
  /// 
  ///  @param op        reference to operator functor
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Mon Feb  3 10:36:04 2014
  virtual void forEachModel(DeviceModelOp &op) const = 0;

  ///
  ///  Executes op on each DeviceInstance pointer of the device
  /// 
  ///  To use this function, create a class which derives from DeviceInstanceOp.  Then implement the
  ///  operator()(DeviceInstance *).  Sample classes are defined in N_DEV_Algorithm.C.
  /// 
  ///  @param op        reference to operator functor
  /// 
  ///  @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  ///  @date   Mon Feb  3 10:36:04 2014
  virtual void forEachInstance(DeviceInstanceOp &op) const = 0;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Device N_DEV_Device;

#endif
