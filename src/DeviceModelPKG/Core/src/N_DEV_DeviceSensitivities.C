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
// Filename       : $RCSfile: N_DEV_DeviceSensitivities.C,v $
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
// Revision Number: $Revision: 1.50.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <algorithm>
#include <vector>
#include <list>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceSensitivities.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>

#include <N_UTL_Misc.h>
#include <N_UTL_OptionBlock.h>

#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_NLS_Manager.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::DeviceSensitivities
// Purpose       : constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
DeviceSensitivities::DeviceSensitivities(
  DeviceOptions & do1,
  ExternData & extData1,
  SolverState & solState1,
  N_LAS_System & lasSys1)

 : devOptions_(do1),
   devMgr_( *(extData1.devMgrPtr)),
   extData_(extData1),
   solState_(solState1),
   entityMapDone_(false),
   numSensParams_(0),
   lasSys_(lasSys1),
   nlsMgr_(devMgr_.getNlsMgrPtr())
{

}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::DeviceSensitivities
// Purpose       : copy constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
DeviceSensitivities::DeviceSensitivities (const DeviceSensitivities &right)
  :devOptions_(right.devOptions_),
  devMgr_(right.devMgr_),
   extData_(right.extData_),
   solState_(right.solState_),
   entityMapDone_(right.entityMapDone_),
   numSensParams_(right.numSensParams_),
   lasSys_(right.lasSys_),
   nlsMgr_(right.nlsMgr_)

{

}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::~DeviceSensitivities
// Purpose       : destructor
// Special Notes : De-allocates all the devices pointed  to by deviceArray
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
DeviceSensitivities::~DeviceSensitivities()
{

}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivites::registerSensParams
//
// Purpose       : This function takes an option block, which contains a
//                 list of user-defined device parameters, and stores
//                 this list in a map.  The map is the "deviceEntityMap",
//                 which maps a parameter string to a device entity (a
//                 device entity is either an instance or a model).
//
// Special Notes : The map is not complete at the end of this function - it
//                 has not filled in the other side(the dePtr side) of the map.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
bool DeviceSensitivities::registerSensParams
  (const N_UTL_OptionBlock & OB)
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.sensDebugLevel > 0)
  {
    cout << "DeviceSensitivites::registerSensParams called!" <<endl;
  }
#endif

  list<N_UTL_Param>::const_iterator iter = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator end   = OB.getParams().end();

  // set up the  nameDevInstMap_ :
  numSensParams_ = 0;
  for ( ; iter !=  end;  ++iter)
  {
    if ( std::string( iter->uTag() ,0,5) == "PARAM") // this is a vector
    {
      ExtendedString tag(iter->sVal());
      tag.toUpper();

#ifdef Xyce_DEBUG_DEVICE
      cout << "name = " << iter->uTag() << "  tag = " << tag << endl;
#endif

      // set up the initial skeleton of the maps:
      // these will be cleaned up later.
      DeviceEntity *dePtr = NULL;
      nameDevEntityMap_[tag] = dePtr;

      sensParamArray_.push_back(tag);

      ++numSensParams_;
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.sensDebugLevel > 0)
  {
    cout << "numSensParams_ = "<< numSensParams_ << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivites::findDeviceEntity_
//
// Purpose       : This function does a manual search of the various
//                 allocated devices, to find a pointer to the device
//                 instance or model which owns the passed parameter.
//
// Special Notes : This is not a map-based search.  It should only be done
//                 once for any one parameter.
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
const DeviceEntity * DeviceSensitivities::findDeviceEntity_ (const std::string & param) const
{
  string tmpName(param);
  string strippedName;
  string entityName;
  stripParamName (tmpName, strippedName, entityName);
  ExtendedString deviceEntityName(entityName);
  deviceEntityName.toUpper ();

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.sensDebugLevel > 0)
  {
    cout << " deviceEntityName = " << deviceEntityName << endl;
  }
#endif

  // now that we have the name, find the actual entity.
  // loop over all the device types.  query  each one.
  const DeviceEntity * dePtr = 0;

  for (int i = 0; i < ModelType::NUMDEV; ++i)
  {
    if ( devMgr_.getDeviceAllocFlag(i))
    {
      dePtr = devMgr_.getDeviceArray(i).findLocalEntity(deviceEntityName);
    }
    if (dePtr != NULL)
    {
      break;
    }
  }

  return dePtr;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
DeviceEntity * DeviceSensitivities::findDeviceEntity_ (const std::string & param)
{
  string tmpName(param);
  string strippedName;
  string entityName;
  stripParamName (tmpName, strippedName, entityName);
  ExtendedString deviceEntityName(entityName);
  deviceEntityName.toUpper ();

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.sensDebugLevel > 0)
  {
    cout << " deviceEntityName = " << deviceEntityName << endl;
  }
#endif

  // now that we have the name, find the actual entity.
  // loop over all the device types.  query  each one.
  DeviceEntity * dePtr = 0;

  for (int i = 0; i < ModelType::NUMDEV; ++i)
  {
    if ( devMgr_.getDeviceAllocFlag(i))
    {
      dePtr = devMgr_.getDeviceArray(i).findLocalEntity(deviceEntityName);
    }
    if (dePtr != NULL)
    {
      break;
    }
  }

  return dePtr;
}
//-----------------------------------------------------------------------------
// Function      : DeviceSensitivites::setUpDeviceEntityMap_
//
// Purpose       : This function completes the setup of the nameDevEntityMap
//                 that was started in the function registerSensOptions.
//
// Special Notes : This function should be called sometime after all the
//                 device entities (instances and models) have been
//                 allocated.
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
void DeviceSensitivities::setUpDeviceEntityMap_ ()
{

#ifdef Xyce_DEBUG_PARAMETER_HANDLING
  if (devOptions_.sensDebugLevel > 0)
  {
    cout << "\nDeviceSensitivites::setUpDeviceEntityMap_ called!" <<endl;
    cout << endl;
  }
#endif

  if (entityMapDone_) return;

  // Note:
  // device models are looked up based on type and level.
  // device instances are looked up based on names.

  // loop over the device entity map

  map<string, DeviceEntity*>::iterator firstDEM=nameDevEntityMap_.begin();
  map<string, DeviceEntity*>::iterator lastDEM =nameDevEntityMap_.end();
  map<string, DeviceEntity*>::iterator iterDEM;

  for (iterDEM=firstDEM; iterDEM != lastDEM; ++iterDEM)
  {
    // find the last ":" in the string, use this information to obtain
    // the name of the device instance or model.

    string tmpName(iterDEM->first);

#ifdef Xyce_DEBUG_PARAMETER_HANDLING
    if (devOptions_.sensDebugLevel > 0)
    {
      cout << " param name = " << tmpName << endl;
    }
#endif

    DeviceEntity * dePtr = findDeviceEntity_(tmpName);
    bool bsuccess = (dePtr!=0);

    if (bsuccess)
    {
      iterDEM->second = dePtr;
#ifdef Xyce_DEBUG_PARAMETER_HANDLING
      if (devOptions_.sensDebugLevel > 0)
      {
	cout << " Found entity, adding it to the DEM map" << endl;
      }
#endif
    }
  }  // end of loop over the device entity map.

  entityMapDone_ = true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivities::stripParamName
//
// Purpose       : This function changes a parameter name from its "full"
//                 name, which includes the name of the device, to just the
//                 parameter name by itself.  For example, the full name
//                 might be:  R1:R, for the resistor instance R1, with a
//                 resistance of R.  The stripped name would be R, and the
//                 entity name would be R1.
//
//                 ERK. 11/11/03.  Changed to allow for default parameters.
//                 If the user does not include a ":" in the fullName, then
//                 it is assumed that the full name is the entity name, and
//                 that the parameter has not been specified.  Certain devices
//                 have default parameters, so they are used if they exist.
//
// Special Notes : I probably need to move this and some of the other
//                 parameter handling capability into another class, which
//                 is devoted to parameter management.
//
//                 Also, it occurs to me that this function might not work
//                 very well for subcircuits, because of multiple ":".
//                 Revisit later.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/25/03
//-----------------------------------------------------------------------------
bool DeviceSensitivities::stripParamName(
  const string & fullName,
  string & strippedName,
  string & entityName) const
{
  // update the parameter:  first extract its name from the string:
  // the name of the parameter will be everything to
  // the right of the ":".
  string tmpName(fullName);
  string::iterator firstCH = tmpName.begin();
  string::iterator lastCH  = tmpName.end();
  string::iterator iterCH;
  string::iterator iterSave;
  string tmpChar(":");

  // find the instance of ":"
  bool foundFlag = false;
  for (iterCH=firstCH; iterCH!=lastCH ; ++iterCH)
  {
    if ( *iterCH == *(tmpChar.begin()) )
    {
      iterSave = iterCH;
      foundFlag = true;
    }
  }

  // error trap.
#ifdef Xyce_DEBUG_PARAMETER_HANDLING
  if (!foundFlag && devOptions_.sensDebugLevel > 1)
  {
    string msg("DeviceSensitivites::stripParamName ");
    msg += "No : character in " + fullName  + ".";
    msg += "  Assuming default parameter.\n\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg);
  }
#endif

  if (foundFlag)
  {
    // Get the entity name, by erasing everything after the last ":"
    string tmpString (tmpName.begin(),iterSave);
    ExtendedString deviceEntityName(tmpString);
    deviceEntityName.toLower();

    entityName = deviceEntityName;

    // get the parameter (stripped) name:
    ++iterSave;
    strippedName = string(iterSave, lastCH);
  }
  else
  {
    entityName = fullName;
    strippedName = "";
  }

#ifdef Xyce_DEBUG_PARAMETER_HANDLING
  if (devOptions_.sensDebugLevel > 1)
  {
    cout << "fullName = " << fullName;
    cout << "  stripped Name = " << strippedName;
    cout << "  entity Name = " << entityName << endl;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSensitivites::getDeviceEntity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
DeviceEntity * DeviceSensitivities::getDeviceEntity (const string & param)
{
  DeviceEntity * dePtr = 0;

  // Is this named parameter in the map?  If not, add.
  if ( nameDevEntityMap_.find(param) == nameDevEntityMap_.end() )
  {
    nameDevEntityMap_[param] = dePtr;
  }

  // Do we need to find the dePtr manually?  If so, do it.
  dePtr = nameDevEntityMap_[param];
  if (dePtr == 0)
  {
    dePtr = findDeviceEntity_(param);
    nameDevEntityMap_[param] = dePtr;
  }

  return dePtr;
}

} // namespace Device
} // namespace Xyce
