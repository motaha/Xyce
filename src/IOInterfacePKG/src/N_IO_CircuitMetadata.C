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
// Filename      : N_IO_CircuitMetadata.C
//
// Purpose       :
//
// Special Notes :
//
// Creator       :
//
// Creation Date :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.164.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#include <iostream>
#include <string>
#include <vector>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif


// ----------   Xyce Includes   ----------

#include <N_DEV_Param.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitMetadata.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_Factory.h>


namespace Xyce {
namespace IO {

enum Line_Types
{
  _TYPE,
  _LEVEL,
  _NUM_NODES,
  _DEF_PARAM,
  _PAR_VAL,
  _MOD_NAMES
};

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::findDeviceMetadata
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 09/28/2003
//----------------------------------------------------------------------------
DeviceMetadata *CircuitMetadata::findDeviceMetadata(
  const string &        deviceTypeIn,
  int                   level)
{
  if (level == -1)
  {
    level = 1;
  }

  string deviceType = deviceTypeIn;
  if (deviceTypeIn == "K")
  {
    deviceType = "L";
  }

  DeviceMetadataIndexMap::iterator it_device_meta_index = deviceMetadataIndex.find(DeviceLevelKey(deviceType, level));
  if (it_device_meta_index != deviceMetadataIndex.end())
  {
    return &(deviceMetadata_[it_device_meta_index->second]);
  }

  // Handle default model:
  DeviceMetadata *DM = createMetadataEntry (deviceType, level);

  if (Xyce::Device::Tom::exists(Xyce::Device::getXyceModelRegistry(), Xyce::Device::Tom::Key(deviceType, level))) {
    Xyce::Device::ParametricData<void> &map = (*Xyce::Device::Tom::getFunction(Xyce::Device::getXyceModelRegistry(), Xyce::Device::Tom::Key(deviceType, level)))();

    Xyce::Device::populateParams(map.getMap(), DM->modelParameters, DM->modelCompositeParameterMap);
  }

  if (Xyce::Device::Tom::exists(Xyce::Device::getXyceInstanceRegistry(), Xyce::Device::Tom::Key(deviceType, level))) {
    Xyce::Device::ParametricData<void> &data = (*Xyce::Device::Tom::getFunction(Xyce::Device::getXyceInstanceRegistry(), Xyce::Device::Tom::Key(deviceType, level)))();

    Xyce::Device::populateParams(data.getMap(), DM->instanceParameters, DM->instanceCompositeParameterMap);
    DM->numOptionalNodes = data.getNumOptionalNodes();
    DM->numFillNodes = data.getNumFillNodes();
    DM->numNodes = data.getNumNodes();
    DM->modelRequired = data.getModelRequired();
    DM->primaryParameter = data.getPrimaryParameter();
    DM->modelTypes = data.getModelTypes();
  }
  for (std::vector<std::string>::iterator it_device_model_type = DM->modelTypes.begin(); it_device_model_type != DM->modelTypes.end(); ++it_device_model_type)
  {
    deviceMetadataIndex.insert(DeviceMetadataIndexMap::value_type(DeviceLevelKey(*it_device_model_type, level), DeviceLevelKey(deviceType, level)));
  }
  deviceMetadataIndex[DeviceLevelKey(deviceType, level)] = DeviceLevelKey(deviceType, level);
  if (deviceType == "L")
    deviceMetadataIndex[DeviceLevelKey(string("K"), level)] = DeviceLevelKey(deviceType, level);

  return DM;
}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::createMetadataEntry
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Dave Shirley
// Creation Date  : 09/16/05
//----------------------------------------------------------------------------
DeviceMetadata *CircuitMetadata::createMetadataEntry(const string & deviceType, int level)
{
  DeviceMetadata DM(deviceType, level);

  std::pair<DeviceMetadataMap::iterator, bool> mdp = deviceMetadata_.insert(DeviceMetadataMap::value_type(DeviceLevelKey(deviceType, level), DM));

  return &(mdp.first->second);
}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::buildMetadata
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 09/25/2003
//----------------------------------------------------------------------------
void CircuitMetadata::buildMetadata()
{
  // Build source function, and options metadata.

  sourceFunctionMetadata();
  optionsMetadata();
}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::optionsMetadata
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Dave Shirley
// Creation Date  : 05/03/06
//----------------------------------------------------------------------------
void CircuitMetadata::optionsMetadata()
{
  string input_line, parameter, value, optionsPackage;
  vector<N_UTL_Param> optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("DEFAS", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEFAD", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEFL", 1.0e-4));
  optionsParameters.push_back(N_UTL_Param("DEFW", 1.0e-4));
  optionsParameters.push_back(N_UTL_Param("GMIN", 1.0e-12));
  optionsParameters.push_back(N_UTL_Param("GMINSCALAR", 1.0e+10));
  optionsParameters.push_back(N_UTL_Param("GMAX", 1.0e20));
  optionsParameters.push_back(N_UTL_Param("TEMP", 27.0));
  optionsParameters.push_back(N_UTL_Param("TNOM", 27.0));
  optionsParameters.push_back(N_UTL_Param("SCALESRC", 0.0));
  optionsParameters.push_back(N_UTL_Param("NUMJAC", 0));
  optionsParameters.push_back(N_UTL_Param("TESTJAC", 0));
  optionsParameters.push_back(N_UTL_Param("TESTJACSTARTSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("TESTJACSTOPSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("TJRELTOL", 0.01));
  optionsParameters.push_back(N_UTL_Param("TJABSTOL", 1.0e-8));
  optionsParameters.push_back(N_UTL_Param("TJSQRTETA", 1.0e-8));
  optionsParameters.push_back(N_UTL_Param("SENSDP", 1.0e-8));
  optionsParameters.push_back(N_UTL_Param("TESTJACWARN", 0));
  optionsParameters.push_back(N_UTL_Param("TESTJACDEVICENAME", ""));
  optionsParameters.push_back(N_UTL_Param("VOLTLIM", 1));
  optionsParameters.push_back(N_UTL_Param("ICFAC", 10000.0));
  optionsParameters.push_back(N_UTL_Param("LAMBERTW", 0));
  optionsParameters.push_back(N_UTL_Param("MAXTIMESTEP", 1.0e99));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIMESTEP", 65536));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIME", 100.0));
  optionsParameters.push_back(N_UTL_Param("VERBOSELEVEL", 0));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0e-12));
  optionsParameters.push_back(N_UTL_Param("CHGTOL", 1.0e-12));
  optionsParameters.push_back(N_UTL_Param("VDSSCALEMIN", 0.3));
  optionsParameters.push_back(N_UTL_Param("VGSTCONST", 4.5));
  optionsParameters.push_back(N_UTL_Param("LENGTH0", 5.0e-6));
  optionsParameters.push_back(N_UTL_Param("WIDTH0", 200.0e-6));
  optionsParameters.push_back(N_UTL_Param("TOX0", 6.0e-8));
  optionsParameters.push_back(N_UTL_Param("MINRES", 0.0));
  optionsParameters.push_back(N_UTL_Param("MINCAP", 0.0));
  optionsParameters.push_back(N_UTL_Param("SENSDEBUGLEVEL", 0));
  optionsParameters.push_back(N_UTL_Param("NUMGAINSCALEBLOCKS", 1));
  optionsParameters.push_back(N_UTL_Param("STAGGERGAINSCALE", 0));
  optionsParameters.push_back(N_UTL_Param("RANDOMIZEVGSTCONST", 0));
  optionsParameters.push_back(N_UTL_Param("NEWEXCESSPHASE", 1));
  optionsParameters.push_back(N_UTL_Param("EXCESSPHASESCALAR1", 1.0));
  optionsParameters.push_back(N_UTL_Param("EXCESSPHASESCALAR2", 1.0));
  optionsParameters.push_back(N_UTL_Param("RANDOMSEED", 0));
  optionsParameters.push_back(N_UTL_Param("TRYTOCOMPACT", false));
  optionsParameters.push_back(N_UTL_Param("CALCULATEALLLEADCURRENTS", false));
  optionsParameters.push_back(N_UTL_Param("NEWMEYER", false ));
  optionsParameters.push_back(N_UTL_Param("ZERORESISTANCETOL", 1.0e-100 ));
  optionsParameters.push_back(N_UTL_Param("CHECKFORZERORESISTANCE", true ));
  optionsParameters.push_back(N_UTL_Param("DETAILED_DEVICE_COUNTS", false ));
  optionsMetadata_[string("DEVICE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("OBJFUNC", "V+1"));
  optionsParameters.push_back(N_UTL_Param("PARAM", "VECTOR"));
  optionsMetadata_[string("SENS")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 0));
  optionsParameters.push_back(N_UTL_Param("ADJOINT", 0));
  optionsParameters.push_back(N_UTL_Param("DIRECT", 0));
  optionsParameters.push_back(N_UTL_Param("DIFFERENCE", 0));
  optionsParameters.push_back(N_UTL_Param("SQRTETA", 1.0e-8));
  optionsMetadata_[string("SENSITIVITY")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("METHOD", 1));
#ifdef Xyce_DEBUG_ANALYSIS
  optionsParameters.push_back(N_UTL_Param("CONSTSTEP", 0));
#endif
  optionsParameters.push_back(N_UTL_Param("USEDEVICEMAX", 1));
  optionsParameters.push_back(N_UTL_Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(N_UTL_Param("RESTARTSTEPSCALE", .005));
//  optionsParameters.push_back(N_UTL_Param("NLNEARCONV", 1));
  optionsParameters.push_back(N_UTL_Param("NLNEARCONV", 0));
  optionsParameters.push_back(N_UTL_Param("NLSMALLUPDATE", 1));
  optionsParameters.push_back(N_UTL_Param("DOUBLEDCOPSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("FIRSTDCOPSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("LASTDCOPSTEP", 1));
  optionsParameters.push_back(N_UTL_Param("RESETTRANNLS", 1));
  optionsParameters.push_back(N_UTL_Param("BPENABLE", 1));
  optionsParameters.push_back(N_UTL_Param("EXITTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("EXITSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("ERROPTION", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 0));
  optionsParameters.push_back(N_UTL_Param("JACLIMITFLAG", 0));
  optionsParameters.push_back(N_UTL_Param("JACLIMIT", 1.0e17));
  optionsParameters.push_back(N_UTL_Param("DAESTATEDERIV", 0));
  optionsParameters.push_back(N_UTL_Param("TESTFIRSTSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DTMIN", 0.0));
  optionsParameters.push_back(N_UTL_Param("NEWBPSTEPPING", 0));
  optionsParameters.push_back(N_UTL_Param("MINTIMESTEPSBP", 10));
  optionsParameters.push_back(N_UTL_Param("NEWLTE", 0));
  optionsParameters.push_back(N_UTL_Param("MAXORD", 5));
  optionsParameters.push_back(N_UTL_Param("MINORD", 1));
  optionsParameters.push_back(N_UTL_Param("OUTPUTINTERPMPDE", 1));
  optionsParameters.push_back(N_UTL_Param("INTERPOUTPUT", 1));
  optionsParameters.push_back(N_UTL_Param("CONDTEST", 0));
  optionsParameters.push_back(N_UTL_Param("CONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(N_UTL_Param("ISOCONDTEST", 0));
  optionsParameters.push_back(N_UTL_Param("ISOCONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(N_UTL_Param("PASSNLSTALL", false));
  optionsParameters.push_back(N_UTL_Param("NLMIN", 3));
  optionsParameters.push_back(N_UTL_Param("NLMAX", 8));
  optionsParameters.push_back(N_UTL_Param("DELMAX", 1.0e+99));
  optionsParameters.push_back(N_UTL_Param("TIMESTEPSREVERSAL", false));
  optionsParameters.push_back(N_UTL_Param("MINTIMESTEPRECOVERY", 0));
  optionsParameters.push_back(N_UTL_Param("FASTTESTS", false));
  optionsParameters.push_back(N_UTL_Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("CURRZEROTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("HISTORYTRACKINGDEPTH", 50));
  optionsMetadata_[string("TIMEINT")] = optionsParameters;

  // Make a copy for MPDE time integration.  This copy will result in MPDE having
  // the same defaults, so if we want to change defaults later, we'll have to
  // set this up explicitly.
  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("METHOD", 1));
#ifdef Xyce_DEBUG_ANALYSIS
  optionsParameters.push_back(N_UTL_Param("CONSTSTEP", 0));
#endif
  optionsParameters.push_back(N_UTL_Param("USEDEVICEMAX", 1));
  optionsParameters.push_back(N_UTL_Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(N_UTL_Param("RESTARTSTEPSCALE", .005));
//  optionsParameters.push_back(N_UTL_Param("NLNEARCONV", 1));
  optionsParameters.push_back(N_UTL_Param("NLNEARCONV", 0));
  optionsParameters.push_back(N_UTL_Param("NLSMALLUPDATE", 1));
  optionsParameters.push_back(N_UTL_Param("DOUBLEDCOPSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("FIRSTDCOPSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("LASTDCOPSTEP", 1));
  optionsParameters.push_back(N_UTL_Param("RESETTRANNLS", 1));
  optionsParameters.push_back(N_UTL_Param("BPENABLE", 1));
  optionsParameters.push_back(N_UTL_Param("EXITTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("EXITSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("ERROPTION", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 0));
  optionsParameters.push_back(N_UTL_Param("JACLIMITFLAG", 0));
  optionsParameters.push_back(N_UTL_Param("JACLIMIT", 1.0e17));
  optionsParameters.push_back(N_UTL_Param("DAESTATEDERIV", 0));
  optionsParameters.push_back(N_UTL_Param("TESTFIRSTSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DTMIN", 0.0));
  optionsParameters.push_back(N_UTL_Param("MAXORD", 5));
  optionsParameters.push_back(N_UTL_Param("MINORD", 1));
  optionsParameters.push_back(N_UTL_Param("OUTPUTINTERPMPDE", 1));
  optionsParameters.push_back(N_UTL_Param("INTERPOUTPUT", 1));
  optionsParameters.push_back(N_UTL_Param("CONDTEST", 0));
  optionsParameters.push_back(N_UTL_Param("CONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(N_UTL_Param("ISOCONDTEST", 0));
  optionsParameters.push_back(N_UTL_Param("ISOCONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(N_UTL_Param("MINTIMESTEPSBP", 10));
  optionsParameters.push_back(N_UTL_Param("NLMIN", 3));
  optionsParameters.push_back(N_UTL_Param("NLMAX", 8));
  optionsParameters.push_back(N_UTL_Param("DELMAX", 1.0e+99));
  optionsParameters.push_back(N_UTL_Param("TIMESTEPSREVERSAL", false));
  optionsMetadata_[string("TIMEINT-MPDE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 1));
  optionsMetadata_[string("CONDUCTANCE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("NLSTRATEGY", 0));
  optionsParameters.push_back(N_UTL_Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(N_UTL_Param("NOX", 1));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("RELTOL", 1.0E-3));
  optionsParameters.push_back(N_UTL_Param("DELTAXTOL", 1.0));
  optionsParameters.push_back(N_UTL_Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("RHSTOL", 1.0E-6));
  optionsParameters.push_back(N_UTL_Param("MAXSTEP", 200));
  optionsParameters.push_back(N_UTL_Param("MAXSEARCHSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("NORMLVL", 2));
  optionsParameters.push_back(N_UTL_Param("LINOPT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMAX", N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMIN", -N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(N_UTL_Param("IN_FORCING", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("DLSDEBUG", 0));
  optionsParameters.push_back(N_UTL_Param("MATRIXMARKET", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(N_UTL_Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(N_UTL_Param("MEMORY", 400));
  optionsParameters.push_back(N_UTL_Param("CONTINUATION", 0));
  optionsParameters.push_back(N_UTL_Param("ENFORCEDEVICECONV", 1));
  optionsParameters.push_back(N_UTL_Param("FASTTESTS", false));
  optionsParameters.push_back(N_UTL_Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[string("NONLIN")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("NLSTRATEGY", 0));
  optionsParameters.push_back(N_UTL_Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(N_UTL_Param("NOX", 1));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(N_UTL_Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(N_UTL_Param("DELTAXTOL", 0.33));
  optionsParameters.push_back(N_UTL_Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("RHSTOL", 1.0E-2));
  optionsParameters.push_back(N_UTL_Param("MAXSTEP", 20));
  optionsParameters.push_back(N_UTL_Param("MAXSEARCHSTEP", 2));
  optionsParameters.push_back(N_UTL_Param("NORMLVL", 2));
  optionsParameters.push_back(N_UTL_Param("LINOPT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMAX", N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMIN", -N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(N_UTL_Param("IN_FORCING", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("DLSDEBUG", 0));
  optionsParameters.push_back(N_UTL_Param("MATRIXMARKET", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(N_UTL_Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(N_UTL_Param("MEMORY", 400));
  optionsParameters.push_back(N_UTL_Param("CONTINUATION", 0));
  optionsParameters.push_back(N_UTL_Param("ENFORCEDEVICECONV", 0));
  optionsParameters.push_back(N_UTL_Param("FASTTESTS", false));
  optionsParameters.push_back(N_UTL_Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[string("NONLIN-TRAN")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("NLSTRATEGY", 0));
  optionsParameters.push_back(N_UTL_Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(N_UTL_Param("NOX", 0));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0E-9));
  optionsParameters.push_back(N_UTL_Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(N_UTL_Param("DELTAXTOL", 0.33));
  optionsParameters.push_back(N_UTL_Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("RHSTOL", 1.0E-6));
  optionsParameters.push_back(N_UTL_Param("MAXSTEP", 200));
  optionsParameters.push_back(N_UTL_Param("MAXSEARCHSTEP", 2));
  optionsParameters.push_back(N_UTL_Param("NORMLVL", 2));
  optionsParameters.push_back(N_UTL_Param("LINOPT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMAX", N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMIN", -N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(N_UTL_Param("IN_FORCING", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("DLSDEBUG", 0));
  optionsParameters.push_back(N_UTL_Param("MATRIXMARKET", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(N_UTL_Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(N_UTL_Param("MEMORY", 400));
  optionsParameters.push_back(N_UTL_Param("CONTINUATION", 0));
  optionsParameters.push_back(N_UTL_Param("ENFORCEDEVICECONV", 0));
  optionsMetadata_[string("NONLIN-HB")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("NLSTRATEGY", 0));
  optionsParameters.push_back(N_UTL_Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(N_UTL_Param("NOX", 1));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("RELTOL", 1.0E-3));
  optionsParameters.push_back(N_UTL_Param("DELTAXTOL", 1.0));
  optionsParameters.push_back(N_UTL_Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("RHSTOL", 1.0E-6));
  optionsParameters.push_back(N_UTL_Param("MAXSTEP", 200));
  optionsParameters.push_back(N_UTL_Param("MAXSEARCHSTEP", 0));
  optionsParameters.push_back(N_UTL_Param("NORMLVL", 2));
  optionsParameters.push_back(N_UTL_Param("LINOPT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMAX", N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMIN", -N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(N_UTL_Param("IN_FORCING", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("DLSDEBUG", 0));
  optionsParameters.push_back(N_UTL_Param("MATRIXMARKET", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(N_UTL_Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(N_UTL_Param("MEMORY", 400));
  optionsParameters.push_back(N_UTL_Param("CONTINUATION", 0));
  optionsParameters.push_back(N_UTL_Param("ENFORCEDEVICECONV", 0));
  optionsParameters.push_back(N_UTL_Param("ALGORITHM", 0));
  optionsParameters.push_back(N_UTL_Param("MAXCONTSTEPS", 0));
  optionsParameters.push_back(N_UTL_Param("CONTINUATIONFLAG", 1));
  optionsParameters.push_back(N_UTL_Param("INNERFAIL", 1));
  optionsParameters.push_back(N_UTL_Param("EXITWITHFAILURE", 1));
  optionsParameters.push_back(N_UTL_Param("FULLNEWTONENFORCE", 1));
  optionsParameters.push_back(N_UTL_Param("CONPARAM", "VA:V0"));
  optionsParameters.push_back(N_UTL_Param("VOLTLIMTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("FASTTESTS", false));
  optionsParameters.push_back(N_UTL_Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[string("NONLIN-TWOLEVEL")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("NLSTRATEGY", 0));
  optionsParameters.push_back(N_UTL_Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(N_UTL_Param("NOX", 1));
  optionsParameters.push_back(N_UTL_Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(N_UTL_Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(N_UTL_Param("DELTAXTOL", 0.33));
  optionsParameters.push_back(N_UTL_Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("RHSTOL", 1.0E-2));
  optionsParameters.push_back(N_UTL_Param("MAXSTEP", 20));
  optionsParameters.push_back(N_UTL_Param("MAXSEARCHSTEP", 2));
  optionsParameters.push_back(N_UTL_Param("NORMLVL", 2));
  optionsParameters.push_back(N_UTL_Param("LINOPT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMAX", N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTMIN", -N_UTL_MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(N_UTL_Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(N_UTL_Param("IN_FORCING", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("DLSDEBUG", 0));
  optionsParameters.push_back(N_UTL_Param("MATRIXMARKET", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(N_UTL_Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(N_UTL_Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(N_UTL_Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(N_UTL_Param("MEMORY", 400));
  optionsParameters.push_back(N_UTL_Param("CONTINUATION", 0));
  optionsParameters.push_back(N_UTL_Param("ENFORCEDEVICECONV", 0));
  optionsParameters.push_back(N_UTL_Param("ALGORITHM", 0));
  optionsParameters.push_back(N_UTL_Param("MAXCONTSTEPS", 0));
  optionsParameters.push_back(N_UTL_Param("CONTINUATIONFLAG", 1));
  optionsParameters.push_back(N_UTL_Param("INNERFAIL", 1));
  optionsParameters.push_back(N_UTL_Param("EXITWITHFAILURE", 1));
  optionsParameters.push_back(N_UTL_Param("FULLNEWTONENFORCE", 1));
  optionsParameters.push_back(N_UTL_Param("CONPARAM", "VA:V0"));
  optionsParameters.push_back(N_UTL_Param("VOLTLIMTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("FASTTESTS", false));
  optionsParameters.push_back(N_UTL_Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(N_UTL_Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[string("NONLIN-TWOLEVEL-TRAN")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("STEPPER", "NATURAL"));
  optionsParameters.push_back(N_UTL_Param("PREDICTOR", "CONSTANT"));
  optionsParameters.push_back(N_UTL_Param("STEPCONTROL", "CONSTANT"));

  optionsParameters.push_back(N_UTL_Param("CONPARAM", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("INITIALVALUE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MAXVALUE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MINVALUE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("INITIALSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MAXSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MINSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("AGGRESSIVENESS", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("BIFPARAM", "VA:V0"));
  optionsParameters.push_back(N_UTL_Param("MAXSTEPS", 20));
  optionsParameters.push_back(N_UTL_Param("MAXNLITERS", 20));
  optionsParameters.push_back(N_UTL_Param("PARAMLIST", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("VOLTAGELIST", "DOFS"));
  optionsParameters.push_back(N_UTL_Param("VOLTAGESCALEFACTOR", 1.0));
  optionsParameters.push_back(N_UTL_Param("RESIDUALCONDUCTANCE", 0.0));
  optionsMetadata_[string("LOCA")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("STEPPER", "NATURAL"));
  optionsParameters.push_back(N_UTL_Param("PREDICTOR", "CONSTANT"));
  optionsParameters.push_back(N_UTL_Param("STEPCONTROL", "CONSTANT"));
  optionsParameters.push_back(N_UTL_Param("CONPARAM", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("INITIALVALUE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MAXVALUE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MINVALUE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("INITIALSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MAXSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("MINSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("AGGRESSIVENESS", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("BIFPARAM", "VA:V0"));
  optionsParameters.push_back(N_UTL_Param("MAXSTEPS", 20));
  optionsParameters.push_back(N_UTL_Param("MAXNLITERS", 20));
  optionsParameters.push_back(N_UTL_Param("PARAMLIST", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("VOLTAGELIST", "DOFS"));
  optionsParameters.push_back(N_UTL_Param("VOLTAGESCALEFACTOR", 1.0));
  optionsParameters.push_back(N_UTL_Param("RESIDUALCONDUCTANCE", 0.0));
  optionsMetadata_[string("TWOLEVEL-LOCA")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("AZ_max_iter", 500));
  optionsParameters.push_back(N_UTL_Param("AZ_precond", 14));
  optionsParameters.push_back(N_UTL_Param("AZ_solver", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_conv", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_pre_calc", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_keep_info", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_orthog", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_subdomain_solve", 9));
  optionsParameters.push_back(N_UTL_Param("AZ_ilut_fill", 3.0));
  optionsParameters.push_back(N_UTL_Param("AZ_drop", 1.0E-3));
  optionsParameters.push_back(N_UTL_Param("AZ_reorder", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_scaling", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_kspace", 500));
  optionsParameters.push_back(N_UTL_Param("AZ_tol", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("AZ_output", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_diagnostics", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_overlap", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_rthresh", 1.0001));
  optionsParameters.push_back(N_UTL_Param("AZ_athresh", 1.0E-4));
  optionsParameters.push_back(N_UTL_Param("AZ_filter", 0.0));
  optionsParameters.push_back(N_UTL_Param("TR_filter", 0.0));
  optionsParameters.push_back(N_UTL_Param("TR_partition", 1));
#ifdef Xyce_SHYLU
  optionsParameters.push_back(N_UTL_Param("ShyLU_rthresh", 1.0E-3));
#endif
#ifdef Xyce_USE_ISORROPIA
  optionsParameters.push_back(N_UTL_Param("TR_partition_type", "GRAPH"));
#endif
  optionsParameters.push_back(N_UTL_Param("TR_reindex", 1));
  optionsParameters.push_back(N_UTL_Param("TR_solvermap", 1));
  optionsParameters.push_back(N_UTL_Param("TR_amd", 1));
  optionsParameters.push_back(N_UTL_Param("TR_btf", 0));
  optionsParameters.push_back(N_UTL_Param("TR_global_btf", 0));
  optionsParameters.push_back(N_UTL_Param("TR_global_btf_droptol", 1.0E-16));
  optionsParameters.push_back(N_UTL_Param("TR_global_btf_verbose", 0));
#ifdef Xyce_TRILINOS_DEV
  optionsParameters.push_back(N_UTL_Param("TR_global_amd", 0));
  optionsParameters.push_back(N_UTL_Param("TR_global_amd_verbose", 0));
#endif
  optionsParameters.push_back(N_UTL_Param("TR_singleton_filter", 0));
  optionsParameters.push_back(N_UTL_Param("SLU_EQUILIBRATE", 1));
  optionsParameters.push_back(N_UTL_Param("SLU_REFACTOR", 1));
  optionsParameters.push_back(N_UTL_Param("SLU_PERMUTE", 2));
  optionsParameters.push_back(N_UTL_Param("SLU_PIVOT_THRESH", -1.0));
  optionsParameters.push_back(N_UTL_Param("SLU_FILL_FAC", -1));
  optionsParameters.push_back(N_UTL_Param("BTF", 0));
  optionsParameters.push_back(N_UTL_Param("BTF_VERBOSE", 0));
  optionsParameters.push_back(N_UTL_Param("BTF_ATHRESH", 0.0));
  optionsParameters.push_back(N_UTL_Param("BTF_RTHRESH", 0.0));
  optionsParameters.push_back(N_UTL_Param("BTF_RTHRESH_INIT", 0.0));
  optionsParameters.push_back(N_UTL_Param("BTF_INIT", 0));
  optionsParameters.push_back(N_UTL_Param("BTF_THRESHOLDING", 0));
  optionsParameters.push_back(N_UTL_Param("BTF_RNTHRESHFAC", 1.0e-3));
  optionsParameters.push_back(N_UTL_Param("adaptive_solve", 0));
  optionsParameters.push_back(N_UTL_Param("use_aztec_precond", 1));
  optionsParameters.push_back(N_UTL_Param("use_ifpack_factory", 0));
  optionsParameters.push_back(N_UTL_Param("ifpack_type", "Amesos"));
  optionsParameters.push_back(N_UTL_Param("diag_perturb", 0.0));
#ifdef Xyce_ML
  optionsParameters.push_back(N_UTL_Param("ML_max_level", 5));
#endif
  optionsParameters.push_back(N_UTL_Param("TR_rcm", 0));
  optionsParameters.push_back(N_UTL_Param("TR_scale", 0));
  optionsParameters.push_back(N_UTL_Param("TR_scale_left", 0));
  optionsParameters.push_back(N_UTL_Param("TR_scale_right", 0));
  optionsParameters.push_back(N_UTL_Param("TR_scale_exp", 1.0));
  optionsParameters.push_back(N_UTL_Param("TR_scale_iter", 0));
  optionsParameters.push_back(N_UTL_Param("TYPE", "DEFAULT"));
  optionsParameters.push_back(N_UTL_Param("PREC_TYPE", "DEFAULT"));
#ifdef Xyce_BELOS
  optionsParameters.push_back(N_UTL_Param("BELOS_SOLVER_TYPE", "Block GMRES"));
#endif
  optionsParameters.push_back(N_UTL_Param("KLU_REPIVOT", 1));
  optionsParameters.push_back(N_UTL_Param("OUTPUT_LS", 1));
  optionsParameters.push_back(N_UTL_Param("OUTPUT_BASE_LS", 1));
  optionsParameters.push_back(N_UTL_Param("OUTPUT_FAILED_LS", 1));
  optionsMetadata_[string("LINSOL")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("AZ_max_iter", 500));
  optionsParameters.push_back(N_UTL_Param("AZ_solver", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_conv", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_pre_calc", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_keep_info", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_orthog", 1));
  optionsParameters.push_back(N_UTL_Param("AZ_reorder", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_scaling", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_kspace", 500));
  optionsParameters.push_back(N_UTL_Param("AZ_tol", 1.0E-12));
  optionsParameters.push_back(N_UTL_Param("AZ_output", 0));
  optionsParameters.push_back(N_UTL_Param("AZ_diagnostics", 0));
  optionsParameters.push_back(N_UTL_Param("TYPE", "DEFAULT"));
  optionsParameters.push_back(N_UTL_Param("PREC_TYPE", "NONE"));
  optionsMetadata_[string("LINSOL-HB")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("PARTITIONER", "CHACO"));
  optionsParameters.push_back(N_UTL_Param("DISTRIBINDSRCNODES", 0.0));
  optionsParameters.push_back(N_UTL_Param("USE_WEIGHTS", 0));
  optionsParameters.push_back(N_UTL_Param("METHOD", "KWAY"));
  optionsParameters.push_back(N_UTL_Param("USE_VWEIGHTS", 0));
  optionsParameters.push_back(N_UTL_Param("USE_EWEIGHTS", 0));
  optionsParameters.push_back(N_UTL_Param("BALANCE_FACTOR", 0));
  optionsParameters.push_back(N_UTL_Param("BISECTIONS", 10));
  optionsParameters.push_back(N_UTL_Param("VERTEX_GROUPING", 1));
  optionsParameters.push_back(N_UTL_Param("REFINEMENT", 1));
  optionsParameters.push_back(N_UTL_Param("VCYCLE", 1));
  optionsParameters.push_back(N_UTL_Param("RECONST", 0));
  optionsParameters.push_back(N_UTL_Param("PREASSIGN", 0));
  optionsParameters.push_back(N_UTL_Param("SEED", -1));
  optionsParameters.push_back(N_UTL_Param("DEBUG", 0));
  optionsMetadata_[string("PARALLEL")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("PACK", 1));
  optionsParameters.push_back(N_UTL_Param("JOB", ""));
  optionsParameters.push_back(N_UTL_Param("START_TIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("FILE", ""));
  optionsParameters.push_back(N_UTL_Param("INITIAL_INTERVAL", 0.0));
  optionsParameters.push_back(N_UTL_Param("TIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("INTERVAL", 0.0));
  optionsMetadata_[string("RESTART")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("INITIAL_INTERVAL", 0.0));
  optionsParameters.push_back(N_UTL_Param("TIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("INTERVAL", 0.0));
  optionsParameters.push_back(N_UTL_Param("HDF5FILENAME", ""));
  optionsParameters.push_back(N_UTL_Param("PRINTENDOFSIMLINE", true));
  optionsMetadata_[string("OUTPUT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("INPUT", ""));
  optionsParameters.push_back(N_UTL_Param("OUTPUT", ""));
  optionsMetadata_[string("OP_IO")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("TIME", 0.0));
  optionsParameters.push_back(N_UTL_Param("INTERVAL", 0.0));
  optionsMetadata_[string("OUTPUT-LINE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("TYPE", "TRAN"));
  optionsParameters.push_back(N_UTL_Param("FILE", ""));
  optionsParameters.push_back(N_UTL_Param("FORMAT", "STD"));
  optionsParameters.push_back(N_UTL_Param("DELIMITER", ""));
  optionsParameters.push_back(N_UTL_Param("WIDTH", 17));
  optionsParameters.push_back(N_UTL_Param("PRECISION", 8));
  optionsParameters.push_back(N_UTL_Param("TIMESCALEFACTOR", 1.0));
  optionsParameters.push_back(N_UTL_Param("FILTER", 0.0));
  optionsMetadata_[string("PRINT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("REPLICATED_CKT", 1));
  optionsParameters.push_back(N_UTL_Param("CHECK_CONNECTIVITY", 0));
  optionsParameters.push_back(N_UTL_Param("SUPERNODE", false));
  optionsMetadata_[string("TOPOLOGY")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("TEST", 0));
  optionsParameters.push_back(N_UTL_Param("OSCSRC", "VECTOR"));
  optionsParameters.push_back(N_UTL_Param("OSCOUT", ""));
  optionsParameters.push_back(N_UTL_Param("AUTON2", false));
  optionsParameters.push_back(N_UTL_Param("AUTON2MAX", 100));
  optionsParameters.push_back(N_UTL_Param("STARTUPPERIODS", 0));
  optionsParameters.push_back(N_UTL_Param("SAVEICDATA", false));
  optionsParameters.push_back(N_UTL_Param("N2", 10));
  optionsParameters.push_back(N_UTL_Param("T2", 0.0));
  optionsParameters.push_back(N_UTL_Param("IC", 0));
  optionsParameters.push_back(N_UTL_Param("ICPER", 10));
  optionsParameters.push_back(N_UTL_Param("DIFF", 0));
  optionsParameters.push_back(N_UTL_Param("DIFFORDER", 1));
  optionsParameters.push_back(N_UTL_Param("FREQDOMAIN", 0));
  optionsParameters.push_back(N_UTL_Param("WAMPDE", 0));
  optionsParameters.push_back(N_UTL_Param("DCOPEXIT", 0));
  optionsParameters.push_back(N_UTL_Param("ICEXIT", 0));
  optionsParameters.push_back(N_UTL_Param("EXITSAWTOOTHSTEP", -1));
  optionsParameters.push_back(N_UTL_Param("PHASE", 0));
  optionsParameters.push_back(N_UTL_Param("PHASECOEFF", 0));
  optionsParameters.push_back(N_UTL_Param("NONLTESTEPS", 10));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 0));
  optionsMetadata_[string("MPDEINT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("TEST", false));
  optionsParameters.push_back(N_UTL_Param("NUMFREQ", 21));
  optionsParameters.push_back(N_UTL_Param("STARTUPPERIODS", 0));
  optionsParameters.push_back(N_UTL_Param("SAVEICDATA", false));
  optionsParameters.push_back(N_UTL_Param("DEBUGLEVEL", 0));
  optionsMetadata_[string("HBINT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(N_UTL_Param("METHOD", "PRIMA"));
  optionsParameters.push_back(N_UTL_Param("SAVEREDSYS", false));
  optionsParameters.push_back(N_UTL_Param("COMPORIGTF", false));
  optionsParameters.push_back(N_UTL_Param("COMPREDTF", false));
  optionsParameters.push_back(N_UTL_Param("COMPTYPE", "DEC"));
  optionsParameters.push_back(N_UTL_Param("COMPNP", 10));
  optionsParameters.push_back(N_UTL_Param("COMPFSTART", 1.0));
  optionsParameters.push_back(N_UTL_Param("COMPFSTOP", 1.0));
  optionsParameters.push_back(N_UTL_Param("EXPPOINT", 0.0));
  optionsParameters.push_back(N_UTL_Param("SCALETYPE", 0));
  optionsParameters.push_back(N_UTL_Param("SCALEFACTOR", 1));
  optionsParameters.push_back(N_UTL_Param("SCALEFACTOR1", 0.01));
  optionsParameters.push_back(N_UTL_Param("SPARSIFICATIONTYPE", 0));
  optionsMetadata_[string("MOR_OPTS")] = optionsParameters;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::sourceFunctionMetadata
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
void CircuitMetadata::sourceFunctionMetadata()
{
  vector<Device::Param> sourceFcnParamList;
  // Source function metadata.
  sourceFcnParamList.resize(7);
  sourceFcnParamList[0].set("V1", 0.0);
  sourceFcnParamList[1].set("V2", 0.0);
  sourceFcnParamList[2].set("TD", 0.0);
  sourceFcnParamList[3].set("TR", 0.0);
  sourceFcnParamList[4].set("TF", 0.0);
  sourceFcnParamList[5].set("PW", 0.0);
  sourceFcnParamList[6].set("PER", 0.0);
  sourceFcnMap[ "PULSE" ] = sourceFcnParamList;
  sourceFcnParamList.resize(6);
  sourceFcnParamList[0].set("V0", 0.0);
  sourceFcnParamList[1].set("VA", 0.0);
  sourceFcnParamList[2].set("FREQ", 0.0);
  sourceFcnParamList[3].set("TD", 0.0);
  sourceFcnParamList[4].set("THETA", 0.0);
  sourceFcnParamList[5].set("PHASE", 0.0);
  sourceFcnMap[ "SIN" ] = sourceFcnParamList;
  sourceFcnParamList[0].set("V1", 0.0);
  sourceFcnParamList[1].set("V2", 0.0);
  sourceFcnParamList[2].set("TD1", 0.0);
  sourceFcnParamList[3].set("TAU1", 0.0);
  sourceFcnParamList[4].set("TD2", 0.0);
  sourceFcnParamList[5].set("TAU2", 0.0);
  sourceFcnMap[ "EXP" ] = sourceFcnParamList;
  sourceFcnParamList.resize(5);
  sourceFcnParamList[0].set("V0", 0.0);
  sourceFcnParamList[1].set("VA", 0.0);
  sourceFcnParamList[2].set("FC", 0.0);
  sourceFcnParamList[3].set("MDI", 0.0);
  sourceFcnParamList[4].set("FS", 0.0);
  sourceFcnMap[ "SFFM" ] = sourceFcnParamList;
  // sourceFcnParamList is empty for source function "PWL" indicating
  // that it needs special handling in NetlistHandler.C
  sourceFcnParamList.clear();
  sourceFcnMap[ "PWL" ] = sourceFcnParamList;
  sourceFcnParamList.resize(8);
  sourceFcnParamList[0].set("V1", 0.0);
  sourceFcnParamList[1].set("V2", 0.0);
  sourceFcnParamList[2].set("TD", 0.0);
  sourceFcnParamList[3].set("TR", 0.0);
  sourceFcnParamList[4].set("TF", 0.0);
  sourceFcnParamList[5].set("PW", 0.0);
  sourceFcnParamList[6].set("PER", 0.0);
  sourceFcnParamList[7].set("SF", 0.0);
  sourceFcnMap[ "SMOOTHPULSE" ] = sourceFcnParamList;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::isModelTypeValid
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
bool CircuitMetadata::isModelTypeValid(
  const string &        deviceType,
  const string &        modelType,
  int                   modelLevel)
{
  DeviceMetadata* DMptr = findDeviceMetadata(deviceType, modelLevel);
  if (DMptr == NULL)
    return false;

  vector<string>::iterator stringIter = find(DMptr->modelTypes.begin(), DMptr->modelTypes.end(), modelType);

  return stringIter != DMptr->modelTypes.end();
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::isDeviceParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
bool CircuitMetadata::isDeviceParameter(
    const string & deviceType,
    int modelLevel,
    const string & parameterName)
{
  DeviceMetadata* DMptr = findDeviceMetadata(deviceType, modelLevel);
  if (DMptr == NULL)
    return false;

  vector<Device::Param>::iterator paramIter =
    find(DMptr->instanceParameters.begin(),
         DMptr->instanceParameters.end(),
         Device::Param(parameterName, ""));

  return paramIter != DMptr->instanceParameters.end();
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getPrimaryParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
string CircuitMetadata::getPrimaryParameter(
  const string &        deviceType,
  int                   modelLevel)
{
  DeviceMetadata* DMptr = findDeviceMetadata(deviceType, modelLevel);
  if (DMptr == NULL)
    return string("");

  return DMptr->primaryParameter;
}


//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getNumberOfNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
int CircuitMetadata::getNumberOfNodes(
  const string &        deviceType,
  int                   modelLevel)
{
  DeviceMetadata* DMptr = findDeviceMetadata(deviceType, modelLevel);
  if (DMptr == NULL)
    return -1;

  return DMptr->numNodes;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getNumberOfOptionalNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
int CircuitMetadata::getNumberOfOptionalNodes(
  const string &        deviceType,
  int                   modelLevel)
{
  DeviceMetadata* DMptr = findDeviceMetadata(deviceType, modelLevel);
  if (DMptr == NULL)
    return -1;

  return DMptr->numOptionalNodes;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getNumberOfFillNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/04/2004
//-----------------------------------------------------------------------------
int CircuitMetadata::getNumberOfFillNodes(
  const string &        deviceType,
  int                   modelLevel)
{
  DeviceMetadata* DMptr = findDeviceMetadata(deviceType, modelLevel);
  if (DMptr == NULL)
    return -1;

  return DMptr->numFillNodes;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::isModelRequired
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
bool CircuitMetadata::isModelRequired(
  const string &        deviceType,
  int                   modelLevel)
{
  DeviceMetadata* DMptr = findDeviceMetadata(deviceType, modelLevel);
  if (DMptr == NULL)
    return false;

  return DMptr->modelRequired;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getSourceFunctionID
// Purpose       : Return the integer value corresponding to the
//                 given source function. Return NUM_SRC_DATA if the given string
//                 is unknown.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
int CircuitMetadata::getSourceFunctionID(const string & sourceFcn)
{
  if (sourceFcn == "PULSE")     return _PULSE_DATA;
  else if (sourceFcn == "SIN")  return _SIN_DATA;
  else if (sourceFcn == "EXP")  return _EXP_DATA;
  else if (sourceFcn == "SFFM") return _SFFM_DATA;
  else if (sourceFcn == "PWL")  return _PWL_DATA;
  else if (sourceFcn == "DC")   return _DC_DATA;
  else if (sourceFcn == "SMOOTHPULSE") return _SMOOTH_PULSE_DATA;
  else if (sourceFcn == "AC")   return _AC_DATA;  //tmei: 05/02
  //else if (sourceFcn == "DISTOF1")   return _DISTOF1_DATA;
  //else if (sourceFcn == "DISTOF2")   return _DISTOF2_DATA;
  else return _NUM_SRC_DATA;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getPtrToInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
vector<Device::Param> * CircuitMetadata::getPtrToInstanceParameters(
  const string &        deviceType,
  int                   modelLevel)
{
  if (modelLevel == -1 || deviceMetadataIndex.find(DeviceLevelKey(deviceType, modelLevel)) == deviceMetadataIndex.end())
    return &(findDeviceMetadata(deviceType, modelLevel)->instanceParameters);
  else
    return &deviceMetadata_[deviceMetadataIndex[DeviceLevelKey(deviceType, modelLevel)]].instanceParameters;
}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::getInstanceCompositeComponents
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 10/07/2003
//----------------------------------------------------------------------------
void CircuitMetadata::getInstanceCompositeComponents(
  const string &                deviceType,
  const string &                parameterName, int modelLevel,
  vector<Device::Param> &       components)
{
  DeviceMetadata *DMptr = findDeviceMetadata(deviceType, modelLevel);
  DeviceParamMap &icpMap = DMptr->instanceCompositeParameterMap;
  DeviceParamMap::iterator iterIcp = icpMap.find(parameterName);

  if ( iterIcp != icpMap.end() )
  {
    components = iterIcp->second;
  }
  else
  {
    std::ostringstream oss;
    oss << "There are no component parameters in metadata for the VECTOR-COMPOSITE parameter: " << parameterName;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str() );
  }
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getPtrToModelParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
vector<Device::Param> * CircuitMetadata::getPtrToModelParameters(
  const string &        modelType,
  int                   modelLevel)
{
  if (modelLevel == -1 || deviceMetadataIndex.find(DeviceLevelKey(modelType, modelLevel)) == deviceMetadataIndex.end())
    return &(findDeviceMetadata(modelType, modelLevel)->modelParameters);
  else
    return &deviceMetadata_[deviceMetadataIndex[DeviceLevelKey(modelType, modelLevel)]].modelParameters;
}


//----------------------------------------------------------------------------
// Function       : CircuitMetadata::getModelCompositeComponents
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Eric Keiter
// Creation Date  : 5/07/2008
//----------------------------------------------------------------------------
void CircuitMetadata::getModelCompositeComponents(
  const string &                modelType,
  const string &                parameterName, int modelLevel,
  vector<Device::Param> &       components)
{
  DeviceMetadata *DMptr = findDeviceMetadata(modelType, modelLevel);
  DeviceParamMap &mcpMap = DMptr->modelCompositeParameterMap;

  DeviceParamMap::iterator iterIcp = mcpMap.find(parameterName);

  if ( iterIcp != mcpMap.end() )
  {
    components = iterIcp->second;
  }
  else
  {
    string msg("There are no component parameters in metadata for");
    msg += "the VECTOR-COMPOSITE parameter: ";
    msg += parameterName + ".\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::isOptionsPackage
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 10/07/2003
//----------------------------------------------------------------------------
bool CircuitMetadata::isOptionsPackage(const string & optionName)
{
  if (optionsMetadata_.count(optionName))
    return true;
  else
    return false;
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getPtrToOptionsParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
vector<N_UTL_Param> * CircuitMetadata::getPtrToOptionsParameters(const string & optionName)
{
  return &(optionsMetadata_[optionName]);
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getSourceFunctionParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
void CircuitMetadata::getSourceFunctionParameters(const string & sourceFcn,
    vector<Device::Param> & sourceFunctionParameters)
{
  sourceFunctionParameters = sourceFcnMap[sourceFcn];
}

} // namespace IO
} // namespace Xyce
