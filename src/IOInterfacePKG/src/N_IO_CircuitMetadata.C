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
// Revision Number: $Revision: 1.189.2.1 $
//
// Revision Date  : $Date: 2014/03/06 00:41:20 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <N_IO_CircuitMetadata.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_SourceData.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace IO {

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::getDeviceMetadata
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 09/28/2003
//----------------------------------------------------------------------------
DeviceMetadata &
CircuitMetadata::getDeviceMetadata(
  const std::string &   deviceTypeIn,
  int                   level) const
{
  if (level == -1)
  {
    level = 1;
  }

  std::string deviceType = deviceTypeIn;
  if (deviceTypeIn == "K")
  {
    deviceType = "L";
  }

  DeviceMetadataIndexMap::const_iterator it_device_meta_index = deviceMetadataIndex.find(NameLevelKey(deviceType, level));
  if (it_device_meta_index != deviceMetadataIndex.end())
  {
    return deviceMetadata_[it_device_meta_index->second];
  }

  // Handle default model:
  DeviceMetadata &device_metadata = deviceMetadata_[NameLevelKey(deviceType, level)];

  const Device::Configuration *configuration = Device::Configuration::findConfiguration(deviceType, level);
  if (configuration) {
    device_metadata.levelValid = true;
    device_metadata.numOptionalNodes = configuration->getNumOptionalNodes();
    device_metadata.numFillNodes = configuration->getNumFillNodes();
    device_metadata.numNodes = configuration->getNumNodes();
    device_metadata.modelRequired = configuration->getModelRequired();
    device_metadata.primaryParameter = configuration->getPrimaryParameter();
    device_metadata.modelTypes = configuration->getModelTypeNames();

    Xyce::Device::populateParams(configuration->getModelParameters().getMap(), device_metadata.modelParameters, device_metadata.modelCompositeParameterMap);
    Xyce::Device::populateParams(configuration->getInstanceParameters().getMap(), device_metadata.instanceParameters, device_metadata.instanceCompositeParameterMap);
  }

  for (std::vector<std::string>::iterator it_device_model_type = device_metadata.modelTypes.begin(); it_device_model_type != device_metadata.modelTypes.end(); ++it_device_model_type)
  {
    deviceMetadataIndex.insert(DeviceMetadataIndexMap::value_type(NameLevelKey(*it_device_model_type, level), NameLevelKey(deviceType, level)));
  }
  deviceMetadataIndex[NameLevelKey(deviceType, level)] = NameLevelKey(deviceType, level);
  if (deviceType == "L")
    deviceMetadataIndex[NameLevelKey(std::string("K"), level)] = NameLevelKey(deviceType, level);

  return device_metadata;
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
  std::string input_line, parameter, value, optionsPackage;
  std::vector<Util::Param> optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("DEFAS", 0.0));
  optionsParameters.push_back(Util::Param("DEFAD", 0.0));
  optionsParameters.push_back(Util::Param("DEFL", 1.0e-4));
  optionsParameters.push_back(Util::Param("DEFW", 1.0e-4));
  optionsParameters.push_back(Util::Param("GMIN", 1.0e-12));
  optionsParameters.push_back(Util::Param("GMINSCALAR", 1.0e+10));
  optionsParameters.push_back(Util::Param("GMAX", 1.0e20));
  optionsParameters.push_back(Util::Param("TEMP", 27.0));
  optionsParameters.push_back(Util::Param("TNOM", 27.0));
  optionsParameters.push_back(Util::Param("SCALESRC", 0.0));
  optionsParameters.push_back(Util::Param("NUMJAC", 0));
  optionsParameters.push_back(Util::Param("TESTJAC", 0));
  optionsParameters.push_back(Util::Param("TESTJACSTARTSTEP", 0));
  optionsParameters.push_back(Util::Param("TESTJACSTOPSTEP", 0));
  optionsParameters.push_back(Util::Param("TJRELTOL", 0.01));
  optionsParameters.push_back(Util::Param("TJABSTOL", 1.0e-8));
  optionsParameters.push_back(Util::Param("TJSQRTETA", 1.0e-8));
  optionsParameters.push_back(Util::Param("SENSDP", 1.0e-8));
  optionsParameters.push_back(Util::Param("TESTJACWARN", 0));
  optionsParameters.push_back(Util::Param("TESTJACDEVICENAME", ""));
  optionsParameters.push_back(Util::Param("VOLTLIM", 1));
  optionsParameters.push_back(Util::Param("ICFAC", 10000.0));
  optionsParameters.push_back(Util::Param("LAMBERTW", 0));
  optionsParameters.push_back(Util::Param("MAXTIMESTEP", 1.0e99));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(Util::Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIMESTEP", 65536));
  optionsParameters.push_back(Util::Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIME", 100.0));
  optionsParameters.push_back(Util::Param("VERBOSELEVEL", 0));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0e-12));
  optionsParameters.push_back(Util::Param("CHGTOL", 1.0e-12));
  optionsParameters.push_back(Util::Param("VDSSCALEMIN", 0.3));
  optionsParameters.push_back(Util::Param("VGSTCONST", 4.5));
  optionsParameters.push_back(Util::Param("LENGTH0", 5.0e-6));
  optionsParameters.push_back(Util::Param("WIDTH0", 200.0e-6));
  optionsParameters.push_back(Util::Param("TOX0", 6.0e-8));
  optionsParameters.push_back(Util::Param("MINRES", 0.0));
  optionsParameters.push_back(Util::Param("MINCAP", 0.0));
  optionsParameters.push_back(Util::Param("SENSDEBUGLEVEL", 0));
  optionsParameters.push_back(Util::Param("NUMGAINSCALEBLOCKS", 1));
  optionsParameters.push_back(Util::Param("STAGGERGAINSCALE", 0));
  optionsParameters.push_back(Util::Param("RANDOMIZEVGSTCONST", 0));
  optionsParameters.push_back(Util::Param("NEWEXCESSPHASE", 1));
  optionsParameters.push_back(Util::Param("EXCESSPHASESCALAR1", 1.0));
  optionsParameters.push_back(Util::Param("EXCESSPHASESCALAR2", 1.0));
  optionsParameters.push_back(Util::Param("RANDOMSEED", 0));
  optionsParameters.push_back(Util::Param("TRYTOCOMPACT", false));
  optionsParameters.push_back(Util::Param("CALCULATEALLLEADCURRENTS", false));
  optionsParameters.push_back(Util::Param("NEWMEYER", false ));
  optionsParameters.push_back(Util::Param("ZERORESISTANCETOL", 1.0e-100 ));
  optionsParameters.push_back(Util::Param("CHECKFORZERORESISTANCE", true ));
  optionsParameters.push_back(Util::Param("DETAILED_DEVICE_COUNTS", false ));
  optionsMetadata_[std::string("DEVICE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("OBJFUNC", "V+1"));
  optionsParameters.push_back(Util::Param("PARAM", "VECTOR"));
  optionsMetadata_[std::string("SENS")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 0));
  optionsParameters.push_back(Util::Param("ADJOINT", 0));
  optionsParameters.push_back(Util::Param("DIRECT", 0));
  optionsParameters.push_back(Util::Param("OUTPUTSCALED", 0));
  optionsParameters.push_back(Util::Param("OUTPUTUNSCALED", 1));
  optionsParameters.push_back(Util::Param("STDOUTPUT", 1));
  optionsParameters.push_back(Util::Param("DIAGNOSTICFILE", 0));
  optionsParameters.push_back(Util::Param("DAKOTAFILE", 0));
  optionsParameters.push_back(Util::Param("DIFFERENCE", 0));
  optionsParameters.push_back(Util::Param("SQRTETA", 1.0e-8));
  optionsMetadata_[std::string("SENSITIVITY")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("METHOD", 1));
#ifdef Xyce_DEBUG_ANALYSIS
  optionsParameters.push_back(Util::Param("CONSTSTEP", 0));
#endif
  optionsParameters.push_back(Util::Param("USEDEVICEMAX", 1));
  optionsParameters.push_back(Util::Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(Util::Param("RESTARTSTEPSCALE", .005));
//  optionsParameters.push_back(Util::Param("NLNEARCONV", 1));
  optionsParameters.push_back(Util::Param("NLNEARCONV", 0));
  optionsParameters.push_back(Util::Param("NLSMALLUPDATE", 1));
  optionsParameters.push_back(Util::Param("DOUBLEDCOPSTEP", 0));
  optionsParameters.push_back(Util::Param("FIRSTDCOPSTEP", 0));
  optionsParameters.push_back(Util::Param("LASTDCOPSTEP", 1));
  optionsParameters.push_back(Util::Param("RESETTRANNLS", 1));
  optionsParameters.push_back(Util::Param("BPENABLE", 1));
  optionsParameters.push_back(Util::Param("EXITTIME", 0.0));
  optionsParameters.push_back(Util::Param("EXITSTEP", 0));
  optionsParameters.push_back(Util::Param("ERROPTION", 0));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 0));
  optionsParameters.push_back(Util::Param("JACLIMITFLAG", 0));
  optionsParameters.push_back(Util::Param("JACLIMIT", 1.0e17));
  optionsParameters.push_back(Util::Param("DAESTATEDERIV", 0));
  optionsParameters.push_back(Util::Param("TESTFIRSTSTEP", 0));
  optionsParameters.push_back(Util::Param("DTMIN", 0.0));
  optionsParameters.push_back(Util::Param("NEWBPSTEPPING", 0));
  optionsParameters.push_back(Util::Param("MINTIMESTEPSBP", 10));
  optionsParameters.push_back(Util::Param("NEWLTE", 0));
  optionsParameters.push_back(Util::Param("MAXORD", 5));
  optionsParameters.push_back(Util::Param("MINORD", 1));
  optionsParameters.push_back(Util::Param("OUTPUTINTERPMPDE", 1));
  optionsParameters.push_back(Util::Param("INTERPOUTPUT", 1));
  optionsParameters.push_back(Util::Param("CONDTEST", 0));
  optionsParameters.push_back(Util::Param("CONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(Util::Param("ISOCONDTEST", 0));
  optionsParameters.push_back(Util::Param("ISOCONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(Util::Param("PASSNLSTALL", false));
  optionsParameters.push_back(Util::Param("NLMIN", 3));
  optionsParameters.push_back(Util::Param("NLMAX", 8));
  optionsParameters.push_back(Util::Param("DELMAX", 1.0e+99));
  optionsParameters.push_back(Util::Param("TIMESTEPSREVERSAL", false));
  optionsParameters.push_back(Util::Param("MINTIMESTEPRECOVERY", 0));
  optionsParameters.push_back(Util::Param("FASTTESTS", false));
  optionsParameters.push_back(Util::Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("CURRZEROTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("HISTORYTRACKINGDEPTH", 50));
  optionsMetadata_[std::string("TIMEINT")] = optionsParameters;

  // Make a copy for MPDE time integration.  This copy will result in MPDE having
  // the same defaults, so if we want to change defaults later, we'll have to
  // set this up explicitly.
  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("METHOD", 1));
#ifdef Xyce_DEBUG_ANALYSIS
  optionsParameters.push_back(Util::Param("CONSTSTEP", 0));
#endif
  optionsParameters.push_back(Util::Param("USEDEVICEMAX", 1));
  optionsParameters.push_back(Util::Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(Util::Param("RESTARTSTEPSCALE", .005));
//  optionsParameters.push_back(Util::Param("NLNEARCONV", 1));
  optionsParameters.push_back(Util::Param("NLNEARCONV", 0));
  optionsParameters.push_back(Util::Param("NLSMALLUPDATE", 1));
  optionsParameters.push_back(Util::Param("DOUBLEDCOPSTEP", 0));
  optionsParameters.push_back(Util::Param("FIRSTDCOPSTEP", 0));
  optionsParameters.push_back(Util::Param("LASTDCOPSTEP", 1));
  optionsParameters.push_back(Util::Param("RESETTRANNLS", 1));
  optionsParameters.push_back(Util::Param("BPENABLE", 1));
  optionsParameters.push_back(Util::Param("EXITTIME", 0.0));
  optionsParameters.push_back(Util::Param("EXITSTEP", 0));
  optionsParameters.push_back(Util::Param("ERROPTION", 0));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 0));
  optionsParameters.push_back(Util::Param("JACLIMITFLAG", 0));
  optionsParameters.push_back(Util::Param("JACLIMIT", 1.0e17));
  optionsParameters.push_back(Util::Param("DAESTATEDERIV", 0));
  optionsParameters.push_back(Util::Param("TESTFIRSTSTEP", 0));
  optionsParameters.push_back(Util::Param("DTMIN", 0.0));
  optionsParameters.push_back(Util::Param("MAXORD", 5));
  optionsParameters.push_back(Util::Param("MINORD", 1));
  optionsParameters.push_back(Util::Param("OUTPUTINTERPMPDE", 1));
  optionsParameters.push_back(Util::Param("INTERPOUTPUT", 1));
  optionsParameters.push_back(Util::Param("CONDTEST", 0));
  optionsParameters.push_back(Util::Param("CONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(Util::Param("ISOCONDTEST", 0));
  optionsParameters.push_back(Util::Param("ISOCONDTESTDEVICENAME", "dev_name"));
  optionsParameters.push_back(Util::Param("MINTIMESTEPSBP", 10));
  optionsParameters.push_back(Util::Param("NLMIN", 3));
  optionsParameters.push_back(Util::Param("NLMAX", 8));
  optionsParameters.push_back(Util::Param("DELMAX", 1.0e+99));
  optionsParameters.push_back(Util::Param("TIMESTEPSREVERSAL", false));
  optionsMetadata_[std::string("TIMEINT-MPDE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 1));
  optionsMetadata_[std::string("CONDUCTANCE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("NLSTRATEGY", 0));
  optionsParameters.push_back(Util::Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(Util::Param("NOX", 1));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0E-12));
  optionsParameters.push_back(Util::Param("RELTOL", 1.0E-3));
  optionsParameters.push_back(Util::Param("DELTAXTOL", 1.0));
  optionsParameters.push_back(Util::Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("RHSTOL", 1.0E-6));
  optionsParameters.push_back(Util::Param("MAXSTEP", 200));
  optionsParameters.push_back(Util::Param("MAXSEARCHSTEP", 0));
  optionsParameters.push_back(Util::Param("NORMLVL", 2));
  optionsParameters.push_back(Util::Param("LINOPT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTMAX", Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTMIN", -Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(Util::Param("IN_FORCING", 0));
  optionsParameters.push_back(Util::Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(Util::Param("DLSDEBUG", 0));
  optionsParameters.push_back(Util::Param("MATRIXMARKET", 0));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(Util::Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(Util::Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(Util::Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(Util::Param("MEMORY", 400));
  optionsParameters.push_back(Util::Param("CONTINUATION", 0));
  optionsParameters.push_back(Util::Param("ENFORCEDEVICECONV", 1));
  optionsParameters.push_back(Util::Param("FASTTESTS", false));
  optionsParameters.push_back(Util::Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[std::string("NONLIN")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("NLSTRATEGY", 0));
  optionsParameters.push_back(Util::Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(Util::Param("NOX", 1));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(Util::Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(Util::Param("DELTAXTOL", 0.33));
  optionsParameters.push_back(Util::Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("RHSTOL", 1.0E-2));
  optionsParameters.push_back(Util::Param("MAXSTEP", 20));
  optionsParameters.push_back(Util::Param("MAXSEARCHSTEP", 2));
  optionsParameters.push_back(Util::Param("NORMLVL", 2));
  optionsParameters.push_back(Util::Param("LINOPT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTMAX", Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTMIN", -Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(Util::Param("IN_FORCING", 0));
  optionsParameters.push_back(Util::Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(Util::Param("DLSDEBUG", 0));
  optionsParameters.push_back(Util::Param("MATRIXMARKET", 0));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(Util::Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(Util::Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(Util::Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(Util::Param("MEMORY", 400));
  optionsParameters.push_back(Util::Param("CONTINUATION", 0));
  optionsParameters.push_back(Util::Param("ENFORCEDEVICECONV", 0));
  optionsParameters.push_back(Util::Param("FASTTESTS", false));
  optionsParameters.push_back(Util::Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[std::string("NONLIN-TRAN")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("NLSTRATEGY", 0));
  optionsParameters.push_back(Util::Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(Util::Param("NOX", 0));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0E-9));
  optionsParameters.push_back(Util::Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(Util::Param("DELTAXTOL", 0.33));
  optionsParameters.push_back(Util::Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("RHSTOL", 1.0E-6));
  optionsParameters.push_back(Util::Param("MAXSTEP", 200));
  optionsParameters.push_back(Util::Param("MAXSEARCHSTEP", 2));
  optionsParameters.push_back(Util::Param("NORMLVL", 2));
  optionsParameters.push_back(Util::Param("LINOPT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTMAX", Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTMIN", -Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(Util::Param("IN_FORCING", 0));
  optionsParameters.push_back(Util::Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(Util::Param("DLSDEBUG", 0));
  optionsParameters.push_back(Util::Param("MATRIXMARKET", 0));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(Util::Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(Util::Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(Util::Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(Util::Param("MEMORY", 400));
  optionsParameters.push_back(Util::Param("CONTINUATION", 0));
  optionsParameters.push_back(Util::Param("ENFORCEDEVICECONV", 0));
  optionsMetadata_[std::string("NONLIN-HB")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("NLSTRATEGY", 0));
  optionsParameters.push_back(Util::Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(Util::Param("NOX", 1));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0E-12));
  optionsParameters.push_back(Util::Param("RELTOL", 1.0E-3));
  optionsParameters.push_back(Util::Param("DELTAXTOL", 1.0));
  optionsParameters.push_back(Util::Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("RHSTOL", 1.0E-6));
  optionsParameters.push_back(Util::Param("MAXSTEP", 200));
  optionsParameters.push_back(Util::Param("MAXSEARCHSTEP", 0));
  optionsParameters.push_back(Util::Param("NORMLVL", 2));
  optionsParameters.push_back(Util::Param("LINOPT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTMAX", Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTMIN", -Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(Util::Param("IN_FORCING", 0));
  optionsParameters.push_back(Util::Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(Util::Param("DLSDEBUG", 0));
  optionsParameters.push_back(Util::Param("MATRIXMARKET", 0));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(Util::Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(Util::Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(Util::Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(Util::Param("MEMORY", 400));
  optionsParameters.push_back(Util::Param("CONTINUATION", 0));
  optionsParameters.push_back(Util::Param("ENFORCEDEVICECONV", 0));
  optionsParameters.push_back(Util::Param("ALGORITHM", 0));
  optionsParameters.push_back(Util::Param("MAXCONTSTEPS", 0));
  optionsParameters.push_back(Util::Param("CONTINUATIONFLAG", 1));
  optionsParameters.push_back(Util::Param("INNERFAIL", 1));
  optionsParameters.push_back(Util::Param("EXITWITHFAILURE", 1));
  optionsParameters.push_back(Util::Param("FULLNEWTONENFORCE", 1));
  optionsParameters.push_back(Util::Param("CONPARAM", "VA:V0"));
  optionsParameters.push_back(Util::Param("VOLTLIMTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("FASTTESTS", false));
  optionsParameters.push_back(Util::Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[std::string("NONLIN-TWOLEVEL")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("NLSTRATEGY", 0));
  optionsParameters.push_back(Util::Param("SEARCHMETHOD", 0));
  optionsParameters.push_back(Util::Param("NOX", 1));
  optionsParameters.push_back(Util::Param("ABSTOL", 1.0E-6));
  optionsParameters.push_back(Util::Param("RELTOL", 1.0E-2));
  optionsParameters.push_back(Util::Param("DELTAXTOL", 0.33));
  optionsParameters.push_back(Util::Param("SMALLUPDATETOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("RHSTOL", 1.0E-2));
  optionsParameters.push_back(Util::Param("MAXSTEP", 20));
  optionsParameters.push_back(Util::Param("MAXSEARCHSTEP", 2));
  optionsParameters.push_back(Util::Param("NORMLVL", 2));
  optionsParameters.push_back(Util::Param("LINOPT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTBT", 0));
  optionsParameters.push_back(Util::Param("CONSTRAINTMAX", Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTMIN", -Util::MachineDependentParams::DoubleMax()));
  optionsParameters.push_back(Util::Param("CONSTRAINTCHANGE", 0.0));
  optionsParameters.push_back(Util::Param("IN_FORCING", 0));
  optionsParameters.push_back(Util::Param("AZ_TOL", 1.0E-12));
  optionsParameters.push_back(Util::Param("DLSDEBUG", 0));
  optionsParameters.push_back(Util::Param("MATRIXMARKET", 0));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 1));
  optionsParameters.push_back(Util::Param("DEBUGMINTIMESTEP", 0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIMESTEP", 99999999));
  optionsParameters.push_back(Util::Param("DEBUGMINTIME", 0.0));
  optionsParameters.push_back(Util::Param("DEBUGMAXTIME", 1.0E99));
  optionsParameters.push_back(Util::Param("SCREENOUTPUT", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEPTYPE", 0));
  optionsParameters.push_back(Util::Param("RECOVERYSTEP", 1.0));
  optionsParameters.push_back(Util::Param("MEMORY", 400));
  optionsParameters.push_back(Util::Param("CONTINUATION", 0));
  optionsParameters.push_back(Util::Param("ENFORCEDEVICECONV", 0));
  optionsParameters.push_back(Util::Param("ALGORITHM", 0));
  optionsParameters.push_back(Util::Param("MAXCONTSTEPS", 0));
  optionsParameters.push_back(Util::Param("CONTINUATIONFLAG", 1));
  optionsParameters.push_back(Util::Param("INNERFAIL", 1));
  optionsParameters.push_back(Util::Param("EXITWITHFAILURE", 1));
  optionsParameters.push_back(Util::Param("FULLNEWTONENFORCE", 1));
  optionsParameters.push_back(Util::Param("CONPARAM", "VA:V0"));
  optionsParameters.push_back(Util::Param("VOLTLIMTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("FASTTESTS", false));
  optionsParameters.push_back(Util::Param("VOLTZEROTOL", 1.0e-6));
  optionsParameters.push_back(Util::Param("CURRZEROTOL", 1.0e-6));
  optionsMetadata_[std::string("NONLIN-TWOLEVEL-TRAN")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("STEPPER", "NATURAL"));
  optionsParameters.push_back(Util::Param("PREDICTOR", "CONSTANT"));
  optionsParameters.push_back(Util::Param("STEPCONTROL", "CONSTANT"));

  optionsParameters.push_back(Util::Param("CONPARAM", "VECTOR"));
  optionsParameters.push_back(Util::Param("INITIALVALUE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MAXVALUE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MINVALUE", "VECTOR"));
  optionsParameters.push_back(Util::Param("INITIALSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MAXSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MINSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(Util::Param("AGGRESSIVENESS", "VECTOR"));
  optionsParameters.push_back(Util::Param("BIFPARAM", "VA:V0"));
  optionsParameters.push_back(Util::Param("MAXSTEPS", 20));
  optionsParameters.push_back(Util::Param("MAXNLITERS", 20));
  optionsParameters.push_back(Util::Param("PARAMLIST", "VECTOR"));
  optionsParameters.push_back(Util::Param("VOLTAGELIST", "DOFS"));
  optionsParameters.push_back(Util::Param("VOLTAGESCALEFACTOR", 1.0));
  optionsParameters.push_back(Util::Param("RESIDUALCONDUCTANCE", 0.0));
  optionsMetadata_[std::string("LOCA")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("STEPPER", "NATURAL"));
  optionsParameters.push_back(Util::Param("PREDICTOR", "CONSTANT"));
  optionsParameters.push_back(Util::Param("STEPCONTROL", "CONSTANT"));
  optionsParameters.push_back(Util::Param("CONPARAM", "VECTOR"));
  optionsParameters.push_back(Util::Param("INITIALVALUE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MAXVALUE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MINVALUE", "VECTOR"));
  optionsParameters.push_back(Util::Param("INITIALSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MAXSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(Util::Param("MINSTEPSIZE", "VECTOR"));
  optionsParameters.push_back(Util::Param("AGGRESSIVENESS", "VECTOR"));
  optionsParameters.push_back(Util::Param("BIFPARAM", "VA:V0"));
  optionsParameters.push_back(Util::Param("MAXSTEPS", 20));
  optionsParameters.push_back(Util::Param("MAXNLITERS", 20));
  optionsParameters.push_back(Util::Param("PARAMLIST", "VECTOR"));
  optionsParameters.push_back(Util::Param("VOLTAGELIST", "DOFS"));
  optionsParameters.push_back(Util::Param("VOLTAGESCALEFACTOR", 1.0));
  optionsParameters.push_back(Util::Param("RESIDUALCONDUCTANCE", 0.0));
  optionsMetadata_[std::string("TWOLEVEL-LOCA")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("AZ_max_iter", 500));
  optionsParameters.push_back(Util::Param("AZ_precond", 14));
  optionsParameters.push_back(Util::Param("AZ_solver", 1));
  optionsParameters.push_back(Util::Param("AZ_conv", 0));
  optionsParameters.push_back(Util::Param("AZ_pre_calc", 1));
  optionsParameters.push_back(Util::Param("AZ_keep_info", 1));
  optionsParameters.push_back(Util::Param("AZ_orthog", 1));
  optionsParameters.push_back(Util::Param("AZ_subdomain_solve", 9));
  optionsParameters.push_back(Util::Param("AZ_ilut_fill", 3.0));
  optionsParameters.push_back(Util::Param("AZ_drop", 1.0E-3));
  optionsParameters.push_back(Util::Param("AZ_reorder", 0));
  optionsParameters.push_back(Util::Param("AZ_scaling", 0));
  optionsParameters.push_back(Util::Param("AZ_kspace", 500));
  optionsParameters.push_back(Util::Param("AZ_tol", 1.0E-12));
  optionsParameters.push_back(Util::Param("AZ_output", 0));
  optionsParameters.push_back(Util::Param("AZ_diagnostics", 0));
  optionsParameters.push_back(Util::Param("AZ_overlap", 0));
  optionsParameters.push_back(Util::Param("AZ_rthresh", 1.0001));
  optionsParameters.push_back(Util::Param("AZ_athresh", 1.0E-4));
  optionsParameters.push_back(Util::Param("AZ_filter", 0.0));
  optionsParameters.push_back(Util::Param("TR_filter", 0.0));
  optionsParameters.push_back(Util::Param("TR_partition", 1));
#ifdef Xyce_SHYLU
  optionsParameters.push_back(Util::Param("ShyLU_rthresh", 1.0E-3));
#endif
#ifdef Xyce_USE_ISORROPIA
  optionsParameters.push_back(Util::Param("TR_partition_type", "GRAPH"));
#endif
  optionsParameters.push_back(Util::Param("TR_reindex", 1));
  optionsParameters.push_back(Util::Param("TR_solvermap", 1));
  optionsParameters.push_back(Util::Param("TR_amd", 1));
  optionsParameters.push_back(Util::Param("TR_btf", 0));
  optionsParameters.push_back(Util::Param("TR_global_btf", 0));
  optionsParameters.push_back(Util::Param("TR_global_btf_droptol", 1.0E-16));
  optionsParameters.push_back(Util::Param("TR_global_btf_verbose", 0));
#ifdef Xyce_TRILINOS_DEV
  optionsParameters.push_back(Util::Param("TR_global_amd", 0));
  optionsParameters.push_back(Util::Param("TR_global_amd_verbose", 0));
#endif
  optionsParameters.push_back(Util::Param("TR_singleton_filter", 0));
  optionsParameters.push_back(Util::Param("SLU_EQUILIBRATE", 1));
  optionsParameters.push_back(Util::Param("SLU_REFACTOR", 1));
  optionsParameters.push_back(Util::Param("SLU_PERMUTE", 2));
  optionsParameters.push_back(Util::Param("SLU_PIVOT_THRESH", -1.0));
  optionsParameters.push_back(Util::Param("SLU_FILL_FAC", -1));
  optionsParameters.push_back(Util::Param("BTF", 0));
  optionsParameters.push_back(Util::Param("BTF_VERBOSE", 0));
  optionsParameters.push_back(Util::Param("BTF_ATHRESH", 0.0));
  optionsParameters.push_back(Util::Param("BTF_RTHRESH", 0.0));
  optionsParameters.push_back(Util::Param("BTF_RTHRESH_INIT", 0.0));
  optionsParameters.push_back(Util::Param("BTF_INIT", 0));
  optionsParameters.push_back(Util::Param("BTF_THRESHOLDING", 0));
  optionsParameters.push_back(Util::Param("BTF_RNTHRESHFAC", 1.0e-3));
  optionsParameters.push_back(Util::Param("adaptive_solve", 0));
  optionsParameters.push_back(Util::Param("use_aztec_precond", 1));
  optionsParameters.push_back(Util::Param("use_ifpack_factory", 0));
  optionsParameters.push_back(Util::Param("ifpack_type", "Amesos"));
  optionsParameters.push_back(Util::Param("diag_perturb", 0.0));
#ifdef Xyce_ML
  optionsParameters.push_back(Util::Param("ML_max_level", 5));
#endif
  optionsParameters.push_back(Util::Param("TR_rcm", 0));
  optionsParameters.push_back(Util::Param("TR_scale", 0));
  optionsParameters.push_back(Util::Param("TR_scale_left", 0));
  optionsParameters.push_back(Util::Param("TR_scale_right", 0));
  optionsParameters.push_back(Util::Param("TR_scale_exp", 1.0));
  optionsParameters.push_back(Util::Param("TR_scale_iter", 0));
  optionsParameters.push_back(Util::Param("TYPE", "DEFAULT"));
  optionsParameters.push_back(Util::Param("PREC_TYPE", "DEFAULT"));
#ifdef Xyce_BELOS
  optionsParameters.push_back(Util::Param("BELOS_SOLVER_TYPE", "Block GMRES"));
#endif
  optionsParameters.push_back(Util::Param("KLU_REPIVOT", 1));
  optionsParameters.push_back(Util::Param("KLU_REINDEX", 0));
  optionsParameters.push_back(Util::Param("OUTPUT_LS", 1));
  optionsParameters.push_back(Util::Param("OUTPUT_BASE_LS", 1));
  optionsParameters.push_back(Util::Param("OUTPUT_FAILED_LS", 1));
  optionsMetadata_[std::string("LINSOL")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("AZ_max_iter", 500));
  optionsParameters.push_back(Util::Param("AZ_solver", 1));
  optionsParameters.push_back(Util::Param("AZ_conv", 0));
  optionsParameters.push_back(Util::Param("AZ_pre_calc", 1));
  optionsParameters.push_back(Util::Param("AZ_keep_info", 1));
  optionsParameters.push_back(Util::Param("AZ_orthog", 1));
  optionsParameters.push_back(Util::Param("AZ_reorder", 0));
  optionsParameters.push_back(Util::Param("AZ_scaling", 0));
  optionsParameters.push_back(Util::Param("AZ_kspace", 500));
  optionsParameters.push_back(Util::Param("AZ_tol", 1.0E-12));
  optionsParameters.push_back(Util::Param("AZ_output", 0));
  optionsParameters.push_back(Util::Param("AZ_diagnostics", 0));
  optionsParameters.push_back(Util::Param("TYPE", "DEFAULT"));
  optionsParameters.push_back(Util::Param("PREC_TYPE", "NONE"));
  optionsMetadata_[std::string("LINSOL-HB")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("PARTITIONER", "CHACO"));
  optionsParameters.push_back(Util::Param("DISTRIBINDSRCNODES", 0.0));
  optionsParameters.push_back(Util::Param("USE_WEIGHTS", 0));
  optionsParameters.push_back(Util::Param("METHOD", "KWAY"));
  optionsParameters.push_back(Util::Param("USE_VWEIGHTS", 0));
  optionsParameters.push_back(Util::Param("USE_EWEIGHTS", 0));
  optionsParameters.push_back(Util::Param("BALANCE_FACTOR", 0));
  optionsParameters.push_back(Util::Param("BISECTIONS", 10));
  optionsParameters.push_back(Util::Param("VERTEX_GROUPING", 1));
  optionsParameters.push_back(Util::Param("REFINEMENT", 1));
  optionsParameters.push_back(Util::Param("VCYCLE", 1));
  optionsParameters.push_back(Util::Param("RECONST", 0));
  optionsParameters.push_back(Util::Param("PREASSIGN", 0));
  optionsParameters.push_back(Util::Param("SEED", -1));
  optionsParameters.push_back(Util::Param("DEBUG", 0));
  optionsMetadata_[std::string("PARALLEL")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("PACK", 1));
  optionsParameters.push_back(Util::Param("JOB", ""));
  optionsParameters.push_back(Util::Param("START_TIME", 0.0));
  optionsParameters.push_back(Util::Param("FILE", ""));
  optionsParameters.push_back(Util::Param("INITIAL_INTERVAL", 0.0));
  optionsParameters.push_back(Util::Param("TIME", 0.0));
  optionsParameters.push_back(Util::Param("INTERVAL", 0.0));
  optionsMetadata_[std::string("RESTART")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("INITIAL_INTERVAL", 0.0));
  optionsParameters.push_back(Util::Param("TIME", 0.0));
  optionsParameters.push_back(Util::Param("INTERVAL", 0.0));
  optionsParameters.push_back(Util::Param("HDF5FILENAME", ""));
  optionsParameters.push_back(Util::Param("PRINTENDOFSIMLINE", true));
  optionsParameters.push_back(Util::Param("OUTPUTVERSIONINRAWFILE", false));
  optionsMetadata_[std::string("OUTPUT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("INPUT", ""));
  optionsParameters.push_back(Util::Param("OUTPUT", ""));
  optionsMetadata_[std::string("OP_IO")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("TIME", 0.0));
  optionsParameters.push_back(Util::Param("INTERVAL", 0.0));
  optionsMetadata_[std::string("OUTPUT-LINE")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("TYPE", "TRAN"));
  optionsParameters.push_back(Util::Param("FILE", ""));
  optionsParameters.push_back(Util::Param("FORMAT", "STD"));
  optionsParameters.push_back(Util::Param("DELIMITER", ""));
  optionsParameters.push_back(Util::Param("WIDTH", 17));
  optionsParameters.push_back(Util::Param("PRECISION", 8));
  optionsParameters.push_back(Util::Param("TIMESCALEFACTOR", 1.0));
  optionsParameters.push_back(Util::Param("FILTER", 0.0));
  optionsMetadata_[std::string("PRINT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("REPLICATED_CKT", 1));
  optionsParameters.push_back(Util::Param("CHECK_CONNECTIVITY", 0));
  optionsParameters.push_back(Util::Param("SUPERNODE", false));
  optionsParameters.push_back(Util::Param("OUTPUTNAMESFILE", false));
  optionsMetadata_[std::string("TOPOLOGY")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("TEST", 0));
  optionsParameters.push_back(Util::Param("OSCSRC", "VECTOR"));
  optionsParameters.push_back(Util::Param("OSCOUT", ""));
  optionsParameters.push_back(Util::Param("AUTON2", false));
  optionsParameters.push_back(Util::Param("AUTON2MAX", 100));
  optionsParameters.push_back(Util::Param("STARTUPPERIODS", 0));
  optionsParameters.push_back(Util::Param("SAVEICDATA", false));
  optionsParameters.push_back(Util::Param("N2", 10));
  optionsParameters.push_back(Util::Param("T2", 0.0));
  optionsParameters.push_back(Util::Param("IC", 0));
  optionsParameters.push_back(Util::Param("ICPER", 10));
  optionsParameters.push_back(Util::Param("DIFF", 0));
  optionsParameters.push_back(Util::Param("DIFFORDER", 1));
  optionsParameters.push_back(Util::Param("FREQDOMAIN", 0));
  optionsParameters.push_back(Util::Param("WAMPDE", 0));
  optionsParameters.push_back(Util::Param("DCOPEXIT", 0));
  optionsParameters.push_back(Util::Param("ICEXIT", 0));
  optionsParameters.push_back(Util::Param("EXITSAWTOOTHSTEP", -1));
  optionsParameters.push_back(Util::Param("PHASE", 0));
  optionsParameters.push_back(Util::Param("PHASECOEFF", 0));
  optionsParameters.push_back(Util::Param("NONLTESTEPS", 10));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 0));
  optionsMetadata_[std::string("MPDEINT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("TEST", false));
  optionsParameters.push_back(Util::Param("NUMFREQ", 21));
  optionsParameters.push_back(Util::Param("STARTUPPERIODS", 0));
  optionsParameters.push_back(Util::Param("SAVEICDATA", false));
  optionsParameters.push_back(Util::Param("DEBUGLEVEL", 0));
  optionsParameters.push_back(Util::Param("TAHB", 1));
  optionsParameters.push_back(Util::Param("VOLTLIM", 1));
  optionsMetadata_[std::string("HBINT")] = optionsParameters;

  optionsParameters.clear();
  optionsParameters.push_back(Util::Param("METHOD", "PRIMA"));
  optionsParameters.push_back(Util::Param("SAVEREDSYS", false));
  optionsParameters.push_back(Util::Param("COMPORIGTF", false));
  optionsParameters.push_back(Util::Param("COMPREDTF", false));
  optionsParameters.push_back(Util::Param("COMPTYPE", "DEC"));
  optionsParameters.push_back(Util::Param("COMPNP", 10));
  optionsParameters.push_back(Util::Param("COMPFSTART", 1.0));
  optionsParameters.push_back(Util::Param("COMPFSTOP", 1.0));
  optionsParameters.push_back(Util::Param("EXPPOINT", 0.0));
  optionsParameters.push_back(Util::Param("SCALETYPE", 0));
  optionsParameters.push_back(Util::Param("SCALEFACTOR", 1));
  optionsParameters.push_back(Util::Param("SCALEFACTOR1", 0.01));
  optionsParameters.push_back(Util::Param("SPARSIFICATIONTYPE", 0));
  optionsMetadata_[std::string("MOR_OPTS")] = optionsParameters;
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
  std::vector<Device::Param> sourceFcnParamList;
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
// Function      : CircuitMetadata::getSourceFunctionID
// Purpose       : Return the integer value corresponding to the
//                 given source function. Return NUM_SRC_DATA if the given string
//                 is unknown.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
int CircuitMetadata::getSourceFunctionID(const std::string & sourceFcn) const
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
// Function      : CircuitMetadata::getInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
std::vector<Device::Param> &CircuitMetadata::getInstanceParameters(
  const std::string &   deviceType,
  int                   modelLevel)
{
  if (modelLevel == -1 || deviceMetadataIndex.find(NameLevelKey(deviceType, modelLevel)) == deviceMetadataIndex.end())
    return getDeviceMetadata(deviceType, modelLevel).instanceParameters;
  else
    return deviceMetadata_[deviceMetadataIndex[NameLevelKey(deviceType, modelLevel)]].instanceParameters;
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
  const std::string &                deviceType,
  const std::string &                parameterName, int modelLevel,
  std::vector<Device::Param> &       components)
{
  DeviceMetadata &device_metadata = getDeviceMetadata(deviceType, modelLevel);
  DeviceParamMap &icpMap = device_metadata.instanceCompositeParameterMap;
  DeviceParamMap::iterator iterIcp = icpMap.find(parameterName);

  if ( iterIcp != icpMap.end() )
  {
    components = iterIcp->second;
  }
  else
  {
    Report::UserError() << "There are no component parameters in metadata for the VECTOR-COMPOSITE parameter: " << parameterName;
  }
}

//-----------------------------------------------------------------------------
// Function      : CircuitMetadata::getModelParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
std::vector<Device::Param> &CircuitMetadata::getModelParameters(
  const std::string &   modelType,
  int                   modelLevel)
{
  if (modelLevel == -1 || deviceMetadataIndex.find(NameLevelKey(modelType, modelLevel)) == deviceMetadataIndex.end())
    return getDeviceMetadata(modelType, modelLevel).modelParameters;
  else
    return deviceMetadata_[deviceMetadataIndex[NameLevelKey(modelType, modelLevel)]].modelParameters;
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
  const std::string &                modelType,
  const std::string &                parameterName, int modelLevel,
  std::vector<Device::Param> &       components)
{
  DeviceMetadata &device_metadata = getDeviceMetadata(modelType, modelLevel);
  DeviceParamMap &mcpMap = device_metadata.modelCompositeParameterMap;

  DeviceParamMap::iterator iterIcp = mcpMap.find(parameterName);

  if ( iterIcp != mcpMap.end() )
  {
    components = iterIcp->second;
  }
  else
  {
    Report::UserError() << "There are no component parameters in metadata for the VECTOR-COMPOSITE parameter " << parameterName;
  }
}

} // namespace IO
} // namespace Xyce
