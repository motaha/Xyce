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
// Filename      : N_IO_CircuitMetadata.h
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 12/11/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.22.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef DeviceMetadata_H
#define DeviceMetadata_H

// ---------- Standard Includes ----------

#include <list>
#include <map>
#include <string>
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceLevelKey.h>
#include <N_DEV_Param.h>
#include <N_UTL_Misc.h>
#include <N_UTL_NoCase.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : DeviceMetadata
// Purpose       :
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/11/00
//-----------------------------------------------------------------------------

class DeviceMetadata
{
public:
  typedef Device::DeviceParamMap DeviceParamMap;

  DeviceMetadata()
    : deviceType(""),
      level(0),
      numNodes(0),
      numOptionalNodes(0),
      numFillNodes(0),
      modelRequired(0),
      primaryParameter("") {};

  DeviceMetadata(const std::string & dtype, int lev)
    : deviceType(dtype),
      level(lev),
      numNodes(0),
      numOptionalNodes(0),
      numFillNodes(0),
      modelRequired(0),
      primaryParameter("") {};

  std::string deviceType;
  int level;
  int numNodes;
  int numOptionalNodes;
  int numFillNodes;
  int modelRequired;
  std::string primaryParameter;
  std::vector<std::string> modelTypes;
  std::vector<Device::Param> instanceParameters;
  std::vector<Device::Param> modelParameters;
  DeviceParamMap instanceCompositeParameterMap;
  DeviceParamMap modelCompositeParameterMap;
};

//-----------------------------------------------------------------------------
// Class         : CircuitMetadata
// Purpose       :
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/11/00
//-----------------------------------------------------------------------------

class CircuitMetadata
{
public:
  typedef Device::DeviceParamMap DeviceParamMap;
  typedef Device::DeviceLevelKey DeviceLevelKey;
  typedef std::map<DeviceLevelKey, DeviceLevelKey, Device::DeviceLevelLess> DeviceMetadataIndexMap;
  typedef std::map<DeviceLevelKey, DeviceMetadata, Device::DeviceLevelLess> DeviceMetadataMap;

  CircuitMetadata()
    : metadataLocation("")
  {};

  // Destructor
  ~CircuitMetadata()
  {};

  // Build device instance and model metadata.
  void buildMetadata( );

  // Find the metadata for a given device. If not found, return NULL.
  DeviceMetadata* findDeviceMetadata(
    const std::string & deviceType,
    int level);

  // Determine if the given model type is valid for the given device.
  bool isModelTypeValid(
    const std::string & deviceType,
    const std::string & modelType,
    int modelLevel);

  // Determine if the given parameter name is a valid parameter for
  // the given device.
  bool isDeviceParameter(const std::string & deviceType, int modelLevel, const std::string & parameterName);

  std::string getPrimaryParameter(const std::string & deviceType, int modelLevel);

  int getNumberOfNodes(const std::string & deviceType, int modelLevel);

  int getNumberOfOptionalNodes(const std::string & deviceType, int modelLevel);

  int getNumberOfFillNodes(const std::string & deviceType, int modelLevel);

  bool isModelRequired(const std::string & deviceType, int modelLevel);

  int getSourceFunctionID(const std::string & sourceFcn);

  std::vector<Device::Param> * getPtrToInstanceParameters(
    const std::string & deviceType,
    int modelLevel);

  void getInstanceCompositeComponents(
    const std::string & deviceType,
    const std::string & parameterName, int modelLevel,
    std::vector<Device::Param> & components);

  void getModelCompositeComponents(
    const std::string & modelType,
    const std::string & parameterName, int modelLevel,
    std::vector<Device::Param> & components);

  std::vector<Device::Param> * getPtrToModelParameters(
    const std::string & modelType,
    int modelLevel);

  bool isOptionsPackage(const std::string & optionName);

  std::vector<N_UTL_Param> * getPtrToOptionsParameters(const std::string & optionName);

  void getSourceFunctionParameters(const std::string & sourceFcn, std::vector<Device::Param> & sourceFunctionParameters);

  N_DEV_DeviceInterface* devIntPtr_;

private:
  DeviceMetadata * createMetadataEntry(const std::string & deviceType, int level);
  void optionsMetadata();
  void sourceFunctionMetadata();

  std::string metadataLocation;

  // The device metadata container which contains information
  // about the recognizable devices in a netlist.
  DeviceMetadataMap deviceMetadata_;

  // An index to translate from model line indexes to deviceMetadata indexes
  // For example: NPN_1 to Q_1 for BJT
  DeviceMetadataIndexMap deviceMetadataIndex;

  // The options metadata container.
  std::map<std::string, std::vector<N_UTL_Param>, LessNoCase> optionsMetadata_;

  // Source function metadata.
  DeviceParamMap sourceFcnMap;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CircuitMetadata N_IO_CircuitMetadata;
typedef Xyce::IO::DeviceMetadata N_IO_DeviceMetadata;

#endif
