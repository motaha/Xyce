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
// Revision Number: $Revision: 1.31.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_CircuitMetadata_h
#define Xyce_N_IO_CircuitMetadata_h

#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_DEV_Param.h>
#include <N_UTL_NoCase.h>

namespace Xyce {
namespace IO {

typedef std::map<std::string, std::vector<Device::Param>, LessNoCase> DeviceParamMap;
typedef std::map<std::string, std::vector<Util::Param>, LessNoCase> UtilParamMap;

//-----------------------------------------------------------------------------
// Class         : DeviceMetadata
// Purpose       :
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/11/00
//-----------------------------------------------------------------------------

struct DeviceMetadata
{
  DeviceMetadata()
    : levelValid(false),
      deviceType(),
      level(0),
      numNodes(0),
      numOptionalNodes(0),
      numFillNodes(0),
      modelRequired(0),
      primaryParameter(),
      modelTypes(),
      instanceParameters(),
      modelParameters(),
      instanceCompositeParameterMap(),
      modelCompositeParameterMap()
  {}

  bool isModelTypeValid(const std::string & modelType) const 
  {
    return std::find_if(modelTypes.begin(), modelTypes.end(), EqualNoCasePred(modelType)) != modelTypes.end();
  }

  bool isModelLevelValid() const 
  {
    return levelValid;
  }

  bool                        levelValid;
  std::string                 deviceType;
  int                         level;
  int                         numNodes;
  int                         numOptionalNodes;
  int                         numFillNodes;
  int                         modelRequired;
  std::string                 primaryParameter;
  std::vector<std::string>    modelTypes;
  std::vector<Device::Param>  instanceParameters;
  std::vector<Device::Param>  modelParameters;
  DeviceParamMap              instanceCompositeParameterMap;
  DeviceParamMap              modelCompositeParameterMap;
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
  typedef std::map<NameLevelKey, NameLevelKey, NameLevelLess> DeviceMetadataIndexMap;
  typedef std::map<NameLevelKey, DeviceMetadata, NameLevelLess> DeviceMetadataMap;

  CircuitMetadata()
    : deviceMetadata_(),
      deviceMetadataIndex(),
      optionsMetadata_(),
      sourceFcnMap()
  {}

  ~CircuitMetadata()
  {}

  // Build device instance and model metadata.
  void buildMetadata();

  // Find the metadata for a given device. If not found, return NULL.
  DeviceMetadata &getDeviceMetadata(const std::string & deviceType, int level) const;

  // Determine if the given model type is valid for the given device.
  bool isModelLevelValid(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).isModelLevelValid();
  }

  // Determine if the given model type is valid for the given device.
  bool isModelTypeValid(const std::string & deviceType, const std::string & modelType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).isModelTypeValid(modelType);
  }

  // Determine if the given parameter name is a valid parameter for the given device.
  bool isDeviceParameter(const std::string & deviceType, int modelLevel, const std::string & parameterName) const 
  {
    DeviceMetadata &device_metadata = getDeviceMetadata(deviceType, modelLevel);

    std::vector<Device::Param>::const_iterator it = std::find(device_metadata.instanceParameters.begin(),
                                                              device_metadata.instanceParameters.end(),
                                                              Device::Param(parameterName, ""));

    return it != device_metadata.instanceParameters.end();
  }


  const std::string &getPrimaryParameter(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).primaryParameter;
  }

  int getNumberOfNodes(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).numNodes;
  }

  int getNumberOfOptionalNodes(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).numOptionalNodes;
  }

  int getNumberOfFillNodes(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).numFillNodes;
  }

  bool isModelRequired(const std::string & deviceType, int modelLevel) const 
  {
    return getDeviceMetadata(deviceType, modelLevel).modelRequired;
  }

  int getSourceFunctionID(const std::string & sourceFcn) const;

  std::vector<Device::Param> &getInstanceParameters(const std::string & deviceType, int modelLevel);

  void getInstanceCompositeComponents(const std::string & deviceType, const std::string & parameterName, int modelLevel, std::vector<Device::Param> & components);

  void getModelCompositeComponents(const std::string & modelType, const std::string & parameterName, int modelLevel, std::vector<Device::Param> & components);

  std::vector<Device::Param> &getModelParameters(const std::string & modelType, int modelLevel);

  bool isOptionsPackage(const std::string & optionName) 
  {
    return optionsMetadata_.find(optionName) != optionsMetadata_.end();
  }

  const std::vector<Util::Param> &getOptionsParameters(const std::string & optionName) const 
  {
    return optionsMetadata_[optionName];
  }

  const std::vector<Device::Param> &getSourceFunctionParameters(const std::string & sourceFcn) const 
  {
    return sourceFcnMap[sourceFcn];
  }

private:
  void optionsMetadata();
  void sourceFunctionMetadata();

private:
  // The device metadata container which contains information
  // about the recognizable devices in a netlist.
  mutable DeviceMetadataMap deviceMetadata_;

  // An index to translate from model line indexes to deviceMetadata indexes
  // For example: NPN_1 to Q_1 for BJT
  mutable DeviceMetadataIndexMap deviceMetadataIndex;

  // The options metadata container.
  mutable UtilParamMap optionsMetadata_;

  // Source function metadata.
  mutable DeviceParamMap sourceFcnMap;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CircuitMetadata N_IO_CircuitMetadata;
typedef Xyce::IO::DeviceMetadata N_IO_DeviceMetadata;

#endif // Xyce_N_IO_CircuitMetadata_h
