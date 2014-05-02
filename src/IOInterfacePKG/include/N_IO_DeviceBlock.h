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
// Filename       : $RCSfile: N_IO_DeviceBlock.h,v $
//
// Purpose        : Declare the N_IO_DeviceBlock class instantiations of which
//                  are associated with netlist device lines.
//
// Special Notes  : ERK.  It seems that the name "N_IO_InstanceBlock" would have been
//                  more appropriate and less confusing, or possibly
//                  N_IO_InstanceParametersBlock.  Calling it "device block"
//                  makes it sound much more general than it really is.
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/09/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.68.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_IO_DEVICEBLOCK_H
#define N_IO_DEVICEBLOCK_H

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

#include <N_DEV_Param.h>

#include <N_TOP_NodeDevBlock.h>

namespace Xyce {
namespace IO {

class DeviceBlock
{
private:
  DeviceBlock();

public:
  DeviceBlock( CircuitContext & cc, CircuitMetadata & md );

  DeviceBlock(
     CircuitContext &                                                cc,
     CircuitMetadata &                                               md,
     const std::string &                                             fileName,
     const std::vector<SpiceSeparatedFieldTool::StringToken> &       parsedInputLine);

  DeviceBlock(const DeviceBlock &rhsDB );

private:
  DeviceBlock &operator=(const DeviceBlock &rhsDB );

public:
  ~DeviceBlock()
  {}

  // Setters and Getters
  void setParsedLine(const std::vector<SpiceSeparatedFieldTool::StringToken> &deviceLine)
  {
    parsedLine_ = deviceLine;
  }

  void setFileName(const std::string &fileName)
  {
    netlistFileName_ = fileName;
  }

  void setName(const std::string &name )
  {
    deviceData_.getDevBlock().setName(name);
    deviceData_.getNodeBlock().set_id( name );
  }

  void setNetlistType( char type )
  {
    netlistType_ = type;
  }

  void setNetlistType(const std::string &type)
  {
    netlistType_ = type;
  }

  void addNodeValue( std::string const& nodeValue )
  {
    deviceData_.getNodeBlock().addNode(tagged_param(nodeValue, 0));
    checkNode(nodeValue);
  }

  void setModelName( std::string const& modelName )
  {
    deviceData_.getDevBlock().setModelName(modelName);
    deviceData_.getDevBlock().modelFlag = (modelName != "");
  }

  void addInstanceParameter(const Device::Param & parameter )
  {
    deviceData_.getDevBlock().params.push_back( parameter );
  }

  void addInstanceParameters(const std::vector<Device::Param> &parameters )
  {
    deviceData_.getDevBlock().params.insert(deviceData_.getDevBlock().params.end(), parameters.begin(), parameters.end());
  }

  void setLineNumber( std::string & netlistFile, int lineNumber )
  {
    deviceData_.getDevBlock().netlistFileName_ = netlistFile;
    deviceData_.getDevBlock().lineNumber_ = lineNumber;
  }

  const std::vector<SpiceSeparatedFieldTool::StringToken> &getParsedLine() const 
  {
    return parsedLine_;
  }

  const std::string &getName() const 
  {
    return deviceData_.getDevBlock().getName();
  }

  const std::string &getModelName() const 
  {
    return deviceData_.getDevBlock().getModelName();
  }

  const std::string getNetlistDeviceType() const 
  {
    return netlistType_;
  }

  int getNumberOfNodes() const 
  {
    return deviceData_.getNodeBlock().get_NodeList().size();
  }

  const std::list<tagged_param> & getNodeValues() const 
  {
    return deviceData_.getNodeBlock().get_NodeList();
  }


  int getNumberOfInstanceParameters() const 
  {
    return deviceData_.getDevBlock().params.size();
  }

  const Device::Param &getInstanceParameter(int i) const 
  {
    return deviceData_.getDevBlock().params[i];
  }

  bool isSubcircuitInstance() const 
  {
    return subcircuitInstance_;
  }

  bool isExtracted() const 
  {
    return extracted_;
  }

  const Topo::NodeDevBlock &getDeviceData() const 
  {
    return deviceData_;
  }

  void setInstanceParameter(int i, Device::Param &parameter );

  Device::Param* findInstanceParameter( Util::Param const& parameter );
  Device::Param* findInstanceParameter( std::string const& parameter );
  void getInstanceParameters(std::vector<Device::Param>& parameters);

  const std::string &getNodeValue(int i) const;

  void getAllNodeNames( std::list< std::string > & nodeNames );

  void setNodeValues( std::list<tagged_param> const& nodeValues );
  void setNodeValue( int const& i, std::string const& nodeValue );


  // Print the details of a device to standard out.
  void print();

  // Clear the device, reset all attributes to their defaults.
  void clear();

  bool extractData();

  // Extract the subcircuit instance data given on a netlist 'X' line.
  bool extractSubcircuitInstanceData();

private:
  void checkNode (const std::string &n);

  bool extractSourceData();

  bool extractMutualInductanceData();

  bool extractBasicDeviceData();

  bool extractBehavioralDeviceData();

  bool extractYDeviceData();

  bool extractUDeviceData();

  bool extractMIDeviceData();

  bool extractSwitchDeviceData();

  bool extractNodes(int modelLevel, int modelNamePosition);

  bool extractModelName( std::string& modelType,
                         int & modelLevel,
                         int & modelNamePosition );

  void extractInstanceParameters( int searchStartPosition,
                                  int & parameterStartPosition,
                                  int & parameterEndPosition,
                                  std::string const& action,
                                  int modelLevel = -1 );

  void addDefaultInstanceParameters(int modelLevel);


  int findSourceFieldPosition( std::string const& fieldToFind,
                               int startPosition );

  bool extractSourceFields( std::vector<std::string> const& fieldNames,
                            std::vector<int> const& fieldPositions );

  // Look for expression valued parameters in the parameter
  // list, evaluate expressions found and reset the parameter
  // value accordingly.
  bool setParameterValues();

  void issueUnrecognizedParameterError(std::string const& parameterName);

  bool isValidDeviceType(const std::string & deviceType);

private:
  CircuitContext &                                    circuitContext_;
  CircuitMetadata &                                   metadata_;

  std::string                                         netlistFileName_;
  std::vector<SpiceSeparatedFieldTool::StringToken>   parsedLine_;

  std::string                                         netlistType_;
  Topo::NodeDevBlock                                  deviceData_;

  bool                                                subcircuitInstance_;
  bool                                                extracted_;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::DeviceBlock N_IO_DeviceBlock;

#endif
