//----------------------------------------------------------------------------
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
// Filename       : $RCSfile: N_IO_DeviceBlock.C,v $
//
// Purpose        : Define the N_IO_DeviceBlock class instantiations of which
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
// Revision Number: $Revision: 1.179.2.1 $
//
// Revision Date  : $Date: 2014/03/10 16:15:21 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <sstream>
#include <iostream>

#include <N_UTL_Misc.h>

#ifdef HAVE_CASSERT
 #include <cassert>
#else
 #include <assert.h>
#endif

#include <N_IO_CircuitContext.h>
#include <N_IO_CircuitMetadata.h>
#include <N_IO_DeviceBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_Report.h>

#include <N_ERH_ErrorMgr.h>
#include <N_TOP_NodeDevBlock.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Expression.h>

#include <N_DEV_SourceData.h>

namespace Xyce {
namespace IO {

//--------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 09/09/2001
//--------------------------------------------------------------------------
DeviceBlock::DeviceBlock(
  CircuitContext & cc,
  CircuitMetadata & md )
  : circuitContext_(cc),
    metadata_(md)
{
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::DeviceBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
DeviceBlock::DeviceBlock(
    CircuitContext & cc,
    CircuitMetadata & md,
    std::string const& fileName,
    std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedInputLine)
: netlistFileName_(fileName),
  parsedLine_(parsedInputLine),
  subcircuitInstance_(false),
  circuitContext_(cc),
  metadata_(md),
  extracted_(false)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::DeviceBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
DeviceBlock::DeviceBlock(DeviceBlock const& rhsDB)
  : netlistFileName_(rhsDB.netlistFileName_),
    parsedLine_(rhsDB.parsedLine_),
    netlistType_(rhsDB.netlistType_),
    deviceData_(rhsDB.deviceData_),
    subcircuitInstance_(rhsDB.subcircuitInstance_),
    extracted_(rhsDB.extracted_),
    circuitContext_(rhsDB.circuitContext_),
    metadata_(rhsDB.metadata_)
{}


//-----------------------------------------------------------------------------
// Function      : DeviceBlock::print
// Purpose       : Output the details of a device block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
void DeviceBlock::print()
{
  Xyce::dout() << std::endl << Xyce::section_divider << std::endl
               << "Device Information" << std::endl
               << "------------------" << std::endl
               << "device line:" << std::endl;
  int numFields = parsedLine_.size();
  int i;
  for ( i = 0; i < numFields; ++i )
  {
    Xyce::dout() << "  " << parsedLine_[i].string_;
  }
  Xyce::dout() << std::endl;
  Xyce::dout() << "  name: " << getName() << std::endl;

  Xyce::dout() << "  nodes: ";
  int numNodes = getNumberOfNodes();
  for ( i = 0; i < numNodes; ++i )
  {
    Xyce::dout() << "Node " << i+1 << ": ";
    Xyce::dout() << getNodeValue(i) << "  ";
  }
  Xyce::dout() << std::endl;

  if ( getModelName() != "" )
  {
    Xyce::dout() << "  model name: " << getModelName() << std::endl;
  }
  Xyce::dout() << std::endl;

  if ( getNumberOfInstanceParameters() > 0 )
  {
    Xyce::dout() << "  Instance Parameters:" << std::endl;
    size_t numParams = getNumberOfInstanceParameters();
    for ( size_t k = 0; k < numParams; ++k )
    {
      Xyce::dout() << "    " << getInstanceParameter(k).uTag();
      Xyce::dout() << "    " << getInstanceParameter(k).stringValue();

      if ( getInstanceParameter(k).given() )
      {
         Xyce::dout() << "    given";
      }
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;
  }

  Xyce::dout() << std::endl << Xyce::section_divider << std::endl;
  Xyce::dout() << std::endl;
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::clear
// Purpose        : Reset all attributes to their defaults.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/06/2003
//----------------------------------------------------------------------------
void DeviceBlock::clear()
{
  netlistFileName_ = "";
  parsedLine_.clear();
  netlistType_ = "";
  subcircuitInstance_ = false;

  setModelName("");
  deviceData_.clear();
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractData
// Purpose       : Extract the device data from parsed line. Use device
//                 metadata to determine device type, number of nodes, etc.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractData()
{
  // Set the device name and netlist type.
  ExtendedString deviceName(parsedLine_[0].string_);
  deviceName.toUpper();
  setName(deviceName);
  setNetlistType(deviceName[0]);
  setLineNumber(netlistFileName_, parsedLine_[0].lineNumber_);

  if (!isValidDeviceType(getNetlistDeviceType()))
  {
    Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "Invalid device type for device " << deviceName;
    return false;
  }

  bool result = false;

  // Invoke appropriate method to handle the rest of the line. Independent
  // sources and mutual inductances need to be handled differently than
  // other devices.
  if ( getNetlistDeviceType() == "V" || getNetlistDeviceType() == "I" )
  {
    result = extractSourceData();
  }
  else if ( getNetlistDeviceType() == "E" || getNetlistDeviceType() == "G" ||
            getNetlistDeviceType() == "F" || getNetlistDeviceType() == "H" )

  {
    result = extractBehavioralDeviceData();
  }
  else if ( getNetlistDeviceType() == "K" )
  {
    result = extractMutualInductanceData();
  }
  else if ( getNetlistDeviceType() == "S" || getNetlistDeviceType() == "W" )
  {
    result = extractSwitchDeviceData();
  }
  else if ( getNetlistDeviceType() == "X" )
  {
    result = extractSubcircuitInstanceData();
  }
  else if ( getNetlistDeviceType() == "Y" )
  {
    result = extractYDeviceData();
  }
  else if ( getNetlistDeviceType() == "U" )
  {
    result = extractUDeviceData();
  } 
  else
  {
    result = extractBasicDeviceData();
  }

  // Check the status and extracting device data
  // and return if something went wrong.
  if ( !result )
  {
    Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "Incomprehensible syntax for device " << deviceName;
    return false;
  }

  //allows testing later to avoid "reextracting" data
  extracted_ = true;

  // Now that the data has been extracted, given the circuit context,
  // the parameter values can be set.
  setParameterValues();

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractSubcircuitInstanceData
// Purpose       : Extract the subcircuit instance data given on a netlist 'X'
//                 line.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/28/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractSubcircuitInstanceData()
{
  // Set the flag indicating this device is a subcircuit instance.
  subcircuitInstance_ = true;

  // Set the device name and netlist type.
  ExtendedString fieldES( parsedLine_[0].string_ );
  fieldES.toUpper();
  setName( fieldES );
  setNetlistType( fieldES[0] );

  // Set the number of fields on the line.
  size_t numFields = parsedLine_.size();

  // Look for "PARAMS:" on the line, set the parameter start position and
  // subcircuit name position accordingly.
  int parameterPosition = 0;
  int subcircuitNamePosition = 0;
  ExtendedString field("");
  for ( size_t i = 1; i < numFields; ++i )
  {
    field = parsedLine_[i].string_;
    field.toUpper();
    if ( field == "PARAMS:" )
    {
      parameterPosition = i+1;
      subcircuitNamePosition = i-1;
      break;
    }
    else if ( i < numFields-1 && parsedLine_[i+1].string_ == "=" )
    {
      parameterPosition = i;
      subcircuitNamePosition = i-1;
      break;
    }
  }

  // If the "PARAMS:" keyword was not found, then set the
  // subcircuitNamePos to the last field on the line.
  if ( parameterPosition == 0 )
  {
    subcircuitNamePosition = numFields - 1;
  }

  // Use the DeviceBlock getModelName() field to hold the subcircuit name field.
  fieldES = parsedLine_[ subcircuitNamePosition ].string_;
  fieldES.toUpper();
  setModelName( fieldES );

  // Set list of subcircuit external nodes, the last node
  // on the line precedes the subcircuit name.
  for ( int i = 1; i < subcircuitNamePosition; ++i )
  {
    addNodeValue( ExtendedString(parsedLine_[i].string_).toUpper() );
  }

  // Collect the parameters on the line if there are any.
  if ( parameterPosition > 0 )
  {
    size_t i = parameterPosition;
    Device::Param parameter("","");
    while ( i < numFields )
    {
      fieldES =  parsedLine_[i].string_;
      fieldES.toUpper();
      if (!fieldES.possibleParam())
      {
        Report::UserError0().at(netlistFileName_, parsedLine_[i].lineNumber_)
          << "Parameter name " << parameter.uTag() << " contains illegal character(s)";
      }
      parameter.setTag( fieldES );

      if ( parsedLine_[i+1].string_ != "=" )
      {
        Report::UserError().at(netlistFileName_, parsedLine_[i].lineNumber_)
          << "Equal sign required between parameter and value in subcircuit parameter list for subcircuit instance " << getName();
      }

      i+=2; // Advance past "=" sign

      fieldES =  parsedLine_[i].string_;
      fieldES.toUpper();
      if (fieldES.possibleParam());
        fieldES = "{" + parsedLine_[i].string_ + "}";
      parameter.setVal(  fieldES );

      addInstanceParameter( parameter );

      // Advance to next field.
      ++i;
    }
  }

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractBasicDeviceData
// Purpose       : Extract the device data from parsedLine for devices other
//                 than independent sources, mutual inductances and "Y" devices.
//                 These devices are exected to have input lines of the
//                 following format:
//
//    dname node node [node]* [Value] [model-name] [param=Value]* [Area]
//    +     [ON|OFF] [IC = value [,value]*] [TEMP=Value]
//
//                  where fields enclosed in [] indicates are optional, *
//                  indicates a field that may be reapeated 0 or more times,
//                  fields given in all lower case a string is expected, fields
//                  starting with a capital letter indicate a numerical value
//                  is expected, and fields in all upper case are expected
//                  literally.
//
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/24/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractBasicDeviceData()
{
  bool result;
  int i, n_start, n_end, n_req, n_opt, n_fill;

  // The device name has been extracted by extractData. The
  // remaining data in parsedLine is extracted here.

  int numFields = parsedLine_.size();

  // Extract the model name from parsedLine if one exists. If
  // a model name was found, find its type.
  int modelLevel, modelNamePosition;
  std::string modelType;
  bool modelFound = extractModelName( modelType, modelLevel, modelNamePosition );
  if (!modelFound)
    modelLevel = -1;

  // Some devices require a model, check that a model was found if it is
  // required for this device.

  if (metadata_.isModelRequired(getNetlistDeviceType(), modelLevel) && !modelFound)
  {
    Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "Did not find required model for device " << getName();
    return false;
  }

  std::string netlistModelType("");
  if ( modelFound && modelNamePosition < numFields)
  {
    // Query the metadata for the validity of the device model.
    if (!metadata_.isModelTypeValid(getNetlistDeviceType(), modelType, modelLevel))
    {
      Report::UserError().at(netlistFileName_, parsedLine_[modelNamePosition].lineNumber_)
        << "Model type \"" << modelType << "\" not valid for device "
        << getName() << " of type " << getNetlistDeviceType();
      return false;
    }

    netlistModelType = modelType;
  }

  if ( netlistModelType == "" )
  {
    netlistModelType = "default";
  }

  // Extract the device nodes from parsedLine_.
  result = extractNodes(modelLevel, modelNamePosition);
  if (!result)
    return result;

  // Process optional nodes
  n_req = metadata_.getNumberOfNodes(getNetlistDeviceType(), modelLevel);
  int nodesFound = n_req;
  n_opt = metadata_.getNumberOfOptionalNodes(getNetlistDeviceType(), modelLevel);
  std::string primaryDeviceParameter(
    metadata_.getPrimaryParameter(getNetlistDeviceType(), modelLevel));

  if (n_opt > 0)
  {
    n_start = n_req+1;
    n_end = n_start;
    n_fill = metadata_.getNumberOfFillNodes(getNetlistDeviceType(), modelLevel);
    if (modelFound && modelNamePosition > n_req)
    {
      n_end = modelNamePosition;
    }
    else if (!modelFound && primaryDeviceParameter == "")
    {
      for (i=n_start ; i<numFields ; ++i)
      {
        if (parsedLine_[i].string_ == "=")
        {
          n_end = i - 1;
          break;
        }
      }
    }
    if (n_end > n_start)
    {
      for (i=n_start ; i<n_end ; ++i)
      {
        addNodeValue( ExtendedString(parsedLine_[i].string_).toUpper());
        ++nodesFound;
      }
    }
    for (i=n_end ; i<=n_req+n_fill ; ++i) {
      addNodeValue( "0" );
    }
  }
  deviceData_.getDevBlock().numExtVars = nodesFound;

  // Check for value field and set position of start of instance
  // parameters on the line.
  Device::Param primaryParameter( "", "" );
  int parameterStartPosition;
  if ( primaryDeviceParameter != "" )
  {
    // test for un-tagged parameter after model name and use as offset
    int offsetUntaggedParamAfterModel = 0;
    if ( modelFound && modelNamePosition + 1 < numFields &&
     !metadata_.isDeviceParameter( getNetlistDeviceType(), modelLevel,
     parsedLine_[ modelNamePosition + 1 ].string_ ) )
    {
      offsetUntaggedParamAfterModel = 1;
    }

    if ( ((getNumberOfNodes() + 1 < numFields) && !modelFound) ||
         (modelFound && (getNumberOfNodes()+1 < modelNamePosition)) ||
         (modelFound && (getNumberOfNodes()+1 == modelNamePosition) &&
         bool(offsetUntaggedParamAfterModel)) )
    {
      if (metadata_.isDeviceParameter( getNetlistDeviceType(), modelLevel,
            parsedLine_[getNumberOfNodes()+1+offsetUntaggedParamAfterModel].string_) &&
          !(numFields < getNumberOfNodes()+3+offsetUntaggedParamAfterModel ||
           parsedLine_[getNumberOfNodes()+2+offsetUntaggedParamAfterModel].string_ != "="))
      {
        // If the field following the node list is in the instance
        // parameter list then it is not the value field.
        parameterStartPosition = getNumberOfNodes() + 1+offsetUntaggedParamAfterModel;
      }
      else
      {
        // The field following the node list is the value field.
        primaryParameter.set(
                       primaryDeviceParameter,
                       parsedLine_[getNumberOfNodes()+1+offsetUntaggedParamAfterModel].string_ );

        // The value field must either be a valid number or an expression.
        if (primaryParameter.getType() == Xyce::Util::STR && !primaryParameter.isNumeric())
        {
          ExtendedString p_orig(primaryParameter.stringValue());
          p_orig.toUpper();
          if (p_orig.possibleParam())
          {
            primaryParameter.setVal(std::string("{" + p_orig + "}"));
            if (!circuitContext_.resolveParameter(primaryParameter))
              primaryParameter.setVal(p_orig);
          }
        }
        if ( !primaryParameter.isNumeric() &&
             !primaryParameter.hasExpressionValue() )
        {
          Report::UserError().at(netlistFileName_, parsedLine_[getNumberOfNodes()+1+offsetUntaggedParamAfterModel].lineNumber_)
            << "Illegal value found for device " <<  getName();
        }

        if ( modelFound )
        {
          parameterStartPosition = getNumberOfNodes() + 3;
        }
        else
        {
          parameterStartPosition = getNumberOfNodes() + 2;
        }
      }
    }
    else if ( modelFound && (getNumberOfNodes()+1 == modelNamePosition) )
    {
      // No device type parameter on line.
      parameterStartPosition = modelNamePosition + 1;
    }
    else if ( getNumberOfNodes()+1 >= numFields )
    {
      Report::UserWarning0().at(netlistFileName_, parsedLine_[getNumberOfNodes()].lineNumber_)
        << "Expected value field for device " << getName() << ", continuing with value of 0";

      // There are no instance parameters given for this device,
      // set this variable for consistency.
      parameterStartPosition = getNumberOfNodes() + 1;
    }
  }
  else
  {
    // The device does not have a value field.
    if ( modelFound )
    {
      parameterStartPosition = modelNamePosition + 1;
    }
    else
    {
      parameterStartPosition = getNumberOfNodes() + 1;
    }
  }

  // Add the device instance parameters and their default values
  // to instanceParameters and check parsedLine for instance
  // parameters.
  int parameterEndPosition;
  extractInstanceParameters( parameterStartPosition,
                             parameterStartPosition,
                             parameterEndPosition, "break",
                             modelLevel);

  // If the primary device parameter was found, reset its value.
  if ( primaryParameter.tag() != "" )
  {
    Device::Param *parameterPtr = findInstanceParameter( primaryParameter );
    if (parameterPtr)
    {
      // parameterPtr->setVal( static_cast<Util::Param &>(primaryParameter));
      setParamValue(*parameterPtr, primaryParameter);
      parameterPtr->setGiven( true );
    }
    else
    {
      if ( getNetlistDeviceType() == "L")
      {
        addInstanceParameter( Device::Param (primaryParameter.tag(), parsedLine_[modelNamePosition-1].string_));
      }
    }
  }

  // Check for AREA, OFF, IC and TEMP parameters.

  // Check for Area parameter.
  int linePosition = parameterEndPosition + 1;
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedLine_[linePosition].string_ );
    field.toUpper();
    if ( field != "ON" && field != "OFF" && field != "IC" && field != "TEMP"  )
    {
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("AREA", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( parsedLine_[linePosition].string_ );
        parameterPtr->setGiven( true );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
         // Either the device does not have an area parameter, or there is
         // some unrecognized parameter on the line. Either way, flag an error
         // on the line and stop.
         issueUnrecognizedParameterError(field);
      }
    }
  }

  // Check for ON parameter.
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedLine_[linePosition].string_ );
    field.toUpper();
    if ( field == "ON" )
    {
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("ON", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( 1.0 );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
        // There is no ON parameter for this device.
        issueUnrecognizedParameterError("ON");
      }
    }
  }

  // Check for OFF parameter.
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedLine_[linePosition].string_ );
    field.toUpper();
    if ( field == "OFF" )
    {
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("OFF", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( 1.0 );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
        // There is no OFF parameter for this device.
        issueUnrecognizedParameterError("OFF");
      }
    }
  }


  // Check for temperature parameter.
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedLine_[linePosition].string_ );
    field.toUpper();
    if ( field == "TEMP" )
    {
      linePosition += 2; // advance past "=" sign
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("TEMP", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( parsedLine_[linePosition].string_ );
        parameterPtr->setGiven( true );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
         // This device does not have a TEMP parameter.
         issueUnrecognizedParameterError("TEMP");
      }
    }
  }

  // Issue fatal error if there are more fields on the line.
  if ( linePosition < numFields )
  {
    Report::UserError0().at(netlistFileName_, parsedLine_[linePosition+1].lineNumber_)
      << "Unrecognized fields for device " << getName();
  }

  return true; // Only get here on success.
}


//----------------------------------------------------------------------------
// Function       : DeviceBlock::extractBehavioralDeviceData
// Purpose        : Extract the device data for a "E", "F", "G", and "H" devices. If
//                  the usage of these devices is as a PSpice style behavioral
//                  device (i.e. a dependent source device), convert the device
//                  to a Xyce B-source device.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 12/09/2002
//----------------------------------------------------------------------------
bool DeviceBlock::extractBehavioralDeviceData()
{
  bool result;

  int numFields = parsedLine_.size();

  // Check the device line for a "VALUE", "TABLE", or "POLY" field. If
  // none of these fields are found and the device is an E or G device,
  // it can be parsed by the usual device rules (extractBasicDeviceData).
  // If there is a "VALUE", "TABLE", or "POLY" field it should be the
  // fourth field on the line.
  ExtendedString typeField("");
  if ( numFields >= 4 )
  {
    typeField = parsedLine_[3].string_;
    typeField.toUpper();
  }

  if ( (getNetlistDeviceType() == "E" || getNetlistDeviceType() == "G") &&
       (typeField != "VALUE" && typeField != "TABLE" && typeField != "POLY") )
  {
    result = extractBasicDeviceData();
    return result;
  }

  // For all other cases, convert the device to a B-source.

  // Get the device type.
  std::string deviceType(getNetlistDeviceType());

  // Change the type of the device to a B-source.
  setNetlistType ( 'B' );
  deviceData_.getDevBlock().bsourceFlag = true;

  // Extract the device nodes from parsedLine_.
  result = extractNodes(-1, 0);
  if (!result)
    return result;

  std::string expression("");
  if ((deviceType == "E" || deviceType == "G") && typeField != "POLY")
  {
    // Get the expression from the line. The expression may need to be
    // transformed by adding braces where Xyce expects them.
    if ( typeField == "VALUE" )
    {
      int exprStart = 4;
      if ( parsedLine_[exprStart].string_ == "=" ) ++exprStart;

      for ( int i = exprStart; i < numFields; ++i )
      {
        expression += parsedLine_[i].string_;
      }
    }
    else if ( typeField == "TABLE" )
    {
      expression = "TABLE";

      // Find the next "=" sign.
      int equalPos;
      for ( equalPos = 5; equalPos < numFields; ++equalPos )
      {
        if ( parsedLine_[equalPos].string_ == "=" ) break;
      }

      if ( equalPos == numFields )
      {
        // Problem, did not find expected "=" sign.
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Required \"=\" sign missing after the expression in the TABLE for device " << getName();
        return false;
      }

      // Get the table expression.
      std::string tableExpression("");
      for ( int i = 4; i < equalPos; ++i )
      {
        tableExpression += parsedLine_[i].string_;
      }

      // Add braces to tableExpression if needed.
      if ( tableExpression[0] != '{' &&
        tableExpression[tableExpression.size()-1] != '}' )
      {
        tableExpression = "{" + tableExpression + "}";
      }

      expression += " " + tableExpression + " = ";

      // Add the table data to the expression.
      for ( int i = equalPos+1; i < numFields; ++i )
      {
        expression += " " + parsedLine_[i].string_;
      }
    }
  }
  else if (typeField == "POLY")
  {
    // Convert an E, F, G, or H device with a POLY function to the equivalent
    // B-source device.

    // Get the dimension of the polynomial.
    if (parsedLine_.size() < 6)
    {
      Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
        << "Not enough fields on input line for device " << getName();
      return false;
    }

    Device::Param dimensionParam("dim", parsedLine_[5].string_);
    int dimension = dimensionParam.getImmutableValue<int>();

    // Build up the B-source POLY expression, some fields need to
    // be manipulated depending on the device type.
    expression = "POLY(";
    expression += dimensionParam.stringValue() + ")";

    int linePosition = 7;
    if (deviceType == "E" || deviceType == "G")
    {
      if (parsedLine_.size() < 7+2*dimension)
      {
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Not enough fields on input line for device " << getName();
        return false;
      }

      // Get the POLY controlling nodes, there should be as many pairs
      // of these as the dimension of the polynomial.
      for (int i = 0; i < dimension; ++i, linePosition += 2)
      {
        std::string cnode1(parsedLine_[linePosition].string_);
        std::string cnode2(parsedLine_[linePosition+1].string_);

        if (cnode2 == "0")
        {
          expression += " V(" + cnode1 + ")";
        }
        else
        {
          expression += " V(" + cnode1 + "," + cnode2 + ")";
        }
      }
    }
    else
    {
      // This is an F or an H device.
      if (parsedLine_.size() < 7+dimension)
      {
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Not enough fields on input line for device " << getName();
        return false;
      }

      // Get the controlling sources, there should be as many of these as
      // as the dimenstion of the polynomial.
      for (int i = 0; i < dimension; ++i, ++linePosition)
      {
        std::string csource(parsedLine_[linePosition].string_);

        expression += " I(" + csource + ")";
      }
    }

    // Add all remaining fields to the expression.
    for ( ;linePosition < parsedLine_.size(); ++linePosition)
    {
      expression += " " + parsedLine_[linePosition].string_;
    }
  }
  else if (deviceType == "F" || deviceType == "H")
  {
    // Linear F or H conversion to B-source.
    if (parsedLine_.size() < 5)
    {
      Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
        << "Not enough fields on input line for device " << getName();
      return false;
    }

    // Get the source and value.
    std::string source(parsedLine_[3].string_);
    std::string value(parsedLine_[4].string_);

    expression = value + " * I(" + source + ")";
  }
  else
  {
    // Have something unrecognizeable on the line.
  }

  // Added braces to expression if needed.
  if ( expression[0] != '{' &&
      expression[expression.size()-1] != '}' )
  {
    expression = "{" + expression + "}";
  }

  // Set the device parameter
  Device::Param parameter( "", "" );
  if ( deviceType == "E" || deviceType == "H" )
  {
    parameter.setTag( "V" );
    parameter.setVal( expression );
    parameter.setGiven( true );
    addInstanceParameter( parameter );

    parameter.setTag( "I" ); // This B-source parameter is required but
                             // won't be used.
    parameter.setVal( "" );
    parameter.setGiven( false );
    addInstanceParameter( parameter );
  }
  else if ( deviceType == "F" || deviceType == "G" )
  {
    parameter.setTag( "I" );
    parameter.setVal( expression );
    parameter.setGiven( true );
    addInstanceParameter( parameter );

    parameter.setTag( "V" ); // This B-source parameter is required but
                             // won't be used.
    parameter.setVal( "" );
    parameter.setGiven( false );
    addInstanceParameter( parameter );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractYDeviceData
// Purpose       : Extract the device data from parsedLine for "Y" devices.
//                 These devices are exected to have input lines of the following
//                 format:
//
//                 Ytype name node [node]* [Value] [model-name] [param=Value]*
//
//                 where fields enclosed in [] indicates are optional, *
//                 indicates a field that may be reapeated 0 or more times,
//                 fields given in all lower case a string is expected, fields
//                 starting with a capital letter indicate a numerical value
//                 is expected, and fields in all upper case are expected
//                 literally.
//
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/24/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractYDeviceData()
{
  bool result;

  int numFields = parsedLine_.size();

  if ( numFields < 2 )
  {
    Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "Not enough fields on input line for device " << getName();
    return false;
  }

  // Reset the device name and type which were set by extractData based on
  // the standard rules for netlist devices.
  ExtendedString deviceName(parsedLine_[1].string_);
  deviceName.toUpper();
  ExtendedString deviceType ( parsedLine_[0].string_.substr(1,parsedLine_[0].string_.size()) );
  deviceType.toUpper();
  setNetlistType( deviceType );
  setName( "Y%" + deviceType + "%" + deviceName );

  // drop the first element of the parsed line vector
  // and treat this like a basic device
  parsedLine_.erase(parsedLine_.begin());

  // special handling for new mutual inductor
  if( deviceType == "MIL" || deviceType == "MIN" )
  {
    return extractMIDeviceData();
  }
  else
  {
    return extractBasicDeviceData();
  }

}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractUDeviceData
// Purpose       : Extract the device data from parsedLine for "U" devices.
//                 These devices are exected to have input lines of the following
//                 format:
//
//                 Uname type(N) node [node]* [Value] [model-name] [param=Value]*
//
//                 where fields enclosed in [] indicates are optional, *
//                 indicates a field that may be reapeated 0 or more times,
//                 fields given in all lower case a string is expected, fields
//                 starting with a capital letter indicate a numerical value
//                 is expected, and fields in all upper case are expected
//                 literally.  Some digital devices (e.g., NAND, AND, OR and
//                 NOR) may have a variable number of inputs.  In that case,
//                 the number of inputs is denoted by (for example) AND(N),
//                 where N is an integer.  Other gates (e.g, INV) have a fixed
//                 number of inputs, and their type is (for example) INV.
//
// Special Notes : The Xyce U device is intended to be similar to the PSpice U 
//                 digital devices, rather than the SPICE3F5 U device which is
//                 the "uniform RC transmission line".  
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 02/14/2014
//-----------------------------------------------------------------------------
bool DeviceBlock::extractUDeviceData()
{
  bool noError = true;

  int numFields = parsedLine_.size();

  if ( numFields < 2 )
  {
    Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "Not enough fields on input line for U device " << getName();
    noError = false;
  }

  // Parse AND(N) syntax where AND is the gate type, and N is the number
  // of inputs.  Error handling for N not being an integer, N missing or
  // trailing ')' missing.
  std::string numInputs = "";
  if ( numFields >=5)
  {
    if ( parsedLine_[2].string_ == "(" )
    {
      if ( parsedLine_[4].string_ == ")" )
      {
        // found possible (N) syntax
        numInputs = parsedLine_[3].string_;
        if (!Util::isInt(numInputs))
	{
          Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
            << "Found (N) syntax in U device " << getName() << " but N not an integer";
          noError = false;
        }
        else
	{
          // found valid (N) syntax.  Remove (N) from parsed line so that
          // extractBasicDeviceData() function works correctly
          parsedLine_.erase(parsedLine_.begin()+2,parsedLine_.begin()+5);
	}
      }
      else
      {
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Apparent error parsing number of inputs, (N) syntax, for U device " << getName();
        noError = false;
      }
    } 
  }

  // Reset the device name and type which were set by extractData based on
  // the standard rules for netlist devices.  Note that the order of the
  // deviceType and deviceName are reversed from the Y-device syntax
  ExtendedString deviceType(parsedLine_[1].string_);
  deviceType.toUpper();
  ExtendedString deviceName( parsedLine_[0].string_.substr(1,parsedLine_[0].string_.size()) );
  deviceName.toUpper();
  // U devices are not allowed to have a '%' character in them because it hoses
  // subsequent processing in the digital-device constructor.
  if ( (deviceType.find_first_of('%') != std::string::npos) || 
       (deviceName.find_first_of('%') != std::string::npos) )
  {
    Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "U device name: " << deviceName << " or type: " << deviceType <<  " should not have % character.";
    noError = false;  
  }
  setNetlistType( deviceType );
  if (numInputs == "")
  {
    // name string format for U devices with a fixed number of inputs
    setName( "U%" + deviceType + "%" + deviceName );
  }
  else
  {
    // name string format for U devices with a variable number of inputs
    setName( "U%" + deviceType + "%" + deviceName + "%" + numInputs);
  }

  // drop the second element of the parsed line vector, which is the
  // gate type, and treat this like a basic device.  Previous statement 
  // already had warning message: "not enough fields on input line for device"
  if ( numFields >= 2 )
  {
    parsedLine_.erase(parsedLine_.begin()+1);
  }

  // remove the U from the device name, before extracting the device data.
  // this mimics the Y device parsing
  parsedLine_[0].string_ = parsedLine_[0].string_.substr(1,parsedLine_[0].string_.size()); 

  return (noError == true) ? extractBasicDeviceData() : false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractSourceData
// Purpose       : Extract the device data from parsedLine for independent
//                 source devices.
//
// Special Notes : Independent sources have the general form:
//
//   VXXXXXXX N+ N- <DC<> DC/TRAN VALUE> <AC <ACMAG <ACPHASE>>>
// + <DISTOF1 <F1MAG <F1PHASE>>> <DISTOF2 <F2MAG <F2PHASE>>>
//
// N+ and N- are the positive and negative nodes, respectively.
//
// DC/TRAN is the dc and transient analysis value of the source.
// The DC keyword is optional.
//
// ACMAG is the AC magnitude and ACPHASE is the AC phase. AC keyword and
// values only present if this is an AC source.
//
// DISTOF1 and DISTOF2 are the keywords that specify that the independent
// source has distortion inputs at the frequencies F1 and F2 respectively.
// These are not supported yet.
//
// Any independent source can be assigned a time-dependent value for
// transient analysis. If a source is assigned a time-dependent value,
// the time-zero value is used for dc analysis. There are five independent
// source functions: PULSE, EXP, SIN, PWL, and single-frequency FM (SFFM).
//
// Examples:
// VCC 10 0 DC 6
// VIN 13 2 0.001 AC 1 SIN(0 1 1MEG)
// ISRC 23 21 AC 0.333 45.0 SFFM(0 1 10K 5 1K)
// VMEAS 12 9
// VCARRIER 1 0 DISTOF1 0.1 -90.0
// VMODULATOR 2 0 DISTOF2 0.01
// IIN1 1 5 AC 1 DISTOF1 DISTOF2 0.001
//
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/25/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractSourceData()
{
  int numFields = parsedLine_.size();
  bool tranGiven=false;
  bool acGiven=false;
  bool dcGiven=true;

  // Extract the source device nodes.
  extractNodes(-1, 0);


  // Set up the "primary" parameter.  (for vsrc's this is V0 by default, and is DC)
  Util::Param primaryParameter( "", "" );
  std::string primaryDeviceParameter(metadata_.getPrimaryParameter(getNetlistDeviceType(), -1));
  primaryParameter.setTag(primaryDeviceParameter);
  primaryParameter.setVal("0");
  int tranSourceType = _DC_DATA;

  // Set position on line to first field after the nodes.
  int linePosition = getNumberOfNodes() + 1;

  // Check for DC field.  If present, then reset DC parameter (default parameter)
  // to the value from the netlist.  The DC value is "special" because
  // it doesn't require a DC keyword in front of it.  Also, it is always the
  // first type specified on the instance line.
  //
  // All other src types require a keyword.  So, the handling is different.
  if ( linePosition < numFields )
  {
    ExtendedString field(parsedLine_[linePosition].string_);
    field.toUpper();
    int sourceID = metadata_.getSourceFunctionID(field);
    // The field is the DC field, advance past "DC" keyword.
    if ( sourceID == _DC_DATA )
    {
      ++linePosition ;

      field = parsedLine_[linePosition].string_;
      field.toUpper();
      sourceID = metadata_.getSourceFunctionID(field);
    }

    // If the field is NOT related to any recognized non-DC
    // source, then proceed.
    if ( sourceID == _NUM_SRC_DATA )  // ie, type not found
    {
      // Get the DC value.
      if ( linePosition < numFields )
      {
        primaryParameter.setVal( parsedLine_[linePosition].string_ );
        ++linePosition;

        // generate error message if DC value is invalid
        if( !( primaryParameter.isNumeric() || primaryParameter.hasExpressionValue() ) )
        {
          Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
            << "Invalid DC value \"" << primaryParameter.stringValue() << "\" for device " << getName();
        }
      }
    }
    else
    {
      // no-op.  primary param was already set to zero.
    }
  }

  // Find the positions of the "AC", "DISTOF1", "DISTOF2", and source
  // function fields in parsedLine_. Save names and positions of the
  // fields that appear in parsedLine in fieldNames and fieldPositions.
  std::vector<std::string> fieldNames;
  std::vector<int> fieldPositions;

  int fieldPosition = findSourceFieldPosition( "AC", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "AC" );
    fieldPositions.push_back( fieldPosition );
    acGiven=true;
//    tranSourceType = metadata_.getSourceFunctionID("AC"); //tmei: 05/02
//    ++linePosition;
  }

  fieldPosition = findSourceFieldPosition( "DISTOF1", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "DISTOF1" );
    fieldPositions.push_back( fieldPosition );
  }

  fieldPosition = findSourceFieldPosition( "DISTOF2", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "DISTOF2" );
    fieldPositions.push_back( fieldPosition );
  }

  fieldPosition = findSourceFieldPosition( "SOURCEFCN", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "SOURCEFCN" );
    fieldPositions.push_back( fieldPosition );
    ExtendedString field( parsedLine_[fieldPosition].string_ );
    field.toUpper();
    tranSourceType = metadata_.getSourceFunctionID(field);
    tranGiven=true;
  }
  else
  {
    tranGiven=false;
  }

  // Set the type(s) of the source.
  Util::Param tranSourceTypeParameter( "TRANSIENTSOURCETYPE", tranSourceType );
  Util::Param acSourceTypeParameter( "ACSOURCETYPE", (int) _AC_DATA);
  Util::Param dcSourceTypeParameter( "DCSOURCETYPE", (int) _DC_DATA);

  // Add the device instance parameters and their default values
  // to instanceParameters and extract instance parameters from
  // parsedLine_.

  int parameterStartPosition;
  int parameterEndPosition;
  extractInstanceParameters( linePosition,
                             parameterStartPosition,
                             parameterEndPosition,
                             "continue");
  fieldNames.push_back( "PARAMS" );
  fieldPositions.push_back( parameterStartPosition );

  // ERK.  Before extractSourceFields can be called, the field vectors
  // need to be in order of ascending field position.
  //
  // This loop below does a quick insertion sort.
  // (see the insort routine from: http://www.yendor.com/programming/sort/)
  //
  // It is quick // as long as the field arrays are less than 15 elements.  I did
  // not use STL sort b/c I am sorting using the contents of fieldPositions,
  // but fieldNames must also be sorted the same way.
  int len = fieldPositions.size();
  for (int i = 1; i < len; i++)
  {
    int j = i;
    int temp = fieldPositions[j];
    std::string tmpName =  fieldNames[j];
    while (j > 0 && (fieldPositions[j-1] > temp))
    {
      fieldPositions[j] = fieldPositions[j-1];
      fieldNames[j] = fieldNames[j-1];
      j--;
    }
    fieldPositions[j] = temp;
    fieldNames[j] = tmpName;
  }

  // Extract the "AC", "DISTOF1", "DISTOF2", and transient
  // source function fields if they appear in parsedLine_.
  bool result = extractSourceFields( fieldNames, fieldPositions );
  if ( !result )
  {
    return false;
  }

  // Set the DC value.
  Device::Param* parameterPtr;
  if ( primaryParameter.tag() != "" )
  {
    parameterPtr = findInstanceParameter( primaryParameter );
    parameterPtr->setVal( primaryParameter );
    parameterPtr->setGiven( true );
  }

  // Set the transient source type.
  parameterPtr = findInstanceParameter( tranSourceTypeParameter );
  parameterPtr->setVal( tranSourceTypeParameter );
  parameterPtr->setGiven( tranGiven );

  // Set the AC source type, if it was given
  parameterPtr = findInstanceParameter( acSourceTypeParameter );
  parameterPtr->setVal( acSourceTypeParameter );
  parameterPtr->setGiven( acGiven );

  // Set the DC source type, if it was given
  parameterPtr = findInstanceParameter( dcSourceTypeParameter );
  parameterPtr->setVal( dcSourceTypeParameter );
  parameterPtr->setGiven( dcGiven );

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::findSourceFieldPosition
// Purpose       : Find the position of a given field of an independent source
//                 in parsedLine_. Return the position if the field is found,
//                 otherwise return 0.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/02/2001
//-----------------------------------------------------------------------------
int DeviceBlock::findSourceFieldPosition( std::string const& fieldToFind,
                                          int startPosition )
{
  size_t numFields = parsedLine_.size();
  int fieldPosition = 0;
  int temp;

  if ( fieldToFind == "SOURCEFCN" ) // any transient source type.
  {
    ExtendedString field("");
    for ( size_t i = startPosition; i < numFields; ++i )
    {
      // The test here checks the field at position i to see if it
      // corresponds to a source function type registered in metadata.
      field = parsedLine_[i].string_;
      field.toUpper();
      if ( metadata_.getSourceFunctionID(field) != _NUM_SRC_DATA && // ie, not found
           metadata_.getSourceFunctionID(field) != _AC_DATA ) // AC
      {
        fieldPosition = i;
        break;
      }
    }
  }
  else
  {
    ExtendedString field("");
    for ( size_t i = startPosition; i < numFields; ++i )
    {
      field = parsedLine_[i].string_;
      field.toUpper();
      if ( field == fieldToFind )
      {
        fieldPosition = i;
        break;
      }
    }
  }

  return fieldPosition;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractSourceFields
//
// Purpose       : This function loops through the fields on a source instance
//                 line, processes them, adding defaults when necessary, and
//                 then adds them to the instance parameter database for
//                 this device.
//
// Special Notes : ERK:  this function does at least one really strange thing.
//                 Instead of plugging values into existing metadata, the
//                 parameters for the given source type are added to the
//                 metadata container.  So, for example, if the user has
//                 set V1, and V1 is already in the metadata, this function
//                 does not set the existing V1.  It adds a new one so that
//                 two instances of V1 are present.
//
//                 In practice, this still works, as the devices (vsrc and
//                 isrc) will use the last value set for any parameter,
//                 and the "added" version is the last value.  However,
//                 for a develoer it is confusing.
//
//                 I think the reason for the "adding" is the PWL source,
//                 which can have an arbitrary number of entries.
//
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/03/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractSourceFields( std::vector<std::string> const& fieldNames,
                                       std::vector<int> const& fieldPositions )
{
  if ( fieldNames[0] == "PARAMS" )
  {
    // No fields to extract.
    return true;
  }

  // Iterate through the fieldNames and extract the corresponding data
  // for the field. The fieldPositions vector gives the position in
  // parsedLine at which the fieldNames appear.
  int i;
  int numFields = fieldNames.size();
  for ( i = 0; i < numFields; ++i )
  {
//    if ( fieldNames[i] == "AC" || fieldNames[i] == "DISTOF1" ||
//         fieldNames[i] == "DISTOF2" )
    if ( fieldNames[i] == "AC")
    {
      // The code here inserts instance parameters for magnitude and phase
      // of AC stimulus.
      Device::Param mag("", "");
      Device::Param phase("", "");

      int numTerms = fieldPositions[i+1] - fieldPositions[i] - 1;

      std::string magName(fieldNames[i] + "MAG");
      std::string phaseName(fieldNames[i] + "PHASE");
      mag.set( magName, "1.0" );
      phase.set( phaseName, "0.0" );

      if ( numTerms > 2 )
      {
        issueUnrecognizedParameterError(fieldNames[i]+" Too Many Terms");
      }

      if ( numTerms >= 1 )
      {
        mag.setVal( parsedLine_[fieldPositions[i]+1].string_ );
        mag.setGiven( true );
      }

      if ( numTerms == 2 )
      {
        phase.setVal( parsedLine_[fieldPositions[i]+2].string_ );
        phase.setGiven( true );
      }

      // Add magnitude and phase parameters to instanceParameters.
      addInstanceParameter( mag );
      addInstanceParameter( phase );

    }
    else if ( fieldNames[i] == "SOURCEFCN" )
    {
      // Determine the source function.
      ExtendedString sourceFunction ( parsedLine_[ fieldPositions[i] ].string_ );
      sourceFunction.toUpper();

      // Check for optional parentheses around source function parameters,
      // and set start and end position of source function parameters on
      // the input line (skip the initial parenthese).
      int sourceFunctionParamStart;
      int sourceFunctionParamEnd;

      // skip over nested parentheses
      int offset=1;
      while( parsedLine_[ fieldPositions[i] + offset ].string_ == "(" )
      {
        ++offset;
      }

      sourceFunctionParamStart = fieldPositions[i] + offset;
      sourceFunctionParamEnd = fieldPositions[i + 1] - offset;

      if ( sourceFunction != "PWL" )
      {
        // Get the source function parameters from metadata, set their
        // values to the given input value and add to the device instance
        // parameter list.

        // Set the number of parameters for this source function on
        // the input line.
        size_t numInputSourceFunctionParams = sourceFunctionParamEnd - sourceFunctionParamStart + 1;

        const std::vector<Device::Param> &sourceFunctionParameters = metadata_.getSourceFunctionParameters(sourceFunction);
        for ( size_t k = 0; k < sourceFunctionParameters.size(); ++k )
        {
          Device::Param parameter("", 0.0);
          parameter.setTag( sourceFunctionParameters[k].uTag() );
          if ( k < numInputSourceFunctionParams )
          {
            parameter.setVal( parsedLine_[sourceFunctionParamStart + k].string_ );
            parameter.setGiven( true );
            addInstanceParameter( parameter );
          }
        }
      }
      else // PWL source
      {
        if (parsedLine_[sourceFunctionParamStart].string_ != "FILE")
        {
          // Xyce expects the "R" parameter to be tagged "REPEATTIME". Find
          // that parameter and rename it.
          Device::Param* parameterPtr =
            findInstanceParameter( Device::Param("R", "") );
          assert(parameterPtr != NULL);
          parameterPtr->setTag( "REPEATTIME" );

          // Locate the start of the (time, value) pairs in the PWL source.
          int timeValPairStart = sourceFunctionParamStart;
          while ( parsedLine_[timeValPairStart].string_ == "R" ||
              parsedLine_[timeValPairStart].string_ == "TD" )
          {
            // If the repeat time appears, reset REPEAT (a boolean indicating
            // the PWL source function is periodic) to 1, aka true.
            if ( parsedLine_[timeValPairStart].string_ == "R" )
            {
              Device::Param* parameterPtr =
                findInstanceParameter( Device::Param("REPEAT", "") );
              assert(parameterPtr != NULL);
              parameterPtr->setVal( true );
              parameterPtr->setGiven( true );
            }

            timeValPairStart += 3;
          }

          // In the PWL source case, sourceFunctionParamEnd could have been
          // miscaluated if the repeat time or time delay appeared with no
          // other tagged parameters (which are required to be at the end of
          // the line. In this case recalculate teh value of
          // sourceFunctionParamEnd.
          if ( sourceFunctionParamEnd < timeValPairStart )
          {
            //sourceFunctionParamEnd = (int) parsedLine_.size() - 1;
            sourceFunctionParamEnd = static_cast<int> (parsedLine_.size()) - 1;
            if ( parsedLine_[sourceFunctionParamEnd].string_ == ")" )
            {
              --sourceFunctionParamEnd;
            }
          }

          // Add (time, value) pairs to the instance parameters. Also,
          // reset the parameter ("NUM") that indicates the number of
          // such pairs and "REPEAT" which is a boolean that indicates
          // that PWL function is periodic.
          Device::Param time( "T", "" );
          Device::Param value( "V", "" );
          int numTimeValuePairs = 0;
          int i = timeValPairStart;
          while ( i < sourceFunctionParamEnd )
          {
            ++numTimeValuePairs;

            time.setVal( parsedLine_[i].string_ );
            time.setGiven( true );
            addInstanceParameter( time );
            ++i;

            value.setVal( parsedLine_[i].string_ );
            value.setGiven( true );
            addInstanceParameter( value );
            ++i;
          }

          if ( numTimeValuePairs > 0 )
          {
            Device::Param* parameterPtr =
              findInstanceParameter( Device::Param("NUM", "") );
            assert(parameterPtr != NULL);
            parameterPtr->setVal( numTimeValuePairs );
            parameterPtr->setGiven( true );
          }
          else
          {
            // no FILE parameter given (that's the next else block)
            // so there should have been some time, value pairs given
            // however, numTimeValueParis == 0.
            Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
              << "Could not parse time/value pairs for piecewise linear function in device: " << getName();
          }
        }
        else
        {
          // PWL time value pairs are given in specified file. Get the
          // file name, open and read (time, value) pairs.
          std::string tvFileNameIN(parsedLine_[sourceFunctionParamStart+1].string_);

          // If the file name is enclosed in double quotes, strip them.
          std::string tvFileName;
          if (tvFileNameIN[0] == '"' &&
              tvFileNameIN[tvFileNameIN.length()-1] =='"')
          {
            tvFileName = tvFileNameIN.substr(1, tvFileNameIN.length()-2);
          }

          std::ifstream tvDataIn;
          tvDataIn.open(tvFileName.c_str(), std::ios::in);
          if ( !tvDataIn.is_open() )
          {
            Report::UserError() << "Could not find file " << tvFileName;
          }

          Device::Param time( "T", "" );
          Device::Param value( "V", "" );

          int numTimeValuePairs = 0;
          double timeIn;
          double valueIn;
          while (tvDataIn >> timeIn)
          {
            char ch;
            tvDataIn.get(ch);

            if (tvDataIn >> valueIn)
            {
              ++numTimeValuePairs;

              time.setVal(timeIn);
              time.setGiven( true );
              addInstanceParameter(time);

              value.setVal(valueIn);
              value.setGiven( true );
              addInstanceParameter(value);
            }
            else
            {
              Report::UserError() << "Problem reading " << tvFileName << std::endl
                                  << "File format must be comma, tab or space separated value. There should be no extra spaces or tabs "
                                  << "around the comma if it is used as the separator.";
            }
          }

          tvDataIn.close();

          if ( numTimeValuePairs > 0 )
          {
            Device::Param* parameterPtr =
              findInstanceParameter( Device::Param("NUM", "") );
            assert(parameterPtr != NULL);
            parameterPtr->setVal( numTimeValuePairs );
            parameterPtr->setGiven( true );
          }

          else
          {
            Report::UserError() << "Failed to successfully read " << tvFileName;
          }
        }
      }
    }
  }

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractMutualInductanceData
// Purpose       : Extract the device data from parsedLine for mutual
//                 inductances.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/25/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractMutualInductanceData( )
{
  int numFields = parsedLine_.size();
  bool kequals = false;

  // Extract the inductor names.
  int numInductors = 0;
  Device::Param parameter( "", "" );
  while ( parsedLine_[numInductors+1].string_[0] == 'L' ||
          parsedLine_[numInductors+1].string_[0] == 'l' ||
          parsedLine_[numInductors+1].string_[0] == 'X' ||
          parsedLine_[numInductors+1].string_[0] == 'x'
          )
  {
    //Have to add some additional checking if the inductor is part of a
    //subcircuit instance.

    if (parsedLine_[numInductors+1].string_[0] == 'X' ||
          parsedLine_[numInductors+1].string_[0] == 'x'
        )
    {
      //For now, we're just not going to allow this.  Fix later with a new
      //parser, hopefullY:
      Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
        << "Subcircuit calls ('X' devices) are not allowed in mutual inductor definitions.";
    }
    else
    {
      ++numInductors;

      parameter.setTag( parsedLine_[numInductors].string_ );
      parameter.setVal( numInductors );
      addInstanceParameter( parameter );
    }
  }

  parameter.setTag( "COUPLING" );

  //Have to check whether coupling is present in the form "k = {coupling}"

  if ((parsedLine_[numInductors + 1].string_ == "K") ||
      ( (parsedLine_[numInductors + 1].string_ == "k")
        && (parsedLine_[numInductors + 2].string_ == "=")
        && (numInductors + 3 < numFields) ) )
  {
    kequals = true;
    parameter.setVal( parsedLine_[numInductors + 3].string_ );
  }
  else
  {
    parameter.setVal( parsedLine_[numInductors + 1].string_ );
  }

  addInstanceParameter( parameter );

  int modelLevel, modelNamePosition;
  std::string modelType;
  bool modelFound = extractModelName( modelType, modelLevel, modelNamePosition );

  // check format for errors; bogus Lnames are handled later
  if ( modelFound )
  {
    if( kequals )
    {
      if ( ( modelNamePosition != ( numInductors + 4 ) ) ||
           ( modelNamePosition != ( numFields - 1 ) ) )
      {
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Malformed line for device " << getName();
      }
    }
    else
    {
      if ( ( modelNamePosition != ( numInductors + 2 ) ) ||
           ( modelNamePosition != ( numFields - 1 ) ) )
      {
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Malformed line for device " << getName();
      }
    }
  }
  else
  {
    if (kequals)
    {
      if ( ( numInductors + 3 ) != ( numFields - 1 ) )
      {
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Specified model not found for device " << getName();
      }
    }
    else
    {
      if ( ( numInductors + 1 ) != ( numFields - 1 ) )
      {
        Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
          << "Specified model not found for device " << getName();
      }
    }
  }

  return true;  // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractSwitchDeviceData
// Purpose       : Convert SPICE switch formats to general expression controlled
//                 switch
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/04
//-----------------------------------------------------------------------------
bool DeviceBlock::extractSwitchDeviceData( )
{
  bool result, modelFound;
  int modelLevel, modelNamePosition, controlPosition;
  std::string modelType, expression;
  int numFields = parsedLine_.size();
  int i, parameterStartPosition, parameterEndPosition;
  bool is_w=false;

  if ( getNetlistDeviceType() == "W" )
  {
    is_w = true;
    setNetlistType ( 'S' );
  }

  modelFound = extractModelName( modelType, modelLevel, modelNamePosition );
  if (!modelFound)
  {
    Report::UserError message;

    message << "No model found in ";
    for (int k=0 ; k<numFields ; ++k)
      message << parsedLine_[k].string_ << " ";
  }
  controlPosition = 0;
  for (i=modelNamePosition+1 ; i<numFields ; ++i)
  {
    if (ExtendedString(parsedLine_[i].string_).toUpper() == "CONTROL")
    {
      controlPosition = i;
      break;
    }
  }

  if ( is_w )
  {
    if (modelNamePosition != 4 || controlPosition != 0)
    {
      Report::UserError msg;
      if (modelNamePosition != 4)
        msg << "Wrong number of nodes in ";
      else
        msg << "Control param in:\n";
      for (int k=0 ; k<numFields ; ++k)
        msg << parsedLine_[k].string_ << " ";
    }
    setNetlistType ( 'S' );
    expression = "{I(" + parsedLine_[modelNamePosition-1].string_ + ")}";
  }
  else
  {
    if ((modelNamePosition != 3 && controlPosition > 0) ||
        (modelNamePosition != 5 && controlPosition == 0))
    {
      Report::UserError msg;
      msg << "Wrong number of nodes in ";
      for (int k=0 ; k<numFields ; ++k)
        msg << parsedLine_[k].string_ << " ";
    }
    if (modelNamePosition == 5)
      expression = "{V(" + parsedLine_[3].string_ + ")-V(" +parsedLine_[4].string_ + ")}";
  }

  addNodeValue( ExtendedString(parsedLine_[1].string_).toUpper());
  addNodeValue( ExtendedString(parsedLine_[2].string_).toUpper());
  parameterStartPosition = modelNamePosition+1;
  parameterEndPosition = numFields;
  extractInstanceParameters( parameterStartPosition,
                             parameterStartPosition,
                             parameterEndPosition,
                             "break");

  if (controlPosition == 0)
  {
    Device::Param parameter;
    int numParameters = getNumberOfInstanceParameters();
    for ( i = 0; i < numParameters; ++i )
    {
      if ( getInstanceParameter(i).uTag() == "CONTROL")
      {
        parameter = getInstanceParameter(i);
        parameter.setVal( expression );
        parameter.setGiven( true );
        setInstanceParameter( i, parameter );
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/27/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractNodes(int modelLevel, int modelNamePosition)
{
  size_t numFields = parsedLine_.size();
  int numNodes;

  numNodes = metadata_.getNumberOfNodes(getNetlistDeviceType(), modelLevel);
  if (numNodes == -1)
    return false;

  const size_t nodeStartPos = 1;
  const size_t nodeEndPos = nodeStartPos + numNodes - 1;

  if ( modelNamePosition > 0 && modelNamePosition <= nodeEndPos)
  {
    if (metadata_.isModelTypeValid(getNetlistDeviceType(), parsedLine_[modelNamePosition].string_, modelLevel))
      Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
        << "Insufficient nodes specified or node name '" << parsedLine_[modelNamePosition].string_ << "' matches one of this device's model name";
    else
      Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "Insufficient number of nodes specified";

    return false;
  }

  if ( numFields < nodeEndPos + 1)
  {
    Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
      << "Not enough fields on input line for device " << getName();
    return false;
  }

  std::list<tagged_param> nodeValues;
  for ( int i = nodeStartPos; i <= nodeEndPos; ++i )
  {
    nodeValues.push_back(
        tagged_param(ExtendedString(parsedLine_[i].string_).toUpper(),0));
  }
  setNodeValues(nodeValues);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractModelName
// Purpose       : Check parsedLine for existance of model name. If a model
//                 name exists, set its position in parsedLine and the device
//                 model name and return true, otherwise return false.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/27/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractModelName(
  std::string & modelType,
  int &         modelLevel,
  int &         modelNamePosition )
{
  modelType = "";
  modelNamePosition = 0;
  modelLevel = 0;

  const std::string &device_type = getNetlistDeviceType();

  const DeviceMetadata &model_group_metadata = metadata_.getDeviceMetadata(device_type, -1);

  size_t model_search_begin_index = model_group_metadata.numNodes + 1;
  // size_t model_search_end_index = std::min(parsedLine_.size(), model_search_begin_index + model_group_metadata.numFillNodes + model_group_metadata.numOptionalNodes + 1);
  // size_t model_search_begin_index = 1;
  size_t model_search_end_index = parsedLine_.size();

  bool modelFound = false;
  for (size_t i = model_search_begin_index; i < model_search_end_index; ++i )
  {
    const std::string &model_name = parsedLine_[i].string_;
    if (device_type == "K") {
      if (model_name == "=" && !equal_nocase(parsedLine_[i - 1].string_, "K"))
        break;
    }
    else {
      if (model_name == "=") // equal sign implies parameter is being parsed, we've gone to far
        break;
    }

    ParameterBlock* modelPtr;
    if (circuitContext_.findModel(model_name, modelPtr))
    {
      setModelName(model_name);
      modelNamePosition = i;

      modelType = modelPtr->getType();
      modelLevel = modelPtr->getLevel();
      if ( device_type != "K")
      {
        const DeviceMetadata &model_metadata = metadata_.getDeviceMetadata(device_type, modelLevel);
        if (!model_metadata.isModelLevelValid() )
        {
          Report::UserError() << "Model type \"" << modelType << "\" does not have level " << modelLevel << " defined";
          return false;
        }
        if (!model_metadata.isModelTypeValid(modelType))
        {
          Report::UserError() << "Model type \"" << modelType << "\" not valid for device " << getName() << " of type " << device_type;
          return false;
        }

        modelPtr->addDefaultModelParameters(metadata_);
      }
      modelFound = true;
      break;
    }
  }

  return modelFound;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/27/2001
//-----------------------------------------------------------------------------
void DeviceBlock::extractInstanceParameters( int searchStartPosition,
                                             int & parameterStartPosition,
                                             int & parameterEndPosition,
                                             std::string const& action,
                                             int modelLevel )
{
  addDefaultInstanceParameters(modelLevel);

  // Check for instance parameters on the device line, reset the parameter
  // value in instanceParameters to the value given.
  int linePosition = searchStartPosition;
  bool foundParameters = false;
  int numFields = parsedLine_.size();
  parameterStartPosition = numFields; // Default in case there are no
                                      // parameters on the line.
  Device::Param paramToFind;
  while (linePosition < numFields)
  {
    // Look for the parameter in the instance parameter list.
    paramToFind.setTag(parsedLine_[linePosition].string_);
    Device::Param* parameterPtr = findInstanceParameter(paramToFind);

    if ( parameterPtr != NULL )
    {

      ExtendedString es(parsedLine_[linePosition].string_);
      es.toUpper();

      if ( !foundParameters )
      {
        foundParameters = true;
        parameterStartPosition = linePosition;
      }

      if ( es == "ON" || es == "OFF")
      {
        // If the ON or OFF field is untagged set its value to 1 (true)
        // and continue to next field.
        if ( !(linePosition + 1 < numFields &&
             parsedLine_[linePosition + 1].string_ == "=") )
        {
          parameterPtr->setVal( 1.0 );
          parameterPtr->setGiven( true );
          ++linePosition;
          continue;
        }
      }

      ++linePosition;    // advance to next field in parsedLine

      if ( linePosition >= numFields )
      {
        // Hit end of line unexpectedly.
        Report::UserError().at(netlistFileName_, parsedLine_[linePosition-1].lineNumber_)
          << "Unexpectedly reached end of line while looking for parameters for device " << getName();
      }

      if ( parsedLine_[linePosition].string_ == "=" )
      {
        ++linePosition; // if field is "=", advance to next field
      }

      if (parameterPtr->stringValue() != "VECTOR" &&
          parameterPtr->stringValue() != "VECTOR-COMPOSITE")
      {
        if ( linePosition >= numFields )
        {
          // Hit end of line unexpectedly.
          Report::UserError().at(netlistFileName_, parsedLine_[linePosition-1].lineNumber_)
            << "Unexpectedly reached end of line while looking for parameters for device " << getName();
        }

        if (DEBUG_IO)
          Xyce::dout() << " Setting parameter " << parameterPtr->uTag()
                       << "to value " << parsedLine_[linePosition].string_ << std::endl;

        // Set the parameter value.
        if (parameterPtr->getType() == Xyce::Util::DBLE)
        {
          std::string & tmpStr = (parsedLine_[linePosition].string_);
          if (Util::possibleParam(tmpStr))
          {
            parameterPtr->setVal( "{" + parsedLine_[linePosition].string_ + "}" );
          }
          else if (Util::isValue(tmpStr))
          {
            parameterPtr->setVal( Util::Value(tmpStr) );
          }
          else
          {
            parameterPtr->setVal( parsedLine_[linePosition].string_ );
          }
        }
        else if (parameterPtr->getType() == Xyce::Util::INT)
        {
          std::string & tmpStr = (parsedLine_[linePosition].string_);

          if (Util::possibleParam(tmpStr))
          {
            parameterPtr->setVal( "{" + parsedLine_[linePosition].string_ + "}");
          }
          else if (Util::isInt(tmpStr))
          {
            parameterPtr->setVal( Util::Ival(tmpStr) );
          }
          else
          {
            parameterPtr->setVal( parsedLine_[linePosition].string_ );
          }
        }
        else
        {
          parameterPtr->setVal( parsedLine_[linePosition].string_ );
        }
        parameterPtr->setGiven( true );

        ++linePosition; // Advance to next field.
      }
      else if (parameterPtr->stringValue() == "VECTOR")
      {
        std::ostringstream paramName;
        Device::Param parameter;
        std::string paramBaseName(parameterPtr->uTag());
        int j = 1;
        paramName << paramBaseName << j;
        parameter.set(paramName.str(), parsedLine_[linePosition].string_);
        ExtendedString p_orig(parameter.stringValue());
        p_orig.toUpper();
        if (p_orig.possibleParam())
        {
          parameter.setVal(std::string("{" + p_orig + "}"));
          if (!circuitContext_.resolveParameter(parameter))
            parameter.setVal(p_orig);
        }
        parameter.setGiven( true );
        addInstanceParameter(parameter);
        ++linePosition;

        while ( (linePosition < numFields) && (parsedLine_[linePosition].string_ == ",") )
        {
          paramName.str("");
          ++j;
          ++linePosition;
          paramName << paramBaseName << j;
          parameter.set(paramName.str(), parsedLine_[linePosition].string_);
          ExtendedString p_orig(parameter.stringValue());
          p_orig.toUpper();
          if (p_orig.possibleParam())
          {
            parameter.setVal(std::string("{" + p_orig + "}"));
            if (!circuitContext_.resolveParameter(parameter))
              parameter.setVal(p_orig);
          }
          parameter.setGiven( true );
          addInstanceParameter(parameter);
          ++linePosition;
        }
      }
      else if (parameterPtr->stringValue() == "VECTOR-COMPOSITE")
      {
        std::ostringstream paramName;
        // Get the components from metadata.
        std::vector<Device::Param> components;
        metadata_.getInstanceCompositeComponents(
            getNetlistDeviceType(),
            parameterPtr->uTag(), modelLevel,
            components);

        if (DEBUG_IO) {
          Xyce::dout() << " Processing composite for parameter " << parameterPtr->uTag() << std::endl;

          for (int ieric=0; ieric < components.size(); ++ieric)
          {
            Xyce::dout() << "tag[" << ieric << "] = " << components[ieric].uTag() << "  val = " << components[ieric].stringValue() << std::endl;
          }
        }

        // Do some error checking.
        // find the right curly brace, if it exists.
        // This is not a perfect test.  It just finds the first one.
        bool foundRightB = false;
        int rightBracePosition = linePosition;
        int ipos1;

        if (DEBUG_IO)
          Xyce::dout() << "    Doing error checking for vector-composite..." << std::endl;

        for (ipos1=linePosition;ipos1<numFields;++ipos1)
        {
          if (DEBUG_IO)
            Xyce::dout() << "    position..." << ipos1 << " string is "
                         << parsedLine_[ipos1].string_ ;

          if (parsedLine_[ipos1].string_ == "}")
          {
            foundRightB = true;
            rightBracePosition = ipos1;
            if (DEBUG_IO)
              Xyce::dout() << "    found right brace!" ;
          }
          if (DEBUG_IO)
            Xyce::dout() << std::endl;
        }

        if (DEBUG_IO)
          Xyce::dout() << "String at line position is " << parsedLine_[linePosition].string_ << std::endl ;

        if (parsedLine_[linePosition].string_ != "{" || !foundRightB)
        {
          // Expect components of VECTOR-COMPOSITE to be enclosed in braces.
          Report::UserError().at(netlistFileName_, parsedLine_[linePosition].lineNumber_)
            << "Composite parameter " << parameterPtr->uTag() << " must be enclosed in braces {}";
        }
        ++linePosition;

        Device::Param parameter( "", "" );
        std::string paramBase(parameterPtr->uTag());
        int numComponents = 0;
        while (parsedLine_[linePosition].string_ != "}")
        {
          ExtendedString component ( parsedLine_[linePosition].string_ );
          component.toUpper();
          std::vector<Device::Param>::iterator paramIter =
            std::find(components.begin(), components.end(), Device::Param(component, ""));
          if (paramIter == components.end())
          {
            Report::UserError().at(netlistFileName_, parsedLine_[linePosition].lineNumber_)
              << "Found unexpected component \"" << component << "\" for composite parameter " << paramBase;
          }

          if (DEBUG_IO)
            Xyce::dout() << " Found component " << paramIter->uTag() << std::endl;

          linePosition += 2;

          // Mark the component as given in components. Later we will
          // add all of the components and their defaults that were
          // not given.
          paramIter->setGiven(true);

          int j = 0;
          bool startVC = true;
          while (startVC || parsedLine_[linePosition].string_ == ",")
          {
            if (!startVC) ++linePosition;
            startVC = false;
            paramName << paramBase << j << "." << component;

            ExtendedString value ( parsedLine_[linePosition].string_ );
            // Commented out by TVR on 9 Nov 2010
            // See similar comment in ParameterBlock::extractModelData
            // It is unreasonable to up-case parameter values here, because
            // the possibility exists that the string parameter value is
            // intended as a file name, and this will break on any system with
            // a case-dependent file system, such as Linux or BSD or any other
            // *nix other than Mac.
            // The responsibility for up-casing a parameter value should be
            // left to the user of the parameter if necessary, never done here.
            //
            //value.toUpper();

            if (value == "DEFAULT")
            {
              // parameter.set(paramName.str(), static_cast<Util::Param &>(*paramIter));
              setParam(parameter, paramName.str(), *paramIter);
            }
            else
            {
              parameter.set(paramName.str(), value);
            }

            parameter.setGiven(true);
            addInstanceParameter(parameter);
            paramName.str("");

            ++linePosition;
            ++j;
          } // end of while loop.

          numComponents = j;
        } // end of while loop.

        // Now find the components not given in the netlist and
        // add them with their defaults.
        parameter.setGiven(false);
        std::vector<Device::Param>::iterator paramIter = components.begin();
        std::vector<Device::Param>::iterator paramEnd = components.end();
        for (; paramIter != paramEnd; ++paramIter)
        {
          if (!paramIter->given())
          {
            for (int j = 0; j < numComponents; ++j)
            {
              paramName << paramBase << j << "." << paramIter->uTag();
              // parameter.set(paramName.str(), static_cast<Util::Param &>(*paramIter));
              setParam(parameter, paramName.str(), *paramIter);
              addInstanceParameter(parameter);
              paramName.str("");
            }
          }
        } // end of paramIter loop
        // Advance past the right bracket, so that it will not be an
        // unrecognized instance parameter.
        ++linePosition;

      } // end of VECTOR-COMPOSITE if statement.
    }
    else //  parameterPtr != NULL
    {
      if ( action == "break" )
      {
        // Break out of loop and return if a field is encountered
        // that is not an instance parameter. Only care about
        // parameterEndPosition in this case.
        parameterEndPosition = linePosition - 1;
        break;
      }
      else if ( action == "continue" )
      {
        // Continue to the next field if a field is encountered that
        // is not an instance parameter.
        ++linePosition;
      }
    }
  }

  // Set parameterEndPosition if we get here with action = "break".
  if ( action == "break" )
  {
    parameterEndPosition = linePosition - 1;
  }
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::addDefaultInstanceParameters
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 08/11/2003
//----------------------------------------------------------------------------
void DeviceBlock::addDefaultInstanceParameters(int modelLevel)
{
  addInstanceParameters(metadata_.getInstanceParameters(getNetlistDeviceType(), modelLevel));
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::issueUnrecognizedParameterError
// Purpose       : Use the Xyce error manager to issue an error message to
//                 the user and exit the program.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/01/2001
//-----------------------------------------------------------------------------
void DeviceBlock::issueUnrecognizedParameterError(
   std::string const& parameterName)
{
  Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
    << "Unrecognized parameter " << parameterName << " for device " << getName();
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::setParameterValues
// Purpose       : Look for expression valued parameters in the parameter
//                 list, evaluate expressions found and reset the parameter
//                 value accordingly.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/31/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::setParameterValues()
{
  Device::Param parameter( "", "" );
  int numParameters = getNumberOfInstanceParameters();
  int i;
  std::vector<std::string> strings;
  std::vector<std::string> funcs;

  for ( i = 0; i < numParameters; ++i )
  {
    parameter = getInstanceParameter(i);
    if ( parameter.hasExpressionValue() || parameter.isQuoted() )
    {
      if (!circuitContext_.resolveParameter(parameter))
      {
        std::string msg("Parameter " + parameter.uTag() + " for device ");
        msg += getName() + " contains unrecognized symbol";
        if (parameter.getType() == Xyce::Util::EXPR)
        {
          Util::Expression &e = parameter.getValue<Util::Expression>();
          std::vector<std::string>::iterator s;

          strings.clear();
          funcs.clear();
          e.get_names(XEXP_STRING, strings);
          e.get_names(XEXP_FUNCTION, funcs);
          if (strings.size() + funcs.size() == 1)
            msg += ":";
          else if (strings.size() + funcs.size() > 1)
            msg += "s:";
          for (s=strings.begin() ; s != strings.end() ; ++s)
            msg += " " + *s;
          for (s=funcs.begin() ; s != funcs.end() ; ++s)
            msg += " " + *s + "()";
        }
        else
        {
          msg += "(s): " + parameter.stringValue();
        }
        if (strings.size() + funcs.size() > 0)
        {
          Report::UserError().at(netlistFileName_, parsedLine_[0].lineNumber_)
            << msg;
        }
      }
      setInstanceParameter( i, parameter );
    }
    else
    {
      if ( parameter.getType() == Xyce::Util::STR && !(parameter.isNumeric()) )
      {
        if (Util::possibleParam(parameter.stringValue()))
        {
          ExtendedString p_orig(parameter.stringValue()); p_orig.toUpper();

          parameter.setVal(std::string("{" + p_orig + "}"));
          if (!circuitContext_.resolveParameter(parameter))
            parameter.setVal(p_orig);
        }
      }
    }
  }

  return true;
}



//-----------------------------------------------------------------------------
// Function      : DeviceBlock::findInstanceParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2001
//-----------------------------------------------------------------------------
Device::Param* DeviceBlock::findInstanceParameter( Util::Param const& parameter )
{
  std::vector<Device::Param>::iterator paramIter;

  paramIter = find( deviceData_.getDevBlock().params.begin(),
                    deviceData_.getDevBlock().params.end(),
                    parameter );
  if ( paramIter != deviceData_.getDevBlock().params.end() )
  {
    return &(*paramIter);
  }
  else
  {
    return NULL;
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceBlock::findInstanceParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2001
//-----------------------------------------------------------------------------
Device::Param* DeviceBlock::findInstanceParameter(
      std::string const& parameterName )
{
  Device::Param parameter( parameterName, "" );

  return findInstanceParameter( parameter );
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::getInstanceParameters
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/10/2003
//----------------------------------------------------------------------------
void DeviceBlock::getInstanceParameters(
    std::vector<Device::Param>& parameters)
{
  parameters.reserve(deviceData_.getDevBlock().params.size());

  parameters.insert(parameters.end(),
                    deviceData_.getDevBlock().params.begin(),
                    deviceData_.getDevBlock().params.end());
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractMIDeviceData
// Purpose       : Extract the device data from parsedLine for YMIL/YMIN device
//
//                  Y%type%name inductor_list coupling_list model params
//
//                  type :- MIL | MIN
//                  for linear and nonlinear coupling respectively
//
//                  name :- K1_K2_...KN
//                  linear couplings with shared inductors have Knames
//                  concatenated w/underscores otherwise the Kname is used
//
//                  inductor_list :- Lname T1 T2 I
//                  a list of one or more inductors comprising the name,
//                  terminals, and inductance
//
//                  coupling_list :- L1 L2 ... LN C
//                  a list of inductor names followed by coupling value
//                  linear couplings with shared inductors will have several
//                  of these
//
//                  model :- the model name for nonlinear coupling
//
//                  params
//
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceBlock::extractMIDeviceData()
{
  // The device name has been extracted by extractData. The
  // remaining data in parsedLine is extracted here.
  int numFields = parsedLine_.size();

  // Extract the model name from parsedLine if one exists. If
  // a model name was found, find its type.
  int modelLevel, modelNamePosition;
  std::string modelType;
  bool modelFound = extractModelName( modelType, modelLevel, modelNamePosition );
  if( !modelFound ) modelLevel = -1;

  // Some devices require a model, check that a model was found if it is
  // required for this device.
  if( metadata_.isModelRequired(getNetlistDeviceType(), modelLevel ) &&
   !modelFound )
  {
    Report::UserError() << "Did not find required model for device " << getName();
    return false;
  }

  if (modelFound)
  {
    // Query the metadata for the validity of the device model.  Note that the
    // YMIL/YMIN model metadata is extracted from the K model metadata
    if( !metadata_.isModelTypeValid(getNetlistDeviceType(), modelType, modelLevel ) )
    {
      Report::UserError() << "Model type \"" << modelType << "\" not valid for device " << getName() << " of type " << getNetlistDeviceType();
      return false;
    }
  }

  std::set< std::string > inductors;
  Device::Param param;

  // parameters names differ between YMIL and YMIN
  std::string tmpDiff;
  if( modelFound )
  {
    param.setTag( "NONLINEARCOUPLING" );
    param.setGiven( true );
    param.setVal( 1 );
    addInstanceParameter( param );
    tmpDiff = "COUPLEDMutIndNonLin";
  }
  else
  {
    tmpDiff = "COUPLEDMutIndLin";
  }

  // Extract the device nodes from parsedLine_.
  std::list<tagged_param> nodeValues;
  int nodesFound = 0, i = 1;

  // stop when first inductor name is seen again == start of coupling list
  std::string couplingListFlag(parsedLine_[i].string_);

  do
  {
    // add inductor to mutual inductance
    param.setTag( tmpDiff );
    param.setVal( parsedLine_[i].string_ );
    param.setGiven(true);
    addInstanceParameter( param );
    inductors.insert( parsedLine_[i].string_ );

    // store the two node values and add to parameter list
    nodeValues.push_back( tagged_param( ExtendedString(
     parsedLine_[i + 1].string_ ).toUpper(), 0 ) );
    param.setTag( "NODE1" );
    param.setVal( parsedLine_[i + 1].string_ );
    param.setGiven(true);
    addInstanceParameter( param );

    nodeValues.push_back( tagged_param( ExtendedString(
     parsedLine_[i + 2].string_ ).toUpper(), 0 ) );
    param.setTag( "NODE2" );
    param.setVal( parsedLine_[i + 2].string_ );
    param.setGiven(true);
    addInstanceParameter( param );

    nodesFound += 2;

    // add inductance value for this inductor
    param.setTag( "COUPLEDINDUCTANCE" );
    param.setVal( atof( ( parsedLine_[i + 3].string_ ).c_str() ) );
    param.setGiven(true);
    addInstanceParameter( param );

    // move to next inductor in the list
    i += 4;
  } while ( ( i < numFields ) && ( parsedLine_[i].string_ != couplingListFlag ) );

  // shouldn't need to do this so I'll comment it out.
  // remove duplicate node names
  // nodeValues.unique();

  // set device node names
  setNodeValues( nodeValues );
  deviceData_.getDevBlock().numExtVars = nodesFound;

  bool moreCouplings = ( parsedLine_[i].string_ == couplingListFlag );
  int curr = i;
  ++i;

  // retrieve coupling list
  while( moreCouplings && ( i < numFields ) )
  {
    if( inductors.find( parsedLine_[i].string_ ) == inductors.end() )
    {
      // add coupling list
      param.setTag( "COUPLING" );
      std::string &  coup = ((parsedLine_[i].string_));
      if (Util::isValue(coup))
      {
        param.setVal(Util::Value(coup));
      }
      else
      {
        std::string exp("{"+ parsedLine_[i].string_ + "}");
        param.setVal(exp);
      }

      param.setGiven(true);
      addInstanceParameter( param );

      for( ; curr < i; ++curr )
      {
        // add inductor to mutual coupling list
        param.setTag( "COUPLEDINDUCTOR" );
        param.setVal( parsedLine_[curr].string_ );
        param.setGiven(true);
        addInstanceParameter( param );
      }

      // move indices to next list coupling value
      ++curr;
      ++i;

      // check for beginning of another coupling list
      moreCouplings = ( ( i + 1 < numFields ) &&
       ( inductors.find( parsedLine_[curr].string_ ) != inductors.end() ) );
    }

    // move to next item on the list
    ++i;
  }


  // Add the device instance parameters and their default values
  // to instanceParameters and check parsedLine for instance
  // parameters.
  int pStart;
  ( modelFound ) ? pStart = i : pStart = curr;

  int pEnd = numFields - 1 ;

  // extract remaining instance parameters on the device line
  extractInstanceParameters( pStart, pStart, pEnd, "break", modelLevel );

  // Issue warning if there are more fields on the line.
  if( pEnd + 1 < numFields )
  {
    Report::UserError() << "Unrecognized fields for device " << getName();
    return false;
  }

  return true;  // Only get here on success.
}


//-----------------------------------------------------------------------------
// Function      : DeviceBlock::isValidDeviceType
// Purpose       : Given a device type (a letter) return false
//                 if that letter represents  an  illegal, unimplemented
//                 device type
// Special Notes : Having added this function, there is now one extra place
//                 that has to be updated when adding a new device letter.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 06/20/2008
//-----------------------------------------------------------------------------
bool DeviceBlock::isValidDeviceType(const std::string & deviceType)
{

  bool retcode;
  if (deviceType == "A" ||
      deviceType == "N" ||
      deviceType == "P" ) 
  {
    retcode=false;
  }
  else
  {
    retcode=true;
  }

  return (retcode);
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::getNodeValue
// Purpose       :
// Special Notes : It is assumed that getNumberOfNodes was called to ensure
//                 that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
const std::string &DeviceBlock::getNodeValue( int i ) const
{
  std::list<tagged_param>::const_iterator paramIter;
  paramIter = (deviceData_.getNodeBlock().get_NodeList()).begin();
  for ( int j = 0; j < i; j++ )
  {
    paramIter++;
  }
  return paramIter->tag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::getAllNodeNames
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceBlock::getAllNodeNames( std::list< std::string > & nodeNames )
{
  std::list<tagged_param>::const_iterator paramIter, paramEnd;

  paramIter = ( deviceData_.getNodeBlock().get_NodeList() ).begin();
  paramEnd  = ( deviceData_.getNodeBlock().get_NodeList() ).end();

  nodeNames.clear();

  for ( ; paramIter != paramEnd; ++paramIter )
  {
    nodeNames.push_back( paramIter->tag );
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::setNodeValue
// Purpose       :
// Special Notes : It is assumed getNumberOfNodes was called to ensure
//                 that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
void DeviceBlock::setNodeValue(int const& i, std::string const& nodeValue )
{
  std::list<tagged_param>::const_iterator paramIter;
  paramIter = deviceData_.getNodeBlock().get_NodeList().begin();
  for ( int j = 0; j < i; j++ )
  {
    paramIter++;
  }
  tagged_param* paramPtr = const_cast<tagged_param*> (&(*paramIter));
  paramPtr->tag = nodeValue;
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::setNodeValues
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/12/2003
//----------------------------------------------------------------------------
void DeviceBlock::setNodeValues(std::list<tagged_param> const& nodeValues)
{
  std::list<tagged_param>::const_iterator n=nodeValues.begin() ;
  std::list<tagged_param>::const_iterator nEnd=nodeValues.end() ;
  for ( ; n != nEnd ; ++n)
  {
    checkNode(n->tag);
  }

  deviceData_.getNodeBlock().set_NodeList(nodeValues);
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::checkNode
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
void
DeviceBlock::checkNode (std::string const& n)
{

/* DEBUG_ELR_PROFILE - FIXME remove this block; mostly done in pass 2

  CircuitContext* myContext = circuitContext_.getCurrentContextPtr();

  if (!(myContext->devMap.empty()))
  {
    if (myContext->devMap.find(n) != myContext->devMap.end())
    {
      std::string msg("Device and node share the same name: " + n);
      if (myContext->getName() != "")
	msg += " in subcircuit: " + myContext->getName();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
    }
  }

  if ( !(myContext->modMap.empty()) )
  {
    if (myContext->modMap.find(n) != myContext->modMap.end())
    {
      std::string msg("Model and node share the same name: " + n);
      if (myContext->getName() != "")
	msg += " in subcircuit: " + myContext->getName();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
    }
  }

*/

}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::setInstanceParamter
// Purpose       :
// Special Notes : It assumed that getNumberOfInstanceParameters was called
//                 to ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
void
DeviceBlock::setInstanceParameter(
  int                   i,
  Device::Param &       parameter )
{
  // deviceData_.getDevBlock().params[i].setVal(static_cast<Util::Param &>(parameter));
  setParamValue(deviceData_.getDevBlock().params[i], parameter);
  deviceData_.getDevBlock().params[i].setGiven( parameter.given() );
}

} // namespace IO
} // namespace Xyce
