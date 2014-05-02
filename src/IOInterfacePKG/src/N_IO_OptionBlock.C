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
// Filename       : $RCSfile: N_IO_OptionBlock.C,v $
//
// Purpose        : Define the N_IO_OptionBlock class an instantiation of
//                  which is associated with a netlist .PARAM or .OPTIONS
//                  line.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/10/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.158.2.4 $
//
// Revision Date  : $Date: 2014/03/13 21:51:25 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#include <iterator>
#include <iostream>
#include <algorithm>
#include <set>

#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_IO_CircuitMetadata.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_Report.h>
#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OptionBlock::OptionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/17/2006
//-----------------------------------------------------------------------------
OptionBlock::OptionBlock( N_IO_CircuitMetadata & md)
: metadata_(md)
{
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::OptionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
OptionBlock::OptionBlock(
    std::string const& fileName,
    std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedInputLine,
    CircuitMetadata & md)
: netlistFileName_(fileName),
  parsedLine(parsedInputLine),
  metadata_(md)
{
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::OptionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
OptionBlock::OptionBlock( std::string const& fileName,
 CircuitMetadata & md )
 : netlistFileName_( fileName ),
   metadata_ (md)
{
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::OptionBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
OptionBlock::OptionBlock(OptionBlock const& rhsOB)
  : netlistFileName_(rhsOB.netlistFileName_),
    parsedLine(rhsOB.parsedLine),
    optionData(rhsOB.optionData),
    metadata_(rhsOB.metadata_)
{
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/17/2006
//-----------------------------------------------------------------------------
OptionBlock & OptionBlock::operator=(const OptionBlock & right)
{
  metadata_ = right.metadata_;

  parsedLine = right.parsedLine;
  optionData = right.optionData;

  netlistFileName_ = right.netlistFileName_;

  return *this;
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::printDiagnostic
// Purpose       : Output the details of a device block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
void OptionBlock::printDiagnostic() const
{
  Xyce::dout() << std::endl;
  Xyce::dout() << "Option Information" << std::endl;
  Xyce::dout() << "------------------" << std::endl;

  if (getName() != "")
  {
    Xyce::dout() << "input line:" << std::endl;
    unsigned int size = parsedLine.size();
    for (unsigned int i = 0; i < size; ++i)
    {
      Xyce::dout() << "  " << parsedLine[i].string_;
    }
    Xyce::dout() << std::endl;
    Xyce::dout() << "  name: " << getName() << std::endl;
  }

  Xyce::dout() << "  parameters: " << std::endl;
  int numParameters = getNumberOfParameters();
  for (int i = 0; i < numParameters; ++i)
  {
    Xyce::dout() << "  " << getParameter(i).tag() << "  ";
    Xyce::dout() << getParameter(i).stringValue() << std::endl;
  }
  Xyce::dout() << std::endl << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractData
// Purpose       : Determine option type and extract the data appropriately.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 05/22/2002
//-----------------------------------------------------------------------------
bool OptionBlock::extractData()
{
  bool result;
  int line;

  ExtendedString lineType ( parsedLine[0].string_ );
  line = parsedLine[0].lineNumber_;
  lineType.toUpper();
  bool optionsForSensitivity=false;

  setNetlistLocation(netlistFileName_, line);

  ExtendedString optionName("");
  if (parsedLine.size() > 1)
  {
    optionName = parsedLine[1].string_ ;
    optionName.toUpper();
  }

  if (lineType == ".TR")
  {
    lineType = ".TRAN"; // TR is a synonym for TRAN
  }

  // Determine the line type and extract data appropriately.
  if ( lineType == ".OP" )
  {
    result = extractOPData();
  }
  else if ( lineType == ".OPTIONS" )
  {
    result = extractOptionsData();
  }
  else if ( lineType == ".SENS" )
  {
    result = extractSENSData();
    optionsForSensitivity=true;
  }
  else if ( lineType == ".DCOP" )
  {
    result = extractDCOPData();
  }
  else if ( lineType == ".OUTPUT" )
  {
    result = extractOutputData();
  }
  else if ( lineType == ".DC" )
  {
    result = extractDCData();
  }
  else if ( lineType == ".STEP" )
  {
    result = extractSTEPData();
  }
  else if ( lineType == ".PARAM" || lineType == ".GLOBAL_PARAM")
  {
    result = extractParamData();
  }
  else if ( lineType == ".PRINT" )
  {
    result = extractPrintData();
  }
  else if ( lineType == ".TRAN" )
  {
    result = extractTRANData();
  }
  else if ( lineType == ".MPDE" )
  {
    result = extractMPDEData();
  }
  else if ( lineType == ".HB" )
  {
    result = extractHBData();
  }
  else if ( lineType == ".AC" )
  {
    result = extractACData();
  }
  else if ( lineType == ".MOR" )
  {
    result = extractMORData();
  }
  else if ( lineType == ".MEASURE" || lineType == ".MEAS" )
  {
    result = extractMEASUREData();
  }
  else if ( lineType == ".FOUR" )
  {
    result = extractFOURIERData();
  }
  else if ( lineType == ".FFT" )
  {
    Report::UserError0().at(netlistFileName_, line)
      << "OptionBlock::extractData():  .FFT option not supported.";
  }
  else if ( lineType == ".RESULT" )
  {
    result = extractRESULTData();
  }
  else if ( lineType == ".OBJECTIVE" )
  {
    result = extractOBJECTIVEData();
  }
  else if ( lineType == ".IC" || lineType == ".DCVOLT" )
  {
    result = extractICData();
  }
  else if ( lineType == ".NODESET" )
  {
    result = extractNodeSetData();
  }
  else if ( lineType == ".SAVE" )
  {
    result = extractSaveData();
  }
  else if ( lineType == ".LOAD" )
  {
    result = extractLoadData();
  }

  if ( !result )
  {
    return false;
  }

  // ERK: This extra bit is needed to restore ".options sens" to work correctly.
  // It was broken by the error trap, below, that was added to patch
  // bug 1017.  So, OPTIONS will still trigger the trap below, unless the
  // specific option name is SENS.
  if (lineType == ".OPTIONS" && optionName=="SENSITIVITY")
  {
    optionsForSensitivity=true;
  }

  // This block of code is to handle bug 1017.  In most
  // cases, attempts to use expressions on .options lines will result in
  // the value being set incorrectly to zero.
  if (lineType != ".PARAM" && lineType != ".GLOBAL_PARAM" &&
      lineType != ".RESULT" && lineType != ".OBJECTIVE" && 
      lineType != ".TRAN" && lineType != ".MEASURE"
      && !(optionsForSensitivity)
     )
  {
    std::list<Util::Param>::const_iterator paramIter, paramEnd;
    paramEnd = optionData.getParams().end();
    paramIter = optionData.getParams().begin();
    for ( ; paramIter != paramEnd ; ++paramIter)
    {
      if ((*paramIter).hasExpressionValue())
      {
        Report::UserError0().at(netlistFileName_, line)
          << "Expressions are not supported for " << lineType;
      }
    }
  }

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractOPData
// Purpose       : Extract the parameters from a netlist .OP line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool OptionBlock::extractOPData()
{
  int numFields = parsedLine.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields > 1 )
  {
    Report::UserWarning0().at(netlistFileName_, parsedLine[0].lineNumber_) << "Ignoring extra fields on .OP line";
  }

  // Set the OptionBlock name.
  setName( "OP" );

  return true; // Only get here on success.
}

//----------------------------------------------------------------------------
// Function       : OptionBlock::extractDCOPData
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 05/08/06
//----------------------------------------------------------------------------
bool OptionBlock::extractDCOPData()
{
  int numFields = parsedLine.size();
  int linePosition = 0;
  std::set<std::string> pars;

  // Set the OptionBlock name.
  setName( "OP_IO" );

  Util::Param parameter("", "");
  ExtendedString ES("");

  ++linePosition;
  while (linePosition < numFields)
  {
    ES = parsedLine[linePosition].string_ ;
    ES.toUpper ();
    if (pars.find(ES) != pars.end())
    {
      Report::UserError().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "Multiple " << ES << " name specifications in .DCOP line";
    }
    if (ES == "INPUT" || ES == "OUTPUT" )
    {
      pars.insert(ES);
      parameter.setTag( ES );
      if (++linePosition >= numFields) return true;
      if (parsedLine[linePosition].string_ == "=")
        if (++linePosition >= numFields) return true;
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;
    }
    else
    {
      if (linePosition == 1 && numFields == 2)
      {
        pars.insert("INPUT");
        parameter.setTag( "INPUT" );
        parameter.setVal( parsedLine[linePosition].string_ );
        addParameter( parameter );
        pars.insert("OUTPUT");
        parameter.setTag( "OUTPUT" );
        parameter.setVal( parsedLine[linePosition].string_ );
        addParameter( parameter );
        ++linePosition;
      }
      else
      {
        Report::UserError().at(netlistFileName_, parsedLine[0].lineNumber_)
          << "Unrecognized field in .DCOP line";
      }
    }
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : OptionBlock::extractOutputData
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/26/2003
//----------------------------------------------------------------------------
bool OptionBlock::extractOutputData()
{
  int numFields = parsedLine.size();

  // Check that the number of fields is as expected.
  if (numFields != 3)
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".OUTPUT line requires exactly two parameters, the beginning time and the output step size.";
  }

  // Set the OptionBlock name.
  setName("OUTPUT-LINE");

  // Extract the output parameters
  Util::Param time( "TIME", "" );
  time.setVal(parsedLine[1].string_);
  addParameter(time);

  Util::Param interval( "INTERVAL", "" );
  interval.setVal(parsedLine[2].string_);
  addParameter(interval);

  // Further processing is needed by the CircuitBlock which initiated
  // this method. For convenience, we add the line number of the .OUTPUT
  // line as a parameter in case the CircuitBlock cannot find the
  // .OPTIONS OUTPUT line that should have preceded this .OUTPUT line.
  Util::Param lineNum( "LINENUMBER", static_cast<int> (parsedLine[0].lineNumber_) );
  addParameter(lineNum);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractParamData
// Purpose       : Extract the parameters from a netlist .PARAM line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool OptionBlock::extractParamData()
{
  int numFields = parsedLine.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 4 )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".param line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.

  Util::Param parameter("", "");
  while ( linePosition < numFields )
  {
    parameter.setTag( parsedLine[linePosition].string_ );

    if ( (parameter.uTag() == "TEMP") || (parameter.uTag() == "VT") ||
         (parameter.uTag() == "GMIN") || (parameter.uTag() == "TIME") )
    {
      Report::UserError0().at(netlistFileName_, parsedLine[linePosition].lineNumber_)
        << "Parameter name " << parameter.uTag() << " is not permitted";
    }

    if ( linePosition + 2 >= numFields )
    {
      Report::UserError0().at(netlistFileName_, parsedLine[linePosition].lineNumber_)
        << "Unexpectedly reached end of line while looking for parameters in .PARAM statement";
    }

    if ( parsedLine[ linePosition+1].string_ != "=" )
    {
      Report::UserError0().at(netlistFileName_, parsedLine[linePosition+1].lineNumber_)
        << "Equal sign (=) required between parameter and value in .PARAM statement";
    }

    linePosition += 2;   // Advance to parameter value field.

    if (parsedLine[linePosition].string_[0] == '{')
    {
      Util::Expression expPtr(parsedLine[linePosition].string_);

      std::vector<std::string> junk;
      std::string msg;
      expPtr.get_names(XEXP_NODE, junk);
      if (junk.size() > 0)
        msg = "Node Voltage";
      expPtr.get_names(XEXP_INSTANCE, junk);
      if (junk.size() > 0)
        msg = "Device Current";
      expPtr.get_names(XEXP_LEAD, junk);
      if (junk.size() > 0)
        msg = "Lead Current";
      if (msg != "")
      {
        Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
          << msg << " may not be used in parameter expression (" << parameter.tag() << ")";
      }
      parameter.setVal(parsedLine[linePosition].string_);
//DNS: It should be this?      parameter.setVal(expPtr);
//        parameter.setVal(expPtr);
    }
    else
    {
      ExtendedString tmp ( parsedLine[linePosition].string_ );
      if (tmp.possibleParam())
        parameter.setVal(std::string("{" + parsedLine[linePosition].string_ + "}"));
//DNS: this expr should be parsed here?
//        parameter.setVal(Util::Expression(std::string("{" + parsedLine[linePosition].string_ + "}")));
      else
        parameter.setVal(parsedLine[linePosition].string_);
    }

    addParameter( parameter );

    ++linePosition;     // Advance to next parameter.

    if ( linePosition < numFields && parsedLine[linePosition].string_ == "," )
    {
      // Skip past comma separator.
      ++linePosition;
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractOptionsData
// Purpose       : Extract the parameters from a netlist .OPTIONS line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool OptionBlock::extractOptionsData()
{
  int numFields = parsedLine.size();

  // The type of options is given by the second field on the .options line
  // UNLESS they are PSPICE style options. Check the second field and set the
  // name attribute appropriately.
  ExtendedString optionName ( parsedLine[1].string_ );
  optionName.toUpper();
  int parameterStartPos = 2;

  if(metadata_.isOptionsPackage(optionName))
  {
    setName( optionName );
  }
  else
  {
    // Check to see if this is PSPICE style option line or if an unrecognized
    // package name was given. Do this by checking for an "=" in either the
    // 3rd or 4th position on the line.
    if (numFields > 2 && parsedLine[2].string_ == "=" )
    {
      // look for TEMP and/or TNOM  If they are the only tags present on the 
      // rest of the line, they try setting the package name to DEVICE 
      //
      // assume we only found TEMP and TNOM for now.
      bool foundOnlyTempAndTnom=true;
      
      // format of parsed line at this point is <tag> = <value> 
      // with the first tag starting at position 1 (hence why 
      // parsedLine[2].string_ == "="  was true to get into this if block.)
      // Check every 3rd item starting from i=1 to only check the tags. 
      for (int i=1; i<numFields; i=i+3 )
      {
        ExtendedString tag ( parsedLine[i].string_ );
        tag.toUpper();
        if( (tag != "TNOM") && (tag != "TEMP") )
        {
          foundOnlyTempAndTnom=false;
        }
      }
      
      if(foundOnlyTempAndTnom)
      {
        setName("DEVICE");
        parameterStartPos = 1;
      }
      else
      {
        setName( "PSPICE" );
        parameterStartPos = 1;

        Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
          << "PSPICE style options are not currently supported in Xyce";

        return false;
      }
    }
    else
    {
      Report::UserError0().at(netlistFileName_, parsedLine[1].lineNumber_)
        << "Unrecognized .OPTIONS package " << optionName;
      return false;
    }
  }

  // Create an option block to temporarily store the default options.
  OptionBlock defaultOptions("", parsedLine, metadata_ );

  // Get the default options from metadata.
  defaultOptions.addDefaultOptionsParameters( getName() );

  // Extract the parameters from parsedLine.
  int parameterEndPos = numFields - 1;
  std::vector<Util::Param> inputParameters;
  Util::Param parameter("", "");
  int intervalParameterStart = -1;
  int i = parameterStartPos;
  std::string paramBaseName;
  while (i <= parameterEndPos-1)
  {
    // Check for equal sign.
    if ( parsedLine[i+1].string_ != "=" )
    {
      if ( optionName != "OUTPUT" && optionName != "RESTART" )
      {
        Report::UserError0().at(netlistFileName_, parsedLine[i].lineNumber_)
           << "Equal sign required between parameter name and value in .OPTIONS " << getName();
        return false;
      }
      else
      {
        // Stop after the tagged parameters have been extracted
        // from a .OPTIONS RESTART or .OPTIONS OUTPUT line, they
        // will be handled later.
        intervalParameterStart = i;
        break;
      }
    }

    // Extract parameter name and value from parsedLine and add to
    // parameter list. Check to see if the parameter is "VECTOR"
    // valued and treat accordingly.
    parameter.set( parsedLine[i].string_, "" );
    Util::Param* parameterPtr = defaultOptions.findParameter(parameter);
    if (parameterPtr == NULL)
    {
      Report::UserWarning0().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "No options parameter " << parameter.tag() << " found, parameter will be ignored.";
      i+= 3;
    }
    else if (parameterPtr->stringValue() != "VECTOR")
    {
      parameter.setVal( parsedLine[i+2].string_ );
      inputParameters.push_back( parameter );
      i+= 3;

      if (i < parameterEndPos-1 && parsedLine[i].string_ == ",")
      {
        Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
          << "Options parameter " << parameter.tag() << " is flagged as not VECTOR, but has comma in value.";
        return false;
      }
    }
    else
    {
      // We have a vector valued parameter.
      // Name the jth component of the parameter of the vector by appending
      // "j" to the parameter name.
      std::ostringstream paramName;
      paramBaseName = ExtendedString(parsedLine[i].string_).toUpper();
      int j = 1;

      paramName << paramBaseName << j;
      i += 2;
      parameter.set(paramName.str(), parsedLine[i].string_);
      addParameter(parameter);

      // This while loop is dangerous if we are near the end of the
      // parsed line because it still assumed that the format is
      // option=value,value and not option=value,value,
      // that is no trailing comma.  It does work if
      // option=value,value is at the end of a line now (RLS 5/07)
      int testSize = parsedLine.size()-1;
      while ((i < testSize) && (parsedLine[i+1].string_ == ",") )
      {
        paramName.str("");
        ++j;
        paramName << paramBaseName << j;
        i += 2;
        parameter.set(paramName.str(), parsedLine[i].string_);
        addParameter(parameter);
      }

      ++i;
    }
  }

  // For each input parameter, check that it is in the default
  // set and if so, set its value in "parameters" to the input
  // value, otherwise flag it as an unknown parameter.
  int numInputParameters = inputParameters.size();
  for ( int k = 0; k < numInputParameters; ++k )
  {
    Util::Param* parameterPtr =
      defaultOptions.findParameter( inputParameters[k] );
    if ( parameterPtr != NULL )
    {
      parameterPtr->setVal( inputParameters[k] );
      addParameter( *parameterPtr );
    }
    else
    {
      Report::UserWarning0().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "No options parameter " << inputParameters[k].tag() << " found, parameter will be ignored.";
    }
  }

  // If this is a ".OPTIONS OUTPUT" line or ".OPTIONS RESTART" line
  // then get the time and interval pairs.
  if ( (optionName == "OUTPUT" || optionName == "RESTART") &&
        intervalParameterStart != -1)
  {
    Util::Param time( "TIME", "" );
    Util::Param interval( "INTERVAL", "" );

    int i = intervalParameterStart;
    while ( i < parameterEndPos )
    {
      time.setVal( parsedLine[i].string_ );
      addParameter( time );
      ++i;

      interval.setVal( parsedLine[i].string_ );
      addParameter( interval );
      ++i;
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::addDefaultOptionsParameters
// Purpose       : Add the default parameters for a model from metadata.
// Special Notes : While not required, it is generally expected that
//                 the parameters given in the input netlist will have
//                 already been extracted via extractModelData().
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/17/2001
//-----------------------------------------------------------------------------
void OptionBlock::addDefaultOptionsParameters(const std::string &optionName )
{
  addParameters(metadata_.getOptionsParameters(optionName));
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::extractDCData
// Purpose       : Determine number of sweep variables on this DC line
//                : and create option blocks for each, storing the blocks
//                : in the referenced vector.
//-----------------------------------------------------------------------------
int OptionBlock::extractDCData( std::vector< OptionBlock > & oBs )
{
  // length of the original .DC line
  int numFields = parsedLine.size();

  // number of sweep sources on this line, and index to current source
  int sourcesFound = 0;

  // start of parameters (skip over the ".DC")
  int linePosition = 1;

  // line can be variable in length, with muliple sources to a line
  while( linePosition < numFields )
  {
    // create a new option block for this source
    oBs.push_back( OptionBlock( netlistFileName_, metadata_ ) );
    oBs[sourcesFound].setName( "DC" );

    Util::Param parameter("", "");

    std::string stringVal = parsedLine[linePosition + 1].string_;
    Util::toUpper(stringVal);

    ExtendedString curr ( parsedLine[linePosition].string_ );
    curr.toUpper();

    if( "LIST" == stringVal || "LIST" == curr )
    {
      // sweep type is LIST
      parameter.setTag( "TYPE" );
      parameter.setVal( "LIST" );
      oBs[sourcesFound].addParameter( parameter );

      // get sweep variable name and move to just before beginning of list
      parameter.setTag( "PARAM" );
      ( "LIST" == curr ) ?
       parameter.setVal( stringVal ) :
       parameter.setVal( curr );
      oBs[sourcesFound].addParameter( parameter );
      ++linePosition;

      // collect values in this list until next sweep variable name is found
      while( ++linePosition < numFields &&
       ( Util::Param( "", parsedLine[linePosition].string_ ) ).isNumeric() )
      {
        parameter.setTag( "VAL" );
        parameter.setVal( parsedLine[linePosition].string_ );
        oBs[sourcesFound].addParameter( parameter );
      }
    }

    else
    {
      std::string sweepStepTag ("STEP");

      // check for non-LIST default (LINear sweep)
      if( ( Util::Param( "", stringVal ) ).isNumeric() )
      {
        stringVal = "LIN";
      }

      // non-LIST sweep type is given
      else
      {
        // get sweep type name and move to sweep variable name
        stringVal = parsedLine[linePosition++].string_;
//        stringVal.toUpper();
        Util::toUpper(stringVal);

        // change sweep step tag
        sweepStepTag = "NUMSTEPS";
      }

      // Add the type (which was determined above) to the parameter list.
      parameter.setTag( "TYPE" );
      parameter.setVal( stringVal );
      oBs[sourcesFound].addParameter( parameter );

      // simple check for expected four items; semantic errors pass through
      if( numFields <= linePosition + 3 )
      {
        Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
          << ".DC line not formatted correctly, found unexpected number of fields";
        linePosition = numFields;
      }
      else {
        parameter.setTag( "PARAM" );
        parameter.setVal( parsedLine[linePosition++].string_ );
        oBs[sourcesFound].addParameter( parameter );

        parameter.setTag( "START" );
        parameter.setVal( parsedLine[linePosition++].string_ );
        oBs[sourcesFound].addParameter( parameter );

        parameter.setTag( "STOP" );
        parameter.setVal( parsedLine[linePosition++].string_ );
        oBs[sourcesFound].addParameter( parameter );

        parameter.setTag( sweepStepTag );
        parameter.setVal( parsedLine[linePosition++].string_ );
        oBs[sourcesFound].addParameter( parameter );
      }
    }

    // record this source (and move on to the next)
    ++sourcesFound;
  }

  return sourcesFound;
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractDCData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool OptionBlock::extractDCData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "DC" );

  // Check that the minimum required number of fields are on the line.
  if ( (numFields-1)%4 != 0 )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".DC line not formatted correctly, found unexpected number of fields";
  }
  else {
    int linePosition = 1;   // Start of parameters on .param line.

    Util::Param parameter("", "");
    while ( linePosition < numFields )
    {
      parameter.setTag( "VSOURCE" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "VSTART" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "VSTOP" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "VSTEP" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.
    }
  }
  
  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractSTEPData
// Purpose       : Extract the parameters from a netlist .STEP line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/30/2003
//-----------------------------------------------------------------------------
bool OptionBlock::extractSTEPData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "STEP" );

  // First check if the type has been explicitly set.
  // If not, set it to the default, LIN.
  int pos1=1;

  bool typeExplicitSetLinDecOct = false;
  bool typeExplicitSetList = false;
  std::string type("LIN");
  int typeIndex = -1;
  while ( pos1 < numFields )
  {
    ExtendedString stringVal ( parsedLine[pos1].string_ );
    stringVal.toUpper ();
    if (stringVal == "LIN" ||
        stringVal == "DEC" ||
        stringVal == "OCT")
    {
      typeExplicitSetLinDecOct = true;
      type = stringVal;
      typeIndex = pos1;
    }
    else if (stringVal == "LIST")
    {
      typeExplicitSetList = true;
      type = stringVal;
      typeIndex = pos1;
    }

    ++pos1;
  }

  // Check that the minimum required number of fields are on the line.
  int offset = 1;
  if (typeExplicitSetLinDecOct)
  {
    offset = 2;
  }

  if (!typeExplicitSetList)// if this is a list, number of fields is arbitrary.
  {
    if ( (numFields-offset)%4 != 0 )
    {
      Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
        << ".STEP line not formatted correctly.";
    }
  }

  int linePosition = 1;   // Start of parameters on .param line.
  Util::Param parameter("", "");

  // Add the type (which was determined above) to the parameter list.
  parameter.setTag( "TYPE" );
  parameter.setVal( type );
  addParameter( parameter );

  if (type=="LIN")
  {
    if (typeExplicitSetLinDecOct) linePosition=2;
    while ( linePosition < numFields )
    {
      parameter.setTag( "PARAM" );
      parameter.setVal( ExtendedString(parsedLine[linePosition].string_).toUpper() );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "START" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "STOP" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "STEP" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.
    }
  }
  else if (type=="DEC")
  {
    if (typeExplicitSetLinDecOct) linePosition=2;

    while ( linePosition < numFields )
    {
      parameter.setTag( "PARAM" );
      parameter.setVal( ExtendedString(parsedLine[linePosition].string_).toUpper() );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "START" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "STOP" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "NUMSTEPS" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.
    }
  }
  else if (type=="OCT")
  {
    if (typeExplicitSetLinDecOct) linePosition=2;

    while ( linePosition < numFields )
    {
      parameter.setTag( "PARAM" );
      parameter.setVal( ExtendedString(parsedLine[linePosition].string_).toUpper() );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "START" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "STOP" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.

      parameter.setTag( "NUMSTEPS" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;     // Advance to next parameter.
    }

  }
  else if (type=="LIST")
  {
    parameter.setTag( "PARAM" );
    parameter.setVal( ExtendedString(parsedLine[1].string_).toUpper() );
    addParameter( parameter );

    int linePosition=3;
    while (linePosition<numFields)
    {
      parameter.setTag( "VAL" );
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;
    }
  }
  else
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".STEP line contains an unrecognized type";
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractMPDEData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool OptionBlock::extractMPDEData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "MPDE" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 || numFields > 6 )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".MPDE line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields - 1;

  Util::Param parameter("", "");

  // TSTEP and TSTOP are required, get them now.
  parameter.setTag( "TSTEP" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );
  ++linePosition;     // Advance to next parameter.

  parameter.setTag( "TSTOP" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );
  ++linePosition;     // Advance to next parameter.

  // Next check last field to see if it is UIC or NOOP.
  parameter.setTag( parsedLine[endPosition].string_ );
  if ( parameter.uTag() == "NOOP" || parameter.uTag() == "UIC" )
  {
    parameter.setVal( "1" );
    addParameter( parameter );
    --endPosition;
  }
  else if ( numFields == 6)
  {
    Report::UserError0().at(netlistFileName_, parsedLine[endPosition].lineNumber_)
      << "expected NOOP/UIC field on .MPDE line but found" << parameter.usVal();
  }

  if ( linePosition <= endPosition )
  {
    parameter.setTag( "TSTART" );
    parameter.setVal( parsedLine[linePosition].string_ );
    addParameter( parameter );
    ++linePosition;     // Advance to next parameter.
  }

  if ( linePosition <= endPosition )
  {
    parameter.setTag( "DTMAX" );
    parameter.setVal( parsedLine[linePosition].string_ );
    addParameter( parameter );
    ++linePosition;     // Advance to next parameter.
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractHBData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool OptionBlock::extractHBData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "HB" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 2 )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".HB line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields;

  Util::Param parameter("", "");

// frequency of oscillation is required
  std::vector<double> freqs(numFields - 1);

  int i = 0;
  while( linePosition < endPosition )
  {
    std::string & value = parsedLine[linePosition].string_;
    if (Util::isValue(value))
    {
      freqs[i] = Util::Value(value);
    }
    else
    {
      Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "Attempt to assign value for FREQ from " << value;
    }
    ++linePosition;
    ++i;
  }

  parameter.setTag( "FREQ" );
  parameter.setVal( freqs );
  addParameter( parameter );
//  ++linePosition;

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractTRANData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool OptionBlock::extractTRANData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "TRAN" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 || numFields > 6 )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".TRAN line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields;

  Util::Param parameter("", "");

  // TSTEP and TSTOP are required, get them now.
  parameter.setTag( "TSTEP" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );
  ++linePosition;     // Advance to next parameter.

  parameter.setTag( "TSTOP" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );
  ++linePosition;     // Advance to next parameter.

  // at this point we can have
  // [ <TSTART> [<DTMAX>] ] [NOOP|UIC] [{expression for time step schedule}]
  //
  // because these are optional, we'll need to test the type of each
  // parameter before deciding how to handle it
  //
  // need these flags to figure out what any untaged doubles are on the .tran line
  bool tstartFound=false;
  bool dtmaxFound=false;
  while( linePosition < endPosition )
  {
    // convert the next item to a parameter and set the
    // parameter's tag and value.  We set both because
    // in some cases we need to change the tag and in
    // others we need to change the value
    parameter.setTag( parsedLine[linePosition].string_ );
    parameter.setVal( parsedLine[linePosition].string_ );

    if(parameter.hasExpressionValue())
    {
      // found a max time step schedule expression
      parameter.setTag("MAXTIMEEXPRESSION");
      addParameter( parameter );
    }
    else
    {
      // could be TSTART DTMAX or NOOP|UIC
      if ( parameter.uTag() == "NOOP" || parameter.uTag() == "UIC" )
      {
        parameter.setVal( "1" );
        addParameter( parameter );
      }
      else
      {
        if( !tstartFound )
        {
          parameter.setTag( "TSTART" );
          addParameter( parameter );
          tstartFound=true;
        }
        else if( !dtmaxFound )
        {
          parameter.setTag( "DTMAX" );
          addParameter( parameter );
          dtmaxFound=true;
        }
        else
        {
          Report::UserError0().at(netlistFileName_, parsedLine[linePosition].lineNumber_)
            << "expected NOOP/UIC field on .TRAN line but found" << parameter.usVal();
        }
      }
    }
    linePosition++;
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractPrintData
// Purpose       : Extract the parameters from a netlist .PRINT line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 05/02/2002
//-----------------------------------------------------------------------------
bool OptionBlock::extractPrintData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "PRINT" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    Xyce::dout() << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << std::endl;
  }
#endif

  // Set the TYPE and add it the parameters.
  std::string tmpString = parsedLine[1].string_;
  Util::toUpper(tmpString);
  if (tmpString == "TR")
  {
    tmpString = "TRAN"; // TR is a synonym for TRAN
  }
  Util::Param typeParameter("TYPE", tmpString);
  if ( typeParameter.usVal() != "DC" &&
       typeParameter.usVal() != "TRAN" &&
       typeParameter.usVal() != "AC" &&
       typeParameter.usVal() != "HB" &&
       typeParameter.usVal() != "MOR" )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[1].lineNumber_)
      << "Invalid type \"" << typeParameter.usVal() << "\" for .PRINT statement";
  }

  // Add the options parameter set to the model.
  addDefaultOptionsParameters( "PRINT" );

  // Reset the default TYPE with the value found.
  Util::Param* parameterPtr = findParameter( typeParameter );

  if( parameterPtr == NULL )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[1].lineNumber_)
      << "Failed to find parameter\"" << typeParameter.usVal() << "\" in optionData for .PRINT statement";
  }

  parameterPtr->setVal( typeParameter.usVal() );

  // Set no default value for the output file.
  // if the user set FILE=vale it will be picked up here.  Otherwise
  // the output manager will use the top level simulation netlist name
  // as the default.
  Util::Param fileParameter("FILE", "");
  parameterPtr = findParameter( fileParameter );

  if( parameterPtr == NULL )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[1].lineNumber_)
      << "Failed to find parameter\"" << typeParameter.usVal() << "\" in optionData for .PRINT statement";
  }

  parameterPtr->setVal( fileParameter.stringValue() );

  //
  // Note: The PRINT control parameters (FORMAT, WIDTH, FILE, ...) must be
  // tagged on the print line. Tagged parameters are "better" since they can
  // be identifed and parsed based soley on metadata with generic code.
  // The code that used to allow untagged parameters has been removed.

  // Check for tagged parameters.
  int parameterStartPos  = 2;
  int position = parameterStartPos;
  if ( numFields > parameterStartPos + 1 && parsedLine[parameterStartPos + 1].string_ == "=" )
  {
     // Tagged parameters found.
     Util::Param* parameterPtr;
     while ( position+1 < parsedLine.size() &&
             parsedLine[position+1].string_ == "=" )
     {

       parameterPtr = findParameter( Util::Param(parsedLine[position].string_, "") );
       if ( parameterPtr != NULL )
       {
         if (parameterPtr->tag() != "FILE")
           parameterPtr->setVal(
               ExtendedString(parsedLine[position+2].string_ ).toUpper());
         else
           parameterPtr->setVal(ExtendedString(parsedLine[position+2].string_));
       }
       else
       {
         Report::UserWarning0().at(netlistFileName_, parsedLine[0].lineNumber_)
           << "No PRINT parameter " << parsedLine[position].string_ << " found, parameter will be ignored.";
         return false;
       }

       position += 3;
     }
  }
  else if (Xyce::DEBUG_IO) // There is no point to this warning unless we're debugging
  {
    Report::UserWarning0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << "No tagged parameters found";
  }

  // Complete the PRINT line parsing.
  // Some of the remaining fields are of the
  // form I(Vname), V(node) or V(node1,node2), and these fields need
  // special treatment.
  Util::Param parameter;
  ExtendedString field("");
  std::ostringstream msg;
  int p_err;
  while ( position < numFields )
  {
    if (position+1 < numFields && parsedLine[position+1].string_ == "(")
    {
      if (parsedLine[position].string_ == "I" || parsedLine[position].string_ == "i"
        || (parsedLine[position].string_.size() == 2 &&
        (parsedLine[position].string_[0] == 'I' || parsedLine[position].string_[0] == 'i')))
      {
        if ((position+3 < numFields && parsedLine[position+3].string_ == ")") ||
            (position+4 < numFields && parsedLine[position+4].string_ == ")"))
        {
          field = parsedLine[position].string_;
          field.toUpper();
          parameter.setTag(field);
          parameter.setVal(1.0);
          addParameter( parameter );

          if (parsedLine[position+3].string_ == ")")
          {
            field = parsedLine[position+2].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );
            addParameter( parameter );

            position += 4;
          }
          else
          {
            field = parsedLine[position+2].string_ + " " + parsedLine[position+3].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );
            addParameter( parameter );

            position += 5;
          }
        }
        else
        {
          msg << "Unrecognized current specification";
          p_err = position;
        }
      }
      else if( parsedLine[position].string_ == "V" || parsedLine[position].string_ == "v"
        || ( (parsedLine[position].string_.size() == 2 || parsedLine[position].string_.size() == 3) &&
             (parsedLine[position].string_[0] == 'V' || parsedLine[position].string_[0] == 'v')))
      {
        if( parsedLine[position+3].string_ == ")" )
        {
          field = parsedLine[position].string_;
          field.toUpper();
          parameter.setTag(field);
          parameter.setVal( 1.0 );
          addParameter( parameter );

          field = parsedLine[position+2].string_;
          field.toUpper();
          parameter.setTag( field );
          parameter.setVal( 0.0 );
          addParameter( parameter );

          position += 4;
        }
        else if( parsedLine[position+5].string_ == ")" )
        {
          field = parsedLine[position].string_;
          field.toUpper();
          parameter.setTag(field);
          parameter.setVal( 2.0 );
          addParameter( parameter );

          field = parsedLine[position+2].string_;
          field.toUpper();
          parameter.setTag( field );
          parameter.setVal( 0.0 );
          addParameter( parameter );

          field = parsedLine[position+4].string_;
          field.toUpper();
          parameter.setTag( field );
          parameter.setVal( 0.0 );
          addParameter( parameter );

          position += 6;
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else if( parsedLine[position].string_ == "N" || parsedLine[position].string_ == "n")
      {
        if( parsedLine[position+3].string_ == ")" )
        {
          parameter.setTag("N");
          parameter.setVal( 1.0 );
          addParameter( parameter );

          field = parsedLine[position+2].string_;
          field.toUpper();
          parameter.setTag( field );
          parameter.setVal( 0.0 );
          addParameter( parameter );

          position += 4;
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else
      {
        msg << "Unrecognized parenthetical specification";
        p_err = position;
      }
    }
    else // This is either an expression, or a STEP parameter.
    {
      if (parsedLine[position].string_ == "(" || parsedLine[position].string_ == ")")
      {
        msg << "Unrecognized parenthesis";
        p_err = position;
      }
      field = parsedLine[position].string_;
      field.toUpper();
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      addParameter ( parameter );
      ++position;
    }
    if (!msg.str().empty())
    {
      msg << " in .print near ";
      position = p_err+4;
      p_err -= 2;
      if (p_err<0)
        p_err = 0;
      if (position >= numFields)
        position = numFields;
      while (p_err < position)
      {
        msg << parsedLine[p_err].string_ + " ";
        ++p_err;
      }
      Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
        << msg.str();
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractRESULTData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 08/29/2004
//-----------------------------------------------------------------------------
bool OptionBlock::extractRESULTData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "RESULT" );

#ifdef Xyce_DEBUG_IO
  // print out the parsed line
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    Xyce::dout() << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << std::endl;
  }
#endif

  int linePosition = 1;   // Start of parameters on .param line.
  Util::Param parameter("", "");

  if (parsedLine.size () <= 1)
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << "Too few fields in the .RESULT line";
    return false;
  }

  parameter.setTag( "EXPRESSION" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );

  if (linePosition != (parsedLine.size () - 1))
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << "Too many fields in the .RESULT line";
    return false;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractACData
// Purpose       : Extract the parameters from a netlist .AC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool OptionBlock::extractACData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "AC" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 4 || numFields > 5)
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".AC line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields - 1;

  Util::Param parameter("", "");

  // type is required
  parameter.setTag( "TYPE" );
//  parameter.setVal( parsedLine[linePosition].string_ );
  ExtendedString stringVal ( parsedLine[linePosition].string_ );
  stringVal.toUpper();
  parameter.setVal(stringVal);
  addParameter( parameter );

  ++linePosition;     // Advance to next parameter.

  // np is required
  parameter.setTag( "NP" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );
  ++linePosition;     // Advance to next parameter.

  // fstart is required
  parameter.setTag( "FSTART" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );
  ++linePosition;     // Advance to next parameter.


  // fstop is required
  parameter.setTag( "FSTOP" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );
//  ++linePosition;     // Advance to next parameter.

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractMORData
// Purpose       : Extract the parameters from a netlist .MOR line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 5/25/12
//-----------------------------------------------------------------------------
bool OptionBlock::extractMORData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "MOR" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 2 )
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << ".MOR line has an unexpected number of fields";
  }

  int linePosition = 1;   // Start of parameters on .param line.

  Util::Param parameter("", "");

  // type is required
  parameter.setTag( "SIZE" );
//  parameter.setVal( parsedLine[linePosition].string_ );
  ExtendedString stringVal ( parsedLine[linePosition].string_ );
  stringVal.toUpper();
  parameter.setVal(stringVal);
  addParameter( parameter );

  std::vector<std::string> portNames(numFields - 2);
  for (int i=0; i<(numFields-2); ++i)
  {
    ++linePosition;     // Advance to next parameter.

    // Insert next port name
    ExtendedString stringVal ( parsedLine[linePosition].string_ );
    stringVal.toUpper();
    portNames[i] = stringVal;
  }

  parameter.setTag( "PORTLIST" );
  parameter.setVal( portNames );
  addParameter( parameter );

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractOBJETIVEData
// Purpose       : Extract the parameters from a netlist .OBJECTIVE line 
//
//-----------------------------------------------------------------------------
bool OptionBlock::extractOBJECTIVEData()
{
  int numFields = parsedLine.size();
  int linePosition = 0;
  std::set<std::string> pars;

  // Set the OptionBlock name.
  setName( "OBJECTIVE" );

  Util::Param parameter("", "");
  ExtendedString ES("");

  ++linePosition;
  while (linePosition < numFields-1)
  {
    ES = parsedLine[linePosition].string_ ;
    ES.toUpper ();
    if (pars.find(ES) != pars.end())
    {
      Report::UserError().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "Multiple " << ES << " name specifications in .OBJECTIVE line";
    }
    if (ES == "FILE" || ES == "FUNCTION" || ES == "VALUE" ||
        ES == "WEIGHT" || ES == "NAME")
    {
      pars.insert(ES);
      parameter.setTag( ES );
      if (++linePosition >= numFields) return true;
      if (parsedLine[linePosition].string_ == "=")
        if (++linePosition >= numFields) return true;
      parameter.setVal( parsedLine[linePosition].string_ );
      if (ES == "VALUE" || ES == "WEIGHT")
      {
        if (!parameter.hasExpressionValue ())
        {
          Report::UserError().at(netlistFileName_, parsedLine[0].lineNumber_)
            << "Non-Expression specified in .OBJECTIVE argument: " << ES << " = " << parsedLine[linePosition].string_;
        }
      }
      addParameter( parameter );
      ++linePosition;
    }
// MATCH was an idea for mapping names between the data file and the netlist.
// This has not seemed useful in initial testing, so it is commented out, and will
// not be implemented unless a need arises
//  else if (ES == "MATCH" || parameter.tag() == "MATCH")
//  {
//    if (ES == "MATCH")
//    {
//      if (++linePosition >= numFields) return true;
//      parameter.setTag( "MATCH" );
//    }
//    std::string::size_type first=parsedLine[linePosition].string_.find_first_of(':');
//    if (first != parsedLine[linePosition].string_.find_last_of(':') ||
//        first == std::string::npos)
//    {
//      std::string msg("Match specification '");
//      msg += parsedLine[linePosition].string_;
//      msg += "' in .OBJECTIVE does not contain a single colon";
//      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
//      netlistFileName_, parsedLine[0].lineNumber_);
//    }
//    parameter.setVal( parsedLine[linePosition].string_ );
//    addParameter( parameter );
//    ++linePosition;
//  }
    else
    {
      if (linePosition == 1)
      {
        pars.insert("NAME");
        parameter.setTag( "NAME" );
        parameter.setVal( parsedLine[linePosition].string_ );
        addParameter( parameter );
        ++linePosition;
      }
      else
      {
        Report::UserError().at(netlistFileName_, parsedLine[0].lineNumber_)
          << "Unrecognized field in .OBJECTIVE line";
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractMEASUREData
// Purpose       : Convert a .measure line to an options block
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool OptionBlock::extractMEASUREData()
{
  // .measure has a lot of optonal and complex syntax.
  //
  // .measure mode name type output_var... options...
  //
  // where mode = AC, DC or TRAN
  //       name = any text string for a name (can't be AC, DC or TRAN)
  //       type = keywords that let us figure out what type of measurement is needed
  //              TRIG and/or TARG = RiseFallDelay
  //              AVG, MAX, MIN, PP, RMS, INTEG = Statistics
  //              FIND/WHEN = FindWhen
  //              PARAM = Equation evaluation
  //              DERIVATIVE or DERIV = Derivative
  //              INTEGRAL = Integral
  //              ERROR = Error
  //              FOUR = Fourier analysis (similar to .FOUR)
  //       output_var = simulation variale to be measured.  This can be any of the following
  //              v(a), v(a,b), v(a)=number v(a)=v(b), i(a) ix(a) or an expression.  
  //              So it really could be one or more
  //              variables and or a real number.  The comparison is always equity.
  //       options = these are keywords=value pairs that set limits on when or how
  //              the measurement is done. Value is usually a number, but can be a string in
  //              at least one case.
  //

  // we will use these sets to test for measure types and keywords.
  std::set<std::string> typeSet;
  typeSet.insert( std::string("TRIG") );
  typeSet.insert( std::string("TARG") );
  typeSet.insert( std::string("AVG") );
  typeSet.insert( std::string("MAX") );
  typeSet.insert( std::string("MIN") );
  typeSet.insert( std::string("PP") );
  typeSet.insert( std::string("RMS") );
  typeSet.insert( std::string("FREQ") );
  typeSet.insert( std::string("FIND") );
  typeSet.insert( std::string("WHEN") );
  typeSet.insert( std::string("PARAM") );
  typeSet.insert( std::string("EQN") );
  typeSet.insert( std::string("DERIVATIVE") );
  typeSet.insert( std::string("DERIV") );
  typeSet.insert( std::string("DUTY") );
  typeSet.insert( std::string("INTEGRAL") );
  typeSet.insert( std::string("INTEG") );
  typeSet.insert( std::string("ERROR") );
  typeSet.insert( std::string("ON_TIME") );
  typeSet.insert( std::string("OFF_TIME") );
  typeSet.insert( std::string("FOUR") );

  std::set<std::string> keywords;
  std::set<std::string> simpleKeywords;
  std::set<std::string> numOrTextKeywords;

  simpleKeywords.insert( std::string("TD") );
  simpleKeywords.insert( std::string("GOAL") );
  simpleKeywords.insert( std::string("WEIGHT") );
  simpleKeywords.insert( std::string("MINVAL") );
  simpleKeywords.insert( std::string("AT") );
  simpleKeywords.insert( std::string("FROM") );
  simpleKeywords.insert( std::string("TO") );
  simpleKeywords.insert( std::string("IGNORE") );
  simpleKeywords.insert( std::string("YMIN") );
  simpleKeywords.insert( std::string("YMAX") );
  simpleKeywords.insert( std::string("ON") );
  simpleKeywords.insert( std::string("OFF") );
  simpleKeywords.insert( std::string("FRAC_MAX") );
  simpleKeywords.insert( std::string("MIN_THRESH") );
  simpleKeywords.insert( std::string("MAX_THRESH") );
  simpleKeywords.insert( std::string("NUMFREQ") );
  simpleKeywords.insert( std::string("GRIDSIZE") );
  simpleKeywords.insert( std::string("DEFAULT_VAL") );
  simpleKeywords.insert( std::string("INDEPVARCOL") ); 
  simpleKeywords.insert( std::string("INDEPVAR2COL") );
  simpleKeywords.insert( std::string("DEPVARCOL") );
 
  numOrTextKeywords.insert( std::string("RISE") );
  numOrTextKeywords.insert( std::string("FALL") );
  numOrTextKeywords.insert( std::string("CROSS") );
  numOrTextKeywords.insert( std::string("FILE") );
  numOrTextKeywords.insert( std::string("COMP_FUNCTION") );
  // make a union for the keywords set
  set_union( simpleKeywords.begin(), simpleKeywords.end(), numOrTextKeywords.begin(),
             numOrTextKeywords.end(), std::inserter<std::set<std::string> >(keywords, keywords.begin()) );

  int numFields = parsedLine.size();

  if( numFields < 4 )
  {
    Report::UserError().at(netlistFileName_,parsedLine[0].lineNumber_)
      << "Too few items on .MEASURE line.  Need at lest .MEASURE <mode> <name> <type>";
  }

  // Set the OptionBlock name.
  setName( "MEASURE" );

  Util::Param parameter;

  // look for the fixed items
  ExtendedString currentWord(parsedLine[1].string_);
  currentWord.toUpper();
  if (currentWord == "TR")
  {
    currentWord = "TRAN"; // TR is a synonym for TRAN
  }

  if( currentWord == "DC" || currentWord == "TRAN" || currentWord == "AC" )
  {
    parameter.set("MODE", currentWord);
    addParameter(parameter);
  }
  else
  {
    Report::UserError().at(netlistFileName_, parsedLine[1].lineNumber_)
      << "Unknown mode in .MEASURE line.  Should be DC, AC or TRAN/TR";
  }

  currentWord = parsedLine[2].string_;
  currentWord.toUpper();
  if (currentWord == "TR")
  {
    currentWord = "TRAN"; // TR is a synonym for TRAN
  }
  if( currentWord != "DC" && currentWord != "TRAN" && currentWord != "AC" )
  {
    parameter.set("NAME", currentWord);
    addParameter(parameter);
  }
  else
  {
    Report::UserError().at(netlistFileName_, parsedLine[2].lineNumber_)
      << "Illegal name in .MEASURE line.  Cannot be AC, DC or TRAN/TR ";
  }

  currentWord = parsedLine[3].string_;
  currentWord.toUpper();
  if( typeSet.find( currentWord ) != typeSet.end() )
  {
    parameter.set("TYPE", currentWord);
    addParameter(parameter);
  }
  else
  {
    Report::UserError().at(netlistFileName_, parsedLine[3].lineNumber_)
      << "Illegal type in .MEASURE line.  Must be one of: TRIG, TARG, AVG, MAX, MIN, PP, RMS, INTEG, FIND, WHEN, PARAM, DERIVATIVE, DERIV, INTEGRAL, ERROR, FOUR";
  }

  // already got MEASURE <DC|AC|TRAN> name TYPE
  // now try and parse of the rest of the line catching keywords when they're found
  int position = 4;
  while ( position < numFields )
  {
    currentWord = parsedLine[position].string_;
    currentWord.toUpper();
    // The type tag, TARG can get repeated, so check for it with other params.
    if( typeSet.find( currentWord ) != typeSet.end() )
    {
      parameter.set("TYPE", currentWord);
      addParameter(parameter);
    }
    else if( numOrTextKeywords.find( currentWord ) != numOrTextKeywords.end() )
    {
      // these can be in the form of TAG=value or TAG=keyword where keyword = LAST
      // and potentially the "=" is optional.
      if( ((position+1) < numFields) || ((parsedLine[(position+1)].string_ == "=") && ((position+2) < numFields)) )
      {
        // value is in next or next + 1 position
        int valPosition = ++position;
        if( parsedLine[valPosition].string_ == "=" )
        {
          valPosition = ++position;
        }
        std::string & value = parsedLine[valPosition].string_;
        if( Util::isInt(value) )
        {
          parameter.set(currentWord, Util::Ival(value) );
        }
        else if( Util::isValue(value) )
        {
          int valAsInt = static_cast<int>(Util::Value(value));
          parameter.set(currentWord, valAsInt);
        }
        else
        {
          parameter.set(currentWord, ExtendedString (value).toUpper());
        }
        addParameter(parameter);
      }
      else
      {
        Report::UserError().at(netlistFileName_, parsedLine[position].lineNumber_)
          << "Incomplete .MEASURE line.  RISE, FALL, CROSS must be followed by a value or LAST";
      }

    }
    else if( simpleKeywords.find( currentWord ) != simpleKeywords.end() )
    {
      // these are in the form of TAG=value
      // and potentially the "=" is optional.
      if( ((position+1) < numFields) || ((parsedLine[(position+1)].string_ == "=") && ((position+2) < numFields)) )
      {
        // value is in next or next + 1 position
        int valPosition = ++position;
        if( parsedLine[valPosition].string_ == "=" )
        {
          valPosition = ++position;
        }
        std::string & value = parsedLine[valPosition].string_;
        if( Util::isValue(value) )
        {
          parameter.set(currentWord, Util::Value(value) );
        }
        else if( Util::isInt(value) )
        {
          parameter.set(currentWord, Util::Ival(value) );
        }
        else
        {
          Report::UserError().at(netlistFileName_, parsedLine[position].lineNumber_)
            << "Incomplete .MEASURE line.  TD, GOAL, WEIGHT, MINVAL, AT, TO, IGNORE, YMIN, YMAX must be followed by a value";
        }
        addParameter(parameter);
      }
      else
      {
        Report::UserError().at(netlistFileName_, parsedLine[position].lineNumber_)
          << "Incomplete .MEASURE line.  TD, GOAL, WEIGHT, MINVAL, AT, TO, IGNORE, YMIN, YMAX must be followed by a value";
      }
    }
    else
    {
      if( (currentWord[0]=='{') && (currentWord[currentWord.size()-1]=='}') ) 
      {
         parameter.set(currentWord, currentWord);
         addParameter(parameter);
      }
      else 
      {
      // the last form of legal syntax is out_var, out_var=val or out_var1=out_var2
      // first figure out how many positions between here and the next legal keyword or
      // the end of the line
      int endPosition=position;
      bool endNotFound=true;
      ExtendedString nextWord("");
      while( endNotFound )
      {
        endPosition++;
        if( endPosition == numFields)
        {
          endNotFound=false;
        }
        else
        {
          nextWord = parsedLine[endPosition].string_;
          nextWord.toUpper();
          if( (keywords.find( nextWord ) != keywords.end()) || (typeSet.find( nextWord ) != typeSet.end()) )
          {
            // found a keyword or type was found so we're at the end of space to
            // use in finding out_var, out_var, out_var=val or out_var1=out_var2
            endNotFound=false;
          }
        }
      }
      // ok, package up the output vars
      // it will be in the form ( [] indicate optional items )
      // V or I nodeName [nodeName] [ val | V or I nodeName [nodeName] ]
      // or { } expression deliminated

      while( position < endPosition )
      {
        nextWord = parsedLine[position].string_;
        nextWord.toUpper();
        // the second part of this if clause it to ensure we don't catch keywords that start 
        // with I, V or N and mistake them for Ixxx( ) or Vxxx() 
        if( (nextWord[0] == 'I' || nextWord[0] == 'V' || nextWord[0] == 'N') && 
            (simpleKeywords.find( nextWord ) == simpleKeywords.end() ) )
        {
          // need to do a bit of look ahead here to see if this is a V(a) or V(a,b)
          int numNodes = 1;
          if( ((position+3) < endPosition) && (parsedLine[position+3].string_ == ",") )
          {
            numNodes = 2;
          }
          parameter.set(nextWord, numNodes);
          addParameter( parameter );

          if( (position+2) < endPosition )
          {
            position+=2;
            nextWord = parsedLine[position].string_;
            nextWord.toUpper();
            parameter.set( nextWord, 0.0 );
            addParameter( parameter );

            if( ((position+2) < endPosition) && (parsedLine[position+1].string_ == ",") )
            {
              // may have a voltage difference request
              position+=2;
              nextWord = parsedLine[position].string_;
              nextWord.toUpper();
              parameter.set( nextWord, 0.0 );
              addParameter( parameter );
            }
          }
          else
          {
            Report::UserError().at(netlistFileName_, parsedLine[position].lineNumber_)
              << "Error in .MEASURE line.  Could not parse voltage/current variable";
          }

        }
        else if( nextWord.isValue() )
        {
          std::string objValue("OBJVAL");
          parameter.set( objValue, nextWord.Value());
          addParameter( parameter );
        }
        else if( nextWord.isInt() )
        {
          std::string objValue( "OBJVAL");
          parameter.set( objValue, nextWord.Ival());
          addParameter( parameter );
        }
        else if( nextWord.possibleParam() )
        {
          // could be a param or measure name that will be evaluated 
          // during the simulation.  
          std::string objValue( "OBJVAL");
          parameter.set( objValue, nextWord);
          addParameter( parameter );
        }
        position++;
      }
      // reset the position indicator to the end - 1 because
      // we're going to increment it at the end of the while loop
      position = endPosition - 1;
      }

    }
    position++;
  }
#ifdef Xyce_DEBUG_IO
  OptionBlock::printDiagnostic();
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractFOURIERData
// Purpose       : Convert a .four line to an options block
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 06/03/2013
//-----------------------------------------------------------------------------
bool OptionBlock::extractFOURIERData()
{
  // Set the OptionBlock name
  setName( "FOUR" );

  Util::Param parameter;
  ExtendedString nextWord("");
  
  if(parsedLine.size() < 3)
  {
    Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
      << "Error: the .FOUR line requires at least 3 arguments '.FOUR freq ov1 <ov2 ov3 ...>'";
  }

  parameter.setTag("FREQ");
  parameter.setVal(parsedLine[1].string_ );
  addParameter( parameter );

  int position = 2;
  int endPosition = parsedLine.size();
  while (position < endPosition)
  {

    nextWord = parsedLine[position].string_;
    nextWord.toUpper();
    if( nextWord == "I" || nextWord == "V" )
    {
      // need to do a bit of look ahead here to see if this is a V(a) or V(a,b)
      int numNodes = 1;
      if( ((position+3) < endPosition) && (parsedLine[position+3].string_ == ",") )
      {
        numNodes = 2;
      }
      parameter.set(nextWord, numNodes);
      addParameter( parameter );

      if( (position+2) < endPosition )
      {
        position+=2;
        nextWord = parsedLine[position].string_;
        nextWord.toUpper();
        parameter.set( nextWord, 0.0 );
        addParameter( parameter );

        if( ((position+2) < endPosition) && (parsedLine[position+1].string_ == ",") )
        {
          // may have a voltage difference request
          position+=2;
          nextWord = parsedLine[position].string_;
          nextWord.toUpper();
          parameter.set( nextWord, 0.0 );
          addParameter( parameter );
        }  
      }
      else
      {
        Report::UserError().at(netlistFileName_, parsedLine[position].lineNumber_)
          << "Could not parse .FOUR variable";
      }
    }
    else
    {
      Report::UserFatal().at(netlistFileName_, parsedLine[position].lineNumber_)
        << "Could not parse .FOUR variable " << nextWord;
    }
    position+=2;

  }

return true;

}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractICData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/07/2007
//-----------------------------------------------------------------------------
bool OptionBlock::extractICData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "IC" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    Xyce::dout() << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << std::endl;
  }
#endif

  std::ostringstream msg;
  Util::Param parameter;
  ExtendedString field("");
  int p_err;
  int parameterStartPos  = 1;
  int position = parameterStartPos;

  // .IC can have two formats:
  //   (1)   .ic  V(a)=1.0   V(2)=2.0
  //   (2)   .ic  a 1.0   2 2.0
  // Check for case 1 by looking for = in the line.
  bool formatOne(false);

  while ( position < numFields )
  {
    if (parsedLine[position].string_ == "=")
    {
      formatOne = true;
    }
    ++position;
  }

  position = parameterStartPos;

  if (formatOne)
  {
    while ( position < numFields )
    {
      if (position+5 < numFields && parsedLine[position+1].string_ == "(")
      {
        // current variables.  ERK: Note, originally I thought to have it
        // be possible to set both voltages and currents.  However, it
        // appears (I think) that other codes just set voltages.
        std::string & posString = parsedLine[position].string_;

        if (posString == "I" || posString == "i"
         || (posString.size() == 2 &&
             (posString[0] == 'I' || posString[0] == 'i')))
        {
          {
            msg << "Unsupported current specification.";
            p_err = position;
          }
        }
        // voltage variables:
        else if( posString == "V" || posString == "v")
        {
          if( parsedLine[position+3].string_ == ")" )
          {
            parameter.setTag("V");
            parameter.setVal( 1.0 );
            addParameter( parameter );

            // node name:
            field = parsedLine[position+2].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );
            addParameter( parameter );

            // value:
            field = parsedLine[position+5].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );
            addParameter( parameter );

            position += 6;
          }
          else if( parsedLine[position+6].string_ == "," )
          {
            msg << "Voltage differences not supported.";
            p_err = position;
          }
          else
          {
            msg << "Unrecognized parenthetical specification.";
            p_err = position;
          }
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else // This is either an expression, or a STEP parameter.
      {
        msg << "Unrecognized .IC and/or .DCVOLT specification";
        p_err = position;
      }

      if (!msg.str().empty())
      {
        msg << " in .IC and/or .DCVOLT near ";
        position = p_err+4;
        p_err -= 2;
        if (p_err<0)
          p_err = 0;
        if (position >= numFields)
          position = numFields;
        while (p_err < position)
        {
          msg << parsedLine[p_err].string_ + " ";
          ++p_err;
        }
        Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
          << msg.str();
      }
    }
  }
  else
  {
    while ( position < numFields )
    {
      parameter.setTag("V");
      parameter.setVal( 1.0 );
      addParameter( parameter );

      // node name:
      field = parsedLine[position].string_;
      field.toUpper();
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      addParameter( parameter );

      // value:
      field = parsedLine[position+1].string_;
      field.toUpper();
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      addParameter( parameter );

      position += 2;
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractNodeSetData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/07/2007
//-----------------------------------------------------------------------------
bool OptionBlock::extractNodeSetData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "NODESET" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    Xyce::dout() << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << std::endl;
  }
#endif

  std::ostringstream msg;

  Util::Param parameter;
  ExtendedString field("");
  int p_err;
  int parameterStartPos  = 1;
  int position = parameterStartPos;

  // .NODESET can have two formats:
  //   (1)   .ic  V(a)=1.0   V(2)=2.0
  //   (2)   .ic  a 1.0   2 2.0
  // Check for case 1 by looking for = in the line.
  bool formatOne(false);

  while ( position < numFields )
  {
    if (parsedLine[position].string_ == "=")
    {
      formatOne = true;
    }
    ++position;
  }


  position = parameterStartPos;

  if (formatOne)
  {
    while ( position < numFields )
    {
      if (position+5 < numFields && parsedLine[position+1].string_ == "(")
      {
        // current variables:
        std::string & posString = parsedLine[position].string_;

        if (posString == "I" || posString == "i"
         || (posString.size() == 2 &&
             (posString[0] == 'I' || posString[0] == 'i')))
        {
          {
            msg << "Unsupported current specification.";
            p_err = position;
          }
        }
        // voltage variables:
        else if( posString == "V" || posString == "v")
        {
          if( parsedLine[position+3].string_ == ")" )
          {
            parameter.setTag("V");
            parameter.setVal( 1.0 );
            addParameter( parameter );

            field = parsedLine[position+2].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );
            addParameter( parameter );

            field = parsedLine[position+5].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );

            addParameter( parameter );

            position += 6;
          }
          else if( parsedLine[position+6].string_ == "," )
          {
            msg << "Voltage differences not supported.";
            p_err = position;
          }
          else
          {
            msg << "Unrecognized parenthetical specification.";
            p_err = position;
          }
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else // This is either an expression, or a STEP parameter.
      {
        msg << "Unrecognized .NODESET specification";
        p_err = position;
      }

      if (!msg.str().empty())
      {
        msg << " in .NODESET near ";
        position = p_err+4;
        p_err -= 2;
        if (p_err<0)
          p_err = 0;
        if (position >= numFields)
          position = numFields;
        while (p_err < position)
        {
          msg << parsedLine[p_err].string_ + " ";
          ++p_err;
        }
        Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
          << msg.str();
      }
    }
  }
  else
  {
    while ( position < numFields )
    {
      parameter.setTag("V");
      parameter.setVal( 1.0 );
      addParameter( parameter );

      // node name:
      field = parsedLine[position].string_;
      field.toUpper();
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      addParameter( parameter );

      // value:
      field = parsedLine[position+1].string_;
      field.toUpper();
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      addParameter( parameter );

      position += 2;
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractSaveData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool OptionBlock::extractSaveData()
{
  int numFields = parsedLine.size();
  int linePosition = 0;

  // Set the OptionBlock name.
  setName( "SAVE" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    Xyce::dout() << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << std::endl;
  }
#endif

  std::string msg("");

  Util::Param parameter("", "");
  ExtendedString ES("");

  ++linePosition;
  while (linePosition < numFields)
  {
    ES = parsedLine[linePosition].string_ ;
    ES.toUpper ();

    if (ES == "TYPE" || ES == "FILE" || ES == "TIME" || ES == "LEVEL" )
    {
      parameter.setTag( ES );
      if (++linePosition >= numFields) return true;
      if (parsedLine[linePosition].string_ == "=")
        if (++linePosition >= numFields) return true;
      parameter.setVal( parsedLine[linePosition].string_ );
      addParameter( parameter );
      ++linePosition;
    }
    else
    {
      Report::UserError().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "Unrecognized field in .SAVE line";
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);


  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractLoadData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool OptionBlock::extractLoadData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "LOAD" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    Xyce::dout() << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << std::endl;
  }
#endif

  std::string msg("");


  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::extractSENSData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool OptionBlock::extractSENSData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "SENS" );

  int parameterStartPos = 1;

  // Create an option block to temporarily store the default options.
  OptionBlock defaultOptions("", parsedLine, metadata_ );

  // Get the default options from metadata.
  defaultOptions.addDefaultOptionsParameters( getName() );

  // Extract the parameters from parsedLine.
  int parameterEndPos = numFields - 1;
  std::vector<Util::Param> inputParameters;
  Util::Param parameter("", "");
  int intervalParameterStart = -1;
  int i = parameterStartPos;
  std::string paramBaseName;
  while (i <= parameterEndPos-1)
  {
    // Check for equal sign.
    if ( parsedLine[i+1].string_ != "=" )
    {
      // Stop after the tagged parameters have been extracted
      // from a .OPTIONS RESTART or .OPTIONS OUTPUT line, they
      // will be handled later.
      intervalParameterStart = i;
      break;
    }

    // Extract parameter name and value from parsedLine and add to
    // parameter list. Check to see if the parameter is "VECTOR"
    // valued and treat accordingly.
    parameter.set( parsedLine[i].string_, "" );
    Util::Param* parameterPtr = defaultOptions.findParameter(parameter);
    if (parameterPtr == NULL)
    {
      Report::UserWarning0().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "No options parameter " << parameter.tag() << " found, parameter will be ignored.";
      i+= 3;
    }
    else if (parameterPtr->stringValue() != "VECTOR")
    {
      parameter.setVal( parsedLine[i+2].string_ );
      inputParameters.push_back( parameter );
      i+= 3;

      if (i < parameterEndPos-1 && parsedLine[i].string_ == ",")
      {
        Report::UserError0().at(netlistFileName_, parsedLine[0].lineNumber_)
          << "Options parameter " << parameter.tag() << " is flagged as not VECTOR, but has comma in value.";
      }
    }
    else
    {
      // We have a vector valued parameter.
      // Name the jth component of the parameter of the vector by appending
      // "j" to the parameter name.
      std::ostringstream paramName;
      paramBaseName = ExtendedString(parsedLine[i].string_).toUpper();
      int j = 1;

      paramName << paramBaseName << j;
      i += 2;
      parameter.set(paramName.str(), parsedLine[i].string_);
      addParameter(parameter);

      // This while loop is dangerous if we are near the end of the
      // parsed line because it still assumed that the format is
      // option=value,value and not option=value,value,
      // that is no trailing comma.  It does work if
      // option=value,value is at the end of a line now (RLS 5/07)
      int testSize = parsedLine.size()-1;
      while ((i < testSize) && (parsedLine[i+1].string_ == ",") )
      {
        paramName.str("");
        ++j;
        paramName << paramBaseName << j;
        i += 2;
        parameter.set(paramName.str(), parsedLine[i].string_);
        addParameter(parameter);
      }

      ++i;
    }
  }

  // For each input parameter, check that it is in the default
  // set and if so, set its value in "parameters" to the input
  // value, otherwise flag it as an unknown parameter.
  int numInputParameters = inputParameters.size();
  for ( int k = 0; k < numInputParameters; ++k )
  {
    Util::Param* parameterPtr =
      defaultOptions.findParameter( inputParameters[k] );
    if ( parameterPtr != NULL )
    {
      parameterPtr->setVal( inputParameters[k] );
      addParameter( *parameterPtr );
    }
    else
    {
      Report::UserWarning0().at(netlistFileName_, parsedLine[0].lineNumber_)
        << "No options parameter " << inputParameters[k].tag() << " found, parameter will be ignored.";
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::setName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void OptionBlock::setName(const std::string &nameIn)
{
  optionData.setName(nameIn);
}

void OptionBlock::setNetlistLocation(const std::string &path, int line_number) 
{
  optionData.setNetlistLocation(NetlistLocation(path, line_number));
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::addParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void OptionBlock::addParameter( Util::Param const& parameter )
{
  optionData.getParams().push_back( parameter );
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::addParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void OptionBlock::addParameters( std::vector<Util::Param> const& parametersIn )
{
  optionData.getParams().insert( optionData.getParams().end(), parametersIn.begin(), parametersIn.end() );
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::clearParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void OptionBlock::clearParameters()
{
  optionData.getParams().clear();
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::setParamter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void OptionBlock::setParameter( int const& i, Util::Param const& parameter )
{
  if ( i < getNumberOfParameters() )
  {
    std::list<Util::Param>::iterator paramIter;
    paramIter = optionData.getParams().begin();
    for ( int j = 0; j < i; ++j )
    {
      ++paramIter;
    }

    paramIter->setVal( parameter.stringValue() );
  }
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
const std::string &OptionBlock::getName() const
{
  return optionData.getName();
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::findParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
Util::Param* OptionBlock::findParameter( Util::Param const& parameter )
{
  std::list<Util::Param>::iterator paramIter = find( optionData.getParams().begin(), optionData.getParams().end(), parameter );
  if ( paramIter != optionData.getParams().end())
  {
    return &(*paramIter);
  }
  else
  {
    return NULL;
  }
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::findParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
Util::Param* OptionBlock::findParameter( std::string const& parameterName )
{
  Util::Param parameter( parameterName, "");

  return findParameter( parameter );
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::getNumberOfParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
int OptionBlock::getNumberOfParameters() const
{
  return optionData.getParams().size();
}

//-----------------------------------------------------------------------------
// Function      : optionBlock::getParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
Util::Param OptionBlock::getParameter( int const& i ) const
{
  if ( i < getNumberOfParameters() )
  {
    std::list<Util::Param>::const_iterator paramIter;
    paramIter = optionData.getParams().begin();
    for ( int j = 0; j < i; ++j )
    {
      ++paramIter;
    }

    return *paramIter;
  }
  else
  {
    return Util::Param( "", "" );
  }
}

} // namespace IO
} // namespace Xyce
