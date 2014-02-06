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
// Revision Number: $Revision: 1.134.2.5 $
//
// Revision Date  : $Date: 2013/12/03 23:30:12 $
//
// Current Owner  : $Author: rlschie $
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
#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Expression.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::N_IO_OptionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/17/2006
//-----------------------------------------------------------------------------
N_IO_OptionBlock::N_IO_OptionBlock( N_IO_CircuitMetadata & md)
: metadata_(md)
{
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::N_IO_OptionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
N_IO_OptionBlock::N_IO_OptionBlock(
    string const& fileName,
    vector<N_IO_SpiceSeparatedFieldTool::StringToken> const& parsedInputLine,
    N_IO_CircuitMetadata & md)
: netlistFileName_(fileName),
  parsedLine(parsedInputLine),
  metadata_(md)
{
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::N_IO_OptionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_IO_OptionBlock::N_IO_OptionBlock( string const& fileName,
 N_IO_CircuitMetadata & md )
 : netlistFileName_( fileName ),
   metadata_ (md)
{
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::N_IO_OptionBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
N_IO_OptionBlock::N_IO_OptionBlock(N_IO_OptionBlock const& rhsOB)
  : netlistFileName_(rhsOB.netlistFileName_),
    parsedLine(rhsOB.parsedLine),
    optionData(rhsOB.optionData),
    metadata_(rhsOB.metadata_)
{
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 01/17/2006
//-----------------------------------------------------------------------------
N_IO_OptionBlock & N_IO_OptionBlock::operator=(const N_IO_OptionBlock & right)
{
  metadata_ = right.metadata_;

  parsedLine = right.parsedLine;
  optionData = right.optionData;

  netlistFileName_ = right.netlistFileName_;

  return *this;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::printDiagnostic
// Purpose       : Output the details of a device block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
void N_IO_OptionBlock::printDiagnostic() const
{
  cout << endl;
  cout << "Option Information" << endl;
  cout << "------------------" << endl;

  if (getName() != "")
  {
    cout << "input line:" << endl;
    unsigned int size = parsedLine.size();
    for (unsigned int i = 0; i < size; ++i)
    {
      cout << "  " << parsedLine[i].string_;
    }
    cout << endl;
    cout << "  name: " << getName() << endl;
  }

  cout << "  parameters: " << endl;
  int numParameters = getNumberOfParameters();
  for (int i = 0; i < numParameters; ++i)
  {
    cout << "  " << getParameter(i).tag() << "  ";
    cout << getParameter(i).sVal() << endl;
  }
  cout << endl;

  cout << "------------------" << endl;
  cout << endl;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractData
// Purpose       : Determine option type and extract the data appropriately.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 05/22/2002
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractData()
{
  bool result;
  int line;

  ExtendedString lineType ( parsedLine[0].string_ );
  line = parsedLine[0].lineNumber_;
  lineType.toUpper();
  bool optionsForSensitivity=false;

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
    string msg("N_IO_OptionBlock::extractData():  .FFT option not supported.");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
       netlistFileName_, line);
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
      lineType != ".RESULT" && lineType != ".OBJECTIVE" && lineType != ".TRAN"
      && !(optionsForSensitivity)
     )
  {
    list<N_UTL_Param>::const_iterator paramIter, paramEnd;
    paramEnd = optionData.getParams().end();
    paramIter = optionData.getParams().begin();
    for ( ; paramIter != paramEnd ; ++paramIter)
    {
      if ((*paramIter).hasExpressionValue())
      {
        string msg("Expressions are not supported for ");
        msg += lineType;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, line);
      }
    }
  }

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractOPData
// Purpose       : Extract the parameters from a netlist .OP line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractOPData()
{
  int numFields = parsedLine.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields > 1 )
  {
    string msg("ignoring extra fields on .OP line\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  // Set the OptionBlock name.
  setName( "OP" );

  return true; // Only get here on success.
}

//----------------------------------------------------------------------------
// Function       : N_IO_OptionBlock::extractDCOPData
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 05/08/06
//----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractDCOPData()
{
  int numFields = parsedLine.size();
  int linePosition = 0;
  set<string> pars;

  // Set the OptionBlock name.
  setName( "OP_IO" );

  N_UTL_Param parameter("", "");
  ExtendedString ES("");

  ++linePosition;
  while (linePosition < numFields)
  {
    ES = parsedLine[linePosition].string_ ;
    ES.toUpper ();
    if (pars.find(ES) != pars.end())
    {
      string msg("Multiple ");
      msg += ES;
      msg += " name specifications in .DCOP line";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
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
        string msg("Unrecognized field in .DCOP line");
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
      }
    }
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : N_IO_OptionBlock::extractOutputData
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/26/2003
//----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractOutputData()
{
  int numFields = parsedLine.size();

  // Check that the number of fields is as expected.
  if (numFields != 3)
  {
    string msg(".OUTPUT line requires exactly two parameters, the\n");
    msg += "beginning time and the output step size.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  // Set the OptionBlock name.
  setName("OUTPUT-LINE");

  // Extract the output parameters
  N_UTL_Param time( "TIME", "" );
  time.setVal(parsedLine[1].string_);
  addParameter(time);

  N_UTL_Param interval( "INTERVAL", "" );
  interval.setVal(parsedLine[2].string_);
  addParameter(interval);

  // Further processing is needed by the N_IO_CircuitBlock which initiated
  // this method. For convenience, we add the line number of the .OUTPUT
  // line as a parameter in case the N_IO_CircuitBlock cannot find the
  // .OPTIONS OUTPUT line that should have preceded this .OUTPUT line.
  N_UTL_Param lineNum( "LINENUMBER", static_cast<int> (parsedLine[0].lineNumber_) );
  addParameter(lineNum);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractParamData
// Purpose       : Extract the parameters from a netlist .PARAM line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractParamData()
{
  int numFields = parsedLine.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 4 )
  {
    string msg(".param line has an unexpected number of fields\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  int linePosition = 1;   // Start of parameters on .param line.

  N_UTL_Param parameter("", "");
  while ( linePosition < numFields )
  {
    parameter.setTag( parsedLine[linePosition].string_ );

    if ( (parameter.uTag() == "TEMP") || (parameter.uTag() == "VT") ||
         (parameter.uTag() == "GMIN") || (parameter.uTag() == "TIME") )
    {
      string msg("Parameter name " + parameter.uTag());
      msg += " is not permitted\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, parsedLine[linePosition].lineNumber_);
    }

    if ( linePosition + 2 >= numFields )
    {
      // Hit end of line unexpectedly.
      string msg("Unexpectedly reached end of line while looking for\n");
      msg += " parameters in .PARAM statement";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, parsedLine[linePosition].lineNumber_);
    }

    if ( parsedLine[ linePosition+1].string_ != "=" )
    {
      string msg("Equal sign (=) required between parameter and value in");
      msg += ".PARAM statement\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, parsedLine[linePosition+1].lineNumber_);
    }

    linePosition += 2;   // Advance to parameter value field.

    if (parsedLine[linePosition].string_[0] == '{')
    {
      N_UTL_Expression expPtr(parsedLine[linePosition].string_);

      vector<string> junk;
      string msg;
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
        msg += " may not be used in parameter expression (";
        msg += parameter.tag();
        msg += ")";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
            netlistFileName_, parsedLine[0].lineNumber_);
      }
      parameter.setVal(parsedLine[linePosition].string_);
//DNS: It should be this?      parameter.setVal(expPtr);
//        parameter.setVal(expPtr);
    }
    else
    {
      ExtendedString tmp ( parsedLine[linePosition].string_ );
      if (tmp.possibleParam())
        parameter.setVal(string("{" + parsedLine[linePosition].string_ + "}"));
//DNS: this expr should be parsed here?
//        parameter.setVal(N_UTL_Expression(string("{" + parsedLine[linePosition].string_ + "}")));
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
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractOptionsData
// Purpose       : Extract the parameters from a netlist .OPTIONS line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractOptionsData()
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
    if ( parsedLine[2].string_ == "=" )
    {
      setName( "PSPICE" );
      parameterStartPos = 1;

      string msg("PSPICE style options are not currently supported in Xyce\n");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
      return false;
    }
    else
    {
      string msg("Unrecognized .OPTIONS package " + optionName + "\n");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, parsedLine[1].lineNumber_);
    }
  }

  // Create an option block to temporarily store the default options.
  N_IO_OptionBlock defaultOptions("", parsedLine, metadata_ );

  // Get the default options from metadata.
  defaultOptions.addDefaultOptionsParameters( getName() );

  // Extract the parameters from parsedLine.
  int parameterEndPos = numFields - 1;
  vector<N_UTL_Param> inputParameters;
  N_UTL_Param parameter("", "");
  int intervalParameterStart = -1;
  int i = parameterStartPos;
  string paramBaseName;
  while (i <= parameterEndPos-1)
  {
    // Check for equal sign.
    if ( parsedLine[i+1].string_ != "=" )
    {
      if ( optionName != "OUTPUT" && optionName != "RESTART" )
      {
        string msg("Equal sign required between parameter name\n");
        msg += "and value in .OPTIONS " + getName() + "\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
            netlistFileName_, parsedLine[i].lineNumber_);
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
    N_UTL_Param* parameterPtr = defaultOptions.findParameter(parameter);
    if (parameterPtr == NULL)
    {
      // Options parameter not found, print out warning.
      string msg("No options parameter " + parameter.tag());
      msg += " found in metadata.\n";
      msg += "This parameter will be ignored.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
      i+= 3;
    }
    else if (parameterPtr->sVal() != "VECTOR")
    {
      parameter.setVal( parsedLine[i+2].string_ );
      inputParameters.push_back( parameter );
      i+= 3;

      if (i < parameterEndPos-1 && parsedLine[i].string_ == ",")
      {
        string msg("Options parameter " + parameter.tag());
        msg += " is flagged as not VECTOR, but has comma in value.\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
                                 netlistFileName_, parsedLine[0].lineNumber_);
      }
    }
    else
    {
      // We have a vector valued parameter.
      // Name the jth component of the parameter of the vector by appending
      // "j" to the parameter name.
      ostringstream paramName;
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
    N_UTL_Param* parameterPtr =
      defaultOptions.findParameter( inputParameters[k] );
    if ( parameterPtr != NULL )
    {
      parameterPtr->setVal( inputParameters[k] );
      addParameter( *parameterPtr );
    }
    else
    {
      // Options parameter not found, print out warning.
      string msg("No options parameter " + inputParameters[k].tag());
      msg += " found in metadata.\n";
      msg += "This parameter will be ignored.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
    }
  }

  // If this is a ".OPTIONS OUTPUT" line or ".OPTIONS RESTART" line
  // then get the time and interval pairs.
  if ( (optionName == "OUTPUT" || optionName == "RESTART") &&
        intervalParameterStart != -1)
  {
    N_UTL_Param time( "TIME", "" );
    N_UTL_Param interval( "INTERVAL", "" );

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
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::addDefaultOptionsParameters
// Purpose       : Add the default parameters for a model from metadata.
// Special Notes : While not required, it is generally expected that
//                 the parameters given in the input netlist will have
//                 already been extracted via extractModelData().
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/17/2001
//-----------------------------------------------------------------------------
void N_IO_OptionBlock::addDefaultOptionsParameters( string const& optionName )
{
  vector<N_UTL_Param> * optionsParameterPtr =
    metadata_.getPtrToOptionsParameters(optionName);

  addParameters(*optionsParameterPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::extractDCData
// Purpose       : Determine number of sweep variables on this DC line
//                : and create option blocks for each, storing the blocks
//                : in the referenced vector.
//-----------------------------------------------------------------------------
int N_IO_OptionBlock::extractDCData( vector< N_IO_OptionBlock > & oBs )
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
    oBs.push_back( N_IO_OptionBlock( netlistFileName_, metadata_ ) );
    oBs[sourcesFound].setName( "DC" );

    N_UTL_Param parameter("", "");
    ExtendedString stringVal ( parsedLine[linePosition + 1].string_ );
    stringVal.toUpper();

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
       ( N_UTL_Param( "", parsedLine[linePosition].string_ ) ).isNumeric() )
      {
        parameter.setTag( "VAL" );
        parameter.setVal( parsedLine[linePosition].string_ );
        oBs[sourcesFound].addParameter( parameter );
      }
    }

    else
    {
      string sweepStepTag ("STEP");

      // check for non-LIST default (LINear sweep)
      if( ( N_UTL_Param( "", stringVal ) ).isNumeric() )
      {
        stringVal = "LIN";
      }

      // non-LIST sweep type is given
      else
      {
        // get sweep type name and move to sweep variable name
        stringVal = parsedLine[linePosition++].string_;
        stringVal.toUpper();

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
        string msg(".DC line not formatted correctly, found unexpected number ");
        msg += "of fields\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
         netlistFileName_, parsedLine[0].lineNumber_);
      }

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

    // record this source (and move on to the next)
    ++sourcesFound;
  }

  return sourcesFound;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractDCData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractDCData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "DC" );

  // Check that the minimum required number of fields are on the line.
  if ( (numFields-1)%4 != 0 )
  {
    string msg(".DC line not formatted correctly, found unexpected number ");
    msg += "of fields\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  int linePosition = 1;   // Start of parameters on .param line.

  N_UTL_Param parameter("", "");
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

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractSTEPData
// Purpose       : Extract the parameters from a netlist .STEP line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/30/2003
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractSTEPData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "STEP" );

  // First check if the type has been explicitly set.
  // If not, set it to the default, LIN.
  int pos1=1;

  bool typeExplicitSetLinDecOct = false;
  bool typeExplicitSetList = false;
  string type("LIN");
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
      string msg(".STEP line not formatted correctly.\n");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
    }
  }

  int linePosition = 1;   // Start of parameters on .param line.
  N_UTL_Param parameter("", "");

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
    string msg(".STEP line contains an unrecognized type");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractMPDEData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractMPDEData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "MPDE" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 || numFields > 6 )
  {
    string msg(".MPDE line has an unexpected number of fields\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields - 1;

  N_UTL_Param parameter("", "");

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
    string msg("expected NOOP/UIC field on .MPDE line but found");
    msg += parameter.usVal() + "\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[endPosition].lineNumber_);
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
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractHBData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractHBData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "HB" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 1 || numFields > 2 )
  {
    string msg(".HB line has an unexpected number of fields\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields - 1;

  N_UTL_Param parameter("", "");

  // frequency of oscillation is required
  parameter.setTag( "FREQ" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractTRANData
// Purpose       : Extract the parameters from a netlist .DC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractTRANData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "TRAN" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 || numFields > 6 )
  {
    string msg(".TRAN line has an unexpected number of fields\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields;

  N_UTL_Param parameter("", "");

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
          string msg("expected NOOP/UIC field on .TRAN line but found");
          msg += parameter.usVal() + "\n";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
              netlistFileName_, parsedLine[linePosition].lineNumber_);
        }
      }
    }
    linePosition++;
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractPrintData
// Purpose       : Extract the parameters from a netlist .PRINT line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 05/02/2002
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractPrintData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "PRINT" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    cout << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << endl;
  }
#endif

  // Set the TYPE and add it the parameters.
  ExtendedString tmpString = parsedLine[1].string_;
  tmpString.toUpper();
  if (tmpString == "TR")
  {
    tmpString = "TRAN"; // TR is a synonym for TRAN
  }
  //N_UTL_Param typeParameter("TYPE", parsedLine[1].string_);
  N_UTL_Param typeParameter("TYPE", tmpString);
  if ( typeParameter.usVal() != "DC" &&
       typeParameter.usVal() != "TRAN" &&
       typeParameter.usVal() != "AC" &&
       typeParameter.usVal() != "HB" &&
       typeParameter.usVal() != "MOR" )
  {
    string msg("Invalid type \"" + typeParameter.usVal());
    msg += "\" for .PRINT statement\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[1].lineNumber_);
  }

  // Add the options parameter set to the model.
  addDefaultOptionsParameters( "PRINT" );

  // Reset the default TYPE with the value found.
  N_UTL_Param* parameterPtr = findParameter( typeParameter );

  if( parameterPtr == NULL )
  {
    string msg("Failed to find parameter\"" + typeParameter.usVal());
    msg += "\" in optionData for .PRINT statement\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
                             netlistFileName_, parsedLine[1].lineNumber_);
  }

  parameterPtr->setVal( typeParameter.usVal() );

  // Set no default value for the output file.
  // if the user set FILE=vale it will be picked up here.  Otherwise
  // the output manager will use the top level simulation netlist name
  // as the default.
  N_UTL_Param fileParameter("FILE", "");
  parameterPtr = findParameter( fileParameter );

  if( parameterPtr == NULL )
  {
    string msg("Failed to find parameter\"" + typeParameter.usVal());
    msg += "\" in optionData for .PRINT statement\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
                             netlistFileName_, parsedLine[1].lineNumber_);
  }

  parameterPtr->setVal( fileParameter.sVal() );

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
     N_UTL_Param* parameterPtr;
     while ( position+1 < parsedLine.size() &&
             parsedLine[position+1].string_ == "=" )
     {

       parameterPtr = findParameter( N_UTL_Param(parsedLine[position].string_, "") );
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
         // Options parameter not found, print out warning.
         string msg("No PRINT parameter " + parsedLine[position].string_);
         msg += " found in metadata.\n";
         msg += "This parameter will be ignored.\n";
         N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg,
             netlistFileName_, parsedLine[position].lineNumber_);
         return false;
       }

       position += 3;
     }
  }
#ifdef Xyce_DEBUG_IO
  // There is no point to this warning unless we're debugging
  else
  {
    string msg("No tagged parameters found\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }
#endif

  // Complete the PRINT line parsing.
  // Some of the remaining fields are of the
  // form I(Vname), V(node) or V(node1,node2), and these fields need
  // special treatment.
  N_UTL_Param parameter;
  ExtendedString field("");
  string msg("");
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
          msg = "Unrecognized current specification";
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
          parameter.setTag("V");
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
          msg = "Unrecognized parenthetical specification";
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
          msg = "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else
      {
        msg = "Unrecognized parenthetical specification";
        p_err = position;
      }
    }
    else // This is either an expression, or a STEP parameter.
    {
      if (parsedLine[position].string_ == "(" || parsedLine[position].string_ == ")")
      {
        msg = "Unrecognized parenthesis";
        p_err = position;
      }
      field = parsedLine[position].string_;
      field.toUpper();
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      addParameter ( parameter );
      ++position;
    }
    if (msg != "")
    {
      msg += " in .print near:\n";
      position = p_err+4;
      p_err -= 2;
      if (p_err<0)
        p_err = 0;
      if (position >= numFields)
        position = numFields;
      while (p_err < position)
      {
        msg += parsedLine[p_err].string_ + " ";
        ++p_err;
      }
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractRESULTData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 08/29/2004
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractRESULTData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "RESULT" );

#ifdef Xyce_DEBUG_IO
  // print out the parsed line
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    cout << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << endl;
  }
#endif

  int linePosition = 1;   // Start of parameters on .param line.
  N_UTL_Param parameter("", "");

  parameter.setTag( "EXPRESSION" );
  parameter.setVal( parsedLine[linePosition].string_ );
  addParameter( parameter );

  if (linePosition != (parsedLine.size () - 1))
  {
    string msg("Too many fields in the .RESULT line\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
    return false;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractACData
// Purpose       : Extract the parameters from a netlist .AC line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractACData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "AC" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 4 || numFields > 5)
  {
    string msg(".AC line has an unexpected number of fields\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  int linePosition = 1;   // Start of parameters on .param line.
  int endPosition = numFields - 1;

  N_UTL_Param parameter("", "");

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
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractMORData
// Purpose       : Extract the parameters from a netlist .MOR line held in
//                 parsedLine.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 5/25/12
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractMORData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "MOR" );

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 2 )
  {
    string msg(".MOR line has an unexpected number of fields\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
  }

  int linePosition = 1;   // Start of parameters on .param line.

  N_UTL_Param parameter("", "");

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
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractOBJETIVEData
// Purpose       : Extract the parameters from a netlist .OBJECTIVE line 
//
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractOBJECTIVEData()
{
  int numFields = parsedLine.size();
  int linePosition = 0;
  set<string> pars;

  // Set the OptionBlock name.
  setName( "OBJECTIVE" );

  N_UTL_Param parameter("", "");
  ExtendedString ES("");

  ++linePosition;
  while (linePosition < numFields-1)
  {
    ES = parsedLine[linePosition].string_ ;
    ES.toUpper ();
    if (pars.find(ES) != pars.end())
    {
      string msg("Multiple ");
      msg += ES;
      msg += " name specifications in .OBJECTIVE line";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
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
          string msg("Non-Expression specified in .OBJECTIVE argument: ");
          msg += ES;
          msg += " = ";
          msg += parsedLine[linePosition].string_;
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
              netlistFileName_, parsedLine[0].lineNumber_);
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
//    string::size_type first=parsedLine[linePosition].string_.find_first_of(':');
//    if (first != parsedLine[linePosition].string_.find_last_of(':') ||
//        first == string::npos)
//    {
//      string msg("Match specification '");
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
        string msg("Unrecognized field in .OBJECTIVE line");
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractMEASUREData
// Purpose       : Convert a .measure line to an options block
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractMEASUREData()
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
  //              v(a), v(a,b), v(a)=number v(a)=v(b).  So it really could be one or more
  //              variables and or a real number.  The comparison is always equity.
  //       options = these are keywords=value pairs that set limits on when or how
  //              the measurement is done. Value is usually a number, but can be a string in
  //              at least one case.
  //

  // we will use these sets to test for measure types and keywords.
  set<string> typeSet;
  typeSet.insert( string("TRIG") );
  typeSet.insert( string("TARG") );
  typeSet.insert( string("AVG") );
  typeSet.insert( string("MAX") );
  typeSet.insert( string("MIN") );
  typeSet.insert( string("PP") );
  typeSet.insert( string("RMS") );
  typeSet.insert( string("FREQ") );
  typeSet.insert( string("FIND") );
  typeSet.insert( string("WHEN") );
  typeSet.insert( string("PARAM") );
  typeSet.insert( string("DERIVATIVE") );
  typeSet.insert( string("DERIV") );
  typeSet.insert( string("DUTY") );
  typeSet.insert( string("INTEGRAL") );
  typeSet.insert( string("INTEG") );
  typeSet.insert( string("ERROR") );
  typeSet.insert( string("ON_TIME") );
  typeSet.insert( string("OFF_TIME") );
  typeSet.insert( string("FOUR") );

  set<string> keywords;
  set<string> simpleKeywords;
  set<string> numOrTextKeywords;

  simpleKeywords.insert( string("TD") );
  simpleKeywords.insert( string("GOAL") );
  simpleKeywords.insert( string("WEIGHT") );
  simpleKeywords.insert( string("MINVAL") );
  simpleKeywords.insert( string("AT") );
  simpleKeywords.insert( string("FROM") );
  simpleKeywords.insert( string("TO") );
  simpleKeywords.insert( string("IGNORE") );
  simpleKeywords.insert( string("YMIN") );
  simpleKeywords.insert( string("YMAX") );
  simpleKeywords.insert( string("ON") );
  simpleKeywords.insert( string("OFF") );
  simpleKeywords.insert( string("FRAC_MAX") );
  simpleKeywords.insert( string("MIN_THRESH") );
  simpleKeywords.insert( string("MAX_THRESH") );
  simpleKeywords.insert( string("NUMFREQ") );
  simpleKeywords.insert( string("GRIDSIZE") );

  numOrTextKeywords.insert( string("RISE") );
  numOrTextKeywords.insert( string("FALL") );
  numOrTextKeywords.insert( string("CROSS") );

  // make a union for the keywords set
  set_union( simpleKeywords.begin(), simpleKeywords.end(), numOrTextKeywords.begin(),
             numOrTextKeywords.end(), inserter<set<string> >(keywords, keywords.begin()) );

  int numFields = parsedLine.size();
  if( numFields < 4 )
  {
    string msg("Too few items on .MEASURE line.  Need at lest .MEASURE <mode> <name> <type>");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
         parsedLine[0].lineNumber_);
  }

  // Set the OptionBlock name.
  setName( "MEASURE" );

  N_UTL_Param parameter;

  // look for the fixed items
  ExtendedString currentWord(parsedLine[1].string_);
  currentWord.toUpper();
  if (currentWord == "TR")
  {
    currentWord = "TRAN"; // TR is a synonym for TRAN
  }

  if( currentWord == "DC" || currentWord == "TRAN")
  {
    parameter.set("MODE", currentWord);
    addParameter(parameter);
  }
  else
  {
    string msg("Unknown mode in .MEASURE line.  Should be DC or TRAN/TR");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
        parsedLine[1].lineNumber_);
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
    string msg("Illegal name in .MEASURE line.  Cannot be AC, DC or TRAN/TR ");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
        parsedLine[2].lineNumber_);
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
    string msg("Illegal type in .MEASURE line.  Must be one of: TRIG, TARG, AVG, ");
    msg += "MAX, MIN, PP, RMS, INTEG, FIND, WHEN, PARAM, DERIVATIVE, DERIV, ";
    msg += "INTEGRAL, ERROR, FOUR";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
        parsedLine[3].lineNumber_);
  }

  // already got MEASURE <DC|TRAN> name TYPE
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
        string & value = parsedLine[valPosition].string_;
        if( N_UTL::isInt(value) )
        {
          parameter.set(currentWord, N_UTL::Ival(value) );
        }
        else if( N_UTL::isValue(value) )
        {
          int valAsInt = static_cast<int>(N_UTL::Value(value));
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
        string msg("Incomplete .MEASURE line.  RISE, FALL, CROSS must be ");
        msg += "followed by a value or LAST";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
            parsedLine[position].lineNumber_);
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
        string & value = parsedLine[valPosition].string_;
        if( N_UTL::isValue(value) )
        {
          parameter.set(currentWord, N_UTL::Value(value) );
        }
        else if( N_UTL::isInt(value) )
        {
          parameter.set(currentWord, N_UTL::Ival(value) );
        }
        else
        {
          string msg("Incomplete .MEASURE line.  TD, GOAL, WEIGHT, MINVAL, AT, ");
          msg += "TO, IGNORE, YMIN, YMAX must be followed by a value";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
              parsedLine[position].lineNumber_);
        }
        addParameter(parameter);
      }
      else
      {
        string msg("Incomplete .MEASURE line.  TD, GOAL, WEIGHT, MINVAL, AT, TO, ");
        msg += "IGNORE, YMIN, YMAX must be followed by a value";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
            parsedLine[position].lineNumber_);
      }
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

      while( position < endPosition )
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
            string msg("Error in .MEASURE line.  Could not parse voltage/current variable");
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
                parsedLine[position].lineNumber_);
          }

        }
        else if( nextWord.isValue() )
        {
          string objValue("OBJVAL");
          parameter.set( objValue, nextWord.Value());
          addParameter( parameter );
        }
        else if( nextWord.isInt() )
        {
          string objValue( "OBJVAL");
          parameter.set( objValue, nextWord.Ival());
          addParameter( parameter );
        }
        position++;
      }
      // reset the position indicator to the end - 1 because
      // we're going to increment it at the end of the while loop
      position = endPosition - 1;

    }
    position++;
  }
#ifdef Xyce_DEBUG_IO
  N_IO_OptionBlock::printDiagnostic();
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractFOURIERData
// Purpose       : Convert a .four line to an options block
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 06/03/2013
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractFOURIERData()
{
  // Set the OptionBlock name
  setName( "FOUR" );

  N_UTL_Param parameter;
  ExtendedString nextWord("");
  
  if(parsedLine.size() < 3)
  {
    string msg("Error: the .FOUR line requires at least 3 arguments '.FOUR freq ov1 <ov2 ov3 ...>'");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg, netlistFileName_,
          parsedLine[0].lineNumber_);
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
        string msg("Error in .FOUR line.  Could not parse variable");
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
              parsedLine[position].lineNumber_);
      }
    }
    else
    {
      string msg = "Error in .FOUR line.  Could not parse variable: " + nextWord;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg, netlistFileName_,
            parsedLine[position].lineNumber_);
    }
    position+=2;

  }

return true;

}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractICData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/07/2007
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractICData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "IC" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    cout << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << endl;
  }
#endif

  string msg("");
  N_UTL_Param parameter;
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
        string & posString = parsedLine[position].string_;

        if (posString == "I" || posString == "i"
         || (posString.size() == 2 &&
             (posString[0] == 'I' || posString[0] == 'i')))
        {
          {
            msg = "Unsupported current specification.";
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
            msg = "Voltage differences not supported.";
            p_err = position;
          }
          else
          {
            msg = "Unrecognized parenthetical specification.";
            p_err = position;
          }
        }
        else
        {
          msg = "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else // This is either an expression, or a STEP parameter.
      {
        msg = "Unrecognized .IC and/or .DCVOLT specification";
        p_err = position;
      }

      if (msg != "")
      {
        msg += " in .IC and/or .DCVOLT near:\n";
        position = p_err+4;
        p_err -= 2;
        if (p_err<0)
          p_err = 0;
        if (position >= numFields)
          position = numFields;
        while (p_err < position)
        {
          msg += parsedLine[p_err].string_ + " ";
          ++p_err;
        }
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
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
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractNodeSetData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/07/2007
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractNodeSetData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "NODESET" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    cout << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << endl;
  }
#endif

  string msg("");

  N_UTL_Param parameter;
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
        string & posString = parsedLine[position].string_;

        if (posString == "I" || posString == "i"
         || (posString.size() == 2 &&
             (posString[0] == 'I' || posString[0] == 'i')))
        {
          {
            msg = "Unsupported current specification.";
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
            msg = "Voltage differences not supported.";
            p_err = position;
          }
          else
          {
            msg = "Unrecognized parenthetical specification.";
            p_err = position;
          }
        }
        else
        {
          msg = "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else // This is either an expression, or a STEP parameter.
      {
        msg = "Unrecognized .NODESET specification";
        p_err = position;
      }

      if (msg != "")
      {
        msg += " in .NODESET near:\n";
        position = p_err+4;
        p_err -= 2;
        if (p_err<0)
          p_err = 0;
        if (position >= numFields)
          position = numFields;
        while (p_err < position)
        {
          msg += parsedLine[p_err].string_ + " ";
          ++p_err;
        }
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
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
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractSaveData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractSaveData()
{
  int numFields = parsedLine.size();
  int linePosition = 0;

  // Set the OptionBlock name.
  setName( "SAVE" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    cout << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << endl;
  }
#endif

  string msg("");

  N_UTL_Param parameter("", "");
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
      string msg("Unrecognized field in .SAVE line");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg,
        netlistFileName_, parsedLine[0].lineNumber_);
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);


  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractLoadData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractLoadData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "LOAD" );

#ifdef Xyce_DEBUG_IO
  for (int ieric=0;ieric<parsedLine.size();++ieric)
  {
    cout << "parsedLine["<<ieric<<"] = " << parsedLine[ieric].string_ << endl;
  }
#endif

  string msg("");


  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::extractSENSData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool N_IO_OptionBlock::extractSENSData()
{
  int numFields = parsedLine.size();

  // Set the OptionBlock name.
  setName( "SENS" );

  int parameterStartPos = 1;

  // Create an option block to temporarily store the default options.
  N_IO_OptionBlock defaultOptions("", parsedLine, metadata_ );

  // Get the default options from metadata.
  defaultOptions.addDefaultOptionsParameters( getName() );

  // Extract the parameters from parsedLine.
  int parameterEndPos = numFields - 1;
  vector<N_UTL_Param> inputParameters;
  N_UTL_Param parameter("", "");
  int intervalParameterStart = -1;
  int i = parameterStartPos;
  string paramBaseName;
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
    N_UTL_Param* parameterPtr = defaultOptions.findParameter(parameter);
    if (parameterPtr == NULL)
    {
      // Options parameter not found, print out warning.
      string msg("No options parameter " + parameter.tag());
      msg += " found in metadata.\n";
      msg += "This parameter will be ignored.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
      i+= 3;
    }
    else if (parameterPtr->sVal() != "VECTOR")
    {
      parameter.setVal( parsedLine[i+2].string_ );
      inputParameters.push_back( parameter );
      i+= 3;

      if (i < parameterEndPos-1 && parsedLine[i].string_ == ",")
      {
        string msg("Options parameter " + parameter.tag());
        msg += " is flagged as not VECTOR, but has comma in value.\n";
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg,
                                 netlistFileName_, parsedLine[0].lineNumber_);
      }
    }
    else
    {
      // We have a vector valued parameter.
      // Name the jth component of the parameter of the vector by appending
      // "j" to the parameter name.
      ostringstream paramName;
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
    N_UTL_Param* parameterPtr =
      defaultOptions.findParameter( inputParameters[k] );
    if ( parameterPtr != NULL )
    {
      parameterPtr->setVal( inputParameters[k] );
      addParameter( *parameterPtr );
    }
    else
    {
      // Options parameter not found, print out warning.
      string msg("No options parameter " + inputParameters[k].tag());
      msg += " found in metadata.\n";
      msg += "This parameter will be ignored.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING_0, msg,
          netlistFileName_, parsedLine[0].lineNumber_);
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  parsedLine.clear();
  vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    (parsedLine).swap(parsedLine);

  return true; // Only get here on success.
}


//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::setName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void N_IO_OptionBlock::setName( string const& nameIn )
{
  optionData.setName(nameIn);
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::addParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void N_IO_OptionBlock::addParameter( N_UTL_Param const& parameter )
{
  optionData.getParams().push_back( parameter );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::addParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void N_IO_OptionBlock::addParameters( vector<N_UTL_Param> const& parametersIn )
{
  optionData.getParams().insert( optionData.getParams().end(), parametersIn.begin(), parametersIn.end() );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::clearParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void N_IO_OptionBlock::clearParameters()
{
  optionData.getParams().clear();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::setParamter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
void N_IO_OptionBlock::setParameter( int const& i, N_UTL_Param const& parameter )
{
  if ( i < getNumberOfParameters() )
  {
    list<N_UTL_Param>::iterator paramIter;
    paramIter = optionData.getParams().begin();
    for ( int j = 0; j < i; ++j )
    {
      ++paramIter;
    }

    paramIter->setVal( parameter.sVal() );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
const string& N_IO_OptionBlock::getName() const
{
  return optionData.getName();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::findParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
N_UTL_Param* N_IO_OptionBlock::findParameter( N_UTL_Param const& parameter )
{
  list<N_UTL_Param>::iterator paramIter = find( optionData.getParams().begin(), optionData.getParams().end(), parameter );
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
// Function      : N_IO_OptionBlock::findParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
N_UTL_Param* N_IO_OptionBlock::findParameter( string const& parameterName )
{
  N_UTL_Param parameter( parameterName, "");

  return findParameter( parameter );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_OptionBlock::getNumberOfParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
int N_IO_OptionBlock::getNumberOfParameters() const
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
N_UTL_Param N_IO_OptionBlock::getParameter( int const& i ) const
{
  if ( i < getNumberOfParameters() )
  {
    list<N_UTL_Param>::const_iterator paramIter;
    paramIter = optionData.getParams().begin();
    for ( int j = 0; j < i; ++j )
    {
      ++paramIter;
    }

    return *paramIter;
  }
  else
  {
    return N_UTL_Param( "", "" );
  }
}
