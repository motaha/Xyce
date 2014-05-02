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
// Filename       : $RCSfile: N_IO_OptionBlock.h,v $
//
// Purpose        : Declare the N_IO_OptionBlock class an instantiation of
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
// Revision Number: $Revision: 1.39.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_IO_OPTIONBLOCK_H
#define N_IO_OPTIONBLOCK_H

#include <string>
#include <vector>
#include <iosfwd>

#include <N_UTL_OptionBlock.h>
#include <N_IO_fwd.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace IO {

class OptionBlock
{
public:
  // Constructor.
  OptionBlock(CircuitMetadata & md);

  // Constructor.
  OptionBlock(
     std::string const& fileName,
     std::vector<SpiceSeparatedFieldTool::StringToken>
     const& parsedInputLine,
     CircuitMetadata & md);

  // ctor
  OptionBlock( std::string const& fileName,
               CircuitMetadata & md );

  // Copy Constructor.
  OptionBlock(OptionBlock const& rhsOB);

  // Destructor.
  ~OptionBlock() {};

  // Assignment operator
  OptionBlock & operator = (const OptionBlock & right);

  // Public data.
  std::vector<SpiceSeparatedFieldTool::StringToken> parsedLine;
  Util::OptionBlock optionData;

  // Public methods.

  // Determine the line type and extract the data appropriately.
  bool extractData();

  // Extract the parameters from a netlist .OPTIONS line held in parsedLine.
  bool extractOPData();

  // Extract the parameters from a netlist .DCOP line held in parsedLine.
  bool extractDCOPData();

  // Extract the parameters from a netlist .OUTPUT line held in parsedLine.
  bool extractOutputData();

  // Extract the parameters from a netlist .PARAM line held in parsedLine.
  bool extractParamData();

  // Extract the parameters from a netlist .OPTIONS line held in parsedLine.
  bool extractOptionsData();

  // Extract the parameters from a netlist .DC line held in parsedLine.
  bool extractDCData();

  int extractDCData( std::vector< OptionBlock > & oBs );

  // Extract the parameters from a netlist .STEP line held in parsedLine.
  bool extractSTEPData();

  // Extract the parameters from a netlist .RESULT line held in parsedLine.
  bool extractRESULTData();

  // Extract the parameters from a netlist .OBJECTIVE line held in parsedLine.
  bool extractOBJECTIVEData();

  // Extract the parameters from a netlist .TRAN line held in parsedLine.
  bool extractTRANData();

  // Extract the parameters from a netlist .MPDE line held in parsedLine.
  bool extractMPDEData();

  // Extract the parameters from a netlist .MEASURE line held in parsedLine.
  bool extractMEASUREData();

  // Extract the parameters from a netlist .FOURIER line held in parsedLine.
  bool extractFOURIERData();

  // Extract the parameters from a netlist .HB line held in parsedLine.
  bool extractHBData();

  // Extract the parameters from a netlist .AC line held in parsedLine.
  bool extractACData();

  // Extract the parameters from a netlist .MOR line held in parsedLine.
  bool extractMORData();

  // Extract the parameters from a netlist .PRINT line held in parsedLine.
  bool extractPrintData();

  // Extract the parameters from a netlist .IC line held in parsedLine.
  bool extractICData();

  // Extract the parameters from a netlist .NODESET line held in parsedLine.
  bool extractNodeSetData();

  // Extract the parameters from a netlist .SAVE line held in parsedLine.
  bool extractSaveData();

  // Extract the parameters from a netlist .LOAD line held in parsedLine.
  bool extractLoadData();

  // Extract the parameters from a netlist .SENS line held in parsedLine.
  bool extractSENSData();

  // Prints the details of an OptionBlock to standard out.
  void printDiagnostic() const;

  void setName( const std::string &name );
  const std::string &getName() const;

  void setNetlistLocation(const std::string &path, int line_number);
    
  void addParameter( Util::Param const& parameter );
  void addParameters( std::vector<Util::Param> const& parameters );
  void setParameter( int const& i, Util::Param const& parameter );
  void clearParameters();

  Util::Param* findParameter( Util::Param const& parameter );
  Util::Param* findParameter( std::string const& parameterName );
  int getNumberOfParameters() const;
  Util::Param getParameter( int const& i ) const;

private:
  OptionBlock();

  // The name of the netlist file containing the .options line (for error
  // reporting.
  std::string netlistFileName_;

  CircuitMetadata & metadata_;

  // Add the parameters for a .options line together with their default
  // values.
  void addDefaultOptionsParameters( std::string const& optionName );
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::OptionBlock N_IO_OptionBlock;

#endif
