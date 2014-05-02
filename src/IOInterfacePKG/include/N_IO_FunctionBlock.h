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
// Filename       : $RCSfile: N_IO_FunctionBlock.h,v $
//
// Purpose        : Declare an N_IO_FunctionBlock instantiations of which are
//                  associated with netlist .FUNC lines.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 12/26/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.22.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_IO_FUNCTIONBLOCK_H
#define N_IO_FUNCTIONBLOCK_H

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_UTL_Packable.h> 

namespace Xyce {
namespace IO {

class FunctionBlock : public Packable
{
public:
  // Constructor.
  FunctionBlock(){};

  // Constructor.
  FunctionBlock(
     std::string const& fileName,
     std::vector<N_IO_SpiceSeparatedFieldTool::StringToken> 
     const& parsedInputLine);

  // Copy Constructor
  FunctionBlock(FunctionBlock const& rhsFB);

  // Destructor
  ~FunctionBlock(){};

  // Public data.
  std::vector<N_IO_SpiceSeparatedFieldTool::StringToken> parsedLine;

  std::string functionName;

  std::string functionNameAndArgs;

  std::vector<std::string> functionArgs;

  std::string functionBody;

  // Public methods.

  void print();
  // Prints the details of an ParameterBlock to standard out.

  bool extractData();
  //- Extract model data from parsedLine using model metadata.

    
  // Packing functionality.
  Packable * instance() const;
  
  // Counts bytes needed to pack block.
  int packedByteCount() const;
  
  // Packs OptionBlock into char buffer using MPI_PACK.
  void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;
  
  // Unpacks OptionBlock from char buffer using MPI_UNPACK.
  void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);

    
private:
  std::string netlistFileName_;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::FunctionBlock FunctionBlock;

#endif

