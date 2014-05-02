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
// Filename       : $RCSfile: N_IO_ParameterBlock.h,v $
//
// Purpose        : Declare the N_IO_ParameterBlock class instantiations of
//                  which are associated with netlist .MODEL lines.
//
// Special Notes  : ERK.  It seems that the name "N_IO_ModelBlock" would have been
//                  more appropriate and less confusing, or possibly
//                  N_IO_ModelParametersBlock.  Calling it "parameters block"
//                  makes it sound much more general than it really is.
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/10/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.41.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_IO_PARAMETERBLOCK_H
#define N_IO_PARAMETERBLOCK_H

// ---------- Standard Includes ----------

#include <string>

#include <vector>


// ----------   Xyce Includes   ----------

#include <N_IO_fwd.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

#include <N_DEV_Param.h>

#include <N_UTL_Packable.h>

#include <N_DEV_DeviceBlock.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : ParameterBlock
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------

class ParameterBlock : public Packable
{
public:
  // Constructor.
  ParameterBlock(){};

  // Constructor.
  ParameterBlock(
     std::string const& fileName,
     std::vector<N_IO_SpiceSeparatedFieldTool::StringToken>
     const& parsedInputLine);

  // Copy Constructor.
  ParameterBlock(ParameterBlock const& rhsPB)
    : netlistFileName_(rhsPB.netlistFileName_),
      expressionValuedParams_(rhsPB.expressionValuedParams_),
      parsedLine(rhsPB.parsedLine),
      modelData(rhsPB.modelData)
  { };

  // Destructor.
  ~ParameterBlock(){};

  // Public data.
  std::vector<N_IO_SpiceSeparatedFieldTool::StringToken> parsedLine;
  Device::ModelBlock modelData;

  std::map< std::string, std::vector<std::vector<Device::Param> > > inputCompositeParamVecMap;

  // Public methods.

  void print();
  // Prints the details of an ParameterBlock to standard out.

  bool extractModelData( N_IO_CircuitMetadata & metadata );
  // Extract model data from parsedLine using model metadata.

  void addDefaultModelParameters( N_IO_CircuitMetadata & metadata );
  // Add the default parameters for a model from metadata.

  bool hasExpressionValuedParams(){return !expressionValuedParams_.empty();};

  void setParameterValues(N_IO_CircuitContext* contextPtr);
  // Check each parameter to see if it is an expression, if
  // so, evaluate the parameter.

  void setName( std::string const& name );
  void setType( std::string const& type );
  void setLevel( std::string const& level );
  void setLineNumber( std::string & netlistFile, int lineNumber );
  void addParameter( Device::Param const& parameter );
  void addParameters( std::vector<Device::Param> const& parameters );
  void setParameter( int const& i, Device::Param const& parameter );

  NetlistLocation netlistLocation() const;
    
  const std::string& getName() const;
  const std::string& getType() const;
  int getLevel() const;
  std::vector<Device::Param> & getParams();

  Device::Param* findParameter( Device::Param const& parameter );
  int getNumberOfParameters() const;
  Device::Param getParameter( int const& i ) const;

  ParameterBlock & operator=(ParameterBlock const& rhsPB)
  {
    netlistFileName_ = rhsPB.netlistFileName_;
    expressionValuedParams_ = rhsPB.expressionValuedParams_;
    parsedLine = rhsPB.parsedLine;
    modelData = rhsPB.modelData;
    return *this;
  };

  bool operator==(std::string const& name);
  bool operator!=(std::string const& name);

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
  int lineNumber_;
  bool defaultApplied_;
  std::vector<Device::Param> expressionValuedParams_;

  void addDefaultCompositeModelParameters_
  ( N_IO_CircuitMetadata & metadata,
    Device::Param & baseParam ,
    std::map<std::string,bool> & paramMetadataExistMap);

};

inline bool ParameterBlock::operator==(std::string const& rhsName)
{
  return (getName() == rhsName);
}

inline bool ParameterBlock::operator!=(std::string const& rhsName)
{
  return (getName() != rhsName);
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setName( std::string const& nameIn )
{
  modelData.name = nameIn;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setType( std::string const& typeIn )
{
  modelData.type = typeIn;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setLevel( std::string const& levelIn )
{
  Device::Param levelParam( "LEVEL", levelIn);
  modelData.level = levelParam.getImmutableValue<int>();
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setLineNumber
// Purpose       : Pass netlist file name and line number to Parameter Block
//               : for error reporting
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/02/2006
//-----------------------------------------------------------------------------
inline void ParameterBlock::setLineNumber( std::string & netlistFile,
                                           int lineNumber )
{
  modelData.netlistFileName_ = netlistFile;
  modelData.lineNumber_ = lineNumber;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline const std::string& ParameterBlock::getName() const
{
  return modelData.name;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline const std::string& ParameterBlock::getType() const
{
  return modelData.type;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline int ParameterBlock::getLevel() const
{
  Device::Param levelParam( "LEVEL", modelData.level );
  return levelParam.getImmutableValue<int>();
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/26/05
//-----------------------------------------------------------------------------
inline std::vector<Device::Param> & ParameterBlock::getParams()
{
  return modelData.params;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::addParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::addParameter( Device::Param const& parameter )
{
  modelData.params.push_back( parameter );
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::addParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::addParameters(
   std::vector<Device::Param> const& parametersIn )
{
  modelData.params.insert( modelData.params.end(),
                           parametersIn.begin(), parametersIn.end() );
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setParamter
// Purpose       :
// Special Notes : It is assumed that getNumberOfParemeters was called to
//                 ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setParameter( int const& i,
                                          Device::Param const& parameter )
{
  modelData.params[i].setVal(parameter.stringValue());
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getNumberOfParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline int ParameterBlock::getNumberOfParameters() const
{
  return modelData.params.size();
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getParameter
// Purpose       :
// Special Notes : It is assumed getNumberOfParameters was called to ensure
//                 that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline Device::Param ParameterBlock::getParameter( int const& i ) const
{
  return modelData.params[i];
}

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::ParameterBlock N_IO_ParameterBlock;

#endif
