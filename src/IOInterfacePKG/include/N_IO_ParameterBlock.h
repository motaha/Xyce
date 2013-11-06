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
// Revision Number: $Revision: 1.32.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
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

// ---------- Forward Declarations ----------

//-----------------------------------------------------------------------------
// Class         : N_IO_ParameterBlock
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------

class N_IO_ParameterBlock : public Packable
{
  public:
    // Constructor.
    N_IO_ParameterBlock(){};

    // Constructor.
    N_IO_ParameterBlock(
        string const& fileName,
        vector<N_IO_SpiceSeparatedFieldTool::StringToken>
        const& parsedInputLine);

    // Copy Constructor.
    N_IO_ParameterBlock(N_IO_ParameterBlock const& rhsPB)
    : netlistFileName_(rhsPB.netlistFileName_),
      expressionValuedParams_(rhsPB.expressionValuedParams_),
      parsedLine(rhsPB.parsedLine),
      modelData(rhsPB.modelData)
      { };

    // Destructor.
    ~N_IO_ParameterBlock(){};

    // Public data.
    vector<N_IO_SpiceSeparatedFieldTool::StringToken> parsedLine;
    N_DEV_ModelBlock modelData;

    map < string, vector <vector<N_DEV_Param> > > inputCompositeParamVecMap;

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

    void setName( string const& name );
    void setType( string const& type );
    void setLevel( string const& level );
    void setLineNumber( string & netlistFile, int lineNumber );
    void addParameter( N_DEV_Param const& parameter );
    void addParameters( vector<N_DEV_Param> const& parameters );
    void setParameter( int const& i, N_DEV_Param const& parameter );

    const string& getName() const;
    const string& getType() const;
    int getLevel() const;
    vector<N_DEV_Param> & getParams();

    N_DEV_Param* findParameter( N_DEV_Param const& parameter );
    int getNumberOfParameters() const;
    N_DEV_Param getParameter( int const& i ) const;

    N_IO_ParameterBlock & operator=(N_IO_ParameterBlock const& rhsPB)
    {
      netlistFileName_ = rhsPB.netlistFileName_;
      expressionValuedParams_ = rhsPB.expressionValuedParams_;
      parsedLine = rhsPB.parsedLine;
      modelData = rhsPB.modelData;
      return *this;
    };

    bool operator==(string const& name);
    bool operator!=(string const& name);

    // Packing functionality.
    Packable * instance() const;

    // Counts bytes needed to pack block.
    int packedByteCount() const;

    // Packs OptionBlock into char buffer using MPI_PACK.
    void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;

    // Unpacks OptionBlock from char buffer using MPI_UNPACK.
    void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);


  private:
    string netlistFileName_;
    int lineNumber_;
    bool defaultApplied_;
    vector<N_DEV_Param> expressionValuedParams_;

    void addDefaultCompositeModelParameters_
       ( N_IO_CircuitMetadata & metadata,
         N_DEV_Param & baseParam ,
        map<string,bool> & paramMetadataExistMap);

};

inline bool N_IO_ParameterBlock::operator==(string const& rhsName)
{
  return (getName() == rhsName);
}

inline bool N_IO_ParameterBlock::operator!=(string const& rhsName)
{
  return (getName() != rhsName);
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::setName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void N_IO_ParameterBlock::setName( string const& nameIn )
{
  modelData.name = nameIn;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::setType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void N_IO_ParameterBlock::setType( string const& typeIn )
{
  modelData.type = typeIn;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::setLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void N_IO_ParameterBlock::setLevel( string const& levelIn )
{
  N_DEV_Param levelParam( "LEVEL", levelIn);
  modelData.level = levelParam.iVal();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::setLineNumber
// Purpose       : Pass netlist file name and line number to Parameter Block
//               : for error reporting
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/02/2006
//-----------------------------------------------------------------------------
inline void N_IO_ParameterBlock::setLineNumber( string & netlistFile,
                                        int lineNumber )
{
  modelData.netlistFileName_ = netlistFile;
  modelData.lineNumber_ = lineNumber;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline const string& N_IO_ParameterBlock::getName() const
{
  return modelData.name;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::getType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline const string& N_IO_ParameterBlock::getType() const
{
  return modelData.type;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::getLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline int N_IO_ParameterBlock::getLevel() const
{
  N_DEV_Param levelParam( "LEVEL", modelData.level );
  return levelParam.iVal();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::getParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/26/05
//-----------------------------------------------------------------------------
inline vector<N_DEV_Param> & N_IO_ParameterBlock::getParams()
{
  return modelData.params;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::addParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void N_IO_ParameterBlock::addParameter( N_DEV_Param const& parameter )
{
  modelData.params.push_back( parameter );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::addParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void N_IO_ParameterBlock::addParameters(
    vector<N_DEV_Param> const& parametersIn )
{
  modelData.params.insert( modelData.params.end(),
                           parametersIn.begin(), parametersIn.end() );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::setParamter
// Purpose       :
// Special Notes : It is assumed that getNumberOfParemeters was called to
//                 ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void N_IO_ParameterBlock::setParameter( int const& i,
                                          N_DEV_Param const& parameter )
{
  modelData.params[i].setVal(parameter.sVal());
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::getNumberOfParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline int N_IO_ParameterBlock::getNumberOfParameters() const
{
  return modelData.params.size();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_ParameterBlock::getParameter
// Purpose       :
// Special Notes : It is assumed getNumberOfParameters was called to ensure
//                 that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline N_DEV_Param N_IO_ParameterBlock::getParameter( int const& i ) const
{
  return modelData.params[i];
}


#endif
