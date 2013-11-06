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
// Revision Number: $Revision: 1.56.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_IO_DEVICEBLOCK_H
#define N_IO_DEVICEBLOCK_H

// ---------- Standard Includes ----------

#include <string>
#include <vector>
#include <map>

// ----------   Xyce Includes   ----------

#include <N_IO_fwd.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Param.h>

#include <N_TOP_NodeDevBlock.h>

class N_IO_CircuitContext;

// ---------- Forward Declarations   ----------

class N_IO_DeviceBlock
{
  public:
    // Constructor.
    N_IO_DeviceBlock( N_IO_CircuitContext & cc, N_IO_CircuitMetadata & md );

    // Constructor.
    N_IO_DeviceBlock(
        N_IO_CircuitContext & cc,
        N_IO_CircuitMetadata & md,
        string const& fileName,
        vector<N_IO_SpiceSeparatedFieldTool::StringToken> const& parsedInputLine);

    N_IO_DeviceBlock( N_IO_DeviceBlock const& rhsDB );
    // Copy Constructor.

    ~N_IO_DeviceBlock();
    // Destructor.


    // Public member data.


    // Public methods.

    // Print the details of a device to standard out.
    void print();

    // Clear the device, reset all attributes to their defaults.
    void clear();

    // Extract the device data from parsed line. Use device metadata
    // to determine device type, number of nodes, etc.
    //bool extractData(N_IO_CircuitContext* circuitContextPtr,
        //N_IO_CircuitMetadata* metadataPtr_);

    bool extractData();

    // Extract the subcircuit instance data given on a netlist 'X' line.
    bool extractSubcircuitInstanceData();

    // Setters and Getters
    void setParsedLine(
      vector<N_IO_SpiceSeparatedFieldTool::StringToken> const& deviceLine);

    void setFileName(string const& fileName){netlistFileName_ = fileName;};

    void setName( string const& name );
    void setNetlistType( char type );
    void setNetlistType( string type );
    void addNodeValue( string const& nodeValue );
    void setNodeValues( list<tagged_param> const& nodeValues );
    void setNodeValue( int const& i, string const& nodeValue );
    void setModelName( string const& modelName );
    void addInstanceParameter( N_DEV_Param const& parameter );
    void addInstanceParameters( vector<N_DEV_Param> const& parameters );
    void setInstanceParameter( int const& i, N_DEV_Param & parameter );
    void setLineNumber( string & netlistFile, int lineNumber );

    vector<N_IO_SpiceSeparatedFieldTool::StringToken> getParsedLine() const;
    const string& getName();
    const string getNetlistDeviceType() const;
    const string& getModelName();
    int getNumberOfNodes();
    string getNodeValue( int const& i );
    void getAllNodeNames( list< string > & nodeNames );
    const list<tagged_param> & getNodeValues();
    N_DEV_Param* findInstanceParameter( N_DEV_Param const& parameter );
    N_DEV_Param* findInstanceParameter( string const& parameter );
    int getNumberOfInstanceParameters();
    N_DEV_Param getInstanceParameter( int const& i );
    void getInstanceParameters(vector<N_DEV_Param>& parameters);
    bool isSubcircuitInstance() const;
    bool isExtracted() const { return extracted_; }

    N_TOP_NodeDevBlock * getDeviceData();

  private:
    N_IO_DeviceBlock();

    // Private data.
    string netlistFileName_;
    vector<N_IO_SpiceSeparatedFieldTool::StringToken> parsedLine_;

    string netlistType_;
    N_TOP_NodeDevBlock deviceData_;

    bool subcircuitInstance_;

    N_IO_CircuitContext & circuitContext_;
    N_IO_CircuitMetadata & metadata_;

    bool extracted_;

    // Private methods.

    void checkNode (string const& n);
    bool extractSourceData();

    bool extractMutualInductanceData();

    bool extractBasicDeviceData();

    bool extractBehavioralDeviceData();

    bool extractYDeviceData();

    bool extractMIDeviceData();

    bool extractSwitchDeviceData();

    bool extractNodes(int modelLevel, int modelNamePosition);

    bool extractModelName( string& modelType,
                           int & modelLevel,
                           int & modelNamePosition );

    void extractInstanceParameters( int searchStartPosition,
                                    int & parameterStartPosition,
                                    int & parameterEndPosition,
                                    string const& action,
                                    int modelLevel = -1 );

    void addDefaultInstanceParameters(int modelLevel);


    int findSourceFieldPosition( string const& fieldToFind,
                                 int startPosition );

    bool extractSourceFields( vector<string> const& fieldNames,
                              vector<int> const& fieldPositions );

    // Look for expression valued parameters in the parameter
    // list, evaluate expressions found and reset the parameter
    // value accordingly.
    bool setParameterValues();

    void issueUnrecognizedParameterError(string const& parameterName);

    bool isValidDeviceType(const string & deviceType);

};

//----------------------------------------------------------------------------
// Function       : N_IO_DeviceBlock::setParsedLine
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/06/2003
//----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setParsedLine(
    vector<N_IO_SpiceSeparatedFieldTool::StringToken> const& deviceLine)
{
  parsedLine_ = deviceLine;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::setName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/31/2001
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setName( string const& nameIn )
{
  deviceData_.getDevBlock().setName(nameIn);
  deviceData_.getNodeBlock().set_id( nameIn );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::setNetlistType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/31/2001
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setNetlistType( char typeIn )
{
  netlistType_ = typeIn;
}

inline void N_IO_DeviceBlock::setNetlistType( string typeIn )
{
  netlistType_ = typeIn;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::addNodeValue
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::addNodeValue( string const& nodeValue )
{
  deviceData_.getNodeBlock().addNode(tagged_param(nodeValue, 0));
  checkNode(nodeValue);

}

//----------------------------------------------------------------------------
// Function       : N_IO_DeviceBlock::setNodeValues
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/12/2003
//----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setNodeValues(
    list<tagged_param> const& nodeValues)
{
  list<tagged_param>::const_iterator n=nodeValues.begin() ;
  list<tagged_param>::const_iterator nEnd=nodeValues.end() ;
  for ( ; n != nEnd ; ++n)
  {
    checkNode(n->tag);
  }

  deviceData_.getNodeBlock().set_NodeList(nodeValues);
}

//----------------------------------------------------------------------------
// Function       : N_IO_DeviceBlock::checkNode
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::checkNode (string const& n)
{

/* DEBUG_ELR_PROFILE - FIXME remove this block; mostly done in pass 2

  N_IO_CircuitContext* myContext = circuitContext_.getCurrentContextPtr();

  if (!(myContext->devMap.empty()))
  {
    if (myContext->devMap.find(n) != myContext->devMap.end())
    {
      string msg("Device and node share the same name: " + n);
      if (myContext->getName() != "")
	msg += " in subcircuit: " + myContext->getName();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
    }
  }

  if ( !(myContext->modMap.empty()) )
  {
    if (myContext->modMap.find(n) != myContext->modMap.end())
    {
      string msg("Model and node share the same name: " + n);
      if (myContext->getName() != "")
	msg += " in subcircuit: " + myContext->getName();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, msg );
    }
  }

*/

}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::setModelName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setModelName( string const& modelNameIn )
{
  deviceData_.getDevBlock().setModelName(modelNameIn);
  if( modelNameIn != "" ) deviceData_.getDevBlock().modelFlag = true;
  else                    deviceData_.getDevBlock().modelFlag = false;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::addInstanceParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::addInstanceParameter(N_DEV_Param const& parameter)
{
  deviceData_.getDevBlock().params.push_back( parameter );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::addInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::addInstanceParameters(
                                    vector<N_DEV_Param> const& parameters )
{
  deviceData_.getDevBlock().params.insert(
      deviceData_.getDevBlock().params.end(),
      parameters.begin(), parameters.end());
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::setInstanceParamter
// Purpose       :
// Special Notes : It assumed that getNumberOfInstanceParameters was called
//                 to ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setInstanceParameter( int const& i,
                                        N_DEV_Param & parameter )
{
  deviceData_.getDevBlock().params[i].setVal(parameter);
  deviceData_.getDevBlock().params[i].setGiven( parameter.given() );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::setLineNumber
// Purpose       : Pass netlist file name and line number to Device Block
//               : for error reporting
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/02/2006
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setLineNumber( string & netlistFile,
                                        int lineNumber )
{
  deviceData_.getDevBlock().netlistFileName_ = netlistFile;
  deviceData_.getDevBlock().lineNumber_ = lineNumber;
}

//----------------------------------------------------------------------------
// Function       : N_IO_DeviceBlock::getParsedLine
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline vector<N_IO_SpiceSeparatedFieldTool::StringToken>
    N_IO_DeviceBlock::getParsedLine() const
{
  return parsedLine_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline const string& N_IO_DeviceBlock::getName()
{
  return deviceData_.getDevBlock().getName();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getNetlistDeviceType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline const string N_IO_DeviceBlock::getNetlistDeviceType() const
{
  return netlistType_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getModelName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline const string& N_IO_DeviceBlock::getModelName()
{
  return deviceData_.getDevBlock().getModelName();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getNumberOfNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline int N_IO_DeviceBlock::getNumberOfNodes()
{
  return (deviceData_.getNodeBlock().get_NodeList()).size();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getNodeValue
// Purpose       :
// Special Notes : It is assumed that getNumberOfNodes was called to ensure
//                 that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline string N_IO_DeviceBlock::getNodeValue( int const& i )
{
  list<tagged_param>::const_iterator paramIter;
  paramIter = (deviceData_.getNodeBlock().get_NodeList()).begin();
  for ( int j = 0; j < i; j++ )
  {
    paramIter++;
  }
  return paramIter->tag;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getAllNodeNames
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::getAllNodeNames( list< string > & nodeNames )
{
  list<tagged_param>::const_iterator paramIter, paramEnd;

  paramIter = ( deviceData_.getNodeBlock().get_NodeList() ).begin();
  paramEnd  = ( deviceData_.getNodeBlock().get_NodeList() ).end();

  nodeNames.clear();

  for ( ; paramIter != paramEnd; ++paramIter )
  {
    nodeNames.push_back( paramIter->tag );
  }
}


//----------------------------------------------------------------------------
// Function       : N_IO_DeviceBlock::getNodeValues
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/13/2003
//----------------------------------------------------------------------------
inline const list<tagged_param> & N_IO_DeviceBlock::getNodeValues()
{
  return deviceData_.getNodeBlock().get_NodeList();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::setNodeValue
// Purpose       :
// Special Notes : It is assumed getNumberOfNodes was called to ensure
//                 that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline void N_IO_DeviceBlock::setNodeValue(
    int const& i, string const& nodeValue )
{
  list<tagged_param>::const_iterator paramIter;
  paramIter = deviceData_.getNodeBlock().get_NodeList().begin();
  for ( int j = 0; j < i; j++ )
  {
    paramIter++;
  }
  tagged_param* paramPtr = const_cast<tagged_param*> (&(*paramIter));
  paramPtr->tag = nodeValue;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getNumberOfInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline int N_IO_DeviceBlock::getNumberOfInstanceParameters()
{
  return deviceData_.getDevBlock().params.size();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getInstanceParameter
// Purpose       :
// Special Notes : It is assumed getNumberOfInstanceParameters was called
//                 to ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline N_DEV_Param N_IO_DeviceBlock::getInstanceParameter( int const& i )
{
  return deviceData_.getDevBlock().params[i];
}

//-----------------------------------------------------------------------------
// Function      : N_IO_DeviceBlock::getDeviceData
// Purpose       :
// Special Notes : It is assumed getNumberOfInstanceParameters was called
//                 to ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
inline N_TOP_NodeDevBlock * N_IO_DeviceBlock::getDeviceData()
{
  return &deviceData_;
}

//----------------------------------------------------------------------------
// Function       : N_IO_DeviceBlock::isSubcircuitInstance
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/10/2003
//----------------------------------------------------------------------------
inline bool N_IO_DeviceBlock::isSubcircuitInstance() const
{
  return subcircuitInstance_;
}

#endif
