//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_CircuitContext.h,v $
//
// Purpose        : Declare the circuit "context".
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 01/21/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revsion$
//
// Revsion Date   : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#ifndef N_IO_CIRCUITCONTEXT_H
#define N_IO_CIRCUITCONTEXT_H

// ---------- Standard Includes ----------
#include <list>
#include <map>
#include <set>
#include <vector>
#include <string>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Param.h>

#include <N_IO_fwd.h>
#include <N_IO_FunctionBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Param.h>
#include <N_UTL_NoCase.h>

#include <N_UTL_Packable.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace IO {

//----------------------------------------------------------------------------
// Class          : CircuitContext
// Purpose        :
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------

class CircuitContext : public Packable
{
public:
  typedef std::map<std::string, N_IO_ParameterBlock *, LessNoCase> ModelMap;

  struct MutualInductance : public Packable
  {
    std::map<std::string,double> inductors;
    std::string coupling;
    std::string model;
    std::string firstInductor;

    // set of inductor names and associated terminals
    std::map< std::string, std::vector< std::string > > terminals;

    // MI name
    std::string name;

    // sharedInductorTable and fastLookupTable key
    int sharedKey;

    MutualInductance() {}

    MutualInductance( N_IO_DeviceBlock & device );

    // Packing functionality.
    Packable * instance() const
    { return new MutualInductance; }

    // Counts bytes needed to pack block.
    int packedByteCount() const;

    // Packs into char buffer using MPI_PACK.
    void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;

    // Unpacks from char buffer using MPI_UNPACK.
    void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);

  }; // end of MutualIndutance struct

  // Constructors.
  CircuitContext(CircuitMetadata & md,
                 std::list<CircuitContext*> & cL,
                 CircuitContext*& ccPtr);

  // Destructor.
  ~CircuitContext();

  CircuitContext & operator=( const CircuitContext & rhs );

  // Begin a subcircuit to context.
  // Note: the current subcircuit context will remain active until
  // an endSubcircuitContext is issued.
  bool beginSubcircuitContext(
     std::string const& netlistFileName,
     std::vector<N_IO_SpiceSeparatedFieldTool::StringToken> & subcircuitLine);

  // End the current context, push it onto the previous contexts
  // list of contexts, and reset the currentContext_ pointer to the
  // previous context.
  void endSubcircuitContext();

  // Add a model to the current context. Note, the circuit context
  // will assume responsibility for model pointers.
  void addModel(N_IO_ParameterBlock * modelPtr);

  // Add a set of .PARAM parameters to the current context.
  void addParams(N_IO_OptionBlock const& param);

  // Do early resolution of quoted parameters
  void resolveQuote (Util::Param & parameter);

  // Add a set of .GLOBAL_PARAM parameters to the current context.
  void addGlobalParams(N_IO_OptionBlock const& param);

  // Add a global node
  void addGlobalNode (std::string &gnode);

  // Add a .FUNC function to the current context.
  void addFunction(N_IO_FunctionBlock const& function);

  // Setters and getters.
  void setName(std::string const& name);
  const std::string& getName() const;
  const std::string& getCurrentContextName() const;
  void setPrefix(std::string const& prefix);
  const std::string& getPrefix() const;
  std::map<std::string, std::string>* getNodeMapPtr() const;
  void setParentContextPtr( CircuitContext * const ptr );
  const CircuitContext * getParentContextPtr() const;
  CircuitContext * getCurrentContextPtr() const;

  // Get the node list for the current context.
  std::vector<std::string> const& getNodeList() const;

  // Increment the device count in the current context.
  void incrementDeviceCount();

  // Add an subcircuit name to the instance list. The instance list
  // will be needed to get a total device count.
  void addInstance(std::string const& subcircuitName,
                   std::string const& fileName, 
                   int const& lineNumber);

  // Resolve the parameters and functions in the curren context. Since this
  // operation is dependent on a given subcircuit instance, the resolution
  // for a given context may occur repeatedly.
  bool resolve(std::vector<Device::Param> const& subcircuitInstanceParams);

  // Set the current context to that corresponding to the given subcircuit
  // name. Save the previous context on the stack for later retrieval.
  // If there is no context corresponding to the subcircuit name, recursively
  // search parent contexts. Return true if found, false otherwise.
  bool setContext(std::string const& subcircuitName,
                  std::string const& subcircuitPrefix = "",
                  std::list<std::string> const& instanceNodes = std::list<std::string>(),
                  CircuitContext* previousContext = NULL) const;
  void setContext(CircuitContext* context) const;

  // Reset the context the context prior to the last invocation of
  // setContext.
  void restorePreviousContext() const;
  bool globalNode (const std::string &nodeName) const;

  // If the input paramter has an expression value, replace replace the
  // parameters and functions in the expression with their actual values.
  bool resolveParameter(Util::Param& parameter,
                        std::vector<std::string> exceptionStrings = std::vector<std::string>());

  // Determine if expressionString has any unresolved strings and
  // resolve appropriately. Return true if all strings are resolved
  // otherwise return false.
  bool resolveStrings(Util::Expression & expression,
                      std::vector<std::string> exceptionStrings = std::vector<std::string>());

  // Determine if expressionString has any unresolved functions and
  // resolve appropriately. Return true if all functions are resolved
  // otherwise return false.
  bool resolveFunctions(Util::Expression & expression);

  // Look for a parameter with tag parameterName in resolvedParams_.
  // Check current context and recursively check parent
  // contexts. Return the parameter if it is found, set the
  // parameter value to the empty string if it is not found.
  bool getResolvedParameter(Util::Param & parameter);

  // Look for a parameter with tag parameterName in resolvedGlobalParams_.
  // Check current context and recursively check parent
  // contexts. Return the parameter if it is found, set the
  // parameter value to the empty string if it is not found.
  bool getResolvedGlobalParameter(Util::Param & parameter);

  // Look for a function with tag functionName in resolvedFunctions_.
  // Check current context and recursively check parent
  // contexts. Return the function (as an Util::Param) if it is found,
  // set the Util::Param value to the empty string if it is not found.
  bool getResolvedFunction(Util::Param & parameter);

  void addMutualInductance( N_IO_DeviceBlock & device )
  {
    currentContextPtr_->mutualInductances_.push_back( MutualInductance( device ) );
  }

  std::vector<MutualInductance> & getMutualInductances()
  {
    return currentContextPtr_->mutualInductances_;
  }

  std::vector< std::set< std::string > > & getSharedInductorTable()
  {
    return currentContextPtr_->sharedInductorTable_;
  }

  std::set< std::string > & getAllCoupledInductors()
  {
    return currentContextPtr_->allCoupledInductors_;
  }

  std::vector< std::vector< int > > & getAllIndexedMIs()
  {
    return currentContextPtr_->allIndexedMIs_;
  }

  int getNumMILines()
  {
    return currentContextPtr_->kLines_.size();
  }

  // Convert all MIs into tokenized device lines
  void bundleMIs();

  // Retrieve one tokenized device line
  std::vector< N_IO_SpiceSeparatedFieldTool::StringToken > & getMILine( int i );

  bool haveMutualInductances()
  {
    return !( currentContextPtr_->mutualInductances_.empty() );
  }

  int totalMutualInductanceCount();

  // Search the models in the current context for the model of the
  // given name. If it is not found, recursively search each parent
  // context. Return a pointer to the parameterBlock for the model
  // if it is found, otherwise return NULL. Also, if the model is found,
  // construct the appropriate model prefix.
  bool findModel(std::string const& modelName,
                 N_IO_ParameterBlock* & modelPtr,
                 std::string& modelPrefix) const;

  bool findModel(std::string const& modelName, N_IO_ParameterBlock* & modelPtr) const;

  // Check whether a subcircuit context is dependent on subcircuit parameters.
  // These are parameters on the .subckt line identified by "params:"
  // keyword. The result should be true if either the current subcircuit
  // context or any subcircuit context in the hierarchy containing the
  // current subcircuit context has subcircuit parameters.
  bool hasSubcircuitParams();

  // Calculate the total number of devices starting at current context
  // and including all subcircuit instances.
  int getTotalDeviceCount();

  // Correct total number of devices after processing K-devices:
  void augmentTotalDeviceCount(int kLineCount, int coupledICount, int YDeviceCount);

  // Packing functionality.
  Packable * instance() const;

  // Counts bytes needed to pack block.
  int packedByteCount() const;

  // Packs OptionBlock into char buffer using MPI_PACK.
  void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;

  // Unpacks OptionBlock from char buffer using MPI_UNPACK.
  void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);

  std::map<std::string,int> devMap;
  ModelMap modMap;

  N_IO_OptionBlock *getGlobals() {return &resolvedGlobalParams_;}

private:
  CircuitContext();

  // Reference to a Pointer to the current context;
  // ERK.  This was once static data, but is now non-static, and ultimately
  // owned by the IO_NetlistImportTool class.  For that reason, it is a
  // reference to a pointer, as it still needs to act static.
  CircuitContext*& currentContextPtr_;

  mutable CircuitContext* parentContextPtr_;

  // Stack of contexts to track context changes.
  std::list<CircuitContext*> & contextList_;

  // Context name, corresponds to subcircuit names.
  std::string name_;

  int deviceCount_;
  std::list<std::string> instanceList_;
  std::map< std::string, NetlistLocation> instanceErrorInfo_; // This data can
  // be dropped at
  // serialization time.

  std::vector<std::string> nodeList_;
  std::vector<N_UTL_Param> subcircuitParameters_;

  std::map< std::string, CircuitContext* > circuitContextTable_;

  ModelMap models_;

  std::vector<N_UTL_Param> unresolvedParams_;
  std::set<std::string> globalNodes_;
  std::vector<N_UTL_Param> unresolvedGlobalParams_;
  std::vector<N_IO_FunctionBlock> unresolvedFunctions_;

  std::vector<MutualInductance> mutualInductances_;

  // lookup tables used to create semiPackedMIs_
  std::vector< std::set< std::string > > sharedInductorTable_;
  std::set< std::string > allCoupledInductors_;
  std::vector< std::vector< int > > allIndexedMIs_;

  // tokenized MIs
  std::vector< std::vector< N_IO_SpiceSeparatedFieldTool::StringToken > >kLines_;

  // Each of the following attributes is not set until pass 2, so they
  // do not need to be serialized.
  std::string subcircuitPrefix_;
  std::map<std::string, std::string> nodeMap_; // note: does not need to be serialized.
  bool resolved_;
  N_IO_OptionBlock resolvedParams_;
  N_IO_OptionBlock resolvedGlobalParams_;
  std::map<std::string, Util::Param> resolvedFunctions_;

  CircuitMetadata & metadata_;
};

//----------------------------------------------------------------------------
// Function       : CircuitContext::setName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline void CircuitContext::setName(std::string const& name)
{
  name_ = name;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline const std::string& CircuitContext::getName() const
{
  return name_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::setPrefix
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline void CircuitContext::setPrefix(std::string const& prefix)
{
  currentContextPtr_->subcircuitPrefix_ = prefix;
}


//----------------------------------------------------------------------------
// Function       : CircuitContext::getPrefix
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline const std::string& CircuitContext::getPrefix() const
{
  return currentContextPtr_->subcircuitPrefix_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getNodeMapPtr
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline std::map<std::string, std::string> * CircuitContext::getNodeMapPtr() const
{
  return &(currentContextPtr_->nodeMap_);
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getCurrentContextName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline const std::string& CircuitContext::getCurrentContextName() const
{
  return currentContextPtr_->name_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getNodeList
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
inline std::vector<std::string> const& CircuitContext::getNodeList() const
{
  return currentContextPtr_->nodeList_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::incrementDeviceCount
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
inline void CircuitContext::incrementDeviceCount()
{
  currentContextPtr_->deviceCount_++;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addInstance
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
inline void CircuitContext::addInstance(std::string const& subcircuitName, std::string const& fileName, int const& lineNumber)
{
  std::string subcircuitNameUpper(ExtendedString(subcircuitName).toUpper());

  currentContextPtr_->instanceList_.push_back(subcircuitNameUpper);
  currentContextPtr_->instanceErrorInfo_[subcircuitNameUpper] = NetlistLocation(fileName, lineNumber);
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::setParentContextPtr
// Purpose        : accessor
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
inline void CircuitContext::setParentContextPtr(
   CircuitContext * const ptr )
{
  parentContextPtr_ = ptr;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getParentContextPtr
// Purpose        : accessor
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
inline const CircuitContext * CircuitContext::getParentContextPtr()
  const
{
  return parentContextPtr_;
}

inline CircuitContext * CircuitContext::getCurrentContextPtr()
  const
{
  return currentContextPtr_;
}

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CircuitContext CircuitContext;

#endif
