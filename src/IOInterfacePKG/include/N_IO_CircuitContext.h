//-------------------------------------------------------------------------
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
// Revsion Date   : $Date: 2013/10/03 17:23:42 $
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

#include <N_UTL_Packable.h>

// ---------- Forward Declarations ----------
class N_IO_ParameterBlock;
class N_IO_DeviceBlock;

//----------------------------------------------------------------------------
// Class          : N_IO_CircuitContext
// Purpose        :
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------

class N_IO_CircuitContext : public Packable
{
  public:

    struct MutualInductance : public Packable
    {
      map<string,double> inductors;
      string coupling;
      string model;
      string firstInductor;

      // set of inductor names and associated terminals
      map< string, vector< string > > terminals;

      // MI name
      string name;

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
    N_IO_CircuitContext(N_IO_CircuitMetadata & md,
                        list<N_IO_CircuitContext*> & cL,
                        N_IO_CircuitContext*& ccPtr);

    // Destructor.
    ~N_IO_CircuitContext();

    N_IO_CircuitContext & operator=( const N_IO_CircuitContext & rhs );

    // Begin a subcircuit to context.
    // Note: the current subcircuit context will remain active until
    // an endSubcircuitContext is issued.
    void beginSubcircuitContext(
        string const& netlistFileName,
        vector<N_IO_SpiceSeparatedFieldTool::StringToken> & subcircuitLine);

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
    void resolveQuote (N_UTL_Param & parameter);

    // Add a set of .GLOBAL_PARAM parameters to the current context.
    void addGlobalParams(N_IO_OptionBlock const& param);

    // Add a global node
    void addGlobalNode (string &gnode);

    // Add a .FUNC function to the current context.
    void addFunction(N_IO_FunctionBlock const& function);

    // Setters and getters.
    void setName(string const& name);
    const string& getName() const;
    const string& getCurrentContextName() const;
    void setPrefix(string const& prefix);
    const string& getPrefix() const;
    map<string, string>* getNodeMapPtr() const;
    void setParentContextPtr( N_IO_CircuitContext * const ptr );
    const N_IO_CircuitContext * getParentContextPtr() const;
    N_IO_CircuitContext * getCurrentContextPtr() const;

    // Get the node list for the current context.
    vector<string> const& getNodeList() const;

    // Increment the device count in the current context.
    void incrementDeviceCount();

    // Add an subcircuit name to the instance list. The instance list
    // will be needed to get a total device count.
    void addInstance(string const& subcircuitName,
                     int const& lineNumber,
                     string const& fileName);

    // Resolve the parameters and functions in the curren context. Since this
    // operation is dependent on a given subcircuit instance, the resolution
    // for a given context may occur repeatedly.
    bool resolve(vector<N_DEV_Param> const& subcircuitInstanceParams);

    // Set the current context to that corresponding to the given subcircuit
    // name. Save the previous context on the stack for later retrieval.
    // If there is no context corresponding to the subcircuit name, recursively
    // search parent contexts. Return true if found, false otherwise.
    bool setContext(string const& subcircuitName,
                    string const& subcircuitPrefix = "",
                    list<string> const& instanceNodes = list<string>(),
                    N_IO_CircuitContext* previousContext = NULL);
    void setContext(N_IO_CircuitContext* context);

    // Reset the context the context prior to the last invocation of
    // setContext.
    void restorePreviousContext();
    bool globalNode (const string &nodeName);

    // If the input paramter has an expression value, replace replace the
    // parameters and functions in the expression with their actual values.
    bool resolveParameter(N_UTL_Param& parameter,
                          vector<string> exceptionStrings = vector<string>());

    // Determine if expressionString has any unresolved strings and
    // resolve appropriately. Return true if all strings are resolved
    // otherwise return false.
    bool resolveStrings(N_UTL_Expression & expression,
                          vector<string> exceptionStrings = vector<string>());

    // Determine if expressionString has any unresolved functions and
    // resolve appropriately. Return true if all functions are resolved
    // otherwise return false.
    bool resolveFunctions(N_UTL_Expression & expression);

    // Look for a parameter with tag parameterName in resolvedParams_.
    // Check current context and recursively check parent
    // contexts. Return the parameter if it is found, set the
    // parameter value to the empty string if it is not found.
    bool getResolvedParameter(N_UTL_Param & parameter);

    // Look for a parameter with tag parameterName in resolvedGlobalParams_.
    // Check current context and recursively check parent
    // contexts. Return the parameter if it is found, set the
    // parameter value to the empty string if it is not found.
    bool getResolvedGlobalParameter(N_UTL_Param & parameter);

    // Look for a function with tag functionName in resolvedFunctions_.
    // Check current context and recursively check parent
    // contexts. Return the function (as an N_UTL_Param) if it is found,
    // set the N_UTL_Param value to the empty string if it is not found.
    bool getResolvedFunction(N_UTL_Param & parameter);

    void addMutualInductance( N_IO_DeviceBlock & device )
    {
      currentContextPtr_->mutualInductances_.push_back( MutualInductance( device ) );
    }

    vector<MutualInductance> & getMutualInductances()
    {
      return currentContextPtr_->mutualInductances_;
    }

    vector< set< string > > & getSharedInductorTable()
    {
      return currentContextPtr_->sharedInductorTable_;
    }

    set< string > & getAllCoupledInductors()
    {
      return currentContextPtr_->allCoupledInductors_;
    }

    vector< vector< int > > & getAllIndexedMIs()
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
    vector< N_IO_SpiceSeparatedFieldTool::StringToken > & getMILine( int i );

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
    bool findModel(string const& modelName,
        N_IO_ParameterBlock* & modelPtr,
        string& modelPrefix);

    bool findModel(string const& modelName, N_IO_ParameterBlock* & modelPtr);

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

    map<string,int> devMap;
    map<string,N_IO_ParameterBlock *> modMap;

    N_IO_OptionBlock *getGlobals() {return &resolvedGlobalParams_;}

  private:
    N_IO_CircuitContext();

    // Reference to a Pointer to the current context;
    // ERK.  This was once static data, but is now non-static, and ultimately
    // owned by the IO_NetlistImportTool class.  For that reason, it is a
    // reference to a pointer, as it still needs to act static.
    N_IO_CircuitContext*& currentContextPtr_;

    N_IO_CircuitContext* parentContextPtr_;

    // Stack of contexts to track context changes.
    list<N_IO_CircuitContext*> & contextList_;

    // Context name, corresponds to subcircuit names.
    string name_;

    int deviceCount_;
    list<string> instanceList_;
    map< string, pair<int, string> > instanceErrorInfo_; // This data can
                                                         // be dropped at
                                                         // serialization time.

    vector<string> nodeList_;
    vector<N_UTL_Param> subcircuitParameters_;

    map< string, N_IO_CircuitContext* > circuitContextTable_;

    map<string, N_IO_ParameterBlock*> models_;

    vector<N_UTL_Param> unresolvedParams_;
    set<string> globalNodes_;
    vector<N_UTL_Param> unresolvedGlobalParams_;
    vector<N_IO_FunctionBlock> unresolvedFunctions_;

    vector<MutualInductance> mutualInductances_;

    // lookup tables used to create semiPackedMIs_
    vector< set< string > > sharedInductorTable_;
    set< string > allCoupledInductors_;
    vector< vector< int > > allIndexedMIs_;

    // tokenized MIs
    vector< vector< N_IO_SpiceSeparatedFieldTool::StringToken > >kLines_;

    // Each of the following attributes is not set until pass 2, so they
    // do not need to be serialized.
    string subcircuitPrefix_;
    map<string, string> nodeMap_; // note: does not need to be serialized.
    bool resolved_;
    N_IO_OptionBlock resolvedParams_;
    N_IO_OptionBlock resolvedGlobalParams_;
    map<string, N_UTL_Param> resolvedFunctions_;

    N_IO_CircuitMetadata & metadata_;
};

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::setName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline void N_IO_CircuitContext::setName(string const& name)
{
  name_ = name;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline const string& N_IO_CircuitContext::getName() const
{
  return name_;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::setPrefix
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline void N_IO_CircuitContext::setPrefix(string const& prefix)
{
  currentContextPtr_->subcircuitPrefix_ = prefix;
}


//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getPrefix
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline const string& N_IO_CircuitContext::getPrefix() const
{
  return currentContextPtr_->subcircuitPrefix_;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getNodeMapPtr
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline map<string, string> * N_IO_CircuitContext::getNodeMapPtr() const
{
  return &(currentContextPtr_->nodeMap_);
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getCurrentContextName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline const string& N_IO_CircuitContext::getCurrentContextName() const
{
  return currentContextPtr_->name_;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getNodeList
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
inline vector<string> const& N_IO_CircuitContext::getNodeList() const
{
  return currentContextPtr_->nodeList_;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::incrementDeviceCount
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
inline void N_IO_CircuitContext::incrementDeviceCount()
{
  currentContextPtr_->deviceCount_++;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::addInstance
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
inline void N_IO_CircuitContext::addInstance(string const& subcircuitName,
    int const& lineNumber, string const& fileName)
{
  string subcircuitNameUpper(ExtendedString(subcircuitName).toUpper());

  currentContextPtr_->instanceList_.push_back(
      subcircuitNameUpper);

  currentContextPtr_->instanceErrorInfo_[subcircuitNameUpper] =
    pair<int, string>(lineNumber, fileName);
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::setParentContextPtr
// Purpose        : accessor
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
inline void N_IO_CircuitContext::setParentContextPtr(
 N_IO_CircuitContext * const ptr )
{
  parentContextPtr_ = ptr;
}

//----------------------------------------------------------------------------
// Function       : N_IO_CircuitContext::getParentContextPtr
// Purpose        : accessor
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
inline const N_IO_CircuitContext * N_IO_CircuitContext::getParentContextPtr()
 const
{
  return parentContextPtr_;
}

inline N_IO_CircuitContext * N_IO_CircuitContext::getCurrentContextPtr()
 const
{
  return currentContextPtr_;
}


#endif
