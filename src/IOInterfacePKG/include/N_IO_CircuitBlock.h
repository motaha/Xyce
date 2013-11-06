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
// Filename       : $RCSfile: N_IO_CircuitBlock.h,v $
//
// Purpose        : Declare the circuit level containers for holding netlist
//                  circuit data and the associated circuit level methods.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/06/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.65.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_IO_CIRCUITBLOCK_H
#define N_IO_CIRCUITBLOCK_H

// ---------- Standard Includes ----------

#include <list>
#include <map>
#include <string>

#include <iostream>
#include <fstream>

#include <vector>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_IO_fwd.h>
#include <N_IO_FunctionBlock.h>

#include <N_IO_OptionBlock.h>
#include <N_IO_CircuitContext.h>

#include <N_TOP_NodeDevBlock.h>
#include <N_DEV_DeviceBlock.h>
#include <N_TOP_NodeBlock.h>
#include <N_UTL_OptionBlock.h>
#include <N_IO_OutputMgr.h>

#include <N_TOP_InsertionTool.h>

#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------

class N_IO_CircuitBlockData;
class N_IO_DeviceBlock;
class N_IO_ParameterBlock;
class N_IO_SpiceSeparatedFieldTool;
class N_IO_DistributionTool;
class N_IO_CmdParse;

typedef pair<ifstream*,N_IO_SpiceSeparatedFieldTool*> FileSSFPair;

//-----------------------------------------------------------------------------
// Class         : N_IO_CircuitBlock
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 09/06/2001
//-----------------------------------------------------------------------------

class N_IO_CircuitBlock
{
  public:
    // Constructor.
    N_IO_CircuitBlock(
        string const& netlistFileName,
        N_IO_CmdParse & cp,
        N_IO_CircuitMetadata & md,
        map<string,int> & mn,
        map<string,FileSSFPair> & ssfm,
        N_IO_CircuitContext & cc,
        N_IO_OutputMgr * outputMgrPtr,
        int & uc,
        bool & gPI,
        map<string,RCP<N_DEV_InstanceBlock> > & dNames,
        set<string> & nNames,
        const vector< pair< string, string> > & externalNetlistParams
      );

    // Constructor.
    N_IO_CircuitBlock(
        string const& fileName,
        vector<N_IO_SpiceSeparatedFieldTool::StringToken>
        const& parsedInputLine,
        N_IO_CmdParse & cp,
        N_IO_CircuitMetadata & md,
        map<string,int> & mn,
        map<string,FileSSFPair> & ssfm,
        N_IO_CircuitContext & cc,
        N_IO_OutputMgr * outputMgrPtr,
        int & uc,
        bool & gPI,
        N_IO_CircuitBlock * mainCircPtr,
        N_IO_DistributionTool* dtPtr,
        Xyce::Topology::InsertionTool* itPtr,
        N_DEV_DeviceInterface* diPtr,
        map<string,RCP<N_DEV_InstanceBlock> > & dNames,
        set<string> & nNames,
        const vector< pair< string, string> > & externalNetlistParams,
        bool removeCvar,
        bool removeDvar,
        bool removeIvar,
        bool removeLvar,
        bool removeMvar,
        bool removeQvar,
        bool removeRvar,
        bool removeVvar,
        bool replgndvar
	);

    // Destructor.
    ~N_IO_CircuitBlock();

    // Member functions.
    const string& getName() const;

    // This function parses the netlist file and fills in the
    // details of the circuit. This is phase 1 of netlist parsing.
    // The devices cannot be completely handled in this phase.
    bool parseNetlistFilePass1();
    bool parseNetlistFilePass1(const string &libSelect, string libInside);

    // Perform the second pass over the netlist, this phase primarily
    // handles devices.
    bool parseNetlistFilePass2();

    // Perform special pass for mutual inductances
    bool parseMutualInductances();

    // Print the contents of CircuitBlock.
    void print();

    // Set data_->ssfPtr_ .
    void setSSFPtr( N_IO_SpiceSeparatedFieldTool* ssfPtr );

    void setStartPosition();
    void setEndPosition();
    void setFilePosition(streampos const& position);
    void setLinePosition(int const& position);
    const streampos getStartPosition() const;
    const streampos getEndPosition() const;
    int getLineStartPosition() const;
    int getLineEndPosition() const;

    // Extract subcircuit data from parsedLine.
    bool extractSubcircuitData(string, int);

    // Instatiate all of the devices in the current (sub)circuit. This
    // method will be invoked during as a part of pass 2 operations. It
    // will be invoked to create instances of the devices in the main
    // circuit and in each subcircuit instance.
    bool instantiateDevices(string libSelect, string libInside);

    void fixupYDeviceNames();

#ifdef Xyce_PARALLEL_MPI
    int getDeviceNames (vector<string> &);
    void checkDeviceNames (vector<string> &);
#endif
    void addTableData( N_IO_DeviceBlock & device );
    //I- device
    //- Add a device to the circuit.

    // Add a model to the circuit.
    void addModel(N_IO_ParameterBlock & model, string const& modelPrefix);

    // Add a set of options corresponding to a .OPTIONS netlist line
    // to the circuit.
    void addOptions( N_IO_OptionBlock const& options );

    // Search the subcircuitInstanceTable of the current circuit block for the
    // subcircuit of the given name. If it is not found, recursively
    // search each parent subcircuit. Return a pointer to the circuit
    // block if it is found, otherwise return NULL.
    N_IO_CircuitBlock* findSubcircuit( string const& subcircuitName );

    void registerDistributionTool(N_IO_DistributionTool* dtPtr);
    void registerInsertionTool(Xyce::Topology::InsertionTool* insertionToolPtr);
    void registerDeviceInterface(N_DEV_DeviceInterface* devIntPtr);

    // Receive the circuit context (from the distribution tool).
    bool receiveCircuitContext(N_IO_CircuitContext&  ccIn);

    // Change netlist file name in slave procs
    void setFileName ( string & );

    // Process a device line on processor zero, or serial.
    bool handleDeviceLine(
        vector<N_IO_SpiceSeparatedFieldTool::StringToken> const& deviceLine,
        const string &libSelect="", const string &libInside="");

    // Member Data.
    string netlistFileName;

    N_IO_CircuitBlock* parentCircuitPtr;  // For subcircuits, points to the
                                     // circuitBlock instance that contains
                                     // this subcircuit. NULL for top level.

    N_IO_CircuitMetadata & metadata;

    // Circuit Context object to hold context information.
    N_IO_CircuitContext & circuitContext;
    inline N_IO_CircuitContext *getCircuitContextPtr() {return &circuitContext;}

    N_IO_OutputMgr * outputMgrPtr_;

    int & useCount;  // Counts copies so deletion of subcircuitInstanceTable can
                     // be done properly in the class destructor.
    vector<N_IO_DeviceBlock> mutualInductors_;
    map<string, N_IO_CircuitBlock*> circuitBlockTable_;
    vector<N_IO_DeviceBlock> subcircuitInstances;

    // The following containers need to be registered with Xyce
    map<string, N_TOP_NodeDevBlock*> deviceTable;
    map<string, N_TOP_NodeBlock*> nodeTable;
    map<string, N_DEV_ModelBlock*> modelTable;
    list<N_UTL_OptionBlock> optionsTable;
    N_UTL_OptionBlock deviceOptions;

    // Lookup table for initial conditions
    map< string, vector< N_IO_SpiceSeparatedFieldTool::StringToken > >
     initCondIndex;

    set<string> & nodeNames_;
    map<string,Teuchos::RCP<N_DEV_InstanceBlock> > & deviceNames_;

    // keep track of K lines that need extracted
    multimap< N_IO_CircuitContext *, N_IO_DeviceBlock > rawMIs;


    // This is a map of node aliases.  Interface nodes to a subcircuit are removed as
    // the subcircuit is expanded into the netlist. We'll store the names of the interface
    // nodes in case the user accesses them elsewhere (as in a print statement)
    // The keys are the alias names and the values are the real cicuit node names
    map<string,string> aliasNodeMap_;
    set< string > aliasNodeMapHelper_;

    // this is the function that does the substitution of aliased nodes for real nodes
    bool substituteNodeAliases();

    // resolve expressions in optionBlocks like .print
    bool resolveExpressionsInOptionBlocks();

  private: // Private attributes

    // Copy Constructor.
    N_IO_CircuitBlock( N_IO_CircuitBlock const& rhsCB );
    N_IO_CircuitBlock& operator=(const N_IO_CircuitBlock& rhsCB);

    N_IO_CircuitBlockData* data_;

    // read a line of input from an istream
    void readline( istream & in, char * line );

    N_IO_CmdParse & commandLine_;

    vector< pair< string, string> > externalNetlistParams_;

    // This is a counter variable that was previously a
    // static variable local to N_IO_CircuitBlock::receiveXMLBuf.
    int numXMLBufReceived;

    // These are previously static variables
    // related to checking the netlist syntax.
    bool cmdChecked;
    bool netlistSave;

    // This is another formerly static variable.
#ifdef Xyce_DEBUG_IO
    int devProcessedNumber;
#endif

    vector<string> nodeList_;         // External nodes of a subcircuit.
    vector<N_IO_DeviceBlock> deviceList_;

};

#endif
