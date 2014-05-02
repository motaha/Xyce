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
// Revision Number: $Revision: 1.76.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_CircuitBlock_h
#define Xyce_N_IO_CircuitBlock_h

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>

#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_UTL_OptionBlock.h>
#include <N_TOP_InsertionTool.h>

namespace Xyce {
namespace IO {

class CircuitBlockData;

typedef std::pair<std::ifstream*,SpiceSeparatedFieldTool*> FileSSFPair;

//-----------------------------------------------------------------------------
// Class         : CircuitBlock
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 09/06/2001
//-----------------------------------------------------------------------------

class CircuitBlock
{
public:
  // Constructor.
  CircuitBlock(
     const std::string & netlistFileName,
     CmdParse & cp,
     CircuitMetadata & md,
     std::map<std::string,int> & mn,
     std::map<std::string,FileSSFPair> & ssfm,
     CircuitContext & cc,
     OutputMgr * outputMgrPtr,
     int & uc,
     bool & gPI,
     std::map<std::string, RCP<Device::InstanceBlock> > & dNames,
     std::set<std::string> & nNames,
     const std::vector< std::pair< std::string, std::string> > & externalNetlistParams
               );

  // Constructor.
  CircuitBlock(
     std::string const& fileName,
     std::vector<SpiceSeparatedFieldTool::StringToken>
     const& parsedInputLine,
     CmdParse & cp,
     CircuitMetadata & md,
     std::map<std::string,int> & mn,
     std::map<std::string,FileSSFPair> & ssfm,
     CircuitContext & cc,
     OutputMgr * outputMgrPtr,
     int & uc,
     bool & gPI,
     CircuitBlock * mainCircPtr,
     DistributionTool* dtPtr,
     Topo::InsertionTool* itPtr,
     Device::DeviceInterface* diPtr,
     std::map<std::string, RCP<Device::InstanceBlock> > & dNames,
     std::set<std::string> & nNames,
     const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
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
  ~CircuitBlock();

  // Member functions.
  const std::string& getName() const;

  // This function parses the netlist file and fills in the
  // details of the circuit. This is phase 1 of netlist parsing.
  // The devices cannot be completely handled in this phase.
  bool parseNetlistFilePass1();
  bool parseNetlistFilePass1(const std::string &libSelect, std::string libInside);

  // Perform the second pass over the netlist, this phase primarily
  // handles devices.
  bool parseNetlistFilePass2();

  // Perform special pass for mutual inductances
  bool parseMutualInductances();

  // Print the contents of CircuitBlock.
  void print();

  // Set data_->ssfPtr_ .
  void setSSFPtr( SpiceSeparatedFieldTool* ssfPtr );

  void setStartPosition();
  void setEndPosition();
  void setFilePosition(std::streampos const& position);
  void setLinePosition(int const& position);
  const std::streampos getStartPosition() const;
  const std::streampos getEndPosition() const;
  int getLineStartPosition() const;
  int getLineEndPosition() const;

  // Extract subcircuit data from parsedLine.
  bool extractSubcircuitData(std::string, int);

  // Instatiate all of the devices in the current (sub)circuit. This
  // method will be invoked during as a part of pass 2 operations. It
  // will be invoked to create instances of the devices in the main
  // circuit and in each subcircuit instance.
  bool instantiateDevices(std::string libSelect, std::string libInside);

  void fixupYDeviceNames();

#ifdef Xyce_PARALLEL_MPI
  int getDeviceNames (std::vector<std::string> &);
  void checkDeviceNames (const std::vector<std::string> &names);
#endif
  void addTableData( DeviceBlock & device );
  //I- device
  //- Add a device to the circuit.

  // Add a model to the circuit.
  void addModel(ParameterBlock & model, std::string const& modelPrefix);

  // Add a set of options corresponding to a .OPTIONS netlist line
  // to the circuit.
  void addOptions(const OptionBlock &options );

  // Search the subcircuitInstanceTable of the current circuit block for the
  // subcircuit of the given name. If it is not found, recursively
  // search each parent subcircuit. Return a pointer to the circuit
  // block if it is found, otherwise return NULL.
  CircuitBlock* findSubcircuit( std::string const& subcircuitName );

  void registerDistributionTool(DistributionTool* dtPtr);
  void registerInsertionTool(Xyce::Topo::InsertionTool* insertionToolPtr);
  void registerDeviceInterface(Device::DeviceInterface* devIntPtr);

  // Receive the circuit context (from the distribution tool).
  bool receiveCircuitContext(CircuitContext&  ccIn);

  // Change netlist file name in slave procs
  void setFileName ( std::string & );

  // Process a device line on processor zero, or serial.
  bool handleDeviceLine(
     std::vector<SpiceSeparatedFieldTool::StringToken> const& deviceLine,
     const std::string &libSelect="", const std::string &libInside="");

  // Member Data.
  std::string netlistFileName;

  CircuitBlock* parentCircuitPtr;  // For subcircuits, points to the
  // circuitBlock instance that contains
  // this subcircuit. NULL for top level.

  CircuitMetadata & metadata;

  // Circuit Context object to hold context information.
  CircuitContext & circuitContext;
  inline CircuitContext *getCircuitContextPtr() {return &circuitContext;}

  OutputMgr * outputMgrPtr_;

  int & useCount;  // Counts copies so deletion of subcircuitInstanceTable can
  // be done properly in the class destructor.
  std::vector<DeviceBlock> mutualInductors_;
  std::map<std::string, CircuitBlock*> circuitBlockTable_;
  std::vector<DeviceBlock> subcircuitInstances;

  // The following containers need to be registered with Xyce
  std::map<std::string, Topo::NodeDevBlock *> deviceTable;
  std::map<std::string, Topo::NodeBlock *> nodeTable;
  std::map<std::string, Device::ModelBlock *> modelTable;
  std::list<Util::OptionBlock> optionsTable;
  Util::OptionBlock deviceOptions;

  // Lookup table for initial conditions
  std::map< std::string, std::vector< SpiceSeparatedFieldTool::StringToken > >
  initCondIndex;

  std::set<std::string> & nodeNames_;
  std::map<std::string, RCP<Device::InstanceBlock> > & deviceNames_;

  // keep track of K lines that need extracted
  std::multimap< CircuitContext *, DeviceBlock > rawMIs;


  // This is a map of node aliases.  Interface nodes to a subcircuit are removed as
  // the subcircuit is expanded into the netlist. We'll store the names of the interface
  // nodes in case the user accesses them elsewhere (as in a print statement)
  // The keys are the alias names and the values are the real cicuit node names
  std::map<std::string,std::string> aliasNodeMap_;
  std::set< std::string > aliasNodeMapHelper_;

  // this is the function that does the substitution of aliased nodes for real nodes
  bool substituteNodeAliases();

  // resolve expressions in optionBlocks like .print
  bool resolveExpressionsInOptionBlocks();

private: // Private attributes

  // Copy Constructor.
  CircuitBlock( CircuitBlock const& rhsCB );
  CircuitBlock& operator=(const CircuitBlock& rhsCB);

  CircuitBlockData * data_;

  // read a line of input from an istream
  void readline( std::istream & in, char * line );

  CmdParse & commandLine_;

  std::vector< std::pair< std::string, std::string> > externalNetlistParams_;

  // This is a counter variable that was previously a
  // static variable local to CircuitBlock::receiveXMLBuf.
  int numXMLBufReceived;

  // Checking the netlist syntax.
  bool cmdChecked;
  bool netlistSave;

  int devProcessedNumber;

  std::vector<std::string> nodeList_;         // External nodes of a subcircuit.
  std::vector<DeviceBlock> deviceList_;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CircuitBlock N_IO_CircuitBlock;

#endif // Xyce_N_IO_CircuitBlock_h
