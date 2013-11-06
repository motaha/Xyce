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

//-------------------------------------------------------------------------
// Filename      :  NetlistImportTool.h
//
// Purpose       : Declare the interface to read and parse a netlist for
//                 an electrical circuit.
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef NetlistImportTool_H
#define NetlistImportTool_H

#include <string>
#include <vector>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_IO_CircuitMetadata.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OutputMgr.h>

// ---------- FORWARD DECLARATIONS -------

class NetlistImportToolData;

class N_IO_DistributionTool;

#ifdef Xyce_PARALLEL_MPI
class N_PDS_Comm;
#endif

class N_IO_CmdParse;
class N_IO_PkgOptionsMgr;

namespace Xyce
{
  namespace Topology
  {
    class Manager;
  }
}

//-----------------------------------------------------------------------------
// Class         : NetlistImportTool
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 07/28/00 ?
//-----------------------------------------------------------------------------
class NetlistImportTool
{
public:

  // Factory to generate singleton of class
  static NetlistImportTool * factory( N_IO_CmdParse & cp, Xyce::Topology::Manager & tm);

  // Registers N_PDS_Comm object for parallel communication and sets up
  // numProc_ and procID_
#ifdef Xyce_PARALLEL_MPI
  bool registerParallelServices(N_PDS_Comm * tmp_pds_ptr);
#endif

  // Destructor
  ~NetlistImportTool();

  // R Result
  // R- The result of constructing a circuit from a netlist.
  // I netlistFile
  // I- The file containing the netlist describing the circuit to be
  // I- constructed.
  // This function performs three basic steps. First, it reads and parses a
  // netlist from a specified file and stores the parsed fields in a structure
  // as strings. Second, it interprets the parsed data and stores the data in
  // an appropriate structure. Finally, it builds the Xyce Circuit.
  int constructCircuitFromNetlist(std::string const & netlistFile,
                                  const vector< pair< string, string> > & externalNetlistParams );

  bool registerDevMgr(N_DEV_DeviceInterface * devPtr);
  bool registerOutputMgr(N_IO_OutputMgr * outputMgrPtr);
  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

protected:

private:

  // constructor
  NetlistImportTool(N_IO_CmdParse & cp, Xyce::Topology::Manager & tm);

  N_IO_CircuitBlock * circuitBlock_;
  N_DEV_DeviceInterface* devIntPtr_;
  N_IO_OutputMgr * outputMgrPtr_;
  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

  bool isSerial_;

  N_IO_CmdParse & commandLine_;
  Xyce::Topology::Manager & topMgr_;

  /////////////////// Beginning of formerly static data //////////////
  list<N_IO_CircuitContext*> contextList_;
  N_IO_CircuitContext * currentContextPtr_;

  // These variables, etc., were originally in the
  // N_IO_CircuitBlock class, as static variables.  When
  // they were made non-static, they were moved up here.
  N_IO_DistributionTool* distToolPtr_;
  N_IO_CircuitMetadata metadata_;
  N_IO_CircuitContext circuitContext_;
  map<string,int> modelNames_;
  int useCount_;  // Counts copies so deletion of subcircuitList can
                  // be done properly in the class destructor.

  // This was originally static data in the N_IO_CircuitBlockData class.
  map<string,FileSSFPair> ssfMap_;

  // This was originally global data in the N_IO_CircuitBlock.C file.
  bool globalParamsInserted_;
  /////////////////////// End of formerly static data /////////////////

#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm * pdsCommPtr_;
#endif
};

#endif // NetlistImportTool_H
