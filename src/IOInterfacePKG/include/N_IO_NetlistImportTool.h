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
// Revision Number: $Revision: 1.38 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef NetlistImportTool_H
#define NetlistImportTool_H

#include <string>
#include <vector>

#include <N_PDS_ParallelMachine.h>
#include <N_PDS_Comm.h>

#include <N_IO_fwd.h>
#include <N_IO_CircuitMetadata.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {

namespace Topo {
class Manager;
}

namespace IO {

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
  static NetlistImportTool * factory( CmdParse & cp, Xyce::Topo::Manager & tm);

  // Registers N_PDS_Comm object for parallel communication and sets up
  // numProc_ and procID_
  bool registerParallelServices(N_PDS_Comm * tmp_pds_ptr);

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
                                  const std::vector< std::pair< std::string, std::string> > & externalNetlistParams );

  bool registerDevMgr(Device::DeviceInterface * devPtr);
  bool registerOutputMgr(OutputMgr * outputMgrPtr);
  // Method to register the package options manager
  bool registerPkgOptionsMgr( PkgOptionsMgr *pkgOptPtr );

protected:

private:

  // constructor
  NetlistImportTool(CmdParse & cp, Xyce::Topo::Manager & tm);

  CircuitBlock * circuitBlock_;
  Device::DeviceInterface* devIntPtr_;
  OutputMgr * outputMgrPtr_;
  // package options manager
  PkgOptionsMgr * pkgOptMgrPtr_;

  CmdParse & commandLine_;
  Xyce::Topo::Manager & topMgr_;

  std::list<CircuitContext*> contextList_;
  CircuitContext * currentContextPtr_;

  // These variables, etc., were originally in the
  // CircuitBlock class, as static variables.  When
  // they were made non-static, they were moved up here.
  DistributionTool* distToolPtr_;
  CircuitMetadata metadata_;
  CircuitContext circuitContext_;
  std::map<std::string,int> modelNames_;
  int useCount_;  // Counts copies so deletion of subcircuitList can
                  // be done properly in the class destructor.

  // This was originally static data in the CircuitBlockData class.
  std::map<std::string,FileSSFPair> ssfMap_;

  // This was originally global data in the CircuitBlock.C file.
  bool globalParamsInserted_;

  N_PDS_Comm * pdsCommPtr_;
  Parallel::Machine   comm_;
};

} // Namespace IO
} // namespace Xyce

#endif // NetlistImportTool_H
