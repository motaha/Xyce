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
// Filename       : $RCSfile: N_TOP_TopologyMgr.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/02/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_Topology_Manager_h
#define Xyce_Topology_Manager_h 1

// ---------- Standard Includes ----------

//#include <iosfwd>
#include <map>
#include <string>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------

class N_TOP_Topology;
class N_IO_PkgOptionsMgr;

/*
class N_PDS_Manager;
class N_ANP_AnalysisInterface;

//class N_LAS_QueryUtil;

namespace Xyce {

namespace Parallel {
 class PartitionTool;
}

namespace InputOutput {
 class RestartDataTool;
}
*/

class N_IO_CmdParse;

namespace Xyce {
namespace Topology {

// class DeviceIface;
 class InsertionTool;
// class NodeTool;
//

//-----------------------------------------------------------------------------
// Class         : Xyce::Topology::Manager
// Purpose       : Management of Topology System construction and usage
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/24/03
//-----------------------------------------------------------------------------
class Manager
{

 public:

  //typedefs
//  typedef map<string,System*> SystemMapType;

  // This used to be a singleton, but not anymore.  ERK. 1/30/2006
  static Manager * instance() { Manager *instPtr = new Manager(); return instPtr; }

  // Destructor
  ~Manager();

  N_TOP_Topology * createTopology(N_IO_CmdParse & cp);

  InsertionTool * getInsertionTool( const string & sysName );
  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

 private:

  Manager()
  : topo_(0), 
  currentDevInsertionTool(NULL)
  {}

  // Don't allow copy construction or assignment
  Manager( const Manager & );
  Manager & operator=( const Manager & );

  N_TOP_Topology * topo_;
  
  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

//  friend ostream & operator<<( ostream & os, const Manager & );

//  getInsertionTool creates an object and returns a pointer to it
//  so we need to track this so we can delete it.
//
  InsertionTool * currentDevInsertionTool;
};

} //namespace Topology
} //namespace Xyce

#endif //Xyce_Topology_Manager_h
